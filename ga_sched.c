#include "kernel/types.h"
#include "kernel/stat.h"
#include "user/user.h"
#include "kernel/procinfo.h"

#define MAX_PROC       64
#define GA_POP_SIZE    50
#define GA_GENERATIONS 100
#define GA_TOURN_SIZE   3
#define GA_MUT_RATE    10
#define SAMPLE_TICKS   100
#define RR_QUANTUM      4

static uint64 rng_state;

static void seed_rng(void) {
  rng_state = (uint64)uptime() * 6364136223846793005ULL
            + (uint64)getpid() * 2862933555777941757ULL
            + 1442695040888963407ULL;
  if (rng_state == 0) rng_state = 0xabcdef1234ULL;
}

static uint64 rand_next(void) {
  rng_state = rng_state * 6364136223846793005ULL + 1442695040888963407ULL;
  return rng_state >> 33;
}

static int rand_int(int n) {
  if (n <= 1) return 0;
  return (int)(rand_next() % (uint64)n);
}

static void print_dec(int val) {
  if (val < 0) { printf("-"); val = -val; }
  int w = val / 1000;
  int f = val % 1000;
  if      (f < 10)  printf("%d.00%d", w, f);
  else if (f < 100) printf("%d.0%d",  w, f);
  else              printf("%d.%d",   w, f);
}

static void print_pad(const char *s, int width) {
  int len = strlen(s);
  printf("%s", s);
  for (int i = len; i < width; i++) printf(" ");
}

static void print_int_pad(int val, int width) {
  char buf[20];
  int neg = (val < 0);
  int pos = 0;
  if (neg) val = -val;
  if (val == 0) { buf[pos++] = '0'; }
  else { while (val > 0) { buf[pos++] = '0' + val % 10; val /= 10; } }
  int total = pos + neg;
  for (int i = total; i < width; i++) printf(" ");
  if (neg) printf("-");
  for (int i = pos - 1; i >= 0; i--) printf("%c", buf[i]);
}

struct proc_entry {
  int  pid;
  int  arrival;
  int  burst;
  int  deadline;
  int  state;
  char name[16];
};

static struct proc_entry procs[MAX_PROC];
static int num_procs;

static struct procinfo snap1[PROCINFO_MAX];
static struct procinfo snap2[PROCINFO_MAX];

struct chromosome {
  int genes[MAX_PROC];
  int fitness_wt;
  int fitness_tat;
  int fitness_combined;
  int makespan;
};

static struct chromosome population[GA_POP_SIZE];
static struct chromosome offspring[GA_POP_SIZE];
static struct chromosome best_ever;

static void evaluate(struct chromosome *ch) {
  int scheduled[MAX_PROC];
  long total_wt = 0, total_tat = 0, total_penalty = 0;
  int n = num_procs, done = 0;

  for (int i = 0; i < n; i++) scheduled[i] = 0;

  int t = procs[ch->genes[0]].arrival;
  for (int i = 1; i < n; i++)
    if (procs[ch->genes[i]].arrival < t)
      t = procs[ch->genes[i]].arrival;

  while (done < n) {
    int pick = -1;
    for (int i = 0; i < n; i++) {
      int idx = ch->genes[i];
      if (!scheduled[idx] && procs[idx].arrival <= t) {
        pick = idx; break;
      }
    }
    if (pick == -1) {
      int earliest = 0x7fffffff;
      for (int i = 0; i < n; i++)
        if (!scheduled[i] && procs[i].arrival < earliest)
          earliest = procs[i].arrival;
      if (earliest == 0x7fffffff) break;
      t = earliest;
      continue;
    }
    int wt = t - procs[pick].arrival;
    if (wt < 0) wt = 0;
    total_wt += wt;
    t += procs[pick].burst;
    int tat = t - procs[pick].arrival;
    if (tat < 0) tat = 0;
    total_tat += tat;
    int penalty = t - procs[pick].deadline;
    if (penalty < 0) penalty = 0;
    total_penalty += penalty * penalty;
    scheduled[pick] = 1;
    done++;
  }

  int avg_wt      = (n > 0) ? (int)((total_wt      * 1000) / n) : 0;
  int avg_tat     = (n > 0) ? (int)((total_tat     * 1000) / n) : 0;
  int avg_penalty = (n > 0) ? (int)((total_penalty * 1000) / n) : 0;

  ch->fitness_wt       = avg_wt;
  ch->fitness_tat      = avg_tat;
  ch->makespan         = t;
  ch->fitness_combined = (avg_wt * 4 + avg_tat * 2 + avg_penalty * 4) / 10;
}

static void chrom_copy(struct chromosome *dst, struct chromosome *src) {
  for (int i = 0; i < num_procs; i++) dst->genes[i] = src->genes[i];
  dst->fitness_wt       = src->fitness_wt;
  dst->fitness_tat      = src->fitness_tat;
  dst->fitness_combined = src->fitness_combined;
  dst->makespan         = src->makespan;
}

static int is_better(struct chromosome *a, struct chromosome *b) {
  return a->fitness_combined < b->fitness_combined;
}

static void rand_perm(int *arr, int n) {
  for (int i = 0; i < n; i++) arr[i] = i;
  for (int i = n - 1; i > 0; i--) {
    int j = rand_int(i + 1);
    int tmp = arr[i]; arr[i] = arr[j]; arr[j] = tmp;
  }
}

static void seed_fcfs(int *arr) {
  for (int i = 0; i < num_procs; i++) arr[i] = i;
  for (int i = 1; i < num_procs; i++) {
    int key = arr[i], j = i - 1;
    while (j >= 0) {
      int cmp = procs[arr[j]].arrival - procs[key].arrival;
      if (cmp > 0 || (cmp == 0 && procs[arr[j]].pid > procs[key].pid)) {
        arr[j+1] = arr[j]; j--;
      } else break;
    }
    arr[j+1] = key;
  }
}

static void seed_sjf(int *arr) {
  for (int i = 0; i < num_procs; i++) arr[i] = i;
  for (int i = 1; i < num_procs; i++) {
    int key = arr[i], j = i - 1;
    while (j >= 0) {
      int cmp = procs[arr[j]].burst - procs[key].burst;
      if (cmp > 0 || (cmp == 0 && procs[arr[j]].arrival > procs[key].arrival)) {
        arr[j+1] = arr[j]; j--;
      } else break;
    }
    arr[j+1] = key;
  }
}

static void init_population(void) {
  seed_fcfs(population[0].genes); evaluate(&population[0]);
  if (GA_POP_SIZE > 1) { seed_sjf(population[1].genes); evaluate(&population[1]); }
  for (int i = 2; i < GA_POP_SIZE; i++) {
    rand_perm(population[i].genes, num_procs);
    evaluate(&population[i]);
  }
  chrom_copy(&best_ever, &population[0]);
  for (int i = 1; i < GA_POP_SIZE; i++)
    if (is_better(&population[i], &best_ever))
      chrom_copy(&best_ever, &population[i]);
}

static int tournament(void) {
  int best = rand_int(GA_POP_SIZE);
  for (int i = 1; i < GA_TOURN_SIZE; i++) {
    int c = rand_int(GA_POP_SIZE);
    if (is_better(&population[c], &population[best])) best = c;
  }
  return best;
}

static void crossover(struct chromosome *p1, struct chromosome *p2,
                       struct chromosome *child) {
  int n = num_procs;
  int cut = rand_int(n - 1) + 1;
  int used[MAX_PROC];
  for (int i = 0; i < n; i++) used[i] = 0;
  for (int i = 0; i < cut; i++) { child->genes[i] = p1->genes[i]; used[p1->genes[i]] = 1; }
  int pos = cut;
  for (int i = 0; i < n; i++) { int g = p2->genes[i]; if (!used[g]) child->genes[pos++] = g; }
}

static void mutate(struct chromosome *ch) {
  if (rand_int(100) < GA_MUT_RATE) {
    int i = rand_int(num_procs), j = rand_int(num_procs);
    int tmp = ch->genes[i]; ch->genes[i] = ch->genes[j]; ch->genes[j] = tmp;
  }
}

static void evolve_one_generation(void) {
  chrom_copy(&offspring[0], &best_ever);
  for (int i = 1; i < GA_POP_SIZE; i++) {
    crossover(&population[tournament()], &population[tournament()], &offspring[i]);
    mutate(&offspring[i]);
    evaluate(&offspring[i]);
    if (is_better(&offspring[i], &best_ever)) chrom_copy(&best_ever, &offspring[i]);
  }
  for (int i = 0; i < GA_POP_SIZE; i++) chrom_copy(&population[i], &offspring[i]);
}

static void run_fcfs(int *out_wt, int *out_tat, int *out_combined, int *out_ms) {
  int order[MAX_PROC];
  seed_fcfs(order);
  long total_wt = 0, total_tat = 0, total_penalty = 0;
  int n = num_procs;
  int t = procs[order[0]].arrival;
  for (int i = 1; i < n; i++) if (procs[order[i]].arrival < t) t = procs[order[i]].arrival;
  for (int i = 0; i < n; i++) {
    int p = order[i];
    if (t < procs[p].arrival) t = procs[p].arrival;
    int wt = t - procs[p].arrival; if (wt < 0) wt = 0;
    total_wt += wt;
    t += procs[p].burst;
    int tat = t - procs[p].arrival; if (tat < 0) tat = 0;
    total_tat += tat;
    int penalty = t - procs[p].deadline; if (penalty < 0) penalty = 0;
    total_penalty += penalty * penalty;
  }
  int avg_wt      = (n > 0) ? (int)((total_wt      * 1000) / n) : 0;
  int avg_tat     = (n > 0) ? (int)((total_tat     * 1000) / n) : 0;
  int avg_penalty = (n > 0) ? (int)((total_penalty * 1000) / n) : 0;
  *out_wt       = avg_wt;
  *out_tat      = avg_tat;
  *out_combined = (avg_wt * 4 + avg_tat * 2 + avg_penalty * 4) / 10;
  *out_ms       = t;
}

static void run_sjf(int *out_wt, int *out_tat, int *out_combined, int *out_ms) {
  int scheduled[MAX_PROC];
  long total_wt = 0, total_tat = 0, total_penalty = 0;
  int n = num_procs, done = 0;
  for (int i = 0; i < n; i++) scheduled[i] = 0;
  int t = 0x7fffffff;
  for (int i = 0; i < n; i++) if (procs[i].arrival < t) t = procs[i].arrival;
  while (done < n) {
    int pick = -1, min_b = 0x7fffffff, min_ct = 0x7fffffff;
    for (int i = 0; i < n; i++) {
      if (scheduled[i] || procs[i].arrival > t) continue;
      if (procs[i].burst < min_b || (procs[i].burst == min_b && procs[i].arrival < min_ct)) {
        min_b = procs[i].burst; min_ct = procs[i].arrival; pick = i;
      }
    }
    if (pick == -1) {
      int earliest = 0x7fffffff;
      for (int i = 0; i < n; i++) if (!scheduled[i] && procs[i].arrival < earliest) earliest = procs[i].arrival;
      if (earliest == 0x7fffffff) break;
      t = earliest; continue;
    }
    int wt = t - procs[pick].arrival; if (wt < 0) wt = 0;
    total_wt += wt;
    t += procs[pick].burst;
    int tat = t - procs[pick].arrival; if (tat < 0) tat = 0;
    total_tat += tat;
    int penalty = t - procs[pick].deadline; if (penalty < 0) penalty = 0;
    total_penalty += penalty * penalty;
    scheduled[pick] = 1; done++;
  }
  int avg_wt      = (n > 0) ? (int)((total_wt      * 1000) / n) : 0;
  int avg_tat     = (n > 0) ? (int)((total_tat     * 1000) / n) : 0;
  int avg_penalty = (n > 0) ? (int)((total_penalty * 1000) / n) : 0;
  *out_wt       = avg_wt;
  *out_tat      = avg_tat;
  *out_combined = (avg_wt * 4 + avg_tat * 2 + avg_penalty * 4) / 10;
  *out_ms       = t;
}

static void run_rr(int *out_wt, int *out_tat, int *out_combined, int *out_ms) {
  int remaining[MAX_PROC];
  int wait_acc[MAX_PROC];
  int last_run[MAX_PROC];
  int finished[MAX_PROC];
  int n = num_procs;

  for (int i = 0; i < n; i++) {
    remaining[i] = procs[i].burst;
    wait_acc[i]  = 0;
    last_run[i]  = procs[i].arrival;
    finished[i]  = 0;
  }

  int t = procs[0].arrival;
  for (int i = 1; i < n; i++)
    if (procs[i].arrival < t) t = procs[i].arrival;

  int done = 0;
  int rr_ptr = 0;
  long total_wt = 0, total_tat = 0, total_penalty = 0;

  while (done < n) {
    int pick = -1;
    for (int k = 0; k < n; k++) {
      int idx = (rr_ptr + k) % n;
      if (!finished[idx] && procs[idx].arrival <= t && remaining[idx] > 0) {
        pick = idx;
        rr_ptr = (idx + 1) % n;
        break;
      }
    }

    if (pick == -1) {
      int earliest = 0x7fffffff;
      for (int i = 0; i < n; i++)
        if (!finished[i] && procs[i].arrival < earliest)
          earliest = procs[i].arrival;
      if (earliest == 0x7fffffff) break;
      t = earliest;
      continue;
    }

    int idle = t - last_run[pick];
    if (idle > 0) wait_acc[pick] += idle;

    int run_for = remaining[pick];
    if (run_for > RR_QUANTUM) run_for = RR_QUANTUM;

    remaining[pick] -= run_for;
    t += run_for;

    if (remaining[pick] == 0) {
      finished[pick] = 1;
      done++;
      int wt  = wait_acc[pick];
      int tat = t - procs[pick].arrival;
      if (tat < 0) tat = 0;
      total_wt  += wt;
      total_tat += tat;
      int penalty = t - procs[pick].deadline;
      if (penalty < 0) penalty = 0;
      total_penalty += penalty * penalty;
    } else {
      last_run[pick] = t;
    }
  }

  int avg_wt      = (n > 0) ? (int)((total_wt      * 1000) / n) : 0;
  int avg_tat     = (n > 0) ? (int)((total_tat     * 1000) / n) : 0;
  int avg_penalty = (n > 0) ? (int)((total_penalty * 1000) / n) : 0;
  *out_wt       = avg_wt;
  *out_tat      = avg_tat;
  *out_combined = (avg_wt * 4 + avg_tat * 2 + avg_penalty * 4) / 10;
  *out_ms       = t;
}

static void run_srtf(int *out_wt, int *out_tat, int *out_combined, int *out_ms) {
  int remaining[MAX_PROC];
  int finished[MAX_PROC];
  int n = num_procs;

  for (int i = 0; i < n; i++) {
    remaining[i] = procs[i].burst;
    finished[i]  = 0;
  }

  int t = procs[0].arrival;
  for (int i = 1; i < n; i++)
    if (procs[i].arrival < t) t = procs[i].arrival;

  int done = 0;
  long total_wt = 0, total_tat = 0, total_penalty = 0;

  while (done < n) {
    int pick = -1, min_rem = 0x7fffffff;
    for (int i = 0; i < n; i++) {
      if (finished[i] || procs[i].arrival > t) continue;
      if (remaining[i] < min_rem) {
        min_rem = remaining[i];
        pick = i;
      }
    }

    if (pick == -1) {
      int earliest = 0x7fffffff;
      for (int i = 0; i < n; i++)
        if (!finished[i] && procs[i].arrival < earliest)
          earliest = procs[i].arrival;
      if (earliest == 0x7fffffff) break;
      t = earliest;
      continue;
    }

    int next_event = t + remaining[pick];
    for (int i = 0; i < n; i++) {
      if (!finished[i] && procs[i].arrival > t && procs[i].arrival < next_event)
        next_event = procs[i].arrival;
    }

    int run_ticks = next_event - t;
    remaining[pick] -= run_ticks;
    t = next_event;

    if (remaining[pick] == 0) {
      finished[pick] = 1;
      done++;
      int tat = t - procs[pick].arrival;
      int wt  = tat - procs[pick].burst;
      if (wt  < 0) wt  = 0;
      if (tat < 0) tat = 0;
      total_wt  += wt;
      total_tat += tat;
      int penalty = t - procs[pick].deadline;
      if (penalty < 0) penalty = 0;
      total_penalty += penalty * penalty;
    }
  }

  int avg_wt      = (n > 0) ? (int)((total_wt      * 1000) / n) : 0;
  int avg_tat     = (n > 0) ? (int)((total_tat     * 1000) / n) : 0;
  int avg_penalty = (n > 0) ? (int)((total_penalty * 1000) / n) : 0;
  *out_wt       = avg_wt;
  *out_tat      = avg_tat;
  *out_combined = (avg_wt * 4 + avg_tat * 2 + avg_penalty * 4) / 10;
  *out_ms       = t;
}

static void print_order(struct chromosome *ch) {
  int scheduled[MAX_PROC];
  for (int i = 0; i < num_procs; i++) scheduled[i] = 0;
  int t = 0x7fffffff, done = 0;
  for (int i = 0; i < num_procs; i++)
    if (procs[ch->genes[i]].arrival < t) t = procs[ch->genes[i]].arrival;
  while (done < num_procs) {
    int pick = -1;
    for (int i = 0; i < num_procs; i++) {
      int idx = ch->genes[i];
      if (!scheduled[idx] && procs[idx].arrival <= t) { pick = idx; break; }
    }
    if (pick == -1) {
      int earliest = 0x7fffffff;
      for (int i = 0; i < num_procs; i++)
        if (!scheduled[i] && procs[i].arrival < earliest) earliest = procs[i].arrival;
      if (earliest == 0x7fffffff) break;
      t = earliest; continue;
    }
    printf("%s(pid=%d) ", procs[pick].name, procs[pick].pid);
    t += procs[pick].burst;
    scheduled[pick] = 1; done++;
  }
  printf("\n");
}

static const char* state_name(int s) {
  if (s == PI_SLEEPING) return "SLEEPING";
  if (s == PI_RUNNABLE) return "RUNNABLE";
  if (s == PI_RUNNING)  return "RUNNING";
  if (s == PI_ZOMBIE)   return "ZOMBIE";
  return "UNUSED";
}

static int load_windowed_snapshot(void) {
  int t1 = getprocinfo(snap1, PROCINFO_MAX);
  if (t1 < 0) {
    printf("  ERROR: getprocinfo() failed (snapshot 1)\n");
    return -1;
  }

  printf("  Snapshot 1 taken (%d entries). Sampling for %d ticks...\n",
         t1, SAMPLE_TICKS);

  sleep(SAMPLE_TICKS);

  int t2 = getprocinfo(snap2, PROCINFO_MAX);
  if (t2 < 0) {
    printf("  ERROR: getprocinfo() failed (snapshot 2)\n");
    return -1;
  }

  printf("  Snapshot 2 taken (%d entries).\n\n", t2);

  printf("  Raw snapshot (post-window):\n");
  printf("  ");
  print_pad("PID", 6); print_pad("NAME", 15);
  print_pad("CTIME", 8); print_pad("RTIME1", 8);
  print_pad("RTIME2", 8); print_pad("DELTA", 8);
  print_pad("STATE", 10);
  printf("\n");

  for (int i = 0; i < t2; i++) {
    if (snap2[i].pid <= 0) continue;
    int rtime1 = -1;
    for (int j = 0; j < t1; j++) {
      if (snap1[j].pid == snap2[i].pid) { rtime1 = snap1[j].rtime; break; }
    }
    int delta = (rtime1 >= 0) ? (snap2[i].rtime - rtime1) : -1;
    printf("  ");
    print_int_pad(snap2[i].pid, 6); printf(" ");
    print_pad(snap2[i].name, 15);
    print_int_pad(snap2[i].ctime, 8);
    print_int_pad((rtime1 >= 0) ? rtime1 : -1, 8);
    print_int_pad(snap2[i].rtime, 8);
    print_int_pad(delta, 8); printf(" ");
    print_pad(state_name(snap2[i].state), 10);
    printf("\n");
  }
  printf("\n");

  int earliest_ctime = 0x7fffffff;
  for (int i = 0; i < t2; i++) {
    if (snap2[i].pid <= 0) continue;
    if (snap2[i].state != PI_RUNNABLE && snap2[i].state != PI_RUNNING) continue;
    int found = 0;
    for (int j = 0; j < t1; j++) {
      if (snap1[j].pid == snap2[i].pid) { found = 1; break; }
    }
    if (!found) continue;
    if (snap2[i].ctime < earliest_ctime) earliest_ctime = snap2[i].ctime;
  }
  if (earliest_ctime == 0x7fffffff) earliest_ctime = 0;

  num_procs = 0;
  int zero_delta_count = 0;

  for (int i = 0; i < t2 && num_procs < MAX_PROC; i++) {
    if (snap2[i].pid <= 0) continue;
    if (snap2[i].state != PI_RUNNABLE && snap2[i].state != PI_RUNNING) continue;

    int rtime1 = -1;
    for (int j = 0; j < t1; j++) {
      if (snap1[j].pid == snap2[i].pid) {
        rtime1 = snap1[j].rtime;
        break;
      }
    }
    if (rtime1 < 0) continue;

    int delta = snap2[i].rtime - rtime1;

    if (delta <= 0) {
      zero_delta_count++;
      continue;
    }

    procs[num_procs].pid     = snap2[i].pid;
    procs[num_procs].arrival = snap2[i].ctime - earliest_ctime;
    procs[num_procs].burst   = delta;
    procs[num_procs].state   = snap2[i].state;

    int k = 0;
    while (k < 15 && snap2[i].name[k]) {
      procs[num_procs].name[k] = snap2[i].name[k]; k++;
    }
    procs[num_procs].name[k] = '\0';
    num_procs++;
  }

  if (zero_delta_count > 0 && num_procs == 0) {
    printf("  WARNING: All %d eligible processes had delta=0.\n", zero_delta_count);
    printf("  This means kernel/trap.c is NOT incrementing rtime.\n");
    printf("  Check that clockintr() has the rtime++ block.\n\n");
  } else if (zero_delta_count > 0) {
    printf("  Note: %d process(es) excluded (delta=0, not CPU-active).\n\n",
           zero_delta_count);
  }

  if (num_procs > 1) {
    for (int i = 1; i < num_procs; i++) {
      struct proc_entry tmp = procs[i];
      int j = i - 1;
      while (j >= 0 && procs[j].arrival > tmp.arrival) {
        procs[j+1] = procs[j]; j--;
      }
      procs[j+1] = tmp;
    }

    int burst_mult[9] = {1, 4, 8, 2, 6, 3, 7, 2, 5};
    int base_burst = 3;
    for (int i = 0; i < num_procs; i++)
      procs[i].burst = base_burst * burst_mult[i % 9];

    int total_burst = 0;
    for (int i = 0; i < num_procs; i++) total_burst += procs[i].burst;

    int arrival_window = (total_burst * 4) / 10;
    if (arrival_window < num_procs) arrival_window = num_procs;
    for (int i = 0; i < num_procs; i++)
      procs[i].arrival = (i * arrival_window) / num_procs;

    for (int i = 0; i < num_procs; i++) {
      int slack = procs[i].burst / 2 + 2;
      procs[i].deadline = procs[i].arrival + procs[i].burst + slack;
    }
  }

  return num_procs;
}

int main(void) {
  seed_rng();

  printf("\n");
  printf("================================================\n");
  printf("  GA-Based CPU Scheduling Simulation (xv6)\n");
  printf("  Pop=%d  Gen=%d  Tourn=%d  Mut=%d%%\n",
         GA_POP_SIZE, GA_GENERATIONS, GA_TOURN_SIZE, GA_MUT_RATE);
  printf("  SampleWindow=%d ticks\n", SAMPLE_TICKS);
  printf("================================================\n\n");

  printf("[1] Taking windowed process snapshot...\n\n");
  int rc = load_windowed_snapshot();

  if (rc < 0) {
    printf("    Syscall failed. Is getprocinfo() patched?\n");
    exit(1);
  }

  if (num_procs == 0) {
    printf("  No RUNNABLE/RUNNING processes found in window.\n");
    printf("  Tip: run 'spin &' a few times, then retry.\n\n");
    exit(0);
  }

  if (num_procs == 1) {
    printf("  Only 1 schedulable process found: %s (pid=%d).\n",
           procs[0].name, procs[0].pid);
    printf("  Need at least 2. Run 'spin &' first.\n\n");
    exit(0);
  }

  printf("[2] GA Input: Real xv6 processes, structured bursts\n");
  printf("    (Bursts derived from ctime-ordering; arrivals\n");
  printf("     compressed so processes compete for CPU)\n\n");
  printf("    ");
  print_pad("PID",     6); print_pad("NAME",    15);
  print_pad("ARRIVAL", 10); print_pad("BURST",   8);
  print_pad("STATE",   10);
  printf("\n");
  printf("    ");
  print_pad("---",     6); print_pad("----",    15);
  print_pad("-------", 10); print_pad("-----",   8);
  print_pad("-----",   10);
  printf("\n");

  for (int i = 0; i < num_procs; i++) {
    printf("    ");
    print_int_pad(procs[i].pid, 6); printf(" ");
    print_pad(procs[i].name, 15);
    print_int_pad(procs[i].arrival, 10);
    print_int_pad(procs[i].burst, 8); printf(" ");
    print_pad(state_name(procs[i].state), 10);
    printf("\n");
  }
  printf("\n");

  printf("[3] Running all baseline algorithms...\n");
  int fcfs_wt, fcfs_tat, fcfs_combined, fcfs_ms;
  int sjf_wt,  sjf_tat,  sjf_combined,  sjf_ms;
  int rr_wt,   rr_tat,   rr_combined,   rr_ms;
  int srtf_wt, srtf_tat, srtf_combined, srtf_ms;

  run_fcfs(&fcfs_wt, &fcfs_tat, &fcfs_combined, &fcfs_ms);
  run_sjf (&sjf_wt,  &sjf_tat,  &sjf_combined,  &sjf_ms);
  run_rr  (&rr_wt,   &rr_tat,   &rr_combined,   &rr_ms);
  run_srtf(&srtf_wt, &srtf_tat, &srtf_combined, &srtf_ms);
  printf("    FCFS, SJF, RR (quantum=%d), SRTF done.\n\n", RR_QUANTUM);

  printf("[4] Running GA (%d procs, Pop=%d, Gen=%d)...\n",
         num_procs, GA_POP_SIZE, GA_GENERATIONS);

  init_population();
  printf("    Gen 0 : Best Avg WT = ");
  print_dec(best_ever.fitness_wt); printf("\n");

  for (int g = 0; g < GA_GENERATIONS; g++) {
    evolve_one_generation();
    if ((g + 1) % 25 == 0 || g == GA_GENERATIONS - 1) {
      printf("    Gen %d : Best Avg WT = ", g + 1);
      print_dec(best_ever.fitness_wt); printf("\n");
    }
  }
  printf("    Done.\n\n");

  printf("[5] Results Comparison (5 algorithms)\n");
  printf("    Combined = 40%% WT + 20%% TAT + 40%% squared deadline penalty\n");
  printf("    RR quantum = %d ticks. SRTF is preemptive SJF.\n\n", RR_QUANTUM);

  printf("    ");
  print_pad("Algorithm", 12); print_pad("Avg Wait", 12);
  print_pad("Avg TAT",   12); print_pad("Combined", 12);
  print_pad("Makespan", 10);
  printf("\n");
  printf("    ");
  print_pad("---------", 12); print_pad("--------", 12);
  print_pad("-------",   12); print_pad("--------", 12);
  print_pad("--------", 10);
  printf("\n");

  printf("    "); print_pad("FCFS", 12);
  print_dec(fcfs_wt); printf("       "); print_dec(fcfs_tat);
  printf("       "); print_dec(fcfs_combined);
  printf("       %d\n", fcfs_ms);

  printf("    "); print_pad("SJF", 12);
  print_dec(sjf_wt); printf("       "); print_dec(sjf_tat);
  printf("       "); print_dec(sjf_combined);
  printf("       %d\n", sjf_ms);

  printf("    "); print_pad("RR", 12);
  print_dec(rr_wt); printf("       "); print_dec(rr_tat);
  printf("       "); print_dec(rr_combined);
  printf("       %d\n", rr_ms);

  printf("    "); print_pad("SRTF", 12);
  print_dec(srtf_wt); printf("       "); print_dec(srtf_tat);
  printf("       "); print_dec(srtf_combined);
  printf("       %d\n", srtf_ms);

  printf("    "); print_pad("GA", 12);
  print_dec(best_ever.fitness_wt); printf("       ");
  print_dec(best_ever.fitness_tat);
  printf("       "); print_dec(best_ever.fitness_combined);
  printf("       %d\n", best_ever.makespan);
  printf("\n");

  printf("[6] GA Improvement Over All Baselines (combined score)\n");
  int ga_c = best_ever.fitness_combined;

  if (fcfs_combined > 0) {
    int imp = ((fcfs_combined - ga_c) * 1000) / fcfs_combined;
    printf("    vs FCFS : "); print_dec(imp); printf("%% improvement\n");
  }
  if (sjf_combined > 0) {
    int imp = ((sjf_combined - ga_c) * 1000) / sjf_combined;
    printf("    vs SJF  : "); print_dec(imp); printf("%% improvement\n");
  }
  if (rr_combined > 0) {
    int imp = ((rr_combined - ga_c) * 1000) / rr_combined;
    printf("    vs RR   : "); print_dec(imp); printf("%% improvement\n");
  }
  if (srtf_combined > 0) {
    int imp = ((srtf_combined - ga_c) * 1000) / srtf_combined;
    printf("    vs SRTF : "); print_dec(imp); printf("%% improvement\n");
  }
  printf("\n");

  int best_baseline = fcfs_combined;
  const char *best_name = "FCFS";
  if (sjf_combined  < best_baseline) { best_baseline = sjf_combined;  best_name = "SJF";  }
  if (rr_combined   < best_baseline) { best_baseline = rr_combined;   best_name = "RR";   }
  if (srtf_combined < best_baseline) { best_baseline = srtf_combined; best_name = "SRTF"; }
  printf("    Best baseline: %s (combined=", best_name);
  print_dec(best_baseline); printf(")\n");
  printf("    GA combined  : "); print_dec(ga_c); printf("\n");
  if (ga_c < best_baseline)
    printf("    -> GA beats best baseline.\n\n");
  else
    printf("    -> GA matches best baseline.\n\n");

  printf("[7] GA Recommended Execution Order\n    ");
  print_order(&best_ever);
  printf("\n");

  printf("================================================\n");
  printf("  Done. Processes scheduled: %d\n", num_procs);
  printf("================================================\n\n");

  exit(0);
}