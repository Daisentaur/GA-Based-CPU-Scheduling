# GA-Based CPU Scheduler for xv6

A Genetic Algorithm CPU scheduler implemented inside xv6, a minimal RISC-V teaching OS. The GA runs against actual live kernel processes — not simulated data — using a multi-objective fitness function that combines waiting time, turnaround time, and a squared deadline-penalty term.

```
Algorithm   Avg Wait   Avg TAT   Combined Cost   Makespan
GA            27.92     40.15       508.92          160
SJF           26.31     38.54       524.45          160   (-2.9% vs GA)
SRTF          24.85     37.08       542.89          160   (-6.3% vs GA)
FCFS          44.77     57.00       942.45          160  (-46.0% vs GA)
RR (q=4)      56.23     68.46      1413.08          160  (-64.0% vs GA)
```
*13-process workload run inside xv6/QEMU. Lower combined cost = better.*

## Why this exists

Every paper we found on GA-based CPU scheduling runs on simulated process data with hardcoded arrival and burst times. We wanted to know what happens when you run a GA against processes that actually exist in a kernel — with real arrival ticks, real CPU consumption, and the constraints that come with integer-only arithmetic and a 4KB kernel stack.

The short answer: it works, and it beats every classical baseline on a multi-objective metric that accounts for deadline compliance. The slightly longer answer is in [our paper](paper/FinalReportOs.pdf).

## What it does, end to end

1. Extends the xv6 kernel with two new fields in `struct proc` (`ctime`, `rtime`) and a new system call (`getprocinfo`) that snapshots the live process table into user space.
2. `ga_sched` calls `getprocinfo` twice, 100 ticks apart. The burst proxy for each process is the delta `rtime₂ − rtime₁`. Because xv6's round-robin scheduler distributes CPU time almost perfectly evenly, a deterministic multiplier pattern is applied to produce a realistic mix of short, medium, and long burst estimates.
3. Arrival times are compressed into a window so processes genuinely compete for the CPU.
4. Per-process soft deadlines are assigned: `Dᵢ = ATᵢ + BTᵢ + ⌊BTᵢ/2⌋ + 2`. Short jobs get tight deadlines; long jobs get slightly more slack.
5. The GA evaluates chromosomes (permutations of process indices) by simulating the non-preemptive schedule forward in time and computing `Cost = 0.4·WT + 0.2·TAT + 0.4·L²`, where `L²` is average squared lateness.
6. The population is seeded with FCFS and SJF orderings so the GA can never start worse than both classical baselines.
7. Tournament selection, single-point OX1 crossover with repair, swap mutation, and elitism evolve the population for 100 generations.
8. Results are printed as a comparison table: FCFS, SJF, RR, SRTF, and GA side by side.

## Why the squared penalty matters

SJF is provably optimal for average waiting time. Any GA using only waiting time as its fitness criterion will converge to SJF at best and can't beat it. The `L²` term is what gives the GA real headroom.

Under a linear lateness penalty, being 100 ticks late costs the same per-tick as being 1 tick late. Under a squared penalty, a job already 50 ticks past its deadline costs 100 units per additional tick (derivative of `x²` at 50 is 100). This makes the GA aggressively prevent catastrophic lateness on long jobs — which is exactly the failure mode SJF and SRTF both have.

Concretely: SRTF has the lowest raw waiting time of all five algorithms, yet its combined cost (542.89) is *higher* than the GA's (508.92). SRTF's aggressive preemption continuously postpones long jobs until every short job has run, causing those long jobs to miss their deadlines by huge margins. The GA learns to pull some long jobs forward, accepting a small wait-time penalty in exchange for a large reduction in squared lateness.

## Project structure and layout

```
ga-cpu-scheduler-xv6/
├── user/
│   ├── ga_sched.c         # main program: reads live kernel procs, runs GA
│   └── spin.c             # helper: infinite-loop process for testing
├── kernel/
│   ├── procinfo.h         # new: shared struct procinfo, state constants
│   └── sysproc.c          # modified: adds sys_getprocinfo()
├── paper/
│   └── FinalReportOs.pdf  # full IEEE-format paper
└── README.md
```

Files you also need to patch in your xv6 tree (changes described below):

```
kernel/proc.h              add ctime, rtime fields to struct proc
kernel/proc.c              initialise ctime = ticks, rtime = 0 in allocproc()
kernel/trap.c              increment rtime in clockintr()
kernel/syscall.h           assign syscall number 23 to getprocinfo
kernel/syscall.c           register the handler
user/usys.pl               add entry("getprocinfo")
user/user.h                add prototype
Makefile                   add $U/_ga_sched and $U/_spin to UPROGS
```

## Setup

### Prerequisites

- Ubuntu 22.04 or 24.04
- Python 3 (for the QEMU launch script)
- About 500MB of disk space for the toolchain

### Install dependencies

```bash
sudo apt-get install -y git qemu-system-misc gcc-riscv64-linux-gnu \
    binutils-riscv64-linux-gnu make
```

### Clone xv6

This repo ships only the modified and added files. You need the base xv6 tree:

```bash
git clone https://github.com/mit-pdos/xv6-riscv.git
cd xv6-riscv
```

### Apply kernel patches

**1. `kernel/proc.h`** — add two fields to `struct proc`, after the `pid` field:

```c
int ctime;   // tick at which process was created
int rtime;   // total ticks spent in RUNNING state
```

**2. `kernel/proc.c`** — in `allocproc()`, after `p->pid = allocpid()`:

```c
p->ctime = ticks;
p->rtime = 0;
```

**3. `kernel/trap.c`** — at the end of `clockintr()`, after `release(&tickslock)`:

```c
struct proc *p = myproc();
if (p && p->state == RUNNING)
    p->rtime++;
```

**4. `kernel/syscall.h`** — add:

```c
#define SYS_getprocinfo 23
```

**5. `kernel/syscall.c`** — in the `syscalls[]` array:

```c
[SYS_getprocinfo] sys_getprocinfo,
```

Also add the extern declaration near the top:

```c
extern uint64 sys_getprocinfo(void);
```

**6. `kernel/procinfo.h`** — copy from this repo's `kernel/procinfo.h`.

**7. `kernel/sysproc.c`** — copy from this repo's `kernel/sysproc.c` (replaces the original; adds `sys_getprocinfo` and includes `procinfo.h`).

**8. `user/usys.pl`** — add:

```perl
entry("getprocinfo");
```

**9. `user/user.h`** — add prototype (before the `uchar` section):

```c
int getprocinfo(struct procinfo*, int);
```

Also add the forward declaration for `struct procinfo` near the top:

```c
struct procinfo;
```

**10. Copy user programs:**

```bash
cp /path/to/this-repo/user/ga_sched.c user/
cp /path/to/this-repo/user/spin.c user/
```

**11. `Makefile`** — in the `UPROGS` list, add:

```
$U/_ga_sched\
$U/_spin\
```

### Build and run

```bash
make TOOLPREFIX=riscv64-linux-gnu- clean
make TOOLPREFIX=riscv64-linux-gnu- qemu
```

If you're on a system where the toolchain is just `riscv64-unknown-elf-`, use that prefix instead.

## Usage

Boot xv6, then spawn background processes so there's something to schedule:

```
$ spin &
$ spin &
$ spin &
$ spin &
$ spin &
$ ga_sched
```

Wait a second or two between launches so the processes get different `ctime` values. The 100-tick sampling window in `ga_sched` will capture the CPU deltas and produce differentiated burst estimates.

To exit QEMU: `Ctrl-A` then `X`.

## How it works in detail

### The getprocinfo syscall

`ga_sched` calls `getprocinfo(buf, max)` which fills a `struct procinfo[]` with a snapshot of the kernel process table. The kernel-side handler:

1. Iterates over `proc[NPROC]`.
2. Acquires each process's spinlock before reading any field — never held across `copyout()`.
3. Skips `UNUSED` slots.
4. Copies to a kernel buffer, then calls `copyout()` to transfer safely to user space.

The `procinfo` struct is minimal by design:

```c
struct procinfo {
    int pid;
    int ctime;   // creation tick
    int rtime;   // CPU ticks consumed
    int state;
    char name[16];
};
```

### Windowed sampling

A single snapshot of `rtime` is useless because xv6's round-robin scheduler gives every runnable process almost exactly the same CPU time per tick. `ga_sched` calls `getprocinfo` twice with a 100-tick sleep in between. The burst estimate for each process is `Δrtime = rtime₂ − rtime₁`. Processes not present in both snapshots (created or terminated mid-window) are excluded.

### Burst augmentation

Even with windowed sampling, the deltas are nearly uniform (6–7 ticks each) since xv6's scheduler is perfectly fair. A deterministic multiplier pattern indexed by each process's creation-time rank produces realistic burst diversity:

```
multipliers = [1, 4, 8, 2, 6, 3, 7, 2, 5, ...]
base unit   = 3 ticks
→ burst times: 3, 12, 24, 6, 18, 9, 21, 6, ...
```

Arrival times are then compressed into a window equal to 40% of total burst time, ensuring the CPU is never idle and processes genuinely compete.

### GA operators

| Component | Implementation |
|-----------|----------------|
| Chromosome | Permutation of process indices; position = scheduling priority |
| Fitness | Non-preemptive schedule simulation → `Cost = 0.4·WT + 0.2·TAT + 0.4·L²` |
| Initialisation | Slot 0 = FCFS order, Slot 1 = SJF order, Slots 2–49 = random |
| Selection | Tournament (k = 3) |
| Crossover | Single-point OX1 with repair (guarantees valid permutation) |
| Mutation | Swap mutation, probability = 0.10 |
| Elitism | `best_ever` always survives into the next generation |
| Complexity | O(P · G · n²); with P=50, G=100, n≤15: ~1.1M operations |

### Why no floating point

xv6 user-space doesn't support floating-point (no FPU context save/restore for user processes in the default build). All fitness values are stored as integers scaled by 1000, and a custom `print_dec()` function handles formatted output like `27.920`.

## Results

Representative run with 13 processes (arrivals 0–58 ticks, bursts 3–24 ticks):

| Algorithm | Avg WT | Avg TAT | Combined Cost | Makespan |
|-----------|--------|---------|---------------|----------|
| GA | 27.92 | 40.15 | **508.92** | 160 |
| SJF | 26.31 | 38.54 | 524.45 | 160 |
| SRTF | 24.85 | 37.08 | 542.89 | 160 |
| FCFS | 44.77 | 57.00 | 942.45 | 160 |
| RR (q=4) | 56.23 | 68.46 | 1413.08 | 160 |

GA recommended execution order:
```
P4(3) → P6(12) → P10(6) → P12(18) → P14(9) → P24(3) →
P20(6) → P22(15) → P32(6) → P28(12) → P8(24) → P18(21) → P30(24)
```
Numbers in parentheses are burst times. The GA inserts P8 (burst=24) before P18 and P30 rather than deferring it until last as SJF would — this is the deadline-protection trade-off in action.

## A note on xv6 and scale

xv6 was used here because of a course requirement, not because it's the right environment to showcase what this algorithm is actually built for. That's worth being explicit about.

A GA's core strength is navigating a combinatorially large search space that greedy algorithms can't meaningfully explore. With 5–13 processes, that search space is tiny — `13!` is only about 6 billion permutations, and greedy heuristics like SJF already get very close to optimal by luck of being the right kind of greedy for small inputs. There's simply not enough room for evolution to demonstrate an advantage on waiting time alone, which is why we introduced the squared deadline penalty to create a harder, multi-dimensional objective.

The seeding strategy — initialising the population with FCFS and SJF orderings — was specifically designed with large workloads in mind. On a small dataset, the GA converges in a handful of generations because the seeded solutions are already nearly optimal and there's nowhere better to go. On a large dataset (50–200 processes), those seeds act as a strong starting point that prevents wasted generations on random garbage, while the remaining 48 random chromosomes provide the diversity needed to explore meaningfully different regions of the schedule space. That's when the generational improvement curve actually looks like what a GA should produce — rapid early gains as good building blocks recombine, then gradual refinement as the population converges on something neither FCFS nor SJF could reach alone.

The real target environment for this algorithm is batch and HPC job scheduling: a queue with 50–200 processes, known or predictable burst times, heterogeneous arrival patterns, and per-job deadlines. In that setting, the combinatorial search space is genuinely intractable for greedy methods, burst-time estimates are available upfront from job submission metadata, and the multi-objective fitness function can be tuned to reflect real SLAs. xv6's process count and its perfectly fair round-robin scheduler are both fundamental mismatches with those conditions — we worked around them with burst augmentation, but that's a patch on the evaluation environment, not a limitation of the algorithm.

If you're building on this and want to see it perform as intended, run it against a synthetic workload of 100+ processes with varied burst times and tight per-process deadlines. That's what it was designed for.

## Limitations

**Burst augmentation is a workaround.** Because xv6's scheduler distributes CPU time almost uniformly, the burst estimates used by the GA are synthetically diversified rather than reflecting true process CPU requirements. A production-grade implementation would use exponential averaging over historical `rtime` samples to produce real burst predictions.

**RR and SRTF baselines are user-space approximations.** They use the same deadline formula and combined-cost metric as the GA for a fair comparison, but they don't capture real context-switch overhead. Raw waiting-time numbers would differ slightly from a true kernel-level implementation.

**No dynamic arrivals.** Processes that appear in only one of the two snapshots are discarded. A production scheduler would need to handle mid-window arrivals and terminations.

**xv6 has very few concurrent processes.** The live-kernel demo typically schedules 5–13 processes. The interesting scheduling complexity only appears when burst augmentation creates enough diversity. On a real Linux system with 50–100 concurrent runnable processes and meaningful burst-time differences, the GA's O(P·G·n²) cost remains tractable (P=50, G=100, n=100 → ~50M operations, sub-second on any modern CPU) and the combinatorial search space grows large enough to make the GA's evolutionary approach genuinely advantageous over greedy heuristics.

## Troubleshooting

| Symptom | Likely cause | Fix |
|---------|-------------|-----|
| `exec ga_sched failed` | Old filesystem image cached | `make clean && make qemu` |
| All processes show `rtime=0` | `clockintr()` patch not applied | Add `p->rtime++` inside `clockintr()` in `trap.c` |
| All processes show `ctime=0` | `allocproc()` patch not applied | Add `p->ctime = ticks` in `allocproc()` in `proc.c` |
| `getprocinfo() returned -1` | `copyout` failed; buffer too small or wrong address | Check `PROCINFO_MAX` in `procinfo.h` matches between kernel and user |
| `No RUNNABLE/RUNNING processes` | spin processes finished before `ga_sched` ran | Make `spin.c` loop forever: `while(1) {}` |
| Only 2 processes scheduled | Didn't launch spin processes | Run `spin &` several times before `ga_sched` |
| Compile error: `stray '@'` | Terminal prompt accidentally pasted into source | Remove the line starting with a username from `sysproc.c` |

## Paper

The full IEEE-format paper is in [`paper/FinalReportOs.pdf`](paper/FinalReportOs.pdf). It covers the formal problem definition, kernel extension details, GA design, the multi-objective fitness function, and complete experimental results.

## References

1. S. R. Sakhare and M. S. Ali, "Genetic Algorithm Based Adaptive Scheduling Algorithm for Real-Time Operating Systems," *IJESA*, Vol. 2, No. 3, 2012.
2. M. González-Rodríguez et al., "Study and evaluation of CPU scheduling algorithms," *Heliyon*, Elsevier, Vol. 10, No. 9, 2024.
3. Elsayed et al., "Optimizing multiprocessor performance in real-time systems using an innovative genetic algorithm approach," *Scientific Reports*, Nature, 2025.
4. D. E. Goldberg, *Genetic Algorithms in Search, Optimization, and Machine Learning*. Addison-Wesley, 1989.
5. R. Cox, F. Kaashoek, R. Morris, "xv6: a simple, Unix-like teaching operating system," MIT CSAIL, 2023.

## License

MIT
