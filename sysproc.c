#include "types.h"
#include "riscv.h"
#include "defs.h"
#include "param.h"
#include "memlayout.h"
#include "spinlock.h"
#include "vm.h"
#include "proc.h"
#include "procinfo.h"

extern struct proc proc[];

uint64 sys_getprocinfo(void)
{
  uint64 buf_addr;
  int max;

  argaddr(0, &buf_addr);
  argint(1, &max);

  if (max <= 0)
    return 0;
  if (max > PROCINFO_MAX)
    max = PROCINFO_MAX;

  struct procinfo kbuf[PROCINFO_MAX];
  int count = 0;

  for (struct proc *p = proc; p < &proc[NPROC] && count < max; p++) {
    acquire(&p->lock);

    if (p->state == UNUSED) {
      release(&p->lock);
      continue;
    }

    kbuf[count].pid   = p->pid;
    kbuf[count].ctime = p->ctime;
    kbuf[count].rtime = p->rtime;
    kbuf[count].state = (int)p->state;

    int k = 0;
    while (k < 15 && p->name[k]) {
      kbuf[count].name[k] = p->name[k];
      k++;
    }
    kbuf[count].name[k] = '\0';

    release(&p->lock);
    count++;
  }

  if (copyout(myproc()->pagetable,
              buf_addr,
              (char *)kbuf,
              count * sizeof(struct procinfo)) < 0)
    return -1;

  return count;
}

uint64
sys_exit(void)
{
  int n;
  argint(0, &n);
  kexit(n);
  return 0;
}

uint64
sys_getpid(void)
{
  return myproc()->pid;
}

uint64
sys_fork(void)
{
  return kfork();
}

uint64
sys_wait(void)
{
  uint64 p;
  argaddr(0, &p);
  return kwait(p);
}

uint64
sys_sbrk(void)
{
  uint64 addr;
  int t;
  int n;

  argint(0, &n);
  argint(1, &t);
  addr = myproc()->sz;

  if(t == SBRK_EAGER || n < 0) {
    if(growproc(n) < 0) {
      return -1;
    }
  } else {
    if(addr + n < addr)
      return -1;
    if(addr + n > TRAPFRAME)
      return -1;
    myproc()->sz += n;
  }
  return addr;
}

uint64
sys_pause(void)
{
  int n;
  uint ticks0;

  argint(0, &n);
  if(n < 0)
    n = 0;
  acquire(&tickslock);
  ticks0 = ticks;
  while(ticks - ticks0 < n){
    if(killed(myproc())){
      release(&tickslock);
      return -1;
    }
    sleep(&ticks, &tickslock);
  }
  release(&tickslock);
  return 0;
}

uint64
sys_kill(void)
{
  int pid;

  argint(0, &pid);
  return kkill(pid);
}

uint64
sys_uptime(void)
{
  uint xticks;

  acquire(&tickslock);
  xticks = ticks;
  release(&tickslock);
  return xticks;
}

uint64
sys_sleep(void)
{
  int n;
  argint(0, &n);

  acquire(&tickslock);

  uint ticks0 = ticks;
  while(ticks - ticks0 < n){
    if(killed(myproc())){
      release(&tickslock);
      return -1;
    }
    sleep(&ticks, &tickslock);
  }

  release(&tickslock);
  return 0;
}
