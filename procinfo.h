#ifndef PROCINFO_H
#define PROCINFO_H

#define PROCINFO_MAX 64

#define PI_UNUSED   0
#define PI_SLEEPING 1
#define PI_RUNNABLE 2
#define PI_RUNNING  3
#define PI_ZOMBIE   4

struct procinfo {
  int pid;
  int ctime;
  int rtime;
  int state;
  char name[16];
};

#endif
