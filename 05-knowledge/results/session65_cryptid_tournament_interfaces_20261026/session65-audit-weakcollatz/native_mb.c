/* Independently authored audit simulator, not copied from tinycryptids.
 * Table: Pascal Michel, Behavior of Busy Beavers, MB winner subsection.
 * No external code is loaded or executed. Bounds and checks survive NDEBUG. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum { CELLS=200001, ORIGIN=100000, HALT=5 };
static unsigned char tape[CELLS];
typedef struct { int write, move, next; } Move;
static const Move table[5][2] = {
  {{1,1,1},{1,-1,2}}, {{1,1,2},{1,1,1}},
  {{1,1,3},{0,-1,4}}, {{1,-1,0},{1,-1,3}},
  {{1,1,HALT},{0,-1,0}}
};
static void check(int ok, const char *reason) {
  if (!ok) { fprintf(stderr,"FAIL: %s\n",reason); exit(1); }
}
static void step(int *state, int *head) {
  check(*state>=0 && *state<5,"state before step");
  check(*head>0 && *head<CELLS-1,"head before step");
  Move t=table[*state][tape[*head]];
  tape[*head]=(unsigned char)t.write;
  *head+=t.move;
  *state=t.next;
}
static int is_config(int state, int head, int n) {
  if (state!=0 || head+n>=CELLS) return 0;
  for (int p=0;p<CELLS;p++)
    if (tape[p]!=(unsigned char)(p>head && p<=head+n)) return 0;
  return 1;
}
static int halt_shape(int state, int head, int k) {
  if (state!=HALT || head+3*k+5>=CELLS || head<1) return 0;
  for (int p=0;p<CELLS;p++) {
    int d=p-head;
    int expected=(d==-1 || d==1 || d==3*k+5 ||
      (d>=4 && d<=3*k+4 && (d-4)%3==0));
    if(tape[p]!=(unsigned char)expected) return 0;
  }
  return 1;
}
int main(void) {
  const int orbit[]={0,6,16,34,64,114,196,334,564,946,1584,2646,4416,7366,12284};
  int state=0,head=ORIGIN;
  long long total=0;
  printf("{\"macro_checkpoints\":[");
  for(int i=0;i<15;i++) {
    int n=orbit[i],k=n/3,r=n%3;
    check(is_config(state,head,n),"native C(n) at macro start");
    long long cost = r==0 ? 5LL*k*k+19*k+15 : r==1 ? 5LL*k*k+25*k+27 : 6LL*k+12;
    for(long long j=0;j<cost;j++) step(&state,&head);
    total+=cost;
    if(r!=2) check(is_config(state,head,orbit[i+1]),"native macro endpoint");
    else check(halt_shape(state,head,k),"native full halting tape pattern");
    printf("%s{\"n\":%d,\"cost\":%lld,\"cumulative\":%lld}",i?",":"",n,cost,total);
  }
  int ones=0,lo=CELLS,hi=-1;
  for(int p=0;p<CELLS;p++) if(tape[p]) {ones++;if(p<lo)lo=p;hi=p;}
  check(state==HALT && total==47176870 && ones==4098,"blank run final controls");
  printf("],\"steps\":%lld,\"ones\":%d,\"nonzero_span\":%d,\"head_relative_to_origin\":%d,",total,ones,hi-lo+1,head-ORIGIN);
  /* Exhaustive finite macro tests, separate fresh native tapes. */
  int cases=0;
  for(int n=0;n<=120;n++) {
    memset(tape,0,sizeof tape); state=0;head=ORIGIN;
    for(int j=1;j<=n;j++) tape[head+j]=1;
    int k=n/3,r=n%3;
    long long cost=r==0?5LL*k*k+19*k+15:r==1?5LL*k*k+25*k+27:6LL*k+12;
    for(long long j=0;j<cost;j++) step(&state,&head);
    check(r==2?halt_shape(state,head,k):is_config(state,head,5*k+(r==0?6:9)),"finite macro test");
    cases++;
  }
  printf("\"fresh_native_macro_cases\":%d,\"status\":\"PASS\"}\n",cases);
  return 0;
}
