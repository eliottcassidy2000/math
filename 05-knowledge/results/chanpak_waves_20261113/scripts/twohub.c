// Two-hub tournaments: spine TT_m (i->j for i<j), hubs u1,u2 with x_i=[u1->i], y_i=[u2->i], e=[u1->u2].
// H = sum over hub orders (h1,h2) of #assignments of spine vertices to increasing segments S0,S1,S2 of the
// path S0 h1 S1 h2 S2 with junction arcs.  Incremental DFS over (x_i,y_i); marks odd H <= LIM.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
static u64 LIM; static int MMAX; static unsigned char *cov;
// state: s0 in {0 none,1 last ok(beats? we store label),...}: we track label of last S0 elem wrt h1: 0 none, 1 last beats h1 (bad if final), 2 last beaten by h1
// simpler: s0: 0 none, 1 last->h1 ok, 2 last not->h1 ; s1: 0 none, 1 started last->h2 ok, 2 started last not->h2 ; s2: 0 none,1 started
typedef struct { u64 c[2][3][3][2]; } St;   // [order][s0][s1][s2]
static void step(const St*in, St*out, int xi, int yi){
  memset(out,0,sizeof(St));
  for(int o=0;o<2;o++){
    // order 0: h1=u1,h2=u2 ; order 1: h1=u2,h2=u1
    int h1beats = o==0? xi: yi;      // h1 -> i ?
    int h2beats = o==0? yi: xi;      // h2 -> i ?
    for(int s0=0;s0<3;s0++)for(int s1=0;s1<3;s1++)for(int s2=0;s2<2;s2++){
      u64 c=in->c[o][s0][s1][s2]; if(!c) continue;
      // put i in S0 (only if S1,S2 not started? no: segments are value-increasing independently; positions interleave freely)
      out->c[o][h1beats?2:1][s1][s2]+=c;            // i becomes last of S0; i->h1 iff h1 does not beat i
      // put i in S1: if S1 not started, need h1 -> i
      if(s1==0){ if(h1beats) out->c[o][s0][h2beats?2:1][s2]+=c; }
      else out->c[o][s0][h2beats?2:1][s2]+=c;
      // put i in S2: if not started need h2 -> i
      if(s2==0){ if(h2beats) out->c[o][s0][s1][1]+=c; }
      else out->c[o][s0][s1][1]+=c;
    }
  }
}
static u64 finish(const St*s, int e){
  u64 H=0;
  for(int o=0;o<2;o++){
    int h1toh2 = o==0? e: !e;
    for(int s0=0;s0<3;s0++)for(int s1=0;s1<3;s1++)for(int s2=0;s2<2;s2++){
      u64 c=s->c[o][s0][s1][s2]; if(!c) continue;
      if(s0==2) continue;                 // last of S0 must beat h1
      if(s1==2) continue;                 // last of S1 must beat h2
      if(s1==0 && !h1toh2) continue;      // S1 empty: need h1 -> h2
      H+=c;
    }
  }
  return H;
}
static long long nodes=0;
static void dfs(int m, St*s){
  nodes++;
  for(int e=0;e<2;e++){ u64 H=finish(s,e); if(H<=LIM && (H&1)) cov[H>>1]=1; }
  if(m>=MMAX) return;
  // prune: if minimum possible H already exceeds LIM? H is not monotone in general; rely on MMAX
  St t;
  for(int xi=0;xi<2;xi++)for(int yi=0;yi<2;yi++){ step(s,&t,xi,yi); dfs(m+1,&t); }
}
int main(int argc,char**argv){
  LIM=strtoull(argv[1],0,10); MMAX=atoi(argv[2]);
  if(argc>3){ // single evaluation: x bits, y bits, e
    char*X=argv[3],*Y=argv[4]; int e=atoi(argv[5]); int m=strlen(X);
    St s; memset(&s,0,sizeof s); s.c[0][0][0][0]=1; s.c[1][0][0][0]=1;
    for(int i=0;i<m;i++){ St t; step(&s,&t,X[i]-'0',Y[i]-'0'); s=t; }
    printf("%llu\n",finish(&s,e)); return 0; }
  cov=calloc(LIM/2+2,1);
  St s; memset(&s,0,sizeof s); s.c[0][0][0][0]=1; s.c[1][0][0][0]=1;
  dfs(0,&s);
  long long miss=0; u64 first=0; int shown=0;
  printf("nodes %lld\nmisses:",nodes);
  for(u64 N=1;N<=LIM;N+=2) if(!cov[N>>1]){ if(!first||1) {} miss++; if(shown<25){printf(" %llu",N);shown++;} }
  printf("\ntotal misses %lld of %llu odd <= %llu\n",miss,(LIM+1)/2,LIM);
  return 0;
}
