// Binary inheritance in the interval-hub family: is every odd N' obtained from a realizer of its odd parent
// n (N' = 2n-1 or 2n+1) by inserting one spine vertex (any position, any x-bit, interval kept an interval)?
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define NS 36
typedef unsigned long long u64;
static int idx(int o,int a,int b,int d){return ((o*3+a)*3+b)*2+d;}
static int T[4][NS][3], Tn[4][NS], fin[2][NS];
static void build(void){
  for(int let=0;let<4;let++){ int xi=let>>1, yi=let&1;
    for(int o=0;o<2;o++)for(int a=0;a<3;a++)for(int b=0;b<3;b++)for(int d=0;d<2;d++){ int s=idx(o,a,b,d), k=0;
      int h1=o?yi:xi, h2=o?xi:yi;
      T[let][s][k++]=idx(o,h1?2:1,b,d);
      if(b||h1) T[let][s][k++]=idx(o,a,h2?2:1,d);
      if(d||h2) T[let][s][k++]=idx(o,a,b,1);
      Tn[let][s]=k; } }
  for(int e=0;e<2;e++) for(int o=0;o<2;o++)for(int a=0;a<3;a++)for(int b=0;b<3;b++)for(int d=0;d<2;d++){
      int h12=o?!e:e; fin[e][idx(o,a,b,d)] = !(a==2||b==2||(b==0&&!h12)); } }
static u64 LIM; static unsigned char *inh, *cov;
int main(int argc,char**argv){
  int M=atoi(argv[1]); LIM=strtoull(argv[2],0,10); build();
  inh=calloc(LIM/2+2,1); cov=calloc(LIM/2+2,1);
  static u64 P[40][NS], S[40][NS]; int let[40];
  for(int m=0;m<=M;m++) for(u64 X=0; X<(1ULL<<m); X++) for(int p=0;p<=m;p++) for(int q=0;p+q<=m;q++){
    if(q==0 && p>0) continue;   // empty interval: one representative
    for(int i=0;i<m;i++) let[i]=(int)((((X>>i)&1)<<1) | (i>=p && i<p+q));
    memset(P[0],0,sizeof P[0]); P[0][idx(0,0,0,0)]=1; P[0][idx(1,0,0,0)]=1;
    for(int i=0;i<m;i++){ memset(P[i+1],0,sizeof P[i+1]); for(int s=0;s<NS;s++) if(P[i][s]) for(int t=0;t<Tn[let[i]][s];t++) P[i+1][T[let[i]][s][t]]+=P[i][s]; }
    for(int e=0;e<2;e++){
      for(int s=0;s<NS;s++) S[m][s]=fin[e][s];
      for(int i=m-1;i>=0;i--) for(int s=0;s<NS;s++){ u64 qq=0; for(int t=0;t<Tn[let[i]][s];t++) qq+=S[i+1][T[let[i]][s][t]]; S[i][s]=qq; }
      u64 n=0; for(int s=0;s<NS;s++) n+=P[0][s]*S[0][s];
      if(n<=LIM) cov[n>>1]=1;
      if(m==M || 2*n+1>LIM) continue;
      // insertion at position i (before old vertex i), new x-bit b, membership y in I allowed if the interval stays an interval
      for(int i=0;i<=m;i++) for(int y=0;y<2;y++){
        int ok;
        if(q==0) ok = 1;                              // empty interval: new vertex may become the interval or stay outside
        else if(y==1) ok = (i>=p && i<=p+q);          // inside: must touch/extend the interval
        else ok = (i<=p || i>=p+q);                   // outside: must not split the interval
        if(!ok) continue;
        for(int b=0;b<2;b++){ int l=(b<<1)|y; u64 h=0;
          for(int s=0;s<NS;s++) if(P[i][s]){ u64 qq=0; for(int t=0;t<Tn[l][s];t++) qq+=S[i][T[l][s][t]]; h+=P[i][s]*qq; }
          if(h==2*n-1 || h==2*n+1) inh[h>>1]=1; }
      }
    }
  }
  u64 X=23; while(X<=LIM && cov[X>>1]) X+=2;
  printf("M=%d: coverage horizon %llu; odd N' in [45, %llu) not inherited from their odd parent:",M,X,X);
  int shown=0; u64 cnt=0, tot=0;
  for(u64 N=45; N<X; N+=2){ tot++; if(!inh[N>>1]){ cnt++; if(shown<40){printf(" %llu",N);shown++;} } }
  printf("\n  not inherited: %llu of %llu\n",cnt,tot);
  return 0;
}
