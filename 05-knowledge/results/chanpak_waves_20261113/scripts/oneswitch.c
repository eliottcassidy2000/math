// Sub-family of two-hub words: hub-1 word x arbitrary, hub-2 word y with one switch (y = 0^p 1^(m-p) or 1^p 0^(m-p)).
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
int main(int argc,char**argv){
  int M=atoi(argv[1]); u64 LIM=strtoull(argv[2],0,10); int fam=atoi(argv[3]); // fam 0: y=0^p1^q ; 1: y=1^p0^q ; 2: both
  build(); unsigned char *cov=calloc(LIM/2+2,1);
  static u64 P[40][NS], S[40][NS];
  for(int m=1;m<=M;m++) for(u64 X=0; X<(1ULL<<m); X++){
    for(int dir=0; dir<2; dir++){ if(fam!=2 && dir!=fam) continue;
      int ylo = dir==0?0:1, yhi = dir==0?1:0;   // y = ylo on [0,p), yhi on [p,m)
      memset(P[0],0,sizeof P[0]); P[0][idx(0,0,0,0)]=1; P[0][idx(1,0,0,0)]=1;
      for(int i=0;i<m;i++){ int let=(int)(((X>>i)&1)<<1)|ylo; memset(P[i+1],0,sizeof P[i+1]); for(int s=0;s<NS;s++) if(P[i][s]) for(int t=0;t<Tn[let][s];t++) P[i+1][T[let][s][t]]+=P[i][s]; }
      for(int e=0;e<2;e++){
        for(int s=0;s<NS;s++) S[m][s]=fin[e][s];
        for(int i=m-1;i>=0;i--){ int let=(int)(((X>>i)&1)<<1)|yhi; for(int s=0;s<NS;s++){ u64 q=0; for(int t=0;t<Tn[let][s];t++) q+=S[i+1][T[let][s][t]]; S[i][s]=q; } }
        for(int p=0;p<=m;p++){ u64 h=0; for(int s=0;s<NS;s++) h+=P[p][s]*S[p][s]; if(h<=LIM) cov[h>>1]=1; }
      }
    }
  }
  int shown=0; u64 cnt=0, first=0; printf("fam %d, m<=%d: misses among odd <= %llu:",fam,M,LIM);
  for(u64 n=1;n<=LIM;n+=2) if(!cov[n>>1]){ cnt++; if(shown<40){printf(" %llu",n);shown++;} }
  printf("\n  total misses %llu of %llu\n",cnt,(LIM+1)/2);
  return 0;
}
