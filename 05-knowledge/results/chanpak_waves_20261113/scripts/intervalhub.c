// Two-hub words with x arbitrary and y having at most two switches: y = A^p B^q A^r (A,B in {0,1}, A != B).
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
  int M=atoi(argv[1]); u64 LIM=strtoull(argv[2],0,10); int A=atoi(argv[3]);   // outer y value A, inner 1-A
  build(); unsigned char *cov=calloc(LIM/2+2,1);
  static u64 P[40][NS], Q[40][NS], S[40][NS];
  for(int m=1;m<=M;m++){ for(u64 X=0; X<(1ULL<<m); X++){
    // prefix vectors with y=A (outer) for all positions; with y=1-A (inner)
    memset(P[0],0,sizeof P[0]); P[0][idx(0,0,0,0)]=1; P[0][idx(1,0,0,0)]=1;
    for(int i=0;i<m;i++){ int let=(int)(((X>>i)&1)<<1)|A; memset(P[i+1],0,sizeof P[i+1]); for(int s=0;s<NS;s++) if(P[i][s]) for(int t=0;t<Tn[let][s];t++) P[i+1][T[let][s][t]]+=P[i][s]; }
    for(int e=0;e<2;e++){
      // suffix responses with outer letters
      for(int s=0;s<NS;s++) S[m][s]=fin[e][s];
      for(int i=m-1;i>=0;i--){ int let=(int)(((X>>i)&1)<<1)|A; for(int s=0;s<NS;s++){ u64 q=0; for(int t=0;t<Tn[let][s];t++) q+=S[i+1][T[let][s][t]]; S[i][s]=q; } }
      for(int p=0;p<=m;p++){
        // inner block starting at p: propagate Q = P[p] through inner letters for q = 0..m-p
        u64 v[NS]; memcpy(v,P[p],sizeof v);
        for(int q=0;p+q<=m;q++){
          u64 h=0; for(int s=0;s<NS;s++) h+=v[s]*S[p+q][s]; if(h<=LIM) cov[h>>1]=1;
          if(p+q<m){ int let=(int)(((X>>(p+q))&1)<<1)|(1-A); u64 w[NS]; memset(w,0,sizeof w); for(int s=0;s<NS;s++) if(v[s]) for(int t=0;t<Tn[let][s];t++) w[T[let][s][t]]+=v[s]; memcpy(v,w,sizeof v); }
        }
      }
    }
  }
  { u64 Xh=23; while(Xh<=LIM && cov[Xh>>1]) Xh+=2; printf("m<=%d: horizon %llu  (2^m=%llu, ratio %.3f, /3^m %.4f)\n",m,Xh,1ULL<<m,(double)Xh/(1ULL<<m),(double)Xh/__builtin_powi(3.0,m)); fflush(stdout);} }
  int shown=0; u64 cnt=0; printf("y = %d^p %d^q %d^r, m<=%d: misses among odd <= %llu:",A,1-A,A,M,LIM);
  for(u64 n=1;n<=LIM;n+=2) if(!cov[n>>1]){ cnt++; if(shown<30){printf(" %llu",n);shown++;} }
  printf("\n  total misses %llu of %llu\n",cnt,(LIM+1)/2);
  return 0;
}
