// Per-value ternary fertility: for each odd n, do some realizer w (two-hub word, length <= LMAX) and some single-letter
// insertion (any position, any letter) give H(w') = 3n + d for each d in {0,2,4}?  Also binary: H(w') = 2n +- 1.
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
  int LMAX=atoi(argv[1]); u64 NMAX=strtoull(argv[2],0,10); build();
  unsigned char *f3=calloc(NMAX/2+2,1), *f2=calloc(NMAX/2+2,1), *seen=calloc(NMAX/2+2,1);
  int w[32]; u64 P[33][NS], S[33][NS];
  for(int m=1;m<=LMAX;m++){
    long nw=1L<<(2*m);
    for(long code=0;code<nw;code++){
      long c=code; for(int i=0;i<m;i++){ w[i]=c&3; c>>=2; }
      memset(P[0],0,sizeof P[0]); P[0][idx(0,0,0,0)]=1; P[0][idx(1,0,0,0)]=1;
      for(int i=0;i<m;i++){ memset(P[i+1],0,sizeof P[i+1]); for(int s=0;s<NS;s++) if(P[i][s]) for(int t=0;t<Tn[w[i]][s];t++) P[i+1][T[w[i]][s][t]]+=P[i][s]; }
      for(int e=0;e<2;e++){
        for(int s=0;s<NS;s++) S[m][s]=fin[e][s];
        for(int i=m-1;i>=0;i--) for(int s=0;s<NS;s++){ u64 q=0; for(int t=0;t<Tn[w[i]][s];t++) q+=S[i+1][T[w[i]][s][t]]; S[i][s]=q; }
        u64 n=0; for(int s=0;s<NS;s++) n+=P[0][s]*S[0][s];
        if(n>NMAX) continue; seen[n>>1]=1;
        for(int p=0;p<=m;p++) for(int let=0;let<4;let++){
          u64 h=0; for(int s=0;s<NS;s++) if(P[p][s]) { u64 q=0; for(int t=0;t<Tn[let][s];t++) q+=S[p][T[let][s][t]]; h+=P[p][s]*q; }
          if(h>=3*n && h<=3*n+4) f3[n>>1]|=1<<((h-3*n)/2);
          if(h==2*n-1) f2[n>>1]|=1; if(h==2*n+1) f2[n>>1]|=2;
        }
      }
    }
  }
  // report
  u64 X=23; while(X<=NMAX && seen[X>>1]) X+=2;
  long nf3=0, nf2=0, tot=0; int shown=0;
  printf("LMAX=%d: words of length <= %d realize every odd n in [23,%llu)\n",LMAX,LMAX,X);
  printf("ternary-infertile odd n in [23,%llu) (missing children 3n+d, d in {0,2,4}):",X/3);
  for(u64 n=23;n<X/3;n+=2){ tot++; if(f3[n>>1]==7) nf3++; else if(shown<40){ printf(" %llu[%d]",n,f3[n>>1]); shown++; } if(f2[n>>1]==3) nf2++; }
  printf("\n  ternary fertile: %ld / %ld;  binary (2n-1 and 2n+1) fertile: %ld / %ld\n",nf3,tot,nf2,tot);
  return 0;
}
