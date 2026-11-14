// For each fixed (m, p, q, e), the value set of H(x; I, e) over x in {0,1}^m: longest run of consecutive odd values.
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
static int cmpu(const void*a,const void*b){ u64 x=*(u64*)a,y=*(u64*)b; return x<y?-1:x>y; }
int main(int argc,char**argv){
  int m=atoi(argv[1]); build(); long n=1L<<m; u64 *vals=malloc(sizeof(u64)*n*2);
  u64 bestlen=0; int bp=0,bq=0,be=0; u64 bA=0;
  for(int p=0;p<=m;p++) for(int q=0;p+q<=m;q++) for(int e=0;e<2;e++){
    long k=0;
    for(long X=0;X<n;X++){ u64 v[NS]; memset(v,0,sizeof v); v[idx(0,0,0,0)]=1; v[idx(1,0,0,0)]=1;
      for(int i=0;i<m;i++){ int let=(int)(((X>>i)&1)<<1)|(i>=p&&i<p+q); u64 w[NS]; memset(w,0,sizeof w);
        for(int s=0;s<NS;s++) if(v[s]) for(int t=0;t<Tn[let][s];t++) w[T[let][s][t]]+=v[s]; memcpy(v,w,sizeof v); }
      u64 h=0; for(int s=0;s<NS;s++) if(fin[e][s]) h+=v[s]; vals[k++]=h; }
    qsort(vals,k,8,cmpu);
    // longest run of consecutive odd values
    u64 runlen=1, runstart=vals[0];
    for(long i=1;i<k;i++){ if(vals[i]==vals[i-1]) continue; if(vals[i]==vals[i-1]+2){ runlen++; } else { if(runlen>bestlen){bestlen=runlen;bp=p;bq=q;be=e;bA=runstart;} runlen=1; runstart=vals[i]; } }
    if(runlen>bestlen){bestlen=runlen;bp=p;bq=q;be=e;bA=runstart;}
  }
  printf("m=%d: longest run of consecutive odd values from ONE interval: %llu values [%llu, %llu] at p=%d q=%d e=%d (2^m=%ld)\n",m,bestlen,bA,bA+2*(bestlen-1),bp,bq,be,n);
  return 0;
}
