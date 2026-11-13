/* fprof.c: full path-cover profiles of tournament classes.
   input: classgen file (first token per line = code; bit k over pairs i<j in order, set => i->j), order n.
   output per class:  n code strong nsink nsource H C F_1 ... F_n
   F_j = number of ordered j-tuples of nonempty vertex-disjoint paths covering V (so Q(t) = sum_j F_j binom(t,j)),
   C = number of directed Hamiltonian cycles, nsink/nsource = numbers of detachable trailing/leading singletons
   (vertices removable one at a time as global sinks / sources). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
#define MAXN 10
static u64 Hs[1<<MAXN];            /* Hamiltonian paths of T[S], Hs[0]=0 here */
static u64 endp[(1<<MAXN)*MAXN];   /* paths covering S ending at v */
static u64 g[2][1<<MAXN];
int main(int argc,char**argv){
  if(argc<3){fprintf(stderr,"usage: fprof clsfile n\n");return 1;}
  FILE*f=fopen(argv[1],"r"); int n=atoi(argv[2]); char line[8192];
  unsigned N=1u<<n, FULL=N-1;
  while(fgets(line,sizeof line,f)){
    unsigned long long code; if(sscanf(line,"%llu",&code)!=1) continue;
    unsigned out[MAXN]={0}; int k=0;
    for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){ if(code>>k&1) out[i]|=1u<<j; else out[j]|=1u<<i; k++; }
    /* Hamiltonian paths of every subset */
    Hs[0]=0;
    for(unsigned S=1;S<N;S++){ u64 tot=0;
      for(int v=0;v<n;v++) if(S>>v&1){ unsigned P=S&~(1u<<v); u64 a=0;
        if(!P) a=1; else for(int w=0;w<n;w++) if((P>>w&1)&&(out[w]>>v&1)) a+=endp[P*n+w];
        endp[S*n+v]=a; tot+=a; }
      Hs[S]=tot; }
    /* Hamiltonian cycles: paths from 0 covering all ending at v with v->0, each cycle counted once */
    u64 C=0;
    if(n>=3){ /* count paths starting at vertex 0: recompute with fixed start */
      static u64 st0[(1<<MAXN)*MAXN];
      memset(st0,0,sizeof(u64)*N*n);
      st0[1u*n+0]=1;
      for(unsigned S=1;S<N;S++){ if(!(S&1)) continue;
        for(int v=0;v<n;v++){ u64 a=st0[S*n+v]; if(!a) continue;
          for(int w=0;w<n;w++) if(!(S>>w&1)&&(out[v]>>w&1)) st0[(S|1u<<w)*n+w]+=a; } }
      for(int v=1;v<n;v++) if(out[v]&1u) C+=st0[FULL*n+v];
    }
    /* strong connectivity */
    int strong=1;
    { unsigned seen=1, fr=1; while(fr){ unsigned nf=0; for(int v=0;v<n;v++) if(fr>>v&1) nf|=out[v]; nf&=~seen; seen|=nf; fr=nf; }
      if(seen!=FULL) strong=0;
      unsigned in[MAXN]={0}; for(int v=0;v<n;v++) for(int w=0;w<n;w++) if(out[v]>>w&1) in[w]|=1u<<v;
      seen=1; fr=1; while(fr){ unsigned nf=0; for(int v=0;v<n;v++) if(fr>>v&1) nf|=in[v]; nf&=~seen; seen|=nf; fr=nf; }
      if(seen!=FULL) strong=0; }
    /* detachable sinks / sources */
    int nsink=0,nsrc=0;
    { unsigned S=FULL; for(;;){ int found=-1; for(int v=0;v<n;v++) if((S>>v&1)&&!(out[v]&S&~(1u<<v))&&__builtin_popcount(S)>1) { found=v; break; }
        if(found<0) break; S&=~(1u<<found); nsink++; } }
    { unsigned S=FULL; for(;;){ int found=-1; for(int v=0;v<n;v++) if((S>>v&1)&&((out[v]&S)==(S&~(1u<<v)))&&__builtin_popcount(S)>1){ found=v; break; }
        if(found<0) break; S&=~(1u<<found); nsrc++; } }
    /* profile */
    printf("%d %llu %d %d %d %llu %llu",n,code,strong,nsink,nsrc,(unsigned long long)Hs[FULL],(unsigned long long)C);
    for(unsigned S=0;S<N;S++) g[0][S]=Hs[S];
    printf(" %llu",(unsigned long long)g[0][FULL]);
    int cur=0;
    for(int j=2;j<=n;j++){ int nx=cur^1;
      for(unsigned S=0;S<N;S++){ u64 s=0;
        for(unsigned B=S; B; B=(B-1)&S){ unsigned A=S&~B; if(!A) continue; s+=g[cur][A]*Hs[B]; }
        g[nx][S]=s; }
      cur=nx; printf(" %llu",(unsigned long long)g[cur][FULL]); }
    printf("\n");
  }
  return 0;
}
