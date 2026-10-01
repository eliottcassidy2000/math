// procgen_edim_20261001_sad.c -- descending simulated annealing for edge-multiset resolving sets of Q_d.
// Additive 64-bit Zobrist keys key(e)=sum_s R[d(e,s)] (false collisions only make it conservative);
// every reported set is re-verified exactly (full histogram comparison) before printing.
// usage: sad d k0 iters_per_level seed T0 T1 kmin
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
static int D,N,E; static uint8_t *dist; static uint64_t R[16]; static uint64_t rs;
static inline uint64_t rng(void){rs^=rs<<13;rs^=rs>>7;rs^=rs<<17;return rs;}
static inline double ur(void){return (rng()>>11)*(1.0/9007199254740992.0);}
static uint64_t *hk; static int *hc,*hu; static int HS,nu;
static long cost(const uint64_t*key){ long c=0; nu=0;
  for(int e=0;e<E;e++){ uint64_t k=key[e]; int h=(int)((k*0x9E3779B97F4A7C15ULL)>>40)&(HS-1);
    while(1){ if(hc[h]==0){hk[h]=k;hc[h]=1;hu[nu++]=h;break;} if(hk[h]==k){c+=hc[h];hc[h]++;break;} h=(h+1)&(HS-1);} }
  for(int i=0;i<nu;i++)hc[hu[i]]=0; return c; }
static int cmp16(const void*a,const void*b){return memcmp(a,b,16);}
static int exact_ok(int*S,int k){ // exact check: sort full histograms
  uint8_t (*H)[16]=calloc(E,16); for(int e=0;e<E;e++){ for(int j=0;j<k;j++) H[e][dist[(size_t)e*N+S[j]]]++; }
  // counts may exceed 255 only if k>255
  qsort(H,E,16,cmp16); int ok=1; for(int e=1;e<E;e++) if(!memcmp(H[e],H[e-1],16)){ok=0;break;} free(H); return ok; }
int main(int argc,char**argv){ if(argc<8){fprintf(stderr,"usage\n");return 1;}
  D=atoi(argv[1]); int K=atoi(argv[2]); long iters=atol(argv[3]); rs=strtoull(argv[4],0,10)*0x9E3779B97F4A7C15ULL+1;
  double T0=atof(argv[5]),T1=atof(argv[6]); int kmin=atoi(argv[7]);
  N=1<<D; E=D*(N/2); dist=malloc((size_t)E*N); HS=1; while(HS<4*E)HS<<=1; hk=malloc(8*HS);hc=calloc(HS,4);hu=malloc(4*HS);
  for(int r=0;r<16;r++) R[r]=rng()|1;
  int e=0; for(int u=0;u<N;u++)for(int i=0;i<D;i++)if(!(u>>i&1)){ for(int w=0;w<N;w++) dist[(size_t)e*N+w]=__builtin_popcount((u^w)&~(1<<i)); e++; }
  uint64_t *key=malloc(8*E),*nk=malloc(8*E); int *S=malloc(4*N),*out=malloc(4*N),*in=calloc(N,4);
  int c=0; while(c<K){int v=rng()%N; if(!in[v]){in[v]=1;S[c++]=v;}}
  int k=K;
  while(k>=kmin){
    int no=0; for(int v=0;v<N;v++) if(!in[v]) out[no++]=v;
    for(int f=0;f<E;f++){uint64_t x=0; for(int j=0;j<k;j++) x+=R[dist[(size_t)f*N+S[j]]]; key[f]=x;}
    long cur=cost(key); long it;
    for(it=0; it<iters && cur>0; it++){
      double T=T0*pow(T1/T0,(double)it/iters); int a=rng()%k,b=rng()%no,s=S[a],t=out[b];
      for(int f=0;f<E;f++) nk[f]=key[f]-R[dist[(size_t)f*N+s]]+R[dist[(size_t)f*N+t]];
      long nc=cost(nk);
      if(nc<=cur || ur()<exp(-(nc-cur)/T)){ uint64_t*tmp=key;key=nk;nk=tmp; cur=nc; S[a]=t; out[b]=s; in[s]=0; in[t]=1; }
    }
    if(cur>0){ printf("d=%d k=%d FAILED (cost %ld after %ld iters)\n",D,k,cur,it); break; }
    if(!exact_ok(S,k)){ printf("d=%d k=%d hash-ok but exact FAIL?!\n",D,k); break; }
    printf("d=%d k=%d FOUND iters=%ld set:",D,k,it); for(int j=0;j<k;j++)printf(" %d",S[j]); printf("\n"); fflush(stdout);
    // remove the vertex whose removal gives least cost
    long bestc=-1; int besta=0;
    for(int a=0;a<k;a++){ for(int f=0;f<E;f++) nk[f]=key[f]-R[dist[(size_t)f*N+S[a]]]; long cc=cost(nk); if(bestc<0||cc<bestc){bestc=cc;besta=a;} }
    in[S[besta]]=0; S[besta]=S[k-1]; k--;
  }
  return 0; }
