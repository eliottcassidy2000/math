// collatz_procgen_20260923_sweep_tightdp.c
//
// Sparse "tight" dynamic programme for record loops and hub cycles of the E-graph
// (HYP-9122: loops through 1; HYP-9125: loops through -1, written in u = -v > 0).
//
// Reverse move (sign +1): x -> y = (2^k x - 1)/3, legal iff y is a positive integer, 3 does not divide y.
// Reverse move (sign -1): u -> y = (2^k u + 1)/3, legal iff y is a positive integer, 3 does not divide y.
//   (u = -v: the negative E-graph v -> 3v+1, v -> v/2 becomes u -> 3u-1, u -> u/2.)
//
// A reverse path from the start h with i moves and K halvings has carry B (2^K h = 3^i x + B for +,
// 3^i u = 2^K h + B for -), with B >= (3^i-1)/2.  Put E = B - (3^i-1)/2 >= 0 and G = E/(2^K h).
// Then along a path  G' = G + (1-2^-k) nu,  nu = (3^i-1)/(2^(K+1) h),  so G is nondecreasing.
//
// TIGHTNESS LEMMA (proved in the note): on every record object (a loop with the extremal number of
// halvings, or a cycle with (q, ceil(q log2 3)) resp. (n, ceil(n log3 2))), every prefix state (i,x)
// carries the extremal K, hence K is a function of (i,x), and G at (i,x) is <= G_final (the value at
// the end).  So a boolean DP over the states with G <= G_final is EXACT for existence; this window
// holds only ~2*G_final*XMAX values per layer.
//
// Usage:
//   tightdp SIGN H LAYERS KTARGET GMAX XMAX [CLIMB] [RECON_T] [CKPT_DIR]
//     SIGN   : +1 or -1
//     H      : start = end value (1 for loops, the hub for cycles)
//     LAYERS : number of reverse moves (multiplications)
//     KTARGET: required total number of halvings at the end
//     GMAX   : G_final (decimal); states with G > GMAX + margin are discarded
//     XMAX   : value cap
//     CLIMB  : 1 = report, for loops through 1, the largest N such that the climb point c_N is
//              visited (c_N = (3^(N+1)-1)/2 for +, (3^N+1)/2 for -) with the final K
//     RECON_T: 0 = no reconstruction; T > 0 = write checkpoints every T layers and reconstruct one
//              object (target, or the deepest climb point if CLIMB) -> prints the move list k_1..k_S
//     CKPT_DIR: directory for checkpoint files
//
// Output lines are machine readable (key=value).
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

typedef unsigned __int128 u128;
typedef struct { uint64_t x; uint32_t K; uint32_t pad; double G, nu; } St;   // 32 bytes
typedef struct { uint64_t x; uint32_t K; } Ck;                                  // 16 bytes in memory

static int SIGN; static uint64_t H; static long LAYERS; static long KTARGET; static double GMAX, MARGIN;
static uint64_t XMAX; static double LOG2H;
static int ESF=0;   // env ES_FILTER=1: E_S (excursion only at v = 6 mod 8) -> M-points must be odd or (sign-1: u = 2 mod 8; sign+1: x = 6 mod 8)
// log2(3) as a double-double: L_HI has 26 significant bits so i*L_HI is exact for i < 2^27.
static const double L_FULL = 1.5849625007211562;
static double L_HI, L_LO;

static int cmp_st(const void*a,const void*b){ uint64_t x=((const St*)a)->x, y=((const St*)b)->x; return x<y?-1:(x>y?1:0); }

// radix sort by x (64-bit keys, but values < 2^40 typically); LSD on 8-bit digits
static void radix_sort(St*a, St*tmp, size_t n, uint64_t maxkey){
  int passes=0; uint64_t m=maxkey; while(m){ passes++; m>>=8; } if(passes==0) passes=1;
  St*src=a,*dst=tmp;
  for(int p=0;p<passes;p++){
    size_t cnt[257]; memset(cnt,0,sizeof cnt);
    int sh=8*p;
    for(size_t i=0;i<n;i++) cnt[((src[i].x>>sh)&255)+1]++;
    for(int b=0;b<256;b++) cnt[b+1]+=cnt[b];
    for(size_t i=0;i<n;i++) dst[cnt[(src[i].x>>sh)&255]++]=src[i];
    St*t=src;src=dst;dst=t;
  }
  if(src!=a) memcpy(a,src,n*sizeof(St));
}

static int better(const St*a,const St*b){ // is a better than b (same x)?
  if(a->K!=b->K) return SIGN>0 ? (a->K<b->K) : (a->K>b->K);
  return a->G < b->G;
}

// expand layer i (cur, n) into children (out), returns count; layer index i = number of moves done
static size_t expand(const St*cur,size_t n,long i,St*out,size_t cap,int*overflow){
  size_t m=0; (void)i;
  for(size_t j=0;j<n;j++){
    uint64_t x=cur[j].x; int r=(int)(x%3);
    int k0 = (SIGN>0) ? (r==1?0:1) : (r==2?0:1);
    for(int k=k0;k<=100;k+=2){
      u128 t=((u128)x)<<k;
      if(SIGN>0){ if(t>(u128)3*XMAX+1) break; if((uint64_t)(t%9)==1) continue; }
      else      { if(t>(u128)3*XMAX-1) break; if((uint64_t)(t%9)==8) continue; }
      uint64_t y = (SIGN>0) ? (uint64_t)((t-1)/3) : (uint64_t)((t+1)/3);
      if(y==0) continue;
      if(ESF && !(y&1) && (y&7)!=(uint64_t)(SIGN>0?6:2)) continue;
      double f = 1.0 - ldexp(1.0,-k);
      double G2 = cur[j].G + f*cur[j].nu;
      if(G2 > GMAX + MARGIN) continue;
      if(m>=cap){ *overflow=1; return m; }
      St s; s.x=y; s.K=cur[j].K+(uint32_t)k; s.pad=0; s.G=G2;
      s.nu = 3.0*ldexp(cur[j].nu,-k) + ldexp(1.0/(double)H, -(int)s.K);
      out[m++]=s;
    }
  }
  return m;
}

static size_t dedup(St*a,size_t n){
  if(n==0) return 0; size_t w=0;
  for(size_t i=1;i<n;i++){
    if(a[i].x==a[w].x){ if(better(&a[i],&a[w])) a[w]=a[i]; }
    else a[++w]=a[i];
  }
  return w+1;
}

static long find(const St*a,size_t n,uint64_t x){
  size_t lo=0,hi=n; while(lo<hi){ size_t mid=(lo+hi)/2; if(a[mid].x<x) lo=mid+1; else hi=mid; }
  if(lo<n && a[lo].x==x) return (long)lo; return -1;
}

// direct (non-incremental) G and nu of a state (i, x, K); used after loading checkpoints
static void direct_G(long i,uint64_t x,uint32_t K,double*G,double*nu){
  double iL_hi=(double)i*L_HI, iL_lo=(double)i*L_LO;
  double big = iL_hi - (double)(K+1);          // exact
  if(SIGN>0){
    double t = big + iL_lo + log2(2.0*(double)x+1.0) - LOG2H;
    *G = 1.0 - exp2(t) + ldexp(1.0/(double)H,-(int)(K+1));
  } else {
    double t = big + iL_lo + log2(2.0*(double)x-1.0) - LOG2H;
    *G = exp2(t) + ldexp(1.0/(double)H,-(int)(K+1)) - 1.0;
  }
  double t2 = big + iL_lo - LOG2H;
  *nu = exp2(t2) - ldexp(1.0/(double)H,-(int)(K+1));
  if(i==0){ *G=0.0; *nu=0.0; }
}

static uint64_t climb_point(int N){ // c_N
  u128 p=1; for(int j=0;j<N+(SIGN>0?1:0);j++) p*=3;
  if(SIGN>0) return (uint64_t)((p-1)/2); else return (uint64_t)((p+1)/2);
}

int main(int argc,char**argv){
  if(argc<7){ fprintf(stderr,"usage: tightdp SIGN H LAYERS KTARGET GMAX XMAX [CLIMB] [RECON_T] [CKPT_DIR]\n"); return 1; }
  SIGN=atoi(argv[1]); H=strtoull(argv[2],0,10); LAYERS=atol(argv[3]); KTARGET=atol(argv[4]);
  GMAX=strtod(argv[5],0); XMAX=strtoull(argv[6],0,10);
  int CLIMB = argc>7 ? atoi(argv[7]) : 0;
  long RT = argc>8 ? atol(argv[8]) : 0;
  const char* CDIR = argc>9 ? argv[9] : ".";
  MARGIN = 1e-9*GMAX + 1e-13;
  ESF = getenv("ES_FILTER") ? atoi(getenv("ES_FILTER")) : 0;
  LOG2H = log2((double)H);
  // log2(3) = L_HI + L_LO with L_HI = floor(2^25 log2 3)/2^25 (exact) and L_LO correctly rounded
  L_HI = 1.58496248722076416015625; L_LO = 1.350039202129749e-08; (void)L_FULL;
  size_t cap = 1<<20; St *cur=malloc(cap*sizeof(St)), *nxt=malloc(cap*sizeof(St)), *tmp=malloc(cap*sizeof(St));
  size_t n=1; cur[0].x=H; cur[0].K=0; cur[0].G=0; cur[0].nu=0; cur[0].pad=0;
  size_t maxn=1, totn=0; long found_layer=-1; int best_climb=-1; long climb_layer=-1;
  uint64_t maxx=H; int climbN_max = 0; while(climb_point(climbN_max+1) <= XMAX && climbN_max<38) climbN_max++;
  FILE*ck=NULL; char ckname[512];
  if(RT>0){ snprintf(ckname,sizeof ckname,"%s/ckpt_%d_%llu_%ld.bin",CDIR,SIGN,(unsigned long long)H,LAYERS); ck=fopen(ckname,"wb"); if(!ck){perror("ckpt");return 1;} }
  long *ckoff = RT>0 ? calloc(LAYERS/RT+2,sizeof(long)) : NULL; size_t *ckn = RT>0 ? calloc(LAYERS/RT+2,sizeof(size_t)) : NULL;
  for(long i=0;i<=LAYERS;i++){
    // checkpoint
    if(RT>0 && i%RT==0){
      long c=i/RT; ckoff[c]=ftell(ck); ckn[c]=n;
      for(size_t j=0;j<n;j++){ fwrite(&cur[j].x,8,1,ck); fwrite(&cur[j].K,4,1,ck); }
    }
    // climb check: layer i = LAYERS - N holds c_N ?
    if(CLIMB){
      long N = LAYERS - i;
      if(N>=1 && N<=climbN_max){
        uint64_t cN=climb_point((int)N); long idx=find(cur,n,cN);
        if(idx>=0 && (long)cur[idx].K==KTARGET){ if((int)N>best_climb){ best_climb=(int)N; climb_layer=i; } }
      }
    }
    if(i==LAYERS){
      long idx=find(cur,n,H);
      if(idx>=0 && (long)cur[idx].K==KTARGET) found_layer=i;
      break;
    }
    // expand
    int of=0; size_t m;
    for(;;){
      of=0; m=expand(cur,n,i,nxt,cap,&of);
      if(!of) break;
      cap*=2; nxt=realloc(nxt,cap*sizeof(St)); tmp=realloc(tmp,cap*sizeof(St)); cur=realloc(cur,cap*sizeof(St));
      if(!nxt||!tmp||!cur){ fprintf(stderr,"oom\n"); return 2; }
    }
    radix_sort(nxt,tmp,m,XMAX);
    m=dedup(nxt,m);
    St*t=cur;cur=nxt;nxt=t; n=m; totn+=n; if(n>maxn) maxn=n;
    if(n && cur[n-1].x>maxx) maxx=cur[n-1].x;
    if(n==0){ printf("status=dead layer=%ld\n",i+1); break; }
    if((i+1)%5000==0){ fprintf(stderr,"layer %ld states %zu\n",i+1,n); }
  }
  printf("sign=%d H=%llu layers=%ld K=%ld Gmax=%.6e xmax=%llu maxstates=%zu totstates=%zu maxx=%llu found=%d climb=%d\n",
    SIGN,(unsigned long long)H,LAYERS,KTARGET,GMAX,(unsigned long long)XMAX,maxn,totn,(unsigned long long)maxx,found_layer>=0,best_climb);
  fflush(stdout);
  if(RT>0){
    fclose(ck);
    long end_layer; uint64_t ex;
    if(CLIMB && best_climb>0){ end_layer=climb_layer; ex=climb_point(best_climb); }
    else if(found_layer>=0){ end_layer=LAYERS; ex=H; }
    else { remove(ckname); return 0; }
    // reconstruction
    int *ks=malloc((LAYERS+2)*sizeof(int)); uint64_t *xs=malloc((LAYERS+2)*sizeof(uint64_t));
    long pos=end_layer; uint64_t cx=ex; uint32_t cK=(uint32_t)KTARGET; xs[pos]=cx;
    FILE*in=fopen(ckname,"rb");
    Ck **segc=calloc(RT+1,sizeof(Ck*)); size_t *segn=calloc(RT+1,sizeof(size_t)); size_t *segcapc=calloc(RT+1,sizeof(size_t));
    size_t wcap=cap; St *wa=malloc(wcap*sizeof(St)), *wb=malloc(wcap*sizeof(St));
    for(long c=(end_layer-1)/RT; c>=0; c--){
      long s0=c*RT, s1=pos; // recompute layers s0..s1-1 (we are at layer s1 = pos)
      fseek(in,ckoff[c],SEEK_SET); size_t n0=ckn[c];
      if(wcap<n0){ wcap=2*n0; wa=realloc(wa,wcap*sizeof(St)); wb=realloc(wb,wcap*sizeof(St)); }
      if(segcapc[0]<n0){ segcapc[0]=n0; segc[0]=realloc(segc[0],n0*sizeof(Ck)); }
      for(size_t j=0;j<n0;j++){ uint64_t x; uint32_t K; if(fread(&x,8,1,in)!=1||fread(&K,4,1,in)!=1){fprintf(stderr,"read\n");return 3;}
        wa[j].x=x; wa[j].K=K; wa[j].pad=0; direct_G(s0,x,K,&wa[j].G,&wa[j].nu); segc[0][j].x=x; segc[0][j].K=K; }
      segn[0]=n0; size_t wn=n0;
      double saveM=MARGIN; MARGIN = 1e-7*GMAX + 1e-11;   // looser margin: superset of true states, all genuine
      for(long l=s0;l<s1-1;l++){
        long a=l-s0; int of; size_t m;
        for(;;){ of=0; m=expand(wa,wn,l,wb,wcap,&of); if(!of) break;
          wcap*=2; wa=realloc(wa,wcap*sizeof(St)); wb=realloc(wb,wcap*sizeof(St)); tmp=realloc(tmp,wcap*sizeof(St)); cap=wcap; }
        if(cap<wcap){ cap=wcap; tmp=realloc(tmp,cap*sizeof(St)); }
        radix_sort(wb,tmp,m,XMAX); m=dedup(wb,m);
        if(segcapc[a+1]<m){ segcapc[a+1]=m+16; segc[a+1]=realloc(segc[a+1],segcapc[a+1]*sizeof(Ck)); }
        for(size_t j=0;j<m;j++){ segc[a+1][j].x=wb[j].x; segc[a+1][j].K=wb[j].K; }
        segn[a+1]=m;
        St*t=wa; wa=wb; wb=t; wn=m;
      }
      MARGIN=saveM;
      // backtrack from layer s1 to s0
      for(long l=s1; l>s0; l--){
        long a=l-1-s0; int ok=0;
        u128 num = (SIGN>0) ? ((u128)3*cx+1) : ((u128)3*cx-1);
        for(int k=0;k<=100 && !ok;k++){
          if(k>0 && (num & (((u128)1<<k)-1))) break;
          u128 xp128 = num>>k; if(xp128==0) break; if(xp128>(u128)XMAX) continue;
          uint64_t xp=(uint64_t)xp128;
          if((long)cK - k < 0) break;
          // binary search in compact layer
          size_t lo=0,hi=segn[a]; while(lo<hi){ size_t mid=(lo+hi)/2; if(segc[a][mid].x<xp) lo=mid+1; else hi=mid; }
          if(lo<segn[a] && segc[a][lo].x==xp && segc[a][lo].K==cK-(uint32_t)k){ ks[l]=k; cx=xp; cK-=(uint32_t)k; xs[l-1]=cx; ok=1; }
        }
        if(!ok){ fprintf(stderr,"backtrack failed at layer %ld x=%llu K=%u\n",l,(unsigned long long)cx,cK); return 5; }
      }
      pos=s0;
    }
    fclose(in); remove(ckname);
    if(pos!=0 || cx!=H || cK!=0){ fprintf(stderr,"reconstruction did not reach the start\n"); return 6; }
    // complete a climb-point object with the k=0 descent
    printf("recon_end_layer=%ld\n",end_layer);
    printf("moves=");
    for(long l=1;l<=end_layer;l++) printf("%d%s",ks[l],l<end_layer?",":"");
    for(long l=end_layer+1;l<=LAYERS;l++) printf(",0");
    printf("\n");
    uint64_t hmax=0; for(long l=0;l<=end_layer;l++) if(xs[l]>hmax) hmax=xs[l];
    printf("height=%llu\n",(unsigned long long)hmax);
  }
  return 0;
}
