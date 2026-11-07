/* rc_depth.c (variant of rc_census.c) - kind-pasteur-2026-11-07-S14 (HYP-9480 (D) extension).
 * Same physical grammar and enumeration order as the S13 packet's root_changes.c:
 * sources t . w_1 (J w_2) ... (J w_(j+1)), j <= K root changes, from the rooted edge t = (D,E,U,d) = (1,0,2,1);
 * O(D,E,U,d) = (D,0,2U,d), Z = (D-E+U,d,2U-d,d), J = (D-E+U,D-E,D+U,D).
 * The final O/Z block is enumerated through Session 205's closed form h + R m - d Q(m) - e (m mod 2) by a DFS over m;
 * intermediate blocks letter by letter, pruned when U > X.
 * rc_depth.c: identical enumeration; up to 7 root-change levels and a 4th argument LOPRINT (print every miss of the
 * top level that is >= LOPRINT). Used for the root-change depth table (exploratory, FINITE-EXACT).
 * New here: OpenMP task split of the block tree, per-level coverage bitmaps (as before) AND a saturating 8-bit
 * multiplicity counter for the union (<= K root changes): the number of grammar words (final block modulo trailing O)
 * with value N.  usage: rc_depth X K Usplit LOPRINT
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>
typedef long long ll; typedef unsigned char u8;
static ll X; static int K; static u8 *hit[8]; static u8 *cnt;
static inline void mark(u8*b,ll v){
  if(v>=1&&v<=X&&(v&1)){
    __atomic_fetch_or(&b[v>>4],(u8)(1u<<((v>>1)&7)),__ATOMIC_RELAXED);
    u8 *c=&cnt[v>>1]; u8 o=__atomic_load_n(c,__ATOMIC_RELAXED);
    while(o<255 && !__atomic_compare_exchange_n(c,&o,(u8)(o+1),1,__ATOMIC_RELAXED,__ATOMIC_RELAXED)){}
  }
}
static inline int get(u8*b,ll v){ return (b[v>>4]>>((v>>1)&7))&1; }
typedef struct { ll h,R,e,d; u8 *t; } tailctx;
static void tdfs(const tailctx*c, ll m, ll q){ // m>=1
  if (c->h+(c->R-c->d)*m > X) return;
  mark(c->t, c->h+c->R*m-c->d*q-c->e*(m&1));
  tdfs(c,2*m,q); tdfs(c,2*m+1,q+m+(m&1));
}
static void tail(u8*b, ll h, ll R, ll e, ll d){ if(h>X) return; tailctx c={h,R,e,d,b}; mark(b,h); tdfs(&c,1,0); }
static void block(ll D, ll E, ll U, ll d, int j){
  if (j < K){ ll D2=D-E+U, E2=D-E, U2=D+U, d2=D; if (D2<=X){ tail(hit[j+1], D2, U2, E2, d2); block(D2,E2,U2,d2,j+1); } }
  if (j < K){ if (2*U <= X) block(D,0,2*U,d,j); if (2*U-d <= X) block(D-E+U,d,2*U-d,d,j); }
}
/* task list: type 0 = tail(level j, seed), type 1 = block subtree from node (D,E,U,d,j) */
typedef struct { int type,j; ll D,E,U,d; } task;
static task *T; static size_t nT=0, capT=0; static ll USPLIT;
static void push(int type,int j,ll D,ll E,ll U,ll d){ if(nT==capT){capT=capT?2*capT:1024; T=realloc(T,capT*sizeof(task));} T[nT++] = (task){type,j,D,E,U,d}; }
static void expand(ll D, ll E, ll U, ll d, int j){ // same semantics as block(), but small-U nodes are split into tasks
  if (U >= USPLIT){ push(1,j,D,E,U,d); return; }
  if (j < K){ ll D2=D-E+U, E2=D-E, U2=D+U, d2=D; if (D2<=X){ push(0,j+1,D2,E2,U2,d2); expand(D2,E2,U2,d2,j+1); } }
  if (j < K){ if (2*U <= X) expand(D,0,2*U,d,j); if (2*U-d <= X) expand(D-E+U,d,2*U-d,d,j); }
}
int main(int argc,char**argv){
  if(argc<3){fprintf(stderr,"usage: rc_census X K [Usplit]\n");return 1;}
  X=atoll(argv[1]); K=atoi(argv[2]); USPLIT = argc>3? atoll(argv[3]) : 4096; ll LOPRINT = argc>4? atoll(argv[4]) : 100001;
  if(K>7){fprintf(stderr,"K must be at most 7\n");return 1;}
  for(int i=0;i<=K;i++) hit[i]=calloc((X>>4)+2,1);
  cnt=calloc((X>>1)+2,1);
  push(0,0,1,0,2,1);                 // zero root changes: tails of t  (seed h=D=1,R=U=2,e=E=0,d=1)
  expand(1,0,2,1,0);
  fprintf(stderr,"tasks %zu\n",nT);
  double t0=omp_get_wtime();
  #pragma omp parallel for schedule(dynamic,1)
  for (ll i=0;i<(ll)nT;i++){ task*k=&T[i];
    if(k->type==0) tail(hit[k->j], k->D, k->U, k->E, k->d); else block(k->D,k->E,k->U,k->d,k->j); }
  fprintf(stderr,"enumeration %.1fs\n",omp_get_wtime()-t0);
  for(int i=1;i<=K;i++) for(ll b=0;b<(X>>4)+2;b++) hit[i][b]|=hit[i-1][b];
  for(int i=0;i<=K;i++){ ll miss=0,last=0; for(ll v=1;v<=X;v+=2) if(!get(hit[i],v)){miss++;last=v;}
    printf("<=%d root changes: odd misses %lld, largest %lld\n",i,miss,last); }
  printf("misses (<=%d) from %lld:",K,LOPRINT); for(ll v=LOPRINT|1;v<=X;v+=2) if(!get(hit[K],v)) printf(" %lld",v); printf("\n");
  for(ll lo=1;lo<X;lo*=2){ ll hi=lo*2>X?X:lo*2; printf("  [%lld,%lld):",lo,hi); for(int i=0;i<=K;i++){ ll c=0; for(ll v=lo|1; v<hi; v+=2) if(!get(hit[i],v)) c++; printf(" k%d:%lld",i,c);} printf("\n"); }
  /* union multiplicity: per dyadic window, min multiplicity, counts of multiplicity 1..8, mean, and per residue mod 3 the min */
  printf("multiplicity (<=%d root changes, words with saturating counter 255):\n",K);
  for(ll lo=1;lo<X;lo*=2){ ll hi=lo*2>X?X:lo*2; ll c[9]={0}; double s=0; ll n=0; int mn=999, mn3[3]={999,999,999}; ll sat=0;
    for(ll v=lo|1; v<hi; v+=2){ int m=cnt[v>>1]; n++; s+=m; if(m<mn)mn=m; if(m<mn3[v%3])mn3[v%3]=m; if(m<=8)c[m]++; if(m==255)sat++; }
    printf("  [%lld,%lld): n=%lld mean=%.2f min=%d min_mod3=(%d,%d,%d) sat255=%lld c0..8=",lo,hi,n,s/n,mn,mn3[0],mn3[1],mn3[2],sat);
    for(int m=0;m<=8;m++) printf("%s%lld",m?",":"",c[m]); printf("\n"); }
  printf("odd N > 10^5 with multiplicity <= 2:"); for(ll v=100001;v<=X;v+=2) if(cnt[v>>1]<=2) printf(" %lld:%d",v,cnt[v>>1]); printf("\n");
  return 0;
}
