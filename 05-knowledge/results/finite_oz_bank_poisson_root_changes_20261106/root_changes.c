/* root_changes.c - kind-pasteur-2026-11-06-S13 (HYP-9480 (D), prediction P4).
 * H-values (Hamiltonian-path counts) of the physical sources t . w_1 (J w_2) ... (J w_(j+1)), j <= K root changes,
 * from the rooted edge t = (D,E,U,d) = (1,0,2,1), all words w_i over {O,Z}.
 * Physical moves: O(D,E,U,d) = (D,0,2U,d); Z = (D-E+U,d,2U-d,d); J = (D-E+U,D-E,D+U,D).
 * The final O/Z block is enumerated through Session 205's closed form h + R m - d Q(m) - e (m mod 2) by a DFS over m
 * (Q updated exactly, pruned by H >= h + (R-d) m); intermediate blocks letter by letter (their trailing O letters
 * matter before a J), pruned when U > X (a later J then gives D >= U > X).
 * usage: root_changes X K   (output deterministic text; P4 used X = 100000000, K = 2)
 * The code below is byte-identical in logic to the exploratory kj.c that produced root_changes_1e8.txt.
 */
// H-values of sources t.w1 (J w)^j for j <= K root changes, from the rooted edge t=(D,E,U,d)=(1,0,2,1).
// Exact physical moves: O(D,E,U,d)=(D,0,2U,d); Z=(D-E+U,d,2U-d,d); J=(D-E+U,D-E,D+U,D).
// The final O/Z block is enumerated through the closed form h+Rm-dQ(m)-e(m mod 2) (S205) by a DFS over m (no arrays).
// Intermediate blocks are enumerated letter by letter (trailing O letters matter before a J); prune when U > X (then J gives D >= U > X).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll; typedef unsigned char u8;
static ll X; static int K; static u8 *hit[4];
static inline void mark(u8*b,ll v){ if(v>=1&&v<=X&&(v&1)) b[v>>4]|=(u8)(1u<<((v>>1)&7)); }
static inline int get(u8*b,ll v){ return (b[v>>4]>>((v>>1)&7))&1; }
static ll H0,R0,E0,D0; static u8* tgt;
static void tdfs(ll m, ll q){ // m>=1
  if (H0+(R0-D0)*m > X) return;
  mark(tgt, H0+R0*m-D0*q-E0*(m&1));
  tdfs(2*m,q); tdfs(2*m+1,q+m+(m&1));
}
static void tail(u8*b, ll h, ll R, ll e, ll d){ if(h>X) return; H0=h;R0=R;E0=e;D0=d;tgt=b; mark(b,h); tdfs(1,0); }
static void block(ll D, ll E, ll U, ll d, int j){ // node inside block j+1 (after j root changes); all its O/Z extensions' values are in tail()
  // root change here (if allowed): new state J
  if (j < K){ ll D2=D-E+U, E2=D-E, U2=D+U, d2=D; if (D2<=X){ tail(hit[j+1], D2, U2, E2, d2); block(D2,E2,U2,d2,j+1); } }
  // extend the current block by one letter, only needed to reach later root changes (U<=X keeps J values <= X possible)
  if (j < K){ if (2*U <= X) block(D,0,2*U,d,j); if (2*U-d <= X) block(D-E+U,d,2*U-d,d,j); }
}
int main(int argc,char**argv){
  X=atoll(argv[1]); K=atoi(argv[2]);
  for(int i=0;i<=K;i++) hit[i]=calloc((X>>4)+2,1);
  tail(hit[0],1,2,0,1);          // zero root changes: tails of t
  block(1,0,2,1,0);
  for(int i=1;i<=K;i++) for(ll b=0;b<(X>>4)+2;b++) hit[i][b]|=hit[i-1][b];
  for(int i=0;i<=K;i++){ ll miss=0,last=0; for(ll v=1;v<=X;v+=2) if(!get(hit[i],v)){miss++;last=v;}
    printf("<=%d root changes: odd misses %lld, largest %lld\n",i,miss,last); }
  printf("misses (<=%d) above 10^5:",K); for(ll v=100001;v<=X;v+=2) if(!get(hit[K],v)) printf(" %lld",v); printf("\n");
  for(ll lo=1;lo<X;lo*=2){ ll hi=lo*2>X?X:lo*2; printf("  [%lld,%lld):",lo,hi); for(int i=0;i<=K;i++){ ll c=0; for(ll v=lo|1; v<hi; v+=2) if(!get(hit[i],v)) c++; printf(" k%d:%lld",i,c);} printf("\n"); }
  return 0;
}
