/* mirror_check.c - opus-2026-11-09-S87.
 * The blue-line H-mod-4 law at odd order n = 2c+1 via the mirror-path decomposition (see mirror_paths.py):
 *   (H(t)-1)/2 = A(y) + sum_i x_i B_i(y) (mod 2), B_i(y) = N(i) + N(n-1-i), N(u) = # mirror-free paths c -> u.
 * For every y (the blue coordinates on the transposed tile pairs; the mirror tiles x do not affect N) this checks
 *   (P2) sum_i B_i(y) = 1 (mod 2)   and   (P3) B(y) != B(y + 1),
 * where y + 1 is the complement tiling's y-part.  It also checks the decomposition against H (Hamiltonian-path DP) on
 * SAMPLES random blue tilings.
 * usage: mirror_check n SLICE NSLICES SAMPLES
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned long long u64;
#define MAXN 15
static int n, c, m, ti[110], tj[110];
static int Nc[MAXN];
static uint16_t ADJ[MAXN];
static void dfs(int v, unsigned used){
  unsigned nb = ADJ[v] & ~used;
  while (nb){ int u = __builtin_ctz(nb); nb &= nb - 1;
    if (u == c) continue;
    if (used >> (n - 1 - u) & 1) continue;          /* mirror-free: never both u and iota(u) */
    Nc[u] ^= 1;                                     /* parities suffice */
    dfs(u, used | (1u << u)); }
}
static void build(u64 lo, u64 hi, uint16_t *adj){       /* tiling bits in two words (m <= 128) */
  memset(adj, 0, sizeof(uint16_t)*MAXN);
  for (int i=0;i+1<n;i++) adj[i] |= 1u<<(i+1);
  for (int b=0;b<m;b++){ int bit = b < 64 ? (lo>>b & 1) : (hi>>(b-64) & 1); if (bit) adj[tj[b]] |= 1u<<ti[b]; else adj[ti[b]] |= 1u<<tj[b]; }
}
static void Bvec(u64 lo, u64 hi, int *B){
  build(lo, hi, ADJ); memset(Nc, 0, sizeof Nc); dfs(c, 1u << c);
  for (int i=0;i<c;i++) B[i] = Nc[i] ^ Nc[n-1-i];
}
static u64 hpaths(const uint16_t *adj){ static u64 f[1<<13][13]; int F=1<<n;
  for (int S=0;S<F;S++) for (int v=0;v<n;v++) f[S][v]=0; for (int v=0;v<n;v++) f[1<<v][v]=1;
  for (int S=1;S<F;S++) for (int v=0;v<n;v++) if (f[S][v]){ u64 x=f[S][v]; uint16_t o=adj[v]&~S; while(o){ int u=__builtin_ctz(o); o&=o-1; f[S|1<<u][u]+=x; } }
  u64 h=0; for (int v=0;v<n;v++) h+=f[F-1][v]; return h; }
int main(int argc, char **argv){
  n = atoi(argv[1]); c = n/2; int slice = atoi(argv[2]), ns = atoi(argv[3]); long samples = atol(argv[4]);
  int idx[MAXN][MAXN]; m=0; for (int i=0;i<n;i++) for (int j=i+2;j<n;j++){ ti[m]=i; tj[m]=j; idx[i][j]=m; m++; }
  int sg[110]; for (int b=0;b<m;b++) sg[b] = idx[n-1-tj[b]][n-1-ti[b]];
  int pa[64], pb[64], np=0, mir[MAXN]; int nm=0;
  for (int b=0;b<m;b++){ if (sg[b]==b) mir[nm++]=b; else if (b < sg[b]){ pa[np]=b; pb[np]=sg[b]; np++; } }
  /* all-ones on the non-mirror tiles */
  long total = 1L << np, bad2 = 0, bad3 = 0, checked = 0; static long joint[2][2];
  for (long y = slice; y < total; y += ns){
    u64 lo=0, hi=0, lo2=0, hi2=0;
    for (int q=0;q<np;q++){ int on = y>>q & 1; int bs[2] = {pa[q], pb[q]};
      for (int z=0;z<2;z++){ int b = bs[z];
        if (on){ if (b<64) lo |= 1ULL<<b; else hi |= 1ULL<<(b-64); }
        else   { if (b<64) lo2 |= 1ULL<<b; else hi2 |= 1ULL<<(b-64); } } }
    int B[MAXN], Bp[MAXN]; Bvec(lo, hi, B); Bvec(lo2, hi2, Bp);
    int s = 0, eq = 1; for (int i=0;i<c;i++){ s ^= B[i]; if (B[i] != Bp[i]) eq = 0; }
    if (s != 1) bad2++; if (eq) bad3++; checked++;
    /* exact joint distribution over the 2^c choices of the mirror tiles x:
       f(t) = A + x.B and f(t+1) = A' + (x+1).B' */
    int A = 0, Ap = 0; { int Bt[MAXN]; Bvec(lo, hi, Bt); for (int i=0;i<c;i++) A ^= Nc[i]; Bvec(lo2, hi2, Bt); for (int i=0;i<c;i++) Ap ^= Nc[i]; }
    int sBp = 0; for (int i=0;i<c;i++) sBp ^= Bp[i];
    for (int x=0; x<(1<<c); x++){ int f1 = A, f2 = Ap ^ sBp; for (int i=0;i<c;i++) if (x>>i & 1){ f1 ^= B[i]; f2 ^= Bp[i]; } joint[f1][f2]++; }
  }
  /* sampled check of the decomposition against H */
  long dbad = 0; srand(12345 + slice);
  for (long k=0; k<samples; k++){
    long y = (((long)rand()<<15) ^ rand()) % total; int x = rand() % (1<<c);
    u64 lo=0, hi=0;
    for (int q=0;q<np;q++) if (y>>q & 1){ int bs[2] = {pa[q], pb[q]}; for (int z=0;z<2;z++){ int b=bs[z]; if (b<64) lo |= 1ULL<<b; else hi |= 1ULL<<(b-64); } }
    int B[MAXN]; Bvec(lo, hi, B); int A = 0; for (int i=0;i<c;i++) A ^= Nc[i];
    int f = A; for (int i=0;i<c;i++) if (x>>i & 1){ f ^= B[i]; int b = mir[0]; for (int z=0; z<nm; z++) if (ti[mir[z]]==i) b = mir[z]; if (b<64) lo |= 1ULL<<b; else hi |= 1ULL<<(b-64); }
    uint16_t adj[MAXN]; build(lo, hi, adj); u64 h = hpaths(adj);
    if ((int)(((h-1)/2) & 1) != f) dbad++;
  }
  printf("n=%d slice %d/%d: y checked %ld (of %ld); (P2) failures %ld; (P3) failures %ld; decomposition sampled %ld, mismatches %ld; joint tilings (f,f') 00 %ld 01 %ld 10 %ld 11 %ld\n",
         n, slice, ns, checked, total, bad2, bad3, samples, dbad, joint[0][0], joint[0][1], joint[1][0], joint[1][1]);
  return 0;
}
