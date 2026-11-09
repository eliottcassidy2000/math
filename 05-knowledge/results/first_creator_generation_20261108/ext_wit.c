/* ext_wit.c (witness variant of ext.c) - kind-pasteur-2026-11-08-S17.  One-vertex extensions of every tournament class of order n, through the
 * quadratic attachment form.
 *
 * For a tournament T on V = {0..n-1} and a new vertex v with x_a = [a -> v]:
 *   H(T + v_x) = H(T) + sum_a x_a L_a - sum_{a<b} x_a x_b M_ab,
 *   W[a][b] = sum over S with a in S, b not in S of f[S][a] g[V\S][b]
 *             (f[S][a] = Hamiltonian paths of T[S] ending at a, g[R][b] = those of T[R] starting at b),
 *   L_a = F_a - G_a + sum_{b != a} W[a][b]   (F_a = f[V][a], G_a = g[V][a]),   M_ab = W[a][b] + W[b][a].
 * (v first: sum_b (1-x_b) G_b; v last: sum_a x_a F_a; v after a prefix on S ending at a, before a suffix starting at b:
 *  x_a (1 - x_b) f[S][a] g[V\S][b].)  All 2^n patterns are visited in Gray order with O(n) updates per step.
 *
 * For first creators T (H(T) not in S(n-1)) it also carries the forms of the n deletions T - u, which give
 * H(T' - u) for T' = T + v_x, and counts the hypo analogues: T' with H(T') in S(n) (not a first creator at n+1)
 * whose every deletion is a first creator at order n (H(T' - v) = H(T) and H(T' - u) not in S(n-1) for all u).
 *
 * usage: ext <n> <fmt: txt|bin> <classfile> <lo> <hi> <S(n-1) values> <S(n) values> <outprefix> [nchk]
 *   txt: lines "code H ..." (classgen.c cls_n.txt);  bin: 34-byte records of dump10.c.
 *   S files: one value per line.  Writes <outprefix>.cnt (u64 counts per H of all (class, pattern) pairs, then of
 *   first-creator pairs) and prints statistics.  nchk > 0: brute-force check of every pattern for the first nchk
 *   classes of the slice (direct subset DP at order n+1, including all deletion counts). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned long long u64;
#define MAXN 12
#define HCAP 100000
static int n;
static unsigned char inPrev[HCAP], inCur[HCAP];        /* S(n-1), S(n) */
static u64 cntAll[HCAP], cntFC[HCAP]; static u64 witCode[HCAP]; static unsigned witX[HCAP]; static unsigned char haveW[HCAP];
static uint32_t f[1<<11][MAXN], g[1<<11][MAXN];

static void decode(u64 c, int m, uint16_t *adj){ memset(adj,0,sizeof(uint16_t)*m); int b=0;
  for (int i=0;i<m;i++) for (int j=i+1;j<m;j++){ if (c>>b & 1) adj[i] |= 1u<<j; else adj[j] |= 1u<<i; b++; } }

static void dps(const uint16_t *adj){
  int F = 1<<n;
  for (int S=0;S<F;S++) for (int v=0;v<n;v++){ f[S][v]=0; g[S][v]=0; }
  for (int v=0;v<n;v++){ f[1<<v][v]=1; g[1<<v][v]=1; }
  uint16_t in[MAXN]; for (int v=0;v<n;v++){ in[v]=0; for (int u=0;u<n;u++) if (adj[u]>>v&1) in[v]|=1u<<u; }
  for (int S=1;S<F;S++) for (int v=0;v<n;v++){
    if (f[S][v]){ uint32_t x=f[S][v]; uint16_t o = adj[v] & ~S; while(o){ int u=__builtin_ctz(o); o&=o-1; f[S|1<<u][u]+=x; } }
    if (g[S][v]){ uint32_t x=g[S][v]; uint16_t o = in[v] & ~S; while(o){ int u=__builtin_ctz(o); o&=o-1; g[S|1<<u][u]+=x; } }
  }
}
/* form of T restricted to the vertex set U (bitmask): h0 = H(T[U]), L[a], M[a][b] for a,b in U */
static int32_t form(int U, int32_t *L, int32_t M[MAXN][MAXN]){
  int32_t W[MAXN][MAXN]; memset(W,0,sizeof W);
  for (int S=(U-1)&U; S>0; S=(S-1)&U){            /* nonempty proper subsets of U */
    int R = U & ~S;
    for (int a=0;a<n;a++) if ((S>>a&1) && f[S][a]){ int32_t fa=f[S][a];
      for (int b=0;b<n;b++) if ((R>>b&1) && g[R][b]) W[a][b] += fa * (int32_t)g[R][b]; }
  }
  int32_t h0=0; for (int b=0;b<n;b++) if (U>>b&1) h0 += g[U][b];
  for (int a=0;a<n;a++) if (U>>a&1){ int32_t s=(int32_t)f[U][a]-(int32_t)g[U][a];
    for (int b=0;b<n;b++) if ((U>>b&1) && b!=a){ s += W[a][b]; M[a][b] = W[a][b]+W[b][a]; }
    L[a]=s; }
  return h0;
}
static u64 hdirect(const uint16_t *adj, int m, u64 *dh){
  static u64 F2[1<<12][MAXN]; int F = 1<<m;
  for (int S=0;S<F;S++) for (int v=0;v<m;v++) F2[S][v]=0;
  for (int v=0;v<m;v++) F2[1<<v][v]=1;
  for (int S=1;S<F;S++) for (int v=0;v<m;v++) if (F2[S][v]){ u64 x=F2[S][v]; uint16_t o = adj[v] & ~S; while(o){ int u=__builtin_ctz(o); o&=o-1; F2[S|1<<u][u]+=x; } }
  u64 h=0; for (int v=0;v<m;v++) h+=F2[F-1][v];
  for (int x=0;x<m;x++){ int S=(F-1)&~(1<<x); u64 s=0; for(int v=0;v<m;v++) s+=F2[S][v]; dh[x]=s; }
  return h;
}
static void loadset(const char *fn, unsigned char *a){ FILE *fp=fopen(fn,"r"); long v; if(!fp){perror(fn); exit(1);} while (fscanf(fp,"%ld",&v)==1) if (v>=0 && v<HCAP) a[v]=1; fclose(fp); }

int main(int argc, char **argv){
  if (argc < 9){ fprintf(stderr,"usage\n"); return 1; }
  n = atoi(argv[1]); int bin = !strcmp(argv[2],"bin"); const char *cf = argv[3];
  long lo = atol(argv[4]), hi = atol(argv[5]); loadset(argv[6], inPrev); loadset(argv[7], inCur);
  const char *outp = argv[8]; long nchk = argc > 9 ? atol(argv[9]) : 0;
  FILE *fp = fopen(cf, bin ? "rb" : "r"); if (!fp){ perror(cf); return 1; }
  long idx = -1, done = 0, nfc = 0, hypo = 0, sterile = 0, chkfail = 0, chkdone = 0;
  long fertile_hist[12] = {0};   /* fertile patterns of first creators, binned by log2 */
  long lowfert = 0, minfert = 1L<<30; u64 minfert_code = 0; uint32_t minfert_H = 0;
  char line[8192];
  int full = (1<<n)-1;
  for (;;){
    u64 code; uint32_t H;
    if (bin){ unsigned char rec[34]; if (fread(rec,34,1,fp)!=1) break; memcpy(&code,rec,8); memcpy(&H,rec+8,4); }
    else { if (!fgets(line,sizeof line,fp)) break; unsigned long long c,h; if (sscanf(line,"%llu %llu",&c,&h)!=2) continue; code=c; H=(uint32_t)h; }
    idx++; if (idx < lo) continue; if (idx >= hi) break;
    uint16_t adj[MAXN]; decode(code, n, adj); dps(adj);
    int32_t L[MAXN], M[MAXN][MAXN]; int32_t h0 = form(full, L, M);
    if ((uint32_t)h0 != H){ fprintf(stderr,"H mismatch at %ld: %d vs %u\n", idx, h0, H); return 2; }
    int fc = !inPrev[H]; if (fc) nfc++;
    /* deletion forms for first creators */
    int32_t Lu[MAXN][MAXN], Mu[MAXN][MAXN][MAXN], hu[MAXN], Du[MAXN][MAXN];
    if (!fc) continue;
    int32_t D[MAXN]; for (int a=0;a<n;a++) D[a]=L[a];
    int32_t h = h0; unsigned x = 0; long fertile = 0;
    int chk = (chkdone < nchk);
    for (unsigned i=0;; i++){
      /* visit pattern x with value h */
      if (h < 0 || h >= HCAP){ fprintf(stderr,"range %d\n", h); return 3; }
      cntAll[h]++;
      if (fc && !inCur[h] && !haveW[h]){ haveW[h]=1; witCode[h]=code; witX[h]=x; }
      if (0){ cntFC[h]++; if (!inCur[h]) fertile++;
        else { int ok = 1; for (int u=0;u<n;u++) if (inPrev[hu[u]]){ ok = 0; break; } if (ok) { hypo++; if (hypo <= 20) printf("HYPO code=%llu x=%u H=%d\n", code, x, h); } } }
      if (chk){
        uint16_t a2[MAXN]; for (int v=0;v<n;v++) a2[v]=adj[v]; a2[n]=0;
        for (int v=0;v<n;v++){ if (x>>v&1) a2[v] |= 1u<<n; else a2[n] |= 1u<<v; }
        u64 dh[MAXN]; u64 hd = hdirect(a2, n+1, dh);
        if (hd != (u64)h) chkfail++;
        if (dh[n] != (u64)h0) chkfail++;
        if (fc) for (int u=0;u<n;u++) if (dh[u] != (u64)hu[u]) chkfail++;
      }
      if (i == (unsigned)full) break;
      int c = __builtin_ctz(i+1);
      if (!(x>>c & 1)){ h += D[c]; x |= 1u<<c; for (int b=0;b<n;b++) if (b!=c) D[b] -= M[b][c];
        }
      else { x &= ~(1u<<c); for (int b=0;b<n;b++) if (b!=c) D[b] += M[b][c]; h -= D[c];
        }
    }
    if (chk) chkdone++;
    if (fc){ if (fertile < (1L<<(n-1))){ lowfert++; if (lowfert <= 50) printf("LOWFERT code=%llu H=%u fertile=%ld\n", code, H, fertile); }
      if (fertile < minfert){ minfert = fertile; minfert_code = code; minfert_H = H; }
      if (fertile == 0) sterile++; int bucket = 0; long t = fertile; while (t){ bucket++; t >>= 1; } if (bucket > 11) bucket = 11; fertile_hist[bucket]++; }
    done++;
  }
  char fn[1024]; snprintf(fn, sizeof fn, "%s.wit", outp); FILE *o = fopen(fn,"w");
  for (int hh=0; hh<HCAP; hh++) if (haveW[hh]) fprintf(o, "%d %llu %u\n", hh, witCode[hh], witX[hh]); fclose(o);
  printf("slice [%ld,%ld): classes %ld, first creators %ld, sterile first creators %ld, hypo pairs %ld, checked classes %ld, check failures %ld\n", lo, hi, done, nfc, sterile, hypo, chkdone, chkfail);
  printf("fertile-pattern histogram (bucket b: 2^(b-1) <= #fertile < 2^b):"); for (int b=0;b<12;b++) printf(" %ld", fertile_hist[b]); printf("\n");
  printf("first creators with fewer than half of the patterns fertile: %ld; minimum fertile %ld (code %llu, H %u)\n", lowfert, minfert, minfert_code, minfert_H);
  return 0;
}
