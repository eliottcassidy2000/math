/* census_audit.c -- independent auditor census of root-change words (audit of kind-pasteur-2026-11-07-S14).
 * Written from the definitions only (no code shared with rc_census.c / rc_depth.c): every word is built letter by
 * letter; no closed form for the final O/Z block is used.
 *
 * Moves on s = (D,E,U,d) from t = (1,0,2,1): O = (D,0,2U,d), Z = (D-E+U,d,2U-d,d), J = (D-E+U,D-E,D+U,D).
 * A prefix w (any last letter) with j = #J(w) letters J contributes its value F(t.w) = D-E+U to
 *   level j coverage (F_<=k = {1} u {F(t.w) : #J(w) <= k}), and
 *   the multiplicity M_K(N) = [N=1] + #{wJ : #J(w) <= K-1, F = N} + #{wZ : #J(w) <= K, F = N}
 *   (+2 if j <= K-1, +1 if j = K); THM-4904 (1) is checked separately by chamber_audit.c.
 * F strictly increases along every letter (checked in chamber_audit.c on 7.6e7 states, and proved), so pruning
 * the DFS at F > HI is exact.
 *
 * The odd N in [1, X] are processed in segments of 2^SEGLOG integers so that the counters stay cache resident;
 * every segment re-traverses the word tree pruned at its upper end and counts only values inside it.
 * Counters saturate at 2^CNTBITS - 1 (CNTBITS = 4, 8, 16 or 32); LEVELS = 1 adds one coverage bitmap per exact j.
 *
 * usage: census_audit X K SEGLOG CNTBITS LEVELS [LISTLO] [USPLIT] [SEG0]  (SEG0: first segment start, a multiple of 2^SEGLOG)
 * output: per-level misses (if LEVELS), union misses (M = 0), window table per dyadic window, low lists.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>

typedef long long ll;
static ll X, LISTLO, USPLIT; static int K, SEGLOG, CNTBITS, LEVELS; static ll CAP;
static ll segA, segHI;                 /* current segment: count values F with segA <= F <= segHI */
static uint8_t *c8; static uint16_t *c16; static uint32_t *c32; static uint8_t *c4;
static uint64_t *lev[8];

static inline void bump(ll i, int inc){
  if (CNTBITS == 4){
    uint8_t *p = c4 + (i >> 1); int sh = (int)(i & 1) * 4;
    uint8_t o = __atomic_load_n(p, __ATOMIC_RELAXED);
    for (;;){ int v = (o >> sh) & 15; if (v == 15) return; int nv = v + inc; if (nv > 15) nv = 15;
      uint8_t n = (uint8_t)((o & ~(15 << sh)) | (nv << sh));
      if (__atomic_compare_exchange_n(p, &o, n, 1, __ATOMIC_RELAXED, __ATOMIC_RELAXED)) return; }
  } else if (CNTBITS == 8){
    uint8_t *p = c8 + i; uint8_t o = __atomic_load_n(p, __ATOMIC_RELAXED);
    for (;;){ if (o == 255) return; int nv = o + inc; if (nv > 255) nv = 255;
      if (__atomic_compare_exchange_n(p, &o, (uint8_t)nv, 1, __ATOMIC_RELAXED, __ATOMIC_RELAXED)) return; }
  } else if (CNTBITS == 16){
    uint16_t *p = c16 + i; uint16_t o = __atomic_load_n(p, __ATOMIC_RELAXED);
    for (;;){ if (o == 65535) return; int nv = o + inc; if (nv > 65535) nv = 65535;
      if (__atomic_compare_exchange_n(p, &o, (uint16_t)nv, 1, __ATOMIC_RELAXED, __ATOMIC_RELAXED)) return; }
  } else {
    __atomic_fetch_add(c32 + i, (uint32_t)inc, __ATOMIC_RELAXED);
  }
}
static inline ll getc_(ll i){
  if (CNTBITS == 4) return (c4[i >> 1] >> ((i & 1) * 4)) & 15;
  if (CNTBITS == 8) return c8[i];
  if (CNTBITS == 16) return c16[i];
  return c32[i];
}
static inline void record(ll F, int j){
  if (F < segA) return;
  ll i = (F - segA) >> 1;      /* segA even, F odd */
  if (LEVELS) __atomic_fetch_or(lev[j] + (i >> 6), 1ULL << (i & 63), __ATOMIC_RELAXED);
  bump(i, j <= K - 1 ? 2 : 1);
}
/* level-K subtree: only O and Z (no root change left). Precondition: F(D,E,U,d) <= segHI.
   F(O s) = D + 2U (= F(J s)), F(Z s) = F + 2U - 2d; children are entered only if their F <= segHI. */
static void dfsK(ll D, ll E, ll U, ll d){
  for (;;){
    ll F = D - E + U;
    if (F >= segA) record(F, K);
    if (D + 2*U <= segHI) dfsK(D, 0, 2*U, d);           /* O */
    if (F + 2*U - 2*d > segHI) return;                    /* Z */
    E = d; U = 2*U - d; D = F;
  }
}
static void dfs(ll D, ll E, ll U, ll d, int j){            /* j < K, F <= segHI */
  for (;;){
    ll F = D - E + U;
    if (F >= segA) record(F, j);
    if (D + 2*U <= segHI){
      if (j + 1 < K) dfs(F, D - E, D + U, D, j + 1); else dfsK(F, D - E, D + U, D);   /* J */
      dfs(D, 0, 2*U, d, j);                                                           /* O */
    }
    if (F + 2*U - 2*d > segHI) return;                                               /* Z */
    E = d; U = 2*U - d; D = F;
  }
}
static void dfs_entry(ll D, ll E, ll U, ll d, int j){
  if (D - E + U > segHI) return;
  if (j == K) dfsK(D, E, U, d); else dfs(D, E, U, d, j);
}
typedef struct { ll D, E, U, d; int j; } node;
static node *T; static ll nT, capT;
static void gen(ll D, ll E, ll U, ll d, int j){
  ll F = D - E + U;
  if (F > segHI) return;
  if (U >= USPLIT){ if (nT == capT){ capT = capT ? 2*capT : 4096; T = realloc(T, capT * sizeof(node)); } T[nT++] = (node){D,E,U,d,j}; return; }
  record(F, j);
  if (j < K) gen(F, D - E, D + U, D, j + 1);
  gen(D, 0, 2*U, d, j);
  gen(F, d, 2*U - d, d, j);
}

/* accumulated statistics */
#define NW 64
static ll w_n[NW], w_min[NW], w_min3[NW][3], w_c[NW][9], w_sat[NW], w_sat255[NW]; static double w_sum255[NW], w_sumx[NW];
static double w_v3s[NW][7], w_v3n[NW][7], w_m9s[NW][9], w_m9n[NW][9];
static ll lev_miss[8], lev_last[8], lev_wmiss[NW][8];
static ll un_miss = 0, un_last = 0;
static ll *L_miss; static ll nL_miss = 0; static ll *L_low2; static int *L_low2v; static ll nL_low2 = 0;
static ll *L_low12; static int *L_low12v; static ll nL_low12 = 0; static ll maxM = 0;
static int wlog(ll N){ int k = 0; while ((2LL << k) <= N) k++; return k; }

int main(int argc, char **argv){
  if (argc < 6){ fprintf(stderr, "usage: census_audit X K SEGLOG CNTBITS LEVELS [LISTLO] [USPLIT]\n"); return 1; }
  X = atoll(argv[1]); K = atoi(argv[2]); SEGLOG = atoi(argv[3]); CNTBITS = atoi(argv[4]); LEVELS = atoi(argv[5]);
  LISTLO = argc > 6 ? atoll(argv[6]) : 100001; USPLIT = argc > 7 ? atoll(argv[7]) : 512; ll SEG0 = argc > 8 ? atoll(argv[8]) : 0;
  if (K < 0 || K > 7){ fprintf(stderr, "K in 0..7\n"); return 1; }
  CAP = CNTBITS == 32 ? (1LL << 32) - 1 : (1LL << CNTBITS) - 1;
  ll S = 1LL << SEGLOG, H = S >> 1;
  if (CNTBITS == 4) c4 = malloc(H / 2 + 1); else if (CNTBITS == 8) c8 = malloc(H); else if (CNTBITS == 16) c16 = malloc(2 * H); else c32 = malloc(4 * H);
  if (LEVELS) for (int j = 0; j <= K; j++) lev[j] = malloc((H / 64 + 1) * 8);
  L_miss = malloc(sizeof(ll) * 1000000); L_low2 = malloc(sizeof(ll) * 1000000); L_low2v = malloc(sizeof(int) * 1000000);
  L_low12 = malloc(sizeof(ll) * 1000000); L_low12v = malloc(sizeof(int) * 1000000);
  for (int w = 0; w < NW; w++){ w_min[w] = 1LL << 40; for (int r = 0; r < 3; r++) w_min3[w][r] = 1LL << 40; }
  double t0 = omp_get_wtime(); ll ntasks_total = 0;
  for (segA = SEG0; segA <= X; segA += S){
    segHI = segA + S - 1; if (segHI > X) segHI = X;
    if (CNTBITS == 4) memset(c4, 0, H / 2 + 1); else if (CNTBITS == 8) memset(c8, 0, H); else if (CNTBITS == 16) memset(c16, 0, 2 * H); else memset(c32, 0, 4 * H);
    if (LEVELS) for (int j = 0; j <= K; j++) memset(lev[j], 0, (H / 64 + 1) * 8);
    if (segA == 0){ bump(0, 1); if (LEVELS) lev[0][0] |= 1; }   /* the empty word: N = 1 at level 0 */
    nT = 0; gen(1, 0, 2, 1, 0); ntasks_total += nT;
    #pragma omp parallel for schedule(dynamic, 1)
    for (ll q = 0; q < nT; q++) dfs_entry(T[q].D, T[q].E, T[q].U, T[q].d, T[q].j);
    /* harvest */
    ll lo = segA | 1;
    for (ll N = lo; N <= segHI; N += 2){
      ll i = (N - segA) >> 1; ll m = getc_(i); int w = wlog(N);
      if (m > maxM) maxM = m;
      w_n[w]++; if (m < w_min[w]) w_min[w] = m; if (m < w_min3[w][N % 3]) w_min3[w][N % 3] = m;
      if (m <= 8) w_c[w][m]++;
      if (m >= CAP) w_sat[w]++;
      if (CAP >= 255){ ll ms = m > 255 ? 255 : m; w_sum255[w] += ms; if (ms == 255) w_sat255[w]++; }
      w_sumx[w] += m;
      if (CNTBITS >= 16){ int k3 = 0; ll v = N; while (v % 3 == 0 && k3 < 6){ v /= 3; k3++; } w_v3s[w][k3] += m; w_v3n[w][k3]++; w_m9s[w][N % 9] += m; w_m9n[w][N % 9]++; }
      if (m == 0){ un_miss++; un_last = N; if (N >= LISTLO && nL_miss < 1000000) L_miss[nL_miss++] = N; }
      if (N > 100000 && m <= 2 && nL_low2 < 1000000){ L_low2[nL_low2] = N; L_low2v[nL_low2++] = (int)m; }
      if (N > 4194304 && m <= 12 && CAP >= 12 && nL_low12 < 1000000){ L_low12[nL_low12] = N; L_low12v[nL_low12++] = (int)m; }
      if (LEVELS){ int cov = 0; for (int j = 0; j <= K; j++){ cov |= (int)((lev[j][i >> 6] >> (i & 63)) & 1);
          if (!cov){ lev_miss[j]++; lev_last[j] = N; lev_wmiss[w][j]++; } } }
    }
    fprintf(stderr, "segment [%lld,%lld] done, tasks %lld, %.1fs\n", segA, segHI, nT, omp_get_wtime() - t0);
  }
  double el = omp_get_wtime() - t0;
  printf("census_audit X=%lld K=%d SEGLOG=%d CNTBITS=%d LEVELS=%d threads=%d tasks=%lld time=%.1fs\n", X, K, SEGLOG, CNTBITS, LEVELS, omp_get_max_threads(), ntasks_total, el);
  if (LEVELS) for (int j = 0; j <= K; j++) printf("<=%d root changes: odd misses %lld, largest %lld\n", j, lev_miss[j], lev_last[j]);
  printf("union (M_%d = 0): odd misses %lld, largest %lld\n", K, un_miss, un_last);
  printf("misses (<=%d) from %lld:", K, LISTLO); for (ll q = 0; q < nL_miss; q++) printf(" %lld", L_miss[q]); printf("\n");
  if (LEVELS) for (int w = 0; w < NW; w++) if (w_n[w]){ printf("  [%lld,%lld):", 1LL << w, (2LL << w) > X ? X : (2LL << w)); for (int j = 0; j <= K; j++) printf(" k%d:%lld", j, lev_wmiss[w][j]); printf("\n"); }
  printf("multiplicity M_%d per window (counter cap %lld; mean255 = mean of min(M,255); max seen %lld):\n", K, CAP, maxM);
  for (int w = 0; w < NW; w++) if (w_n[w]){
    printf("  [%lld,%lld): n=%lld min=%lld min_mod3=(%lld,%lld,%lld) satcap=%lld", 1LL << w, (2LL << w) > X ? X : (2LL << w), w_n[w], w_min[w],
           w_min3[w][0] == (1LL<<40) ? 999 : w_min3[w][0], w_min3[w][1] == (1LL<<40) ? 999 : w_min3[w][1], w_min3[w][2] == (1LL<<40) ? 999 : w_min3[w][2], w_sat[w]);
    if (CAP >= 255) printf(" mean255=%.2f sat255=%lld", w_sum255[w] / w_n[w], w_sat255[w]);
    if (CNTBITS >= 16) printf(" exactmean=%.2f", w_sumx[w] / w_n[w]);
    printf(" c0..8="); for (int m = 0; m <= 8; m++) printf("%s%lld", m ? "," : "", w_c[w][m]); printf("\n");
  }
  if (CNTBITS >= 16) for (int w = 20; w < NW; w++) if (w_n[w]){
    printf("  [%lld,%lld) mean by v3(N)=0..5,>=6:", 1LL << w, 2LL << w); for (int k = 0; k < 7; k++) printf(" %.1f", w_v3n[w][k] ? w_v3s[w][k] / w_v3n[w][k] : 0.0);
    double s3[3] = {0}, n3[3] = {0}; for (int r = 0; r < 9; r++){ s3[r % 3] += w_m9s[w][r]; n3[r % 3] += w_m9n[w][r]; }
    printf(" | mod 3: %.2f %.2f %.2f | mod 9:", s3[0] / n3[0], s3[1] / n3[1], s3[2] / n3[2]); for (int r = 0; r < 9; r++) printf(" %.1f", w_m9s[w][r] / w_m9n[w][r]); printf("\n");
  }
  printf("odd N > 10^5 with M <= 2:"); for (ll q = 0; q < nL_low2; q++) printf(" %lld:%d", L_low2[q], L_low2v[q]); printf("\n");
  if (CAP >= 12){ printf("odd N > 2^22 with M <= 12 (%lld values):", nL_low12); for (ll q = 0; q < nL_low12; q++) printf(" %lld:%d", L_low12[q], L_low12v[q]); printf("\n"); }
  return 0;
}
