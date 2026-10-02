// anneal.c -- simulated annealing for edge-multiset resolving k-sets of Q_d (d <= 7, k <= 31).
// Exact keys (histogram levels 0..d-2 in 5-bit fields, as in edimsearch.c); cost = number of colliding
// edge pairs = sum over histogram classes of C(m,2).  Moves: swap a landmark s in S with a vertex t not
// in S; with probability pfocus the move is "focused": pick a random colliding edge pair (e,f) and a
// vertex t with d(e,t) != d(f,t) (t separates e and f).  Metropolis acceptance, geometric cooling per
// restart.  Any cost-0 set is re-verified from scratch (sorting full histograms) before it is printed.
// usage: anneal d k seconds seed [iters_per_restart T0 T1 pfocus]
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#define MAXE 448
static int D, N, N1, E, K;
static uint32_t C[128][MAXE];
static uint8_t DIST[MAXE][128];
static int EU[MAXE], EI[MAXE];
static uint64_t rs;
static inline uint64_t rng(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static inline double ur(void) { return (rng() >> 11) * (1.0 / 9007199254740992.0); }
// cost via hash table with counts; also records up to 64 colliding pairs for focused moves
static uint32_t hk[2048], hs[2048], hc[2048]; static uint16_t he[2048]; static uint32_t gen = 0;
static int cp_e[64], cp_f[64], ncp;
static long cost(const uint32_t *key) {
  if (++gen == 0) { memset(hs, 0, sizeof hs); gen = 1; }
  long c = 0; ncp = 0;
  for (int e = 0; e < E; e++) {
    uint32_t k = key[e], h = (k * 2654435761u) >> 21;
    while (hs[h] == gen && hk[h] != k) h = (h + 1) & 2047;
    if (hs[h] == gen) { c += hc[h]; hc[h]++; if (ncp < 64) { cp_e[ncp] = he[h]; cp_f[ncp] = e; ncp++; } }
    else { hs[h] = gen; hk[h] = k; hc[h] = 1; he[h] = (uint16_t)e; }
  }
  return c;
}
static int cmp8(const void *a, const void *b) { return memcmp(a, b, 8); }
static int exact_resolving(const int *S, int k) {   // from scratch: d(e,s) = min over endpoints, full histograms sorted
  uint8_t (*H)[8] = calloc(E, 8);
  for (int e = 0; e < E; e++) { int u = EU[e], v = EU[e] | (1 << EI[e]);
    for (int j = 0; j < k; j++) { int du = __builtin_popcount(u ^ S[j]), dv = __builtin_popcount(v ^ S[j]); H[e][du < dv ? du : dv]++; } }
  qsort(H, E, 8, cmp8); int ok = 1; for (int e = 1; e < E; e++) if (!memcmp(H[e], H[e - 1], 8)) { ok = 0; break; }
  free(H); return ok;
}
int main(int argc, char **argv) {
  if (argc < 5) { fprintf(stderr, "usage: anneal d k seconds seed [iters_per_restart T0 T1 pfocus]\n"); return 1; }
  D = atoi(argv[1]); K = atoi(argv[2]); double secs = atof(argv[3]); rs = strtoull(argv[4], 0, 10) * 0x9E3779B97F4A7C15ULL + 12345;
  long iters = argc > 5 ? atol(argv[5]) : 2000000; double T0 = argc > 6 ? atof(argv[6]) : 2.0, T1 = argc > 7 ? atof(argv[7]) : 0.05;
  double pfocus = argc > 8 ? atof(argv[8]) : 0.5;
  if (D < 2 || D > 7 || K < 1 || K > 31) return 1;
  N = 1 << D; N1 = N / 2; E = 0;
  for (int u = 0; u < N; u++) for (int i = 0; i < D; i++) if (!(u >> i & 1)) { EU[E] = u; EI[E] = i; E++; }
  for (int e = 0; e < E; e++) for (int s = 0; s < N; s++) { int d = __builtin_popcount((EU[e] ^ s) & ~(1 << EI[e])); DIST[e][s] = d; C[s][e] = d <= D - 2 ? 1u << (5 * d) : 0; }
  struct timespec t0, t1; clock_gettime(CLOCK_MONOTONIC, &t0);
  int S[64], in[128], best_S[64]; long best = 1L << 60; long restarts = 0, totmoves = 0;
  uint32_t key[MAXE], nk[MAXE];
  for (;;) {
    clock_gettime(CLOCK_MONOTONIC, &t1);
    double el = (t1.tv_sec - t0.tv_sec) + 1e-9 * (t1.tv_nsec - t0.tv_nsec);
    if (el > secs) break;
    memset(in, 0, sizeof in); int c = 0;
    while (c < K) { int v = rng() % N; if (!in[v]) { in[v] = 1; S[c++] = v; } }
    for (int e = 0; e < E; e++) { uint32_t x = 0; for (int j = 0; j < K; j++) x += C[S[j]][e]; key[e] = x; }
    long cur = cost(key), rbest = cur;
    int lcp_e[64], lcp_f[64], lncp = ncp; memcpy(lcp_e, cp_e, sizeof cp_e); memcpy(lcp_f, cp_f, sizeof cp_f);
    for (long it = 0; it < iters && cur > 0; it++) {
      double T = T0 * pow(T1 / T0, (double)it / iters);
      int a = rng() % K, s = S[a], t;
      if (lncp > 0 && ur() < pfocus) {   // focused: t separates a random colliding pair
        int q = rng() % lncp, e = lcp_e[q], f = lcp_f[q], tries = 0;
        do { t = rng() % N; tries++; } while ((in[t] || DIST[e][t] == DIST[f][t]) && tries < 200);
        if (in[t] || DIST[e][t] == DIST[f][t]) continue;
      } else { do { t = rng() % N; } while (in[t]); }
      for (int e = 0; e < E; e++) nk[e] = key[e] - C[s][e] + C[t][e];
      long nc = cost(nk); totmoves++;
      if (nc <= cur || ur() < exp(-(double)(nc - cur) / T)) {
        memcpy(key, nk, sizeof(uint32_t) * E); cur = nc; S[a] = t; in[s] = 0; in[t] = 1;
        lncp = ncp; memcpy(lcp_e, cp_e, sizeof cp_e); memcpy(lcp_f, cp_f, sizeof cp_f);
        if (cur < rbest) rbest = cur;
        if (cur < best) { best = cur; memcpy(best_S, S, sizeof(int) * K); }
      }
    }
    restarts++;
    if (cur < best) { best = cur; memcpy(best_S, S, sizeof(int) * K); }
    if (cur == 0) {
      if (!exact_resolving(S, K)) { printf("FATAL: cost 0 but exact check fails\n"); return 9; }
      printf("RESOLVING d=%d k=%d set=", D, K); for (int j = 0; j < K; j++) printf(j ? ",%d" : "%d", S[j]); printf("\n"); fflush(stdout);
    }
    printf("restart %ld: final cost %ld, best in run %ld, overall best %ld\n", restarts, cur, rbest, best); fflush(stdout);
  }
  printf("SUMMARY d=%d k=%d restarts=%ld moves=%ld best_colliding_pairs=%ld best_set=", D, K, restarts, totmoves, best);
  for (int j = 0; j < K; j++) printf(j ? ",%d" : "%d", best_S[j]);
  printf("\n");
  return 0;
}
