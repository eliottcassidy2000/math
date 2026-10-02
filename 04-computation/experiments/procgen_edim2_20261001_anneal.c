// procgen_edim2_20261001_anneal.c -- memory-light descending annealing for edge-multiset resolving sets of Q_d.
// Distances are computed on the fly by the projection lemma: d(e,w) = popcount((u ^ w) & ~(1<<i)) for
// e = {u, u + e_i}.  Histograms are hashed additively, key(e) = sum_{s in S} R[d(e,s)] (64-bit Zobrist);
// false hash collisions only make the search conservative.  Every reported set is re-verified EXACTLY
// (full histograms, sorted and compared) before it is printed.
// usage: anneal d k0 iters_per_level seed T0 T1 kmin
// output lines: "d=.. k=.. FOUND iters=.. set: v1 v2 ..." (exactly verified) or "... FAILED ...".
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
static int D, N, E;
static uint32_t *EU; static uint8_t *EI;          // edge e = (EU[e], EI[e])
static uint64_t R[64]; static uint64_t rs;
static inline uint64_t rng(void){ rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static inline double ur(void){ return (rng() >> 11) * (1.0 / 9007199254740992.0); }
static inline int dist(int e, int w){ return __builtin_popcount((EU[e] ^ (uint32_t)w) & ~(1u << EI[e])); }
static uint64_t *hk; static int *hc, *hu; static int HS, nu;
static long cost(const uint64_t *key){
  long c = 0; nu = 0;
  for (int e = 0; e < E; e++){
    uint64_t k = key[e]; int h = (int)((k * 0x9E3779B97F4A7C15ULL) >> 40) & (HS - 1);
    while (1){
      if (hc[h] == 0){ hk[h] = k; hc[h] = 1; hu[nu++] = h; break; }
      if (hk[h] == k){ c += hc[h]; hc[h]++; break; }
      h = (h + 1) & (HS - 1);
    }
  }
  for (int i = 0; i < nu; i++) hc[hu[i]] = 0;
  return c;
}
static int HW; // histogram width in 16-bit counters
static int cmpH(const void *a, const void *b){ return memcmp(a, b, (size_t)HW * 2); }
static int exact_ok(const int *S, int k){
  uint16_t *H = calloc((size_t)E * HW, 2);
  for (int e = 0; e < E; e++) for (int j = 0; j < k; j++) H[(size_t)e * HW + dist(e, S[j])]++;
  // store big-endian-ish ordering irrelevant: only equality matters after sorting
  qsort(H, E, (size_t)HW * 2, cmpH);
  int ok = 1;
  for (int e = 1; e < E; e++) if (!memcmp(H + (size_t)e * HW, H + (size_t)(e - 1) * HW, (size_t)HW * 2)){ ok = 0; break; }
  free(H); return ok;
}
int main(int argc, char **argv){
  if (argc < 8){ fprintf(stderr, "usage: anneal d k0 iters seed T0 T1 kmin\n"); return 1; }
  D = atoi(argv[1]); int K = atoi(argv[2]); long iters = atol(argv[3]);
  rs = strtoull(argv[4], 0, 10) * 0x9E3779B97F4A7C15ULL + 1;
  double T0 = atof(argv[5]), T1 = atof(argv[6]); int kmin = atoi(argv[7]);
  N = 1 << D; E = D * (N / 2); HW = D;
  EU = malloc(4 * (size_t)E); EI = malloc((size_t)E);
  int e = 0;
  for (int u = 0; u < N; u++) for (int i = 0; i < D; i++) if (!(u >> i & 1)){ EU[e] = u; EI[e] = i; e++; }
  HS = 1; while (HS < 4 * E) HS <<= 1;
  hk = malloc(8 * (size_t)HS); hc = calloc(HS, 4); hu = malloc(4 * (size_t)HS);
  for (int r = 0; r < 64; r++) R[r] = rng() | 1;
  uint64_t *key = malloc(8 * (size_t)E), *nk = malloc(8 * (size_t)E);
  int *S = malloc(4 * (size_t)N), *out = malloc(4 * (size_t)N), *in = calloc(N, 4);
  int c = 0; while (c < K){ int v = rng() % N; if (!in[v]){ in[v] = 1; S[c++] = v; } }
  int k = K;
  while (k >= kmin){
    int no = 0; for (int v = 0; v < N; v++) if (!in[v]) out[no++] = v;
    for (int f = 0; f < E; f++){ uint64_t x = 0; for (int j = 0; j < k; j++) x += R[dist(f, S[j])]; key[f] = x; }
    long cur = cost(key); long it;
    for (it = 0; it < iters && cur > 0; it++){
      double T = T0 * pow(T1 / T0, (double)it / iters);
      int a = rng() % k, b = rng() % no, s = S[a], t = out[b];
      for (int f = 0; f < E; f++) nk[f] = key[f] - R[dist(f, s)] + R[dist(f, t)];
      long nc = cost(nk);
      if (nc <= cur || ur() < exp(-(nc - cur) / T)){
        uint64_t *tmp = key; key = nk; nk = tmp; cur = nc; S[a] = t; out[b] = s; in[s] = 0; in[t] = 1;
      }
    }
    if (cur > 0){ printf("d=%d k=%d FAILED (cost %ld after %ld iters)\n", D, k, cur, it); fflush(stdout); break; }
    if (!exact_ok(S, k)){ printf("d=%d k=%d hash-ok but exact FAIL\n", D, k); fflush(stdout); break; }
    printf("d=%d k=%d FOUND iters=%ld set:", D, k, it);
    for (int j = 0; j < k; j++) printf(" %d", S[j]);
    printf("\n"); fflush(stdout);
    // remove the landmark whose removal leaves the fewest collisions
    long bestc = -1; int besta = 0;
    for (int a = 0; a < k; a++){
      for (int f = 0; f < E; f++) nk[f] = key[f] - R[dist(f, S[a])];
      long cc = cost(nk); if (bestc < 0 || cc < bestc){ bestc = cc; besta = a; }
    }
    in[S[besta]] = 0; S[besta] = S[k - 1]; k--;
  }
  return 0;
}
