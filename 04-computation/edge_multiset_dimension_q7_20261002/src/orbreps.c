// orbreps.c -- Aut(Q_n)-orbit representatives of a-subsets of V(Q_n) (n <= 6) by orderly generation.
//
// Canonical form of a set T: the lexicographically smallest SORTED tuple in its Aut(Q_n)-orbit
// (Aut(Q_n) = { x -> pi(x xor t) }, order 2^n n!).  This canonical form is hereditary: if T is
// canonical then T minus its largest element is canonical (proof in q7_report.md), so every canonical
// (a+1)-tuple is a canonical a-tuple plus one element larger than its maximum.  Level a+1 is therefore
// obtained from level a by trying all such one-element extensions and keeping the canonical ones;
// every orbit is produced exactly once.
//
// Canonicity test of T (mask M): T must contain 0; for every s in T translate s to 0 and minimise over
// all n! coordinate permutations (enumerated by the Steinhaus-Johnson-Trotter sequence of adjacent
// transpositions, each applied to the 2^n-bit mask as one delta swap).  Lex order of sorted tuples of
// equal size = "lowest bit of X xor Y lies in X".  Shortcut (exact): with t2 the 2nd smallest element of
// T and w_s the distance from s to its nearest other element of T, translations s with 2^w_s - 1 > t2
// cannot produce a smaller image, and 2^w_s - 1 < t2 proves T is not canonical.
//
// usage: orbreps n amax outdir [selftest]
//   writes outdir/reps_n<n>_a<a>.bin (little-endian uint64 masks, lexicographic tuple order), a=0..amax,
//   and prints the number of representatives per level.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

static int n, N, nperm;
static int sjt_pos[720];
static uint64_t swapmask[6], FULL;
static const uint64_t TM[6] = {0x5555555555555555ULL, 0x3333333333333333ULL, 0x0F0F0F0F0F0F0F0FULL,
                               0x00FF00FF00FF00FFULL, 0x0000FFFF0000FFFFULL, 0x00000000FFFFFFFFULL};

static inline uint64_t translate(uint64_t m, int t) {     // bit x -> bit x^t
  for (int j = 0; j < n; j++) if (t >> j & 1) { int s = 1 << j; m = ((m & TM[j]) << s) | ((m >> s) & TM[j]); }
  return m & FULL;
}
static inline uint64_t adjswap(uint64_t m, int i) {      // exchange coordinates i and i+1
  int s = 1 << i; uint64_t t = ((m >> s) ^ m) & swapmask[i]; return m ^ t ^ (t << s);
}
static inline int lexless(uint64_t X, uint64_t Y) {      // sorted tuple of X < sorted tuple of Y (|X|=|Y|)
  uint64_t z = X ^ Y; return z && ((X & z & (~z + 1)) != 0);
}
static void build_sjt(void) {
  int p[6], dir[6]; nperm = 1;
  for (int i = 0; i < n; i++) { p[i] = i; dir[i] = -1; nperm *= (i + 1); }
  int cnt = 0;
  for (;;) {
    int mi = -1;
    for (int i = 0; i < n; i++) {
      int j = i + dir[p[i]];
      if (j >= 0 && j < n && p[j] < p[i] && (mi < 0 || p[i] > p[mi])) mi = i;
    }
    if (mi < 0) break;
    int j = mi + dir[p[mi]], v = p[mi];
    p[mi] = p[j]; p[j] = v; sjt_pos[cnt++] = mi < j ? mi : j;
    for (int i = 0; i < n; i++) if (p[i] > v) dir[p[i]] = -dir[p[i]];
  }
  if (cnt != nperm - 1) { fprintf(stderr, "SJT length %d != %d\n", cnt, nperm - 1); exit(2); }
  for (int i = 0; i + 1 < n; i++) {
    swapmask[i] = 0;
    for (int x = 0; x < N; x++) if ((x >> i & 1) && !(x >> (i + 1) & 1)) swapmask[i] |= 1ULL << x;
  }
  // self-check: the cumulative products visit nperm distinct coordinate permutations
  // (probe: images of the n unit-vector singletons)
  uint64_t *seen = malloc(sizeof(uint64_t) * nperm * 6); int ns = 0;
  uint64_t pr[6]; for (int c = 0; c < n; c++) pr[c] = 1ULL << (1 << c);
  for (int q = 0; q < nperm; q++) {
    if (q > 0) for (int c = 0; c < n; c++) pr[c] = adjswap(pr[c], sjt_pos[q - 1]);
    for (int r = 0; r < ns; r++) if (!memcmp(seen + 6 * r, pr, sizeof(uint64_t) * n)) { fprintf(stderr, "SJT repeat\n"); exit(2); }
    memcpy(seen + 6 * ns, pr, sizeof(uint64_t) * n); ns++;
    for (int c = 0; c < n; c++) if (__builtin_popcountll(pr[c]) != 1 || __builtin_popcountll(__builtin_ctzll(pr[c])) != 1) { fprintf(stderr, "SJT not a coordinate perm\n"); exit(2); }
  }
  free(seen);
}
static int is_canon(uint64_t M) {
  int a = __builtin_popcountll(M);
  if (a == 0) return 1;
  if (!(M & 1)) return 0;
  if (a == 1) return 1;
  int T[64], c = 0; for (uint64_t x = M; x; x &= x - 1) T[c++] = __builtin_ctzll(x);
  int t2 = T[1];
  for (int j = 0; j < a; j++) {
    int s = T[j], ws = 99;
    for (int l = 0; l < a; l++) if (l != j) { int w = __builtin_popcount(T[l] ^ s); if (w < ws) ws = w; }
    int lo = (1 << ws) - 1;
    if (lo > t2) continue;
    if (lo < t2) return 0;
    uint64_t m = translate(M, s);
    if (lexless(m, M)) return 0;
    for (int q = 0; q < nperm - 1; q++) { m = adjswap(m, sjt_pos[q]); if (lexless(m, M)) return 0; }
  }
  return 1;
}
// brute-force canonical form via explicit vertex permutations (self-test only)
static uint64_t canon_brute(uint64_t M) {
  int perms[720][6]; int np = 0, q[6];
  // all permutations of n by recursion-free counting
  int idx[6] = {0};
  for (;;) {
    int used = 0, ok = 1;
    for (int i = 0; i < n; i++) { q[i] = idx[i]; if (used >> q[i] & 1) ok = 0; used |= 1 << q[i]; }
    if (ok) { memcpy(perms[np], q, sizeof q); np++; }
    int i = 0; while (i < n && ++idx[i] == n) { idx[i] = 0; i++; } if (i == n) break;
  }
  uint64_t best = 0; int have = 0;
  for (int p = 0; p < np; p++) for (int t = 0; t < N; t++) {
    uint64_t img = 0;
    for (uint64_t x = M; x; x &= x - 1) { int v = __builtin_ctzll(x) ^ t, u = 0; for (int i = 0; i < n; i++) if (v >> i & 1) u |= 1 << perms[p][i]; img |= 1ULL << u; }
    if (!have || lexless(img, best)) { best = img; have = 1; }
  }
  return best;
}
static uint64_t rng_s = 88172645463325252ULL;
static uint64_t rng(void) { rng_s ^= rng_s << 13; rng_s ^= rng_s >> 7; rng_s ^= rng_s << 17; return rng_s; }

int main(int argc, char **argv) {
  if (argc < 4) { fprintf(stderr, "usage: orbreps n amax outdir [selftest]\n"); return 1; }
  n = atoi(argv[1]); int amax = atoi(argv[2]); const char *outdir = argv[3];
  if (n < 1 || n > 6) return 1;
  N = 1 << n; FULL = (N == 64) ? ~0ULL : ((1ULL << N) - 1);
  build_sjt();
  if (argc > 4) {   // self-test: is_canon(M) == (canon_brute(M) == M) on random sets of each size
    int bad = 0, tests = atoi(argv[4]);
    for (int t = 0; t < tests; t++) {
      int a = 1 + rng() % (N - 1); uint64_t M = 0;
      while (__builtin_popcountll(M) < a) M |= 1ULL << (rng() % N);
      uint64_t c = canon_brute(M);
      if (is_canon(M) != (c == M)) bad++;
      if (!is_canon(c)) bad++;
    }
    printf("selftest n=%d: %d random sets, %d disagreements\n", n, tests, bad);
    return bad != 0;
  }
  uint64_t *cur = malloc(8), *nxt; size_t ncur = 1; cur[0] = 0;
  for (int a = 0; a <= amax; a++) {
    char fn[512]; snprintf(fn, sizeof fn, "%s/reps_n%d_a%d.bin", outdir, n, a);
    FILE *f = fopen(fn, "wb"); if (!f) { perror(fn); return 3; }
    if (fwrite(cur, 8, ncur, f) != ncur) { perror("write"); return 3; }
    fclose(f);
    printf("n=%d a=%d reps=%zu  (t=%.1fs)\n", n, a, ncur, (double)clock() / CLOCKS_PER_SEC); fflush(stdout);
    if (a == amax) break;
    size_t cap = 1024, nn = 0; nxt = malloc(8 * cap);
    for (size_t r = 0; r < ncur; r++) {
      uint64_t M = cur[r]; int top = M ? 63 - __builtin_clzll(M) : -1;
      for (int x = top + 1; x < N; x++) {
        uint64_t M2 = M | (1ULL << x);
        if (is_canon(M2)) { if (nn == cap) { cap *= 2; nxt = realloc(nxt, 8 * cap); } nxt[nn++] = M2; }
      }
    }
    free(cur); cur = nxt; ncur = nn;
  }
  return 0;
}
