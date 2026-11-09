/* S19 audit: core statistics of every ported tournament of small order, generated independently.
 *
 * mode "classes NMAX": isomorphism classes of tournaments of order 1..NMAX (NMAX <= 8) by one-vertex extension of the
 *   previous level's representatives plus a brute-force canonical form (minimum arc code over all n! relabellings).
 *   Prints "#count n classes" lines, then one line per class representative and port:
 *     n cls u  h2 R2 D2 hu2  h3 R3 D3 hu3
 *   with h = I(Omega(T),x), R = E + O, D = E - O, hu = I(Omega(T - u), x) from tcount.h.
 * mode "labelled N": every labelled tournament of order N with port 0 (every ported tournament is isomorphic to one
 *   with port 0); prints the distinct tuples "x h R D hu" for x = 2 and x = 3 (no isomorphism testing at all).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tcount.h"

static int PERM[40320][8], NPERM;
static void gen_perms(int n) {
  int p[8]; for (int i = 0; i < n; i++) p[i] = i; NPERM = 0;
  /* lexicographic next_permutation */
  while (1) {
    memcpy(PERM[NPERM++], p, sizeof p);
    int i = n - 2; while (i >= 0 && p[i] > p[i + 1]) i--; if (i < 0) break;
    int j = n - 1; while (p[j] < p[i]) j--; int t = p[i]; p[i] = p[j]; p[j] = t;
    for (int a = i + 1, b = n - 1; a < b; a++, b--) { t = p[a]; p[a] = p[b]; p[b] = t; }
  }
}
/* arc code: bit index over pairs (a<b) in row order; bit = 1 iff a -> b */
static unsigned long long code_of(int n, const uint32_t *out, const int *p) {
  /* p[v] = new label of old vertex v; build new out masks */
  uint32_t no[8] = {0};
  for (int v = 0; v < n; v++) for (int w = 0; w < n; w++) if (out[v] >> w & 1) no[p[v]] |= 1u << p[w];
  unsigned long long c = 0; int b = 0;
  for (int a = 0; a < n; a++) for (int d = a + 1; d < n; d++, b++) if (no[a] >> d & 1) c |= 1ULL << b;
  return c;
}
static void decode(int n, unsigned long long c, uint32_t *out) {
  for (int v = 0; v < n; v++) out[v] = 0; int b = 0;
  for (int a = 0; a < n; a++) for (int d = a + 1; d < n; d++, b++) { if (c >> b & 1) out[a] |= 1u << d; else out[d] |= 1u << a; }
}
static int cmpull(const void *a, const void *b) { unsigned long long x = *(const unsigned long long *)a, y = *(const unsigned long long *)b; return x < y ? -1 : x > y; }

static void stats_line(int n, const uint32_t *out, int cls) {
  static ll I2[1 << MAXN], I3[1 << MAXN];
  TN = n; tc_alloc(n); memcpy(TOUT, out, n * sizeof(uint32_t));
  tc_cycles(); tc_packing(2, I2); tc_packing(3, I3);
  size_t FULL = ((size_t)1 << n) - 1;
  for (int u = 0; u < n; u++) {
    ll E2, O2, E3, O3; tc_port(u, I2, &E2, &O2); tc_port(u, I3, &E3, &O3);
    printf("%d %d %d  %lld %lld %lld %lld  %lld %lld %lld %lld\n", n, cls, u,
           I2[FULL], E2 + O2, E2 - O2, I2[FULL & ~((size_t)1 << u)],
           I3[FULL], E3 + O3, E3 - O3, I3[FULL & ~((size_t)1 << u)]);
  }
}

/* small open-addressing set of 4-tuples */
typedef struct { ll k[4]; int used; } slot;
static slot *ST[2]; static size_t SSZ = 1 << 16, SCNT[2];
static void set_add(int which, ll a, ll b, ll c, ll d) {
  if (!ST[which]) ST[which] = calloc(SSZ, sizeof(slot));
  unsigned long long h = (unsigned long long)a * 1000003ULL ^ (unsigned long long)b * 998244353ULL ^ (unsigned long long)(c + 1000000) * 19260817ULL ^ (unsigned long long)d * 7919ULL;
  size_t i = (h ^ (h >> 29)) & (SSZ - 1);
  while (ST[which][i].used) {
    slot *s = &ST[which][i];
    if (s->k[0] == a && s->k[1] == b && s->k[2] == c && s->k[3] == d) return;
    i = (i + 1) & (SSZ - 1);
  }
  ST[which][i].used = 1; ST[which][i].k[0] = a; ST[which][i].k[1] = b; ST[which][i].k[2] = c; ST[which][i].k[3] = d;
  if (++SCNT[which] * 2 > SSZ) { fprintf(stderr, "set too full\n"); exit(3); }
}

int main(int argc, char **argv) {
  if (argc >= 3 && !strcmp(argv[1], "classes")) {
    int NMAX = atoi(argv[2]);
    unsigned long long *prev = malloc(sizeof(unsigned long long)); prev[0] = 0; size_t nprev = 1; /* order 1 */
    printf("#count 1 1\n");
    { uint32_t out[1] = {0}; stats_line(1, out, 0); }
    for (int n = 2; n <= NMAX; n++) {
      gen_perms(n);
      size_t cap = nprev << (n - 1); unsigned long long *cand = malloc(cap * sizeof(unsigned long long)); size_t nc = 0;
      for (size_t i = 0; i < nprev; i++) {
        uint32_t base[8]; decode(n - 1, prev[i], base);
        for (unsigned m = 0; m < (1u << (n - 1)); m++) {
          uint32_t out[8];
          for (int v = 0; v < n - 1; v++) out[v] = base[v] | ((m >> v & 1) ? 0u : (1u << (n - 1)));
          out[n - 1] = 0; for (int v = 0; v < n - 1; v++) if (m >> v & 1) out[n - 1] |= 1u << v; /* new vertex beats v iff bit */
          unsigned long long best = ~0ULL;
          for (int k = 0; k < NPERM; k++) { unsigned long long c = code_of(n, out, PERM[k]); if (c < best) best = c; }
          cand[nc++] = best;
        }
      }
      qsort(cand, nc, sizeof *cand, cmpull);
      size_t nu = 0; for (size_t i = 0; i < nc; i++) if (!nu || cand[i] != cand[nu - 1]) cand[nu++] = cand[i];
      printf("#count %d %zu\n", n, nu);
      for (size_t i = 0; i < nu; i++) { uint32_t out[8]; decode(n, cand[i], out); stats_line(n, out, (int)i); }
      free(prev); prev = cand; nprev = nu;
      fflush(stdout);
    }
    return 0;
  }
  if (argc >= 3 && !strcmp(argv[1], "labelled")) {
    int n = atoi(argv[2]); int npairs = n * (n - 1) / 2;
    static ll I2[1 << MAXN], I3[1 << MAXN];
    size_t FULL = ((size_t)1 << n) - 1;
    for (unsigned long long c = 0; c < (1ULL << npairs); c++) {
      uint32_t out[8]; decode(n, c, out);
      TN = n; tc_alloc(n); memcpy(TOUT, out, n * sizeof(uint32_t));
      tc_cycles(); tc_packing(2, I2); tc_packing(3, I3);
      ll E2, O2, E3, O3; tc_port(0, I2, &E2, &O2); tc_port(0, I3, &E3, &O3);
      set_add(0, I2[FULL], E2 + O2, E2 - O2, I2[FULL & ~1]);
      set_add(1, I3[FULL], E3 + O3, E3 - O3, I3[FULL & ~1]);
    }
    for (int w = 0; w < 2; w++)
      for (size_t i = 0; i < SSZ; i++) if (ST[w] && ST[w][i].used)
        printf("%d %lld %lld %lld %lld\n", w + 2, ST[w][i].k[0], ST[w][i].k[1], ST[w][i].k[2], ST[w][i].k[3]);
    fprintf(stderr, "order %d: %llu labelled tournaments, distinct port states x=2: %zu, x=3: %zu\n", n, 1ULL << npairs, SCNT[0], SCNT[1]);
    return 0;
  }
  fprintf(stderr, "usage: cores classes NMAX | cores labelled N\n");
  return 1;
}
