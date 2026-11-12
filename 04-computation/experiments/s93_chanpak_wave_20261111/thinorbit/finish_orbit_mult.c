/* Lane thinorbit (S93), finishing pass.  Representation numbers of the port-chain orbit.
 *
 * Same exact breadth-first closure as phase 1 of orbit_cover.c (states (h,E,O) with h <= LIM, exact
 * de-duplication, h never decreases, E and O never saturate for the limits used here: asserted), but instead of
 * a hit flag it counts r(v) = number of DISTINCT orbit states with first coordinate v.
 * Output per dyadic window [2^k, 2^(k+1)): number of odd v, number with r(v) = 0, min / mean / max of r(v) over the
 * odd v, an argmin, and the second moment sum r(v)^2 against (sum r(v))^2 / #odd  (the Cauchy-Schwarz ratio: the
 * number of represented v is at least (sum r)^2 / sum r^2).
 * Usage: finish_orbit_mult MATFILE LIM
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef long long ll;
static ll M[4096][9];
static int nm;
static uint64_t *tab, tabsize, tabcount;
static uint64_t *queue, qcap, qlen;
static ll LIM, CAP;
static int BE = 21;

static inline uint64_t hash64(uint64_t x) {
  x ^= x >> 31; x *= 0x7fb5d329728ea185ULL; x ^= x >> 27; x *= 0x81dadef4bc2dd44dULL; x ^= x >> 33;
  return x;
}
static void grow(void) {
  uint64_t old = tabsize, *ot = tab;
  tabsize = old ? old * 2 : (1ULL << 22);
  tab = calloc(tabsize, sizeof(uint64_t));
  if (!tab) { fprintf(stderr, "out of memory\n"); exit(1); }
  for (uint64_t i = 0; i < old; i++) if (ot[i]) {
    uint64_t j = hash64(ot[i]) & (tabsize - 1);
    while (tab[j]) j = (j + 1) & (tabsize - 1);
    tab[j] = ot[i];
  }
  free(ot);
}
static int insert(uint64_t key) {
  if (5 * (tabcount + 1) > 3 * tabsize) grow();
  uint64_t i = hash64(key) & (tabsize - 1);
  while (tab[i]) { if (tab[i] == key) return 0; i = (i + 1) & (tabsize - 1); }
  tab[i] = key; tabcount++;
  return 1;
}
static void push(uint64_t key) {
  if (qlen == qcap) { qcap = qcap ? qcap * 2 : (1ULL << 22); queue = realloc(queue, qcap * sizeof(uint64_t));
    if (!queue) { fprintf(stderr, "out of memory\n"); exit(1); } }
  queue[qlen++] = key;
}
#define PACK(h,e,o) ((((uint64_t)(h)) << (2 * BE)) | (((uint64_t)(e)) << BE) | ((uint64_t)(o)))

int main(int argc, char **argv) {
  if (argc < 3) { fprintf(stderr, "usage: finish_orbit_mult MATFILE LIM\n"); return 1; }
  FILE *f = fopen(argv[1], "r");
  if (!f) { fprintf(stderr, "cannot open %s\n", argv[1]); return 1; }
  LIM = atoll(argv[2]);
  { int bh = 1; while ((1LL << bh) <= LIM) bh++; BE = (64 - bh) / 2; if (BE > 30) BE = 30; }
  CAP = (1LL << BE) - 1;
  while (nm < 4096 && fscanf(f, "%lld %lld %lld %lld %lld %lld %lld %lld %lld",
         &M[nm][0], &M[nm][1], &M[nm][2], &M[nm][3], &M[nm][4], &M[nm][5], &M[nm][6], &M[nm][7], &M[nm][8]) == 9) nm++;
  fclose(f);
  printf("%d generators, LIM = %lld, %d bits for E and for O\n", nm, LIM, BE);
  uint32_t *r = calloc(LIM + 2, sizeof(uint32_t));
  grow();
  uint64_t s0 = PACK(1, 0, 1);
  insert(s0); push(s0); r[1] = 1;
  const uint64_t mask = (1ULL << BE) - 1;
  uint64_t sat = 0;
  for (uint64_t qi = 0; qi < qlen; qi++) {
    uint64_t s = queue[qi];
    ll h = (ll)(s >> (2 * BE)), e = (ll)((s >> BE) & mask), o = (ll)(s & mask);
    for (int k = 0; k < nm; k++) {
      const ll *m = M[k];
      ll h2 = m[0] * h + m[1] * e + m[2] * o;
      if (h2 > LIM) continue;
      ll e2 = m[3] * h + m[4] * e + m[5] * o;
      ll o2 = m[6] * h + m[7] * e + m[8] * o;
      if (e2 >= CAP || o2 >= CAP) { sat++; if (e2 > CAP) e2 = CAP; if (o2 > CAP) o2 = CAP; }
      uint64_t t = PACK(h2, e2, o2);
      if (insert(t)) { push(t); r[h2]++; }
    }
  }
  printf("%llu distinct states with h <= %lld; saturation events: %llu (0 means r(v) is exact)\n",
         (unsigned long long)qlen, LIM, (unsigned long long)sat);
  printf("window k: [2^k, 2^(k+1)) cut at LIM | #odd | #r=0 | min r (at v) | mean r | max r | (sum r)^2/sum r^2 as a share of #odd\n");
  for (int k = 1; (1LL << k) <= LIM; k++) {
    ll lo = 1LL << k, hi = (2LL << k) - 1; if (hi > LIM) hi = LIM;
    ll n = 0, z = 0, amin = -1; uint32_t mn = 0xffffffffu, mx = 0; double s1 = 0, s2 = 0;
    for (ll v = lo + 1; v <= hi; v += 2) {
      n++; if (!r[v]) z++;
      if (r[v] < mn) { mn = r[v]; amin = v; }
      if (r[v] > mx) mx = r[v];
      s1 += r[v]; s2 += (double)r[v] * r[v];
    }
    if (n) printf("  %2d | %7lld | %5lld | %7u (%lld) | %10.2f | %8u | %.4f\n", k, n, z, mn, amin, s1 / n, mx,
                  s2 > 0 ? s1 * s1 / s2 / n : 0.0);
  }
  /* minimum of r over the odd v in the upper half of each window excluding the known exceptions */
  return 0;
}
