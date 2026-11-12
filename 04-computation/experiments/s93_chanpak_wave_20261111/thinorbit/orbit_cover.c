/* Lane thinorbit (S93).  Coverage and growth of the port-chain orbit of (h,E,O) = (1,0,1) under a set of
 * nonnegative 3 x 3 integer matrices (HYP-9640 (B)).  Written independently of chain_cover.c.
 *
 * Phase 1: breadth-first closure with exact de-duplication of all states with h <= LSTORE.  E and O are
 *          saturated at CAP = 2^BE - 1 >= LIM/2 + 1 (BE = number of bits left for E and for O in a 64-bit key).
 *          The first row of every generator is (eta, 2 q_o, 2 q_e), so a coordinate >= LIM/2 + 1 with a positive
 *          coefficient forces h > LIM, and min(., CAP) commutes with nonnegative integer rows: the set of
 *          reached h <= LIM is exact.  The number of saturated stored states is printed (0 means that the
 *          state counts are exact counts of orbit points).
 * Phase 2: from every stored state, depth-first search WITHOUT de-duplication through the states with
 *          LSTORE < h <= LIM, marking hit[h].  h never decreases (eta >= 1), so every state with h <= LIM is
 *          reached from a stored state through states with h > LSTORE only; the search is complete when no
 *          depth limit is given.  A child equal to its parent (the fixed points (h,0,h) of the rank-one
 *          generator (1,0,1)(1,2,0)) is skipped; the depth cap 400 is asserted never to be reached.
 * Output: the odd v <= LIM not reached; the number of stored states with h <= 2^k (growth of the orbit);
 *         the number of phase-2 visits.
 * Usage: orbit_cover MATFILE LIM LSTORE [MAXDEPTH]
 *        MAXDEPTH (optional) limits the number of phase-2 steps; then the miss list is only an upper bound
 *        for the true misses (every reported hit is still a true orbit value).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>

typedef long long ll;
static ll M[4096][9];
static int nm;
static uint64_t *tab, tabsize, tabcount;
static uint64_t *queue, qcap, qlen;
static ll LIM, LSTORE, CAP;
static int BE = 21;
static int MAXDEPTH = 400;
static unsigned char *hit;
static int deep_flag = 0;

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

static uint64_t dfs(ll h, ll e, ll o, int depth) {
  uint64_t visits = 1;
  hit[h] = 1;
  if (depth >= MAXDEPTH) { if (MAXDEPTH >= 400) deep_flag = 1; return visits; }
  for (int k = 0; k < nm; k++) {
    const ll *m = M[k];
    ll h2 = m[0] * h + m[1] * e + m[2] * o;
    if (h2 > LIM) continue;
    ll e2 = m[3] * h + m[4] * e + m[5] * o;
    ll o2 = m[6] * h + m[7] * e + m[8] * o;
    if (e2 > CAP) e2 = CAP;
    if (o2 > CAP) o2 = CAP;
    if (h2 == h && e2 == e && o2 == o) continue;
    visits += dfs(h2, e2, o2, depth + 1);
  }
  return visits;
}

int main(int argc, char **argv) {
  if (argc < 4) { fprintf(stderr, "usage: orbit_cover MATFILE LIM LSTORE [MAXDEPTH]\n"); return 1; }
  FILE *f = fopen(argv[1], "r");
  if (!f) { fprintf(stderr, "cannot open %s\n", argv[1]); return 1; }
  LIM = atoll(argv[2]); LSTORE = atoll(argv[3]);
  if (LSTORE > LIM) LSTORE = LIM;
  if (argc > 4) MAXDEPTH = atoi(argv[4]);
  { int bh = 1; while ((1LL << bh) <= LIM) bh++; BE = (64 - bh) / 2; if (BE > 30) BE = 30; }
  CAP = (1LL << BE) - 1;
  if (CAP < LIM / 2 + 1) { fprintf(stderr, "LIM too large for the packing\n"); return 1; }
  while (nm < 4096 && fscanf(f, "%lld %lld %lld %lld %lld %lld %lld %lld %lld",
         &M[nm][0], &M[nm][1], &M[nm][2], &M[nm][3], &M[nm][4], &M[nm][5], &M[nm][6], &M[nm][7], &M[nm][8]) == 9) nm++;
  fclose(f);
  printf("%d generators, LIM = %lld, LSTORE = %lld, CAP = %lld (%d bits for E and for O), phase-2 depth limit %d\n",
         nm, LIM, LSTORE, CAP, BE, MAXDEPTH);
  hit = calloc(LIM + 2, 1);
  grow();
  uint64_t s0 = PACK(1, 0, 1);
  insert(s0); push(s0); hit[1] = 1;
  const uint64_t mask = (1ULL << BE) - 1;
  for (uint64_t qi = 0; qi < qlen; qi++) {
    uint64_t s = queue[qi];
    ll h = (ll)(s >> (2 * BE)), e = (ll)((s >> BE) & mask), o = (ll)(s & mask);
    for (int k = 0; k < nm; k++) {
      const ll *m = M[k];
      ll h2 = m[0] * h + m[1] * e + m[2] * o;
      if (h2 > LSTORE) continue;
      ll e2 = m[3] * h + m[4] * e + m[5] * o;
      ll o2 = m[6] * h + m[7] * e + m[8] * o;
      if (e2 > CAP) e2 = CAP;
      if (o2 > CAP) o2 = CAP;
      uint64_t t = PACK(h2, e2, o2);
      if (insert(t)) { push(t); hit[h2] = 1; }
    }
  }
  printf("phase 1: %llu distinct states with h <= %lld\n", (unsigned long long)qlen, LSTORE);
  {
    uint64_t cnt[70]; memset(cnt, 0, sizeof(cnt));
    uint64_t sat = 0;
    double maxe = 0, maxo = 0, mino = 1e18;
    for (uint64_t qi = 0; qi < qlen; qi++) {
      ll h = (ll)(queue[qi] >> (2 * BE));
      ll e = (ll)((queue[qi] >> BE) & mask), o = (ll)(queue[qi] & mask);
      int k = 0; while ((1LL << k) < h) k++;
      cnt[k]++;
      if (e >= CAP || o >= CAP) { sat++; continue; }
      if ((double)e / h > maxe) maxe = (double)e / h;
      if ((double)o / h > maxo) maxo = (double)o / h;
      if ((double)o / h < mino) mino = (double)o / h;
    }
    uint64_t cum = 0;
    printf("saturated stored states: %llu; over the others: max E/h = %.6f, max O/h = %.6f, min O/h = %.6f\n",
           (unsigned long long)sat, maxe, maxo, mino);
    printf("growth: k, #distinct states with h <= 2^k (meaningful for 2^k <= LSTORE)\n");
    for (int k = 0; k < 60 && (1LL << k) < 2 * LSTORE; k++) { cum += cnt[k]; printf("  %2d %llu\n", k, (unsigned long long)cum); }
  }
  uint64_t visits = 0;
  if (LSTORE < LIM) {
    double t0 = omp_get_wtime();
    #pragma omp parallel for schedule(dynamic, 2048) reduction(+:visits)
    for (long long qi = 0; qi < (long long)qlen; qi++) {
      uint64_t s = queue[qi];
      ll h = (ll)(s >> (2 * BE)), e = (ll)((s >> BE) & mask), o = (ll)(s & mask);
      for (int k = 0; k < nm; k++) {
        const ll *m = M[k];
        ll h2 = m[0] * h + m[1] * e + m[2] * o;
        if (h2 <= LSTORE || h2 > LIM) continue;
        ll e2 = m[3] * h + m[4] * e + m[5] * o;
        ll o2 = m[6] * h + m[7] * e + m[8] * o;
        if (e2 > CAP) e2 = CAP;
        if (o2 > CAP) o2 = CAP;
        visits += dfs(h2, e2, o2, 1);
      }
    }
    printf("phase 2: %llu visits (states with LSTORE < h <= LIM, with multiplicity), %.1f s, depth cap reached: %d\n",
           (unsigned long long)visits, omp_get_wtime() - t0, deep_flag);
  }
  ll nmiss = 0, lastmiss = -1;
  printf("odd v <= %lld not reached:", LIM);
  for (ll v = 1; v <= LIM; v += 2) if (!hit[v]) { nmiss++; lastmiss = v; if (nmiss <= 130) printf(" %lld", v); }
  if (nmiss > 130) printf(" ...");
  printf("\ncount of misses %lld, largest miss %lld\n", nmiss, lastmiss);
  printf("reached odd values by dyadic window [2^k, 2^(k+1)) (cut at LIM): k reached total\n");
  for (int k = 1; (1LL << k) <= LIM; k++) {
    ll lo = 1LL << k, hi = (2LL << k) - 1; if (hi > LIM) hi = LIM;
    ll c = 0, tot = 0; for (ll v = lo + 1; v <= hi; v += 2) { c += hit[v]; tot++; }
    printf("  %2d %lld %lld\n", k, c, tot);
  }
  return 0;
}
