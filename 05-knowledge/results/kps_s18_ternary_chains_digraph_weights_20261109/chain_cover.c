/* Port-chain coverage (THM-5070): breadth-first closure of the start state (h,E,O) = (1,0,1) under the
 * nonnegative block matrices written by chain_blocks.py, keeping states with h <= LIM.
 * E and O are saturated at CAP = LIM/X + 1: every block adds X*(q_o E + q_e O) to h with integer q's,
 * so a coordinate >= CAP with a positive coefficient forces h > LIM; nonnegative matrices make the
 * saturation exact (min(.,CAP) is a semiring homomorphism on the reachable values).
 * Usage: chain_cover MATFILE X LIM [MODCLASS]
 * Prints the values v <= LIM with v = 1 (mod X) (X = 2: odd v) that are not reached.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef long long ll;
static ll (*M)[9];
static int nm;
static uint64_t *tab;
static uint64_t tabsize, tabcount;
static uint64_t *queue;
static uint64_t qcap, qlen;

static uint64_t hash64(uint64_t x) {
  x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL; x ^= x >> 33;
  return x;
}
static void rehash(void);
static int insert(uint64_t key) { /* key != 0 */
  if (2 * (tabcount + 1) > tabsize) rehash();
  uint64_t i = hash64(key) & (tabsize - 1);
  while (tab[i]) {
    if (tab[i] == key) return 0;
    i = (i + 1) & (tabsize - 1);
  }
  tab[i] = key; tabcount++;
  return 1;
}
static void rehash(void) {
  uint64_t old = tabsize, *ot = tab;
  tabsize = old ? old * 2 : (1ULL << 20);
  tab = calloc(tabsize, sizeof(uint64_t));
  if (!tab) { fprintf(stderr, "out of memory\n"); exit(1); }
  tabcount = 0;
  for (uint64_t i = 0; i < old; i++) if (ot[i]) {
    uint64_t j = hash64(ot[i]) & (tabsize - 1);
    while (tab[j]) j = (j + 1) & (tabsize - 1);
    tab[j] = ot[i]; tabcount++;
  }
  free(ot);
}
static void push(uint64_t key) {
  if (qlen == qcap) { qcap = qcap ? qcap * 2 : (1ULL << 20); queue = realloc(queue, qcap * sizeof(uint64_t));
    if (!queue) { fprintf(stderr, "out of memory\n"); exit(1); } }
  queue[qlen++] = key;
}

int main(int argc, char **argv) {
  if (argc < 4) { fprintf(stderr, "usage: chain_cover MATFILE X LIM\n"); return 1; }
  FILE *f = fopen(argv[1], "r");
  ll X = atoll(argv[2]), LIM = atoll(argv[3]);
  ll CAP = LIM / X + 1;
  int B = 0; while ((1LL << B) <= (LIM > CAP ? LIM : CAP)) B++;
  if (3 * B > 63) { fprintf(stderr, "LIM too large for packing\n"); return 1; }
  M = malloc(sizeof(*M) * 100000);
  while (nm < 100000 && fscanf(f, "%lld %lld %lld %lld %lld %lld %lld %lld %lld",
         &M[nm][0], &M[nm][1], &M[nm][2], &M[nm][3], &M[nm][4], &M[nm][5], &M[nm][6], &M[nm][7], &M[nm][8]) == 9) nm++;
  fclose(f);
  fprintf(stderr, "%d matrices, X=%lld LIM=%lld CAP=%lld bits=%d\n", nm, X, LIM, CAP, B);
  unsigned char *hit = calloc(LIM + 1, 1);
  uint64_t mask = (1ULL << B) - 1;
#define PACK(h,e,o) ((((uint64_t)(h)) << (2*B)) | (((uint64_t)(e)) << B) | ((uint64_t)(o)))
  rehash();
  uint64_t s0 = PACK(1, 0, 1);
  insert(s0); push(s0); hit[1] = 1;
  for (uint64_t qi = 0; qi < qlen; qi++) {
    uint64_t s = queue[qi];
    ll h = (ll)(s >> (2 * B)), e = (ll)((s >> B) & mask), o = (ll)(s & mask);
    for (int k = 0; k < nm; k++) {
      ll *m = M[k];
      ll h2 = m[0] * h + m[1] * e + m[2] * o;
      if (h2 > LIM) continue;
      ll e2 = m[3] * h + m[4] * e + m[5] * o;
      ll o2 = m[6] * h + m[7] * e + m[8] * o;
      if (e2 > CAP) e2 = CAP;
      if (o2 > CAP) o2 = CAP;
      uint64_t t = PACK(h2, e2, o2);
      if (insert(t)) { push(t); hit[h2] = 1; }
    }
    if ((qi & ((1ULL << 24) - 1)) == 0 && qi) fprintf(stderr, "processed %llu of %llu states\n",
        (unsigned long long)qi, (unsigned long long)qlen);
  }
  fprintf(stderr, "states %llu\n", (unsigned long long)qlen);
  ll nmiss = 0;
  printf("misses (v = 1 mod %lld, v <= %lld):", X, LIM);
  for (ll v = 1; v <= LIM; v++) {
    if (v % X != 1 % X) continue;
    if (X == 2 && v % 2 == 0) continue;
    if (!hit[v]) { nmiss++; if (nmiss <= 200) printf(" %lld", v); }
  }
  printf("\ncount %lld states %llu\n", nmiss, (unsigned long long)qlen);
  return 0;
}
