/* AUDIT (independent): orbit of (1,0,1) under nonnegative 3x3 port-chain matrices, values h <= LIM.
 * Usage: a4_cover MATFILE X LIM CAPMODE
 *   MATFILE: one matrix per line, 9 integers (rows h', E', O').
 *   CAPMODE: 0 -> E, O saturated at floor(LIM/X)+1 (the packet's cap); c>0 -> saturated at c (any c >= the
 *            packet's cap must give the same reached h-set if saturation is exact); -1 -> no saturation at all
 *            (exact states; aborts if a coordinate would not fit).
 * Every reached h-value is then CERTIFIED: the word of matrices leading to its first state is rebuilt from parent
 * pointers and multiplied out exactly (128-bit, no saturation); the exact first coordinate must equal h.
 * Prints the complete list of unreached v <= LIM with v = 1 (mod X) (odd v for X = 2), state count, certificate count.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef long long ll;
typedef __int128 i128;
static ll Mt[400][9];
static int nm;
static uint64_t *keys, *slot;          /* queue of packed states; hash slots */
static uint32_t *par;                  /* parent index */
static uint16_t *via;                  /* matrix index */
static uint64_t qn, qcap, hsize, hcount;
static int shift;

static uint64_t hidx(uint64_t k) { return (k * 0x9E3779B97F4A7C15ULL) >> shift; }
static void grow(void) {
  uint64_t old = hsize, *os = slot;
  hsize = old ? old * 2 : (1ULL << 22);
  shift = 64 - __builtin_ctzll(hsize);
  slot = calloc(hsize, 8);
  if (!slot) { fprintf(stderr, "oom\n"); exit(2); }
  for (uint64_t i = 0; i < old; i++) if (os[i]) {
    uint64_t j = hidx(os[i]);
    while (slot[j]) j = (j + 1) & (hsize - 1);
    slot[j] = os[i];
  }
  free(os);
}
static int add(uint64_t k) {
  if (4 * (hcount + 1) > 3 * hsize) grow();
  uint64_t j = hidx(k);
  while (slot[j]) { if (slot[j] == k) return 0; j = (j + 1) & (hsize - 1); }
  slot[j] = k; hcount++;
  return 1;
}

int main(int argc, char **argv) {
  if (argc < 5) { fprintf(stderr, "usage\n"); return 1; }
  FILE *f = fopen(argv[1], "r");
  ll X = atoll(argv[2]), LIM = atoll(argv[3]), capmode = atoll(argv[4]);
  while (nm < 400 && fscanf(f, "%lld %lld %lld %lld %lld %lld %lld %lld %lld", &Mt[nm][0], &Mt[nm][1], &Mt[nm][2],
         &Mt[nm][3], &Mt[nm][4], &Mt[nm][5], &Mt[nm][6], &Mt[nm][7], &Mt[nm][8]) == 9) nm++;
  fclose(f);
  for (int k = 0; k < nm; k++) for (int t = 0; t < 9; t++) if (Mt[k][t] < 0) { fprintf(stderr, "negative entry\n"); return 1; }
  ll CAP = capmode == 0 ? LIM / X + 1 : (capmode > 0 ? capmode : 0);
  int B = 1; while ((1LL << B) <= LIM) B++;
  int BE = capmode < 0 ? (63 - B) / 2 : 1;
  if (capmode >= 0) { while ((1LL << BE) <= CAP) BE++; }
  if (B + 2 * BE > 63) { fprintf(stderr, "packing too wide\n"); return 1; }
  ll EMAX = (1LL << BE) - 1;
  fprintf(stderr, "%d matrices X=%lld LIM=%lld CAP=%lld (mode %lld) bits h=%d e/o=%d\n", nm, X, LIM, CAP, capmode, B, BE);
  qcap = 1 << 20; qn = 0;
  keys = malloc(qcap * 8); par = malloc(qcap * 4); via = malloc(qcap * 2);
  grow();
  uint32_t *first = malloc((LIM + 1) * sizeof(uint32_t));
  memset(first, 0xff, (LIM + 1) * sizeof(uint32_t));
#define PK(h, e, o) (((uint64_t)(h) << (2 * BE)) | ((uint64_t)(e) << BE) | (uint64_t)(o))
  uint64_t s0 = PK(1, 0, 1);
  add(s0); keys[0] = s0; par[0] = 0; via[0] = 0xffff; qn = 1; first[1] = 0;
  for (uint64_t qi = 0; qi < qn; qi++) {
    uint64_t s = keys[qi];
    ll h = (ll)(s >> (2 * BE)), e = (ll)((s >> BE) & EMAX), o = (ll)(s & EMAX);
    for (int k = 0; k < nm; k++) {
      ll *m = Mt[k];
      ll h2 = m[0] * h + m[1] * e + m[2] * o;
      if (h2 > LIM) continue;
      ll e2 = m[3] * h + m[4] * e + m[5] * o, o2 = m[6] * h + m[7] * e + m[8] * o;
      if (capmode >= 0) { if (e2 > CAP) e2 = CAP; if (o2 > CAP) o2 = CAP; }
      else if (e2 > EMAX || o2 > EMAX) { fprintf(stderr, "exact mode overflow\n"); return 3; }
      uint64_t t = PK(h2, e2, o2);
      if (add(t)) {
        if (qn == qcap) { qcap *= 2; keys = realloc(keys, qcap * 8); par = realloc(par, qcap * 4); via = realloc(via, qcap * 2);
          if (!keys || !par || !via) { fprintf(stderr, "oom\n"); return 2; } }
        keys[qn] = t; par[qn] = (uint32_t)qi; via[qn] = (uint16_t)k;
        if (first[h2] == 0xffffffffu) first[h2] = (uint32_t)qn;
        qn++;
      }
    }
  }
  /* certificates */
  ll ncert = 0, nbad = 0, maxlen = 0;
  int *word = malloc(sizeof(int) * 100000);
  for (ll v = 1; v <= LIM; v++) {
    if (first[v] == 0xffffffffu) continue;
    int len = 0; uint32_t i = first[v];
    while (i != 0) { word[len++] = via[i]; i = par[i]; }
    i128 h = 1, e = 0, o = 1;
    for (int t = len - 1; t >= 0; t--) {
      ll *m = Mt[word[t]];
      i128 h2 = m[0] * h + m[1] * e + m[2] * o, e2 = m[3] * h + m[4] * e + m[5] * o, o2 = m[6] * h + m[7] * e + m[8] * o;
      h = h2; e = e2; o = o2;
    }
    if (h == v) ncert++; else nbad++;
    if (len > maxlen) maxlen = len;
  }
  ll nmiss = 0;
  printf("misses (v = 1 mod %lld, v <= %lld):", X, LIM);
  for (ll v = 1; v <= LIM; v++) {
    if (X == 2 ? (v % 2 == 0) : (v % X != 1)) continue;
    if (first[v] == 0xffffffffu) { nmiss++; printf(" %lld", v); }
  }
  printf("\ncount %lld; states %llu; reached h-values certified by exact word products: %lld, failures %lld;"
         " longest certificate word %lld\n", nmiss, (unsigned long long)qn, ncert, nbad, maxlen);
  return 0;
}
