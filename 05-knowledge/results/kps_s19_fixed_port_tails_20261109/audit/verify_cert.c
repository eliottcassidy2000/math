/* S19 audit, Q4: re-verify the LP dual certificate with per-core STATE CLOSURE (not the closed form used by q4.c).
 * Usage: verify_cert CORES X LIM CERT
 *   CERT lines: "target num/den" (y_t > 0), first line a comment.
 * For every core state c: run the S/B closure to LIM (as cover.c), and add y_t for each certificate target it reaches.
 * Feasibility of the dual (every core column sum <= 1) proves: every cover of the certificate targets, hence of the
 * whole interval, uses at least ceil(sum y_t) core states.  All arithmetic is in integers (units of 1/10^9).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef long long ll;
static uint64_t *HT; static size_t HSZ, HCNT;
static uint64_t hmix(uint64_t k) { k ^= k >> 31; k *= 0x7fb5d329728ea185ULL; k ^= k >> 27; k *= 0x81dadef4bc2dd44dULL; k ^= k >> 33; return k; }
static void ht_init(void) { if (!HT) { HSZ = 1 << 16; HT = calloc(HSZ, 8); } else memset(HT, 0, HSZ * 8); HCNT = 0; }
static int ht_add(uint64_t key) {
  if ((HCNT + 1) * 10 > HSZ * 7) { size_t o = HSZ; uint64_t *old = HT; HSZ *= 2; HT = calloc(HSZ, 8);
    for (size_t i = 0; i < o; i++) if (old[i]) { size_t j = hmix(old[i]) & (HSZ - 1); while (HT[j]) j = (j + 1) & (HSZ - 1); HT[j] = old[i]; } free(old); }
  size_t j = hmix(key) & (HSZ - 1); while (HT[j]) { if (HT[j] == key) return 0; j = (j + 1) & (HSZ - 1); } HT[j] = key; HCNT++; return 1; }
typedef struct { ll h, R, D; } st; static st *S; static size_t SC, SL;
static void push(ll h, ll R, ll D) { if (SL == SC) { SC = SC ? 2 * SC : 1 << 16; S = realloc(S, SC * sizeof(st)); } S[SL].h = h; S[SL].R = R; S[SL].D = D; SL++; }
int main(int argc, char **argv) {
  FILE *f = fopen(argv[1], "r"); ll X = atoll(argv[2]), LIM = atoll(argv[3]); FILE *g = fopen(argv[4], "r");
  ll *yv = calloc(LIM + 1, sizeof(ll)); char line[512]; ll ntar = 0, total = 0;
  while (fgets(line, sizeof line, g)) {
    if (line[0] == '#') continue; ll t, num, den;
    if (sscanf(line, "%lld %lld/%lld", &t, &num, &den) != 3) { if (sscanf(line, "%lld %lld", &t, &num) == 2) den = 1; else continue; }
    if (1000000000LL % den) { fprintf(stderr, "unexpected denominator %lld\n", den); return 2; }
    yv[t] = num * (1000000000LL / den); total += yv[t]; ntar++;
  }
  unsigned char *seen = calloc(LIM + 1, 1);
  ll h0, R0, D0, hu, nc = 0, maxcol = 0, cmax = -1;
  while (fscanf(f, "%lld %lld %lld %lld", &h0, &R0, &D0, &hu) == 4) {
    memset(seen, 0, LIM + 1); ht_init(); SL = 0;
    if (h0 <= LIM) push(h0, R0, D0);
    while (SL) { st s = S[--SL]; uint64_t tag = s.D == -hu ? 1 : (s.D == 0 ? 2 : 3);
      if (!ht_add(((uint64_t)s.h << 42) | ((uint64_t)s.R << 2) | tag)) continue;
      seen[s.h] = 1;
      ll hb = s.h + X * (s.R + s.D) / 2; if (hb <= LIM) push(hb, 2 * s.R - hu, -hu);
      if (s.h + X * s.R <= LIM) push(s.h, 2 * s.R, 0); }
    ll col = 0; for (ll v = 0; v <= LIM; v++) if (seen[v] && yv[v]) col += yv[v];
    if (col > maxcol) { maxcol = col; cmax = nc; }
    nc++;
  }
  printf("certificate targets %lld, sum y = %lld/1e9; cores %lld; max core column sum = %lld/1e9 (core #%lld)\n", ntar, total, nc, maxcol, cmax);
  if (maxcol <= 1000000000LL) {
    ll lb = (total + 1000000000LL - 1) / 1000000000LL;
    printf("dual feasible: every cover needs at least %lld core states\n", lb);
  } else printf("dual NOT feasible as given\n");
  return 0;
}
