/* S19 audit: values of fixed-port tails over a list of core states, written independently of fp_cover.c.
 * Usage: cover CORES X LIM [OUT_HITS]
 *   CORES: lines "h0 R0 D0 hu" (distinct core states, from dedup_cores.py)
 *   Letters: S (h,R,D) -> (h,2R,0);  B (h,R,D) -> (h + X(R+D)/2, 2R - hu, -hu).  Only S and B are used.
 * Exactness of the pruning:
 *   - h never decreases (B adds X*E >= 0), so a state with h > LIM is useless;
 *   - after k >= 1 letters S the next B adds X 2^(k-1) R >= X R, so S is useful only if h + X R <= LIM;
 *   - words ending in S have the value of the word without its trailing S's, so every value is the h of a visited state.
 * States are deduplicated per core on (h, R, D) (D is one of D0, 0, -hu, encoded as a 2-bit tag).
 * Prints the misses (odd v for X = 2; v = 1 mod 3 for X = 3; v = 1 mod X in general), the number of states,
 * and optionally writes the hit bitmap.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef long long ll;

static uint64_t *HT; static size_t HSZ, HCNT;
static void ht_init(void) { if (!HT) { HSZ = 1 << 16; HT = calloc(HSZ, 8); } else memset(HT, 0, HSZ * 8); HCNT = 0; }
static uint64_t hmix(uint64_t k) { k ^= k >> 31; k *= 0x7fb5d329728ea185ULL; k ^= k >> 27; k *= 0x81dadef4bc2dd44dULL; k ^= k >> 33; return k; }
static int ht_add(uint64_t key) { /* key != 0 ; returns 1 if new */
  if ((HCNT + 1) * 10 > HSZ * 7) {
    size_t osz = HSZ; uint64_t *old = HT; HSZ *= 2; HT = calloc(HSZ, 8);
    for (size_t i = 0; i < osz; i++) if (old[i]) { size_t j = hmix(old[i]) & (HSZ - 1); while (HT[j]) j = (j + 1) & (HSZ - 1); HT[j] = old[i]; }
    free(old);
  }
  size_t j = hmix(key) & (HSZ - 1);
  while (HT[j]) { if (HT[j] == key) return 0; j = (j + 1) & (HSZ - 1); }
  HT[j] = key; HCNT++; return 1;
}
typedef struct { ll h, R, D; } state;
static state *STK; static size_t SCAP, SLEN;
static void push(ll h, ll R, ll D) { if (SLEN == SCAP) { SCAP = SCAP ? SCAP * 2 : 1 << 16; STK = realloc(STK, SCAP * sizeof(state)); } STK[SLEN].h = h; STK[SLEN].R = R; STK[SLEN].D = D; SLEN++; }

int main(int argc, char **argv) {
  if (argc < 4) { fprintf(stderr, "usage: cover CORES X LIM [OUT_HITS]\n"); return 1; }
  FILE *f = fopen(argv[1], "r"); ll X = atoll(argv[2]), LIM = atoll(argv[3]);
  unsigned char *hit = calloc(LIM + 1, 1);
  ll h0, R0, D0, hu; ll ncores = 0, nstates = 0;
  while (fscanf(f, "%lld %lld %lld %lld", &h0, &R0, &D0, &hu) == 4) {
    ncores++;
    if (h0 > LIM) continue;
    ht_init(); SLEN = 0; push(h0, R0, D0);
    while (SLEN) {
      state s = STK[--SLEN];
      uint64_t tag = (s.D == -hu) ? 1 : (s.D == 0 ? 2 : 3);
      if (s.R < 0 || s.R >= (1LL << 40)) { fprintf(stderr, "R out of range\n"); return 2; }
      uint64_t key = ((uint64_t)s.h << 42) | ((uint64_t)s.R << 2) | tag;
      if (!ht_add(key)) continue;
      nstates++;
      hit[s.h] = 1;
      if ((s.R + s.D) % 2 != 0 && X % 2 != 0) { fprintf(stderr, "odd 2E?\n"); return 3; }
      ll hb = s.h + X * (s.R + s.D) / 2;              /* B */
      if (hb <= LIM) push(hb, 2 * s.R - hu, -hu);
      if (s.h + X * s.R <= LIM) push(s.h, 2 * s.R, 0); /* S */
    }
  }
  ll nmiss = 0;
  printf("x=%lld LIM=%lld cores=%lld states=%lld misses:", X, LIM, ncores, nstates);
  for (ll v = 1; v <= LIM; v++) {
    if (v % X != 1 % X) continue;           /* admissible residue: I(Omega,x) = 1 mod x */
    if (X == 2 && v % 2 == 0) continue;
    if (!hit[v]) { nmiss++; if (nmiss <= 400) printf(" %lld", v); }
  }
  printf("\ncount %lld\n", nmiss);
  /* sanity: no hit outside the admissible residue class */
  ll bad = 0; for (ll v = 1; v <= LIM; v++) if (hit[v] && (v % X != 1 % X)) bad++;
  printf("hits outside residue class 1 mod x: %lld\n", bad);
  if (argc > 4) { FILE *g = fopen(argv[4], "wb"); fwrite(hit, 1, LIM + 1, g); fclose(g); }
  return 0;
}
