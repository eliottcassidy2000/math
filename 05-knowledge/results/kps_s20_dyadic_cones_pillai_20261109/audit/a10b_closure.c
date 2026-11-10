/* a10b_closure.c (S20 audit; a10_closure.c with a three-word key for h, e > 2^36): bounded odd-part closures C(r0) of HYP-9656, written independently.
   Start (1,1,1); moves: every matrix of the library (exit inside the gadget) and sink padding (h,e,R) -> (h,0,2R);
   keep a triple iff oddpart(R) <= r0.
   mode 0 (exact): store every reachable triple with h <= HMAX (any R) and record its value h; expand a stored triple
                   only if R <= HMAX.  Justification: every library matrix has E_u >= 1 and 0 <= s_u <= h_G (checked in
                   a09), so on realisable triples a gadget child has h' >= E_u R >= R; a sink child keeps h.  Hence
                   triples with R > HMAX cannot produce new values <= HMAX, and the value set below HMAX is exact.
   mode 1 (author emulation): drop triples with R > 2*HMAX + 2 (as rclosure3/4.c do); expand everything stored.
   usage: a10_closure MATFILE HMAX r0 mode log2hash [valuesfile]
   prints: triples, values, max value, values per decade, complete prefix (least odd n not 7, 21 missing). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll; typedef unsigned long long u64; typedef __int128 i128;
static u64 *K1, *K2, *K3; static u64 HS, MASK, ntrip = 0;
static ll *stk; static size_t sp = 0, scap = 0;
static ll HMAX, r0; static int mode;
static ll *vals; static size_t nvals = 0, vcap = 0;
static inline u64 hsh(u64 a, u64 b){ u64 x = a * 0x9FB21C651E98DF25ULL + b * 0xD6E8FEB86659FD93ULL; x ^= x >> 29; x *= 0xBF58476D1CE4E5B9ULL; x ^= x >> 32; return x; }
static void push(ll h, ll e, ll R){
  if (sp + 3 > scap) { scap = scap ? 2 * scap : (1 << 24); stk = realloc(stk, scap * sizeof(ll)); if (!stk) { fprintf(stderr, "oom\n"); exit(3); } }
  stk[sp++] = h; stk[sp++] = e; stk[sp++] = R;
}
static void add(i128 h, i128 e, i128 R){
  if (h < 1 || h > HMAX) return;
  if (R < 1 || e < 0 || e > h) { fprintf(stderr, "non-realisable triple produced\n"); exit(5); }
  if (mode == 1 && R > (i128)(2 * HMAX + 2)) return;
  if (R >= ((i128)1 << 62)) { fprintf(stderr, "R overflow\n"); exit(6); }
  u64 r = (u64)R; int k = __builtin_ctzll(r); u64 o = r >> k;
  if ((ll)o > r0) return;
  u64 a = (u64)h, b = (u64)e, c3 = (u64)o | ((u64)k << 32);
  if ((u64)e >= (1ULL << 36) || o >= (1ULL << 16)) { fprintf(stderr, "pack overflow\n"); exit(7); }
  u64 i = hsh(a, b ^ (c3 * 0x94D049BB133111EBULL)) & MASK;
  while (K1[i]) { if (K1[i] == a && K2[i] == b && K3[i] == c3) return; i = (i + 1) & MASK; }
  K1[i] = a; K2[i] = b; K3[i] = c3; ntrip++;
  if (ntrip > HS / 10 * 8) { fprintf(stderr, "hash too full\n"); exit(2); }
  if (nvals == vcap) { vcap = vcap ? 2 * vcap : (1 << 20); vals = realloc(vals, vcap * sizeof(ll)); }
  vals[nvals++] = (ll)h;
  if (mode == 1 || R <= (i128)HMAX) push((ll)h, (ll)e, (ll)R);
}
static int cmpll(const void *x, const void *y){ ll a = *(const ll *)x, b = *(const ll *)y; return a < b ? -1 : a > b; }
int main(int argc, char **argv){
  if (argc < 6) { fprintf(stderr, "usage\n"); return 1; }
  FILE *f = fopen(argv[1], "r"); HMAX = atoll(argv[2]); r0 = atoll(argv[3]); mode = atoi(argv[4]); int lg = atoi(argv[5]);
  ll (*M)[9] = malloc(sizeof(ll[9]) * 50000); int nm = 0;
  while (1) { ll m[9]; int ok = 1; for (int i = 0; i < 9; i++) if (fscanf(f, "%lld", &m[i]) != 1) { ok = 0; break; } if (!ok) break; memcpy(M[nm++], m, sizeof m); }
  HS = 1ULL << lg; MASK = HS - 1; K1 = calloc(HS, sizeof(u64)); K2 = calloc(HS, sizeof(u64)); K3 = calloc(HS, sizeof(u64));
  if (!K1 || !K2) { fprintf(stderr, "oom hash\n"); return 3; }
  add(1, 1, 1);
  while (sp) {
    ll R = stk[--sp], e = stk[--sp], h = stk[--sp];
    add(h, 0, (i128)2 * R);
    for (int i = 0; i < nm; i++) {
      ll *q = M[i];
      add((i128)q[0] * h + (i128)q[1] * e + (i128)q[2] * R, (i128)q[3] * h + (i128)q[4] * e + (i128)q[5] * R,
          (i128)q[6] * h + (i128)q[7] * e + (i128)q[8] * R);
    }
  }
  qsort(vals, nvals, sizeof(ll), cmpll);
  size_t d = 0; for (size_t i = 0; i < nvals; i++) if (i == 0 || vals[i] != vals[i - 1]) vals[d++] = vals[i];
  printf("matrices=%d HMAX=%lld r0=%lld mode=%s triples=%llu values=%zu max=%lld\n", nm, HMAX, r0, mode ? "author-Rcap" : "exact",
         (unsigned long long)ntrip, d, d ? vals[d - 1] : 0);
  printf("  per decade:");
  for (ll lim = 10; lim <= HMAX && lim > 0; lim *= 10) { size_t c = 0; for (size_t i = 0; i < d; i++) if (vals[i] <= lim) c++; printf(" <=%lld:%zu", lim, c); }
  printf("\n");
  ll n = 1; size_t idx = 0;
  for (;; n += 2) { while (idx < d && vals[idx] < n) idx++; if (n == 7 || n == 21) continue; if (idx < d && vals[idx] == n) continue; break; }
  printf("  complete prefix: every odd n < %lld except 7, 21 is a value\n", n);
  if (argc > 6) { FILE *g = fopen(argv[6], "w"); for (size_t i = 0; i < d; i++) fprintf(g, "%lld\n", vals[i]); fclose(g); }
  return 0;
}
