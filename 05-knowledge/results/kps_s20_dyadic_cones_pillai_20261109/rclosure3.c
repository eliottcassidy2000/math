/* rclosure3: closure of triples (h,e,R) with oddpart(R) <= r0 under one-arc gadget matrices (exit inside the
   gadget) and sink padding (h,e,R) -> (h,0,2R).  Along every step h does not decrease (THM-4875 (C) on realisable
   triples), so the closure below a bound HMAX is exact.  Wide (64-bit) storage of h, e, R.
   usage: rclosure3 GADGETFILE GMAX HMAX r0 [log2 hash size]
   output: number of distinct matrices, triples, distinct values, values by decade, all values when <= 600. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll; typedef unsigned long long ull;
static ll *kh, *ke, *kR; static char *ku; static long HS, ntrip = 0;
static ll *sh, *se, *sR; static long sp = 0, scap;
static ll HMAX, RMAX, r0;
static ull mix(ull x){ x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL; x ^= x >> 33; return x; }
static void add(ll h, ll e, ll R){
  if (h < 1 || h > HMAX || R < 1 || R > RMAX || e < 0 || e > h) return;
  ll o = R; while (!(o & 1)) o >>= 1; if (o > r0) return;
  ull i = mix((ull)h * 0x9E3779B97F4A7C15ULL ^ (ull)e * 0xC2B2AE3D27D4EB4FULL ^ (ull)R) & (HS - 1);
  while (ku[i]) { if (kh[i] == h && ke[i] == e && kR[i] == R) return; i = (i + 1) & (HS - 1); }
  ku[i] = 1; kh[i] = h; ke[i] = e; kR[i] = R; ntrip++;
  if (ntrip > HS * 0.7) { fprintf(stderr, "hash full\n"); exit(2); }
  if (sp >= scap) { scap *= 2; sh = realloc(sh, scap * sizeof(ll)); se = realloc(se, scap * sizeof(ll)); sR = realloc(sR, scap * sizeof(ll));
    if (!sh || !se || !sR) { fprintf(stderr, "oom\n"); exit(3); } }
  sh[sp] = h; se[sp] = e; sR[sp] = R; sp++;
}
static int cmpll(const void *a, const void *b){ ll x = *(const ll *)a, y = *(const ll *)b; return x < y ? -1 : x > y; }
int main(int argc, char **argv){
  FILE *f = fopen(argv[1], "r"); int GMAX = atoi(argv[2]); HMAX = atoll(argv[3]); r0 = atoll(argv[4]); RMAX = 2 * HMAX + 2;
  HS = 1L << atoi(argc > 5 ? argv[5] : "22");
  ll (*M)[9] = malloc(sizeof(ll[9]) * 40000); int nm = 0; int g, u, z; char codebuf[64]; ll m[9];
  while (fscanf(f, "%d %63s %d %d", &g, codebuf, &u, &z) == 4) {
    for (int i = 0; i < 9; i++) if (fscanf(f, "%lld", &m[i]) != 1) return 1;
    if (g > GMAX) continue;
    int dup = 0; for (int i = 0; i < nm; i++) if (!memcmp(M[i], m, sizeof m)) { dup = 1; break; }
    if (!dup) { memcpy(M[nm], m, sizeof m); nm++; }
  }
  kh = calloc(HS, sizeof(ll)); ke = calloc(HS, sizeof(ll)); kR = calloc(HS, sizeof(ll)); ku = calloc(HS, 1);
  scap = 1 << 20; sh = malloc(scap * sizeof(ll)); se = malloc(scap * sizeof(ll)); sR = malloc(scap * sizeof(ll));
  add(1, 1, 1);
  while (sp) {
    sp--; ll h = sh[sp], e = se[sp], R = sR[sp];
    add(h, 0, 2 * R);
    for (int i = 0; i < nm; i++) { ll *q = M[i]; add(q[0]*h + q[1]*e + q[2]*R, q[3]*h + q[4]*e + q[5]*R, q[6]*h + q[7]*e + q[8]*R); }
  }
  ll *vals = malloc(sizeof(ll) * (ntrip + 1)); long nv = 0;
  for (long i = 0; i < HS; i++) if (ku[i]) vals[nv++] = kh[i];
  qsort(vals, nv, sizeof(ll), cmpll); long d = 0;
  for (long i = 0; i < nv; i++) if (i == 0 || vals[i] != vals[i - 1]) vals[d++] = vals[i];
  printf("GMAX=%d r0=%lld HMAX=%lld distinct-matrices=%d triples=%ld values=%ld max=%lld\n", GMAX, r0, HMAX, nm, ntrip, d, d ? vals[d - 1] : 0);
  for (ll lim = 10; lim <= HMAX; lim *= 10) { long c = 0; for (long i = 0; i < d; i++) if (vals[i] <= lim) c++; printf("  values<=%lld: %ld\n", lim, c); }
  if (d <= 600) { printf("  all:"); for (long i = 0; i < d; i++) printf(" %lld", vals[i]); printf("\n"); }
  return 0;
}
