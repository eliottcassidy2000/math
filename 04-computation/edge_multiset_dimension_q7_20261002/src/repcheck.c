// repcheck.c -- independent brute-force check of an orbit-representative file (a-subsets of Q_n, n <= 6).
// For every listed set T: |T| = a, and T is the lexicographically smallest sorted tuple among its images
// under the automorphisms x -> pi(x xor t), computed from explicit vertex tables (no SJT, no shortcuts,
// no code shared with orbreps.c).  The sets must be pairwise distinct.  Lex-min sorted tuples are a
// canonical form, so a list of distinct canonical sets is pairwise inequivalent; if moreover its length
// equals the Burnside count, it contains exactly one representative of every orbit.
//   mode full : all 2^n n! group elements are tried.
//   mode fast : (default) first require 0 in T (otherwise translating an element of T to 0 gives a sorted
//               tuple starting with 0 < min T); then only the images that contain 0 can be lex-smaller than
//               T, and pi(x xor t) = 0 iff x = t, so only the translations t in T (with all n! pi) are tried.
// Lex order of equal-size sets X, Y as masks: X <_lex Y iff the lowest element of X xor Y lies in X.
// usage: repcheck n a repfile [first last] [full]
//        prints: REPCHECK n= a= mode= file_reps= checked= noncanonical= wrongsize= duplicates=
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
static int n, N, nperm, img[720][64];
static int cmpu64(const void *x, const void *y) { uint64_t a = *(const uint64_t *)x, b = *(const uint64_t *)y; return a < b ? -1 : a > b; }
int main(int argc, char **argv) {
  if (argc < 4) { fprintf(stderr, "usage: repcheck n a repfile [first last] [full]\n"); return 1; }
  n = atoi(argv[1]); int a = atoi(argv[2]); N = 1 << n;
  if (n < 1 || n > 6) return 1;
  int full = !strcmp(argv[argc - 1], "full");
  int idx[6] = {0}; nperm = 0;
  for (;;) {   // all maps i -> idx[i]; keep the bijections
    int used = 0, ok = 1;
    for (int i = 0; i < n; i++) { if (used >> idx[i] & 1) ok = 0; used |= 1 << idx[i]; }
    if (ok) { for (int v = 0; v < N; v++) { int u = 0; for (int i = 0; i < n; i++) if (v >> i & 1) u |= 1 << idx[i]; img[nperm][v] = u; } nperm++; }
    int i = 0; while (i < n && ++idx[i] == n) { idx[i] = 0; i++; } if (i == n) break;
  }
  { long f = 1; for (int i = 2; i <= n; i++) f *= i; if (nperm != f) { printf("FATAL nperm\n"); return 9; } }
  FILE *f = fopen(argv[3], "rb"); if (!f) return 2;
  fseek(f, 0, SEEK_END); long nrep = ftell(f) / 8; fseek(f, 0, SEEK_SET);
  uint64_t *R = malloc(8 * (nrep ? nrep : 1)); if (fread(R, 8, nrep, f) != (size_t)nrep) return 3; fclose(f);
  int haverange = argc >= 6 && strcmp(argv[4], "full");
  long first = haverange ? atol(argv[4]) : 0, last = haverange ? atol(argv[5]) : nrep; if (last > nrep) last = nrep;
  long bad = 0, badsize = 0;
  for (long r = first; r < last; r++) {
    uint64_t M = R[r]; if (__builtin_popcountll(M) != a || (N < 64 && (M >> N))) { badsize++; continue; }
    int T[64], k = 0; for (uint64_t x = M; x; x &= x - 1) T[k++] = __builtin_ctzll(x);
    int canon = 1;
    if (!full && a > 0 && !(M & 1)) canon = 0;
    int nt = full ? N : k;
    for (int ti = 0; ti < nt && canon; ti++) {
      int t = full ? ti : T[ti];
      for (int p = 0; p < nperm; p++) {
        uint64_t m = 0; for (int j = 0; j < k; j++) m |= 1ULL << img[p][T[j] ^ t];
        uint64_t z = m ^ M; if (z && (m & z & (~z + 1))) { canon = 0; break; }   // image is lex-smaller
      }
    }
    if (!canon) bad++;
  }
  qsort(R, nrep, 8, cmpu64); long dup = 0; for (long r = 1; r < nrep; r++) dup += R[r] == R[r - 1];
  printf("REPCHECK n=%d a=%d mode=%s file_reps=%ld checked=%ld noncanonical=%ld wrongsize=%ld duplicates=%ld\n",
         n, a, full ? "full" : "fast", nrep, last - first, bad, badsize, dup);
  return 0;
}
