// canon.c -- canonical form of vertex subsets of Q_n (n <= 7) under the full group Aut(Q_n)
// (x -> pi(x xor t), 2^n n! elements, enumerated explicitly).  Canonical form = the minimum mask
// (bit v = vertex v; for n = 7 the 128-bit mask is compared as an unsigned 128-bit integer).
// Also prints the order of the setwise stabilizer.  Independent of the search code.
// usage: canon n < sets      (each input line: a set as a list of vertex numbers, or "...set=v1,v2,...")
// output per line: <canonical mask hex (32 hex digits for n=7, 16 otherwise)> <stabilizer order> <size>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned __int128 u128;
static int n, N, nperm, perms[5040][7];
static void gen_perms(void) {
  int idx[7] = {0}; nperm = 0;
  for (;;) {
    int used = 0, ok = 1;
    for (int i = 0; i < n; i++) { if (used >> idx[i] & 1) ok = 0; used |= 1 << idx[i]; }
    if (ok) { memcpy(perms[nperm], idx, sizeof idx); nperm++; }
    int i = 0; while (i < n && ++idx[i] == n) { idx[i] = 0; i++; } if (i == n) break;
  }
}
int main(int argc, char **argv) {
  if (argc < 2) return 1;
  n = atoi(argv[1]); if (n < 1 || n > 7) return 1; N = 1 << n;
  gen_perms();
  static int img[5040][128];
  for (int p = 0; p < nperm; p++) for (int v = 0; v < N; v++) { int u = 0; for (int i = 0; i < n; i++) if (v >> i & 1) u |= 1 << perms[p][i]; img[p][v] = u; }
  char line[1 << 16];
  while (fgets(line, sizeof line, stdin)) {
    char *p = strstr(line, "set="); p = p ? p + 4 : line;
    int S[128], k = 0;
    while (*p) { while (*p && (*p < '0' || *p > '9')) p++; if (!*p) break; S[k++] = (int)strtol(p, &p, 10); if (k >= 128) break; }
    if (!k) continue;
    u128 orig = 0; for (int j = 0; j < k; j++) { if (S[j] < 0 || S[j] >= N) { fprintf(stderr, "bad vertex\n"); return 2; } orig |= (u128)1 << S[j]; }
    u128 best = ~(u128)0; long stab = 0;
    for (int pp = 0; pp < nperm; pp++) for (int t = 0; t < N; t++) {
      u128 m = 0; for (int j = 0; j < k; j++) m |= (u128)1 << img[pp][S[j] ^ t];
      if (m < best) best = m;
      if (m == orig) stab++;
    }
    if (n == 7) printf("%016llx%016llx %ld %d\n", (unsigned long long)(best >> 64), (unsigned long long)best, stab, k);
    else printf("%016llx %ld %d\n", (unsigned long long)best, stab, k);
  }
  return 0;
}
