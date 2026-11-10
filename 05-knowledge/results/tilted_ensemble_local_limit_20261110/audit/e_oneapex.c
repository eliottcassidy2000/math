// e_oneapex.c -- AUDIT check (e): the remark of section 4.2.  Written from scratch by the auditor.
// Builds the ACTUAL one-apex tournaments for the word  1^(d-1) 0 b 1^d 0  (length L = 2d+2, plus the apex) and counts
// Hamiltonian paths by a plain subset DP (no recurrence, no cut formula), for both ways of reading the word:
//   reading S (source first): letter t (t = 1..L, left to right) is the chain vertex that beats all later letters' vertices;
//   reading K (sink first)  : letter t is the chain vertex that loses to all later letters' vertices.
// In both readings the apex beats exactly the chain vertices whose letter is 1.
// Also evaluates, for every word, the closed form derived by the auditor from (4.1) for a transitive base (source-first):
//   H = [w_1 = 1] + [w_L = 0] + sum_{i<j, w_i=1, w_j=0} 2^(j-i-1) + #{j : w_j = 0, w_{j+1} = 1}.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
static u64 ham(int n, const int *outm) {           /* paths ending at v on mask */
  size_t sz = (size_t)1 << n; u64 *dp = calloc(sz * n, sizeof(u64)); if (!dp) { fprintf(stderr, "no memory\n"); exit(1); }
  int inm[32]; for (int v = 0; v < n; v++) { inm[v] = 0; for (int u = 0; u < n; u++) if (outm[u] >> v & 1) inm[v] |= 1 << u; }
  for (size_t S = 1; S < sz; S++) { unsigned vv = (unsigned)S; while (vv) { int v = __builtin_ctz(vv); vv &= vv - 1; unsigned R = (unsigned)S ^ (1u << v);
      if (!R) { dp[S * n + v] = 1; continue; } u64 t = 0; unsigned pp = R & (unsigned)inm[v]; while (pp) { int u = __builtin_ctz(pp); pp &= pp - 1; t += dp[(size_t)R * n + u]; } dp[S * n + v] = t; } }
  u64 H = 0; for (int v = 0; v < n; v++) H += dp[(sz - 1) * n + v]; free(dp); return H;
}
static u64 build_and_count(const char *w, int sinkfirst) {
  int L = (int)strlen(w), n = L + 1, outm[32]; memset(outm, 0, sizeof outm);
  /* chain vertex of letter t (0-based) is t; source-first: t beats t' for t < t'; sink-first: t beats t' for t > t' */
  for (int s = 0; s < L; s++) for (int t = s + 1; t < L; t++) { if (!sinkfirst) outm[s] |= 1 << t; else outm[t] |= 1 << s; }
  for (int t = 0; t < L; t++) { if (w[t] == '1') outm[L] |= 1 << t; else outm[t] |= 1 << L; }
  return ham(n, outm);
}
static u64 closed(const char *w) {
  int L = (int)strlen(w); u64 H = (w[0] == '1') + (w[L - 1] == '0');
  for (int i = 0; i < L; i++) for (int j = i + 1; j < L; j++) if (w[i] == '1' && w[j] == '0') H += 1ULL << (j - i - 1);
  for (int j = 0; j + 1 < L; j++) if (w[j] == '0' && w[j + 1] == '1') H++;
  return H;
}
int main(int argc, char **argv) {
  int dmax = argc > 1 ? atoi(argv[1]) : 8, bad = 0;
  /* control: closed form against the brute-force DP on all words of length <= 9 (source-first) */
  for (int L = 1; L <= 9; L++) for (int x = 0; x < (1 << L); x++) { char w[16]; for (int t = 0; t < L; t++) w[t] = (x >> t & 1) ? '1' : '0'; w[L] = 0;
      if (closed(w) != build_and_count(w, 0)) bad++; }
  printf("control: closed form vs brute-force DP on all %d words of length <= 9 (source-first reading): %d mismatches\n", (1 << 10) - 2, bad);
  printf(" d  b  word                     claimed    H(source-first)  H(sink-first)   closed form\n");
  for (int d = 2; d <= dmax; d++) for (int b = 0; b <= 1; b++) {
    char w[64]; int k = 0; for (int t = 0; t < d - 1; t++) w[k++] = '1'; w[k++] = '0'; w[k++] = b ? '1' : '0'; for (int t = 0; t < d; t++) w[k++] = '1'; w[k++] = '0'; w[k] = 0;
    long long lo = (1LL << (2 * d + 1)) - 3 * (1LL << (d - 1)) - 1, hi = lo + 2;
    u64 hs = build_and_count(w, 0), hk = build_and_count(w, 1);
    printf("%2d  %d  %-22s %5lld/%-5lld %12llu %14llu %12llu   %s\n", d, b, w, lo, hi, hs, hk, closed(w),
           (hs == (u64)(b ? hi : lo)) ? "source-first matches (b=0 -> minus, b=1 -> plus)" : (hs == (u64)(b ? lo : hi)) ? "source-first matches with b swapped" : "source-first does NOT match");
  }
  /* closed form only, larger d (source-first) */
  int ok = 1;
  for (int d = 2; d <= 28; d++) for (int b = 0; b <= 1; b++) {
    char w[80]; int k = 0; for (int t = 0; t < d - 1; t++) w[k++] = '1'; w[k++] = '0'; w[k++] = b ? '1' : '0'; for (int t = 0; t < d; t++) w[k++] = '1'; w[k++] = '0'; w[k] = 0;
    u64 want = (1ULL << (2 * d + 1)) - 3 * (1ULL << (d - 1)) + (b ? 1 : -1); if (closed(w) != want) ok = 0; }
  printf("closed form, d = 2..28, source-first: b=0 gives 2^(2d+1)-3*2^(d-1)-1 and b=1 gives +1: %s\n", ok ? "yes" : "NO");
  return 0;
}
