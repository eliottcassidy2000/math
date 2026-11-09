/* AUDIT: complete census of I(Omega(T), b), b = 2..8, over all labelled tournaments of order <= 7, by the same
 * cycle-DP / packing-DP routine as a1_brute.c. Prints, per base, the alpha-sums m <= 200 NOT attained at order <= 7
 * (m = (I - 1)/b). Used to spot-check the attainment side of THM-4940 (6) for b = 4..8 (m <= 80).
 * Optional: "a9_small B TARGET NMIN NMAX ITERS SEED" hill-climbs for I(Omega, B) = TARGET (orders NMIN..NMAX <= 12)
 * and prints the tournament found (rows), to supply witnesses for values outside the order-7 census.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef long long ll;
static int n, A[12][12];
static ll cyc[1 << 12], F[1 << 12], dp[1 << 12][12];

static void cycles(void) {
  int N = 1 << n;
  memset(cyc, 0, sizeof(ll) * N);
  for (int s = 0; s < n; s++) {
    for (int S = 0; S < N; S++) if ((S >> s & 1) && !(S & ((1 << s) - 1))) for (int v = 0; v < n; v++) dp[S][v] = 0;
    dp[1 << s][s] = 1;
    for (int S = 1 << s; S < N; S++) {
      if (!(S >> s & 1) || (S & ((1 << s) - 1))) continue;
      int pc = __builtin_popcount(S);
      for (int v = 0; v < n; v++) { ll c = dp[S][v]; if (!c) continue;
        if (pc >= 3 && (pc & 1) && A[v][s]) cyc[S] += c;
        for (int t = s + 1; t < n; t++) if (!(S >> t & 1) && A[v][t]) dp[S | 1 << t][t] += c; }
    }
  }
}
static ll Ix(ll x) {
  int N = 1 << n; F[0] = 1;
  for (int S = 1; S < N; S++) { int m = __builtin_ctz(S), rest = S & ~(1 << m); ll val = F[rest];
    for (int T = rest;; T = (T - 1) & rest) { int U = T | (1 << m); if (cyc[U]) val += x * cyc[U] * F[S & ~U]; if (!T) break; }
    F[S] = val; }
  return F[N - 1];
}
static unsigned long long st = 0x2545F4914F6CDD1DULL;
static unsigned long long rnd(void) { unsigned long long z = (st += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL; z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL; return z ^ (z >> 31); }

int main(int argc, char **argv) {
  if (argc > 1) {
    ll B = atoll(argv[1]), target = atoll(argv[2]); int nmin = atoi(argv[3]), nmax = atoi(argv[4]);
    ll iters = atoll(argv[5]); st ^= strtoull(argv[6], 0, 10) * 0x9E3779B97F4A7C15ULL;
    for (ll it = 0; it < iters; it++) {
      n = nmin + rnd() % (nmax - nmin + 1);
      for (int i = 0; i < n; i++) { A[i][i] = 0; for (int j = i + 1; j < n; j++) { int b = rnd() & 1; A[i][j] = b; A[j][i] = !b; } }
      cycles(); ll cur = Ix(B), d = llabs(cur - target);
      for (int step = 0; step < 600 && d; step++) {
        int i = rnd() % n, j = rnd() % n; if (i == j) continue;
        int t = A[i][j]; A[i][j] = A[j][i]; A[j][i] = t;
        cycles(); ll v = Ix(B), d2 = llabs(v - target);
        if (d2 <= d || rnd() % 25 == 0) d = d2; else { t = A[i][j]; A[i][j] = A[j][i]; A[j][i] = t; }
      }
      if (!d) { printf("FOUND base %lld target %lld n=%d rows", B, target, n);
        for (int i = 0; i < n; i++) { printf(" "); for (int j = 0; j < n; j++) printf("%d", A[i][j]); } printf("\n"); return 0; }
    }
    printf("not found base %lld target %lld\n", B, target);
    return 0;
  }
  static char att[9][201];
  for (n = 1; n <= 7; n++) {
    int m = n * (n - 1) / 2;
    for (ll code = 0; code < (1LL << m); code++) {
      int k = 0; memset(A, 0, sizeof A);
      for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++, k++) { if (code >> k & 1) A[i][j] = 1; else A[j][i] = 1; }
      cycles();
      for (ll b = 2; b <= 8; b++) { ll v = Ix(b); ll a = (v - 1) / b; if ((v - 1) % b) { printf("NOT 1 mod b!\n"); return 1; } if (a <= 200) att[b][a] = 1; }
    }
  }
  for (int b = 2; b <= 8; b++) {
    printf("base %d: alpha-sums m <= 200 not attained by any tournament of order <= 7:", b);
    for (int a = 0; a <= 200; a++) if (!att[b][a]) printf(" %d", a);
    printf("\n");
  }
  return 0;
}
