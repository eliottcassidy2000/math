/* port statistics of every tournament class and port: for x = 2, 3:
 * h = I(Omega(T), x), hu = I(Omega(T - u), x), E/O = sum over paths P from u of even/odd order of I(Omega(T - V(P)), x).
 * input: classgen file (first token = canonical code), n.  output: code u h2 E2 O2 hu2 h3 E3 O3 hu3
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll;
static int n, A[12][12];
static ll cyc[1 << 12], I2[1 << 12], I3[1 << 12], P[1 << 12][12], D[1 << 12][12];
static void cycles(void) {
  int N = 1 << n;
  for (int S = 0; S < N; S++) cyc[S] = 0;
  for (int s = 0; s < n; s++) {
    int base = 1 << s;
    for (int S = 0; S < N; S++) if ((S & base) && !(S & (base - 1))) for (int v = 0; v < n; v++) P[S][v] = 0;
    P[base][s] = 1;
    for (int S = 0; S < N; S++) {
      if (!(S & base) || (S & (base - 1))) continue;
      for (int v = 0; v < n; v++) { ll c = P[S][v]; if (!c) continue;
        for (int w = s + 1; w < n; w++) if (!(S >> w & 1) && A[v][w]) P[S | 1 << w][w] += c;
        if (A[v][s] && __builtin_popcount(S) >= 3) cyc[S] += c; }
    }
  }
}
static void packing(ll x, ll *F) {
  int N = 1 << n; F[0] = 1;
  for (int S = 1; S < N; S++) {
    int m = __builtin_ctz(S); ll val = F[S & ~(1 << m)]; int rest = S & ~(1 << m);
    for (int T = rest;; T = (T - 1) & rest) {
      int U = T | (1 << m); int pc = __builtin_popcount(U);
      if (pc % 2 == 1 && pc >= 3 && cyc[U]) val += x * cyc[U] * F[S & ~U];
      if (T == 0) break;
    }
    F[S] = val;
  }
}
int main(int argc, char **argv) {
  FILE *f = fopen(argv[1], "r"); n = atoi(argv[2]);
  char line[4096]; unsigned long long code;
  int N = 1 << n, FULL = N - 1;
  while (fgets(line, sizeof line, f)) {
    if (sscanf(line, "%llu", &code) != 1) continue;
    memset(A, 0, sizeof A);
    int b = 0;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (code >> b & 1) A[i][j] = 1; else A[j][i] = 1; b++; }
    cycles(); packing(2, I2); packing(3, I3);
    for (int u = 0; u < n; u++) {
      /* D[S][v]: number of paths starting at u with vertex set S ending at v */
      for (int S = 0; S < N; S++) for (int v = 0; v < n; v++) D[S][v] = 0;
      D[1 << u][u] = 1;
      for (int S = 0; S < N; S++) { if (!(S >> u & 1)) continue;
        for (int v = 0; v < n; v++) { ll c = D[S][v]; if (!c) continue;
          for (int w = 0; w < n; w++) if (!(S >> w & 1) && A[v][w]) D[S | 1 << w][w] += c; } }
      ll E2 = 0, O2 = 0, E3 = 0, O3 = 0;
      for (int S = 0; S < N; S++) { if (!(S >> u & 1)) continue;
        ll pc = 0; for (int v = 0; v < n; v++) pc += D[S][v];
        if (!pc) continue;
        int rest = FULL & ~S;
        if (__builtin_popcount(S) % 2 == 0) { E2 += pc * I2[rest]; E3 += pc * I3[rest]; }
        else { O2 += pc * I2[rest]; O3 += pc * I3[rest]; } }
      int Ru = FULL & ~(1 << u);
      printf("%llu %d %lld %lld %lld %lld %lld %lld %lld %lld\n", code, u, I2[FULL], E2, O2, I2[Ru], I3[FULL], E3, O3, I3[Ru]);
    }
  }
  return 0;
}
