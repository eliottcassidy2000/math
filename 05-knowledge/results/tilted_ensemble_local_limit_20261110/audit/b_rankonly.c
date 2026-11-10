// b_rankonly.c -- AUDIT (beyond the README; covers the packet's PREREGISTRATION.md P1/P2 at order 10):
// rank distribution of beta_R = (mu_R/2 mod 2) over the classes read from stdin, with class and labelled weights.
// mu_R(i,j) = M[i][j] + M[j][i], M[j][i] = sum_{S: j in S, i notin S} E[S][j] B[V-S][i]  (auditor's own code).
// usage: gentourng -q m [res/mod] | ./b_rankonly m      -> one line of raw counts (parts can be added)
#include <stdio.h>
#include <stdlib.h>
typedef long long i64; typedef unsigned long long u64;
#define MAXM 10
static int m, full, adj[MAXM][MAXM], inm[MAXM], outm[MAXM], sc[MAXM], perm_[MAXM], used_;
static i64 E[1 << MAXM][MAXM], B[1 << MAXM][MAXM], Mx[MAXM][MAXM];
static i64 autcount(int k) {
  if (k == m) return 1; i64 c = 0;
  for (int w = 0; w < m; w++) { if (used_ >> w & 1) continue; if (sc[w] != sc[k]) continue; int ok = 1;
    for (int j = 0; j < k && ok; j++) if (adj[j][k] != adj[perm_[j]][w]) ok = 0;
    if (ok) { perm_[k] = w; used_ |= 1 << w; c += autcount(k + 1); used_ &= ~(1 << w); } }
  return c;
}
int main(int argc, char **argv) {
  m = atoi(argv[1]); full = (1 << m) - 1; i64 mf = 1; for (int i = 2; i <= m; i++) mf *= i;
  char line[256]; u64 nt = 0, odd = 0, cls[MAXM + 2] = {0}; i64 lab[MAXM + 2] = {0}, minmu = 1LL << 60, lt = 0;
  while (fgets(line, sizeof line, stdin)) {
    if (line[0] != '0' && line[0] != '1') continue; int k = 0;
    for (int i = 0; i < m; i++) { inm[i] = outm[i] = 0; sc[i] = 0; adj[i][i] = 0; }
    for (int i = 0; i < m; i++) for (int j = i + 1; j < m; j++) { int b = line[k++] == '1'; adj[i][j] = b; adj[j][i] = !b;
      if (b) { outm[i] |= 1 << j; inm[j] |= 1 << i; sc[i]++; } else { outm[j] |= 1 << i; inm[i] |= 1 << j; sc[j]++; } }
    for (int S = 1; S <= full; S++) { int vv = S; while (vv) { int v = __builtin_ctz(vv); vv &= vv - 1; int R = S ^ (1 << v);
        if (!R) { E[S][v] = B[S][v] = 1; continue; } i64 e = 0, b = 0; int pp = R & inm[v]; while (pp) { int u = __builtin_ctz(pp); pp &= pp - 1; e += E[R][u]; }
        pp = R & outm[v]; while (pp) { int u = __builtin_ctz(pp); pp &= pp - 1; b += B[R][u]; } E[S][v] = e; B[S][v] = b; } }
    for (int j = 0; j < m; j++) for (int i = 0; i < m; i++) { Mx[j][i] = 0; if (i == j) continue; int rest = full & ~(1 << i) & ~(1 << j); i64 s = 0;
      for (int S0 = rest;; S0 = (S0 - 1) & rest) { int S = S0 | (1 << j); s += E[S][j] * B[full ^ S][i]; if (!S0) break; } Mx[j][i] = s; }
    int rows[MAXM];
    for (int i = 0; i < m; i++) { rows[i] = 0; for (int j = 0; j < m; j++) { if (i == j) continue; i64 mu = Mx[i][j] + Mx[j][i]; if (mu & 1) odd++; if (mu < minmu) minmu = mu; if ((mu >> 1) & 1) rows[i] |= 1 << j; } }
    int rk = 0; for (int c = 0; c < m; c++) { int p = -1; for (int i = rk; i < m; i++) if (rows[i] >> c & 1) { p = i; break; } if (p < 0) continue;
      int t = rows[p]; rows[p] = rows[rk]; rows[rk] = t; for (int i = 0; i < m; i++) if (i != rk && (rows[i] >> c & 1)) rows[i] ^= rows[rk]; rk++; }
    used_ = 0; i64 au = autcount(0); cls[rk]++; lab[rk] += mf / au; lt += mf / au; nt++;
  }
  printf("order %d classes %llu labelled %lld oddmu %llu minmu %lld", m, nt, lt, odd, minmu);
  for (int r = 0; r <= m; r++) printf(" r%d %llu %lld", r, cls[r], lab[r]);
  printf("\n"); return 0;
}
