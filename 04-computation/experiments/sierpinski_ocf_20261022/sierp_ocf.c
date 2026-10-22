/* sierp_ocf.c -- per-tournament OCF and Steiner-quasigroup data, for gentourng output on stdin (n <= 10).
 *   H        : Hamiltonian paths (Held-Karp);
 *   alpha_k  : sets of k pairwise vertex-disjoint directed odd cycles (OCF: H = sum_k 2^k alpha_k, checked);
 *   Im1      : I(Omega,-1) = sum (-1)^k alpha_k  (= minus the reduced Euler characteristic of Ind(Omega));
 *   J        : sum_k J_k alpha_k, J_k = (2^k - (-1)^k)/3 the Jacobsthal numbers, so that H = 3 J + Im1;
 *   S0,S1,S2 : Steiner signature: label chi(v) = outdeg(v) mod 3, reduce the label string of each Hamiltonian
 *              path by x*y = -x-y mod 3, i.e. R(P) = (-1)^(n-1) sum_i binom(n-1,i) chi(v_i) mod 3; S_r = #{P : R = r};
 *   aut      : |Aut(T)| by brute force over permutations (only if n <= 8, else 0).
 * mac-mini-2026-10-22-sierp.   usage: gentourng N | ./sierp_ocf  (one line out per tournament) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef long long ll; typedef uint64_t u64;
static int n; static uint32_t out[16];
static int binom3[16];
static u64 cyc[1 << 10]; static u64 F[1 << 10][5];
static u64 dpH[1 << 10][10];
static u64 dpS[1 << 10][10][3];
static int perm[10], used[10]; static long autc;
static void perms(int k){
  if (k == n){ for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (i != j && ((out[i] >> j) & 1) != ((out[perm[i]] >> perm[j]) & 1)) return; autc++; return; }
  for (int v = 0; v < n; v++) if (!used[v]){ used[v] = 1; perm[k] = v; perms(k + 1); used[v] = 0; }
}
int main(void){
  char line[256];
  while (fgets(line, sizeof line, stdin)){
    if (line[0] == '>') continue;
    int L = (int)strlen(line); while (L && (line[L-1] == '\n' || line[L-1] == '\r')) line[--L] = 0;
    n = 0; while (n * (n - 1) / 2 < L) n++;
    if (n * (n - 1) / 2 != L || n > 10) { fprintf(stderr, "bad line\n"); return 1; }
    for (int i = 0; i < n; i++) out[i] = 0;
    int p = 0;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++, p++){ if (line[p] == '1') out[i] |= 1u << j; else out[j] |= 1u << i; }
    uint32_t N = 1u << n;
    /* binomials mod 3 with sign (-1)^(n-1) */
    for (int i = 0; i < n; i++){ ll b = 1; for (int t = 0; t < i; t++) b = b * (n - 1 - t) / (t + 1); int s = (int)(b % 3); if ((n - 1) & 1) s = (3 - s) % 3; binom3[i] = s; }
    int chi[16]; for (int v = 0; v < n; v++) chi[v] = __builtin_popcount(out[v]) % 3;
    /* H and Steiner signature by DP over (mask, last, residue) */
    memset(dpH, 0, sizeof(u64) * N * 10); memset(dpS, 0, sizeof(u64) * N * 10 * 3);
    for (int v = 0; v < n; v++){ dpH[1u << v][v] = 1; dpS[1u << v][v][(binom3[0] * chi[v]) % 3] = 1; }
    for (uint32_t m = 1; m < N; m++){
      int pos = __builtin_popcount(m);
      for (int v = 0; v < n; v++) if ((m >> v) & 1){
        u64 h = dpH[m][v]; if (!h) continue;
        uint32_t nxt = out[v] & ~m;
        for (uint32_t q = nxt; q; q &= q - 1){ int w = __builtin_ctz(q); uint32_t m2 = m | (1u << w);
          dpH[m2][w] += h;
          int add = (binom3[pos] * chi[w]) % 3;
          for (int r = 0; r < 3; r++) dpS[m2][w][(r + add) % 3] += dpS[m][v][r]; }
      }
    }
    u64 H = 0, S[3] = {0, 0, 0};
    for (int v = 0; v < n; v++){ H += dpH[N - 1][v]; for (int r = 0; r < 3; r++) S[r] += dpS[N - 1][v][r]; }
    /* odd cycles per vertex set: Hamiltonian cycles of T[U], start at min vertex */
    static u64 dp[1 << 10][10];
    memset(dp, 0, sizeof dp); memset(cyc, 0, sizeof(u64) * N);
    for (uint32_t m = 1; m < N; m++){
      int v0 = __builtin_ctz(m);
      if (m == (1u << v0)) { dp[m][v0] = 1; continue; }
      uint32_t rest = m & ~(1u << v0);
      for (uint32_t r = rest; r; r &= r - 1){ int u = __builtin_ctz(r); uint32_t pm = m & ~(1u << u); u64 s = 0;
        for (uint32_t q = pm; q; q &= q - 1){ int w = __builtin_ctz(q); if ((out[w] >> u) & 1) s += dp[pm][w]; }
        dp[m][u] = s; }
      int sz = __builtin_popcount(m);
      if (sz >= 3 && (sz & 1)){ u64 c = 0; for (uint32_t r = rest; r; r &= r - 1){ int u = __builtin_ctz(r); if ((out[u] >> v0) & 1) c += dp[m][u]; } cyc[m] = c; }
    }
    int K = n / 3; memset(F, 0, sizeof(u64) * N * 5); F[0][0] = 1;
    u64 alpha[5] = {1, 0, 0, 0, 0};
    for (uint32_t W = 1; W < N; W++){
      int v0 = __builtin_ctz(W); uint32_t Wp = W & ~(1u << v0);
      for (uint32_t sub = Wp; ; sub = (sub - 1) & Wp){
        uint32_t U = sub | (1u << v0); u64 c = cyc[U];
        if (c) for (int k = 1; k <= K; k++) F[W][k] += c * F[W & ~U][k - 1];
        if (!sub) break;
      }
      for (int k = 1; k <= K; k++) alpha[k] += F[W][k];
    }
    ll Hocf = 0, Im1 = 0, J = 0; ll Jk[5] = {0, 1, 1, 3, 5};
    for (int k = 0; k <= K; k++){ Hocf += (ll)alpha[k] << k; Im1 += (k & 1) ? -(ll)alpha[k] : (ll)alpha[k]; J += Jk[k] * (ll)alpha[k]; }
    if ((u64)Hocf != H){ fprintf(stderr, "OCF mismatch %s\n", line); return 1; }
    if (3 * J + Im1 != Hocf){ fprintf(stderr, "Jacobsthal mismatch\n"); return 1; }
    autc = 0; if (n <= 8){ memset(used, 0, sizeof used); perms(0); }
    char scores[32]; for (int v = 0; v < n; v++) scores[v] = '0' + __builtin_popcount(out[v]); scores[n] = 0;
    printf("%s H=%llu a=%llu,%llu,%llu,%llu Im1=%lld J=%lld S=%llu,%llu,%llu aut=%ld sc=%s\n", line, (unsigned long long)H,
           (unsigned long long)alpha[1], (unsigned long long)alpha[2], (unsigned long long)alpha[3], (unsigned long long)alpha[4],
           Im1, J, (unsigned long long)S[0], (unsigned long long)S[1], (unsigned long long)S[2], autc, scores);
  }
  return 0;
}
