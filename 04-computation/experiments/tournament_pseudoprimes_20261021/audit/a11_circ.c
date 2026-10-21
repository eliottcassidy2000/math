/* AUDIT a11 (independent auditor, 2026-10-21): circulant tournaments Cay(Z_n, S), n odd.
 * For each connection set (or each multiplier class) compute
 *   H    = # Hamiltonian paths  (full subset DP over all start vertices when n <= 17, cross-checked with n * paths-from-0;
 *          for n > 17 only n * paths-from-0, via a popcount-layered DP that keeps two layers in memory);
 *   C    = # directed Hamiltonian cycles = sum over v with v->0 of #(paths from 0 ending at v);
 *   rho0 = # Hamiltonian paths P with P reversed = -P (reflection x -> -x), by a DP over the (n-1)/2 pairs {x,-x};
 *   mu   = |{a in Z_n^* : aS = S}|,  A = |S cap Z_n^*|.
 * Usage: ./a11_circ n all        (one line per multiplier class)
 *        ./a11_circ n set s1,s2,...  |  ./a11_circ n paley  |  ./a11_circ n rot
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
static int n;
static int inS[64];
static int gcdi(int a, int b) { while (b) { int t = a % b; a = b; b = t; } return a; }
static int arc(int x, int y) { return inS[((y - x) % n + n) % n]; }
/* full DP: dp[mask][v] = # paths covering mask ending at v (any start) */
static u64 H_full(void) {
    size_t N = (size_t)1 << n; u64 *dp = calloc(N * n, sizeof(u64));
    for (int v = 0; v < n; v++) dp[((size_t)1 << v) * n + v] = 1;
    for (size_t m = 1; m < N; m++)
        for (int v = 0; v < n; v++) { u64 x = dp[m * n + v]; if (!x) continue;
            for (int w = 0; w < n; w++) if (!(m >> w & 1) && arc(v, w)) dp[(m | ((size_t)1 << w)) * n + w] += x; }
    u64 H = 0; for (int v = 0; v < n; v++) H += dp[(N - 1) * n + v];
    free(dp); return H;
}
/* binomial table for ranking */
static u64 Bn[64][64];
static void binom(void) { for (int i = 0; i < 64; i++) { Bn[i][0] = 1; for (int j = 1; j <= i; j++) Bn[i][j] = Bn[i - 1][j - 1] + (j <= i - 1 ? Bn[i - 1][j] : 0); } }
/* colex rank of a k-subset of {0..K-1} given as bitmask */
static u64 rank_mask(u64 m) { u64 r = 0; int i = 0; while (m) { int b = __builtin_ctzll(m); r += (b >= i + 1) ? Bn[b][i + 1] : 0; i++; m &= m - 1; } return r; }
/* next k-subset in colex order (Gosper) */
static u64 gosper(u64 x) { u64 c = x & -x, r = x + c; return (((r ^ x) >> 2) / c) | r; }
/* paths starting at 0: vertices 1..n-1 mapped to bits 0..n-2. Layered DP. Returns #paths from 0 (all ends) and C. */
static void from0(u64 *P0, u64 *Cyc) {
    int K = n - 1; /* bits */
    u64 maxl = Bn[K][K / 2];
    u64 *cur = calloc(maxl * K, sizeof(u64)), *nxt = calloc(maxl * K, sizeof(u64));
    /* layer 1: single vertex w with 0 -> w */
    for (int b = 0; b < K; b++) if (arc(0, b + 1)) cur[rank_mask(1ULL << b) * K + b] = 1;
    for (int k = 1; k < K; k++) {
        u64 cnt = Bn[K][k + 1]; memset(nxt, 0, cnt * K * sizeof(u64));
        u64 x = (1ULL << (k + 1)) - 1;
        for (u64 idx = 0; idx < cnt; idx++, x = gosper(x)) {
            /* target mask x of size k+1; pull from x - w for each w in x */
            int bits[64], nb = 0; u64 t = x; while (t) { bits[nb++] = __builtin_ctzll(t); t &= t - 1; }
            /* prefix/suffix ranks: rank(x - bits[j]) = sum_{i<j} B[b_i][i+1] + sum_{i>j} B[b_i][i] */
            u64 pre[65], suf[65]; pre[0] = 0; for (int i = 0; i < nb; i++) pre[i + 1] = pre[i] + ((bits[i] >= i + 1) ? Bn[bits[i]][i + 1] : 0);
            suf[nb] = 0; for (int i = nb - 1; i >= 0; i--) suf[i] = suf[i + 1] + ((bits[i] >= i) ? Bn[bits[i]][i] : 0);
            for (int j = 0; j < nb; j++) {
                int w = bits[j]; u64 r = pre[j] + suf[j + 1]; u64 acc = 0; u64 *src = cur + r * K;
                for (int i = 0; i < nb; i++) { if (i == j) continue; int v = bits[i]; if (arc(v + 1, w + 1)) acc += src[v]; }
                nxt[idx * K + w] = acc;
            }
        }
        u64 *tmp = cur; cur = nxt; nxt = tmp;
    }
    u64 tot = 0, cyc = 0;
    for (int v = 0; v < K; v++) { tot += cur[v]; if (arc(v + 1, 0)) cyc += cur[v]; }
    *P0 = tot; *Cyc = cyc; free(cur); free(nxt);
}
/* rho0: sequences v_1..v_m (m=(n-1)/2), one from each pair {x,-x}, v_i -> v_{i+1}, v_m -> 0 */
static u64 rho0(void) {
    int m = (n - 1) / 2; size_t N = (size_t)1 << m;
    /* endpoint encoded as 2*pair + sign; pair j <-> {j+1, n-j-1} */
    u64 *dp = calloc(N * 2 * m, sizeof(u64));
    for (int j = 0; j < m; j++) for (int sg = 0; sg < 2; sg++) dp[((size_t)1 << j) * 2 * m + 2 * j + sg] = 1;
    for (size_t M = 1; M < N; M++) for (int e = 0; e < 2 * m; e++) { u64 x = dp[M * 2 * m + e]; if (!x) continue;
        int j = e / 2, v = (e & 1) ? n - j - 1 : j + 1;
        for (int j2 = 0; j2 < m; j2++) if (!(M >> j2 & 1)) for (int sg = 0; sg < 2; sg++) { int w = sg ? n - j2 - 1 : j2 + 1;
            if (arc(v, w)) dp[(M | ((size_t)1 << j2)) * 2 * m + 2 * j2 + sg] += x; } }
    u64 r = 0; for (int e = 0; e < 2 * m; e++) { int j = e / 2, v = (e & 1) ? n - j - 1 : j + 1; if (arc(v, 0)) r += dp[(N - 1) * 2 * m + e]; }
    free(dp); return r;
}
static void report(void) {
    int mu = 0, A = 0;
    for (int a = 1; a < n; a++) { if (gcdi(a, n) != 1) continue; int ok = 1; for (int s = 1; s < n; s++) if (inS[s] && !inS[(a * s) % n]) { ok = 0; break; } mu += ok; }
    for (int s = 1; s < n; s++) if (inS[s] && gcdi(s, n) == 1) A++;
    u64 P0, C; from0(&P0, &C);
    u64 Hf = (n <= 17) ? H_full() : 0;
    u64 r = rho0();
    printf("n=%d S={", n); int first = 1; for (int s = 1; s < n; s++) if (inS[s]) { printf(first ? "%d" : ",%d", s); first = 0; }
    printf("} mu=%d A=%d H=%llu", mu, A, (u64)n * P0);
    if (n <= 17) printf(" Hfull=%llu%s", Hf, Hf == (u64)n * P0 ? "" : " MISMATCH");
    printf(" C=%llu rho0=%llu\n", C, r); fflush(stdout);
}
int main(int argc, char **argv) {
    n = atoi(argv[1]); binom(); int m = (n - 1) / 2;
    if (!strcmp(argv[2], "all")) {
        int nc = 0;
        unsigned char *seen = calloc(1u << m, 1);
        for (unsigned code = 0; code < (1u << m); code++) {
            if (seen[code]) continue; nc++;
            memset(inS, 0, sizeof inS);
            for (int j = 0; j < m; j++) { int s = (code >> j & 1) ? n - (j + 1) : j + 1; inS[s] = 1; }
            /* mark the multiplier orbit */
            for (int a = 1; a < n; a++) { if (gcdi(a, n) != 1) continue; unsigned c2 = 0;
                for (int j = 0; j < m; j++) { int s = (a * (j + 1)) % n; /* image of pair j's +-element */ }
                for (int s = 1; s < n; s++) if (inS[s]) { int t = (a * s) % n; int j = (t <= m) ? t - 1 : n - t - 1; if (t > m) c2 |= 1u << j; }
                seen[c2] = 1; }
            report();
        }
        printf("n=%d classes=%d\n", n, nc);
    } else {
        memset(inS, 0, sizeof inS);
        if (!strcmp(argv[2], "paley")) { for (int x = 1; x < n; x++) inS[(x * x) % n] = 1; }
        else if (!strcmp(argv[2], "rot")) { for (int s = 1; s <= m; s++) inS[s] = 1; }
        else { char *t = strtok(argv[3], ","); while (t) { inS[atoi(t) % n] = 1; t = strtok(0, ","); } }
        int cnt = 0; for (int s = 1; s < n; s++) { cnt += inS[s]; if (inS[s] && inS[n - s]) { printf("not a tournament\n"); return 1; } }
        if (cnt != m) { printf("not a tournament (|S|=%d)\n", cnt); return 1; }
        report();
    }
    return 0;
}
