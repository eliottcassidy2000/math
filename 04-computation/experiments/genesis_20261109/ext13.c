/* ext13.c - one-vertex extensions T = E + v of order-12 parents E, searching order-13 escapers with
 * minimum semidegree >= 4 (equivalently lambda >= 4 at n = 13 <= 4*4+1, THM-4811 (B))   (genesis lane, Q1-Q2).
 *
 * For a parent E (order 12) we evaluate, for every pattern x (x_a = [a -> v]):
 *   H(E + v)        by the attachment form of E,
 *   H((E - u) + v)  for every u in E, by the attachment form of E - u (order 11),
 * all updated incrementally along a Gray code.  H(T - v) = H(E) is fixed.
 * A pattern is reported when: min semidegree of T >= 4; H(T - w) in S(11) for all 13 w (exact, klein's census via
 * unattained_11.txt); and H(T) <= 545835 (= max H(12), proxy; S(12) membership is certified afterwards by wit.c).
 * The reported H values are re-checked by a direct subset DP before printing.
 * usage: ext13 < parents   (lines "12 code"); prints "C H code13 Hdel_0..Hdel_12"
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned long long u64;
typedef __int128 i128;
static unsigned char hole11[95096];

static void form(const uint32_t *o, int k, i128 *H0, i128 *L, i128 (*M)[16]) {
    size_t S = (size_t)1 << k;
    u64 *F = calloc(S * k, sizeof(u64)), *G = calloc(S * k, sizeof(u64));
    uint32_t in[16] = {0};
    for (int a = 0; a < k; a++) for (int b = 0; b < k; b++) if (o[a] >> b & 1) in[b] |= 1u << a;
    for (int a = 0; a < k; a++) { F[((size_t)1 << a) * k + a] = 1; G[((size_t)1 << a) * k + a] = 1; }
    for (size_t s = 1; s < S; s++) for (int a = 0; a < k; a++) {
        u64 c = F[s * k + a];
        if (c) { uint32_t nx = o[a] & ~(uint32_t)s; while (nx) { int b = __builtin_ctz(nx); nx &= nx - 1; F[(s | ((size_t)1 << b)) * k + b] += c; } }
        u64 d = G[s * k + a];
        if (d) { uint32_t nx = in[a] & ~(uint32_t)s; while (nx) { int b = __builtin_ctz(nx); nx &= nx - 1; G[(s | ((size_t)1 << b)) * k + b] += d; } }
    }
    size_t full = S - 1;
    u64 H = 0; for (int a = 0; a < k; a++) H += F[full * k + a];
    u64 W[16][16]; memset(W, 0, sizeof W);
    for (size_t s = 1; s < full; s++) { size_t r = full & ~s;
        for (int a = 0; a < k; a++) { u64 fa = F[s * k + a]; if (!fa) continue;
            for (int b = 0; b < k; b++) { u64 gb = G[r * k + b]; if (gb) W[a][b] += fa * gb; } } }
    *H0 = H;
    for (int a = 0; a < k; a++) {
        i128 l = (i128)F[full * k + a] - (i128)G[full * k + a];
        for (int b = 0; b < k; b++) l += W[a][b];
        L[a] = l;
        for (int b = 0; b < k; b++) M[a][b] = (a == b) ? 0 : (i128)W[a][b] + (i128)W[b][a];
    }
    free(F); free(G);
}
static u64 hdirect(const uint32_t *o, int n) {
    size_t S = (size_t)1 << n; u64 *dp = calloc(S * n, sizeof(u64));
    for (int a = 0; a < n; a++) dp[((size_t)1 << a) * n + a] = 1;
    for (size_t m = 1; m < S; m++) for (int a = 0; a < n; a++) { u64 c = dp[m * n + a]; if (!c) continue;
        uint32_t nx = o[a] & ~(uint32_t)m; while (nx) { int b = __builtin_ctz(nx); nx &= nx - 1; dp[(m | ((size_t)1 << b)) * n + b] += c; } }
    u64 h = 0; for (int a = 0; a < n; a++) h += dp[(S - 1) * n + a]; free(dp); return h;
}
static int inS11(i128 h) { return h > 0 && h <= 95095 && ((long)h & 1) && h != 7 && h != 21 && !hole11[(long)h]; }

int main(void) {
    FILE *f = fopen("unattained_11.txt", "r"); char line[512];
    while (fgets(line, sizeof line, f)) { if (line[0] == '#') continue; long h = atol(line); if (h > 0 && h <= 95095) hole11[h] = 1; }
    fclose(f);
    int n; char code[256]; long parents = 0, reported = 0;
    while (scanf("%d %255s", &n, code) == 2) {
        if (n != 12) continue;
        parents++;
        uint32_t E[16] = {0};
        int p = 0;
        for (int i = 0; i < 12; i++) for (int j = i + 1; j < 12; j++) { if (code[p] == '1') E[i] |= 1u << j; else E[j] |= 1u << i; p++; }
        i128 H0, L[16], M[16][16];
        form(E, 12, &H0, L, M);
        if (!inS11(H0)) continue;   /* T - v = E must be a non-creator */
        /* forms of E - u */
        i128 Hu[12], Lu[12][16], Mu[12][16][16];
        int map[12][12];  /* map[u][a] = index of a in E-u, or -1 */
        for (int u = 0; u < 12; u++) {
            uint32_t o[16] = {0}; int idx = 0, inv[12];
            for (int a = 0; a < 12; a++) { if (a == u) { map[u][a] = -1; continue; } map[u][a] = idx; inv[idx++] = a; }
            for (int i = 0; i < 11; i++) for (int j = 0; j < 11; j++) if (E[inv[i]] >> inv[j] & 1) o[i] |= 1u << j;
            form(o, 11, &Hu[u], Lu[u], Mu[u]);
        }
        /* degrees in E */
        int outE[12];
        for (int a = 0; a < 12; a++) outE[a] = __builtin_popcount(E[a]);
        /* Gray code over x */
        i128 val = H0, sx[16] = {0};
        i128 vu[12], su[12][16];
        for (int u = 0; u < 12; u++) { vu[u] = Hu[u]; for (int a = 0; a < 16; a++) su[u][a] = 0; }
        uint32_t x = 0;
        for (u64 t = 0; t < (1ULL << 12); t++) {
            if (t) {
                int c = __builtin_ctzll(t);
                int on = !(x >> c & 1);
                if (on) { val += L[c] - sx[c]; x |= 1u << c; for (int b = 0; b < 12; b++) sx[b] += M[b][c]; }
                else { val -= L[c] - sx[c]; x &= ~(1u << c); for (int b = 0; b < 12; b++) sx[b] -= M[b][c]; }
                for (int u = 0; u < 12; u++) {
                    int cc = map[u][c]; if (cc < 0) continue;
                    if (on) { vu[u] += Lu[u][cc] - su[u][cc]; for (int b = 0; b < 11; b++) su[u][b] += Mu[u][b][cc]; }
                    else { vu[u] -= Lu[u][cc] - su[u][cc]; for (int b = 0; b < 11; b++) su[u][b] -= Mu[u][b][cc]; }
                }
            }
            /* min semidegree >= 4 in T */
            int xv = __builtin_popcount(x);            /* in-degree of v = #a with a -> v */
            if (xv < 4 || 12 - xv < 4) continue;
            int ok = 1;
            for (int a = 0; a < 12 && ok; a++) { int o = outE[a] + (x >> a & 1); if (o < 4 || 12 - o < 4) ok = 0; }
            if (!ok) continue;
            if (val > 545835) continue;
            for (int u = 0; u < 12 && ok; u++) if (!inS11(vu[u])) ok = 0;
            if (!ok) continue;
            /* verify directly */
            uint32_t T[16]; memcpy(T, E, sizeof(uint32_t) * 12); T[12] = 0;
            for (int a = 0; a < 12; a++) { if (x >> a & 1) T[a] |= 1u << 12; else T[12] |= 1u << a; }
            u64 hT = hdirect(T, 13);
            if ((i128)hT != val) { fprintf(stderr, "MISMATCH T\n"); continue; }
            printf("C %llu ", hT);
            for (int i = 0; i < 13; i++) for (int j = i + 1; j < 13; j++) putchar(T[i] >> j & 1 ? '1' : '0');
            for (int u = 0; u < 12; u++) printf(" %lld", (long long)vu[u]);
            printf(" %lld\n", (long long)H0);
            fflush(stdout); reported++;
        }
    }
    fprintf(stderr, "parents %ld reported %ld\n", parents, reported);
    return 0;
}
