/* wit.c - certify membership h in S(m) by an explicit witness (genesis lane, mac-mini-2026-11-09).
 *
 * For random base tournaments T of order m-1 (arc i->j for i<j in a random vertex order with probability p, else j->i),
 * compute the attachment form of a new vertex v (KNOWN: HYP-2879, THM-4097, THM-4114; restated in THM-4942 (1)):
 *     H(T + v_x) = H(T) + sum_a x_a L_a - sum_{a<b} x_a x_b M_ab,    x_a = [a -> v],
 *     L_a = f[V][a] - g[V][a] + sum_b W[a][b],  M_ab = W[a][b] + W[b][a],
 *     W[a][b] = #{(path ending at a, path starting at b) partitioning V}.
 * All 2^(m-1) patterns are scanned in Gray-code order; any value in the target list is reported with the base code and
 * pattern, and RE-VERIFIED by a direct subset DP of the m-tournament before printing.
 * usage: wit m p seed maxbases target1 target2 ...      (m <= 17)
 * output: "W h code_m"  (code_m = upper-triangle code of the m-tournament, '1' iff i -> j)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned long long u64;
typedef __int128 i128;
#define MAXM 17
static int k;                 /* base order = m - 1 */
static uint32_t out[MAXM];
static u64 *F, *G;            /* F[S*k+a]: paths covering S ending at a; G[S*k+a]: covering S starting at a */
static u64 rng = 0x9E3779B97F4A7C15ULL;
static inline u64 xr(void) { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return rng; }

static u64 hdirect(const uint32_t *o, int n) {
    size_t S = (size_t)1 << n;
    u64 *dp = calloc(S * n, sizeof(u64));
    for (int a = 0; a < n; a++) dp[((size_t)1 << a) * n + a] = 1;
    for (size_t m = 1; m < S; m++) for (int a = 0; a < n; a++) {
        u64 c = dp[m * n + a]; if (!c) continue;
        uint32_t nx = o[a] & ~(uint32_t)m;
        while (nx) { int b = __builtin_ctz(nx); nx &= nx - 1; dp[(m | ((size_t)1 << b)) * n + b] += c; }
    }
    u64 h = 0; for (int a = 0; a < n; a++) h += dp[(S - 1) * n + a];
    free(dp); return h;
}

int main(int argc, char **argv) {
    if (argc < 6) { fprintf(stderr, "usage: wit m p seed maxbases targets...\n"); return 1; }
    int m = atoi(argv[1]); double p = atof(argv[2]); rng ^= (u64)atoll(argv[3]) * 0xD1B54A32D192ED03ULL;
    long maxbases = atol(argv[4]);
    int nt = argc - 5; u64 *tg = malloc(nt * sizeof(u64)); int *done = calloc(nt, sizeof(int));
    for (int i = 0; i < nt; i++) tg[i] = strtoull(argv[5 + i], 0, 10);
    for (int i = 0; i < 20; i++) xr();
    k = m - 1;
    size_t S = (size_t)1 << k;
    F = malloc(S * k * sizeof(u64)); G = malloc(S * k * sizeof(u64));
    int remaining = nt;
    for (long base = 0; base < maxbases && remaining > 0; base++) {
        int perm[MAXM];
        for (int i = 0; i < k; i++) perm[i] = i;
        for (int i = k - 1; i > 0; i--) { int j = xr() % (i + 1); int t = perm[i]; perm[i] = perm[j]; perm[j] = t; }
        memset(out, 0, sizeof out);
        for (int i = 0; i < k; i++) for (int j = i + 1; j < k; j++) {
            int a = perm[i], b = perm[j];
            if ((double)(xr() % 1000000) / 1e6 < p) out[a] |= 1u << b; else out[b] |= 1u << a;
        }
        /* forward table F and backward table G */
        memset(F, 0, S * k * sizeof(u64)); memset(G, 0, S * k * sizeof(u64));
        uint32_t in[MAXM] = {0};
        for (int a = 0; a < k; a++) for (int b = 0; b < k; b++) if (out[a] >> b & 1) in[b] |= 1u << a;
        for (int a = 0; a < k; a++) { F[((size_t)1 << a) * k + a] = 1; G[((size_t)1 << a) * k + a] = 1; }
        for (size_t s = 1; s < S; s++) for (int a = 0; a < k; a++) {
            u64 c = F[s * k + a];
            if (c) { uint32_t nx = out[a] & ~(uint32_t)s; while (nx) { int b = __builtin_ctz(nx); nx &= nx - 1; F[(s | ((size_t)1 << b)) * k + b] += c; } }
            u64 d = G[s * k + a];
            if (d) { uint32_t nx = in[a] & ~(uint32_t)s; while (nx) { int b = __builtin_ctz(nx); nx &= nx - 1; G[(s | ((size_t)1 << b)) * k + b] += d; } }
        }
        size_t full = S - 1;
        u64 H = 0; for (int a = 0; a < k; a++) H += F[full * k + a];
        static u64 W[MAXM][MAXM];
        memset(W, 0, sizeof W);
        for (size_t s = 1; s < full; s++) {
            size_t r = full & ~s;
            for (int a = 0; a < k; a++) { u64 fa = F[s * k + a]; if (!fa) continue;
                for (int b = 0; b < k; b++) { u64 gb = G[r * k + b]; if (gb) W[a][b] += fa * gb; } }
        }
        i128 L[MAXM], M[MAXM][MAXM];
        for (int a = 0; a < k; a++) {
            i128 l = (i128)F[full * k + a] - (i128)G[full * k + a];
            for (int b = 0; b < k; b++) l += W[a][b];
            L[a] = l;
            for (int b = 0; b < k; b++) M[a][b] = (a == b) ? 0 : (i128)W[a][b] + (i128)W[b][a];
        }
        /* Gray code over x in {0,1}^k */
        i128 val = H, sx[MAXM];
        for (int a = 0; a < k; a++) sx[a] = 0;   /* sx[a] = sum_b x_b M[a][b] */
        uint32_t x = 0;
        for (u64 t = 0; t < ((u64)1 << k); t++) {
            if (t) {
                int c = __builtin_ctzll(t);
                if (x >> c & 1) { val -= L[c] - sx[c]; x &= ~(1u << c); for (int b = 0; b < k; b++) sx[b] -= M[b][c]; }
                else { val += L[c] - sx[c]; x |= 1u << c; for (int b = 0; b < k; b++) sx[b] += M[b][c]; }
            }
            for (int i = 0; i < nt; i++) if (!done[i] && val == (i128)tg[i]) {
                uint32_t o2[MAXM]; memcpy(o2, out, sizeof o2); o2[k] = 0;
                for (int a = 0; a < k; a++) { if (x >> a & 1) o2[a] |= 1u << k; else o2[k] |= 1u << a; }
                u64 hv = hdirect(o2, m);
                if (hv != tg[i]) { fprintf(stderr, "MISMATCH form %llu direct %llu\n", tg[i], hv); continue; }
                printf("W %llu ", tg[i]);
                for (int a = 0; a < m; a++) for (int b = a + 1; b < m; b++) putchar(o2[a] >> b & 1 ? '1' : '0');
                printf("\n"); fflush(stdout);
                done[i] = 1; remaining--;
            }
        }
    }
    fprintf(stderr, "remaining %d\n", remaining);
    return 0;
}
