/* minH.c - simulated annealing for upper bounds on the connectivity floor f_k(n) = min H over n-tournaments with
 * minimum semidegree >= k (= lambda >= k when n <= 4k + 1, THM-4811 (B)).  EXPLORATORY (upper bounds only).
 * Calibration: S17's exact f_3(n) = 171, 621, 2055, 5753, 11475 at n = 7..11 (THM-4942 (4)).
 * usage: minH n k seed restarts steps_per_restart
 * output: best H found and its code
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
static int N, K;
static uint32_t out[20];
static uint64_t rng = 0x243F6A8885A308D3ULL;
static inline uint64_t xr(void) { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return rng; }
static uint64_t *dp;
static uint64_t H(void) {
    size_t S = (size_t)1 << N;
    memset(dp, 0, S * N * sizeof(uint64_t));
    for (int a = 0; a < N; a++) dp[((size_t)1 << a) * N + a] = 1;
    for (size_t m = 1; m < S; m++) for (int a = 0; a < N; a++) {
        uint64_t c = dp[m * N + a]; if (!c) continue;
        uint32_t nx = out[a] & ~(uint32_t)m;
        while (nx) { int b = __builtin_ctz(nx); nx &= nx - 1; dp[(m | ((size_t)1 << b)) * N + b] += c; }
    }
    uint64_t h = 0; for (int a = 0; a < N; a++) h += dp[(S - 1) * N + a];
    return h;
}
static int minsemi(void) { int m = 99; for (int v = 0; v < N; v++) { int o = __builtin_popcount(out[v]), i = N - 1 - o; if (o < m) m = o; if (i < m) m = i; } return m; }
int main(int argc, char **argv) {
    N = atoi(argv[1]); K = atoi(argv[2]); rng ^= (uint64_t)atoll(argv[3]) * 0x9E3779B97F4A7C15ULL;
    int R = atoi(argv[4]); long steps = atol(argv[5]);
    for (int i = 0; i < 20; i++) xr();
    dp = malloc(((size_t)1 << N) * N * sizeof(uint64_t));
    uint64_t best = UINT64_MAX; uint32_t bestout[20];
    for (int r = 0; r < R; r++) {
        do { memset(out, 0, sizeof out);
             for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++) { if (xr() & 1) out[i] |= 1u << j; else out[j] |= 1u << i; }
        } while (minsemi() < K);
        uint64_t h = H();
        for (long t = 0; t < steps; t++) {
            double temp = 0.3 * (1.0 - (double)t / steps) + 0.005;
            int a = xr() % N, b = xr() % N; if (a == b) continue;
            out[a] ^= 1u << b; out[b] ^= 1u << a;
            if (minsemi() < K) { out[a] ^= 1u << b; out[b] ^= 1u << a; continue; }
            uint64_t h2 = H();
            if (h2 <= h || exp(-(log((double)h2) - log((double)h)) / temp) > (double)(xr() % 1000000) / 1e6) {
                h = h2;
                if (h < best) { best = h; memcpy(bestout, out, sizeof out); }
            } else { out[a] ^= 1u << b; out[b] ^= 1u << a; }
        }
    }
    printf("n=%d k=%d best H=%llu code ", N, K, (unsigned long long)best);
    for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++) putchar(bestout[i] >> j & 1 ? '1' : '0');
    printf("\n");
    return 0;
}
