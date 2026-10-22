/* aud_fwdpoly.c -- AUDIT (independent of forest_numbers.c).  Forward-step polynomial
 *     F(T, x) = sum over orderings sigma of [n] of x^(#i : sigma_i -> sigma_(i+1) is an arc of T)  (mod p)
 * by a LAYERED "pull" dynamic program: layer k holds g[S][j] = sum over orderings of the k-set S ending at its j-th
 * smallest element.  Sets of one size are ranked by the combinatorial number system (colex = Gosper order), so only
 * two layers are in memory (about 250 MB at n = 23, uint64).  Transition (explicit weights, no (x-1) trick):
 *     g[S][w] = sum_(v in S \ {w}) g[S \ {w}][v] * (x if v -> w else 1).
 * Input: one upper-triangle code on stdin ('1' at the k-th pair (i<j), row-major, iff i -> j).
 * Usage: aud_fwdpoly p x1 [x2 ...]   (p < 2^62; accumulation in unsigned __int128).  Prints "x F(T,x) mod p". */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef unsigned __int128 u128;
static uint64_t binom[33][33];

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: aud_fwdpoly p x1 [x2 ...] < code\n"); return 1; }
    uint64_t p = strtoull(argv[1], 0, 10);
    char *line = malloc(1 << 16);
    if (!fgets(line, 1 << 16, stdin)) return 1;
    int len = (int)strlen(line);
    while (len && (line[len - 1] == '\n' || line[len - 1] == '\r' || line[len - 1] == ' ')) line[--len] = 0;
    int n = 1; while (n * (n - 1) / 2 < len) n++;
    if (n * (n - 1) / 2 != len || n > 30) { fprintf(stderr, "bad code length %d\n", len); return 1; }
    static unsigned char beats[32][32];
    int q = 0;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++, q++) {
        if (line[q] == '1') { beats[i][j] = 1; beats[j][i] = 0; } else { beats[i][j] = 0; beats[j][i] = 1; }
    }
    for (int a = 0; a <= 32; a++) { binom[a][0] = 1; for (int b = 1; b <= a; b++) binom[a][b] = binom[a - 1][b - 1] + (b <= a - 1 ? binom[a - 1][b] : 0); }
    size_t maxent = 0;
    for (int k = 1; k <= n; k++) if (binom[n][k] * (uint64_t)k > maxent) maxent = binom[n][k] * (uint64_t)k;
    uint64_t *old = malloc(maxent * sizeof(uint64_t)), *cur = malloc(maxent * sizeof(uint64_t));
    if (!old || !cur) { fprintf(stderr, "alloc failed\n"); return 1; }
    for (int ai = 2; ai < argc; ai++) {
        uint64_t x = strtoull(argv[ai], 0, 10) % p;
        uint64_t wt[32][32];
        for (int v = 0; v < n; v++) for (int w = 0; w < n; w++) wt[v][w] = beats[v][w] ? x : 1 % p;
        /* layer 1: S = {v}, rank v (colex), one entry */
        for (int v = 0; v < n; v++) old[v] = 1 % p;
        for (int k = 2; k <= n; k++) {
            uint32_t S = (1u << k) - 1, lim = (n == 32) ? 0 : (1u << n);
            uint64_t r = 0;
            int el[32];
            uint64_t A[33], B[33];
            while (S < lim) {
                int c = 0;
                for (uint32_t t = S; t; t &= t - 1) el[c++] = __builtin_ctz(t);
                /* prefix sums of C(s_i, i+1) and suffix sums of C(s_i, i) */
                A[0] = 0; for (int i = 0; i < k; i++) A[i + 1] = A[i] + binom[el[i]][i + 1];
                B[k] = 0; for (int i = k - 1; i >= 0; i--) B[i] = B[i + 1] + binom[el[i]][i];
                uint64_t *dst = cur + r * k;
                for (int j = 0; j < k; j++) {
                    int w = el[j];
                    uint64_t rs = A[j] + B[j + 1];            /* rank of S \ {w} among (k-1)-sets */
                    const uint64_t *src = old + rs * (k - 1);
                    u128 acc = 0;
                    for (int i = 0; i < j; i++) acc += (u128)src[i] * wt[el[i]][w];
                    for (int i = j + 1; i < k; i++) acc += (u128)src[i - 1] * wt[el[i]][w];
                    dst[j] = (uint64_t)(acc % p);
                }
                r++;
                /* Gosper: next k-subset in increasing order */
                uint32_t cc = S & -S, rr = S + cc;
                S = (((rr ^ S) >> 2) / cc) | rr;
                if (rr == 0) break;
            }
            if (r != binom[n][k]) { fprintf(stderr, "rank count mismatch k=%d: %llu vs %llu\n", k, (unsigned long long)r, (unsigned long long)binom[n][k]); return 1; }
            uint64_t *t = old; old = cur; cur = t;
        }
        u128 tot = 0;
        for (int j = 0; j < n; j++) tot += old[j];
        printf("%llu %llu\n", (unsigned long long)x, (unsigned long long)(uint64_t)(tot % p));
        fflush(stdout);
    }
    free(old); free(cur); free(line);
    return 0;
}
