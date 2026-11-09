/* tprop.c - structural profile of tournaments (genesis lane, mac-mini-2026-11-09).
 *
 * Input lines:  n  code
 *   code = upper triangle in row-major order (pairs (i,j), i<j, i=0..n-2, j=i+1..n-1),
 *          character '1' means i -> j, '0' means j -> i.
 * Output line per tournament:
 *   n H strong lambda prime minsemideg | H(T-v) for v=0..n-1 | code
 *   strong  : 1 if strongly connected
 *   lambda  : arc-strong connectivity = min over proper nonempty S of #arcs S -> V\S
 *   prime   : 1 if T has no module M with 2 <= |M| <= n-1 (modules = intervals)
 * H is computed by the standard subset DP (exact, 64-bit).  n <= 22.
 * Compile: cc -O2 -o tprop tprop.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int n;
static uint32_t out[32];

static uint64_t hcount(uint32_t verts) {
    /* number of Hamiltonian paths of T[verts] */
    int idx[32], k = 0;
    for (int v = 0; v < n; v++) if (verts >> v & 1) idx[k++] = v;
    if (k == 0) return 1;
    if (k == 1) return 1;
    uint32_t lo[32];
    for (int a = 0; a < k; a++) {
        lo[a] = 0;
        for (int b = 0; b < k; b++) if (out[idx[a]] >> idx[b] & 1) lo[a] |= 1u << b;
    }
    size_t S = (size_t)1 << k;
    uint64_t *dp = calloc(S * k, sizeof(uint64_t));
    for (int a = 0; a < k; a++) dp[((size_t)1 << a) * k + a] = 1;
    for (size_t m = 1; m < S; m++) {
        for (int a = 0; a < k; a++) {
            uint64_t c = dp[m * k + a];
            if (!c) continue;
            uint32_t nx = lo[a] & ~(uint32_t)m;
            while (nx) {
                int b = __builtin_ctz(nx); nx &= nx - 1;
                dp[(m | (1u << b)) * k + b] += c;
            }
        }
    }
    uint64_t h = 0;
    for (int a = 0; a < k; a++) h += dp[(S - 1) * k + a];
    free(dp);
    return h;
}

static int lambda_arc(void) {
    uint32_t full = (n == 32) ? 0xffffffffu : ((1u << n) - 1);
    int best = 1 << 30;
    for (uint32_t S = 1; S < full; S++) {
        int c = 0;
        uint32_t r = S;
        while (r) { int v = __builtin_ctz(r); r &= r - 1; c += __builtin_popcount(out[v] & ~S & full); }
        if (c < best) best = c;
        if (best == 0) break;
    }
    return best;
}

static int is_prime(void) {
    uint32_t full = (1u << n) - 1;
    /* module check over all subsets of size 2..n-1 */
    for (uint32_t M = 1; M < full; M++) {
        int sz = __builtin_popcount(M);
        if (sz < 2) continue;
        int ok = 1;
        for (int w = 0; w < n && ok; w++) {
            if (M >> w & 1) continue;
            uint32_t o = out[w] & M;
            if (o != 0 && o != M) ok = 0;
        }
        if (ok) return 0;
    }
    return 1;
}

int main(void) {
    char code[1024];
    while (scanf("%d %1023s", &n, code) == 2) {
        memset(out, 0, sizeof out);
        int p = 0;
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++) {
                if (code[p] == '1') out[i] |= 1u << j; else out[j] |= 1u << i;
                p++;
            }
        uint32_t full = (1u << n) - 1;
        uint64_t H = hcount(full);
        int lam = lambda_arc();
        int strong = lam > 0;
        int pr = is_prime();
        int msd = 1 << 30;
        for (int v = 0; v < n; v++) {
            int o = __builtin_popcount(out[v]);
            int in = n - 1 - o;
            if (o < msd) msd = o;
            if (in < msd) msd = in;
        }
        printf("%d %llu %d %d %d %d |", n, (unsigned long long)H, strong, lam, pr, msd);
        for (int v = 0; v < n; v++) printf(" %llu", (unsigned long long)hcount(full & ~(1u << v)));
        printf(" | %s\n", code);
        fflush(stdout);
    }
    return 0;
}
