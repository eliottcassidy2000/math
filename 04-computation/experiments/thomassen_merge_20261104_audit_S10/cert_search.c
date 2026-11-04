/* Audit of THM-4870 (B): find, for every odd H <= LIM (except 3, 7, 21), an explicit tournament with
 * exactly H Hamiltonian paths and trivial automorphism group.  Random near-transitive tournaments of
 * order 4..12 (transitive order with a random set of reversed arcs).  Each hit is printed as
 *   H n row_0 ... row_{n-1}      (row_v = hex bitmask of out-neighbours of v)
 * and is re-verified independently by verify_certs.py.  Usage: cs LIM SAMPLES_PER_ORDER SEED
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef uint64_t u64; typedef uint32_t u32;
#define MAXN 12
static u64 DP[(1 << MAXN) * MAXN];
static u64 rng_s;
static u64 rng(void) { rng_s ^= rng_s << 13; rng_s ^= rng_s >> 7; rng_s ^= rng_s << 17; return rng_s; }
static u64 hpaths(int n, const u32 *out) {
    int full = (1 << n) - 1;
    memset(DP, 0, sizeof(u64) * (size_t)(1 << n) * n);
    for (int v = 0; v < n; v++) DP[(1 << v) * n + v] = 1;
    for (int m = 1; m <= full; m++) for (int v = 0; v < n; v++) {
        u64 c = DP[m * n + v]; if (!c) continue;
        u32 av = out[v] & ~(u32)m;
        while (av) { int w = __builtin_ctz(av); av &= av - 1; DP[(m | (1 << w)) * n + w] += c; }
    }
    u64 s = 0; for (int v = 0; v < n; v++) s += DP[full * n + v]; return s;
}
/* count automorphisms, stop at 2 */
static int mp[MAXN], cnt; static u32 used;
static void bt(int n, const u32 *out, int i) {
    if (cnt >= 2) return;
    if (i == n) { cnt++; return; }
    int d = __builtin_popcount(out[i]);
    for (int c = 0; c < n && cnt < 2; c++) {
        if (used >> c & 1) continue;
        if (__builtin_popcount(out[c]) != d) continue;
        int ok = 1;
        for (int j = 0; j < i; j++) if (((out[i] >> j) & 1) != ((out[c] >> mp[j]) & 1)) { ok = 0; break; }
        if (!ok) continue;
        mp[i] = c; used |= 1u << c; bt(n, out, i + 1); used &= ~(1u << c);
    }
}
static int asymmetric(int n, const u32 *out) { cnt = 0; used = 0; bt(n, out, 0); return cnt == 1; }
int main(int argc, char **argv) {
    int LIM = atoi(argv[1]); long S = atol(argv[2]); rng_s = strtoull(argv[3], 0, 10) | 1;
    char *found = calloc(LIM + 1, 1);
    int nfound = 0;
    for (int n = 4; n <= 12; n++) {
        for (long s = 0; s < S; s++) {
            u32 out[MAXN] = {0};
            for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) out[i] |= 1u << j; /* transitive */
            /* reversal probability: mixture, favouring few reversals at larger n */
            double p = (double)(rng() % 1000) / 1000.0 * (n <= 7 ? 0.6 : 3.0 / n);
            for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++)
                if ((double)(rng() % 1000000) / 1e6 < p) { out[i] &= ~(1u << j); out[j] |= 1u << i; }
            u64 H = hpaths(n, out);
            if (H > (u64)LIM || found[H]) continue;
            if (!asymmetric(n, out)) continue;
            found[H] = 1; nfound++;
            printf("%llu %d", (unsigned long long)H, n);
            for (int v = 0; v < n; v++) printf(" %x", out[v]);
            printf("\n");
        }
        fprintf(stderr, "after order %d: %d values found\n", n, nfound);
    }
    fprintf(stderr, "odd H <= %d without an asymmetric certificate:", LIM);
    for (int h = 1; h <= LIM; h += 2) if (!found[h]) fprintf(stderr, " %d", h);
    fprintf(stderr, "\n");
    return 0;
}
