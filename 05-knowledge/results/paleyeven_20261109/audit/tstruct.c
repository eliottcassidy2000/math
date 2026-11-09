/* tstruct.c - audit P1 (independent): module defect, twin distance, arc-strong connectivity, min semidegree.
 * Input lines: name n out_0 ... out_{n-1}   (out_i = decimal bitmask of out-neighbours of i), n <= 24.
 * md     = min over M, 2 <= |M| <= n-1, of sum_{z not in M} min(|N+(z) & M|, |N-(z) & M|)
 * twin   = the same minimum over |M| = 2
 * lambda = min over nonempty proper U of #arcs U -> V \ U
 * delta0 = min over v of min(d+(v), d-(v))
 * Also reports how many M attain md, and the minimisers when there are at most 40 of them.
 * Compile: cc -O2 -o tstruct tstruct.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

int main(void) {
    char name[64];
    int n;
    while (scanf("%63s %d", name, &n) == 2) {
        uint32_t out[32], in[32];
        memset(in, 0, sizeof in);
        for (int i = 0; i < n; i++) { unsigned long long m; if (scanf("%llu", &m) != 1) return 1; out[i] = (uint32_t)m; }
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (out[i] >> j & 1) in[j] |= 1u << i;
        /* tournament sanity */
        for (int i = 0; i < n; i++) {
            if (out[i] >> i & 1) { printf("%s: loop\n", name); return 1; }
            for (int j = i + 1; j < n; j++) if (((out[i] >> j) & 1) == ((out[j] >> i) & 1)) { printf("%s: not a tournament\n", name); return 1; }
        }
        uint32_t full = (n == 32) ? 0xffffffffu : ((1u << n) - 1);
        int delta0 = 99;
        for (int v = 0; v < n; v++) {
            int a = __builtin_popcount(out[v]), b = __builtin_popcount(in[v]);
            if (a < delta0) delta0 = a;
            if (b < delta0) delta0 = b;
        }
        int md = 1 << 30, twin = 1 << 30, lam = 1 << 30;
        long long nmd = 0;
        uint32_t mins[40]; int nm = 0;
        for (uint32_t M = 1; M < full; M++) {
            /* cut M -> complement */
            int cut = 0;
            uint32_t s = M;
            while (s) { int z = __builtin_ctz(s); s &= s - 1; cut += __builtin_popcount(out[z] & ~M & full); }
            if (cut < lam) lam = cut;
            int sz = __builtin_popcount(M);
            if (sz < 2) continue;
            int d = 0;
            uint32_t c = ~M & full;
            while (c) {
                int z = __builtin_ctz(c); c &= c - 1;
                int a = __builtin_popcount(out[z] & M), b = sz - a;   /* tournament: N-(z) & M = M \ N+(z) */
                d += (a < b) ? a : b;
            }
            if (sz == 2 && d < twin) twin = d;
            if (d < md) { md = d; nmd = 1; nm = 0; mins[nm++] = M; }
            else if (d == md) { nmd++; if (nm < 40) mins[nm++] = M; }
        }
        printf("%s n=%d md=%d twin=%d lambda=%d delta0=%d n_md_sets=%lld", name, n, md, twin, lam, delta0, nmd);
        if (nmd <= 40) {
            /* classify minimisers by size */
            int bysz[33] = {0};
            for (int i = 0; i < nm; i++) bysz[__builtin_popcount(mins[i])]++;
            printf(" sizes:");
            for (int s2 = 2; s2 < n; s2++) if (bysz[s2]) printf(" %dx%d", s2, bysz[s2]);
            if (nmd <= 40 && n <= 14) {
                printf(" pairs:");
                for (int i = 0; i < nm; i++) if (__builtin_popcount(mins[i]) == 2) {
                    int u = __builtin_ctz(mins[i]); int v = 31 - __builtin_clz(mins[i]);
                    printf(" {%d,%d}", u, v);
                }
            }
        }
        printf("\n");
        fflush(stdout);
    }
    return 0;
}
