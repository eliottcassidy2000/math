/* mdefect.c - exact module defect, twin distance, arc-strong connectivity and semidegrees (paleyeven lane,
 * mac-mini-2026-11-09).
 *
 * Input lines: "n code [tag]"  (code = upper triangle row-major, '1' iff i -> j), n <= 26.
 * For every vertex set M with 2 <= |M| <= n-1:
 *     def(M) = sum_{z not in M} min(|N+(z) & M|, |N-(z) & M|).
 * md = min def(M) (= THM-4985's decomposability defect for n >= 3); twin = min over |M| = 2.
 * lambda = min over nonempty proper S of #arcs S -> V\S (exact, all subsets).
 * Output: tag n md twin lambda minsemideg nmin minsize_profile(md by |M|) example_argmin(bitmask)
 * Compile: cc -O2 -o mdefect mdefect.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int n;
static uint32_t out[32], in_[32];

int main(void) {
    char code[512], tag[256];
    char line[1024];
    while (fgets(line, sizeof line, stdin)) {
        tag[0] = 0;
        if (sscanf(line, "%d %511s %255s", &n, code, tag) < 2) continue;
        if (n > 26 || (int)strlen(code) != n * (n - 1) / 2) { fprintf(stderr, "bad line\n"); return 1; }
        memset(out, 0, sizeof out); memset(in_, 0, sizeof in_);
        int k = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
            if (code[k++] == '1') { out[i] |= 1u << j; in_[j] |= 1u << i; }
            else { out[j] |= 1u << i; in_[i] |= 1u << j; }
        }
        int msd = n;
        for (int v = 0; v < n; v++) {
            int a = __builtin_popcount(out[v]), b = __builtin_popcount(in_[v]);
            int m = a < b ? a : b; if (m < msd) msd = m;
        }
        uint32_t full = (n == 32) ? 0xffffffffu : ((1u << n) - 1);
        int md = 1 << 30, twin = 1 << 30, lam = 1 << 30;
        long long nmin = 0; uint32_t arg = 0;
        int bysize[33]; for (int s = 0; s <= 32; s++) bysize[s] = 1 << 30;
        for (uint32_t M = 1; M < full; M++) {
            /* lambda: arcs from M to complement */
            int cut = 0;
            for (int z = 0; z < n; z++) if (M >> z & 1) cut += __builtin_popcount(out[z] & ~M & full);
            if (cut < lam) lam = cut;
            int sz = __builtin_popcount(M);
            if (sz < 2) continue;
            int d = 0;
            uint32_t C = ~M & full;
            while (C) {
                int z = __builtin_ctz(C); C &= C - 1;
                int a = __builtin_popcount(out[z] & M), b = __builtin_popcount(in_[z] & M);
                d += a < b ? a : b;
            }
            if (d < bysize[sz]) bysize[sz] = d;
            if (sz == 2 && d < twin) twin = d;
            if (d < md) { md = d; nmin = 1; arg = M; }
            else if (d == md) nmin++;
        }
        printf("%s n=%d md=%d twin=%d lambda=%d minsemideg=%d nmin=%lld bysize=", tag[0] ? tag : "-", n, md, twin, lam,
               msd, nmin);
        for (int s = 2; s <= n - 1; s++) printf("%d%s", bysize[s], s < n - 1 ? "," : "");
        printf(" argmin=0x%x\n", arg);
        fflush(stdout);
    }
    return 0;
}
