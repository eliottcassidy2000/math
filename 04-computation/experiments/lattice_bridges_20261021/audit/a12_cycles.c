/* a12_cycles.c -- AUDIT pass 2 (independent of circ_counts.c and alpha_decomp.c): directed cycle counts of a circulant
 * tournament Cay(Z_n, S).
 *
 *   mode hk  : Held-Karp from vertex 0 over subsets of Z_n \ {0}; dp[mask][v] = # paths 0 -> ... -> v with vertex set
 *              {0} u mask.  By rotation invariance c_k = (n/k) * #(k-cycles through 0), and #(k-cycles through 0) =
 *              sum_{|mask| = k-1} sum_{v -> 0} dp[mask][v].  Also H = n * sum_v dp[full][v].  Optional: dump
 *              g(Z) = # directed Hamiltonian cycles of T[Z] for every Z containing 0 (index = mask of Z \ {0}) to a
 *              binary file of uint64 (for a13_alpha.py).
 *   mode iec : inclusion-exclusion for the Hamiltonian cycles only: C = sum_{X containing 0} (-1)^(n-|X|) (A_X^n)_{00}
 *              (closed walks from 0 of length n inside X), modulo 2^128 (exact).
 * stdin: lines "n s_1 ... s_m"; stdout: one line per set.   Usage: ./a12_cycles hk [dumpprefix] | iec
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned __int128 u128;
static void p128(u128 x){ char b[64]; int i = 63; b[i] = 0; if (!x){ printf("0"); return; } while (x){ b[--i] = (char)('0' + (int)(x % 10)); x /= 10; } printf("%s", b + i); }
static int n; static uint32_t outm[33];
static void build(const int *S, int ms){
    int in[64] = {0};
    for (int t = 0; t < ms; t++) in[((S[t] % n) + n) % n] = 1;
    for (int i = 0; i < n; i++){ outm[i] = 0; for (int d = 1; d < n; d++) if (in[d]) outm[i] |= 1u << ((i + d) % n); }
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++)
        if ((((outm[i] >> j) & 1) + ((outm[j] >> i) & 1)) != 1){ fprintf(stderr, "not a tournament\n"); exit(2); }
}
int main(int argc, char **argv){
    if (argc < 2){ fprintf(stderr, "usage\n"); return 1; }
    int hk = !strcmp(argv[1], "hk"); const char *dump = (hk && argc > 2) ? argv[2] : NULL;
    char line[4096];
    while (fgets(line, sizeof line, stdin)){
        int vals[64], k = 0; char *p = line, *e;
        for (;;){ long x = strtol(p, &e, 10); if (e == p) break; vals[k++] = (int)x; p = e; }
        if (k < 1) continue;
        n = vals[0]; build(vals + 1, k - 1);
        for (int t = 0; t < k; t++) printf(t ? " %d" : "%d", vals[t]);
        if (hk){
            int N1 = n - 1; size_t full = (size_t)1 << N1;
            uint64_t *f = calloc(full * (size_t)N1, sizeof(uint64_t)); if (!f){ fprintf(stderr, "oom\n"); return 3; }
            uint64_t *g = dump ? calloc(full, sizeof(uint64_t)) : NULL;
            uint32_t o[32], into0 = 0;
            for (int v = 1; v < n; v++){ o[v - 1] = outm[v] >> 1; if ((outm[v] >> 0) & 1) into0 |= 1u << (v - 1); }
            uint32_t o0 = outm[0] >> 1;
            for (int b = 0; b < N1; b++) if ((o0 >> b) & 1) f[((size_t)1 << b) * N1 + b] = 1;
            u128 thru0[34]; memset(thru0, 0, sizeof thru0);
            for (size_t mask = 1; mask < full; mask++){
                uint64_t *row = f + mask * N1; uint32_t mm = (uint32_t)mask; u128 closing = 0;
                while (mm){ int v = __builtin_ctz(mm); mm &= mm - 1; uint64_t c = row[v]; if (!c) continue;
                    if ((into0 >> v) & 1) closing += c;
                    uint32_t nx = o[v] & ~(uint32_t)mask;
                    while (nx){ int w = __builtin_ctz(nx); nx &= nx - 1; f[(mask | ((size_t)1 << w)) * N1 + w] += c; } }
                int sz = __builtin_popcount((uint32_t)mask) + 1;     /* cycle length */
                if (sz >= 3) thru0[sz] += closing;
                if (g) g[mask] = (uint64_t)closing;
            }
            u128 P0 = 0; uint64_t *last = f + (full - 1) * N1; for (int v = 0; v < N1; v++) P0 += last[v];
            free(f);
            printf(" : H="); p128(P0 * (u128)n);
            for (int L = 3; L <= n; L++){ u128 c = thru0[L] * (u128)n; if (c % (u128)L){ fprintf(stderr, "non-integral c_%d\n", L); return 4; }
                printf(" c%d=", L); p128(c / (u128)L); }
            if (g){ char fn[512]; snprintf(fn, sizeof fn, "%s_n%d.bin", dump, n); FILE *fp = fopen(fn, "wb"); fwrite(g, sizeof(uint64_t), full, fp); fclose(fp); free(g); }
        } else {
            /* closed walks from 0 of length n inside X, X ranging over subsets containing 0 */
            u128 acc = 0; uint32_t restfull = ((n == 32) ? 0xFFFFFFFFu : ((1u << n) - 1)) & ~1u;
            for (uint32_t r = 0; ; r = (r - restfull) & restfull){
                uint32_t X = r | 1u; int sz = __builtin_popcount(X);
                u128 v[32], w[32]; for (int i = 0; i < n; i++) v[i] = 0; v[0] = 1;
                for (int step = 0; step < n; step++){
                    for (int i = 0; i < n; i++) w[i] = 0;
                    for (uint32_t xs = X; xs; xs &= xs - 1){ int i = __builtin_ctz(xs); if (!v[i]) continue;
                        uint32_t nx = outm[i] & X; while (nx){ int j = __builtin_ctz(nx); nx &= nx - 1; w[j] += v[i]; } }
                    for (int i = 0; i < n; i++) v[i] = w[i];
                }
                if ((n - sz) & 1) acc -= v[0]; else acc += v[0];
                if (r == restfull) break;
            }
            printf(" : C(iec)="); p128(acc);
        }
        printf("\n"); fflush(stdout);
    }
    return 0;
}
