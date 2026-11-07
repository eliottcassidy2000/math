/* Independent audit of THM-4890 (D): the randomised-micro model.

   n >= 12, m = floor(n/2) - 2, C = 2^m, P = 2^(n-m-1), B = P/C.  Words v = h 2^m + l, h in [P, 2P), l in [0, C).
   c = C_m(l) (reversed complement), e(l) = [l < 2^(m-1)] = [c odd].
   g_F(v) = g(h) + h c + F(l) + [h odd] e(l), F(l) i.i.d. uniform on {0..C-1}; F = g - 1 recovers g.

   Computes EXACTLY (no sampling):
     type 1 (diagonal) = P C;
     type 2 (h = h', l != l') minimum separation (must be > 0);
     type 3 (h != h', l = l') = deterministic count;
     type 4 (h != h', l != l') expectation = sum over quadruples of P(F(l') - F(l) = Delta) = max(0, C - |Delta|)/C^2;
     the deterministic (DC) count (type 4 with F = g - 1), and E_n = PC + type3 + DC, cross-checked against a sort;
     the strip totals sum_(h != h') Z(h, h') with W = C (Z = #{(c,c') in [0,C)^2 : |h c - h' c' - D| <= C});
   the same expectation for F i.i.d. with the EMPIRICAL law of g(l) - 1 on [0, C) (rho(t) = sum_N p(N) p(N+t));
   plus Monte Carlo realisations of the random model (energy of g_F by sorting) for a few seeds.
   Usage: aud_randmodel n [mc_seeds]
*/
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

static uint32_t rev_w(uint32_t x, int w) { uint32_t r = 0; for (int i = 0; i < w; i++) { r = (r << 1) | (x & 1); x >>= 1; } return r; }
static uint32_t Cw(uint32_t x, int w) { return ((1u << w) - 1) - rev_w(x, w); }

static int cmp_u64(const void *a, const void *b) { uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b; return (x > y) - (x < y); }

static uint64_t rng_state = 88172645463325252ull;
static inline uint64_t xorshift64(void) { rng_state ^= rng_state << 13; rng_state ^= rng_state >> 7; rng_state ^= rng_state << 17; return rng_state; }

static long long energy_of(uint64_t *vals, size_t N) {
    qsort(vals, N, sizeof(uint64_t), cmp_u64);
    long long E = 0; size_t i = 0;
    while (i < N) { size_t k = i + 1; while (k < N && vals[k] == vals[i]) k++; E += (long long)(k - i) * (k - i); i = k; }
    return E;
}

int main(int argc, char **argv) {
    int n = atoi(argv[1]);
    int seeds = argc > 2 ? atoi(argv[2]) : 0;
    int m = n / 2 - 2;
    long long C = 1ll << m, P = 1ll << (n - m - 1), B = P / C;
    uint32_t *G = (uint32_t *)malloc(sizeof(uint32_t) << n);
    G[0] = 1;
    for (uint64_t v = 1; v < (1ull << n); v++) { uint64_t u = v >> 1; G[v] = (v & 1) ? G[u] : (uint32_t)(G[u] + u + (u & 1)); }
    long long *cl = malloc(sizeof(long long) * C), *lc = malloc(sizeof(long long) * C);
    for (long long l = 0; l < C; l++) { cl[l] = Cw((uint32_t)l, m); lc[cl[l]] = l; }
    /* identity check: g(h 2^m + l) = g(h) + h c + (g(l) - 1) + [h odd][l < 2^(m-1)] and 0 <= g(l) - 1 <= C - 1 */
    long long bad_id = 0, bad_F = 0;
    for (long long l = 0; l < C; l++) if (!(G[l] >= 1 && G[l] - 1 <= C - 1)) bad_F++;
    for (long long h = P; h < 2 * P; h++) for (long long l = 0; l < C; l++) {
        long long v = h * C + l;
        long long rhs = (long long)G[h] + h * cl[l] + ((long long)G[l] - 1) + ((h & 1) && (l < C / 2) ? 1 : 0);
        if (rhs != (long long)G[v]) bad_id++;
    }
    /* type 2: minimum over h and c != c' of h|c - c'| - (C - 1) - 1 */
    long long type2_min_sep = P * 1 - (C - 1) - 1;
    /* types 3, 4 */
    long long type3 = 0, dc = 0, Ztot = 0;
    double type4 = 0.0, type4emp = 0.0;
    /* empirical-law model: F(l) i.i.d. with the law of g(l) - 1 on [0, C) */
    double *pN = calloc(C, sizeof(double)), *rho = calloc(2 * C, sizeof(double));
    for (long long l = 0; l < C; l++) pN[G[l] - 1] += 1.0 / C;
    for (long long t = -(C - 1); t <= C - 1; t++) { double s = 0; for (long long N = 0; N < C; N++) { long long M2 = N + t; if (M2 >= 0 && M2 < C) s += pN[N] * pN[M2]; } rho[t + C] = s; }
    for (long long h = P; h < 2 * P; h++) {
        long long gh = G[h], ho = h & 1;
        for (long long hp = P; hp < 2 * P; hp++) {
            if (hp == h) continue;
            long long ghp = G[hp], hpo = hp & 1;
            long long D = ghp - gh;
            /* type 3: (h - h') c = D + e([h' odd] - [h odd]), with e = [c odd] */
            for (int e = 0; e < 2; e++) {
                long long num = D + e * (hpo - ho), den = h - hp;
                if (num % den == 0) { long long c = num / den; if (c >= 0 && c < C && (c & 1) == e) type3++; }
            }
            for (long long c = 0; c < C; c++) {
                long long K = gh - ghp + h * c + ho * (c & 1);       /* Delta = K - h' c' - [h' odd][c' odd] */
                long long c0 = (K - C) / hp - 2; if (c0 < 0) c0 = 0;
                long long c1 = (K + C) / hp + 2; if (c1 > C - 1) c1 = C - 1;
                for (long long cp = c0; cp <= c1; cp++) {
                    long long hcD = h * c - hp * cp - D;               /* strip variable */
                    if (hcD >= -C && hcD <= C) Ztot++;
                    if (cp == c) continue;                             /* l = l' handled by type 3 */
                    long long Delta = K - hp * cp - hpo * (cp & 1);
                    long long ad = Delta < 0 ? -Delta : Delta;
                    if (ad <= C - 1) {
                        type4 += (double)(C - ad) / (double)(C * C);
                        type4emp += rho[Delta + C];
                        long long l = lc[c], lp = lc[cp];
                        if ((long long)G[lp] - (long long)G[l] == Delta) dc++;
                    }
                }
            }
        }
    }
    double PC = (double)(P * C);
    /* deterministic energy by sorting, as a cross-check of PC + type3 + dc */
    size_t N = (size_t)(P * C);
    uint64_t *vals = malloc(sizeof(uint64_t) * N);
    for (size_t z = 0; z < N; z++) vals[z] = G[(1ull << (n - 1)) + z];
    long long En = energy_of(vals, N);
    printf("n=%d m=%d C=%lld P=%lld B=%lld |X_n|=PC=%lld\n", n, m, C, P, B, P * C);
    printf("  identity g = g(h)+hc+(g(l)-1)+[h odd][l<2^(m-1)] violations %lld; 0<=g(l)-1<=C-1 violations %lld\n", bad_id, bad_F);
    printf("  type 2 minimum separation P-C (must be >0): %lld\n", type2_min_sep);
    printf("  type 3 (h!=h', l=l') count %lld = %.4f PC   (bound 2P^2 = %.1f PC)\n", type3, type3 / PC, 2.0 * P * P / PC);
    printf("  type 4 expectation (random F) %.1f = %.4f PC   (bound 76 PC)\n", type4, type4 / PC);
    printf("  strip total sum Z(h,h') (W = C) %lld; /C = %.4f PC\n", Ztot, Ztot / (double)C / PC);
    printf("  empirical-law model (F i.i.d. with the law of g(l)-1 on [0,C)): type 4 expectation = %.4f PC; DC / uniform-model = %.3f; DC / empirical-model = %.3f\n",
           type4emp / PC, dc / type4, dc / type4emp);
    printf("  random-model expected energy E_F = PC + type3 + type4 = %.4f |X_n|   (claimed bound 109 |X_n|)\n", 1.0 + type3 / PC + type4 / PC);
    printf("  deterministic: DC count %lld = %.4f PC; PC + type3 + DC = %.4f |X_n|; sorted energy E_n = %.4f |X_n| (match %s)\n",
           dc, dc / PC, 1.0 + type3 / PC + dc / PC, En / PC, (En == P * C + type3 + dc) ? "yes" : "NO");
    for (int sd = 1; sd <= seeds; sd++) {
        rng_state = 0x9E3779B97F4A7C15ull * (uint64_t)(sd + 1000 * n);
        long long *F = malloc(sizeof(long long) * C);
        for (long long l = 0; l < C; l++) F[l] = (long long)(xorshift64() % (uint64_t)C);
        for (long long h = P, z = 0; h < 2 * P; h++) for (long long l = 0; l < C; l++, z++)
            vals[z] = (uint64_t)((long long)G[h] + h * cl[l] + F[l] + ((h & 1) && (l < C / 2) ? 1 : 0));
        long long EF = energy_of(vals, N);
        printf("  Monte Carlo seed %d: energy of g_F = %.4f |X_n|\n", sd, EF / PC);
        free(F);
    }
    fflush(stdout);
    return 0;
}
