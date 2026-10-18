/* Count Petersen copies in G_q = Cay(F_(q^2), mu_(q+1)), q an odd prime (mac-mini-2026-10-18-unitpath).
   Same search as kuratowski_petersen_20261018/petersen_census.py: injective homomorphisms with vertex 0 -> 0 and
   vertex 1 -> 1 (complete up to the edge-transitive symmetry). F_(q^2) = F_q[t]/(t^2 - n), n a non-residue.
   Petersen: outer 0-1-2-3-4-0, spokes i-(i+5), inner 5-7-9-6-8-5.
   Usage: ./petersen_count QMIN QMAX  -> one line per prime: q, count, count/q^2 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static int q, nres, V, M;
static int *mux, *muy;            /* unit vectors */
static int *isunit;               /* isunit[a*q+b] = 1 if N(a + b t) = 1 */
static inline int add(int u, int v) { int a = u / q + v / q, b = u % q + v % q; if (a >= q) a -= q; if (b >= q) b -= q; return a * q + b; }
static inline int sub(int u, int v) { int a = u / q - v / q, b = u % q - v % q; if (a < 0) a += q; if (b < 0) b += q; return a * q + b; }
static inline int adj(int u, int v) { return isunit[sub(u, v)]; }
static int isprime(int n) { if (n < 2) return 0; for (int d = 2; d * d <= n; d++) if (n % d == 0) return 0; return 1; }
static long long powmod(long long a, long long e, long long m) { long long r = 1; a %= m; while (e) { if (e & 1) r = r * a % m; a = a * a % m; e >>= 1; } return r; }
int main(int argc, char **argv) {
    int qmin = atoi(argv[1]), qmax = atoi(argv[2]);
    for (q = qmin; q <= qmax; q++) {
        if (!isprime(q) || q == 2) continue;
        for (nres = 2; nres < q; nres++) if (powmod(nres, (q - 1) / 2, q) == q - 1) break;
        V = q * q;
        isunit = calloc(V, sizeof(int)); mux = malloc(sizeof(int) * (q + 1)); muy = malloc(sizeof(int) * (q + 1)); M = 0;
        for (int a = 0; a < q; a++) for (int b = 0; b < q; b++) {
            long long nn = ((long long)a * a - (long long)nres * b % q * b) % q; if (nn < 0) nn += q;
            if (nn == 1) { isunit[a * q + b] = 1; mux[M] = a; muy[M] = b; M++; }
        }
        int *mu = malloc(sizeof(int) * M); for (int i = 0; i < M; i++) mu[i] = mux[i] * q + muy[i];
        int z0 = 0, one = 1 * q + 0;
        /* CN[d]: common neighbours of 0 and d (at most two, by the two-circle lemma) */
        int *cn = malloc(sizeof(int) * 2 * V), *ncn = calloc(V, sizeof(int));
        for (int i = 0; i < M; i++) for (int k = 0; k < M; k++) {
            int d = sub(mu[i], mu[k]);           /* w = mu[i] is adjacent to 0 and to d iff w - d = mu[k] */
            if (ncn[d] < 2) cn[2 * d + ncn[d]] = mu[i];
            ncn[d]++;
        }
        for (int d = 1; d < V; d++) if (ncn[d] > 2) { fprintf(stderr, "two-circle lemma violated at q=%d\n", q); return 1; }
        long long count = 0;
        for (int i2 = 0; i2 < M; i2++) {
            int v2 = add(one, mu[i2]); if (v2 == z0) continue;
            for (int i3 = 0; i3 < M; i3++) {
                int v3 = add(v2, mu[i3]); if (v3 == z0 || v3 == one) continue;
                int d4 = v3;                                   /* v4 adjacent to v3 and to 0 */
                for (int a = 0; a < ncn[d4]; a++) {
                    int v4 = cn[2 * d4 + a]; if (v4 == one || v4 == v2) continue;
                    for (int i5 = 0; i5 < M; i5++) {
                        int v5 = mu[i5]; if (v5 == one || v5 == v2 || v5 == v3 || v5 == v4) continue;
                        int d7 = sub(v5, v2);
                        for (int b7 = 0; b7 < ncn[d7]; b7++) {
                            int v7 = add(v2, cn[2 * d7 + b7]);
                            if (v7 == z0 || v7 == one || v7 == v3 || v7 == v4 || v7 == v5) continue;
                            int d9 = sub(v7, v4);
                            for (int b9 = 0; b9 < ncn[d9]; b9++) {
                                int v9 = add(v4, cn[2 * d9 + b9]);
                                if (v9 == z0 || v9 == one || v9 == v2 || v9 == v3 || v9 == v5 || v9 == v7) continue;
                                int d6 = sub(v9, one);
                                for (int b6 = 0; b6 < ncn[d6]; b6++) {
                                    int v6 = add(one, cn[2 * d6 + b6]);
                                    if (v6 == z0 || v6 == v2 || v6 == v3 || v6 == v4 || v6 == v5 || v6 == v7 || v6 == v9) continue;
                                    int d8 = sub(v6, v3);
                                    for (int b8 = 0; b8 < ncn[d8]; b8++) {
                                        int v8 = add(v3, cn[2 * d8 + b8]);
                                        if (!adj(v8, v5)) continue;
                                        if (v8 == z0 || v8 == one || v8 == v2 || v8 == v4 || v8 == v5 || v8 == v6 || v8 == v7 || v8 == v9) continue;
                                        count++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        free(cn); free(ncn);
        printf("q = %4d  q mod 20 = %2d  q mod 60 = %2d  count = %8lld  count/q^2 = %.4f\n", q, q % 20, q % 60, count, (double)count / ((double)q * q));
        fflush(stdout);
        free(isunit); free(mux); free(muy); free(mu);
    }
    return 0;
}
