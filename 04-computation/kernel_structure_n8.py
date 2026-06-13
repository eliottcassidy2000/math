#!/usr/bin/env python3
"""
Detailed structural analysis of the 5 SC+SF kernel classes at n=8.
Uses the verified masks from verify_n8_kernel.py.
Instance: opus-2026-03-06-S8
"""

import subprocess, tempfile, os, time

C_CODE = r"""
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define N 8
#define M 21
#define NMASK (1 << M)

static int tile_x[M], tile_y[M], trans_map[M];

void setup_tiles() {
    int idx = 0;
    for (int y = 1; y <= N-2; y++)
        for (int x = N; x > y+1; x--)
            { tile_x[idx] = x; tile_y[idx] = y; idx++; }
    int tile_lookup[N+1][N+1];
    memset(tile_lookup, -1, sizeof(tile_lookup));
    for (int i = 0; i < M; i++) tile_lookup[tile_x[i]][tile_y[i]] = i;
    for (int i = 0; i < M; i++)
        trans_map[i] = tile_lookup[N - tile_y[i] + 1][N - tile_x[i] + 1];
}

void mask_to_adj(int mask, int A[N][N]) {
    memset(A, 0, sizeof(int)*N*N);
    for (int k = 0; k < N-1; k++) A[k][k+1] = 1;
    for (int i = 0; i < M; i++) {
        int xi = N - tile_x[i];
        int yi = N - tile_y[i];
        if ((mask >> i) & 1) A[yi][xi] = 1;
        else A[xi][yi] = 1;
    }
}

long long count_ham_dp(const int A[N][N]) {
    static long long dp[1<<N][N];
    memset(dp, 0, sizeof(dp));
    for (int v = 0; v < N; v++) dp[1<<v][v] = 1;
    for (int mask = 1; mask < (1<<N); mask++)
        for (int v = 0; v < N; v++) {
            if (!dp[mask][v] || !(mask & (1<<v))) continue;
            for (int u = 0; u < N; u++) {
                if (mask & (1<<u)) continue;
                if (A[v][u]) dp[mask|(1<<u)][u] += dp[mask][v];
            }
        }
    long long total = 0;
    for (int v = 0; v < N; v++) total += dp[(1<<N)-1][v];
    return total;
}

int count_directed_kcycles(const int A[N][N], int k) {
    /* Count directed k-cycles using DP on bitmask */
    /* A k-cycle visits exactly k distinct vertices and returns to start */
    if (k < 3 || k > N) return 0;
    int count = 0;
    /* For each starting vertex s, count paths of length k-1 starting at s
       that end at some vertex with edge back to s */
    for (int s = 0; s < N; s++) {
        long long dp[1<<N][N];
        memset(dp, 0, sizeof(dp));
        dp[1<<s][s] = 1;
        for (int step = 0; step < k-1; step++) {
            long long ndp[1<<N][N];
            memset(ndp, 0, sizeof(ndp));
            for (int mask = 0; mask < (1<<N); mask++)
                for (int v = 0; v < N; v++) {
                    if (!dp[mask][v]) continue;
                    for (int u = 0; u < N; u++) {
                        if (u == s && step < k-2) continue; /* don't return early */
                        if (mask & (1<<u)) continue;
                        if (A[v][u]) ndp[mask|(1<<u)][u] += dp[mask][v];
                    }
                }
            memcpy(dp, ndp, sizeof(dp));
        }
        /* Count paths ending at vertex with edge to s */
        for (int mask = 0; mask < (1<<N); mask++)
            for (int v = 0; v < N; v++) {
                if (!dp[mask][v]) continue;
                if (A[v][s]) count += dp[mask][v];
            }
    }
    /* Each k-cycle is counted k times (once for each starting vertex) */
    return count / k;
}

/* I(Omega(T), 2): independence polynomial of odd cycle conflict graph at x=2 */
/* Omega(T) has vertices = directed odd cycles, edges = sharing a vertex */
/* I(Omega, 2) = sum over independent sets S of 2^|S| */
/* For n=8 this could have many cycles, let's just use OCF: H = I(Omega, 2) */

int main() {
    setup_tiles();

    /* Verified kernel masks from brute-force analysis */
    int kernel_masks[] = {125298, 307609, 334741, 450778, 477910};
    long long kernel_H[] = {641, 653, 637, 621, 657};
    int kernel_aut[] = {1, 1, 1, 3, 1};
    int nkernel = 5;

    /* Also the H-maximizer */
    int hmax_mask = 6603;

    printf("=== DETAILED KERNEL CLASS STRUCTURE (n=8) ===\n\n");

    for (int k = 0; k < nkernel; k++) {
        int A[N][N];
        mask_to_adj(kernel_masks[k], A);

        printf("--- Kernel class %d (mask=%d, H=%lld, |Aut|=%d) ---\n",
               k+1, kernel_masks[k], kernel_H[k], kernel_aut[k]);

        /* Print adjacency matrix */
        printf("Adjacency matrix:\n");
        for (int i = 0; i < N; i++) {
            printf("  ");
            for (int j = 0; j < N; j++) printf("%d", A[i][j]);
            printf("\n");
        }

        /* Out-degrees (unsorted) */
        printf("Out-degrees: ");
        for (int i = 0; i < N; i++) {
            int d = 0;
            for (int j = 0; j < N; j++) d += A[i][j];
            printf("%d ", d);
        }
        printf("\n");

        /* Directed cycle counts */
        int c3 = count_directed_kcycles(A, 3);
        int c5 = count_directed_kcycles(A, 5);
        int c7 = count_directed_kcycles(A, 7);
        printf("Directed cycles: c3=%d, c5=%d, c7=%d\n", c3, c5, c7);

        /* Check OCF: H should equal I(Omega, 2) */
        /* We know H from verification */
        /* alpha_1 = c3+c5+c7 (total odd cycles) */
        int alpha1 = c3 + c5 + c7;
        printf("alpha_1 (total odd cycles) = %d\n", alpha1);

        /* Count vertex-disjoint 3-cycle pairs (alpha_2) */
        /* Brute force: for each pair of 3-cycles, check vertex-disjoint */
        int triples[100][3];
        int ntriples = 0;
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++) {
                if (j == i) continue;
                for (int kk = 0; kk < N; kk++) {
                    if (kk == i || kk == j) continue;
                    if (A[i][j] && A[j][kk] && A[kk][i]) {
                        /* Directed 3-cycle i->j->k->i. Store canonical (min first) */
                        if (i < j && i < kk) {
                            triples[ntriples][0] = i;
                            triples[ntriples][1] = j;
                            triples[ntriples][2] = kk;
                            ntriples++;
                        }
                    }
                }
            }
        int alpha2 = 0;
        for (int a = 0; a < ntriples; a++)
            for (int b = a+1; b < ntriples; b++) {
                int disjoint = 1;
                for (int x = 0; x < 3 && disjoint; x++)
                    for (int y = 0; y < 3 && disjoint; y++)
                        if (triples[a][x] == triples[b][y]) disjoint = 0;
                if (disjoint) alpha2++;
            }
        printf("alpha_2 (disjoint 3-cycle pairs) = %d\n", alpha2);

        /* OCF check: H = 1 + 2*alpha_1 + 4*alpha_2 + ... */
        /* For n=8, there can be higher terms (triples of disjoint 3-cycles, 5+3 pairs, etc.) */
        long long ocf_approx = 1 + 2*alpha1 + 4*alpha2;
        printf("OCF approx (1 + 2*a1 + 4*a2) = %lld (actual H=%lld, diff=%lld)\n",
               ocf_approx, kernel_H[k], kernel_H[k] - ocf_approx);

        /* The difference H - (1+2a1+4a2) accounts for:
           - vertex-disjoint pairs involving 5-cycles or 7-cycles
           - triples of disjoint 3-cycles (contributing 8)
           - etc. */

        /* Compare with H-maximizer */
        printf("\n");
    }

    /* H-maximizer for comparison */
    printf("--- H-MAXIMIZER (mask=%d, H=661) ---\n", hmax_mask);
    {
        int A[N][N];
        mask_to_adj(hmax_mask, A);
        printf("Adjacency matrix:\n");
        for (int i = 0; i < N; i++) {
            printf("  ");
            for (int j = 0; j < N; j++) printf("%d", A[i][j]);
            printf("\n");
        }
        printf("Out-degrees: ");
        for (int i = 0; i < N; i++) {
            int d = 0;
            for (int j = 0; j < N; j++) d += A[i][j];
            printf("%d ", d);
        }
        printf("\n");

        int c3 = count_directed_kcycles(A, 3);
        int c5 = count_directed_kcycles(A, 5);
        int c7 = count_directed_kcycles(A, 7);
        printf("Directed cycles: c3=%d, c5=%d, c7=%d\n", c3, c5, c7);
        printf("alpha_1 = %d\n", c3+c5+c7);

        int triples[100][3];
        int ntriples = 0;
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++) {
                if (j == i) continue;
                for (int kk = 0; kk < N; kk++) {
                    if (kk == i || kk == j) continue;
                    if (A[i][j] && A[j][kk] && A[kk][i]) {
                        if (i < j && i < kk) {
                            triples[ntriples][0] = i;
                            triples[ntriples][1] = j;
                            triples[ntriples][2] = kk;
                            ntriples++;
                        }
                    }
                }
            }
        int alpha2 = 0;
        for (int a = 0; a < ntriples; a++)
            for (int b = a+1; b < ntriples; b++) {
                int disjoint = 1;
                for (int x = 0; x < 3 && disjoint; x++)
                    for (int y = 0; y < 3 && disjoint; y++)
                        if (triples[a][x] == triples[b][y]) disjoint = 0;
                if (disjoint) alpha2++;
            }
        printf("alpha_2 = %d\n", alpha2);
        long long ocf_approx = 1 + 2*(c3+c5+c7) + 4*alpha2;
        printf("OCF approx = %lld (actual H=661, diff=%lld)\n",
               ocf_approx, 661 - ocf_approx);
    }

    /* Summary table */
    printf("\n=== SUMMARY TABLE ===\n");
    printf("Class | H   | c3 | c5 | c7 | a1  | a2 | 1+2a1+4a2 | residual\n");
    printf("------+-----+----+----+----+-----+----+-----------+---------\n");

    for (int k = 0; k < nkernel; k++) {
        int A[N][N];
        mask_to_adj(kernel_masks[k], A);
        int c3 = count_directed_kcycles(A, 3);
        int c5 = count_directed_kcycles(A, 5);
        int c7 = count_directed_kcycles(A, 7);

        int triples[100][3]; int ntriples = 0;
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++) {
                if (j == i) continue;
                for (int kk = 0; kk < N; kk++) {
                    if (kk == i || kk == j) continue;
                    if (A[i][j] && A[j][kk] && A[kk][i] && i<j && i<kk) {
                        triples[ntriples][0]=i; triples[ntriples][1]=j;
                        triples[ntriples][2]=kk; ntriples++;
                    }
                }
            }
        int a2 = 0;
        for (int a = 0; a < ntriples; a++)
            for (int b = a+1; b < ntriples; b++) {
                int d=1;
                for (int x=0;x<3&&d;x++) for (int y=0;y<3&&d;y++)
                    if(triples[a][x]==triples[b][y]) d=0;
                if(d) a2++;
            }
        int a1 = c3+c5+c7;
        long long base = 1+2*a1+4*a2;
        printf("K%d    | %3lld | %2d | %2d | %2d | %3d | %2d | %9lld | %3lld\n",
               k+1, kernel_H[k], c3, c5, c7, a1, a2, base, kernel_H[k]-base);
    }

    /* H-max */
    {
        int A[N][N]; mask_to_adj(hmax_mask, A);
        int c3 = count_directed_kcycles(A, 3);
        int c5 = count_directed_kcycles(A, 5);
        int c7 = count_directed_kcycles(A, 7);
        int triples[100][3]; int ntriples = 0;
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++) {
                if (j == i) continue;
                for (int kk = 0; kk < N; kk++) {
                    if (kk == i || kk == j) continue;
                    if (A[i][j] && A[j][kk] && A[kk][i] && i<j && i<kk) {
                        triples[ntriples][0]=i; triples[ntriples][1]=j;
                        triples[ntriples][2]=kk; ntriples++;
                    }
                }
            }
        int a2 = 0;
        for (int a = 0; a < ntriples; a++)
            for (int b = a+1; b < ntriples; b++) {
                int d=1;
                for (int x=0;x<3&&d;x++) for (int y=0;y<3&&d;y++)
                    if(triples[a][x]==triples[b][y]) d=0;
                if(d) a2++;
            }
        int a1 = c3+c5+c7;
        long long base = 1+2*a1+4*a2;
        printf("Hmax  | 661 | %2d | %2d | %2d | %3d | %2d | %9lld | %3lld\n",
               c3, c5, c7, a1, a2, base, 661-base);
    }

    return 0;
}
"""

def run():
    with tempfile.NamedTemporaryFile(mode='w', suffix='.c', delete=False) as f:
        f.write(C_CODE)
        c_path = f.name
    exe_path = c_path.replace('.c', '')
    try:
        result = subprocess.run(['gcc', '-O2', '-o', exe_path, c_path],
                              capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Compilation failed: {result.stderr}")
            return
        t0 = time.time()
        result = subprocess.run([exe_path], capture_output=True, text=True, timeout=600)
        print(result.stdout)
        if result.stderr:
            print(f"(stderr: {result.stderr.strip()[-200:]})")
        print(f"Time: {time.time()-t0:.1f}s")
    finally:
        os.unlink(c_path)
        if os.path.exists(exe_path): os.unlink(exe_path)

if __name__ == '__main__':
    run()
