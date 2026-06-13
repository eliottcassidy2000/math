#!/usr/bin/env python3
"""
n=7 tiling isomorphism analysis using a C extension for canonicalization.
Falls back to pure Python if C compilation fails.

Key idea: use nauty-style certificate = lexicographically smallest adjacency
matrix under all permutations. Precompute all 5040 permutations.

Instance: opus-2026-03-06-S7
"""

import ctypes
import tempfile
import os
import subprocess
import time
import sys
from collections import defaultdict

C_CODE = r"""
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define N 7
#define NPERM 5040

static int perms[NPERM][N];
static int nperm = 0;

void gen_perms(int *arr, int l, int r) {
    if (l == r) {
        memcpy(perms[nperm++], arr, N*sizeof(int));
        return;
    }
    for (int i = l; i <= r; i++) {
        int tmp = arr[l]; arr[l] = arr[i]; arr[i] = tmp;
        gen_perms(arr, l+1, r);
        tmp = arr[l]; arr[l] = arr[i]; arr[i] = tmp;
    }
}

void init_perms() {
    int arr[N];
    for (int i = 0; i < N; i++) arr[i] = i;
    nperm = 0;
    gen_perms(arr, 0, N-1);
}

/* Compute canonical form and store in out. Returns score signature. */
void canonicalize(const int A[N][N], int out[N*N]) {
    int best_set = 0;
    for (int pi = 0; pi < nperm; pi++) {
        const int *p = perms[pi];
        int better = 0;
        for (int i = 0; i < N && !better; i++) {
            for (int j = 0; j < N && !better; j++) {
                int val = A[p[i]][p[j]];
                int idx = i*N + j;
                if (!best_set || val < out[idx]) {
                    better = 1;  /* new best */
                } else if (val > out[idx]) {
                    better = -1; /* worse */
                }
            }
        }
        if (better == 1 || !best_set) {
            best_set = 1;
            for (int i = 0; i < N; i++)
                for (int j = 0; j < N; j++)
                    out[i*N+j] = A[p[i]][p[j]];
        }
    }
}

/* Count automorphisms */
int count_aut(const int A[N][N]) {
    int count = 0;
    for (int pi = 0; pi < nperm; pi++) {
        const int *p = perms[pi];
        int ok = 1;
        for (int i = 0; i < N && ok; i++)
            for (int j = 0; j < N && ok; j++)
                if (A[p[i]][p[j]] != A[i][j]) ok = 0;
        if (ok) count++;
    }
    return count;
}

/* Count anti-automorphisms */
int count_anti_aut(const int A[N][N]) {
    int count = 0;
    for (int pi = 0; pi < nperm; pi++) {
        const int *p = perms[pi];
        int ok = 1;
        for (int i = 0; i < N && ok; i++)
            for (int j = 0; j < N && ok; j++)
                if (A[p[i]][p[j]] != A[j][i]) ok = 0;
        if (ok) count++;
    }
    return count;
}

/* Count Hamiltonian paths using bitmask DP */
long long count_ham_dp(const int A[N][N]) {
    /* dp[mask][v] */
    static long long dp[1<<N][N];
    memset(dp, 0, sizeof(dp));
    for (int v = 0; v < N; v++) dp[1<<v][v] = 1;
    for (int mask = 1; mask < (1<<N); mask++) {
        for (int v = 0; v < N; v++) {
            if (!dp[mask][v]) continue;
            if (!(mask & (1<<v))) continue;
            for (int u = 0; u < N; u++) {
                if (mask & (1<<u)) continue;
                if (A[v][u]) dp[mask|(1<<u)][u] += dp[mask][v];
            }
        }
    }
    long long total = 0;
    int full = (1<<N)-1;
    for (int v = 0; v < N; v++) total += dp[full][v];
    return total;
}

/* Process all 2^15 tilings */
/* Output: for each tiling, its canonical form hash (index into sorted list) */
/* tile_coords: array of (x,y) pairs, m tiles */
/* trans_map: the transpose permutation on tiles */

struct class_info {
    int canon[N*N];  /* canonical adjacency matrix */
    int count;       /* number of tilings */
    long long H;     /* Hamiltonian path count */
    int aut;         /* automorphism count */
    int anti_aut;    /* anti-automorphism count */
    int gs_count;    /* grid-symmetric tiling count */
    int is_sc;       /* self-converse? */
    int is_sf;       /* self-flip? */
    int scores[N];   /* sorted out-degrees */
};

#define MAX_CLASSES 1000
static struct class_info classes[MAX_CLASSES];
static int nclasses = 0;

/* Tiles for n=7 */
#define M 15
static int tile_x[M], tile_y[M];
static int trans_map_t[M];
static int verts[N]; /* [7,6,5,4,3,2,1] */

void setup_tiles() {
    int idx = 0;
    for (int y = 1; y <= N-2; y++)
        for (int x = N; x > y+1; x--)
        { tile_x[idx] = x; tile_y[idx] = y; idx++; }

    for (int i = 0; i < N; i++) verts[i] = N - i;

    /* Build trans_map */
    int tile_lookup[N+1][N+1];
    memset(tile_lookup, -1, sizeof(tile_lookup));
    for (int i = 0; i < M; i++) tile_lookup[tile_x[i]][tile_y[i]] = i;
    for (int i = 0; i < M; i++) {
        int nx = N - tile_y[i] + 1;
        int ny = N - tile_x[i] + 1;
        trans_map_t[i] = tile_lookup[nx][ny];
    }
}

int vert_idx(int v) {
    for (int i = 0; i < N; i++) if (verts[i] == v) return i;
    return -1;
}

void mask_to_adj(int mask, int A[N][N]) {
    memset(A, 0, sizeof(int)*N*N);
    for (int k = 0; k < N-1; k++) A[k][k+1] = 1;
    for (int i = 0; i < M; i++) {
        int xi = vert_idx(tile_x[i]);
        int yi = vert_idx(tile_y[i]);
        if ((mask >> i) & 1) A[yi][xi] = 1;
        else A[xi][yi] = 1;
    }
}

int is_grid_sym(int mask) {
    for (int i = 0; i < M; i++) {
        if (((mask>>i)&1) != ((mask>>trans_map_t[i])&1)) return 0;
    }
    return 1;
}

int find_class(int canon[N*N]) {
    for (int c = 0; c < nclasses; c++) {
        if (memcmp(classes[c].canon, canon, N*N*sizeof(int)) == 0) return c;
    }
    return -1;
}

int main() {
    init_perms();
    setup_tiles();

    printf("n=7, m=%d, 2^m=%d\n", M, 1<<M);
    printf("Processing %d tilings...\n", 1<<M);

    int tiling_class[1<<M];
    int gs_flags[1<<M];

    for (int mask = 0; mask < (1<<M); mask++) {
        int A[N][N];
        mask_to_adj(mask, A);

        int canon[N*N];
        canonicalize(A, canon);

        gs_flags[mask] = is_grid_sym(mask);

        int ci = find_class(canon);
        if (ci < 0) {
            ci = nclasses++;
            memcpy(classes[ci].canon, canon, N*N*sizeof(int));
            classes[ci].count = 0;
            classes[ci].H = count_ham_dp(A);
            classes[ci].aut = count_aut(A);
            classes[ci].anti_aut = count_anti_aut(A);
            classes[ci].gs_count = 0;
            classes[ci].is_sf = 0;

            /* Self-converse check */
            int A_op[N][N];
            for (int i = 0; i < N; i++)
                for (int j = 0; j < N; j++)
                    A_op[i][j] = A[j][i];
            int canon_op[N*N];
            canonicalize(A_op, canon_op);
            classes[ci].is_sc = (memcmp(canon, canon_op, N*N*sizeof(int)) == 0);

            /* Sorted scores */
            int degs[N];
            for (int i = 0; i < N; i++) {
                degs[i] = 0;
                for (int j = 0; j < N; j++) degs[i] += A[i][j];
            }
            /* Sort descending */
            for (int i = 0; i < N-1; i++)
                for (int j = i+1; j < N; j++)
                    if (degs[j] > degs[i]) { int t=degs[i]; degs[i]=degs[j]; degs[j]=t; }
            memcpy(classes[ci].scores, degs, N*sizeof(int));
        }
        tiling_class[mask] = ci;
        classes[ci].count++;
        if (gs_flags[mask]) classes[ci].gs_count++;

        if ((mask & 0x1FFF) == 0 && mask > 0) {
            fprintf(stderr, "  %d/%d (%.0f%%), %d classes so far\r",
                    mask, 1<<M, 100.0*mask/(1<<M), nclasses);
        }
    }
    fprintf(stderr, "\n");

    /* Self-flip check */
    for (int mask = 0; mask < (1<<M); mask++) {
        int flip = mask ^ ((1<<M)-1);
        if (tiling_class[mask] == tiling_class[flip]) {
            classes[tiling_class[mask]].is_sf = 1;
        }
    }

    printf("\n%d isomorphism classes\n\n", nclasses);

    /* Print results */
    int sc_count = 0, sf_count = 0, both_count = 0;
    long long max_h = 0;
    for (int c = 0; c < nclasses; c++) {
        if (classes[c].is_sc) sc_count++;
        if (classes[c].is_sf) sf_count++;
        if (classes[c].is_sc && classes[c].is_sf) both_count++;
        if (classes[c].H > max_h) max_h = classes[c].H;
    }

    printf("Self-converse: %d\n", sc_count);
    printf("Self-flip: %d\n", sf_count);
    printf("Both SC+SF: %d\n", both_count);
    printf("Max H: %lld\n\n", max_h);

    printf("--- SC+SF KERNEL CLASSES ---\n");
    for (int c = 0; c < nclasses; c++) {
        if (classes[c].is_sc && classes[c].is_sf) {
            printf("  #%d: H=%lld%s, |Aut|=%d, |AntiAut|=%d, GS=%d/%d (%.3f), scores=(",
                   c, classes[c].H, (classes[c].H == max_h ? "*MAX*" : ""),
                   classes[c].aut, classes[c].anti_aut,
                   classes[c].gs_count, classes[c].count,
                   (double)classes[c].gs_count/classes[c].count);
            for (int i = 0; i < N; i++) printf("%d%s", classes[c].scores[i], i<N-1?",":"");
            printf(")\n");
        }
    }

    printf("\n--- H-MAXIMIZERS (H=%lld) ---\n", max_h);
    for (int c = 0; c < nclasses; c++) {
        if (classes[c].H == max_h) {
            printf("  #%d: |Aut|=%d, |AntiAut|=%d, SC=%d, SF=%d, GS=%d/%d, scores=(",
                   c, classes[c].aut, classes[c].anti_aut,
                   classes[c].is_sc, classes[c].is_sf,
                   classes[c].gs_count, classes[c].count);
            for (int i = 0; i < N; i++) printf("%d%s", classes[c].scores[i], i<N-1?",":"");
            printf(")\n");
        }
    }

    printf("\n--- TOP 15 BY H ---\n");
    /* Sort classes by H descending */
    int order[MAX_CLASSES];
    for (int i = 0; i < nclasses; i++) order[i] = i;
    for (int i = 0; i < nclasses-1; i++)
        for (int j = i+1; j < nclasses; j++)
            if (classes[order[j]].H > classes[order[i]].H)
                { int t=order[i]; order[i]=order[j]; order[j]=t; }

    for (int k = 0; k < 15 && k < nclasses; k++) {
        int c = order[k];
        printf("  #%d: H=%lld, |Aut|=%d, |AntiAut|=%d, SC=%d, SF=%d, GS=%d/%d, scores=(",
               c, classes[c].H, classes[c].aut, classes[c].anti_aut,
               classes[c].is_sc, classes[c].is_sf,
               classes[c].gs_count, classes[c].count);
        for (int i = 0; i < N; i++) printf("%d%s", classes[c].scores[i], i<N-1?",":"");
        printf(")\n");
    }

    printf("\n--- GS DISTRIBUTION (top 15) ---\n");
    /* Sort by gs_count descending */
    for (int i = 0; i < nclasses-1; i++)
        for (int j = i+1; j < nclasses; j++)
            if (classes[order[j]].gs_count > classes[order[i]].gs_count)
                { int t=order[i]; order[i]=order[j]; order[j]=t; }

    int total_gs = 0;
    int gs_class_count = 0;
    for (int c = 0; c < nclasses; c++) {
        total_gs += classes[c].gs_count;
        if (classes[c].gs_count > 0) gs_class_count++;
    }
    printf("Total GS: %d, Classes with GS: %d/%d\n", total_gs, gs_class_count, nclasses);

    for (int k = 0; k < 15 && k < nclasses; k++) {
        int c = order[k];
        if (classes[c].gs_count == 0) break;
        printf("  #%d: GS=%d/%d (%.3f), H=%lld, SC=%d, SF=%d, |Aut|=%d\n",
               c, classes[c].gs_count, classes[c].count,
               (double)classes[c].gs_count/classes[c].count,
               classes[c].H, classes[c].is_sc, classes[c].is_sf, classes[c].aut);
    }

    /* SF-only classes */
    printf("\n--- SF-ONLY CLASSES ---\n");
    for (int c = 0; c < nclasses; c++) {
        if (classes[c].is_sf && !classes[c].is_sc) {
            printf("  #%d: H=%lld, |Aut|=%d, GS=%d/%d, scores=(",
                   c, classes[c].H, classes[c].aut,
                   classes[c].gs_count, classes[c].count);
            for (int i = 0; i < N; i++) printf("%d%s", classes[c].scores[i], i<N-1?",":"");
            printf(")\n");
        }
    }

    return 0;
}
"""

def run_c_analysis():
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

        print("Compiled successfully. Running...")
        t0 = time.time()
        result = subprocess.run([exe_path], capture_output=True, text=True, timeout=600)
        print(result.stdout)
        if result.stderr:
            print(f"(stderr: {result.stderr.strip()})")
        print(f"\nTotal time: {time.time()-t0:.1f}s")
    finally:
        os.unlink(c_path)
        if os.path.exists(exe_path):
            os.unlink(exe_path)

if __name__ == '__main__':
    run_c_analysis()
