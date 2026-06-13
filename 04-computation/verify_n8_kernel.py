#!/usr/bin/env python3
"""
Verify n=8 SC+SF kernel: precompute H for all masks, then targeted brute-force iso checks.
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

static int tile_x[M], tile_y[M], trans_map[M], verts[N];

void setup_tiles() {
    int idx = 0;
    for (int y = 1; y <= N-2; y++)
        for (int x = N; x > y+1; x--)
            { tile_x[idx] = x; tile_y[idx] = y; idx++; }
    for (int i = 0; i < N; i++) verts[i] = N - i;
    int tile_lookup[N+1][N+1];
    memset(tile_lookup, -1, sizeof(tile_lookup));
    for (int i = 0; i < M; i++) tile_lookup[tile_x[i]][tile_y[i]] = i;
    for (int i = 0; i < M; i++) {
        trans_map[i] = tile_lookup[N - tile_y[i] + 1][N - tile_x[i] + 1];
    }
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

int count_3cycles(const int A[N][N]) {
    int c = 0;
    for (int i = 0; i < N; i++)
        for (int j = i+1; j < N; j++)
            for (int k = j+1; k < N; k++) {
                if (A[i][j] && A[j][k] && A[k][i]) c++;
                if (A[i][k] && A[k][j] && A[j][i]) c++;
            }
    return c;
}

int is_grid_sym(int mask) {
    for (int i = 0; i < M; i++)
        if (((mask>>i)&1) != ((mask>>trans_map[i])&1)) return 0;
    return 1;
}

/* Brute force permutations */
static int perms[40320][N];
static int nperm = 0;

void gen_perms(int *arr, int l, int r) {
    if (l == r) { memcpy(perms[nperm++], arr, N*sizeof(int)); return; }
    for (int i = l; i <= r; i++) {
        int tmp = arr[l]; arr[l] = arr[i]; arr[i] = tmp;
        gen_perms(arr, l+1, r);
        tmp = arr[l]; arr[l] = arr[i]; arr[i] = tmp;
    }
}

int brute_is_iso(const int A[N][N], const int B[N][N]) {
    for (int pi = 0; pi < nperm; pi++) {
        const int *p = perms[pi];
        int ok = 1;
        for (int i = 0; i < N && ok; i++)
            for (int j = 0; j < N && ok; j++)
                if (A[p[i]][p[j]] != B[i][j]) ok = 0;
        if (ok) return 1;
    }
    return 0;
}

int brute_count_aut(const int A[N][N]) {
    int c = 0;
    for (int pi = 0; pi < nperm; pi++) {
        const int *p = perms[pi];
        int ok = 1;
        for (int i = 0; i < N && ok; i++)
            for (int j = 0; j < N && ok; j++)
                if (A[p[i]][p[j]] != A[i][j]) ok = 0;
        if (ok) c++;
    }
    return c;
}

int brute_count_anti_aut(const int A[N][N]) {
    int c = 0;
    for (int pi = 0; pi < nperm; pi++) {
        const int *p = perms[pi];
        int ok = 1;
        for (int i = 0; i < N && ok; i++)
            for (int j = 0; j < N && ok; j++)
                if (A[p[i]][p[j]] != A[j][i]) ok = 0;
        if (ok) c++;
    }
    return c;
}

/* Precomputed data */
static long long *H_all;    /* H value per mask */
static int *score_hash;     /* hash of sorted scores per mask */
static int *c3_all;         /* 3-cycle count per mask */

void compute_scores_sorted(const int A[N][N], int degs[N]) {
    for (int i = 0; i < N; i++) {
        degs[i] = 0;
        for (int j = 0; j < N; j++) degs[i] += A[i][j];
    }
    for (int i = 0; i < N-1; i++)
        for (int j = i+1; j < N; j++)
            if (degs[j] > degs[i]) { int t=degs[i]; degs[i]=degs[j]; degs[j]=t; }
}

int hash_scores(const int degs[N]) {
    int h = 0;
    for (int i = 0; i < N; i++) h = h * 10 + degs[i];
    return h;
}

int is_palindromic(int sh) {
    /* Decode score hash and check d[i]+d[N-1-i]=N-1 */
    int degs[N];
    int tmp = sh;
    for (int i = N-1; i >= 0; i--) { degs[i] = tmp % 10; tmp /= 10; }
    for (int i = 0; i < N/2; i++)
        if (degs[i] + degs[N-1-i] != N-1) return 0;
    return 1;
}

int main() {
    int arr[N];
    for (int i = 0; i < N; i++) arr[i] = i;
    nperm = 0;
    gen_perms(arr, 0, N-1);
    setup_tiles();

    int all_ones = NMASK - 1;

    /* Phase 1: Precompute H, scores, c3 for ALL masks */
    fprintf(stderr, "Phase 1: Precomputing H, scores, c3 for all %d masks...\n", NMASK);
    H_all = (long long *)malloc(NMASK * sizeof(long long));
    score_hash = (int *)malloc(NMASK * sizeof(int));
    c3_all = (int *)malloc(NMASK * sizeof(int));

    long long max_H = 0;
    int max_mask = -1;

    for (int mask = 0; mask < NMASK; mask++) {
        int A[N][N];
        mask_to_adj(mask, A);
        H_all[mask] = count_ham_dp(A);
        int degs[N];
        compute_scores_sorted(A, degs);
        score_hash[mask] = hash_scores(degs);
        c3_all[mask] = count_3cycles(A);

        if (H_all[mask] > max_H) { max_H = H_all[mask]; max_mask = mask; }
        if ((mask & 0xFFFF) == 0)
            fprintf(stderr, "  %d/%d (%.1f%%)\r", mask, NMASK, 100.0*mask/NMASK);
    }
    fprintf(stderr, "  Done. H-max = %lld at mask %d\n", max_H, max_mask);

    /* Phase 2: Analyze H-maximizer */
    printf("=== H-MAXIMIZER (mask %d) ===\n", max_mask);
    {
        int A[N][N], At[N][N], Af[N][N];
        mask_to_adj(max_mask, A);
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                At[i][j] = A[j][i];
        mask_to_adj(max_mask ^ all_ones, Af);

        printf("H = %lld\n", max_H);
        printf("|Aut| = %d\n", brute_count_aut(A));
        printf("|AntiAut| = %d\n", brute_count_anti_aut(A));
        printf("SC = %d\n", brute_is_iso(A, At));
        printf("SF = %d\n", brute_is_iso(A, Af));
        printf("H(flip) = %lld\n", H_all[max_mask ^ all_ones]);

        int degs[N];
        for (int i = 0; i < N; i++) {
            degs[i] = 0;
            for (int j = 0; j < N; j++) degs[i] += A[i][j];
        }
        printf("Scores (unsorted): ");
        for (int i = 0; i < N; i++) printf("%d ", degs[i]);
        printf("\nc3 = %d\n\n", c3_all[max_mask]);
    }

    /* Phase 3: Find ALL SC+SF masks using cheap filters + brute iso */
    fprintf(stderr, "Phase 3: Scanning for SC+SF tournaments...\n");

    /* Kernel classes: representatives and their properties */
    typedef struct {
        int mask;
        long long H;
    } Rep;
    Rep reps[200];
    int nreps = 0;

    int total_scsf_masks = 0;

    for (int mask = 0; mask < NMASK; mask++) {
        /* Filter 1: palindromic scores (necessary for SC) */
        if (!is_palindromic(score_hash[mask])) continue;

        /* Filter 2: same H as flip (necessary for SF) */
        int flip = mask ^ all_ones;
        if (H_all[mask] != H_all[flip]) continue;

        /* Filter 3: same scores as flip */
        if (score_hash[mask] != score_hash[flip]) continue;

        /* Filter 4: same c3 as flip */
        if (c3_all[mask] != c3_all[flip]) continue;

        /* Expensive check: is SF? */
        int A[N][N], Af[N][N];
        mask_to_adj(mask, A);
        mask_to_adj(flip, Af);
        if (!brute_is_iso(A, Af)) continue;

        /* Expensive check: is SC? */
        int At[N][N];
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                At[i][j] = A[j][i];
        if (!brute_is_iso(A, At)) continue;

        total_scsf_masks++;

        /* Check if isomorphic to existing representative */
        int found = 0;
        for (int r = 0; r < nreps; r++) {
            int B[N][N];
            mask_to_adj(reps[r].mask, B);
            if (brute_is_iso(A, B)) { found = 1; break; }
        }
        if (!found && nreps < 200) {
            reps[nreps].mask = mask;
            reps[nreps].H = H_all[mask];
            nreps++;
            fprintf(stderr, "  Found kernel class %d: mask=%d, H=%lld\n",
                    nreps, mask, H_all[mask]);
        }

        if ((mask & 0x3FFF) == 0)
            fprintf(stderr, "  %d/%d, %d classes, %d masks\r",
                    mask, NMASK, nreps, total_scsf_masks);
    }
    fprintf(stderr, "\n");

    /* Phase 4: Detailed analysis of each kernel class */
    printf("=== SC+SF KERNEL (n=8): %d iso classes from %d masks ===\n\n",
           nreps, total_scsf_masks);

    for (int r = 0; r < nreps; r++) {
        int A[N][N];
        mask_to_adj(reps[r].mask, A);
        int aut = brute_count_aut(A);
        int aaut = brute_count_anti_aut(A);

        /* Count class size and GS by scanning masks with matching H */
        int class_size = 0, gs_count = 0;
        for (int m2 = 0; m2 < NMASK; m2++) {
            if (H_all[m2] != reps[r].H) continue;
            if (score_hash[m2] != score_hash[reps[r].mask]) continue;
            if (c3_all[m2] != c3_all[reps[r].mask]) continue;
            int B[N][N];
            mask_to_adj(m2, B);
            if (brute_is_iso(A, B)) {
                class_size++;
                if (is_grid_sym(m2)) gs_count++;
            }
        }

        printf("Kernel class %d:\n", r+1);
        printf("  mask=%d, H=%lld\n", reps[r].mask, reps[r].H);
        printf("  |Aut|=%d, |AntiAut|=%d\n", aut, aaut);
        printf("  class_size=%d (expected H/|Aut|=%lld/%d=%lld)\n",
               class_size, reps[r].H, aut, reps[r].H / aut);
        printf("  GS=%d/%d (%.4f)\n", gs_count, class_size,
               (double)gs_count / class_size);

        int degs[N];
        for (int i = 0; i < N; i++) {
            degs[i] = 0;
            for (int j = 0; j < N; j++) degs[i] += A[i][j];
        }
        printf("  scores (unsorted): ");
        for (int i = 0; i < N; i++) printf("%d ", degs[i]);
        printf("\n  c3=%d\n\n", c3_all[reps[r].mask]);
    }

    /* Count H=661 iso classes with SC and/or SF */
    printf("=== H=661 CLASSES: SC/SF BREAKDOWN ===\n");
    Rep h661_reps[100];
    int nh661 = 0;

    for (int mask = 0; mask < NMASK; mask++) {
        if (H_all[mask] != 661) continue;

        int A[N][N];
        mask_to_adj(mask, A);

        /* Check if new iso class */
        int found = 0;
        for (int r = 0; r < nh661; r++) {
            int B[N][N];
            mask_to_adj(h661_reps[r].mask, B);
            if (brute_is_iso(A, B)) { found = 1; break; }
        }
        if (!found && nh661 < 100) {
            h661_reps[nh661].mask = mask;
            nh661++;
        }
    }

    printf("Total iso classes with H=661: %d\n", nh661);
    for (int r = 0; r < nh661; r++) {
        int A[N][N], At[N][N], Af[N][N];
        mask_to_adj(h661_reps[r].mask, A);
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                At[i][j] = A[j][i];
        mask_to_adj(h661_reps[r].mask ^ all_ones, Af);

        int sc = brute_is_iso(A, At);
        int sf = brute_is_iso(A, Af);
        int aut = brute_count_aut(A);

        printf("  H=661 class %d: mask=%d, SC=%d, SF=%d, |Aut|=%d\n",
               r+1, h661_reps[r].mask, sc, sf, aut);
    }

    free(H_all);
    free(score_hash);
    free(c3_all);
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

        print("Compiled. Running...")
        t0 = time.time()
        result = subprocess.run([exe_path], capture_output=True, text=True, timeout=7200)
        print(result.stdout)
        if result.stderr:
            lines = result.stderr.strip().split('\n')
            print(f"\n(stderr tail: {chr(10).join(lines[-5:])})")
        print(f"\nTotal time: {time.time()-t0:.1f}s")
    finally:
        os.unlink(c_path)
        if os.path.exists(exe_path):
            os.unlink(exe_path)

if __name__ == '__main__':
    run()
