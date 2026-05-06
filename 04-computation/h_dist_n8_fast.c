/*
 * H-distribution for n=8 via tiling enumeration.
 *
 * Formula: |{labeled T: H=h}| = 8!/h * |{tilings: H=h}|
 * Tilings = 2^{C(7,2)} = 2^21 = 2,097,152 configs.
 * Each H computation: bitmask DP O(2^8 * 8^2) = O(16384).
 *
 * Total: ~34 billion ops. In C: ~30 seconds.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define N 8
#define NTILES 21   /* C(7,2) = 21 tiles */
#define TILE_SPACE (1 << NTILES)  /* 2^21 = 2097152 */
#define FULL_MASK ((1 << N) - 1)
#define HDIST_SIZE 8192   /* H values up to 8191 for n=8 */

/* Tiles: (x,y) with y in {0,...,N-2}, x in {y+2,...,N-1} */
/* Listed in order: y=0: x=2,3,...,7; y=1: x=3,...,7; ... */
static int tile_x[NTILES], tile_y[NTILES];
static int ntiles_actual;

static int h_dist_tiling[HDIST_SIZE] = {0};
static long long h_dist_labeled[HDIST_SIZE] = {0};

/* H computation via bitmask DP */
static int hamiltonian_count(int adj[N][N]) {
    /* dp[mask][v] = paths using exactly vertices in mask, ending at v */
    static int dp[1 << N][N];
    int mask, v, u;

    memset(dp, 0, sizeof(dp));
    for (v = 0; v < N; v++)
        dp[1 << v][v] = 1;

    for (mask = 1; mask <= FULL_MASK; mask++) {
        for (v = 0; v < N; v++) {
            if (!(mask & (1 << v))) continue;
            if (!dp[mask][v]) continue;
            int cnt = dp[mask][v];
            for (u = 0; u < N; u++) {
                if (mask & (1 << u)) continue;
                if (adj[v][u])
                    dp[mask | (1 << u)][u] += cnt;
            }
        }
    }

    int total = 0;
    for (v = 0; v < N; v++)
        total += dp[FULL_MASK][v];
    return total;
}

int main(void) {
    int n_computed = 0;
    clock_t t_start, t_end;

    /* Build tile list */
    ntiles_actual = 0;
    for (int y = 0; y < N-1; y++) {
        for (int x = y+2; x < N; x++) {
            tile_x[ntiles_actual] = x;
            tile_y[ntiles_actual] = y;
            ntiles_actual++;
        }
    }
    if (ntiles_actual != NTILES) {
        fprintf(stderr, "Tile count mismatch: got %d, expected %d\n",
                ntiles_actual, NTILES);
        return 1;
    }

    printf("n=%d: %d tiles, %d tilings\n", N, NTILES, TILE_SPACE);
    printf("Enumerating all tilings...\n");
    fflush(stdout);

    t_start = clock();

    int adj[N][N];
    for (int tmask = 0; tmask < TILE_SPACE; tmask++) {
        /* Build adjacency matrix */
        memset(adj, 0, sizeof(adj));
        /* Base path: k -> k-1 for k=1..N-1 */
        for (int k = 1; k < N; k++)
            adj[k][k-1] = 1;
        /* Tiles */
        for (int bi = 0; bi < NTILES; bi++) {
            int x = tile_x[bi], y = tile_y[bi];
            if ((tmask >> bi) & 1)
                adj[y][x] = 1;  /* reversed: y→x */
            else
                adj[x][y] = 1;  /* default: x→y */
        }
        int h = hamiltonian_count(adj);
        if (h < HDIST_SIZE)
            h_dist_tiling[h]++;
        else
            fprintf(stderr, "WARNING: H=%d exceeds HDIST_SIZE at tmask=%d\n", h, tmask);

        if ((tmask & 0xFFFFF) == 0) {
            t_end = clock();
            double elapsed = (double)(t_end - t_start) / CLOCKS_PER_SEC;
            printf("  progress: %d/%d (%.1f%%) in %.1fs\n",
                   tmask, TILE_SPACE, 100.0*tmask/TILE_SPACE, elapsed);
            fflush(stdout);
        }
    }

    t_end = clock();
    double elapsed = (double)(t_end - t_start) / CLOCKS_PER_SEC;

    printf("\nDone in %.2fs\n\n", elapsed);

    /* Compute labeled counts */
    long long nfact = 40320LL;  /* 8! */
    long long total_tiling = 0, total_labeled = 0;
    int h_min = 0, h_max = 0;

    /* Find range */
    for (int h = 1; h < HDIST_SIZE; h++) {
        if (h_dist_tiling[h] > 0) {
            if (h_min == 0) h_min = h;
            h_max = h;
        }
    }

    printf("H-distribution for n=%d tilings:\n", N);
    printf("  (Tiling space = 2^%d = %d)\n", NTILES, TILE_SPACE);
    printf("\n");
    printf("  H  | #tilings | #labeled_T (= 8!/H * #tilings)\n");
    printf("  ---+----------+----------------------------------\n");

    /* Check forbidden values */
    printf("  Forbidden check: H=7 %s, H=21 %s\n",
           (h_dist_tiling[7] == 0 && 7 <= h_max) ? "ABSENT ✓" : (7 > h_max ? "above range" : "PRESENT!"),
           (h_dist_tiling[21] == 0 && 21 <= h_max) ? "ABSENT ✓" : (21 > h_max ? "above range" : "PRESENT!"));
    printf("\n");

    for (int h = h_min; h <= h_max; h++) {
        if (h_dist_tiling[h] == 0) continue;
        int tc = h_dist_tiling[h];
        long long lc = nfact * (long long)tc / h;
        if (nfact * (long long)tc % h != 0) {
            printf("  WARNING: non-integer labeled count at H=%d (count=%d)!\n", h, tc);
        }
        total_tiling += tc;
        total_labeled += lc;
        printf("  %4d | %8d | %lld\n", h, tc, lc);
    }

    printf("\n  Total tilings:  %lld (expected %d)\n",
           total_tiling, TILE_SPACE);
    printf("  Total labeled:  %lld (expected 2^%d = %lld)\n",
           total_labeled, N*(N-1)/2,
           (long long)1 << (N*(N-1)/2));

    /* ΣH over tilings */
    long long sigma_H = 0;
    for (int h = h_min; h <= h_max; h++)
        sigma_H += (long long)h * h_dist_tiling[h];
    printf("  ΣH (tilings):   %lld\n", sigma_H);
    printf("  H_max:          %d\n", h_max);
    printf("  Distinct H:     ");
    int ndist = 0;
    for (int h = h_min; h <= h_max; h++) if (h_dist_tiling[h]>0) ndist++;
    printf("%d values\n", ndist);

    /* Spectrum */
    printf("\nH-spectrum (achievable values):\n  ");
    int col = 0;
    for (int h = h_min; h <= h_max; h++) {
        if (h_dist_tiling[h] == 0) continue;
        printf("%d ", h);
        col++;
        if (col % 20 == 0) printf("\n  ");
    }
    printf("\n");

    /* Forbidden odd values in range */
    printf("\nForbidden odd values in [1, H_max]:\n  ");
    col = 0;
    for (int h = 1; h <= h_max; h += 2) {
        if (h_dist_tiling[h] == 0) {
            printf("%d ", h);
            col++;
            if (col % 20 == 0) printf("\n  ");
        }
    }
    printf("\n");

    return 0;
}
