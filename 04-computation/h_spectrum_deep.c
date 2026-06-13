/*
 * h_spectrum_deep.c — kind-pasteur-2026-03-21-S17c
 *
 * Exhaustive H-spectrum computation for n=3..8.
 * For each n: compute H(T) for all 2^{C(n,2)} tournaments.
 * Output: complete frequency table of H values.
 *
 * n=7: 2^21 = 2M tournaments (~10s)
 * n=8: 2^28 = 268M tournaments (~55min)
 *
 * Compile: gcc -O3 -o h_spectrum_deep h_spectrum_deep.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAX_N 9
#define MAX_H 1000

typedef unsigned long long u64;

static int adj[MAX_N][MAX_N];

static int ham_paths(int n) {
    int dp[1 << MAX_N][MAX_N];
    int full = (1 << n) - 1;
    memset(dp, 0, sizeof(dp));
    for (int v = 0; v < n; v++) dp[1 << v][v] = 1;
    for (int mask = 1; mask < (1 << n); mask++)
        for (int v = 0; v < n; v++) {
            if (!(mask & (1 << v))) continue;
            int c = dp[mask][v];
            if (c == 0) continue;
            for (int w = 0; w < n; w++) {
                if (mask & (1 << w)) continue;
                if (adj[v][w]) dp[mask | (1 << w)][w] += c;
            }
        }
    int total = 0;
    for (int v = 0; v < n; v++) total += dp[full][v];
    return total;
}

int main(int argc, char **argv) {
    int max_n = 7;
    if (argc > 1) max_n = atoi(argv[1]);

    for (int n = 3; n <= max_n; n++) {
        int m = n * (n-1) / 2;
        u64 total = 1ULL << m;

        /* H frequency table */
        int *freq = calloc(MAX_H + 1, sizeof(int));
        int max_h_seen = 0;

        /* Arc map */
        int arc_i[50], arc_j[50];
        int pos = 0;
        for (int i = 0; i < n; i++)
            for (int j = i+1; j < n; j++) {
                arc_i[pos] = i; arc_j[pos] = j; pos++;
            }

        time_t t0 = time(NULL);

        for (u64 bits = 0; bits < total; bits++) {
            memset(adj, 0, sizeof(adj));
            for (int p = 0; p < m; p++) {
                if (bits & (1ULL << p)) adj[arc_i[p]][arc_j[p]] = 1;
                else adj[arc_j[p]][arc_i[p]] = 1;
            }

            int H = ham_paths(n);
            if (H <= MAX_H) freq[H]++;
            if (H > max_h_seen) max_h_seen = H;

            if (n >= 7 && (bits & 0xFFFFF) == 0 && bits > 0) {
                fprintf(stderr, "  n=%d: %llu/%llu (%.1f%%)\r",
                        n, (unsigned long long)bits, (unsigned long long)total,
                        100.0*bits/total);
            }
        }

        time_t t1 = time(NULL);

        printf("n=%d (2^%d = %llu tournaments, %lds):\n", n, m,
               (unsigned long long)total, (long)(t1-t0));

        /* Print frequency table */
        int n_distinct = 0;
        int n_gaps = 0;
        for (int h = 1; h <= max_h_seen; h += 2) {
            if (freq[h] > 0) n_distinct++;
            else n_gaps++;
        }

        printf("  max_H = %d, distinct = %d, gaps = %d\n", max_h_seen, n_distinct, n_gaps);

        /* Print all H values and their frequencies */
        printf("  H | freq | freq/total\n");
        for (int h = 1; h <= max_h_seen; h += 2) {
            if (freq[h] > 0) {
                printf("  %3d | %10d | %.8f\n", h, freq[h], (double)freq[h]/total);
            }
        }

        /* Print gaps */
        printf("  Gaps: ");
        for (int h = 1; h <= max_h_seen; h += 2) {
            if (freq[h] == 0) printf("%d ", h);
        }
        printf("\n");

        /* H mod 4 distribution */
        int mod4[4] = {0};
        long long sum_mod4[4] = {0};
        for (int h = 1; h <= max_h_seen; h += 2) {
            mod4[h % 4] += (freq[h] > 0) ? 1 : 0;
            sum_mod4[h % 4] += freq[h];
        }
        printf("  H mod 4: {1: %d values, %lld tours} {3: %d values, %lld tours}\n",
               mod4[1], sum_mod4[1], mod4[3], sum_mod4[3]);

        /* H mod 8 distribution */
        printf("  H mod 8: ");
        for (int r = 1; r < 8; r += 2) {
            int cnt = 0;
            for (int h = r; h <= max_h_seen; h += 8)
                if (freq[h] > 0) cnt++;
            printf("{%d: %d} ", r, cnt);
        }
        printf("\n\n");

        free(freq);
        fflush(stdout);
    }
    return 0;
}
