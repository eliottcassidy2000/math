/*
 * h_spectrum_n9_sample.c — Sample H-spectrum at n=9.
 * Quick check: what H values appear? Is H=21 ever seen?
 *
 * Author: opus-2026-03-07-S42
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define N 9

int adj[N][N];

void random_tournament() {
    for (int i = 0; i < N; i++) {
        adj[i][i] = 0;
        for (int j = i+1; j < N; j++) {
            if (rand() & 1) {
                adj[i][j] = 1; adj[j][i] = 0;
            } else {
                adj[i][j] = 0; adj[j][i] = 1;
            }
        }
    }
}

long long hamiltonian_paths() {
    long long dp[512][9];
    memset(dp, 0, sizeof(dp));
    for (int v = 0; v < N; v++) dp[1<<v][v] = 1;
    for (int mask = 1; mask < 512; mask++)
        for (int v = 0; v < N; v++) {
            if (!(mask & (1<<v)) || dp[mask][v] == 0) continue;
            for (int u = 0; u < N; u++) {
                if (mask & (1<<u)) continue;
                if (adj[v][u])
                    dp[mask|(1<<u)][u] += dp[mask][v];
            }
        }
    long long total = 0;
    for (int v = 0; v < N; v++) total += dp[511][v];
    return total;
}

int main() {
    srand(time(NULL));
    setvbuf(stdout, NULL, _IONBF, 0);

    long long TRIALS = 2000000;
    /* Track which odd values in [1..200] are seen */
    int seen[201] = {0};

    for (long long trial = 0; trial < TRIALS; trial++) {
        random_tournament();
        long long H = hamiltonian_paths();
        if (H <= 200) seen[H] = 1;

        if ((trial & 0x1FFFF) == 0x1FFFF) {
            fprintf(stderr, "Trial %lld\n", trial);
        }
    }

    printf("H-spectrum at n=9 (%lld samples):\n", TRIALS);
    printf("Missing odd values in [1..200]: ");
    int count = 0;
    for (int h = 1; h <= 200; h += 2) {
        if (!seen[h]) {
            printf("%d ", h);
            count++;
        }
    }
    printf("\n(%d missing)\n", count);

    printf("\nPresent odd values in [1..50]: ");
    for (int h = 1; h <= 50; h += 2)
        if (seen[h]) printf("%d ", h);
    printf("\n");

    return 0;
}
