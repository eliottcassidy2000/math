/*
 * h21_bottleneck_v2.c — FIXED: use unsigned short for n=9 bitmasks.
 * Verify bottleneck lemmas and analyze no-good-deletion tournaments.
 *
 * Author: opus-2026-03-07-S43
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define N 9

static unsigned int rng_state = 54321;
static inline unsigned int xorshift() {
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 17;
    rng_state ^= rng_state << 5;
    return rng_state;
}

static inline void random_tournament(unsigned short out[N]) {
    memset(out, 0, N * sizeof(unsigned short));
    for (int i = 0; i < N; i++)
        for (int j = i+1; j < N; j++) {
            if (xorshift() & 1)
                out[i] |= (1 << j);
            else
                out[j] |= (1 << i);
        }
}

static inline int has_source_sink(const unsigned short out[N]) {
    for (int i = 0; i < N; i++) {
        int s = __builtin_popcount(out[i]);
        if (s == 0 || s == N-1) return 1;
    }
    return 0;
}

static inline int vertex_in_3cycle(const unsigned short out[N], int v) {
    unsigned short all = (1 << N) - 1;
    unsigned short ov = out[v];
    unsigned short iv = (~ov) & all & ~(1 << v);
    unsigned short ov_copy = ov;
    while (ov_copy) {
        int a = __builtin_ctz(ov_copy);
        ov_copy &= ov_copy - 1;
        if (out[a] & iv & ~(1 << a))
            return 1;
    }
    return 0;
}

static inline int is_cycle_rich(const unsigned short out[N]) {
    if (has_source_sink(out)) return 0;
    for (int v = 0; v < N; v++)
        if (!vertex_in_3cycle(out, v)) return 0;
    return 1;
}

/* Check if u has a 3-cycle NOT through v */
static inline int has_3cycle_avoiding(const unsigned short out[N], int u, int v) {
    unsigned short all = (1 << N) - 1;
    unsigned short mask = all & ~(1 << v) & ~(1 << u);
    unsigned short ou = out[u] & mask;
    unsigned short iu = (~out[u]) & all & mask;

    unsigned short ou_copy = ou;
    while (ou_copy) {
        int a = __builtin_ctz(ou_copy);
        ou_copy &= ou_copy - 1;
        if (out[a] & iu)
            return 1;
    }
    return 0;
}

static inline int deletion_is_cycle_rich(const unsigned short out[N], int v) {
    unsigned short all = (1 << N) - 1;
    unsigned short mask = all & ~(1 << v);

    for (int u = 0; u < N; u++) {
        if (u == v) continue;
        int s = __builtin_popcount(out[u] & mask);
        if (s == 0 || s == N-2) return 0;
    }

    for (int u = 0; u < N; u++) {
        if (u == v) continue;
        if (!has_3cycle_avoiding(out, u, v)) return 0;
    }

    return 1;
}

static int count_3cycles(const unsigned short out[N]) {
    int count = 0;
    for (int a = 0; a < N; a++)
        for (int b = a+1; b < N; b++)
            for (int c = b+1; c < N; c++) {
                int ab = (out[a] >> b) & 1;
                int bc = (out[b] >> c) & 1;
                int ca = (out[c] >> a) & 1;
                if ((ab && bc && ca) || (!ab && !bc && !ca))
                    count++;
            }
    return count;
}

static int max_matching(const unsigned short out[N]) {
    int triples[100][3];
    int nt = 0;
    for (int a = 0; a < N; a++)
        for (int b = a+1; b < N; b++)
            for (int c = b+1; c < N; c++) {
                int ab = (out[a] >> b) & 1;
                int bc = (out[b] >> c) & 1;
                int ca = (out[c] >> a) & 1;
                if ((ab && bc && ca) || (!ab && !bc && !ca)) {
                    if (nt < 100) {
                        triples[nt][0] = a; triples[nt][1] = b; triples[nt][2] = c;
                        nt++;
                    }
                }
            }
    int used = 0, mm = 0;
    for (int i = 0; i < nt; i++) {
        int m = (1 << triples[i][0]) | (1 << triples[i][1]) | (1 << triples[i][2]);
        if (!(used & m)) { used |= m; mm++; }
    }
    return mm;
}

static int hamiltonian_paths(const unsigned short out[N]) {
    int dp[1 << N][N];
    memset(dp, 0, sizeof(dp));
    for (int v = 0; v < N; v++)
        dp[1 << v][v] = 1;
    int full = (1 << N) - 1;
    for (int mask = 1; mask <= full; mask++)
        for (int v = 0; v < N; v++) {
            if (!(mask & (1 << v)) || dp[mask][v] == 0) continue;
            for (int u = 0; u < N; u++) {
                if (mask & (1 << u)) continue;
                if ((out[v] >> u) & 1)
                    dp[mask | (1 << u)][u] += dp[mask][v];
            }
        }
    int total = 0;
    for (int v = 0; v < N; v++) total += dp[full][v];
    return total;
}

int main() {
    setvbuf(stdout, NULL, _IONBF, 0);
    setvbuf(stderr, NULL, _IONBF, 0);

    long long trials = 100000000LL;
    long long cycle_rich = 0;
    long long no_good_del = 0;
    int bottleneck_total_size[20]; /* total sum |B(v)| for no-good cases */
    int bottleneck_centers[20];   /* number of v with B(v) != empty */
    memset(bottleneck_total_size, 0, sizeof(bottleneck_total_size));
    memset(bottleneck_centers, 0, sizeof(bottleneck_centers));

    unsigned short out[N];

    for (long long t = 0; t < trials; t++) {
        random_tournament(out);
        if (!is_cycle_rich(out)) continue;
        cycle_rich++;

        int has_good = 0;
        for (int v = 0; v < N; v++) {
            if (deletion_is_cycle_rich(out, v)) { has_good = 1; break; }
        }
        if (has_good) continue;

        no_good_del++;
        int t3 = count_3cycles(out);
        int mm = max_matching(out);
        int H = hamiltonian_paths(out);

        /* Compute bottleneck structure */
        int B[N]; memset(B, 0, sizeof(B));
        int bn_of[N]; memset(bn_of, -1, sizeof(bn_of));

        for (int u = 0; u < N; u++) {
            for (int v = 0; v < N; v++) {
                if (v == u) continue;
                if (!has_3cycle_avoiding(out, u, v)) {
                    bn_of[u] = v;
                    B[v] |= (1 << u);
                    break;
                }
            }
        }

        int total_bn = 0, centers = 0;
        for (int v = 0; v < N; v++) {
            int sz = __builtin_popcount(B[v]);
            total_bn += sz;
            if (sz > 0) centers++;
        }

        if (total_bn < 20) bottleneck_total_size[total_bn]++;
        if (centers < 20) bottleneck_centers[centers]++;

        if (no_good_del <= 20) {
            int sc[N];
            for (int i = 0; i < N; i++) sc[i] = __builtin_popcount(out[i]);
            printf("\nNo-good #%lld: t3=%d mm=%d H=%d\n", no_good_del, t3, mm, H);
            printf("  Scores: ");
            for (int i = 0; i < N; i++) printf("%d ", sc[i]);
            printf("\n");
            for (int v = 0; v < N; v++) {
                if (B[v]) {
                    printf("  B(%d) = {", v);
                    for (int u = 0; u < N; u++)
                        if (B[v] & (1 << u)) printf("%d ", u);
                    printf("}\n");
                }
            }
            /* Also show score obstructions */
            unsigned short all = (1 << N) - 1;
            for (int v = 0; v < N; v++) {
                unsigned short mask = all & ~(1 << v);
                int obs = 0;
                for (int u = 0; u < N; u++) {
                    if (u == v) continue;
                    int s = __builtin_popcount(out[u] & mask);
                    if (s == 0) { printf("  Del %d: %d->source\n", v, u); obs = 1; }
                    if (s == N-2) { printf("  Del %d: %d->sink\n", v, u); obs = 1; }
                }
                if (!obs && !B[v] && !(B[0] & (1<<v)) && !(B[1] & (1<<v)) && !(B[2] & (1<<v)) &&
                    !(B[3] & (1<<v)) && !(B[4] & (1<<v)) && !(B[5] & (1<<v)) && !(B[6] & (1<<v)) &&
                    !(B[7] & (1<<v)) && !(B[8] & (1<<v))) {
                    /* No score obstruction and B(v)=empty, but still no good — why? */
                    printf("  Del %d: UNEXPLAINED FAILURE\n", v);
                }
            }
        }

        if ((t + 1) % 20000000 == 0 || no_good_del == 1) {
            fprintf(stderr, "t=%lldM, cr=%lld, no_good=%lld\n",
                    (t+1)/1000000, cycle_rich, no_good_del);
        }
    }

    printf("\n=== RESULTS ===\n");
    printf("Trials: %lld\nCycle-rich: %lld\nNo good deletion: %lld\n",
           trials, cycle_rich, no_good_del);
    printf("Total bottleneck size distribution:\n");
    for (int i = 0; i < 20; i++)
        if (bottleneck_total_size[i] > 0)
            printf("  sum|B(v)|=%d: %d cases\n", i, bottleneck_total_size[i]);
    printf("Number of bottleneck centers:\n");
    for (int i = 0; i < 20; i++)
        if (bottleneck_centers[i] > 0)
            printf("  #centers=%d: %d cases\n", i, bottleneck_centers[i]);

    return 0;
}
