/*
 * h21_exhaustive_n8_v3.c — Ultra-fast exhaustive H=21 check at n=8
 *
 * Key optimizations:
 * 1. Precompute out-neighbor bitmasks per vertex (fast score check)
 * 2. Use popcount for scores
 * 3. Early exit on source/sink
 * 4. Bit-parallel 3-cycle counting
 * 5. Only run Held-Karp on tournaments passing all filters
 *
 * Author: opus-2026-03-07-S41
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define N 8
#define NBITS (N*(N-1)/2)
#define TOTAL (1U << NBITS)

/* Precompute out-neighbor masks from tournament encoding */
static inline void compute_out_masks(unsigned int T, unsigned char out[N]) {
    memset(out, 0, N);
    int pos = 0;
    for (int i = 0; i < N; i++) {
        for (int j = i+1; j < N; j++) {
            if ((T >> pos) & 1) {
                out[i] |= (1 << j);  /* i -> j */
            } else {
                out[j] |= (1 << i);  /* j -> i */
            }
            pos++;
        }
    }
}

/* Check source/sink using popcount on out-masks */
static inline int has_source_sink(const unsigned char out[N]) {
    for (int i = 0; i < N; i++) {
        int s = __builtin_popcount(out[i]);
        if (s == 0 || s == N-1) return 1;
    }
    return 0;
}

/* Count 3-cycles using bitmasks:
 * {a,b,c} forms a 3-cycle iff either a->b, b->c, c->a or a->c, c->b, b->a.
 * Equivalently, the subtournament on {a,b,c} has all scores = 1.
 * Score of a in {a,b,c} = (a->b) + (a->c). 3-cycle iff all three scores are 1.
 *
 * Faster: t3 = C(n,3) - sum_v C(s_v, 2) where s_v = score of v.
 * This is because t3 = (C(n,3) - sum_v C(s_v,2)) / ... no that's not right.
 *
 * Actually: sum of C(s_v, 2) over all v counts ordered pairs (j,k) with j<k
 * both in out-neighborhood of v. Each such triple {v,j,k} where v beats both
 * is counted once. The number of "transitive triples" is sum_v C(s_v, 2).
 * And t3 = C(n,3) - sum_v C(s_v, 2).
 */
static inline int count_3cycles(const unsigned char out[N]) {
    int sum_choose2 = 0;
    for (int i = 0; i < N; i++) {
        int s = __builtin_popcount(out[i]);
        sum_choose2 += s * (s - 1) / 2;
    }
    return 56 - sum_choose2;  /* C(8,3) = 56 */
}

/* Check if all vertices are in at least one 3-cycle.
 * Vertex v is in a 3-cycle iff there exist a, b with v->a->b->v or v->b->a->v.
 * v->a means a in out[v]. a->b means b in out[a]. b->v means v in out[b],
 * i.e., !(v in out[v]) which is always false... let me think again.
 * b->v means (out[b] >> v) & 1.
 *
 * More efficient: v is in a 3-cycle iff out[v] & in[v] overlap when considering
 * paths of length 2. Specifically, there exists a with v->a (a in out[v]) and
 * a->b->v. That means b in in[v] (b->v) and a->b (b in out[a]).
 * So we need: exists a in out[v] such that out[a] & in[v] has a bit set (not v, not a).
 * in[v] = ~out[v] & ((1<<N)-1) & ~(1<<v).
 */
static inline int all_in_3cycle(const unsigned char out[N]) {
    unsigned char in_mask[N];
    unsigned char all = (1 << N) - 1;
    for (int v = 0; v < N; v++)
        in_mask[v] = (~out[v]) & all & ~(1 << v);

    for (int v = 0; v < N; v++) {
        int found = 0;
        /* v is in a 3-cycle if exists a in out[v] with out[a] & in_mask[v] & ~(1<<a) != 0 */
        unsigned char ov = out[v];
        unsigned char iv = in_mask[v];
        while (ov && !found) {
            int a = __builtin_ctz(ov);
            ov &= ov - 1;  /* clear lowest bit */
            if (out[a] & iv & ~(1 << a))
                found = 1;
        }
        if (!found) return 0;
    }
    return 1;
}

/* Held-Karp for Hamiltonian path count */
static long long hamiltonian_paths(const unsigned char out[N]) {
    static long long dp[1 << N][N];
    memset(dp, 0, sizeof(dp));
    for (int v = 0; v < N; v++)
        dp[1 << v][v] = 1;
    for (unsigned int mask = 1; mask < (1U << N); mask++) {
        for (int v = 0; v < N; v++) {
            if (!(mask & (1 << v))) continue;
            if (dp[mask][v] == 0) continue;
            /* Extend to neighbors of v */
            unsigned char targets = out[v] & ~mask;
            unsigned char t = targets;
            while (t) {
                int u = __builtin_ctz(t);
                t &= t - 1;
                dp[mask | (1 << u)][u] += dp[mask][v];
            }
        }
    }
    long long total = 0;
    unsigned int full = (1U << N) - 1;
    for (int v = 0; v < N; v++) total += dp[full][v];
    return total;
}

int main() {
    setvbuf(stderr, NULL, _IONBF, 0);
    setvbuf(stdout, NULL, _IONBF, 0);
    long long h21_count = 0;
    long long skipped_src_sink = 0;
    long long skipped_high_t3 = 0;
    long long skipped_free_vertex = 0;
    long long checked_hk = 0;
    unsigned char out[N];

    for (unsigned int T = 0; T < TOTAL; T++) {
        compute_out_masks(T, out);

        /* Filter 1: source/sink */
        if (has_source_sink(out)) {
            skipped_src_sink++;
            continue;
        }

        /* Filter 2: count 3-cycles via score formula */
        int t3 = count_3cycles(out);
        if (t3 > 10) {
            skipped_high_t3++;
            continue;
        }

        /* Filter 3: all vertices in some 3-cycle */
        if (!all_in_3cycle(out)) {
            skipped_free_vertex++;
            continue;
        }

        /* Held-Karp */
        checked_hk++;
        long long H = hamiltonian_paths(out);

        if (H == 21) {
            h21_count++;
            printf("H=21 FOUND! T=%u H=%lld t3=%d\n", T, H, t3);
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    if (i == j) printf(".");
                    else printf("%d", (out[i] >> j) & 1);
                }
                printf("\n");
            }
            fflush(stdout);
        }

        if ((T & 0xFFFFFF) == 0xFFFFFF) {
            fprintf(stderr, "Progress: %u/%u (%.1f%%). src=%lld t3=%lld free=%lld hk=%lld h21=%lld\n",
                    T + 1, TOTAL, 100.0 * (T + 1) / TOTAL,
                    skipped_src_sink, skipped_high_t3, skipped_free_vertex, checked_hk, h21_count);
        }
    }

    printf("\n=== EXHAUSTIVE n=8 RESULTS ===\n");
    printf("Total tournaments: %u\n", TOTAL);
    printf("Skipped (source/sink): %lld\n", skipped_src_sink);
    printf("Skipped (t3 > 10): %lld\n", skipped_high_t3);
    printf("Skipped (free vertex): %lld\n", skipped_free_vertex);
    printf("Checked via Held-Karp: %lld\n", checked_hk);
    printf("H=21 found: %lld\n", h21_count);

    return 0;
}
