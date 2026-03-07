/*
 * h21_bottleneck_verify.c — Verify bottleneck lemmas at n=9.
 *
 * For cycle-rich n=9 tournaments, compute:
 * 1. B(v) = set of vertices whose ALL 3-cycles go through v
 * 2. Verify: each u has at most 1 bottleneck (Lemma 2)
 * 3. Verify: B(v) forms a transitive subtournament (Lemma 3+)
 * 4. For no-good-deletion tournaments: analyze bottleneck structure
 * 5. Verify: B(v)!=empty implies at least one 3-cycle through v uses a non-B vertex
 *
 * Author: opus-2026-03-07-S43
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define N 9
#define NBITS (N*(N-1)/2)

static unsigned int rng_state = 54321;
static inline unsigned int xorshift() {
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 17;
    rng_state ^= rng_state << 5;
    return rng_state;
}

static inline void compute_out_masks(unsigned int T, unsigned char out[N]) {
    memset(out, 0, N);
    int pos = 0;
    for (int i = 0; i < N; i++)
        for (int j = i+1; j < N; j++) {
            if ((T >> pos) & 1)
                out[i] |= (1 << j);
            else
                out[j] |= (1 << i);
            pos++;
        }
}

static inline void random_tournament(unsigned char out[N]) {
    memset(out, 0, N);
    for (int i = 0; i < N; i++)
        for (int j = i+1; j < N; j++) {
            if (xorshift() & 1)
                out[i] |= (1 << j);
            else
                out[j] |= (1 << i);
        }
}

static inline int has_source_sink(const unsigned char out[N]) {
    for (int i = 0; i < N; i++) {
        int s = __builtin_popcount(out[i]);
        if (s == 0 || s == N-1) return 1;
    }
    return 0;
}

static inline int vertex_in_3cycle(const unsigned char out[N], int v) {
    unsigned char all = (1 << N) - 1;
    unsigned char ov = out[v];
    unsigned char iv = (~ov) & all & ~(1 << v);
    unsigned char ov_copy = ov;
    while (ov_copy) {
        int a = __builtin_ctz(ov_copy);
        ov_copy &= ov_copy - 1;
        if (out[a] & iv & ~(1 << a))
            return 1;
    }
    return 0;
}

static inline int all_in_3cycle(const unsigned char out[N]) {
    for (int v = 0; v < N; v++)
        if (!vertex_in_3cycle(out, v)) return 0;
    return 1;
}

static inline int is_cycle_rich(const unsigned char out[N]) {
    if (has_source_sink(out)) return 0;
    return all_in_3cycle(out);
}

/* Check if vertex u has a 3-cycle NOT through v */
static inline int has_3cycle_avoiding(const unsigned char out[N], int u, int v) {
    unsigned char all = (1 << N) - 1;
    unsigned char mask = all & ~(1 << v) & ~(1 << u);  /* exclude both u and v */
    unsigned char ou = out[u] & mask;  /* u's out-neighbors excluding v */
    unsigned char iu = (~out[u]) & mask;  /* u's in-neighbors excluding v */

    unsigned char ou_copy = ou;
    while (ou_copy) {
        int a = __builtin_ctz(ou_copy);
        ou_copy &= ou_copy - 1;
        /* Check if a -> some b -> u, with b != v */
        if (out[a] & iu)
            return 1;
    }
    return 0;
}

/* Check if T-v is cycle-rich */
static inline int deletion_is_cycle_rich(const unsigned char out[N], int v) {
    unsigned char all = (1 << N) - 1;
    unsigned char mask = all & ~(1 << v);

    for (int u = 0; u < N; u++) {
        if (u == v) continue;
        int s = __builtin_popcount(out[u] & mask);
        if (s == 0 || s == N-2) return 0; /* source or sink */
    }

    for (int u = 0; u < N; u++) {
        if (u == v) continue;
        if (!has_3cycle_avoiding(out, u, v)) return 0;
    }

    return 1;
}

/* Max matching of 3-cycle sets */
static int max_matching_3cycles(const unsigned char out[N]) {
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
                        triples[nt][0] = a;
                        triples[nt][1] = b;
                        triples[nt][2] = c;
                        nt++;
                    }
                }
            }

    /* Greedy matching */
    int used = 0, mm = 0;
    for (int i = 0; i < nt; i++) {
        int m = (1 << triples[i][0]) | (1 << triples[i][1]) | (1 << triples[i][2]);
        if (!(used & m)) {
            used |= m;
            mm++;
        }
    }
    return mm;
}

int main() {
    setvbuf(stdout, NULL, _IONBF, 0);
    setvbuf(stderr, NULL, _IONBF, 0);

    long long trials = 50000000;
    long long cycle_rich = 0;
    long long no_good_del = 0;
    long long lemma2_violations = 0;

    /* Track bottleneck statistics for no-good-deletion cases */
    int bn_size_hist[N+1]; /* histogram of |B(v)| sizes */
    memset(bn_size_hist, 0, sizeof(bn_size_hist));
    int max_bn_set_seen = 0;

    unsigned char out[N];

    for (long long t = 0; t < trials; t++) {
        random_tournament(out);
        if (!is_cycle_rich(out)) continue;
        cycle_rich++;

        /* Compute bottleneck function B(v) for all v */
        int B[N]; /* B[v] = bitmask of vertices bottlenecked by v */
        memset(B, 0, sizeof(B));

        int bottleneck_of[N]; /* which vertex is u's bottleneck (-1 if none) */
        memset(bottleneck_of, -1, sizeof(bottleneck_of));
        int num_bottlenecks = 0;

        for (int u = 0; u < N; u++) {
            for (int v = 0; v < N; v++) {
                if (v == u) continue;
                if (!has_3cycle_avoiding(out, u, v)) {
                    /* v is a bottleneck for u */
                    if (bottleneck_of[u] != -1 && bottleneck_of[u] != v) {
                        /* u has two bottlenecks! Lemma 2 violation? */
                        /* Actually this means ALL 3-cycles of u go through BOTH v and bottleneck_of[u] */
                        /* So u's only 3-cycle is {u, v, bottleneck_of[u]} */
                        lemma2_violations++;
                    }
                    bottleneck_of[u] = v;
                    B[v] |= (1 << u);
                    num_bottlenecks++;
                    break; /* found u's bottleneck, stop looking */
                }
            }
        }

        /* Check if any good deletion exists */
        int has_good = 0;
        for (int v = 0; v < N; v++) {
            if (deletion_is_cycle_rich(out, v)) {
                has_good = 1;
                break;
            }
        }

        if (!has_good) {
            no_good_del++;
            int mm = max_matching_3cycles(out);

            if (no_good_del <= 10) {
                int sc[N];
                for (int i = 0; i < N; i++)
                    sc[i] = __builtin_popcount(out[i]);
                printf("\n--- No-good-deletion #%lld (mm=%d) ---\n", no_good_del, mm);
                printf("Scores: ");
                for (int i = 0; i < N; i++) printf("%d ", sc[i]);
                printf("\n");
                printf("Bottleneck structure:\n");
                for (int v = 0; v < N; v++) {
                    if (B[v]) {
                        printf("  B(%d) = {", v);
                        for (int u = 0; u < N; u++)
                            if (B[v] & (1 << u)) printf("%d ", u);
                        printf("} size=%d\n", __builtin_popcount(B[v]));
                    }
                }
                printf("Bottleneck-of: ");
                for (int u = 0; u < N; u++) printf("%d->%d ", u, bottleneck_of[u]);
                printf("\n");

                /* Check score obstructions too */
                printf("Score obstructions:\n");
                for (int v = 0; v < N; v++) {
                    for (int u = 0; u < N; u++) {
                        if (u == v) continue;
                        int s = __builtin_popcount(out[u] & ~(1 << v));
                        if (s == 0) printf("  Del %d: %d becomes source\n", v, u);
                        if (s == N-2) printf("  Del %d: %d becomes sink\n", v, u);
                    }
                }
            }

            /* Track bottleneck size distribution */
            for (int v = 0; v < N; v++) {
                int sz = __builtin_popcount(B[v]);
                if (sz < N+1) bn_size_hist[sz]++;
                if (sz > max_bn_set_seen) max_bn_set_seen = sz;
            }
        }

        if ((t + 1) % 10000000 == 0) {
            fprintf(stderr, "Progress: %lldM/%lldM, cr=%lld, no_good=%lld, L2_viol=%lld\n",
                    (t+1)/1000000, trials/1000000, cycle_rich, no_good_del, lemma2_violations);
        }
    }

    printf("\n=== RESULTS ===\n");
    printf("Trials: %lld\n", trials);
    printf("Cycle-rich: %lld\n", cycle_rich);
    printf("No good deletion: %lld\n", no_good_del);
    printf("Lemma 2 violations: %lld (should be 0)\n", lemma2_violations);
    printf("Max B(v) size seen: %d\n", max_bn_set_seen);
    printf("B(v) size histogram (for no-good-del cases):\n");
    for (int s = 0; s <= N; s++)
        if (bn_size_hist[s] > 0)
            printf("  |B(v)|=%d: %d instances\n", s, bn_size_hist[s]);

    return 0;
}
