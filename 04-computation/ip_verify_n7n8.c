/*
 * EXHAUSTIVE VERIFICATION: H = IP(OddCycleDisjointnessGraph, 2)
 * opus-2026-03-14-S89b
 *
 * n=7: H = 1 + 2(t3+t5+t7) + 4*d33
 * n=8: H = 1 + 2(t3+t5+t7) + 4*(d33+d35) + ???
 *
 * Uses DP bitmask for H, direct counting for cycles.
 * Parallelized with OpenMP.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXN 8

/* Tournament adjacency: adj[i][j] = 1 if i->j */
static int adj[MAXN][MAXN];
static int N;

/* Build tournament from bit vector */
static void build_tournament(int n, uint32_t bits) {
    int idx = 0;
    memset(adj, 0, sizeof(adj));
    for (int i = 0; i < n; i++)
        for (int j = i+1; j < n; j++) {
            if ((bits >> idx) & 1) {
                adj[i][j] = 1;
            } else {
                adj[j][i] = 1;
            }
            idx++;
        }
}

/* DP bitmask Hamiltonian path count */
static int64_t compute_H(int n) {
    /* dp[mask][v] = number of Hamiltonian paths ending at v using vertices in mask */
    int states = 1 << n;
    int64_t *dp = calloc((size_t)states * n, sizeof(int64_t));
    if (!dp) { fprintf(stderr, "OOM\n"); exit(1); }

    for (int v = 0; v < n; v++)
        dp[(1 << v) * n + v] = 1;

    for (int mask = 1; mask < states; mask++) {
        for (int v = 0; v < n; v++) {
            if (!(mask & (1 << v))) continue;
            int64_t val = dp[mask * n + v];
            if (val == 0) continue;
            for (int u = 0; u < n; u++) {
                if (mask & (1 << u)) continue;
                if (adj[v][u]) {
                    dp[(mask | (1 << u)) * n + u] += val;
                }
            }
        }
    }

    int64_t total = 0;
    int full_mask = states - 1;
    for (int v = 0; v < n; v++)
        total += dp[full_mask * n + v];

    free(dp);
    return total;
}

/* Count 3-cycles; also store their vertex sets */
#define MAX_TRIPLES 100
static int triple_verts[MAX_TRIPLES][3];
static int n_triples;

static int count_3cycles(int n) {
    n_triples = 0;
    int count = 0;
    for (int i = 0; i < n; i++)
        for (int j = i+1; j < n; j++)
            for (int k = j+1; k < n; k++) {
                if ((adj[i][j] && adj[j][k] && adj[k][i]) ||
                    (adj[i][k] && adj[k][j] && adj[j][i])) {
                    if (n_triples < MAX_TRIPLES) {
                        triple_verts[n_triples][0] = i;
                        triple_verts[n_triples][1] = j;
                        triple_verts[n_triples][2] = k;
                    }
                    n_triples++;
                    count++;
                }
            }
    return count;
}

/* Count directed 5-cycles on a given 5-subset */
static int count_5cycles_subset(int v[5]) {
    /* Try all 24 permutations of 4 elements (fix v[0]) */
    int count = 0;
    int perm[5];
    perm[0] = v[0];
    int rest[4] = {v[1], v[2], v[3], v[4]};

    /* Generate all 24 permutations of rest */
    for (int a = 0; a < 4; a++)
    for (int b = 0; b < 4; b++) { if (b==a) continue;
    for (int c = 0; c < 4; c++) { if (c==a||c==b) continue;
    for (int d = 0; d < 4; d++) { if (d==a||d==b||d==c) continue;
        perm[1] = rest[a]; perm[2] = rest[b];
        perm[3] = rest[c]; perm[4] = rest[d];
        if (adj[perm[0]][perm[1]] &&
            adj[perm[1]][perm[2]] &&
            adj[perm[2]][perm[3]] &&
            adj[perm[3]][perm[4]] &&
            adj[perm[4]][perm[0]])
            count++;
    }}}
    return count; /* Each directed 5-cycle counted once (start fixed at v[0]) */
}

/* Count all directed 5-cycles in tournament */
static int count_5cycles(int n) {
    int total = 0;
    int v[5];
    for (v[0]=0; v[0]<n; v[0]++)
    for (v[1]=v[0]+1; v[1]<n; v[1]++)
    for (v[2]=v[1]+1; v[2]<n; v[2]++)
    for (v[3]=v[2]+1; v[3]<n; v[3]++)
    for (v[4]=v[3]+1; v[4]<n; v[4]++) {
        /* Count directed 5-cycles on this subset */
        int sub_count = count_5cycles_subset(v);
        total += sub_count;
    }
    /* Each directed 5-cycle on a 5-set is counted once (start fixed).
     * But across subsets, each is counted once since we enumerate each subset once.
     * Wait — fixing the start at v[0] means we count each cycle starting at the
     * smallest vertex. But a 5-cycle v0->v1->v2->v3->v4->v0 might not have
     * v0 as the smallest. We're iterating over all 5-subsets and for each one
     * fixing v[0] as the smallest element. So each directed cycle is counted
     * exactly once (the unique start at the smallest vertex in the cycle).
     * This IS the count of directed 5-cycles. */
    return total;
}

/* Count directed 7-cycles (Hamiltonian cycles for n=7, or subsets for n=8) */
static int count_7cycles_subset(int v[7]) {
    /* Fix v[0], permute rest (6! = 720) */
    int count = 0;
    int rest[6];
    for (int i = 0; i < 6; i++) rest[i] = v[i+1];

    /* Generate all 720 permutations of 6 elements */
    for (int a=0; a<6; a++) {
    for (int b=0; b<6; b++) { if(b==a) continue;
    for (int c=0; c<6; c++) { if(c==a||c==b) continue;
    for (int d=0; d<6; d++) { if(d==a||d==b||d==c) continue;
    for (int e=0; e<6; e++) { if(e==a||e==b||e==c||e==d) continue;
    for (int f=0; f<6; f++) { if(f==a||f==b||f==c||f==d||f==e) continue;
        if (adj[v[0]][rest[a]] &&
            adj[rest[a]][rest[b]] &&
            adj[rest[b]][rest[c]] &&
            adj[rest[c]][rest[d]] &&
            adj[rest[d]][rest[e]] &&
            adj[rest[e]][rest[f]] &&
            adj[rest[f]][v[0]])
            count++;
    }}}}}
    }
    return count;
}

static int count_7cycles(int n) {
    if (n < 7) return 0;
    int total = 0;
    int v[7];
    for (v[0]=0; v[0]<n; v[0]++)
    for (v[1]=v[0]+1; v[1]<n; v[1]++)
    for (v[2]=v[1]+1; v[2]<n; v[2]++)
    for (v[3]=v[2]+1; v[3]<n; v[3]++)
    for (v[4]=v[3]+1; v[4]<n; v[4]++)
    for (v[5]=v[4]+1; v[5]<n; v[5]++)
    for (v[6]=v[5]+1; v[6]<n; v[6]++) {
        total += count_7cycles_subset(v);
    }
    return total;
}

/* Count disjoint 3-cycle pairs */
static int count_d33(void) {
    int count = 0;
    for (int i = 0; i < n_triples; i++)
        for (int j = i+1; j < n_triples; j++) {
            int disjoint = 1;
            for (int a = 0; a < 3 && disjoint; a++)
                for (int b = 0; b < 3 && disjoint; b++)
                    if (triple_verts[i][a] == triple_verts[j][b])
                        disjoint = 0;
            count += disjoint;
        }
    return count;
}

/* Count disjoint (3-cycle, 5-cycle) pairs — for n=8 */
/* A 3-cycle on {a,b,c} and a 5-cycle on 5 vertices disjoint from {a,b,c} */
static int count_d35(int n) {
    if (n < 8) return 0;
    int count = 0;
    /* For each 3-cycle triple */
    for (int t = 0; t < n_triples; t++) {
        int used[MAXN] = {0};
        used[triple_verts[t][0]] = 1;
        used[triple_verts[t][1]] = 1;
        used[triple_verts[t][2]] = 1;

        /* Collect remaining vertices */
        int rem[MAXN], nr = 0;
        for (int v = 0; v < n; v++)
            if (!used[v]) rem[nr++] = v;

        /* Count 5-cycles on subsets of remaining vertices */
        if (nr < 5) continue;
        int v5[5];
        for (int a=0; a<nr; a++)
        for (int b=a+1; b<nr; b++)
        for (int c=b+1; c<nr; c++)
        for (int d=c+1; d<nr; d++)
        for (int e=d+1; e<nr; e++) {
            v5[0]=rem[a]; v5[1]=rem[b]; v5[2]=rem[c]; v5[3]=rem[d]; v5[4]=rem[e];
            count += count_5cycles_subset(v5);
        }
    }
    return count;
}

int main(int argc, char *argv[]) {
    int n = 7;
    if (argc > 1) n = atoi(argv[1]);
    if (n < 3 || n > 8) {
        fprintf(stderr, "Usage: %s [n=7|8]\n", argv[0]);
        return 1;
    }

    int m = n*(n-1)/2;
    uint64_t total = 1ULL << m;

    printf("EXHAUSTIVE VERIFICATION: H = IP(OddCycGraph, 2)\n");
    printf("n=%d, m=%d, total=%llu tournaments\n", n, m, (unsigned long long)total);
    printf("Formula: H = 1 + 2(t3+t5+t7) + 4*d33");
    if (n >= 8) printf(" + 4*d35");
    printf("\n\n");

    int64_t errors = 0;
    int64_t max_error = 0;
    uint64_t progress_step = total / 20;
    if (progress_step == 0) progress_step = 1;

    for (uint64_t bits = 0; bits < total; bits++) {
        if (bits % progress_step == 0) {
            fprintf(stderr, "  %llu%% done (errors: %lld)\n",
                    (unsigned long long)(100*bits/total), (long long)errors);
        }

        build_tournament(n, (uint32_t)bits);
        int64_t H = compute_H(n);

        int t3 = count_3cycles(n);
        int t5 = count_5cycles(n);
        int t7 = count_7cycles(n);
        int d33 = count_d33();

        int64_t predicted = 1 + 2*(t3 + t5 + t7) + 4*d33;

        if (n >= 8) {
            int d35 = count_d35(n);
            predicted += 4 * d35;
        }

        if (predicted != H) {
            errors++;
            int64_t err = H - predicted;
            if (err < 0) err = -err;
            if (err > max_error) max_error = err;
            if (errors <= 10) {
                printf("  MISMATCH #%lld: bits=%llu, H=%lld, pred=%lld, "
                       "t3=%d, t5=%d, t7=%d, d33=%d, diff=%lld\n",
                       (long long)errors, (unsigned long long)bits,
                       (long long)H, (long long)predicted,
                       t3, t5, t7, d33, (long long)(H - predicted));
            }
        }
    }

    printf("\n==============================\n");
    if (errors == 0) {
        printf("*** PERFECT: Formula EXACT for ALL %llu n=%d tournaments! ***\n",
               (unsigned long long)total, n);
    } else {
        printf("FAILED: %lld/%llu mismatches, max error = %lld\n",
               (long long)errors, (unsigned long long)total, (long long)max_error);
    }
    printf("==============================\n");

    return 0;
}
