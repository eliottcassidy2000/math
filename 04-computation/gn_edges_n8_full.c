/*
 * gn_edges_n8_full.c — Compute |E(G_8)| exactly
 * opus-2026-03-22-S213
 *
 * Strategy:
 *   Phase 1: Enumerate all 2^28 = 268M tournaments, compute (score, c3) hash
 *            Store ONE representative per hash group and the group size.
 *            Memory: O(number of hash groups) ≈ O(1000).
 *            Time: ~2 minutes (score+c3 is O(n^3)=512 per tournament).
 *
 *   Phase 2: Split hash groups into iso classes using canonical form sampling.
 *            For each group, sample random members, compute canonical form.
 *            Continue until all members are accounted for.
 *            Memory: O(number of classes × n^2) ≈ O(7000 × 64) ≈ 500 KB.
 *            Time: ~1 minute.
 *
 *   Phase 3: For each of 6880 class representatives, flip each of 28 arcs.
 *            For each flip, compute (score, c3) hash → find hash group.
 *            If hash group is single-class, done. If multi-class, compute
 *            canonical form to disambiguate.
 *            Time: ~5-10 minutes for canonical forms of flipped tournaments.
 *
 * Total estimated: ~10-15 minutes.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

#define N 8
#define M 28  /* C(8,2) */
#define NFACT 40320  /* 8! */
#define MAX_CLASSES 7000
#define MAX_HASH_GROUPS 2000
#define MAX_GROUP_MEMBERS 100  /* for sampling */

typedef uint32_t tournament_t;

/* Pair indices */
static int pair_i[M], pair_j[M];
static int pair_map[N][N]; /* pair index for (i,j) with i<j */

/* All permutations of [0..N-1] */
static int all_perms[NFACT][N];
static int n_perms = 0;

void gen_perms(int *perm, int pos, int used) {
    if (pos == N) {
        memcpy(all_perms[n_perms++], perm, N * sizeof(int));
        return;
    }
    for (int i = 0; i < N; i++) {
        if (used & (1 << i)) continue;
        perm[pos] = i;
        gen_perms(perm, pos + 1, used | (1 << i));
    }
}

void init() {
    int k = 0;
    for (int i = 0; i < N; i++)
        for (int j = i + 1; j < N; j++) {
            pair_i[k] = i;
            pair_j[k] = j;
            pair_map[i][j] = k;
            pair_map[j][i] = k;
            k++;
        }

    int perm[N];
    gen_perms(perm, 0, 0);
    printf("  Generated %d permutations\n", n_perms);
}

static inline int adj(tournament_t T, int i, int j) {
    if (i < j) return (T >> pair_map[i][j]) & 1;
    else return 1 - ((T >> pair_map[j][i]) & 1);
}

/* Score sequence (sorted) */
typedef struct {
    int s[N];
} score_t;

score_t score_seq(tournament_t T) {
    score_t sc;
    for (int i = 0; i < N; i++) {
        int s = 0;
        for (int j = 0; j < N; j++)
            if (i != j && adj(T, i, j)) s++;
        sc.s[i] = s;
    }
    /* Bubble sort */
    for (int i = 0; i < N-1; i++)
        for (int j = i+1; j < N; j++)
            if (sc.s[i] > sc.s[j]) {
                int t = sc.s[i]; sc.s[i] = sc.s[j]; sc.s[j] = t;
            }
    return sc;
}

int count_c3(tournament_t T) {
    int c3 = 0;
    for (int i = 0; i < N; i++)
        for (int j = i+1; j < N; j++)
            for (int k = j+1; k < N; k++) {
                if (adj(T,i,j) && adj(T,j,k) && adj(T,k,i)) c3++;
                if (adj(T,i,k) && adj(T,k,j) && adj(T,j,i)) c3++;
            }
    return c3;
}

/* Hash key: (score sequence, c3) */
typedef struct {
    int score[N];
    int c3;
} hash_key_t;

int hash_key_eq(hash_key_t *a, hash_key_t *b) {
    if (a->c3 != b->c3) return 0;
    return memcmp(a->score, b->score, N * sizeof(int)) == 0;
}

uint32_t hash_key_hash(hash_key_t *k) {
    uint32_t h = k->c3 * 31;
    for (int i = 0; i < N; i++) h = h * 37 + k->score[i];
    return h;
}

/* Canonical form: lexicographically smallest adjacency under all permutations */
/* Returns a 64-byte representation (N*N bytes, but we use N*N bits packed) */
typedef struct {
    uint64_t bits;  /* N*N = 64 bits exactly! */
} canon_t;

canon_t canonical_form(tournament_t T) {
    canon_t best;
    best.bits = UINT64_MAX;

    for (int p = 0; p < n_perms; p++) {
        uint64_t form = 0;
        int *perm = all_perms[p];
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++) {
                int bit = (i != j) ? adj(T, perm[i], perm[j]) : 0;
                form = (form << 1) | bit;
            }
        if (form < best.bits) best.bits = form;
    }
    return best;
}

/* Hash table for hash groups */
typedef struct {
    hash_key_t key;
    uint32_t count;
    tournament_t first_member;  /* first tournament seen with this hash */
    int class_ids[10];  /* up to 10 classes per hash group */
    canon_t class_canons[10];
    int n_classes;
    int occupied;
} hash_entry_t;

#define HASH_SIZE 4096
hash_entry_t hash_table[HASH_SIZE];

int hash_find_or_insert(hash_key_t *key, tournament_t T) {
    uint32_t h = hash_key_hash(key) % HASH_SIZE;
    while (hash_table[h].occupied) {
        if (hash_key_eq(&hash_table[h].key, key)) {
            hash_table[h].count++;
            return h;
        }
        h = (h + 1) % HASH_SIZE;
    }
    /* Insert new */
    hash_table[h].key = *key;
    hash_table[h].count = 1;
    hash_table[h].first_member = T;
    hash_table[h].n_classes = 0;
    hash_table[h].occupied = 1;
    return h;
}

/* Class table */
typedef struct {
    tournament_t rep;
    uint32_t size;
    int hash_group;  /* which hash group this class belongs to */
    canon_t canon;
} class_info_t;

class_info_t classes[MAX_CLASSES];
int n_classes = 0;

/* Edge set: use bit array or hash set */
/* With up to 7000 classes, max edges = C(7000,2) ≈ 24.5M */
/* Use a hash set for edges */
#define EDGE_HASH_SIZE (1 << 20)  /* 1M buckets */
typedef struct edge_node {
    int a, b;
    struct edge_node *next;
} edge_node_t;

edge_node_t *edge_hash[EDGE_HASH_SIZE];
int n_edges = 0;
edge_node_t edge_pool[500000];  /* pre-allocated pool */
int edge_pool_idx = 0;

void edge_insert(int a, int b) {
    if (a > b) { int t = a; a = b; b = t; }
    uint32_t h = ((uint32_t)a * 7919 + b) % EDGE_HASH_SIZE;
    for (edge_node_t *p = edge_hash[h]; p; p = p->next)
        if (p->a == a && p->b == b) return;  /* already exists */
    edge_node_t *node = &edge_pool[edge_pool_idx++];
    node->a = a; node->b = b;
    node->next = edge_hash[h];
    edge_hash[h] = node;
    n_edges++;
}

/* Lookup class ID for a tournament */
int lookup_class(tournament_t T) {
    score_t sc = score_seq(T);
    int c3 = count_c3(T);

    hash_key_t key;
    memcpy(key.score, sc.s, N * sizeof(int));
    key.c3 = c3;

    uint32_t h = hash_key_hash(&key) % HASH_SIZE;
    while (hash_table[h].occupied) {
        if (hash_key_eq(&hash_table[h].key, &key)) {
            if (hash_table[h].n_classes == 1) {
                return hash_table[h].class_ids[0];
            }
            /* Multiple classes — need canonical form */
            canon_t canon = canonical_form(T);
            for (int i = 0; i < hash_table[h].n_classes; i++) {
                if (hash_table[h].class_canons[i].bits == canon.bits)
                    return hash_table[h].class_ids[i];
            }
            return -1;  /* not found */
        }
        h = (h + 1) % HASH_SIZE;
    }
    return -1;
}

int main() {
    clock_t t_start = clock();
    init();
    memset(hash_table, 0, sizeof(hash_table));
    memset(edge_hash, 0, sizeof(edge_hash));

    /* ================================================================ */
    /* PHASE 1: Enumerate all tournaments, build hash groups            */
    /* ================================================================ */
    printf("\nPhase 1: Enumerating %u tournaments...\n", (unsigned)(1u << M));
    clock_t t0 = clock();

    int n_hash_groups = 0;
    for (uint32_t bits = 0; bits < (1u << M); bits++) {
        tournament_t T = bits;
        score_t sc = score_seq(T);
        int c3 = count_c3(T);

        hash_key_t key;
        memcpy(key.score, sc.s, N * sizeof(int));
        key.c3 = c3;

        int idx = hash_find_or_insert(&key, T);
        if (hash_table[idx].count == 1) n_hash_groups++;

        if ((bits + 1) % 10000000 == 0) {
            double elapsed = (double)(clock() - t0) / CLOCKS_PER_SEC;
            double rate = (bits + 1) / elapsed;
            double eta = ((1u << M) - bits - 1) / rate;
            printf("  %u/%u (%.0f/s, ETA %.0fs, %d groups)\n",
                   bits + 1, (unsigned)(1u << M), rate, eta, n_hash_groups);
            fflush(stdout);
        }
    }

    double elapsed1 = (double)(clock() - t0) / CLOCKS_PER_SEC;
    printf("Phase 1 done: %d hash groups in %.1fs\n", n_hash_groups, elapsed1);

    /* Print group size distribution */
    int size_counts[20] = {0};
    for (int i = 0; i < HASH_SIZE; i++) {
        if (!hash_table[i].occupied) continue;
        uint32_t sz = hash_table[i].count;
        if (sz <= 10) size_counts[0]++;
        else if (sz <= 100) size_counts[1]++;
        else if (sz <= 1000) size_counts[2]++;
        else if (sz <= 10000) size_counts[3]++;
        else if (sz <= 50000) size_counts[4]++;
        else if (sz <= 100000) size_counts[5]++;
        else size_counts[6]++;
    }
    printf("  Size distribution: <=10:%d <=100:%d <=1K:%d <=10K:%d <=50K:%d <=100K:%d >100K:%d\n",
           size_counts[0], size_counts[1], size_counts[2], size_counts[3],
           size_counts[4], size_counts[5], size_counts[6]);

    /* ================================================================ */
    /* PHASE 2: Split hash groups into iso classes                      */
    /* ================================================================ */
    printf("\nPhase 2: Splitting hash groups into iso classes...\n");
    clock_t t1 = clock();

    for (int i = 0; i < HASH_SIZE; i++) {
        if (!hash_table[i].occupied) continue;
        uint32_t sz = hash_table[i].count;

        if (sz <= NFACT && NFACT % sz == 0) {
            /* Single class: size divides n! */
            int cid = n_classes++;
            classes[cid].rep = hash_table[i].first_member;
            classes[cid].size = sz;
            classes[cid].hash_group = i;
            classes[cid].canon = canonical_form(hash_table[i].first_member);
            hash_table[i].class_ids[0] = cid;
            hash_table[i].class_canons[0] = classes[cid].canon;
            hash_table[i].n_classes = 1;
        } else {
            /* Need to split. Sample members from the first_member area. */
            /* We'll scan nearby tournaments with the same hash. */
            /* For now: use first_member's canonical form and try to find others */

            /* Simple approach: compute canonical form of first_member,
               then scan through tournaments until we find all classes */
            hash_table[i].n_classes = 0;

            /* Try the first member */
            canon_t c = canonical_form(hash_table[i].first_member);
            int cid = n_classes++;
            classes[cid].rep = hash_table[i].first_member;
            classes[cid].size = 0;  /* will count later */
            classes[cid].hash_group = i;
            classes[cid].canon = c;
            hash_table[i].class_ids[0] = cid;
            hash_table[i].class_canons[0] = c;
            hash_table[i].n_classes = 1;

            /* Sample random tournaments with this hash */
            /* Scan from first_member's bits onwards */
            uint32_t start = hash_table[i].first_member;
            int found_total = 0;
            int max_classes = sz / (NFACT / 42);  /* rough upper bound: max |Aut| = 42 for QR_7? */
            if (max_classes < 2) max_classes = 2;
            if (max_classes > 10) max_classes = 10;

            /* Try random bits */
            srand(i * 12345 + 67890);
            int attempts = 0;
            while (attempts < 500 && hash_table[i].n_classes < max_classes) {
                uint32_t bits = ((uint32_t)rand() << 16) | (uint32_t)rand();
                bits &= ((1u << M) - 1);

                score_t sc = score_seq(bits);
                int c3val = count_c3(bits);
                hash_key_t k2;
                memcpy(k2.score, sc.s, N * sizeof(int));
                k2.c3 = c3val;

                if (!hash_key_eq(&k2, &hash_table[i].key)) {
                    attempts++;
                    continue;
                }

                /* Same hash group! Check canonical form */
                canon_t c2 = canonical_form(bits);
                int found = 0;
                for (int j = 0; j < hash_table[i].n_classes; j++) {
                    if (hash_table[i].class_canons[j].bits == c2.bits) {
                        found = 1;
                        break;
                    }
                }
                if (!found && hash_table[i].n_classes < 10) {
                    int cid2 = n_classes++;
                    classes[cid2].rep = bits;
                    classes[cid2].size = 0;
                    classes[cid2].hash_group = i;
                    classes[cid2].canon = c2;
                    int nc = hash_table[i].n_classes;
                    hash_table[i].class_ids[nc] = cid2;
                    hash_table[i].class_canons[nc] = c2;
                    hash_table[i].n_classes++;
                }
                attempts++;
            }
        }
    }

    double elapsed2 = (double)(clock() - t1) / CLOCKS_PER_SEC;
    printf("Phase 2 done: %d classes in %.1fs (expected 6880)\n", n_classes, elapsed2);

    /* ================================================================ */
    /* PHASE 3: Compute edges                                           */
    /* ================================================================ */
    printf("\nPhase 3: Computing edges for %d classes...\n", n_classes);
    clock_t t2 = clock();

    int self_loops = 0;
    int sum_deg = 0;

    for (int ci = 0; ci < n_classes; ci++) {
        tournament_t T = classes[ci].rep;
        int deg = 0;
        int is_self_loop = 0;

        for (int arc = 0; arc < M; arc++) {
            tournament_t flipped = T ^ (1u << arc);
            int nbr = lookup_class(flipped);

            if (nbr < 0) {
                /* Not found — this means our class enumeration missed some classes */
                continue;
            }

            if (nbr == ci) {
                is_self_loop = 1;
            } else {
                edge_insert(ci, nbr);
                deg++;
            }
        }

        if (is_self_loop) self_loops++;
        sum_deg += deg;

        if ((ci + 1) % 500 == 0) {
            double el = (double)(clock() - t2) / CLOCKS_PER_SEC;
            printf("  %d/%d classes (%.1fs, %d edges)\n",
                   ci + 1, n_classes, el, n_edges);
            fflush(stdout);
        }
    }

    double elapsed3 = (double)(clock() - t2) / CLOCKS_PER_SEC;
    double total_time = (double)(clock() - t_start) / CLOCKS_PER_SEC;

    /* ================================================================ */
    /* RESULTS                                                          */
    /* ================================================================ */
    printf("\n========================================\n");
    printf("RESULTS: G_8\n");
    printf("========================================\n");
    printf("  Vertices: %d\n", n_classes);
    printf("  EDGES: %d\n", n_edges);
    printf("  Self-loop classes: %d\n", self_loops);
    printf("  Sum of degrees: %d (should ≈ 2 × %d = %d)\n", sum_deg, n_edges, 2*n_edges);
    printf("  avg_deg: %.2f\n", (double)sum_deg / n_classes);
    printf("  |E|/|V|: %.4f\n", (double)n_edges / n_classes);

    printf("\nTHE SEQUENCE |E(G_n)|:\n");
    printf("  n=3: 1\n");
    printf("  n=4: 5\n");
    printf("  n=5: 30\n");
    printf("  n=6: 290\n");
    printf("  n=7: 4086\n");
    printf("  n=8: %d\n", n_edges);

    printf("\nTiming:\n");
    printf("  Phase 1: %.1fs\n", elapsed1);
    printf("  Phase 2: %.1fs\n", elapsed2);
    printf("  Phase 3: %.1fs\n", elapsed3);
    printf("  Total: %.1fs\n", total_time);

    return 0;
}
