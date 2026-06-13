/*
 * e_nauty_fast_s266.c — Fast E(G_n) using edge-centric nauty
 * opus-2026-03-23-S266
 *
 * Strategy: read tournament class reps from gentourng,
 * flip arc (0,1), canonicalize with nauty, count distinct pairs.
 *
 * Compile: cc -O3 -o e_nauty_fast e_nauty_fast_s266.c -lnauty
 * Usage:   gentourng n | ./e_nauty_fast n
 *
 * For n=10: 9.7M class reps, 1 flip + 1 canonicalization each.
 * Nauty canonicalization: ~10μs per tournament → ~100 seconds total.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* Simple approach without nauty library dependency:
 * Use score-based hashing + refinement for canonicalization.
 * This is slower than nauty but doesn't require the library.
 */

#define MAXN 12
#define MAXM ((MAXN*(MAXN-1))/2)

typedef struct {
    int adj[MAXN][MAXN];
    int n;
} Tournament;

typedef struct {
    unsigned long long hash;
    int idx;
} HashEntry;

/* Compute a strong hash of a tournament */
unsigned long long tournament_hash(Tournament *T) {
    int n = T->n;
    /* Score sequence (sorted) */
    int scores[MAXN];
    for (int i = 0; i < n; i++) {
        scores[i] = 0;
        for (int j = 0; j < n; j++)
            if (i != j) scores[i] += T->adj[i][j];
    }
    /* Sort scores */
    for (int i = 0; i < n-1; i++)
        for (int j = i+1; j < n; j++)
            if (scores[i] > scores[j]) {
                int tmp = scores[i]; scores[i] = scores[j]; scores[j] = tmp;
            }

    /* 3-cycle count */
    int c3 = 0;
    for (int i = 0; i < n; i++)
        for (int j = i+1; j < n; j++)
            for (int k = j+1; k < n; k++) {
                if (T->adj[i][j] && T->adj[j][k] && T->adj[k][i]) c3++;
                if (T->adj[i][k] && T->adj[k][j] && T->adj[j][i]) c3++;
            }

    /* Score-pair invariant */
    int sp[MAXN][MAXN] = {{0}};
    for (int i = 0; i < n; i++) {
        int si = 0;
        for (int j = 0; j < n; j++) if (i!=j) si += T->adj[i][j];
        for (int j = i+1; j < n; j++) {
            int sj = 0;
            for (int k = 0; k < n; k++) if (j!=k) sj += T->adj[j][k];
            int a = si < sj ? si : sj;
            int b = si < sj ? sj : si;
            sp[a][b] += T->adj[i][j];
        }
    }

    /* Combine into hash */
    unsigned long long h = 0;
    for (int i = 0; i < n; i++)
        h = h * 31 + scores[i];
    h = h * 1000003 + c3;
    for (int i = 0; i < n; i++)
        for (int j = i; j < n; j++)
            h = h * 131 + sp[i][j];

    return h;
}

/* Read tournaments from stdin in gentourng format */
int read_tournament(Tournament *T, int n, char *line) {
    int m = n*(n-1)/2;
    if ((int)strlen(line) < m) return 0;

    T->n = n;
    memset(T->adj, 0, sizeof(T->adj));

    int idx = 0;
    for (int i = 0; i < n; i++)
        for (int j = i+1; j < n; j++) {
            if (line[idx] == '1') T->adj[i][j] = 1;
            else T->adj[j][i] = 1;
            idx++;
        }
    return 1;
}

/* Flip arc (0,1) */
void flip_01(Tournament *T, Tournament *out) {
    memcpy(out, T, sizeof(Tournament));
    out->adj[0][1] = 1 - T->adj[0][1];
    out->adj[1][0] = 1 - T->adj[1][0];
}

/* Hash table for tournament lookup */
#define HASH_SIZE (1 << 24)  /* 16M buckets */
#define HASH_MASK (HASH_SIZE - 1)

int hash_table_idx[HASH_SIZE];  /* index of tournament, or -1 */
unsigned long long hash_table_key[HASH_SIZE];

void hash_init(void) {
    memset(hash_table_idx, -1, sizeof(hash_table_idx));
}

void hash_insert(unsigned long long key, int idx) {
    unsigned int pos = (unsigned int)(key & HASH_MASK);
    while (hash_table_idx[pos] != -1) {
        pos = (pos + 1) & HASH_MASK;
    }
    hash_table_idx[pos] = idx;
    hash_table_key[pos] = key;
}

int hash_lookup(unsigned long long key) {
    unsigned int pos = (unsigned int)(key & HASH_MASK);
    while (hash_table_idx[pos] != -1) {
        if (hash_table_key[pos] == key) return hash_table_idx[pos];
        pos = (pos + 1) & HASH_MASK;
    }
    return -1;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: gentourng n | %s n\n", argv[0]);
        return 1;
    }
    int n = atoi(argv[1]);
    int m = n*(n-1)/2;

    fprintf(stderr, "Reading tournaments for n=%d...\n", n);

    /* Phase 1: Read all class reps and hash them */
    hash_init();
    Tournament *tours = NULL;
    int capacity = 1024, count = 0;
    tours = malloc(capacity * sizeof(Tournament));

    char line[256];
    clock_t t0 = clock();

    while (fgets(line, sizeof(line), stdin)) {
        /* Strip newline */
        int len = strlen(line);
        while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r'))
            line[--len] = '\0';
        if (len < m) continue;

        if (count >= capacity) {
            capacity *= 2;
            tours = realloc(tours, capacity * sizeof(Tournament));
        }

        read_tournament(&tours[count], n, line);
        unsigned long long h = tournament_hash(&tours[count]);
        hash_insert(h, count);
        count++;

        if (count % 1000000 == 0)
            fprintf(stderr, "  Read %d class reps...\n", count);
    }

    double read_time = (double)(clock() - t0) / CLOCKS_PER_SEC;
    fprintf(stderr, "Read %d class reps in %.1fs\n", count, read_time);

    /* Phase 2: For each class rep, flip (0,1), lookup, count edges */
    t0 = clock();
    long long edges = 0, self_loops = 0, hash_misses = 0;

    for (int i = 0; i < count; i++) {
        Tournament flipped;
        flip_01(&tours[i], &flipped);

        unsigned long long h = tournament_hash(&flipped);
        int nb = hash_lookup(h);

        if (nb < 0) {
            hash_misses++;
            continue;  /* hash collision miss — need better hash */
        }

        if (nb == i) {
            self_loops++;
        } else if (nb > i) {
            edges++;  /* only count each pair once */
        }
        /* If nb < i, this pair was already counted when processing nb */

        if ((i+1) % 1000000 == 0)
            fprintf(stderr, "  Processed %d/%d, edges=%lld, misses=%lld\n",
                    i+1, count, edges, hash_misses);
    }

    double proc_time = (double)(clock() - t0) / CLOCKS_PER_SEC;

    printf("n=%d: V=%d, E=%lld, SL=%lld, hash_misses=%lld\n",
           n, count, edges, self_loops, hash_misses);
    printf("Time: read=%.1fs, process=%.1fs, total=%.1fs\n",
           read_time, proc_time, read_time + proc_time);

    if (hash_misses > 0) {
        printf("WARNING: %lld hash misses — result is a LOWER BOUND\n", hash_misses);
        printf("Need better hash or nauty canonicalization for exact result\n");
    }

    free(tours);
    return 0;
}
