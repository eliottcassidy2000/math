/*
 * wn_compact.c — Compute W(n) with minimal memory using compact storage
 * opus-2026-03-15-S89c
 *
 * Key insight: at each popcount level, most dp entries are zero.
 * Only store nonzero (mask, v, value) triples using a flat array.
 * Sort by mask for cache efficiency.
 *
 * For the next level, use a hash map with open addressing keyed by
 * (mask, v) packed into a single 64-bit key: mask in upper bits, v in lower.
 *
 * Compile: gcc -O3 -o wn_compact wn_compact.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef unsigned __int128 u128;
typedef unsigned long long u64;

void print_u128(u128 x) {
    if (x == 0) { printf("0"); return; }
    char buf[50];
    int pos = 0;
    while (x > 0) {
        buf[pos++] = '0' + (int)(x % 10);
        x /= 10;
    }
    for (int i = pos - 1; i >= 0; i--)
        putchar(buf[i]);
}

/* Entry in the sparse DP table */
typedef struct {
    unsigned int mask;
    unsigned char v;
    u128 val;
} entry_t;

/* Hash table for accumulating next level */
typedef struct {
    u64 key;     /* (mask << 5) | v, or 0xFFFFFFFFFFFFFFFF for empty */
    u128 val;
} ht_entry_t;

#define HT_EMPTY 0xFFFFFFFFFFFFFFFFULL

static inline u64 make_key(unsigned int mask, int v) {
    return ((u64)mask << 5) | (unsigned)v;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <n>\n", argv[0]);
        return 1;
    }
    int n = atoi(argv[1]);
    if (n < 2 || n > 30) {
        fprintf(stderr, "n must be 2..30\n");
        return 1;
    }

    time_t t0 = time(NULL);

    /* Start: level 1 entries */
    long long curr_count = n;
    entry_t *curr = (entry_t *)malloc(n * sizeof(entry_t));
    for (int v = 0; v < n; v++) {
        curr[v].mask = 1u << v;
        curr[v].v = v;
        curr[v].val = 1;
    }

    for (int p = 1; p < n; p++) {
        /* Estimate next level size: each entry generates up to (n-p) new entries,
         * but many collide. Use generous hash table.
         * Estimated unique entries: curr_count * (n-p) / 2 (rough) */
        long long est_next = curr_count * (n - p);
        if (est_next > 500000000LL) est_next = 500000000LL; /* cap at 500M */

        long long ht_size = est_next * 2;
        if (ht_size < 1024) ht_size = 1024;
        /* Round up to power of 2 for fast modulo */
        long long ht_pow = 1;
        while (ht_pow < ht_size) ht_pow <<= 1;
        ht_size = ht_pow;
        u64 ht_mask_bits = ht_size - 1;

        fprintf(stderr, "  level %d->%d: %lld entries, ht_size=%lld (%.1f GB) (%lds)...",
                p, p+1, curr_count, ht_size,
                (double)ht_size * sizeof(ht_entry_t) / 1e9,
                time(NULL) - t0);
        fflush(stderr);

        ht_entry_t *ht = (ht_entry_t *)malloc(ht_size * sizeof(ht_entry_t));
        if (!ht) {
            fprintf(stderr, "\nFailed to allocate hash table (%.1f GB)\n",
                    (double)ht_size * sizeof(ht_entry_t) / 1e9);
            free(curr);
            return 1;
        }
        /* Initialize all keys to empty */
        memset(ht, 0xFF, ht_size * sizeof(ht_entry_t));
        /* But values need to be 0 — handle in insert */

        long long next_count = 0;

        for (long long i = 0; i < curr_count; i++) {
            unsigned int mask = curr[i].mask;
            int v = curr[i].v;
            u128 cnt = curr[i].val;

            for (int u = 0; u < n; u++) {
                if (mask & (1u << u)) continue;
                if (u == v - 1) continue;

                u128 weight = (u == v + 1) ? 2 * cnt : cnt;
                unsigned int new_mask = mask | (1u << u);
                u64 key = make_key(new_mask, u);

                /* Insert into hash table */
                u64 h = (key * 0x9E3779B97F4A7C15ULL) >> (64 - 30);
                h &= ht_mask_bits;
                while (1) {
                    if (ht[h].key == HT_EMPTY) {
                        ht[h].key = key;
                        ht[h].val = weight;
                        next_count++;
                        break;
                    } else if (ht[h].key == key) {
                        ht[h].val += weight;
                        break;
                    }
                    h = (h + 1) & ht_mask_bits;
                }
            }
        }

        fprintf(stderr, " %lld entries\n", next_count);

        /* Convert hash table to sorted entry array */
        free(curr);
        curr = (entry_t *)malloc(next_count * sizeof(entry_t));
        if (!curr) {
            fprintf(stderr, "Failed to allocate entry array (%lld entries)\n", next_count);
            free(ht);
            return 1;
        }

        long long idx = 0;
        for (long long h = 0; h < ht_size; h++) {
            if (ht[h].key != HT_EMPTY) {
                u64 key = ht[h].key;
                curr[idx].mask = (unsigned int)(key >> 5);
                curr[idx].v = (unsigned char)(key & 0x1F);
                curr[idx].val = ht[h].val;
                idx++;
            }
        }
        curr_count = next_count;
        free(ht);
    }

    /* Sum all values at the final level (should be full mask) */
    u128 W = 0;
    for (long long i = 0; i < curr_count; i++) {
        W += curr[i].val;
    }

    printf("W(%d) = ", n);
    print_u128(W);
    printf("\n");

    fprintf(stderr, "Total time: %lds\n", time(NULL) - t0);

    free(curr);
    return 0;
}
