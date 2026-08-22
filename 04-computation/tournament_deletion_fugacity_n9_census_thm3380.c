/*
 * Exact full order-nine accelerator for proposed THM-3380.
 *
 * Input: one 36-character upper-triangle 0/1 encoding per tournament
 * isomorphism class, with bit 1 meaning i->j for i<j.  The frozen audit used
 * the 191,536 default-output representatives of nauty-gentourng 2.7r3+ds-1;
 * input SHA-256 4f7d6c43cfed87e1e5293dc751736efe2d7bc1554946cdc83f4026a575fbbbf8.
 *
 * It computes Dham from Hamiltonian-path DP on all induced subsets, counts
 * directed Hamiltonian 9-cycles by an independent anchored DP, counts
 * unordered directed-triangle factors, checks Dham(0)=2*a9+8*b333 on every
 * class, and compares the Dham and (Dham,b333) collision partitions.
 *
 * Reproduce:
 *   gcc -O3 -std=c11 -Wall -Wextra -Werror SOURCE.c -o census
 *   ./census gentourng9.txt
 * Normal frozen transcript SHA-256:
 *   7a5a3605c6a0d8f47c9119564181f5ab603f9cf2eb0a13c244b8a13c78eda80f
 * The transcript also byte-matches builds at -O0 and -fsanitize=undefined.
 *
 * This accelerator establishes a FINITE-EXACT census, not an all-order
 * reconstruction theorem.  The generator-free Python companion independently
 * checks the frozen split witness by direct directed-cycle enumeration.
 */

#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define N 9
#define FULL ((1U << N) - 1U)
#define TABLE_SIZE (1U << 19)

typedef struct {
    uint8_t used;
    uint16_t min_b, max_b;
    uint32_t first_index, count;
    uint32_t first_a;
    int32_t d[N + 1];
    char first_bits[37];
} DEntry;

typedef struct {
    uint8_t used;
    uint16_t b;
    uint32_t count;
    int32_t d[N + 1];
} DBEntry;

static DEntry *d_table;
static DBEntry *db_table;
static uint32_t choose_table[N + 1][N + 1];

static uint64_t mix64(uint64_t x) {
    x ^= x >> 30;
    x *= UINT64_C(0xbf58476d1ce4e5b9);
    x ^= x >> 27;
    x *= UINT64_C(0x94d049bb133111eb);
    x ^= x >> 31;
    return x;
}

static uint64_t hash_d(const int32_t d[N + 1], uint16_t b, int include_b) {
    uint64_t h = UINT64_C(0x6a09e667f3bcc909);
    for (int i = 0; i <= N; ++i)
        h = mix64(h ^ (uint32_t)d[i] ^ ((uint64_t)i << 48));
    if (include_b) h = mix64(h ^ b ^ UINT64_C(0x9e3779b97f4a7c15));
    return h;
}

static int equal_d(const int32_t left[N + 1], const int32_t right[N + 1]) {
    return memcmp(left, right, sizeof(int32_t) * (N + 1)) == 0;
}

static void parse_tournament(const char *line, uint16_t out[N]) {
    memset(out, 0, sizeof(uint16_t) * N);
    int edge = 0;
    for (int left = 0; left < N; ++left) {
        for (int right = left + 1; right < N; ++right) {
            if (line[edge++] == '1') out[left] |= (uint16_t)(1U << right);
            else out[right] |= (uint16_t)(1U << left);
        }
    }
}

static int cyclic_triangle(const uint16_t out[N], uint16_t mask) {
    for (int vertex = 0; vertex < N; ++vertex) {
        if (!(mask & (1U << vertex))) continue;
        if (__builtin_popcount((unsigned)(out[vertex] & mask)) != 1) return 0;
    }
    return 1;
}

static uint16_t triangle_factor_count(const uint16_t out[N]) {
    uint16_t triangles[84];
    int number = 0;
    for (int a = 0; a < N; ++a)
        for (int b = a + 1; b < N; ++b)
            for (int c = b + 1; c < N; ++c) {
                uint16_t mask = (uint16_t)((1U << a) | (1U << b) | (1U << c));
                if (cyclic_triangle(out, mask)) triangles[number++] = mask;
            }
    uint32_t ordered_factors = 0;
    for (int i = 0; i < number; ++i) {
        for (int j = 0; j < number; ++j) {
            if (i == j) continue;
            if (triangles[i] & triangles[j]) continue;
            uint16_t third = (uint16_t)(FULL ^ triangles[i] ^ triangles[j]);
            if (cyclic_triangle(out, third)) ++ordered_factors;
        }
    }
    if (ordered_factors % 6U) {
        fprintf(stderr, "triangle factor orbit failure: %u\n", ordered_factors);
        exit(6);
    }
    return (uint16_t)(ordered_factors / 6U);
}

static void tournament_packet(
    const uint16_t out[N], int32_t diagonal[N + 1],
    uint32_t *cycle9, uint16_t *factor333
) {
    uint32_t dp[1U << N][N];
    uint32_t anchored[1U << N][N];
    uint32_t hamiltonian[1U << N];
    memset(dp, 0, sizeof(dp));
    memset(anchored, 0, sizeof(anchored));
    memset(hamiltonian, 0, sizeof(hamiltonian));
    for (int vertex = 0; vertex < N; ++vertex) dp[1U << vertex][vertex] = 1;
    anchored[1][0] = 1;

    for (uint16_t mask = 1; mask <= FULL; ++mask) {
        for (int endpoint = 0; endpoint < N; ++endpoint) {
            uint16_t available = (uint16_t)(out[endpoint] & (FULL ^ mask));
            uint32_t count = dp[mask][endpoint];
            uint32_t anchor_count = anchored[mask][endpoint];
            while (available) {
                uint16_t bit = (uint16_t)(available & (uint16_t)(-available));
                available ^= bit;
                int successor = __builtin_ctz((unsigned)bit);
                dp[mask | bit][successor] += count;
                anchored[mask | bit][successor] += anchor_count;
            }
        }
        for (int endpoint = 0; endpoint < N; ++endpoint)
            hamiltonian[mask] += dp[mask][endpoint];
    }
    hamiltonian[0] = 1;

    int64_t accumulator[N + 1] = {0};
    for (uint16_t remaining = 0; remaining <= FULL; ++remaining) {
        int deleted = N - __builtin_popcount((unsigned)remaining);
        for (int degree = 0; degree <= deleted; ++degree) {
            int sign = ((deleted - degree) & 1) ? -1 : 1;
            accumulator[degree] +=
                (int64_t)sign * hamiltonian[remaining] * choose_table[deleted][degree];
        }
    }
    for (int degree = 0; degree <= N; ++degree) {
        if (accumulator[degree] < INT32_MIN || accumulator[degree] > INT32_MAX) {
            fprintf(stderr, "coefficient overflow\n");
            exit(2);
        }
        diagonal[degree] = (int32_t)accumulator[degree];
    }

    *cycle9 = 0;
    for (int endpoint = 1; endpoint < N; ++endpoint)
        if (out[endpoint] & 1U) *cycle9 += anchored[FULL][endpoint];
    *factor333 = triangle_factor_count(out);
    if (diagonal[0] != (int32_t)(2U * *cycle9 + 8U * *factor333)) {
        fprintf(stderr, "full-cover identity failed: %d %u %u\n",
                diagonal[0], *cycle9, *factor333);
        exit(3);
    }
}

static DEntry *insert_d(
    const int32_t d[N + 1], uint16_t b, uint32_t a,
    uint32_t index, const char *bits, int *is_new
) {
    uint32_t slot = (uint32_t)hash_d(d, 0, 0) & (TABLE_SIZE - 1);
    while (d_table[slot].used && !equal_d(d_table[slot].d, d))
        slot = (slot + 1) & (TABLE_SIZE - 1);
    DEntry *entry = &d_table[slot];
    if (!entry->used) {
        entry->used = 1;
        entry->min_b = entry->max_b = b;
        entry->first_index = index;
        entry->first_a = a;
        entry->count = 1;
        memcpy(entry->d, d, sizeof(entry->d));
        memcpy(entry->first_bits, bits, 37);
        *is_new = 1;
    } else {
        if (b < entry->min_b) entry->min_b = b;
        if (b > entry->max_b) entry->max_b = b;
        ++entry->count;
        *is_new = 0;
    }
    return entry;
}

static DBEntry *insert_db(const int32_t d[N + 1], uint16_t b, int *is_new) {
    uint32_t slot = (uint32_t)hash_d(d, b, 1) & (TABLE_SIZE - 1);
    while (db_table[slot].used &&
           (db_table[slot].b != b || !equal_d(db_table[slot].d, d)))
        slot = (slot + 1) & (TABLE_SIZE - 1);
    DBEntry *entry = &db_table[slot];
    if (!entry->used) {
        entry->used = 1;
        entry->b = b;
        entry->count = 1;
        memcpy(entry->d, d, sizeof(entry->d));
        *is_new = 1;
    } else {
        ++entry->count;
        *is_new = 0;
    }
    return entry;
}

int main(int argc, char **argv) {
    if (argc != 2) {
        fprintf(stderr, "usage: %s gentourng9.txt\n", argv[0]);
        return 1;
    }
    for (int n = 0; n <= N; ++n) {
        choose_table[n][0] = choose_table[n][n] = 1;
        for (int k = 1; k < n; ++k)
            choose_table[n][k] = choose_table[n - 1][k - 1] + choose_table[n - 1][k];
    }
    d_table = calloc(TABLE_SIZE, sizeof(*d_table));
    db_table = calloc(TABLE_SIZE, sizeof(*db_table));
    if (!d_table || !db_table) return 4;

    FILE *input = fopen(argv[1], "r");
    if (!input) return 5;
    char bits[128];
    uint32_t index = 0, distinct_d = 0, distinct_db = 0;
    int have_witness = 0;
    uint32_t witness_left = 0, witness_right = 0;
    uint32_t witness_left_a = 0, witness_right_a = 0;
    uint16_t witness_left_b = 0, witness_right_b = 0;
    int32_t witness_d[N + 1] = {0};
    char witness_left_bits[37] = {0}, witness_right_bits[37] = {0};
    uint64_t stream_digest = UINT64_C(0xcbf29ce484222325);

    while (fgets(bits, sizeof(bits), input)) {
        bits[strcspn(bits, "\r\n")] = '\0';
        if (strlen(bits) != 36) continue;
        uint16_t out[N];
        int32_t d[N + 1];
        uint32_t a;
        uint16_t b;
        parse_tournament(bits, out);
        tournament_packet(out, d, &a, &b);
        int new_d, new_db;
        DEntry *entry = insert_d(d, b, a, index, bits, &new_d);
        insert_db(d, b, &new_db);
        distinct_d += (uint32_t)new_d;
        distinct_db += (uint32_t)new_db;
        if (!have_witness && !new_d &&
            (b < entry->max_b || b > entry->min_b)) {
            have_witness = 1;
            witness_left = entry->first_index;
            witness_right = index;
            witness_left_a = entry->first_a;
            witness_right_a = a;
            witness_left_b =
                (uint16_t)((d[0] - 2U * entry->first_a) / 8U);
            witness_right_b = b;
            memcpy(witness_d, d, sizeof(d));
            memcpy(witness_left_bits, entry->first_bits, 37);
            memcpy(witness_right_bits, bits, 37);
        }
        stream_digest = mix64(stream_digest ^ hash_d(d, b, 1) ^ a ^ index);
        ++index;
    }
    fclose(input);

    uint32_t d_nontrivial = 0, d_colliding = 0, d_largest = 1, d_split = 0;
    uint32_t db_nontrivial = 0, db_colliding = 0, db_largest = 1;
    for (uint32_t slot = 0; slot < TABLE_SIZE; ++slot) {
        if (d_table[slot].used) {
            uint32_t count = d_table[slot].count;
            if (count > 1) { ++d_nontrivial; d_colliding += count; }
            if (count > d_largest) d_largest = count;
            if (d_table[slot].min_b != d_table[slot].max_b) ++d_split;
        }
        if (db_table[slot].used) {
            uint32_t count = db_table[slot].count;
            if (count > 1) { ++db_nontrivial; db_colliding += count; }
            if (count > db_largest) db_largest = count;
        }
    }
    printf("classes=%u\n", index);
    printf("D: distinct=%u nontrivial_fibers=%u colliding_classes=%u largest=%u\n",
           distinct_d, d_nontrivial, d_colliding, d_largest);
    printf("D_b333: distinct=%u nontrivial_fibers=%u colliding_classes=%u largest=%u\n",
           distinct_db, db_nontrivial, db_colliding, db_largest);
    printf("D_fibers_split_by_b333=%u\n", d_split);
    printf("stream_digest_fnv_mix=%016" PRIx64 "\n", stream_digest);
    if (have_witness) {
        printf("first_split_indices=%u,%u\n", witness_left, witness_right);
        printf("left_bits=%s full=%u*x+%u*x^3\n",
               witness_left_bits, witness_left_a, witness_left_b);
        printf("right_bits=%s full=%u*x+%u*x^3\n",
               witness_right_bits, witness_right_a, witness_right_b);
        printf("common_D=(");
        for (int i = 0; i <= N; ++i)
            printf("%s%d", i ? "," : "", witness_d[i]);
        printf(")\n");
    } else {
        printf("first_split_indices=NONE\n");
    }
    free(d_table);
    free(db_table);
    return 0;
}
