/*
 * h21_no_good_del_n10.c — Check for no-good-deletion cycle-rich tournaments at n=10.
 * Key question: at n >= 10, does every cycle-rich tournament have a good deletion?
 *
 * Author: opus-2026-03-07-S43
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define N 10

static unsigned int rng_state = 99999;
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

static inline int is_cycle_rich(const unsigned short out[N]) {
    unsigned short all = (1 << N) - 1;
    for (int i = 0; i < N; i++) {
        int s = __builtin_popcount(out[i]);
        if (s == 0 || s == N-1) return 0;
    }
    for (int v = 0; v < N; v++) {
        unsigned short ov = out[v];
        unsigned short iv = (~ov) & all & ~(1 << v);
        int found = 0;
        unsigned short ov_copy = ov;
        while (ov_copy && !found) {
            int a = __builtin_ctz(ov_copy);
            ov_copy &= ov_copy - 1;
            if (out[a] & iv & ~(1 << a))
                found = 1;
        }
        if (!found) return 0;
    }
    return 1;
}

static inline int has_3cycle_avoiding(const unsigned short out[N], int u, int v) {
    unsigned short all = (1 << N) - 1;
    unsigned short mask = all & ~(1 << v) & ~(1 << u);
    unsigned short ou = out[u] & mask;
    unsigned short iu = (~out[u]) & all & mask;
    unsigned short ou_copy = ou;
    while (ou_copy) {
        int a = __builtin_ctz(ou_copy);
        ou_copy &= ou_copy - 1;
        if (out[a] & iu) return 1;
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

int main() {
    setvbuf(stdout, NULL, _IONBF, 0);
    setvbuf(stderr, NULL, _IONBF, 0);

    long long trials = 200000000LL;
    long long cycle_rich = 0;
    long long no_good_del = 0;

    unsigned short out[N];

    for (long long t = 0; t < trials; t++) {
        random_tournament(out);
        if (!is_cycle_rich(out)) continue;
        cycle_rich++;

        int has_good = 0;
        for (int v = 0; v < N; v++) {
            if (deletion_is_cycle_rich(out, v)) { has_good = 1; break; }
        }

        if (!has_good) {
            no_good_del++;
            int sc[N];
            for (int i = 0; i < N; i++) sc[i] = __builtin_popcount(out[i]);
            printf("NO GOOD DEL #%lld: scores=", no_good_del);
            for (int i = 0; i < N; i++) printf("%d ", sc[i]);
            printf("\n");
        }

        if ((t + 1) % 50000000 == 0) {
            fprintf(stderr, "t=%lldM, cr=%lld, no_good=%lld\n",
                    (t+1)/1000000, cycle_rich, no_good_del);
        }
    }

    printf("\n=== n=%d RESULTS ===\n", N);
    printf("Trials: %lld\nCycle-rich: %lld\nNo good deletion: %lld\n",
           trials, cycle_rich, no_good_del);

    return 0;
}
