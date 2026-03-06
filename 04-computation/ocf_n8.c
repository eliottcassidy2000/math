/*
 * n=8 OCF verification: H(T) = I(Omega(T), 2) for all 2^27 tournaments.
 *
 * Compile: gcc -O3 -o ocf_n8 ocf_n8.c -lpthread
 * Run: ./ocf_n8
 *
 * Instance: opus-2026-03-05-S4
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <pthread.h>
#include <unistd.h>

#define N 8
#define FULL ((1 << N) - 1)
#define N_ARC_PAIRS 27
#define TOTAL_CONFIGS (1 << 27)

/* Arc pairs (a,b) with a<b, excluding (0,1) */
static int arc_a[N_ARC_PAIRS], arc_b[N_ARC_PAIRS];

void init_arc_pairs() {
    int idx = 0;
    for (int a = 0; a < N; a++)
        for (int b = a + 1; b < N; b++)
            if (!(a == 0 && b == 1))
                { arc_a[idx] = a; arc_b[idx] = b; idx++; }
}

/* Build tournament from bitmask */
void build_T(int mask, int T[N][N]) {
    memset(T, 0, sizeof(int) * N * N);
    T[0][1] = 1;
    for (int bit = 0; bit < N_ARC_PAIRS; bit++) {
        int a = arc_a[bit], b = arc_b[bit];
        if (mask & (1 << bit)) {
            T[a][b] = 1;
        } else {
            T[b][a] = 1;
        }
    }
}

/* Hamiltonian path count via bitmask DP */
long long ham_count(int T[N][N]) {
    /* dp[mask][last] = # Ham paths on mask ending at last */
    long long dp[1 << N][N];
    memset(dp, 0, sizeof(dp));
    for (int v = 0; v < N; v++)
        dp[1 << v][v] = 1;

    for (int mask = 1; mask <= FULL; mask++) {
        for (int last = 0; last < N; last++) {
            long long c = dp[mask][last];
            if (c == 0) continue;
            for (int nxt = 0; nxt < N; nxt++) {
                if (mask & (1 << nxt)) continue;
                if (T[last][nxt])
                    dp[mask | (1 << nxt)][nxt] += c;
            }
        }
    }

    long long total = 0;
    for (int v = 0; v < N; v++)
        total += dp[FULL][v];
    return total;
}

/*
 * Count directed Hamiltonian cycles in T[combo] where combo has m vertices.
 * combo[0..m-1] are the vertex indices.
 * Fix combo[0] as start, count directed paths through all others back to start.
 */
long long count_cycles(int T[N][N], int combo[], int m) {
    if (m < 3) return 0;
    int start = combo[0];
    int others[N];
    int mo = m - 1;
    for (int i = 0; i < mo; i++) others[i] = combo[i + 1];
    int ofull = (1 << mo) - 1;

    /* dp[mask][last_idx] = # paths from start through subset of others */
    long long dp[1 << 7][7];  /* max mo = N-1 = 7 */
    memset(dp, 0, sizeof(dp));

    for (int i = 0; i < mo; i++)
        if (T[start][others[i]])
            dp[1 << i][i] = 1;

    for (int mask = 1; mask <= ofull; mask++) {
        for (int li = 0; li < mo; li++) {
            long long c = dp[mask][li];
            if (c == 0) continue;
            for (int ni = 0; ni < mo; ni++) {
                if (mask & (1 << ni)) continue;
                if (T[others[li]][others[ni]])
                    dp[mask | (1 << ni)][ni] += c;
            }
        }
    }

    long long total = 0;
    for (int li = 0; li < mo; li++)
        if (T[others[li]][start])
            total += dp[ofull][li];
    return total;
}

/*
 * Compute I(Omega(T), 2).
 * At n=8, max independent set size = 2.
 * I = 1 + 2*total_cycles + 4*vd_pairs
 */
long long compute_ocf(int T[N][N]) {
    /* Find all odd-length vertex subsets with directed cycles */
    int vmasks[256];
    long long vcounts[256];
    int nsets = 0;

    int combo[N];
    /* Enumerate odd-length subsets */
    for (int length = 3; length <= N; length += 2) {
        /* Generate all C(N, length) subsets */
        /* Use bitmask enumeration */
        for (int smask = 0; smask <= FULL; smask++) {
            if (__builtin_popcount(smask) != length) continue;

            int m = 0;
            for (int v = 0; v < N; v++)
                if (smask & (1 << v))
                    combo[m++] = v;

            long long cnt = count_cycles(T, combo, m);
            if (cnt > 0) {
                vmasks[nsets] = smask;
                vcounts[nsets] = cnt;
                nsets++;
            }
        }
    }

    long long total_cycles = 0;
    for (int i = 0; i < nsets; i++)
        total_cycles += vcounts[i];

    long long vd_pairs = 0;
    for (int i = 0; i < nsets; i++)
        for (int j = i + 1; j < nsets; j++)
            if ((vmasks[i] & vmasks[j]) == 0)
                vd_pairs += vcounts[i] * vcounts[j];

    return 1 + 2 * total_cycles + 4 * vd_pairs;
}

/* Thread data */
typedef struct {
    int start_mask, end_mask;
    long long checked, fails;
    int thread_id;
} thread_data_t;

void *verify_chunk(void *arg) {
    thread_data_t *td = (thread_data_t *)arg;
    td->checked = 0;
    td->fails = 0;

    for (int mask = td->start_mask; mask < td->end_mask; mask++) {
        int T[N][N];
        build_T(mask, T);

        long long ht = ham_count(T);
        long long it = compute_ocf(T);

        if (ht != it) {
            td->fails++;
            if (td->fails <= 3)
                printf("  FAIL: thread %d, mask %d, H=%lld, I=%lld\n",
                       td->thread_id, mask, ht, it);
        }
        td->checked++;
    }
    return NULL;
}

int main() {
    init_arc_pairs();

    int n_threads = 4;
    int total = TOTAL_CONFIGS;
    int chunk_per_thread = total / n_threads;

    printf("=== n=8 OCF Exhaustive Proof (C, %d threads) ===\n", n_threads);
    printf("Total: %d configs\n", total);

    /* Quick test first */
    printf("Quick test (1000 configs)...\n");
    clock_t t0 = clock();
    long long test_fails = 0;
    for (int mask = 0; mask < 1000; mask++) {
        int T[N][N];
        build_T(mask, T);
        long long ht = ham_count(T);
        long long it = compute_ocf(T);
        if (ht != it) test_fails++;
    }
    double test_time = (double)(clock() - t0) / CLOCKS_PER_SEC;
    printf("  1000 configs in %.2fs (%.0f/s), %lld fails\n",
           test_time, 1000.0/test_time, test_fails);
    double est = (double)total / (1000.0/test_time);
    printf("  Estimated total: %.0fs = %.1fhr (single thread)\n", est, est/3600);
    printf("  Estimated total: %.0fs = %.1fhr (%d threads)\n",
           est/n_threads, est/n_threads/3600, n_threads);

    if (test_fails > 0) {
        printf("Failures in quick test, aborting.\n");
        return 1;
    }

    printf("\nStarting full verification...\n");
    time_t wall_start = time(NULL);

    pthread_t threads[n_threads];
    thread_data_t td[n_threads];

    for (int i = 0; i < n_threads; i++) {
        td[i].thread_id = i;
        td[i].start_mask = i * chunk_per_thread;
        td[i].end_mask = (i == n_threads - 1) ? total : (i + 1) * chunk_per_thread;
        pthread_create(&threads[i], NULL, verify_chunk, &td[i]);
    }

    /* Monitor progress */
    long long total_checked = 0, total_fails = 0;
    while (total_checked < total) {
        sleep(30);
        total_checked = 0;
        total_fails = 0;
        for (int i = 0; i < n_threads; i++) {
            total_checked += td[i].checked;
            total_fails += td[i].fails;
        }
        double elapsed = difftime(time(NULL), wall_start);
        double rate = total_checked / elapsed;
        double eta = (total - total_checked) / rate;
        printf("  [%5.1f%%] %12lld/%d | %lld fails | %.0f/s | ETA %.2fhr\n",
               100.0 * total_checked / total, total_checked, total,
               total_fails, rate, eta / 3600);

        if (total_fails > 10) {
            printf("Too many failures, aborting.\n");
            for (int i = 0; i < n_threads; i++)
                pthread_cancel(threads[i]);
            return 1;
        }
    }

    for (int i = 0; i < n_threads; i++)
        pthread_join(threads[i], NULL);

    total_checked = 0;
    total_fails = 0;
    for (int i = 0; i < n_threads; i++) {
        total_checked += td[i].checked;
        total_fails += td[i].fails;
    }

    double wall_time = difftime(time(NULL), wall_start);
    printf("\n============================================================\n");
    printf("TOTAL: %lld checked, %lld failures, %.0fs\n",
           total_checked, total_fails, wall_time);

    if (total_fails == 0 && total_checked == total) {
        printf("\n*** PROVED: H(T) = I(Omega(T), 2) for ALL n=8 tournaments ***\n");
    }

    return 0;
}
