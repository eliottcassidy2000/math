/*
 * n=8 OCF exhaustive proof: H(T) = I(Omega(T), 2) for all 2^27 tournaments.
 * 4-threaded version based on verified single-threaded code.
 *
 * Compile: gcc -O3 -o ocf_n8_full ocf_n8_full.c -lpthread
 * Instance: opus-2026-03-05-S3
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <pthread.h>
#include <unistd.h>

#define N 8
#define FULL ((1 << N) - 1)
#define N_ARCS 27
#define TOTAL (1 << 27)

static int arc_a[N_ARCS], arc_b[N_ARCS];

void init_arcs() {
    int idx = 0;
    for (int a = 0; a < N; a++)
        for (int b = a+1; b < N; b++)
            if (!(a==0 && b==1))
                { arc_a[idx]=a; arc_b[idx]=b; idx++; }
}

long long ham_count(int T[N][N]) {
    long long dp[256][8];
    memset(dp, 0, sizeof(dp));
    for (int v = 0; v < N; v++) dp[1<<v][v] = 1;
    for (int mask = 1; mask <= FULL; mask++)
        for (int last = 0; last < N; last++) {
            long long c = dp[mask][last];
            if (!c) continue;
            for (int nxt = 0; nxt < N; nxt++) {
                if (mask & (1<<nxt)) continue;
                if (T[last][nxt]) dp[mask|(1<<nxt)][nxt] += c;
            }
        }
    long long t = 0;
    for (int v = 0; v < N; v++) t += dp[FULL][v];
    return t;
}

long long count_cycles(int T[N][N], int combo[], int m) {
    if (m < 3) return 0;
    int start = combo[0], mo = m-1;
    int others[7];
    for (int i = 0; i < mo; i++) others[i] = combo[i+1];
    int ofull = (1<<mo)-1;
    long long dp[128][7];
    memset(dp, 0, sizeof(dp));
    for (int i = 0; i < mo; i++)
        if (T[start][others[i]]) dp[1<<i][i] = 1;
    for (int mask = 1; mask <= ofull; mask++)
        for (int li = 0; li < mo; li++) {
            long long c = dp[mask][li];
            if (!c) continue;
            for (int ni = 0; ni < mo; ni++) {
                if (mask & (1<<ni)) continue;
                if (T[others[li]][others[ni]]) dp[mask|(1<<ni)][ni] += c;
            }
        }
    long long t = 0;
    for (int li = 0; li < mo; li++)
        if (T[others[li]][start]) t += dp[ofull][li];
    return t;
}

long long compute_ocf(int T[N][N]) {
    int vmasks[256]; long long vcounts[256]; int nsets = 0;
    int combo[N];
    for (int smask = 0; smask <= FULL; smask++) {
        int pc = __builtin_popcount(smask);
        if (pc < 3 || pc % 2 == 0) continue;
        int m = 0;
        for (int v = 0; v < N; v++)
            if (smask & (1<<v)) combo[m++] = v;
        long long cnt = count_cycles(T, combo, m);
        if (cnt > 0) { vmasks[nsets]=smask; vcounts[nsets]=cnt; nsets++; }
    }
    long long tc = 0;
    for (int i = 0; i < nsets; i++) tc += vcounts[i];
    long long vd = 0;
    for (int i = 0; i < nsets; i++)
        for (int j = i+1; j < nsets; j++)
            if ((vmasks[i] & vmasks[j]) == 0) vd += vcounts[i]*vcounts[j];
    return 1 + 2*tc + 4*vd;
}

typedef struct {
    int start_mask, end_mask;
    volatile long long checked, fails;
    int thread_id;
} thread_data_t;

void *verify_chunk(void *arg) {
    thread_data_t *td = (thread_data_t *)arg;
    td->checked = 0;
    td->fails = 0;
    for (int mask = td->start_mask; mask < td->end_mask; mask++) {
        int T[N][N];
        memset(T, 0, sizeof(T));
        T[0][1] = 1;
        for (int bit = 0; bit < N_ARCS; bit++) {
            if (mask & (1<<bit)) T[arc_a[bit]][arc_b[bit]] = 1;
            else T[arc_b[bit]][arc_a[bit]] = 1;
        }
        long long h = ham_count(T);
        long long i = compute_ocf(T);
        if (h != i) {
            td->fails++;
            if (td->fails <= 3)
                printf("  FAIL thread %d mask %d: H=%lld I=%lld\n", td->thread_id, mask, h, i);
        }
        td->checked++;
    }
    return NULL;
}

int main() {
    init_arcs();
    int n_threads = 4;
    int chunk = TOTAL / n_threads;

    printf("=== n=8 OCF Exhaustive Proof (%d threads) ===\n", n_threads);
    printf("Total: %d configs (2^27)\n", TOTAL);
    fflush(stdout);

    time_t wall_start = time(NULL);
    pthread_t threads[4];
    thread_data_t td[4];

    for (int i = 0; i < n_threads; i++) {
        td[i].thread_id = i;
        td[i].start_mask = i * chunk;
        td[i].end_mask = (i == n_threads-1) ? TOTAL : (i+1)*chunk;
        pthread_create(&threads[i], NULL, verify_chunk, &td[i]);
    }

    /* Monitor */
    long long total_checked = 0;
    while (total_checked < TOTAL) {
        sleep(30);
        total_checked = 0;
        long long total_fails = 0;
        for (int i = 0; i < n_threads; i++) {
            total_checked += td[i].checked;
            total_fails += td[i].fails;
        }
        double elapsed = difftime(time(NULL), wall_start);
        double rate = elapsed > 0 ? total_checked / elapsed : 0;
        double eta = rate > 0 ? (TOTAL - total_checked) / rate : 0;
        printf("  [%5.1f%%] %12lld/%d | %lld fails | %.0f/s | ETA %.2fhr\n",
               100.0 * total_checked / TOTAL, total_checked, TOTAL,
               total_fails, rate, eta / 3600);
        fflush(stdout);

        if (total_fails > 10) {
            printf("Too many failures, aborting.\n");
            return 1;
        }
    }

    for (int i = 0; i < n_threads; i++)
        pthread_join(threads[i], NULL);

    long long final_checked = 0, final_fails = 0;
    for (int i = 0; i < n_threads; i++) {
        final_checked += td[i].checked;
        final_fails += td[i].fails;
    }

    double wall_time = difftime(time(NULL), wall_start);
    printf("\n============================================================\n");
    printf("TOTAL: %lld checked, %lld failures, %.0fs (%.2fhr)\n",
           final_checked, final_fails, wall_time, wall_time/3600);

    if (final_fails == 0 && final_checked == TOTAL)
        printf("\n*** PROVED: H(T) = I(Omega(T), 2) for ALL n=8 tournaments ***\n");

    return final_fails > 0 ? 1 : 0;
}
