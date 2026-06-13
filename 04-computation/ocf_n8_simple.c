/* Simple single-threaded n=8 OCF test */
#include <stdio.h>
#include <string.h>
#include <time.h>

#define N 8
#define FULL ((1 << N) - 1)

int arc_a[27], arc_b[27];

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

int main() {
    init_arcs();
    int test_count = 100000;
    printf("Testing %d configs...\n", test_count);
    fflush(stdout);
    clock_t t0 = clock();
    int fails = 0;
    for (int mask = 0; mask < test_count; mask++) {
        int T[N][N];
        memset(T, 0, sizeof(T));
        T[0][1] = 1;
        for (int bit = 0; bit < 27; bit++) {
            if (mask & (1<<bit)) T[arc_a[bit]][arc_b[bit]] = 1;
            else T[arc_b[bit]][arc_a[bit]] = 1;
        }
        long long h = ham_count(T);
        long long i = compute_ocf(T);
        if (h != i) {
            fails++;
            if (fails <= 3) printf("FAIL mask %d: H=%lld I=%lld\n", mask, h, i);
        }
    }
    double elapsed = (double)(clock() - t0) / CLOCKS_PER_SEC;
    double rate = test_count / elapsed;
    long long total = 1LL << 27;
    printf("Result: %d/%d fails in %.2fs (%.0f/s)\n", fails, test_count, elapsed, rate);
    printf("Full 2^27 estimate: %.0fs = %.2fhr (1 thread)\n", total/rate, total/rate/3600);
    printf("Full 2^27 estimate: %.0fs = %.2fhr (4 threads)\n", total/rate/4, total/rate/4/3600);
    return fails > 0 ? 1 : 0;
}
