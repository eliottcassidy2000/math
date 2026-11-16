/* verify the k-apex transfer against direct subset DP on random tapes */
#include <stdio.h>
#include <stdlib.h>
#include "w44_kapex.h"
static ull dpH(int n, unsigned *out) {
    static ull dp[1 << 20][20]; int N = 1 << n;
    for (int S = 0; S < N; S++) memset(dp[S], 0, sizeof(ull) * n);
    for (int v = 0; v < n; v++) dp[1 << v][v] = 1;
    for (int S = 1; S < N; S++) for (int v = 0; v < n; v++) { ull c = dp[S][v]; if (!c) continue;
        unsigned t = out[v] & ~S; while (t) { int u = __builtin_ctz(t); t &= t - 1; dp[S | 1u << u][u] += c; } }
    ull s = 0; for (int v = 0; v < n; v++) s += dp[N - 1][v]; return s;
}
int main(void) {
    srand(44); int bad = 0, tests = 0;
    for (int k = 1; k <= 4; k++) {
        K_ = k; kapex_init_perms(k);
        for (int trial = 0; trial < 60; trial++) {
            int m = 2 + rand() % (k == 4 ? 12 : 14); M_ = m; int n = m + k; unsigned out[24] = {0};
            for (int i = 0; i < m; i++) for (int j = i + 1; j < m; j++) out[i] |= 1u << j;
            for (int i = 0; i < m; i++) { beat_[i] = rand() & ((1 << k) - 1);
                for (int a = 0; a < k; a++) { if (beat_[i] >> a & 1) out[i] |= 1u << (m + a); else out[m + a] |= 1u << i; } }
            for (int a = 0; a < k; a++) for (int b = a + 1; b < k; b++) { int o = rand() & 1; arc_[a][b] = o; arc_[b][a] = !o;
                if (o) out[m + a] |= 1u << (m + b); else out[m + b] |= 1u << (m + a); }
            ull h1 = kapex_count(), h2 = dpH(n, out); tests++;
            if (h1 != h2) { bad++; if (bad < 5) printf("MISMATCH k=%d m=%d transfer %llu dp %llu\n", k, m, h1, h2); }
        }
        printf("k = %d apexes: %d permutations of the visiting order\n", k, nperm_);
    }
    printf("k-apex transfer vs subset DP: %d tests, %d mismatches\n", tests, bad);
    return bad != 0;
}
