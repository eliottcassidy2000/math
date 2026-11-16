/* digits.c <n>: for every tournament class of order n (gentourng, stdin) and every vertex t: h = H(T), e_t = Hamiltonian
   paths ending at t, B_t = sum over A subset of V - t of h(T[A]) s_t(T[V - A]) (s_t = Hamiltonian paths of the subset starting
   at t), digit sigma_t = B_t - h (the height-two tower has count 3h + 2 sigma_t, THM-4865), charge a_t = B_t - e_t - h + 1
   (the rooted append has count 2h - 1 + a_t); also E_t (even root-start sum) to check B_t - e_t = 2 E_t.
   Output per count h: whether digits 0, 1, 2 and charge 0 occur, the set of small digits and charges (|.| <= 8), counts. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64; typedef long long i64;
static int n; static u64 adj[12];
static i64 Hm[1<<10]; static i64 dp[1<<10][10], dps[1<<10][10];
#define MAXH 4000
static int seen[MAXH][3], seenA0[MAXH], attained[MAXH], seenReady[MAXH], seenReadyA0[MAXH]; static long digitsets[MAXH][17], chargesets[MAXH][17]; static long mism = 0;
int main(int argc, char **argv) {
    int nn = atoi(argv[1]); char line[256]; long classes = 0;
    while (fgets(line, sizeof line, stdin)) {
        int len = strlen(line); while (len && (line[len-1]=='\n' || line[len-1]=='\r')) line[--len] = 0; if (!len) continue; classes++;
        n = nn; for (int i = 0; i < n; i++) adj[i] = 0; int k = 0;
        for (int i = 0; i < n; i++) for (int j = i+1; j < n; j++) { if (line[k++] == '1') adj[i] |= 1ULL<<j; else adj[j] |= 1ULL<<i; }
        int N = 1 << n;
        for (int m = 0; m < N; m++) { Hm[m] = 0; for (int l = 0; l < n; l++) dp[m][l] = 0; }
        for (int v = 0; v < n; v++) dp[1<<v][v] = 1;
        for (int m = 1; m < N; m++) for (int l = 0; l < n; l++) { i64 a = dp[m][l]; if (!a) continue; u64 nx = adj[l] & ~(u64)m; while (nx) { int v = __builtin_ctzll(nx); nx &= nx-1; dp[m|(1<<v)][v] += a; } }
        for (int m = 1; m < N; m++) { i64 s = 0; for (int l = 0; l < n; l++) s += dp[m][l]; Hm[m] = s; } Hm[0] = 1;
        i64 h = Hm[N-1]; if (h >= MAXH) continue; attained[h] = 1;
        for (int t = 0; t < n; t++) {
            /* paths starting at t: dps[m][l] = paths covering m from t ending at l */
            for (int m = 0; m < N; m++) for (int l = 0; l < n; l++) dps[m][l] = 0;
            dps[1<<t][t] = 1;
            for (int m = 1; m < N; m++) { if (!((m>>t)&1)) continue; for (int l = 0; l < n; l++) { i64 a = dps[m][l]; if (!a) continue; u64 nx = adj[l] & ~(u64)m; while (nx) { int v = __builtin_ctzll(nx); nx &= nx-1; dps[m|(1<<v)][v] += a; } } }
            i64 e = dp[N-1][t];     /* paths ending at t */
            i64 B = 0, E = 0, O = 0;
            for (int m = 1; m < N; m++) if ((m>>t)&1) { i64 st = 0; for (int l = 0; l < n; l++) st += dps[m][l]; i64 w = st * Hm[(N-1) & ~m]; B += w; if (__builtin_popcount(m) % 2 == 0) E += w; else O += w; }
            /* B_t = sum over A subset V-t of h(A) s_t(V-A) = sum over m containing t of s_t(m) h(complement) -- same as above */
            i64 sigma = B - h, a = B - e - h + 1; i64 g = Hm[(N-1) & ~(1<<t)]; int ready = (O - E - g == 0); if (ready) seenReady[h] = 1; if (ready && a == 0) seenReadyA0[h] = 1;
            if (B - e != 2*E) mism++;
            if (sigma >= 0 && sigma <= 2) seen[h][sigma] = 1;
            if (a == 0) seenA0[h] = 1;
            if (sigma >= -8 && sigma <= 8) digitsets[h][sigma+8]++;
            if (a >= -8 && a <= 8) chargesets[h][a+8]++;
        }
    }
    printf("order %d: %ld classes; mismatches of B - e = 2E: %ld\n", nn, classes, mism);
    printf("h: digit0 digit1 digit2 charge0 ready readyCharge0 | small digits (-8..8) | small charges\n");
    for (int h = 1; h < MAXH; h++) if (attained[h]) {
        printf("  %d: %d %d %d %d %d %d |", h, seen[h][0], seen[h][1], seen[h][2], seenA0[h], seenReady[h], seenReadyA0[h]);
        for (int d = -8; d <= 8; d++) if (digitsets[h][d+8]) printf(" %d", d);
        printf(" |");
        for (int d = -8; d <= 8; d++) if (chargesets[h][d+8]) printf(" %d", d);
        printf("\n");
    }
    return 0;
}
