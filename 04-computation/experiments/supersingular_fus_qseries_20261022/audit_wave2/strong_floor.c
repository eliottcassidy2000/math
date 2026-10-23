/* audit_wave2: exhaustive strong floor f(m) = min H(T) over strong labelled tournaments on m vertices,
   m = 3..7, and the set of H values attained by ALL tournaments on m <= 7 vertices (checks THM-1370's
   f(3..7) = 3, 5, 9, 15, 25 and that 7, 21 are never attained up to 7 vertices).
   H(T) = number of Hamiltonian paths (Held-Karp DP over subsets). Strong = every vertex reaches every vertex.
   Build: cc -O2 -o strong_floor strong_floor.c ; Run: ./strong_floor > strong_floor.out */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int adj[8];                /* adj[u] = bitmask of v with u -> v */
static long long dp[1 << 8][8];

static long long hcount(int m) {
    int full = (1 << m) - 1;
    for (int s = 0; s <= full; s++) for (int v = 0; v < m; v++) dp[s][v] = 0;
    for (int v = 0; v < m; v++) dp[1 << v][v] = 1;
    for (int s = 1; s <= full; s++)
        for (int v = 0; v < m; v++) {
            long long c = dp[s][v];
            if (!c) continue;
            int out = adj[v] & ~s;
            while (out) { int w = __builtin_ctz(out); out &= out - 1; dp[s | (1 << w)][w] += c; }
        }
    long long h = 0;
    for (int v = 0; v < m; v++) h += dp[full][v];
    return h;
}

static int strong(int m) {
    int full = (1 << m) - 1;
    for (int st = 0; st < m; st++) {
        int seen = 1 << st, fr = seen;
        while (fr) { int nx = 0, f = fr; while (f) { int u = __builtin_ctz(f); f &= f - 1; nx |= adj[u]; } nx &= ~seen; seen |= nx; fr = nx; }
        if (seen != full) return 0;
    }
    return 1;
}

int main(void) {
    for (int m = 3; m <= 7; m++) {
        int pairs[28][2], np = 0;
        for (int u = 0; u < m; u++) for (int v = u + 1; v < m; v++) { pairs[np][0] = u; pairs[np][1] = v; np++; }
        long long fmin = -1, hmax = 0;
        static char seenH[2000];
        memset(seenH, 0, sizeof seenH);
        for (long long code = 0; code < (1LL << np); code++) {
            for (int u = 0; u < m; u++) adj[u] = 0;
            for (int i = 0; i < np; i++) {
                int u = pairs[i][0], v = pairs[i][1];
                if (code >> i & 1) adj[u] |= 1 << v; else adj[v] |= 1 << u;
            }
            long long h = hcount(m);
            if (h < 2000) seenH[h] = 1;
            if (h > hmax) hmax = h;
            if (strong(m) && (fmin < 0 || h < fmin)) fmin = h;
        }
        printf("m = %d: strong floor f(m) = %lld; max H = %lld; H = 7 attained: %s; H = 21 attained: %s\n",
               m, fmin, hmax, seenH[7] ? "YES" : "no", seenH[21] ? "YES" : "no");
        printf("   odd values <= 61 attained:");
        for (int h = 1; h <= 61; h += 2) if (seenH[h]) printf(" %d", h);
        printf("\n");
    }
    return 0;
}
