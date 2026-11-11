/* AUDIT program (independent): the WHOLE switching class of (source) + T for a tournament T on n <= 15 vertices
 * (N = n + 1 <= 16): for every switching set X (a subset of the vertices of T; 2^n members, Gray code) the
 * numbers of directed cycles of every length (rooted subset programme, exact below 2^32) and the number of
 * Hamiltonian paths (exact, 64 bits). Tallies:
 *   - the parity vector of the odd cycle counts (c_3, c_5, ..., c_(N-1)) against the prediction of (bK),
 *     [N = 4 mod 8] C(N/2 - 1, k) mod 2;
 *   - H mod 4 against 1 + N/2 (bS), H mod 8, 1 + 2 c_odd = H (mod 4);
 *   - the parity vectors of the even cycle counts (how many different ones occur);
 *   - S^2 = -(N-1) I for every member (checked).
 * usage: aud_class < input      (input: n, then n rows of 0/1 for T)
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define MAXN 16
static int N;
static int adj[MAXN][MAXN];

static uint32_t off[(1 << 15) + 1];
static uint32_t dp[15 * (1 << 14)];
static uint64_t hp[1 << MAXN][MAXN];

static void cycles(uint64_t *cyc) {
    int r, i, j, l;
    for (l = 0; l <= N; l++) cyc[l] = 0;
    for (r = 0; r <= N - 3; r++) {
        int m = N - 1 - r;
        uint32_t o[MAXN], start = 0, back = 0, full = (1u << m) - 1, mask;
        for (i = 0; i < m; i++) {
            o[i] = 0;
            for (j = 0; j < m; j++) if (adj[r + 1 + i][r + 1 + j]) o[i] |= 1u << j;
            if (adj[r][r + 1 + i]) start |= 1u << i;
            if (adj[r + 1 + i][r]) back |= 1u << i;
        }
        memset(dp, 0, (size_t)off[full + 1] * sizeof(uint32_t));
        for (i = 0; i < m; i++) if (start >> i & 1) dp[off[1u << i]] = 1;
        for (mask = 1; mask <= full; mask++) {
            uint32_t rem = ~mask & full, mm, idx = 0;
            const uint32_t *rowp = dp + off[mask];
            int k = __builtin_popcount(mask);
            for (mm = mask; mm; mm &= mm - 1, idx++) {
                uint32_t x = rowp[idx];
                if (!x) continue;
                int u = __builtin_ctz(mm);
                if (k >= 2 && (back >> u & 1)) cyc[k + 1] += x;
                uint32_t t;
                for (t = o[u] & rem; t; t &= t - 1) {
                    int w = __builtin_ctz(t);
                    uint32_t nm = mask | (1u << w);
                    dp[off[nm] + (uint32_t)__builtin_popcount(nm & ((1u << w) - 1))] += x;
                }
            }
        }
    }
}

static uint64_t hampaths(void) {
    uint32_t full = (1u << N) - 1, s, outm[MAXN];
    int v, w;
    for (v = 0; v < N; v++) { outm[v] = 0; for (w = 0; w < N; w++) if (adj[v][w]) outm[v] |= 1u << w; }
    for (s = 0; s <= full; s++) for (v = 0; v < N; v++) hp[s][v] = 0;
    for (v = 0; v < N; v++) hp[1u << v][v] = 1;
    for (s = 1; s <= full; s++) {
        uint32_t mm;
        for (mm = s; mm; mm &= mm - 1) {
            v = __builtin_ctz(mm);
            uint64_t x = hp[s][v];
            if (!x) continue;
            uint32_t t;
            for (t = outm[v] & ~s & full; t; t &= t - 1) { w = __builtin_ctz(t); hp[s | (1u << w)][w] += x; }
        }
    }
    uint64_t h = 0;
    for (v = 0; v < N; v++) h += hp[full][v];
    return h;
}

int main(void) {
    int n, i, j, k;
    char row[64];
    if (scanf("%d", &n) != 1 || n < 3 || n > 15) { fprintf(stderr, "bad n\n"); return 1; }
    N = n + 1;
    memset(adj, 0, sizeof adj);
    for (i = 0; i < n; i++) {
        if (scanf("%63s", row) != 1) return 1;
        for (j = 0; j < n; j++) adj[i + 1][j + 1] = row[j] == '1';
    }
    for (j = 1; j < N; j++) adj[0][j] = 1;      /* vertex 0 is the source */
    uint32_t mask;
    off[0] = 0;
    for (mask = 0; mask < (1u << (N - 1)); mask++) off[mask + 1] = off[mask] + (uint32_t)__builtin_popcount(mask);
    int m = N / 2 - 1;
    /* prediction of (bK) */
    long long binom[20][20];
    for (i = 0; i < 20; i++) { binom[i][0] = 1; for (j = 1; j < 20; j++) binom[i][j] = i ? binom[i-1][j-1] + binom[i-1][j] : 0; }
    uint32_t pred = 0;
    for (k = 1; k <= m; k++) if ((N % 8 == 4) && (binom[m][k] & 1)) pred |= 1u << (k - 1);
    uint64_t members = (uint64_t)1 << n, g;
    long long badK = 0, badS = 0, badHc = 0, notSH = 0;
    long long h8[8];
    memset(h8, 0, sizeof h8);
    /* distinct even parity vectors: at most 2^7 */
    static long long evencount[1 << 8];
    static long long oddcount[1 << 8];
    memset(evencount, 0, sizeof evencount);
    memset(oddcount, 0, sizeof oddcount);
    for (g = 0; g < members; g++) {
        if (g) {
            int w = __builtin_ctzll(g) + 1;      /* switch at the single vertex w */
            for (j = 0; j < N; j++) if (j != w) { adj[w][j] ^= 1; adj[j][w] ^= 1; }
        }
        /* skew Hadamard check */
        int ok = 1;
        for (i = 0; i < N && ok; i++) for (j = i + 1; j < N && ok; j++) {
            int d = 0;
            for (k = 0; k < N; k++) if (k != i && k != j) d += (adj[i][k] == adj[j][k]) ? 1 : -1;
            if (d != 0) ok = 0;
        }
        if (!ok) notSH++;
        uint64_t cyc[MAXN + 2];
        cycles(cyc);
        uint32_t vec = 0, evec = 0;
        uint64_t codd = 0;
        for (k = 1; k <= m; k++) { if (cyc[2 * k + 1] & 1) vec |= 1u << (k - 1); codd += cyc[2 * k + 1]; }
        for (k = 2; 2 * k <= N; k++) if (cyc[2 * k] & 1) evec |= 1u << (k - 2);
        oddcount[vec]++;
        evencount[evec]++;
        if (vec != pred) badK++;
        uint64_t H = hampaths();
        h8[H & 7]++;
        if ((H & 3) != (uint64_t)((1 + N / 2) & 3)) badS++;
        if ((H & 3) != ((1 + 2 * codd) & 3)) badHc++;
    }
    printf("N = %d: switching class of (source) + T, T of order %d: %llu switching sets\n", N, n, (unsigned long long)members);
    printf("  members that are not skew Hadamard: %lld\n", notSH);
    printf("  (bK) predicted parity vector (c_3, c_5, ...):");
    for (k = 1; k <= m; k++) printf("%d", (pred >> (k - 1)) & 1);
    printf("; members with another vector: %lld\n", badK);
    printf("  odd parity vectors found:");
    for (mask = 0; mask < (1u << m); mask++) if (oddcount[mask]) {
        printf(" ");
        for (k = 1; k <= m; k++) printf("%d", (mask >> (k - 1)) & 1);
        printf(":%lld", oddcount[mask]);
    }
    printf("\n  (bS) H = 1 + N/2 = %d (mod 4): violations %lld; H = 1 + 2 c_odd (mod 4) violations: %lld\n", (1 + N / 2) & 3, badS, badHc);
    printf("  H mod 8:");
    for (k = 0; k < 8; k++) if (h8[k]) printf(" %d:%lld", k, h8[k]);
    int distinct = 0;
    for (mask = 0; mask < (1u << 8); mask++) if (evencount[mask]) distinct++;
    printf("\n  distinct parity vectors of the even cycle counts (c_4, c_6, ..., c_N): %d\n", distinct);
    return 0;
}
