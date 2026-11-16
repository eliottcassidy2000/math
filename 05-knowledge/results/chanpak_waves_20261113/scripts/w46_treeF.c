/* Tree-reversal tournaments.  Vertices 0..n-1 in build order; parent p[j] < j for j >= 1.
   Arcs: i -> j for i < j, except j -> p[j] (each new vertex beats only its parent).
   Mode "all N": enumerate all recursive trees with n <= N vertices, record the set of H values.
   Mode "read": read lines "n p1 p2 ... p_{n-1}" (parents of vertices 1..n-1, 0-based), print H. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long ull;
static ull *dp;
static ull Hcount(int n, const unsigned *out) {
    int N = 1 << n;
    memset(dp, 0, sizeof(ull) * (size_t)N * n);
    for (int v = 0; v < n; v++) dp[(size_t)(1u << v) * n + v] = 1;
    for (int S = 1; S < N; S++) {
        ull *row = dp + (size_t)S * n;
        for (int v = 0; v < n; v++) { ull c = row[v]; if (!c) continue;
            unsigned t = out[v] & ~(unsigned)S;
            while (t) { int u = __builtin_ctz(t); t &= t - 1; dp[(size_t)(S | 1u << u) * n + u] += c; } }
    }
    ull s = 0; for (int v = 0; v < n; v++) s += dp[(size_t)(N - 1) * n + v]; return s;
}
static void build(int n, const int *p, unsigned *out) {
    for (int i = 0; i < n; i++) out[i] = 0;
    for (int j = 1; j < n; j++) for (int i = 0; i < j; i++) {
        if (p[j] == i) out[j] |= 1u << i; else out[i] |= 1u << j; }
}
#define HMAX 4000000
static unsigned char seen[HMAX];
static int p[32]; static unsigned out[32]; static long cnt; static ull hmaxv;
static int MAXCH = 99, CATER = 0; static int nch[32];
static int okcat(int n) { /* caterpillar: removing leaves leaves a path */
    int deg[32] = {0}; for (int j = 1; j < n; j++) { deg[j]++; deg[p[j]]++; }
    int spine = 0, sp[32]; for (int v = 0; v < n; v++) if (deg[v] > 1 || n <= 2) sp[spine++] = v;
    for (int i = 0; i < spine; i++) { int d = 0; for (int k = 0; k < spine; k++) if (k != i && (p[sp[k]] == sp[i] || p[sp[i]] == sp[k])) d++; if (d > 2) return 0; }
    return 1; }
static void rec(int j, int n) {
    if (j == n) { if (CATER && !okcat(n)) return; build(n, p, out); ull h = Hcount(n, out); cnt++; if (h < HMAX) seen[h] = 1; if (h > hmaxv) hmaxv = h; return; }
    for (int q = 0; q < j; q++) { if (nch[q] >= MAXCH) continue; nch[q]++; p[j] = q; rec(j + 1, n); nch[q]--; }
}
int main(int argc, char **argv) {
    if (argc > 1 && !strcmp(argv[1], "all")) {
        int NM = atoi(argv[2]); if (argc > 3) MAXCH = atoi(argv[3]); if (argc > 4) CATER = atoi(argv[4]); dp = malloc(sizeof(ull) * ((size_t)1 << NM) * NM);
        for (int n = 1; n <= NM; n++) {
            cnt = 0; hmaxv = 0; rec(1, n);
            int miss = 0, firstmiss[40], nm = 0, distinct = 0; ull bound = 0;
            for (ull h = 1; h < HMAX && h <= hmaxv; h += 2) if (seen[h]) distinct++;
            /* largest B such that every odd h <= B except 7, 21 is seen */
            for (ull h = 1; h < HMAX; h += 2) { if (!seen[h] && h != 7 && h != 21) { bound = h - 2; break; } }
            for (ull h = 1; h < HMAX && nm < 12 && h <= hmaxv; h += 2) if (!seen[h]) firstmiss[nm++] = (int)h;
            printf("n<=%d: trees %ld (at n), max H %llu, odd values seen so far %d; complete (except 7,21) up to %llu; first misses:", n, cnt, hmaxv, distinct, bound);
            for (int i = 0; i < nm; i++) printf(" %d", firstmiss[i]); printf("\n"); fflush(stdout);
        }
        return 0;
    }
    dp = malloc(sizeof(ull) * ((size_t)1 << 21) * 21);
    int n; while (scanf("%d", &n) == 1) { p[0] = -1; for (int j = 1; j < n; j++) scanf("%d", &p[j]);
        build(n, p, out); printf("%llu\n", Hcount(n, out)); fflush(stdout); }
    return 0;
}
