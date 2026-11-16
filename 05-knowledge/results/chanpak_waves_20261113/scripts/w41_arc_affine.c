/* # Wave 41 (mac-mini-2026-11-16-dials). Run in a directory holding lit_base.pkl from w39_lit_setup.py (slowtree) / masks files from gentourng (arc_affine). */
/* Wave 41: additive (affine) arc cubes.  For every class representative of order n and every k-subset of arcs
   (k = 2 or 3), compute the 2^k counts H(sigma) and the multilinear coefficients
       H(sigma) = c0 + sum c_i s_i + sum c_ij s_i s_j + c_123 s1 s2 s3   (s_i = 1 iff arc i is flipped).
   A cube is AFFINE if every interaction coefficient vanishes.  Print the histogram of the sorted absolute single-flip
   gaps (|c_1|, ..., |c_k|) over affine cubes, and of the interaction pattern over all cubes.
   usage: arc_affine n k < masks */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll;
static int n; static int out[10];
static ll dp[1 << 9][9];
static ll ham(void) {
    int N = 1 << n;
    for (int S = 0; S < N; S++) memset(dp[S], 0, sizeof(ll) * n);
    for (int v = 0; v < n; v++) dp[1 << v][v] = 1;
    for (int S = 1; S < N; S++)
        for (int v = 0; v < n; v++) {
            ll c = dp[S][v];
            if (c) { int t = out[v] & ~S; while (t) { int b = t & -t; dp[S | b][__builtin_ctz(b)] += c; t ^= b; } }
        }
    ll s = 0; for (int v = 0; v < n; v++) s += dp[N - 1][v]; return s;
}
static void flip(int a, int b) { if ((out[a] >> b) & 1) { out[a] &= ~(1 << b); out[b] |= 1 << a; } else { out[b] &= ~(1 << a); out[a] |= 1 << b; } }
#define MAXG 64
static long hist[MAXG][MAXG][MAXG];
int main(int argc, char **argv) {
    n = atoi(argv[1]); int k = atoi(argv[2]);
    int pa[40], pb[40], np = 0;
    for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) { pa[np] = a; pb[np] = b; np++; }
    char line[256]; long cubes = 0, affine = 0, pairfree = 0; int idx[3];
    while (fgets(line, sizeof line, stdin)) {
        char *p = line; int ok = 1;
        for (int i = 0; i < n; i++) { char *e; long v = strtol(p, &e, 10); if (e == p) { ok = 0; break; } out[i] = (int)v; p = e; }
        if (!ok) continue;
        int base[10]; memcpy(base, out, sizeof base);
        for (int i = 0; i < k; i++) idx[i] = i;
        while (1) {
            int K = 1 << k; ll val[8];
            for (int o = 0; o < K; o++) { memcpy(out, base, sizeof base); for (int i = 0; i < k; i++) if ((o >> i) & 1) flip(pa[idx[i]], pb[idx[i]]); val[o] = ham(); }
            cubes++;
            ll c[3]; for (int i = 0; i < k; i++) c[i] = val[1 << i] - val[0];
            int aff = 1, pf = 1;
            for (int i = 0; i < k; i++) for (int j = i + 1; j < k; j++) { ll cij = val[(1 << i) | (1 << j)] - val[1 << i] - val[1 << j] + val[0]; if (cij) { aff = 0; pf = 0; } }
            if (k == 3) { ll c3 = val[7] - val[3] - val[5] - val[6] + val[1] + val[2] + val[4] - val[0]; if (c3) aff = 0; }
            pairfree += pf;
            if (aff) {
                affine++;
                ll g[3]; for (int i = 0; i < k; i++) g[i] = llabs(c[i]) / 2;
                for (int i = 0; i < k; i++) for (int j = i + 1; j < k; j++) if (g[j] < g[i]) { ll t = g[i]; g[i] = g[j]; g[j] = t; }
                if (g[k - 1] < MAXG) hist[g[0]][k > 1 ? g[1] : 0][k > 2 ? g[2] : 0]++;
            }
            int i = k - 1; while (i >= 0 && idx[i] == np - k + i) i--; if (i < 0) break; idx[i]++; for (int j = i + 1; j < k; j++) idx[j] = idx[j - 1] + 1;
        }
    }
    printf("# order %d, k = %d: cubes %ld, pair-interaction-free %ld, affine %ld\n", n, k, cubes, pairfree, affine);
    printf("# affine cubes by sorted half-gaps (|c_i|/2): count\n");
    for (int a = 0; a < MAXG; a++) for (int b = 0; b < MAXG; b++) for (int d = 0; d < MAXG; d++) if (hist[a][b][d]) printf("%d,%d,%d: %ld\n", a, b, d, hist[a][b][d]);
    return 0;
}
