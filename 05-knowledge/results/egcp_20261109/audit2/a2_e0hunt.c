/* AUDIT2: targeted local search for even graphs G on exactly n vertices with no universal vertex and d(G) = target.
 * Works in the complement H = comp(G): d(G) = 2^n - D(H), D(H) = #sets contained in some open neighbourhood of H
 * (empty set included; THM-5040 Lemma 0).  No universal vertex in G <=> no isolated vertex in H.
 * Simulated annealing on |D(H) - (2^n - target)| with single-edge flips; every exact hit is tested for evenness of G
 * by complete automorphism enumeration (sign on the edges of G); d(G) is then recomputed by the subset DP on G.
 * usage: a2_e0hunt n seed restarts steps target1 [target2 ...]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static int n;
static unsigned long long rs;
static unsigned long long rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static double urand(void) { return (rnd() >> 11) * (1.0 / 9007199254740992.0); }

static unsigned char mark[1 << 16];
static int Dcount(const int *h) {
    int full = 1 << n; memset(mark, 0, full);
    int cnt = 0;
    for (int v = 0; v < n; v++) {
        int N = h[v], sub = N;
        while (1) { if (!mark[sub]) { mark[sub] = 1; cnt++; } if (!sub) break; sub = (sub - 1) & N; }
    }
    if (!mark[0]) cnt++;
    return cnt;
}
static int dcount(const int *g) {
    int full = (1 << n) - 1, c = 0;
    static int dom[1 << 16]; dom[0] = 0;
    for (int S = 1; S <= full; S++) { int v = __builtin_ctz(S); dom[S] = dom[S & (S - 1)] | g[v] | (1 << v); if (dom[S] == full) c++; }
    return c;
}
/* evenness of G */
static int G[16], deg[16], ord[16], perm[16], used[16], oddf; static long long naut;
static void bt(int i) {
    if (oddf) return;
    if (i == n) { naut++; int rev = 0; for (int u = 0; u < n; u++) for (int v = u + 1; v < n; v++) if (((G[u] >> v) & 1) && perm[u] > perm[v]) rev++; if (rev & 1) oddf = 1; return; }
    int v = ord[i];
    for (int w = 0; w < n; w++) {
        if (used[w] || deg[w] != deg[v]) continue;
        int ok = 1; for (int j = 0; j < i; j++) { int u = ord[j]; if (((G[v] >> u) & 1) != ((G[w] >> perm[u]) & 1)) { ok = 0; break; } }
        if (!ok) continue;
        perm[v] = w; used[w] = 1; bt(i + 1); used[w] = 0; perm[v] = -1; if (oddf) return;
    }
}
static int is_even(const int *g) {
    memcpy(G, g, sizeof(int) * n);
    for (int v = 0; v < n; v++) { deg[v] = __builtin_popcount(G[v]); ord[v] = v; used[v] = 0; perm[v] = -1; }
    for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) if (deg[ord[b]] > deg[ord[a]]) { int t = ord[a]; ord[a] = ord[b]; ord[b] = t; }
    oddf = 0; naut = 0; bt(0); return !oddf;
}

int main(int argc, char **argv) {
    n = atoi(argv[1]); rs = 88172645463325252ULL ^ (unsigned long long)atoll(argv[2]) * 2654435761ULL;
    long restarts = atol(argv[3]), steps = atol(argv[4]);
    int nt = argc - 5; int tg[64], found[64] = {0}; long hits[64] = {0}, oddhits[64] = {0};
    for (int i = 0; i < nt; i++) tg[i] = atoi(argv[5 + i]);
    int full = (1 << n) - 1;
    for (long r = 0; r < restarts; r++) {
        int ti = r % nt; if (found[ti]) { int all = 1; for (int i = 0; i < nt; i++) all &= found[i]; if (all) break; continue; }
        int target = (1 << n) - tg[ti];             /* wanted D(H) */
        int h[16]; memset(h, 0, sizeof(h));
        double p = 0.15 + 0.5 * urand();
        for (int u = 0; u < n; u++) for (int v = u + 1; v < n; v++) if (urand() < p) { h[u] |= 1 << v; h[v] |= 1 << u; }
        int cur = abs(Dcount(h) - target);
        double T = 2.0 + 30.0 * urand();
        for (long s = 0; s < steps; s++) {
            int u = rnd() % n, v = rnd() % n; if (u == v) continue;
            h[u] ^= 1 << v; h[v] ^= 1 << u;
            int iso = 0; for (int w = 0; w < n; w++) if (!h[w]) iso = 1;
            int nd = iso ? 1 << 30 : abs(Dcount(h) - target);
            if (nd <= cur || urand() < exp((cur - nd) / T)) cur = nd; else { h[u] ^= 1 << v; h[v] ^= 1 << u; }
            T *= 0.9995; if (T < 0.05) T = 0.05;
            if (cur == 0) {
                int g[16]; for (int w = 0; w < n; w++) g[w] = (full ^ h[w]) & ~(1 << w);
                hits[ti]++;
                if (is_even(g)) {
                    int d = dcount(g); int univ = 0; for (int w = 0; w < n; w++) if (g[w] == (full ^ (1 << w))) univ = 1;
                    if (d == tg[ti] && !univ) {
                        found[ti] = 1;
                        printf("FOUND n=%d d=%d naut=%lld G-masks:", n, d, naut); for (int w = 0; w < n; w++) printf(" %d", g[w]); printf("\n"); fflush(stdout);
                        break;
                    }
                } else oddhits[ti]++;
                /* keep walking: perturb */
                int a = rnd() % n, b = rnd() % n; if (a != b) { h[a] ^= 1 << b; h[b] ^= 1 << a; }
                int iso2 = 0; for (int w = 0; w < n; w++) if (!h[w]) iso2 = 1;
                if (iso2) { h[a] ^= 1 << b; h[b] ^= 1 << a; }
                cur = abs(Dcount(h) - target);
            }
        }
    }
    for (int i = 0; i < nt; i++) printf("target %d: %s (exact-D hits %ld, of which odd %ld)\n", tg[i], found[i] ? "REALISED" : "not found", hits[i], oddhits[i]);
    return 0;
}
