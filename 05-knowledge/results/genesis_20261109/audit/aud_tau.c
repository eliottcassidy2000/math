/* AUDIT G1 (own code): spanning-tree spectra of labelled simple graphs on m <= MMAX vertices.
 *  - tau by Kirchhoff (fraction-free Bareiss on the reduced Laplacian), validated against brute-force
 *    enumeration of (m-1)-edge subsets for m <= 5;
 *  - per m: connected labelled count, |S_tau(m)|, max, new values (not in S_tau(m-1)), labelled share of new;
 *  - alpha(N) = least m with N in S_tau(m), for N <= 300;
 *  - Lemma T1 on every 2-connected graph with m <= 7: tau >= m, equality iff cycle;
 *  - vertex insertion for m <= 6: tau(G) >= deg(v) tau(G-v) for every non-cut v, equality iff deg v = 1.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll;
#define MMAX 8
static int m, E, eu[40], ev[40];
static unsigned adj[MMAX];

static void build(unsigned long long mask) {
    for (int i = 0; i < m; i++) adj[i] = 0;
    for (int e = 0; e < E; e++) if (mask >> e & 1) { adj[eu[e]] |= 1u << ev[e]; adj[ev[e]] |= 1u << eu[e]; }
}
static int connected_on(unsigned vs) {   /* induced subgraph on vertex set vs connected? (vs nonempty) */
    if (!vs) return 1;
    unsigned seen = vs & -vs, frontier = seen;
    while (frontier) {
        unsigned nf = 0;
        for (unsigned f = frontier; f; f &= f - 1) { int v = __builtin_ctz(f); nf |= adj[v] & vs; }
        nf &= ~seen;
        seen |= nf;
        frontier = nf;
    }
    return seen == vs;
}
static ll tau_on(unsigned vs) {           /* spanning trees of induced subgraph on vs (Kirchhoff + Bareiss) */
    int idx[MMAX], k = 0;
    for (int v = 0; v < m; v++) if (vs >> v & 1) idx[k++] = v;
    if (k <= 1) return 1;
    ll a[MMAX][MMAX];
    int d = k - 1;
    for (int i = 0; i < d; i++)
        for (int j = 0; j < d; j++) {
            int u = idx[i + 1], w = idx[j + 1];
            if (i == j) a[i][j] = __builtin_popcount(adj[u] & vs);
            else a[i][j] = (adj[u] >> w & 1) ? -1 : 0;
        }
    ll prev = 1;
    int sign = 1;
    for (int c = 0; c < d; c++) {
        if (a[c][c] == 0) {
            int r = c + 1;
            while (r < d && a[r][c] == 0) r++;
            if (r == d) return 0;
            for (int j = 0; j < d; j++) { ll t = a[c][j]; a[c][j] = a[r][j]; a[r][j] = t; }
            sign = -sign;
        }
        for (int i = c + 1; i < d; i++)
            for (int j = c + 1; j < d; j++)
                a[i][j] = (a[i][j] * a[c][c] - a[i][c] * a[c][j]) / prev;
        prev = a[c][c];
    }
    return sign * a[d - 1][d - 1];
}
static ll brute_tau(unsigned long long mask) {   /* count (m-1)-subsets of edges forming a spanning tree */
    int es[40], ne = 0;
    for (int e = 0; e < E; e++) if (mask >> e & 1) es[ne++] = e;
    ll cnt = 0;
    for (unsigned long long s = 0; s < (1ULL << ne); s++) {
        if (__builtin_popcountll(s) != m - 1) continue;
        int par[MMAX];
        for (int i = 0; i < m; i++) par[i] = i;
        int ok = 1;
        for (int t = 0; t < ne && ok; t++) if (s >> t & 1) {
            int x = eu[es[t]], y = ev[es[t]];
            while (par[x] != x) x = par[x];
            while (par[y] != y) y = par[y];
            if (x == y) ok = 0; else par[x] = y;
        }
        cnt += ok;
    }
    return cnt;
}
int main(int argc, char **argv) {
    int MM = argc > 1 ? atoi(argv[1]) : MMAX;
    static unsigned char inprev[300000], incur[300000];
    static unsigned long long cnt[300000];
    int alpha[301];
    for (int i = 0; i <= 300; i++) alpha[i] = 0;
    memset(inprev, 0, sizeof inprev);
    inprev[1] = 1; /* m = 1: single vertex, tau = 1 */
    long t1bad = 0, t1tests = 0, insbad = 0, instests = 0, brbad = 0, brtests = 0;
    for (m = 2; m <= MM; m++) {
        E = 0;
        for (int i = 0; i < m; i++) for (int j = i + 1; j < m; j++) { eu[E] = i; ev[E] = j; E++; }
        memset(incur, 0, sizeof incur);
        memset(cnt, 0, sizeof cnt);
        unsigned full = (1u << m) - 1;
        unsigned long long nconn = 0;
        ll mx = 0;
        for (unsigned long long mask = 0; mask < (1ULL << E); mask++) {
            build(mask);
            if (!connected_on(full)) continue;
            nconn++;
            ll t = tau_on(full);
            if (t <= 0 || t >= 300000) { fprintf(stderr, "tau out of range %lld\n", t); return 1; }
            incur[t] = 1;
            cnt[t]++;
            if (t > mx) mx = t;
            if (m <= 5 && (mask % 7 == 0)) { brtests++; if (brute_tau(mask) != t) brbad++; }
            if (m <= 7 && m >= 3) {
                int two = 1;
                for (int v = 0; v < m && two; v++) if (!connected_on(full & ~(1u << v))) two = 0;
                if (two) {
                    t1tests++;
                    int cyc = 1;
                    for (int v = 0; v < m; v++) if (__builtin_popcount(adj[v]) != 2) cyc = 0;
                    if (t < m || ((t == m) != cyc)) t1bad++;
                }
            }
            if (m <= 6) {
                for (int v = 0; v < m; v++) {
                    unsigned rest = full & ~(1u << v);
                    if (!connected_on(rest)) continue;          /* v is a cut vertex */
                    instests++;
                    ll tv = tau_on(rest);
                    int dg = __builtin_popcount(adj[v]);
                    if (t < dg * tv) insbad++;
                    if ((t == dg * tv) != (dg == 1)) insbad++;
                    if ((t == tv) != (dg == 1)) insbad++;
                }
            }
        }
        long nvals = 0, newv = 0;
        unsigned long long lnew = 0;
        for (int t = 1; t < 300000; t++) if (incur[t]) {
            nvals++;
            if (!inprev[t]) { newv++; lnew += cnt[t]; if (t <= 300 && !alpha[t]) alpha[t] = m; }
        }
        printf("m=%d connected=%llu |S_tau|=%ld max=%lld new_values=%ld labelled_new=%llu share=%.4f two_in_S=%d\n",
               m, nconn, nvals, mx, newv, lnew, (double)lnew / nconn, incur[2]);
        fflush(stdout);
        /* S(m-1) <= S(m) check, then roll */
        for (int t = 1; t < 300000; t++) if (inprev[t] && !incur[t]) printf("  monotonicity FAIL at %d\n", t);
        memcpy(inprev, incur, sizeof incur);
    }
    printf("brute-force tau checks: %ld, mismatches %ld\n", brtests, brbad);
    printf("Lemma T1: %ld 2-connected labelled graphs (m=3..7), failures %ld\n", t1tests, t1bad);
    printf("vertex insertion: %ld (G, non-cut v) pairs (m<=6), failures %ld\n", instests, insbad);
    printf("alpha(N) for N<=60 (0 = not attained up to m=%d):\n", MM);
    for (int N = 1; N <= 60; N++) printf("%d:%d%s", N, alpha[N], N % 15 ? " " : "\n");
    printf("fixed points alpha(N)=N with N<=%d:", MM);
    for (int N = 3; N <= MM; N++) if (alpha[N] == N) printf(" %d", N);
    printf("\nnot attained up to m=%d among {9,10,13,18,21,22}:", MM);
    int chk[6] = {9, 10, 13, 18, 21, 22};
    for (int i = 0; i < 6; i++) printf(" %d:%s", chk[i], alpha[chk[i]] ? "attained" : "NO");
    printf("\n");
    return 0;
}
