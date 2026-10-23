/* Audit (HYP-9371 section 8): first-principles cluster kernel of a core with prescribed path lengths, written from the
 * THM-4757 section 4 definition only (no code from decomp.c / loop_engine*.py):
 *
 *   K(core, lambda) = sum over run decompositions of the subdivided core (core edge e = a path of lambda_e arcs, with
 *   its reference orientation tail -> head) of  mu * prod_runs ((-1)^alpha + (-1)^(l - alpha)),
 *
 *   * core vertex of degree d with a set of k disjoint half-edge pairs (pass-throughs): mu factor (-1)^(s-1) (s-1)!,
 *     s = d - k (block size);
 *   * internal (degree-2) node: pass-through (factor 1) or junction of two run ends (block of size 2, factor -1);
 *   * runs = components of the arc-end link graph (arc links + pass-through links); a closed component excludes the
 *     configuration; l = #arcs of the run, alpha = #arcs traversed against the reference in one direction.
 *
 * Then [z^L] X_excess = sum_cores (1/|Aut|) sum_(lambda: sum = L) K(core, lambda) hom(core; Shat^lambda), the hom number
 * taken with each core edge in its reference orientation (computed elsewhere).
 * Input (stdin): V E, E lines "a b", m, then m lines of E path lengths.  Output: one line per length vector, K.
 * Build: cc -O2 -o a5_kernel a5_kernel.c */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXV 8
#define MAXE 16
#define MAXARC 64
#define MAXEND (2 * MAXARC)

static int V, E, ea[MAXE], eb[MAXE], lam[MAXE];
static int narc, abase[MAXE];
/* arc-end id = 2*arc + side (side 0 at the arc's reference start, side 1 at its end) */
static int vend[MAXV][2 * MAXE], nvend[MAXV];
/* partial matchings per core vertex: stored as flat arrays of pairs */
typedef struct { int np; int p[MAXE][2]; } Match;
static Match *mt[MAXV]; static int nmt[MAXV], cap[MAXV];

static void addm(int v, int np, int p[][2]) {
    if (nmt[v] == cap[v]) { cap[v] = cap[v] ? 2 * cap[v] : 64; mt[v] = realloc(mt[v], cap[v] * sizeof(Match)); }
    mt[v][nmt[v]].np = np; memcpy(mt[v][nmt[v]].p, p, sizeof(int) * 2 * np); nmt[v]++;
}
static void genm(int v, int *it, int k, int np, int p[][2]) {
    if (k == 0) { addm(v, np, p); return; }
    int h = it[0], rest[2 * MAXE];
    for (int i = 1; i < k; i++) rest[i - 1] = it[i];
    genm(v, rest, k - 1, np, p);
    for (int i = 0; i < k - 1; i++) {
        int r2[2 * MAXE], m = 0;
        for (int j = 0; j < k - 1; j++) if (j != i) r2[m++] = rest[j];
        p[np][0] = h; p[np][1] = rest[i];
        genm(v, r2, m, np + 1, p);
    }
}
static int64_t mob(int s) { int64_t f = 1; for (int i = 2; i <= s - 1; i++) f *= i; return ((s - 1) % 2) ? -f : f; }

int main(void) {
    if (scanf("%d %d", &V, &E) != 2) return 1;
    for (int e = 0; e < E; e++) if (scanf("%d %d", &ea[e], &eb[e]) != 2) return 1;
    int m; if (scanf("%d", &m) != 1) return 1;
    for (int q = 0; q < m; q++) {
        for (int e = 0; e < E; e++) if (scanf("%d", &lam[e]) != 1) return 1;
        narc = 0;
        for (int e = 0; e < E; e++) { abase[e] = narc; narc += lam[e]; }
        if (narc > MAXARC) { fprintf(stderr, "too many arcs\n"); return 1; }
        for (int v = 0; v < V; v++) { nvend[v] = 0; nmt[v] = 0; }
        for (int e = 0; e < E; e++) {
            vend[ea[e]][nvend[ea[e]]++] = 2 * abase[e];                       /* first arc, start side */
            vend[eb[e]][nvend[eb[e]]++] = 2 * (abase[e] + lam[e] - 1) + 1;    /* last arc, end side */
        }
        int p[MAXE][2];
        for (int v = 0; v < V; v++) genm(v, vend[v], nvend[v], 0, p);
        /* internal nodes: list of (end side of arc j-1, start side of arc j) */
        int nint = 0, ia[MAXARC], ib[MAXARC];
        for (int e = 0; e < E; e++)
            for (int j = 1; j < lam[e]; j++) { ia[nint] = 2 * (abase[e] + j - 1) + 1; ib[nint] = 2 * (abase[e] + j); nint++; }
        int nends = 2 * narc;
        int64_t K = 0;
        int idx[MAXV]; memset(idx, 0, sizeof idx);
        int partner[MAXEND], vis[MAXEND];
        for (;;) {
            int64_t mu = 1;
            for (int x = 0; x < nends; x++) partner[x] = -1;
            for (int v = 0; v < V; v++) {
                Match *M = &mt[v][idx[v]];
                mu *= mob(nvend[v] - M->np);
                for (int j = 0; j < M->np; j++) { partner[M->p[j][0]] = M->p[j][1]; partner[M->p[j][1]] = M->p[j][0]; }
            }
            for (int t = 0; t < (1 << nint); t++) {          /* bit i set: internal node i is a junction */
                int64_t w = mu;
                for (int i = 0; i < nint; i++) {
                    if ((t >> i) & 1) { partner[ia[i]] = -1; partner[ib[i]] = -1; w = -w; }
                    else { partner[ia[i]] = ib[i]; partner[ib[i]] = ia[i]; }
                }
                for (int x = 0; x < nends; x++) vis[x] = 0;
                int dead = 0;
                for (int x0 = 0; x0 < nends && !dead; x0++) {
                    if (vis[x0] || partner[x0] >= 0) continue;   /* start only at run ends */
                    int x = x0, l = 0, al = 0;
                    for (;;) {
                        vis[x] = 1;
                        int y = x ^ 1;                            /* traverse the arc */
                        vis[y] = 1; l++;
                        if ((x & 1) == 1) al++;                   /* entered at the end side: against the reference */
                        if (partner[y] < 0) break;                /* run ends */
                        x = partner[y];
                    }
                    if (l & 1) { dead = 1; break; }               /* (-1)^al + (-1)^(l-al) = 0 */
                    w *= (al & 1) ? -2 : 2;
                }
                if (dead) continue;
                int closed = 0;
                for (int x = 0; x < nends; x++) if (!vis[x]) { closed = 1; break; }
                if (closed) continue;
                K += w;
            }
            int v = 0;
            while (v < V && ++idx[v] == nmt[v]) { idx[v] = 0; v++; }
            if (v == V) break;
        }
        printf("%lld\n", (long long)K);
    }
    return 0;
}
