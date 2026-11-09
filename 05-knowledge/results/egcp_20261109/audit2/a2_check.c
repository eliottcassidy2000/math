/* AUDIT2 graph checker (independent code).  Reads graphs from stdin, one per line:
 *     n q a_0 a_1 ... a_{n-1}
 * where a_v is the neighbour bitmask of v (decimal) and q is a vertex to be fixed (-1: none).
 * Prints per line:  d  sign_status  q_status  has_universal  min_degree  max_degree
 *   d            = number of dominating sets (subset DP, exact)
 *   sign_status  = ODD if some automorphism reverses an odd number of edges (reference orientation u<v), else EVEN
 *   q_status     = QODD if some automorphism fixing q is odd, QEVEN if none, '-' if q = -1
 * Automorphisms are enumerated by plain backtracking with a degree filter; the search stops at the first odd
 * automorphism (separately for the q-stabiliser), so EVEN answers come from a complete enumeration.
 * Also checks the input is a simple undirected graph (symmetric, loopless).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXN 24
static int n, adj[MAXN], deg[MAXN], ord[MAXN], perm[MAXN], used[MAXN], fixq;
static int found_odd;
static long long naut;

static void bt(int i) {
    if (found_odd) return;
    if (i == n) {
        naut++;
        int rev = 0;
        for (int u = 0; u < n; u++) for (int v = u + 1; v < n; v++)
            if (((adj[u] >> v) & 1) && perm[u] > perm[v]) rev++;
        if (rev & 1) found_odd = 1;
        return;
    }
    int v = ord[i];
    for (int w = 0; w < n; w++) {
        if (used[w] || deg[w] != deg[v]) continue;
        if (fixq >= 0 && ((v == fixq) != (w == fixq))) continue;
        int ok = 1;
        for (int j = 0; j < i; j++) {
            int u = ord[j];
            if (((adj[v] >> u) & 1) != ((adj[w] >> perm[u]) & 1)) { ok = 0; break; }
        }
        if (!ok) continue;
        perm[v] = w; used[w] = 1;
        bt(i + 1);
        used[w] = 0; perm[v] = -1;
        if (found_odd) return;
    }
}

static int odd_search(int q) {
    fixq = q; found_odd = 0; naut = 0;
    for (int v = 0; v < n; v++) { used[v] = 0; perm[v] = -1; ord[v] = v; }
    /* order: q first if fixed, then by degree descending, then BFS-ish by adjacency to earlier vertices */
    int placed = 0, inord[MAXN] = {0};
    if (q >= 0) { ord[placed++] = q; inord[q] = 1; }
    while (placed < n) {
        int best = -1, bscore = -1;
        for (int v = 0; v < n; v++) if (!inord[v]) {
            int conn = 0; for (int j = 0; j < placed; j++) if ((adj[v] >> ord[j]) & 1) conn++;
            int score = conn * 64 + deg[v];
            if (score > bscore) { bscore = score; best = v; }
        }
        ord[placed++] = best; inord[best] = 1;
    }
    bt(0);
    return found_odd;
}

int main(void) {
    int q;
    while (scanf("%d %d", &n, &q) == 2) {
        for (int v = 0; v < n; v++) scanf("%d", &adj[v]);
        int bad = 0;
        for (int v = 0; v < n; v++) {
            if ((adj[v] >> v) & 1) bad = 1;
            for (int u = 0; u < n; u++) if (((adj[v] >> u) & 1) != ((adj[u] >> v) & 1)) bad = 1;
            if (adj[v] >> n) bad = 1;
            deg[v] = __builtin_popcount(adj[v]);
        }
        int full = (1 << n) - 1;
        long long d = 0;
        int *dom = malloc(sizeof(int) << n);
        dom[0] = 0;
        for (int S = 1; S <= full; S++) {
            int v = __builtin_ctz(S);
            dom[S] = dom[S & (S - 1)] | adj[v] | (1 << v);
            if (dom[S] == full) d++;
        }
        free(dom);
        int univ = 0, mind = 99, maxd = -1;
        for (int v = 0; v < n; v++) { if (adj[v] == (full ^ (1 << v))) univ = 1; if (deg[v] < mind) mind = deg[v]; if (deg[v] > maxd) maxd = deg[v]; }
        int odd = odd_search(-1);
        const char *qs = "-";
        if (q >= 0) qs = odd_search(q) ? "QODD" : "QEVEN";
        printf("%lld %s %s %d %d %d%s\n", d, odd ? "ODD" : "EVEN", qs, univ, mind, maxd, bad ? " BADINPUT" : "");
        fflush(stdout);
    }
    return 0;
}
