/* procgen_selfie_20261001_euler.c  (procgen selfie lane, 2026-10-01)
 *
 * OPEN-Q-060 test: count unlabeled Euler (even) graphs F on n vertices whose automorphism group
 * preserves edge-orientation parity, i.e. eps_F(s) = (-1)^{#edges {a<b} of F with s(a) > s(b)} = +1
 * for every automorphism s.  Claim (twisted Brauer lemma): this count equals A049313(n).
 * Also checks, for every automorphism enumerated, the cycle form
 *     eps_F(s) = (-1)^{#even cycles of s (length 2k) whose antipodal pairs {x, s^k x} are edges of F}.
 *
 * usage: geng -q n | procgen_selfie_euler n      (reads graph6, prints a summary line to stdout)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXN 16
static int n, adj[MAXN][MAXN], deg[MAXN], sig[MAXN];
static int perm[MAXN], used[MAXN];
static long long naut, nneg, cycle_form_fail;
static int stop_early;

static int eps_of_perm(void) {
    int inv = 0;
    for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) if (adj[a][b] && perm[a] > perm[b]) inv++;
    return (inv & 1) ? -1 : 1;
}
static int eps_cycle_form(void) {
    int seen[MAXN] = {0}, cnt = 0;
    for (int s = 0; s < n; s++) {
        if (seen[s]) continue;
        int len = 0, x = s;
        do { seen[x] = 1; x = perm[x]; len++; } while (x != s);
        if (len % 2 == 0) {
            int y = s; for (int k = 0; k < len / 2; k++) y = perm[y];
            if (adj[s][y]) cnt++;
        }
    }
    return (cnt & 1) ? -1 : 1;
}
static void bt(int i) {
    if (stop_early) return;
    if (i == n) {
        naut++;
        int e = eps_of_perm();
        if (e != eps_cycle_form()) cycle_form_fail++;
        if (e < 0) { nneg++; stop_early = 1; }
        return;
    }
    for (int w = 0; w < n; w++) {
        if (used[w] || deg[w] != deg[i] || sig[w] != sig[i]) continue;
        int ok = 1;
        for (int j = 0; j < i && ok; j++) if (adj[i][j] != adj[w][perm[j]]) ok = 0;
        if (!ok) continue;
        used[w] = 1; perm[i] = w; bt(i + 1); used[w] = 0;
        if (stop_early) return;
    }
}

int main(int argc, char **argv) {
    n = atoi(argv[1]);
    char line[512];
    long long total = 0, even = 0, untwisted = 0;
    long long hist_edges_untw[200] = {0};
    while (fgets(line, sizeof line, stdin)) {
        line[strcspn(line, "\r\n")] = 0;
        if (!line[0]) continue;
        int nn = line[0] - 63;
        if (nn != n) continue;
        total++;
        memset(adj, 0, sizeof adj);
        int k = 0, pos = 1, bitsleft = 0, cur = 0;
        for (int j = 1; j < n; j++) for (int i = 0; i < j; i++) {
            if (bitsleft == 0) { cur = line[pos++] - 63; bitsleft = 6; }
            int b = (cur >> (bitsleft - 1)) & 1; bitsleft--;
            adj[i][j] = adj[j][i] = b; k++;
        }
        int ok = 1, m = 0;
        for (int v = 0; v < n; v++) { deg[v] = 0; for (int w = 0; w < n; w++) deg[v] += adj[v][w]; if (deg[v] & 1) ok = 0; m += deg[v]; }
        if (!ok) continue;
        even++;
        for (int v = 0; v < n; v++) { sig[v] = 0; for (int w = 0; w < n; w++) if (adj[v][w]) sig[v] += 1 << (deg[w] & 15); }
        naut = 0; nneg = 0; stop_early = 0;
        memset(used, 0, sizeof used);
        bt(0);
        if (nneg == 0) { untwisted++; hist_edges_untw[m / 2]++; }
    }
    printf("EULER n=%d graphs_read=%lld euler_graphs=%lld untwisted=%lld cycle_form_failures=%lld\n",
           n, total, even, untwisted, cycle_form_fail);
    printf("  untwisted by edge count:");
    for (int e = 0; e < 200; e++) if (hist_edges_untw[e]) printf(" %d:%lld", e, hist_edges_untw[e]);
    printf("\n");
    return 0;
}
