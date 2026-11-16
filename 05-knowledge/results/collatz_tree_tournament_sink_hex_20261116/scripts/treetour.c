/* treetour.c Nmax order rule [inflate]: the Collatz inverse tree (root 1; children of n: 2n, and (n-1)/3 when n = 4 mod 6,
   (n-1)/3 > 1), built one node at a time in a chosen order; the tree tournament T_N on the first N nodes: comparable pairs
   (ancestor, descendant) oriented ancestor -> descendant; incomparable pairs oriented by the rule:
     rule 0: later-built beats earlier-built;  rule 1: earlier beats later (transitive if the order is a linear extension);
     rule 2: the numerically larger label beats the smaller; rule 3: the smaller beats the larger.
   order 0: BFS by height, evens before odds within a level, then by value; order 1: BFS, odds first; order 2: by value
   (not a linear extension); order 3: DFS (even child first).  Prints H(T_N) for N = 1..Nmax (Nmax <= 20).
   With inflate = c > 1, every node is replaced by a transitive block of c vertices (block vertices in a chain) and H of the
   inflated tournament on the first N nodes (c N <= 20) is printed too. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64; typedef long long i64;
static int n; static u64 adj[24];
static i64 *dp; 
static i64 hamcount(void) {
    int N = 1 << n; dp = calloc((size_t)N * n, sizeof(i64));
    for (int v = 0; v < n; v++) dp[(size_t)(1<<v)*n + v] = 1;
    for (int m = 1; m < N; m++) for (int l = 0; l < n; l++) { i64 a = dp[(size_t)m*n + l]; if (!a) continue; u64 nx = adj[l] & ~(u64)m; while (nx) { int v = __builtin_ctzll(nx); nx &= nx-1; dp[(size_t)(m|(1<<v))*n + v] += a; } }
    i64 s = 0; for (int l = 0; l < n; l++) s += dp[(size_t)(N-1)*n + l]; free(dp); return s;
}
#define MAXN 64
static long long lab[MAXN]; static int par[MAXN], depth[MAXN], cnt = 0;
static int anc(int a, int b) { /* is a an ancestor of b */ while (b >= 0) { if (b == a) return 1; b = par[b]; } return 0; }
int main(int argc, char **argv) {
    int Nmax = atoi(argv[1]), order = atoi(argv[2]), rule = atoi(argv[3]); int inflate = argc > 4 ? atoi(argv[4]) : 1;
    /* generate the tree by BFS with the chosen tie-break, enough nodes */
    long long q[4096]; int qp[4096], qd[4096]; int qh = 0, qt = 0;
    q[qt] = 1; qp[qt] = -1; qd[qt] = 0; qt++;
    int want = Nmax + 2;
    while (qh < qt && cnt < want) {
        long long v = q[qh]; int p = qp[qh], d = qd[qh]; qh++;
        lab[cnt] = v; par[cnt] = p; depth[cnt] = d; int me = cnt; cnt++;
        long long kids[2]; int nk = 0;
        long long even = 2*v, odd = (v % 6 == 4 && (v-1)/3 > 1) ? (v-1)/3 : 0;
        if (order == 1) { if (odd) kids[nk++] = odd; kids[nk++] = even; } else { kids[nk++] = even; if (odd) kids[nk++] = odd; }
        for (int i = 0; i < nk; i++) { q[qt] = kids[i]; qp[qt] = me; qd[qt] = d+1; qt++; }
    }
    /* BFS by level with tie-break "evens first then value": our queue order is BFS with children appended in order; within
       a level nodes appear in the order their parents were processed, which is not exactly evens-first by level.  Re-sort
       each level: stable sort by (parity preference, value). */
    int idx[MAXN]; for (int i = 0; i < cnt; i++) idx[i] = i;
    if (order == 0 || order == 1) {
        for (int i = 0; i < cnt; i++) for (int j = i+1; j < cnt; j++) {
            int a = idx[i], b = idx[j]; int swap = 0;
            if (depth[a] != depth[b]) swap = depth[a] > depth[b];
            else { int pa = lab[a] & 1, pb = lab[b] & 1; if (order == 0) { if (pa != pb) swap = pa > pb; else swap = lab[a] > lab[b]; } else { if (pa != pb) swap = pa < pb; else swap = lab[a] > lab[b]; } }
            if (swap) { int t = idx[i]; idx[i] = idx[j]; idx[j] = t; }
        }
    } else if (order == 2) {
        for (int i = 0; i < cnt; i++) for (int j = i+1; j < cnt; j++) if (lab[idx[i]] > lab[idx[j]]) { int t = idx[i]; idx[i] = idx[j]; idx[j] = t; }
    } else if (order == 3) { /* DFS even child first */
        int st[MAXN], sp = 0, out = 0; st[sp++] = 0;
        while (sp) { int v = st[--sp]; idx[out++] = v; int kids[2], nk = 0; for (int u = cnt-1; u >= 0; u--) if (par[u] == v) kids[nk++] = u; /* push odd first so even pops first */
            for (int i = 0; i < nk; i++) st[sp++] = kids[i]; }
        /* nodes outside the generated set are ignored */ 
    }
    printf("# Collatz tree tournament: order=%d rule=%d inflate=%d; nodes in build order:", order, rule, inflate);
    for (int i = 0; i < Nmax; i++) printf(" %lld", lab[idx[i]]); printf("\n");
    for (int N = 1; N <= Nmax; N++) {
        n = N * inflate; if (n > 22) break;
        for (int i = 0; i < n; i++) adj[i] = 0;
        for (int i = 0; i < N; i++) for (int j = i+1; j < N; j++) {
            int a = idx[i], b = idx[j]; int a_beats_b;
            if (anc(a, b)) a_beats_b = 1; else if (anc(b, a)) a_beats_b = 0;
            else if (rule == 0) a_beats_b = 0; else if (rule == 1) a_beats_b = 1; else if (rule == 2) a_beats_b = lab[a] > lab[b]; else a_beats_b = lab[a] < lab[b];
            for (int s = 0; s < inflate; s++) for (int t = 0; t < inflate; t++) { int x = i*inflate+s, y = j*inflate+t; if (a_beats_b) adj[x] |= 1ULL<<y; else adj[y] |= 1ULL<<x; }
        }
        for (int i = 0; i < N; i++) for (int s = 0; s < inflate; s++) for (int t = s+1; t < inflate; t++) adj[i*inflate+s] |= 1ULL<<(i*inflate+t);
        printf("N=%d (node %lld) H=%lld\n", N, lab[idx[N-1]], hamcount());
    }
    return 0;
}
