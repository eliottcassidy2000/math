/* movegraph.c -- the native-move graph on isomorphism classes of n-tournaments (n <= 9).
 * Reads gentourng's text format on stdin.  For every class T and every directed Hamiltonian path
 * Q = (q_0..q_{n-1}) of T, the native move tau_Q(T) keeps the arcs of Q and reverses all other pairs;
 * Q is blue iff the positional reversal q_i -> q_{n-1-i} is an anti-automorphism of T.
 * Output: V lines  "V idx key scores strong sc aut H"  and  E lines  "E src dst blue count".
 * Build: cc -O2 -I/opt/homebrew/include/nauty movegraph.c -L/opt/homebrew/lib -lnauty -o movegraph */
#define MAXN 64
#include "nauty.h"
#include "nautinv.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
static int N;
static uint64_t canon_key(const uint32_t *q, int k, double *grp) {
    graph g[MAXN], cg[MAXN]; int lab[MAXN], ptn[MAXN], orbits[MAXN];
    DEFAULTOPTIONS_DIGRAPH(options); statsblk stats; options.getcanon = TRUE;
    EMPTYGRAPH(g, 1, k);
    for (int i = 0; i < k; i++) for (int j = 0; j < k; j++) if (q[i] >> j & 1) g[i] |= bit[j];
    densenauty(g, lab, ptn, orbits, &options, &stats, 1, k, cg);
    uint64_t key = 0; int bp = 0;
    for (int i = 0; i < k; i++) for (int j = i + 1; j < k; j++) { if (cg[i] & bit[j]) key |= (1ULL << bp); bp++; }
    if (grp) { double s = stats.grpsize1; for (int e = 0; e < stats.grpsize2; e++) s *= 10.0; *grp = s; }
    return key;
}
#define HS (1<<22)
static uint64_t hk[HS]; static int hv[HS]; static unsigned char hu[HS];
static int hfind(uint64_t v) { uint32_t h = (uint32_t)((v * 0x9E3779B97F4A7C15ULL) >> 42) & (HS-1); while (hu[h]) { if (hk[h] == v) return hv[h]; h = (h + 1) & (HS-1); } return -1; }
static void hput(uint64_t v, int idx) { uint32_t h = (uint32_t)((v * 0x9E3779B97F4A7C15ULL) >> 42) & (HS-1); while (hu[h]) h = (h + 1) & (HS-1); hu[h] = 1; hk[h] = v; hv[h] = idx; }
static uint32_t (*adj)[16]; static int nclasses = 0, cap = 0;
static int scc_count(const uint32_t *o, int n) { /* number of strong components: Landau: strong iff no proper initial segment of sorted scores satisfies sum = C(k,2) */
    int s[16]; for (int i = 0; i < n; i++) s[i] = __builtin_popcount(o[i]);
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) if (s[j] < s[i]) { int t = s[i]; s[i] = s[j]; s[j] = t; }
    int cnt = 0, sum = 0; for (int k = 1; k < n; k++) { sum += s[k-1]; if (sum == k*(k-1)/2) cnt++; } return cnt + 1;
}
/* per-class edge aggregation */
typedef struct { int dst; int blue; int count; } edge_t;
static edge_t edges[1<<16]; static int nedges;
static void add_edge(int dst, int blue) { for (int i = 0; i < nedges; i++) if (edges[i].dst == dst && edges[i].blue == blue) { edges[i].count++; return; } edges[nedges].dst = dst; edges[nedges].blue = blue; edges[nedges].count = 1; nedges++; }
static uint32_t cur[16]; static int path[16]; static long Hcount;
static void visit(int depth, uint32_t used) {
    if (depth == N) {
        Hcount++;
        uint32_t t[16]; for (int i = 0; i < N; i++) t[i] = 0;
        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) if (i != j && !(cur[i] >> j & 1)) t[i] |= 1u << j;   /* converse */
        for (int i = 0; i + 1 < N; i++) { int a = path[i], b = path[i+1]; t[a] |= 1u << b; t[b] &= ~(1u << a); }  /* restore Q */
        int pos[16]; for (int i = 0; i < N; i++) pos[path[i]] = i;
        int blue = 1;
        for (int u = 0; u < N && blue; u++) for (int v = u + 1; v < N; v++) {
            int ru = path[N-1-pos[u]], rv = path[N-1-pos[v]];
            int uv = (cur[u] >> v) & 1, rvu = (cur[rv] >> ru) & 1;
            if (uv != rvu) { blue = 0; break; }
        }
        uint64_t k = canon_key(t, N, NULL); int d = hfind(k); if (d < 0) { fprintf(stderr, "unknown class\n"); exit(3); }
        add_edge(d, blue); return;
    }
    int last = path[depth-1];
    for (int v = 0; v < N; v++) if (!(used >> v & 1) && (cur[last] >> v & 1)) { path[depth] = v; visit(depth + 1, used | (1u << v)); }
}
int main(int argc, char **argv) {
    N = atoi(argv[1]); char line[256];
    cap = 1 << 18; adj = malloc(sizeof(uint32_t[16]) * cap);
    static uint64_t keys[1<<18]; static double grps[1<<18];
    while (fgets(line, sizeof line, stdin)) {
        int L = strlen(line); while (L && (line[L-1] == '\n' || line[L-1] == '\r')) line[--L] = 0;
        if (L != N * (N - 1) / 2) continue;
        uint32_t o[16]; for (int i = 0; i < N; i++) o[i] = 0;
        int p = 0; for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++) { if (line[p++] == '1') o[i] |= 1u << j; else o[j] |= 1u << i; }
        double grp; uint64_t key = canon_key(o, N, &grp);
        if (hfind(key) >= 0) { fprintf(stderr, "duplicate input class\n"); exit(2); }
        hput(key, nclasses); memcpy(adj[nclasses], o, sizeof(uint32_t[16])); keys[nclasses] = key; grps[nclasses] = grp; nclasses++;
    }
    for (int c = 0; c < nclasses; c++) {
        memcpy(cur, adj[c], sizeof cur);
        uint32_t conv[16]; for (int i = 0; i < N; i++) { conv[i] = 0; for (int j = 0; j < N; j++) if (i != j && !(cur[i] >> j & 1)) conv[i] |= 1u << j; }
        int sc = (canon_key(conv, N, NULL) == keys[c]);
        nedges = 0; Hcount = 0;
        for (int v = 0; v < N; v++) { path[0] = v; visit(1, 1u << v); }
        printf("V %d %llx", c, (unsigned long long)keys[c]);
        for (int i = 0; i < N; i++) printf("%c%d", i ? ',' : ' ', __builtin_popcount(cur[i]));
        printf(" %d %d %.0f %ld\n", scc_count(cur, N) == 1, sc, grps[c], Hcount);
        for (int i = 0; i < nedges; i++) printf("E %d %d %d %d\n", c, edges[i].dst, edges[i].blue, edges[i].count);
    }
    fprintf(stderr, "n=%d classes=%d\n", N, nclasses);
    return 0;
}
