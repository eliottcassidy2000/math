/* gen_classes.c (AUDIT, independent) - isomorphism classes of tournaments of orders 1..maxn, generated from scratch
 * by one-vertex extension and a canonical form (colour refinement, then the lexicographically least arc sequence
 * over all orderings that respect the colour classes, with prefix pruning).
 * A wrong canonical form can only INCREASE the number of classes (two non-isomorphic tournaments never share a
 * code, because a code is a relabelling of the tournament itself), so the known counts 1, 1, 2, 4, 12, 56, 456,
 * 6880, 191536 (OEIS A000568) certify the lists.
 * Output files <prefix>_<n>.txt:  code H Pi2   (code: bit b, LSB first, for the pair (i<j) in lexicographic
 * order, 1 iff i -> j; H = Hamiltonian paths; Pi2 = ordered pairs of disjoint possibly empty paths covering T).
 * usage: gen_classes maxn prefix */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
#define MAXN 10
static int n, adj[MAXN][MAXN];
static int col[MAXN], order_[MAXN], cls_start[MAXN + 1], ncls;
static int best[MAXN], cur[MAXN], used[MAXN], havebest;
static unsigned char bestseq[64];
static void refine(void) {
    int sig[MAXN][MAXN + 1];
    for (int v = 0; v < n; v++) { int s = 0; for (int u = 0; u < n; u++) s += adj[v][u]; col[v] = s; }
    /* normalise colours to ranks */
    int nc = 0;
    for (int round = 0; round < 2 * n + 2; round++) {
        /* signature: own colour, then sorted colours of out-neighbours */
        for (int v = 0; v < n; v++) { int k = 0; sig[v][k++] = col[v]; int tmp[MAXN], t = 0;
            for (int u = 0; u < n; u++) if (adj[v][u]) tmp[t++] = col[u];
            for (int a = 1; a < t; a++) { int x = tmp[a], b = a - 1; while (b >= 0 && tmp[b] > x) { tmp[b + 1] = tmp[b]; b--; } tmp[b + 1] = x; }
            sig[v][k++] = t; for (int a = 0; a < t; a++) sig[v][k++] = tmp[a]; for (; k <= MAXN; k++) sig[v][k] = -1; }
        /* rank the signatures */
        int idx[MAXN]; for (int v = 0; v < n; v++) idx[v] = v;
        for (int a = 1; a < n; a++) { int x = idx[a], b = a - 1; while (b >= 0 && memcmp(sig[idx[b]], sig[x], sizeof sig[0]) > 0) { idx[b + 1] = idx[b]; b--; } idx[b + 1] = x; }
        int newcol[MAXN], c = 0; newcol[idx[0]] = 0;
        for (int a = 1; a < n; a++) { if (memcmp(sig[idx[a]], sig[idx[a - 1]], sizeof sig[0]) != 0) c++; newcol[idx[a]] = c; }
        for (int v = 0; v < n; v++) col[v] = newcol[v];
        if (c + 1 == nc) break; nc = c + 1;
    }
    /* positions: classes in colour order */
    int pos = 0; ncls = 0;
    for (int c = 0; c < n; c++) { int any = 0; for (int v = 0; v < n; v++) if (col[v] == c) { if (!any) { cls_start[ncls++] = pos; any = 1; } order_[pos++] = v; } }
    cls_start[ncls] = n;
}
static int poscls[MAXN];
/* sequence order: for j = 1..n-1, for i < j: bit [cur[i] -> cur[j]]; position of (i, j) = j (j - 1) / 2 + i */
/* compare the prefix cur[0..k] (all columns j <= k) with the best sequence found so far: -1, 0, +1.
 * Recomputed from scratch at every node (no flags are carried down, so an update of 'best' cannot leave stale state). */
static int cmp_prefix(int k) {
    for (int j = 1; j <= k; j++) { int base = j * (j - 1) / 2;
        for (int i = 0; i < j; i++) { int b = adj[cur[i]][cur[j]], bb = bestseq[base + i]; if (b < bb) return -1; if (b > bb) return 1; } }
    return 0; }
static void search(int k, int smaller) {
    (void)smaller;
    if (k == n) { if (!havebest || cmp_prefix(n - 1) < 0) { memcpy(best, cur, sizeof(int) * n); havebest = 1;
            for (int j = 1; j < n; j++) for (int i = 0; i < j; i++) bestseq[j * (j - 1) / 2 + i] = (unsigned char)adj[cur[i]][cur[j]]; }
        return; }
    int c = poscls[k];
    for (int a = cls_start[c]; a < cls_start[c + 1]; a++) { int v = order_[a]; if (used[v]) continue;
        cur[k] = v;
        if (havebest && cmp_prefix(k) > 0) continue;
        used[v] = 1; search(k + 1, 0); used[v] = 0; }
}
static u64 canon(void) {
    refine();
    for (int c = 0; c < ncls; c++) for (int p = cls_start[c]; p < cls_start[c + 1]; p++) poscls[p] = c;
    havebest = 0; memset(used, 0, sizeof used); search(0, 0);
    u64 code = 0; int b = 0;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (adj[best[i]][best[j]]) code |= (u64)1 << b; b++; }
    return code;
}
static void decode(u64 code, int m) { memset(adj, 0, sizeof adj); int b = 0;
    for (int i = 0; i < m; i++) for (int j = i + 1; j < m; j++) { if (code >> b & 1) adj[i][j] = 1; else adj[j][i] = 1; b++; } }
#define TSIZE (1u << 22)
static u64 *table; static unsigned nkeys;
static int insert(u64 key) { u64 h = (key + 1) * 0x9E3779B97F4A7C15ULL; unsigned i = (unsigned)(h >> 42) & (TSIZE - 1);
    while (table[i] != ~0ULL) { if (table[i] == key) return 0; i = (i + 1) & (TSIZE - 1); } table[i] = key; nkeys++; return 1; }
static int cmpu64(const void *a, const void *b) { u64 x = *(const u64 *)a, y = *(const u64 *)b; return x < y ? -1 : x > y; }
static void hp(int m, u64 *H, u64 *Pi2) {
    static u64 f[1 << MAXN][MAXN][2]; int full = (1 << m) - 1; memset(f, 0, sizeof(u64) * (size_t)(1 << m) * MAXN * 2);
    for (int v = 0; v < m; v++) f[1 << v][v][0] = 1;
    for (int mask = 1; mask <= full; mask++) for (int v = 0; v < m; v++) for (int b = 0; b < 2; b++) { u64 c = f[mask][v][b]; if (!c) continue;
        for (int u = 0; u < m; u++) { if (mask >> u & 1) continue;
            if (adj[v][u]) f[mask | 1 << u][u][b] += c;           /* continue the path */
            if (b == 0) f[mask | 1 << u][u][1] += c; } }           /* start the second path at u */
    u64 h = 0, p2 = 0; for (int v = 0; v < m; v++) { h += f[full][v][0]; p2 += f[full][v][1]; }
    *H = h; *Pi2 = p2 + 2 * h; }
int main(int argc, char **argv) {
    int maxn = atoi(argv[1]); const char *prefix = argv[2];
    table = (u64 *)malloc(sizeof(u64) * TSIZE);
    u64 *prev = (u64 *)malloc(sizeof(u64)); prev[0] = 0; unsigned nprev = 1;
    for (int m = 1; m <= maxn; m++) {
        u64 *codes; unsigned nc;
        if (m == 1) { codes = prev; nc = 1; }
        else {
            for (unsigned i = 0; i < TSIZE; i++) table[i] = ~0ULL; nkeys = 0;
            for (unsigned c = 0; c < nprev; c++) for (unsigned p = 0; p < (1u << (m - 1)); p++) {
                decode(prev[c], m - 1); n = m;
                for (int i = 0; i < m - 1; i++) { if (p >> i & 1) { adj[i][m - 1] = 1; adj[m - 1][i] = 0; } else { adj[m - 1][i] = 1; adj[i][m - 1] = 0; } }
                insert(canon()); }
            codes = (u64 *)malloc(sizeof(u64) * nkeys); nc = 0;
            for (unsigned i = 0; i < TSIZE; i++) if (table[i] != ~0ULL) codes[nc++] = table[i];
            qsort(codes, nc, sizeof(u64), cmpu64);
        }
        char fn[512]; sprintf(fn, "%s_%d.txt", prefix, m); FILE *fp = fopen(fn, "w");
        for (unsigned i = 0; i < nc; i++) { decode(codes[i], m); u64 H, P2; hp(m, &H, &P2); fprintf(fp, "%llu %llu %llu\n", codes[i], H, P2); }
        fclose(fp); printf("order %d: %u classes\n", m, nc); fflush(stdout);
        prev = codes; nprev = nc;
    }
    return 0;
}
