/* prescribed_isolation.c -- mac-mini-2026-10-27-slab
 * For a regular tournament T on n = 2m+1 vertices, a prescribed vertex a and a role (source or sink): is there a word
 * of exactly m BLACK native moves after which a is a source (resp. sink) and T - a is strong?
 * Dipole law (THM-4811 (A)): d(a) = s(a) - m goes to -d(a) + [a = start] - [a = end], so |d(a)| = m after m moves forces
 * a to be an endpoint of every path with alternating roles; the first role is fixed by the target role and m:
 *   first move with a as START: d = 1, -2, 3, ...  (source iff m odd);  as END: d = -1, 2, -3, ... (source iff m even).
 * Depth-first search over black Hamiltonian paths with a at the required end (random order, memo of failed states by a
 * canonical form with a marked), final check that T - a is strong.  Reads gentourng output (regular classes, e.g.
 * gentourng -d3D3 7); tests every vertex and both roles of every class.
 * Build: cc -O3 -I/opt/homebrew/opt/nauty/include prescribed_isolation.c /opt/homebrew/opt/nauty/lib/libnauty.a -o prescribed_isolation
 * Run:   gentourng -q -d<m>D<m> n | ./prescribed_isolation n [node_budget]
 */
#define MAXN 32
#include <nauty/nauty.h>
#include <nauty/nautinv.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
static int n, m; static uint64_t rs = 88172645463325252ULL;
static uint64_t rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static int strongMask(const uint32_t *o, uint32_t vs) {
    if (!vs) return 1; int start = __builtin_ctz(vs); uint32_t seen = 1u << start, fr = seen;
    while (fr) { uint32_t nx = 0; for (uint32_t f = fr; f; f &= f - 1) nx |= o[__builtin_ctz(f)] & vs; nx &= ~seen; seen |= nx; fr = nx; }
    if (seen != vs) return 0; seen = 1u << start; fr = seen;
    while (fr) { uint32_t nx = 0; for (uint32_t f = fr; f; f &= f - 1) { int v = __builtin_ctz(f); for (uint32_t w = vs; w; w &= w - 1) { int u = __builtin_ctz(w); if (o[u] >> v & 1) nx |= 1u << u; } } nx &= ~seen; seen |= nx; fr = nx; }
    return seen == vs;
}
static int isBlack(const uint32_t *o, const int *Q) { int r[32]; for (int i = 0; i < n; i++) r[Q[i]] = Q[n - 1 - i];
    for (int u = 0; u < n; u++) for (uint32_t x = o[u]; x; x &= x - 1) { int v = __builtin_ctz(x); if (!((o[r[v]] >> r[u]) & 1)) return 1; } return 0; }
static void apply_move(const uint32_t *o, const int *Q, uint32_t *r) { uint32_t full = (1u << n) - 1;
    for (int u = 0; u < n; u++) r[u] = (~o[u]) & full & ~(1u << u);
    for (int i = 0; i + 1 < n; i++) { int u = Q[i], v = Q[i + 1]; r[u] |= 1u << v; r[v] &= ~(1u << u); } }
typedef struct { uint64_t w[2]; } Key;
static Key canon_marked(const uint32_t *q, int a) {   /* canonical form with vertex a in its own colour class */
    graph g[MAXN], cg[MAXN]; int lab[MAXN], ptn[MAXN], orbits[MAXN];
    DEFAULTOPTIONS_DIGRAPH(options); statsblk stats; options.getcanon = TRUE; options.defaultptn = FALSE;
    EMPTYGRAPH(g, 1, n);
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (q[i] >> j & 1) ADDONEARC(g, i, j, 1);
    int k = 0; lab[k] = a; ptn[k] = 0; k++;
    for (int v = 0; v < n; v++) if (v != a) { lab[k] = v; ptn[k] = 1; k++; }
    ptn[n - 1] = 0;
    densenauty(g, lab, ptn, orbits, &options, &stats, 1, n, cg);
    Key key = {{0, 0}}; int pos = 0;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (ISELEMENT(GRAPHROW(cg, i, 1), j)) key.w[pos >> 6] |= 1ULL << (pos & 63); pos++; }
    return key;
}
#define HC (1 << 20)
static Key hk[HC]; static int hd[HC]; static unsigned char hu[HC];
static int truncated;
static int seen_fail(Key k, int depth) { uint64_t h = ((k.w[0] * 0x9E3779B97F4A7C15ULL) ^ (k.w[1] * 0xC2B2AE3D27D4EB4FULL) ^ (uint64_t)depth * 0x165667B19E3779F9ULL) >> 44;
    while (hu[h]) { if (hk[h].w[0] == k.w[0] && hk[h].w[1] == k.w[1] && hd[h] == depth) return 1; h = (h + 1) & (HC - 1); } return 0; }
static void add_fail(Key k, int depth) { uint64_t h = ((k.w[0] * 0x9E3779B97F4A7C15ULL) ^ (k.w[1] * 0xC2B2AE3D27D4EB4FULL) ^ (uint64_t)depth * 0x165667B19E3779F9ULL) >> 44;
    while (hu[h]) { if (hk[h].w[0] == k.w[0] && hk[h].w[1] == k.w[1] && hd[h] == depth) return; h = (h + 1) & (HC - 1); } hu[h] = 1; hk[h] = k; hd[h] = depth; }
static long nodes, budget;
/* enumerate Hamiltonian paths with a at the start (asStart) or at the end; collect up to cap of them, black only */
static int collect(const uint32_t *o, int a, int asStart, int (*Qs)[32], int cap) {
    int cnt = 0; uint32_t full = (1u << n) - 1;
    /* DFS over paths; for 'end' role, enumerate paths in the converse from a and reverse them */
    uint32_t oo[32]; if (asStart) memcpy(oo, o, n * sizeof(uint32_t)); else for (int u = 0; u < n; u++) oo[u] = (~o[u]) & full & ~(1u << u);
    int Q[32]; uint32_t cand[32]; int depth = 0; Q[0] = a; uint32_t used = 1u << a; cand[0] = oo[a] & ~used;
    while (depth >= 0) {
        if (depth == n - 1) {
            int P[32]; if (asStart) memcpy(P, Q, n * sizeof(int)); else for (int i = 0; i < n; i++) P[i] = Q[n - 1 - i];
            if (isBlack(o, P)) {   /* reservoir sampling to keep a uniform random subset of size cap */
                if (cnt < cap) memcpy(Qs[cnt], P, sizeof P); else { long r = rnd() % (cnt + 1); if (r < cap) memcpy(Qs[r], P, sizeof P); }
                cnt++;
            }
            used &= ~(1u << Q[depth]); depth--; continue;
        }
        if (!cand[depth]) { used &= ~(1u << Q[depth]); depth--; continue; }
        int v = __builtin_ctz(cand[depth]); cand[depth] &= cand[depth] - 1;
        depth++; Q[depth] = v; used |= 1u << v; cand[depth] = oo[v] & ~used;
    }
    if (cnt > cap) truncated = 1;
    return cnt < cap ? cnt : cap;
}
static int search(const uint32_t *o, int a, int k, int asStart) {   /* k moves left; asStart = role of a in the next move */
    uint32_t full = (1u << n) - 1;
    if (k == 0) { int s = __builtin_popcount(o[a]); return (s == 0 || s == n - 1) && strongMask(o, full & ~(1u << a)); }
    if (++nodes > budget) return -1;
    Key key = canon_marked(o, a); if (seen_fail(key, k * 2 + asStart)) return 0;
    static int Qs[8][16384][32]; int depthIdx = k;   /* per-depth buffers (k <= 7) */
    int nq = collect(o, a, asStart, Qs[depthIdx], 16384);
    for (int i = 0; i < nq; i++) {
        uint32_t r[32]; apply_move(o, Qs[depthIdx][i], r);
        int res = search(r, a, k - 1, !asStart);
        if (res != 0) return res;
    }
    add_fail(key, k * 2 + asStart); return 0;
}
int main(int argc, char **argv) {
    n = atoi(argv[1]); m = (n - 1) / 2; budget = argc > 2 ? atol(argv[2]) : 200000;
    char line[256]; long classes = 0, inst = 0, ok = 0, undec = 0, fail = 0;
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] != '0' && line[0] != '1') continue;
        classes++;
        uint32_t o[32] = {0}; int k = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (line[k] == '1') o[i] |= 1u << j; else o[j] |= 1u << i; k++; }
        for (int a = 0; a < n; a++) for (int role = 0; role < 2; role++) {   /* role 0: source, 1: sink */
            int firstStart = ((m % 2 == 1) == (role == 0));   /* first move with a as start gives source iff m odd */
            inst++; nodes = 0; truncated = 0; memset(hu, 0, sizeof hu);
            int res = search(o, a, m, firstStart);
            if (res == 1) ok++; else if (res < 0) { undec++; printf("UNDECIDED class %ld vertex %d role %s\n", classes, a, role ? "sink" : "source"); }
            else if (truncated) { undec++; printf("UNDECIDED (path sample truncated) class %ld vertex %d role %s\n", classes, a, role ? "sink" : "source"); }
            else { fail++; printf("NO WORD class %ld vertex %d role %s: %s", classes, a, role ? "sink" : "source", line); }
        }
    }
    printf("n=%d regular classes %ld; (vertex, role) instances %ld; isolated in exactly m = %d black moves with strong core: %ld; none exists: %ld; undecided: %ld\n", n, classes, inst, m, ok, fail, undec);
    return 0;
}
