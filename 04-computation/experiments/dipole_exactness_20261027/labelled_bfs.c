/* labelled_bfs.c -- exact labelled native-move distances at order n <= 8 (2^C(n,2) states).
 * Start: the Paley tournament P_7 on Z_7 (i -> i+1, i+2, i+4) for n = 7, or the regular / any tournament given as an
 * arc string.  Moves: tau_Q(T) = T^op xor E(Q) over all Hamiltonian paths Q.  Reports
 *   (a) the distance to the nearest transitive tournament (any labelling: the transitive class),
 *   (b) the distance to the fixed transitive tournament TT_id (i -> j iff i < j),
 *   (c) an explicit geodesic word for (a), every step re-verified,
 *   (d) the parity of the number of forward arcs (THM-4810 (C)) along the way.
 * Build: cc -O2 labelled_bfs.c -o labelled_bfs ;  Run: ./labelled_bfs 7 [arcstring]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
static int n, NP;
static int pi_[64], pj_[64], pidx[16][16];
static void decode(uint32_t code, uint32_t *out) { for (int i = 0; i < n; i++) out[i] = 0; for (int k = 0; k < NP; k++) { int i = pi_[k], j = pj_[k]; if (code >> k & 1) out[i] |= 1u << j; else out[j] |= 1u << i; } }
static uint32_t encode(const uint32_t *out) { uint32_t c = 0; for (int k = 0; k < NP; k++) if (out[pi_[k]] >> pj_[k] & 1) c |= 1u << k; return c; }
static int is_transitive(const uint32_t *out) { int seen = 0; for (int v = 0; v < n; v++) { int s = __builtin_popcount(out[v]); if (seen >> s & 1) return 0; seen |= 1 << s; } return 1; }
static uint32_t TTid;
static int Qs[50000][16]; static int nQ;
static void collect_paths(const uint32_t *out) {
    nQ = 0;
    for (int s = 0; s < n; s++) {
        int Q[16]; uint32_t cand[16]; int depth = 0; Q[0] = s; uint32_t used = 1u << s; cand[0] = out[s] & ~used;
        while (depth >= 0) {
            if (depth == n - 1) { memcpy(Qs[nQ++], Q, sizeof Q); used &= ~(1u << Q[depth]); depth--; continue; }
            if (!cand[depth]) { used &= ~(1u << Q[depth]); depth--; continue; }
            int v = __builtin_ctz(cand[depth]); cand[depth] &= cand[depth] - 1;
            depth++; Q[depth] = v; used |= 1u << v; cand[depth] = out[v] & ~used;
        }
    }
}
static uint32_t move(const uint32_t *out, const int *Q) {
    uint32_t r[16]; uint32_t full = (1u << n) - 1;
    for (int u = 0; u < n; u++) r[u] = (~out[u]) & full & ~(1u << u);
    for (int i = 0; i + 1 < n; i++) { int u = Q[i], v = Q[i + 1]; r[u] |= 1u << v; r[v] &= ~(1u << u); }
    return encode(r);
}
int main(int argc, char **argv) {
    n = atoi(argv[1]); NP = n * (n - 1) / 2;
    int k = 0; for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { pi_[k] = i; pj_[k] = j; pidx[i][j] = k; k++; }
    uint32_t start[16];
    if (argc > 2) { const char *s = argv[2]; uint32_t c = 0; for (int q = 0; q < NP; q++) if (s[q] == '1') c |= 1u << q; decode(c, start); }
    else { for (int i = 0; i < n; i++) { start[i] = 0; for (int j = 0; j < n; j++) { int d = ((j - i) % n + n) % n; if (d == 1 || d == 2 || d == 4) start[i] |= 1u << j; } } }
    uint32_t tt[16]; for (int i = 0; i < n; i++) { tt[i] = 0; for (int j = i + 1; j < n; j++) tt[i] |= 1u << j; } TTid = encode(tt);
    uint32_t S = 1u << NP;
    int8_t *dist = malloc(S); uint32_t *parent = malloc(S * sizeof(uint32_t)); uint32_t *queue = malloc(S * sizeof(uint32_t));
    memset(dist, -1, S);
    uint32_t s0 = encode(start); dist[s0] = 0; parent[s0] = s0; long qh = 0, qt = 0; queue[qt++] = s0;
    int dTrans = -1, dTTid = -1; uint32_t hitTrans = 0;
    while (qh < qt) {
        uint32_t c = queue[qh++]; uint32_t out[16]; decode(c, out);
        if (dTrans < 0 && is_transitive(out)) { dTrans = dist[c]; hitTrans = c; }
        if (c == TTid) { dTTid = dist[c]; }
        if (dTrans >= 0 && dTTid >= 0) break;
        collect_paths(out);
        for (int q = 0; q < nQ; q++) { uint32_t t = move(out, Qs[q]); if (dist[t] < 0) { dist[t] = dist[c] + 1; parent[t] = c; queue[qt++] = t; } }
    }
    printf("n=%d start forward-arc parity %d; distance to the transitive class %d; distance to TT_id %d; states visited %ld\n",
           n, __builtin_popcount(s0) & 1, dTrans, dTTid, qt);
    /* geodesic to the transitive class, re-verified step by step */
    uint32_t chain[64]; int len = 0; for (uint32_t c = hitTrans; ; c = parent[c]) { chain[len++] = c; if (c == s0) break; }
    printf("geodesic (forward-arc count, scores):\n");
    for (int i = len - 1; i >= 0; i--) {
        uint32_t out[16]; decode(chain[i], out);
        int sc[16]; for (int v = 0; v < n; v++) sc[v] = __builtin_popcount(out[v]);
        printf("  step %d: F=%d scores", len - 1 - i, __builtin_popcount(chain[i]));
        for (int v = 0; v < n; v++) printf(" %d", sc[v]);
        if (i > 0) {   /* find and print the path realising the next step */
            collect_paths(out); int found = -1;
            for (int q = 0; q < nQ; q++) if (move(out, Qs[q]) == chain[i - 1]) { found = q; break; }
            if (found < 0) { printf("  VERIFY FAILED\n"); return 1; }
            printf("  path");
            for (int t = 0; t < n; t++) printf(" %d", Qs[found][t]);
        }
        printf("\n");
    }
    return 0;
}
