/* tcgen.c (klein-2026-11-04) -- a provably complete one-vertex extension generator for the tournaments with exactly one or
   exactly two Hamiltonian cycles.
   Theory.  Moon (1966): every strong tournament of order n >= 4 has a vertex v with T - v strong.  Extension lemma (opus-S74
   section 2 = klein's block lemma): if T - v is strong then every Hamiltonian cycle C' of T - v extends through v (v has an
   in -> out switch along C'), so hc(T) >= hc(T - v), and the out-neighbourhood O of v has at most hc(T) in -> out switches
   along C', i.e. O is a union of at most hc(T) cyclic intervals of C'.  Hence for n >= 5 every tournament with hc in {1, 2}
   is an extension of a tournament of order n - 1 with hc in {1, 2} by a vertex whose out-neighbourhood is a union of at most
   two proper cyclic intervals of one Hamiltonian cycle of the parent, and hc(child) = sum_{w in O, u in I} P(w -> u) with
   P the Hamiltonian-path counts of the parent.  Conversely every such extension with hc in {1, 2} is in the class.
   stdin: parent digraph6 strings (one per line, order m, each with hc in {1, 2}).  stdout: "hc digraph6" for every
   extension with hc in {1, 2} (not canonical; canonicalise with labelg and dedupe). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
static int m; static u64 adj[32]; static u64 *dp; static u64 P[32][32];
static int parse_d6(const char *s, u64 *a) {
    if (s[0] != '&') return -1; int n = s[1] - 63; int bits = n * n; const char *p = s + 2;
    for (int i = 0; i < n; i++) a[i] = 0;
    for (int k = 0; k < bits; k++) { int ch = p[k / 6] - 63; int b = (ch >> (5 - k % 6)) & 1; if (b) a[k / n] |= 1ULL << (k % n); }
    return n;
}
static void emit_d6(int n, u64 *a, char *out) {
    out[0] = '&'; out[1] = n + 63; int bits = n * n; int nch = (bits + 5) / 6;
    for (int c = 0; c < nch; c++) { int v = 0; for (int b = 0; b < 6; b++) { int k = c * 6 + b; int bit = 0; if (k < bits) bit = (a[k / n] >> (k % n)) & 1; v = (v << 1) | bit; } out[2 + c] = v + 63; }
    out[2 + nch] = 0;
}
/* paths from start s: dp[S][v] = # paths starting at s, visiting exactly S (contains s), ending at v */
static void paths_from(int s) {
    u64 N = 1ULL << m; memset(dp, 0, sizeof(u64) * N * m); dp[(1ULL << s) * m + s] = 1;
    for (u64 S = 1; S < N; S++) { if (!(S >> s & 1)) continue; for (int v = 0; v < m; v++) { u64 c = dp[S * m + v]; if (!c) continue;
        u64 nx = adj[v] & ~S; while (nx) { int u = __builtin_ctzll(nx); nx &= nx - 1; dp[(S | (1ULL << u)) * m + u] += c; } } }
}
int main(int argc, char **argv) {
    m = atoi(argv[1]); int K = argc > 2 ? atoi(argv[2]) : 2;   /* keep hc <= K; O with <= K intervals */ u64 N = 1ULL << m; dp = malloc(sizeof(u64) * N * m);
    char line[512], out[512]; u64 parents = 0, emitted = 0, skipped = 0;
    int cyc[32];
    while (fgets(line, sizeof line, stdin)) {
        int L = strlen(line); while (L && (line[L - 1] == '\n' || line[L - 1] == '\r')) line[--L] = 0; if (!L) continue;
        if (parse_d6(line, adj) != m) { skipped++; continue; }
        parents++;
        /* all-pairs Hamiltonian path counts P[w][u] and one Hamiltonian cycle (from the DP with start 0) */
        for (int w = 0; w < m; w++) { paths_from(w); for (int u = 0; u < m; u++) P[w][u] = dp[(N - 1) * m + u]; }
        paths_from(0); u64 hc = 0; int endv = -1; for (int v = 1; v < m; v++) if (adj[v] & 1ULL) { hc += dp[(N - 1) * m + v]; if (dp[(N - 1) * m + v]) endv = v; }
        if (hc < 1 || hc > (u64)K) { fprintf(stderr, "parent with hc = %llu skipped: %s\n", hc, line); continue; }
        /* reconstruct one cycle 0 = cyc[0], cyc[1], ..., cyc[m-1] = endv (path from 0 ending at endv, endv -> 0) */
        u64 S = N - 1; int cur = endv; cyc[m - 1] = cur;
        for (int pos = m - 1; pos > 0; pos--) { u64 Sp = S & ~(1ULL << cur); int found = -1;
            for (int i = 0; i < m; i++) if ((Sp >> i & 1) && dp[Sp * m + i] && (adj[i] >> cur & 1)) { found = i; break; }
            cur = found; S = Sp; cyc[pos - 1] = cur; }
        /* patterns: every proper subset O of the cycle with exactly 2 or 4 in/out switches, enumerated once each by its
           sorted switch positions (p < q: O = [p,q) or its complement; p < q < r < s: O = [p,q) u [r,s) or its complement) */
        int pos[4];
        for (int k = 2; k <= 2 * K; k += 2) { if (k > m) break;
            int idx[4]; for (int i = 0; i < k; i++) idx[i] = i;
            while (1) {
                for (int phase = 0; phase < 2; phase++) {
                    u64 O = 0;
                    /* membership of cyc[i] in O: count switches at positions <= i; phase chooses the state at position 0 */
                    for (int i = 0, st = phase, j = 0; i < m; i++) { while (j < k && idx[j] == i) { st ^= 1; j++; } if (st) O |= 1ULL << cyc[i]; }
                    u64 I = (N - 1) & ~O; if (!O || !I) continue; u64 h = 0;
                    for (u64 x = O; x; x &= x - 1) { int w = __builtin_ctzll(x); for (u64 y = I; y; y &= y - 1) { int u = __builtin_ctzll(y); h += P[w][u]; if (h > (u64)K) break; } if (h > (u64)K) break; }
                    if (h >= 1 && h <= (u64)K) { u64 ch[32]; for (int i = 0; i < m; i++) ch[i] = adj[i] | ((I >> i & 1) ? (1ULL << m) : 0); ch[m] = O; emit_d6(m + 1, ch, out); printf("%llu %s %llu\n", h, out, parents); emitted++; }
                }
                /* next k-subset of {0..m-1} in lexicographic order */
                int i = k - 1; while (i >= 0 && idx[i] == m - k + i) i--; if (i < 0) break; idx[i]++; for (int j = i + 1; j < k; j++) idx[j] = idx[j - 1] + 1;
            }
        }
        (void)pos;
    }
    fprintf(stderr, "order %d parents %llu (skipped %llu), extensions with hc in {1,2}: %llu (labelled by parent and pattern, before canonical dedupe)\n", m, parents, skipped, emitted);
    return 0;
}
