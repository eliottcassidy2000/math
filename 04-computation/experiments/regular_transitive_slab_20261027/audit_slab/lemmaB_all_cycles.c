/* lemmaB_all_cycles.c -- independent audit of THM-4819 Lemma B (odd n = 2m+1, T regular) and
 * Lemma B' (even n = 2m, T balanced: scores m-1, m).  Written from the definitions only; it shares
 * no code with the producer's scripts.
 *
 * Input: gentourng upper-triangle lines on stdin (char k for the pair (i,j), i<j, in row order:
 * '1' means i -> j).  For every tournament T and EVERY vertex v it enumerates ALL directed Hamilton
 * cycles c_0 -> ... -> c_{M-1} -> c_0 of T - v (M = n-1; c_0 = the least vertex of T - v, so each
 * directed cycle is met once), and for every rotation P = (c_t, c_{t+1}, ..., c_{t-1}, v) with
 * c_{t-1} -> v it tests whether P is BLUE (positional reversal p_i <-> p_{n-1-i} is an
 * anti-automorphism: x -> y iff sig(y) -> sig(x)).  A cycle all of whose rotations are blue is a
 * counterexample to the strong form of the lemma ("EVERY Hamilton cycle of T - v has a black
 * rotation closed at v").  For every blue rotation the identity
 *     (*)  for x != t:  c_t -> c_x  iff  (2t - x + 1 mod M) in U,   U = {t : c_{t-1} -> v}
 * is checked as well (it must hold whenever sig_t is an anti-automorphism).
 * Vertices are grouped by score: odd n -- all scores are m; even n -- score m (the lemma's case) and
 * score m-1 (not claimed; reported for information).
 *
 * Build: cc -O2 -o lemmaB_all_cycles lemmaB_all_cycles.c      Run: gentourng -q -d4D4 9 | ./lemmaB_all_cycles 9
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int n, M;
static uint32_t out_[32];
static inline int arc(int a, int b) { return (int)((out_[a] >> b) & 1u); }

static int is_anti(const int *sig) {
    for (int x = 0; x < n; x++)
        for (int y = x + 1; y < n; y++)
            if (arc(x, y) != arc(sig[y], sig[x])) return 0;
    return 1;
}

/* statistics, index g = 0 for "score m" (odd: every vertex), g = 1 for "score m-1" (even only) */
static long long ncycles[2], nrot[2], nblue[2], allblue[2], starfail, hist[2][40], noCycleVertex[2];
static long long bluestart_scorem[2];
static int vtx, cyc[32], used[32];
static long long cyc_this_vertex;
static int grp;
static int scorev[32];
static int printed = 0;

static void process_cycle(void) {
    int inU[32];
    for (int t = 0; t < M; t++) inU[t] = arc(cyc[(t - 1 + M) % M], vtx);
    int rots = 0, blue = 0;
    for (int t = 0; t < M; t++) {
        if (!inU[t]) continue;
        rots++;
        int P[32], sig[32];
        for (int i = 0; i < M; i++) P[i] = cyc[(t + i) % M];
        P[M] = vtx;
        for (int i = 0; i < n; i++) sig[P[i]] = P[n - 1 - i];
        for (int i = 0; i + 1 < n; i++) if (!arc(P[i], P[i + 1])) { fprintf(stderr, "not a path\n"); exit(3); }
        if (is_anti(sig)) {
            blue++;
            if (scorev[P[0]] + scorev[P[n - 1]] != n - 1) bluestart_scorem[grp]++;  /* impossible for blue */
            for (int x = 0; x < M; x++) {
                if (x == t) continue;
                int w = ((2 * t - x + 1) % M + 2 * M) % M;
                if (arc(cyc[t], cyc[x]) != inU[w]) starfail++;
            }
        }
    }
    ncycles[grp]++; nrot[grp] += rots; nblue[grp] += blue; hist[grp][blue]++;
    cyc_this_vertex++;
    if (blue == rots) {
        allblue[grp]++;
        if (printed < 12) {
            printed++;
            printf("  ALL-BLUE cycle: v=%d (score %d), cycle", vtx, scorev[vtx]);
            for (int i = 0; i < M; i++) printf(" %d", cyc[i]);
            printf(" (%d rotations)\n", rots);
        }
    }
}

static void dfs(int depth) {
    if (depth == M) { if (arc(cyc[M - 1], cyc[0])) process_cycle(); return; }
    int last = cyc[depth - 1];
    for (int w = 0; w < n; w++) {
        if (used[w] || !arc(last, w)) continue;
        used[w] = 1; cyc[depth] = w; dfs(depth + 1); used[w] = 0;
    }
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s n\n", argv[0]); return 1; }
    n = atoi(argv[1]); M = n - 1;
    int m = n / 2;            /* odd n = 2m+1: m = (n-1)/2; even n = 2m: m = n/2 */
    int odd = n & 1;
    char line[1024]; long ntour = 0;
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] != '0' && line[0] != '1') continue;
        memset(out_, 0, sizeof out_);
        int k = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
            if (line[k] == '1') out_[i] |= 1u << j; else out_[j] |= 1u << i;
            k++;
        }
        for (int v = 0; v < n; v++) scorev[v] = __builtin_popcount(out_[v]);
        for (int v = 0; v < n; v++) {
            if (odd && scorev[v] != m) { fprintf(stderr, "not regular\n"); return 2; }
            if (!odd && scorev[v] != m && scorev[v] != m - 1) { fprintf(stderr, "not balanced\n"); return 2; }
        }
        ntour++;
        for (int v = 0; v < n; v++) {
            vtx = v; grp = (odd || scorev[v] == m) ? 0 : 1;
            memset(used, 0, sizeof used);
            used[v] = 1;
            int c0 = (v == 0) ? 1 : 0;
            used[c0] = 1; cyc[0] = c0; cyc_this_vertex = 0;
            dfs(1);
            if (cyc_this_vertex == 0) noCycleVertex[grp]++;
        }
    }
    printf("n=%d: tournaments %ld\n", n, ntour);
    for (int g = 0; g < (odd ? 1 : 2); g++) {
        printf("  %s: (T, v, Hamilton cycle of T-v) instances %lld; rotations %lld; blue rotations %lld; "
               "cycles with ALL rotations blue %lld; vertices v with T-v non-Hamiltonian %lld; blue rotations with endpoint sum != n-1 %lld\n",
               odd ? "every vertex (score m)" : (g == 0 ? "v of score m (Lemma B')" : "v of score m-1 (not claimed)"),
               ncycles[g], nrot[g], nblue[g], allblue[g], noCycleVertex[g], bluestart_scorem[g]);
        printf("    histogram of #blue rotations per cycle:");
        for (int b = 0; b < 40; b++) if (hist[g][b]) printf(" %d:%lld", b, hist[g][b]);
        printf("\n");
    }
    printf("  identity (*) failures over all blue rotations: %lld\n", starfail);
    return 0;
}
