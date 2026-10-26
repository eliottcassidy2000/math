// legality_strong_lemma_20261026.cpp -- opus-2026-10-26-S61
// Strong form of the single-pair lemma (THM-4806 R (i)):
//   for every tournament T and every dark arc e = u -> v of T (an arc on no Hamiltonian path), the tournament T^(uv)
//   has a Hamiltonian path through v -> u that uses no other dark arc of T.
// Why it matters: if it held for all STRONG tournaments it would imply all-order injectivity of T -> L(T): given
// L(T) = L(T'), the two tournaments have the same strong components (THM-4806 R (ii)); inside a component C the
// differing pairs D are dark in C; for e in D the avoiding path of C^(e) uses e reversed and otherwise arcs of C
// outside D, hence is a Hamiltonian path of C' through a D-pair, which extends to one of T' - a contradiction.
// This program tests the strong form exhaustively (all labelled tournaments, n <= maxA; all of X_8) and on random
// samples of X_9, X_10, reporting violations overall and among strong tournaments.  Arc masks are 128-bit (index i*n+j).
// Build: g++ -O2 -std=c++17 -o lsl legality_strong_lemma_20261026.cpp ; run: ./lsl [maxA=7] [samples9=30000] [samples10=3000]
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <vector>
#include <algorithm>
using namespace std;
typedef unsigned long long u64; typedef unsigned __int128 u128;

static int n; static int adj[11][11];
static u128 arcsUsed; static int order[11];
static u128 darkMask; static int flipU, flipV; static bool found;
static inline u128 bit(int i, int j) { return (u128)1 << (i * n + j); }

static void dfsAll(int pos, int used) {
    if (pos == n) { for (int i = 0; i + 1 < n; i++) arcsUsed |= bit(order[i], order[i + 1]); return; }
    int last = order[pos - 1];
    for (int v = 0; v < n; v++) if (!(used >> v & 1) && adj[last][v]) { order[pos] = v; dfsAll(pos + 1, used | 1 << v); }
}
static void usedArcs() { arcsUsed = 0; for (int v = 0; v < n; v++) { order[0] = v; dfsAll(1, 1 << v); } }
static void dfsAvoid(int pos, int used, bool hasFlip) {
    if (found) return;
    if (pos == n) { if (hasFlip) found = true; return; }
    int last = order[pos - 1];
    for (int v = 0; v < n; v++) {
        if ((used >> v & 1) || !adj[last][v]) continue;
        bool isFlip = (last == flipV && v == flipU);
        if (!isFlip && (darkMask & bit(last, v))) continue;
        order[pos] = v; dfsAvoid(pos + 1, used | 1 << v, hasFlip || isFlip);
    }
}
static bool avoidingPathExists() { found = false; for (int v = 0; v < n && !found; v++) { order[0] = v; dfsAvoid(1, 1 << v, false); } return found; }
static bool strongT() {
    int r = 1, prev = 0; while (r != prev) { prev = r; for (int i = 0; i < n; i++) if (r >> i & 1) for (int j = 0; j < n; j++) if (adj[i][j]) r |= 1 << j; }
    if (r != (1 << n) - 1) return false;
    r = 1; prev = 0; while (r != prev) { prev = r; for (int i = 0; i < n; i++) if (r >> i & 1) for (int j = 0; j < n; j++) if (adj[j][i]) r |= 1 << j; }
    return r == (1 << n) - 1;
}
struct Stat { u64 tourn = 0, withDark = 0, darkArcs = 0, violations = 0, strongCnt = 0, strongDark = 0, strongViol = 0; };
static void testTournament(Stat& st, u64 code) {
    usedArcs(); u128 used = arcsUsed; u128 dark = 0; bool stg = strongT();
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (adj[i][j] && !(used & bit(i, j))) dark |= bit(i, j);
    st.tourn++; if (stg) st.strongCnt++; if (!dark) return; st.withDark++;
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (dark & bit(i, j)) {
        st.darkArcs++; if (stg) st.strongDark++;
        flipU = i; flipV = j; darkMask = dark;
        adj[i][j] = 0; adj[j][i] = 1;
        bool ok = avoidingPathExists();
        adj[i][j] = 1; adj[j][i] = 0;
        if (!ok) { st.violations++; if (stg) { st.strongViol++; if (st.strongViol <= 4) printf("  strong violation: n = %d code %llu dark arc %d->%d\n", n, code, i, j); } }
    }
}
static u64 rng = 0x9E3779B97F4A7C15ULL;
static u64 next64() { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return rng * 0x2545F4914F6CDD1DULL; }
static void report(const char* label, const Stat& st) {
    printf("%s: tournaments %llu, with dark arcs %llu, dark arcs %llu, violations %llu; strong %llu with %llu dark arcs, violations among strong %llu\n",
           label, st.tourn, st.withDark, st.darkArcs, st.violations, st.strongCnt, st.strongDark, st.strongViol);
}

int main(int argc, char** argv) {
    int maxA = argc > 1 ? atoi(argv[1]) : 7; int s9 = argc > 2 ? atoi(argv[2]) : 30000; int s10 = argc > 3 ? atoi(argv[3]) : 3000;
    vector<pair<int,int>> pairs; char label[64];
    for (n = 3; n <= maxA; n++) {
        pairs.clear(); for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) pairs.push_back({i, j});
        int M = (int)pairs.size(); Stat st;
        for (u64 code = 0; code < (1ULL << M); code++) {
            for (int k = 0; k < M; k++) { int i = pairs[k].first, j = pairs[k].second; adj[i][j] = (code >> k) & 1; adj[j][i] = 1 - adj[i][j]; }
            testTournament(st, code);
        }
        snprintf(label, sizeof label, "n = %d (all labelled)", n); report(label, st);
    }
    for (auto [nn, samples] : vector<pair<int,int>>{{8, 0}, {9, s9}, {10, s10}}) {
        n = nn; pairs.clear(); for (int i = 0; i < n; i++) for (int j = i + 2; j < n; j++) pairs.push_back({i, j}); int M = (int)pairs.size(); Stat st;
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) adj[i][j] = 0;
        for (int i = 0; i + 1 < n; i++) { adj[i][i + 1] = 1; adj[i + 1][i] = 0; }
        if (samples == 0) {
            for (u64 code = 0; code < (1ULL << M); code++) {
                for (int k = 0; k < M; k++) { int i = pairs[k].first, j = pairs[k].second; adj[i][j] = (code >> k) & 1; adj[j][i] = 1 - adj[i][j]; }
                testTournament(st, code);
            }
            snprintf(label, sizeof label, "X_%d (exhaustive)", n);
        } else {
            for (int t = 0; t < samples; t++) {
                u64 code = next64() & ((1ULL << M) - 1);
                for (int k = 0; k < M; k++) { int i = pairs[k].first, j = pairs[k].second; adj[i][j] = (code >> k) & 1; adj[j][i] = 1 - adj[i][j]; }
                testTournament(st, code);
            }
            snprintf(label, sizeof label, "X_%d (random %d)", n, samples);
        }
        report(label, st);
    }
    return 0;
}
