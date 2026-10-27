// black_distance.cpp -- exact black distance to the non-strong locus by labelled BFS (black moves only).
// A black move from T is K_Q(T) for a BLACK Hamiltonian path Q of T: keep the arcs of Q, reverse every other pair.
// Usage: echo CODE | ./black_distance n [maxdepth]
// Prints, level by level, the number of distinct labelled tournaments reached, their lambda profile, and stops at
// the first level containing a non-strong tournament.  Also independently re-verifies the cut identity
// f_{K_Q T}(V\X) = f_T(X) - [q1 in X] + [qn in X] on every move from the root.
#include "tour.hpp"
#include <unordered_set>
#include <map>
static Tour kmove(const Tour &T, const int *Q) {
    Tour R; R.n = T.n; uint32_t full = (1u << T.n) - 1;
    for (int u = 0; u < T.n; u++) R.out[u] = (~T.out[u]) & full & ~(1u << u);
    for (int i = 0; i + 1 < T.n; i++) { int u = Q[i], v = Q[i + 1]; R.out[u] |= 1u << v; R.out[v] &= ~(1u << u); }
    R.finish(); return R;
}
template <class F> static void allPaths(const Tour &T, F &&f) {
    int n = T.n;
    for (int s0 = 0; s0 < n; s0++) {
        int Q[20]; uint32_t cand[20]; int depth = 0; Q[0] = s0; uint32_t used = 1u << s0; cand[0] = T.out[s0] & ~used;
        while (depth >= 0) {
            if (depth == n - 1) { f(Q); used &= ~(1u << Q[depth]); depth--; continue; }
            if (!cand[depth]) { used &= ~(1u << Q[depth]); depth--; continue; }
            int v = __builtin_ctz(cand[depth]); cand[depth] &= cand[depth] - 1;
            depth++; Q[depth] = v; used |= 1u << v; cand[depth] = T.out[v] & ~used;
        }
    }
}
static string key(const Tour &T) { string s((size_t)T.n * 4, 0); memcpy(&s[0], T.out, (size_t)T.n * 4); return s; }
int main(int argc, char **argv) {
    int n = argc > 1 ? atoi(argv[1]) : 10; int maxd = argc > 2 ? atoi(argv[2]) : 3;
    char line[512]; Tour T0;
    while (fgets(line, sizeof line, stdin)) {
        if (!parseLine(line, n, T0)) continue;
        printf("root %s  lambda=%d\n", encode(T0).c_str(), ({ vector<uint32_t> c; minCuts(T0, c); }));
        // cut identity check on all moves from the root
        long long idChecks = 0, idFails = 0; uint32_t full = (1u << n) - 1;
        allPaths(T0, [&](const int *Q) {
            Tour R = kmove(T0, Q);
            for (uint32_t X = 1; X < full; X++) {
                int lhs = fOut(R, full & ~X), rhs = fOut(T0, X) - (int)(X >> Q[0] & 1) + (int)(X >> Q[n - 1] & 1);
                idChecks++; if (lhs != rhs) idFails++;
            }
        });
        printf("  cut identity checked on %lld (path, X) pairs: %lld failures\n", idChecks, idFails);
        unordered_set<string> seen; vector<Tour> frontier{T0}; seen.insert(key(T0));
        bool done = false;
        for (int d = 1; d <= maxd && !done; d++) {
            vector<Tour> next; long long moves = 0;
            for (const Tour &T : frontier) {
                allPaths(T, [&](const int *Q) {
                    if (isBlue(T, Q)) return;
                    moves++;
                    Tour R = kmove(T, Q);
                    string k = key(R);
                    if (seen.insert(k).second) next.push_back(R);
                });
            }
            map<int, long long> lamHist; long long nonstrong = 0;
            for (const Tour &R : next) { vector<uint32_t> c; int l = minCuts(R, c); lamHist[l]++; if (l == 0) nonstrong++; }
            printf("  level %d: black moves examined %lld, new labelled tournaments %zu, lambda histogram:", d, moves, next.size());
            for (auto &p : lamHist) printf(" %d:%lld", p.first, p.second);
            printf("\n");
            if (nonstrong) { printf("  => black distance to the non-strong locus = %d\n", d); done = true; }
            frontier.swap(next);
            if (frontier.size() > 3000000) { printf("  frontier too large, stopping\n"); break; }
        }
        if (!done) printf("  => black distance > %d (or search stopped)\n", maxd);
    }
    return 0;
}
