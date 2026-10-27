// black_descent_profile.cpp -- mac-mini-2026-10-27-potentials
// Profiles the black descents of black_descent.cpp to guide the all-order proof:
//   per strong class: number of descent candidates (legal Q with start in some minimum out-cut X and end outside),
//   number of BLACK candidates, and min over minimum out-cuts X of the number of black candidates through X;
//   also whether every minimum out-cut X admits a black candidate that is a rotation of a Hamiltonian cycle
//   (closed word: end(Q) -> start(Q)), and whether the class is self-converse-capable (has a blue path at all).
// Build: g++ -O2 -std=c++17 -o black_descent_profile black_descent_profile.cpp ; run: gentourng -cq n | ./black_descent_profile n
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
using namespace std;
static int n;
struct Tour { uint16_t out[16]; };
static int lambdaOf(const Tour &T, vector<uint16_t> &cuts) {
    int full = (1 << n) - 1, best = 1 << 20;
    for (int X = 1; X < full; X++) {
        int f = 0; for (int x = 0; x < n; x++) if (X >> x & 1) f += __builtin_popcount(T.out[x] & ~X & full);
        if (f < best) { best = f; cuts.clear(); }
        if (f == best) cuts.push_back((uint16_t)X);
    }
    return best;
}
static inline bool isBlack(const Tour &T, const int *Q) {
    int r[16]; for (int i = 0; i < n; i++) r[Q[i]] = Q[n - 1 - i];
    for (int u = 0; u < n; u++) { uint16_t o = T.out[u]; while (o) { int v = __builtin_ctz(o); o &= o - 1; if (!((T.out[r[v]] >> r[u]) & 1)) return true; } }
    return false;
}
template <class F> static void allPaths(const Tour &T, F &&f) {
    for (int s = 0; s < n; s++) {
        int Q[16]; uint16_t cand[16]; int depth = 0; Q[0] = s; uint16_t used = (uint16_t)(1u << s); cand[0] = T.out[s] & ~used;
        while (depth >= 0) {
            if (depth == n - 1) { f(Q); used &= (uint16_t)~(1u << Q[depth]); depth--; continue; }
            if (!cand[depth]) { used &= (uint16_t)~(1u << Q[depth]); depth--; continue; }
            int v = __builtin_ctz(cand[depth]); cand[depth] &= cand[depth] - 1;
            depth++; Q[depth] = v; used |= (uint16_t)(1u << v); cand[depth] = T.out[v] & ~used;
        }
    }
}
int main(int argc, char **argv) {
    n = argc > 1 ? atoi(argv[1]) : 7;
    char line[256]; long long classes = 0;
    long long everyCutBlack = 0, everyCutClosedBlack = 0, hasBlue = 0;
    int minBlackCand = 1 << 30, minBlackPerCut = 1 << 30; string argMinCand, argMinCut;
    map<int, long long> histMinPerCut;
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] != '0' && line[0] != '1') continue;
        Tour T; memset(&T, 0, sizeof T); int k = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (line[k] == '1') T.out[i] |= (uint16_t)(1u << j); else T.out[j] |= (uint16_t)(1u << i); k++; }
        classes++;
        vector<uint16_t> cuts; int lam = lambdaOf(T, cuts); (void)lam;
        int C = (int)cuts.size();
        vector<int> blackThrough(C, 0), closedBlackThrough(C, 0);
        int cand = 0, blackCand = 0; bool blue = false;
        allPaths(T, [&](const int *Q) {
            bool b = isBlack(T, Q); if (!b) blue = true;
            bool closed = (T.out[Q[n - 1]] >> Q[0]) & 1;
            bool isCand = false;
            for (int c = 0; c < C; c++) if ((cuts[c] >> Q[0] & 1) && !(cuts[c] >> Q[n - 1] & 1)) {
                isCand = true; if (b) { blackThrough[c]++; if (closed) closedBlackThrough[c]++; }
            }
            if (isCand) { cand++; if (b) blackCand++; }
        });
        if (blue) hasBlue++;
        int mpc = *min_element(blackThrough.begin(), blackThrough.end());
        int mcpc = *min_element(closedBlackThrough.begin(), closedBlackThrough.end());
        histMinPerCut[min(mpc, 10)]++;
        if (mpc > 0) everyCutBlack++;
        if (mcpc > 0) everyCutClosedBlack++;
        string code(line, strlen(line) - 1);
        if (blackCand < minBlackCand) { minBlackCand = blackCand; argMinCand = code; }
        if (mpc < minBlackPerCut) { minBlackPerCut = mpc; argMinCut = code; }
    }
    printf("n=%d strong classes=%lld  with a blue path=%lld\n", n, classes, hasBlue);
    printf("  every minimum out-cut has a black candidate: %lld / %lld\n", everyCutBlack, classes);
    printf("  every minimum out-cut has a CLOSED black candidate (rotation of a Ham cycle): %lld / %lld\n", everyCutClosedBlack, classes);
    printf("  min #black candidates over classes = %d at %s\n", minBlackCand, argMinCand.c_str());
    printf("  min over (class, cut) of #black candidates through the cut = %d at %s\n", minBlackPerCut, argMinCut.c_str());
    printf("  histogram of min-per-cut black candidates (capped at 10):");
    for (auto &p : histMinPerCut) printf(" %d:%lld", p.first, p.second);
    printf("\n");
    return 0;
}
