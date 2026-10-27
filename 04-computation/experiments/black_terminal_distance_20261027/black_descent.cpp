// black_descent.cpp -- mac-mini-2026-10-27-potentials
//
// Arc-strong connectivity as the black distance to the non-strong locus.
//
// Conventions (blue_line_black_detours_20261020.cpp, endpoint packets): a legal word Q is a directed
// Hamiltonian path of T; the native complement K_Q(T) = T^op xor E(Q) keeps the arcs of Q and reverses
// every other pair; Q is BLUE when its positional reversal v_i -> v_(n+1-i) is an anti-automorphism of T,
// BLACK otherwise.
//
// lambda(T) = min over nonempty proper X of e(X -> X^c)  (arc-strong connectivity; 0 iff T is not strong).
// Lemma (proved in the note): for every legal Q, lambda(K_Q T) >= lambda(T) - 1, with equality iff some
// minimum out-cut X (e(X -> X^c) = lambda(T)) contains start(Q) and misses end(Q).
//
// For every strong class of order n (read from `gentourng -c n`, upper triangle row by row, '1' = i -> j)
// this program records lambda, whether a BLACK descent (a black Q with lambda(K_Q T) = lambda - 1) exists,
// and whether any descent exists.  With "check" it also verifies the Lemma literally on every legal path
// (lambda of the complemented tournament computed from scratch) -- intended for n <= 8.
//
// Build: g++ -O2 -std=c++17 -o black_descent black_descent.cpp
// Run:   gentourng -cq n | ./black_descent n [check]
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <map>
#include <string>
using namespace std;

static int n;

struct Tour { uint16_t out[16]; };

static int lambdaOf(const Tour &T, vector<uint16_t> *mincuts) {
    int full = (1 << n) - 1, best = 1 << 20;
    for (int X = 1; X < full; X++) {
        int f = 0;
        for (int x = 0; x < n; x++) if (X >> x & 1) f += __builtin_popcount(T.out[x] & ~X & full);
        if (f < best) { best = f; if (mincuts) mincuts->clear(); }
        if (f == best && mincuts) mincuts->push_back((uint16_t)X);
    }
    return best;
}

static inline bool isBlack(const Tour &T, const int *Q) {
    int r[16];
    for (int i = 0; i < n; i++) r[Q[i]] = Q[n - 1 - i];
    for (int u = 0; u < n; u++) {
        uint16_t o = T.out[u];
        while (o) { int v = __builtin_ctz(o); o &= o - 1; if (!((T.out[r[v]] >> r[u]) & 1)) return true; }
    }
    return false;
}

static inline Tour kmove(const Tour &T, const int *Q) {
    Tour R; uint16_t full = (uint16_t)((1u << n) - 1);
    for (int u = 0; u < n; u++) R.out[u] = (uint16_t)((~T.out[u]) & full & ~(1u << u));
    for (int i = 0; i + 1 < n; i++) { int u = Q[i], v = Q[i + 1]; R.out[u] |= (uint16_t)(1u << v); R.out[v] &= (uint16_t)~(1u << u); }
    return R;
}

// DFS over directed Hamiltonian paths starting at s; callback returns true to stop.
template <class F> static bool pathsFrom(const Tour &T, int s, F &&f) {
    int Q[16]; uint16_t cand[16]; int depth = 0; Q[0] = s; uint16_t used = (uint16_t)(1u << s);
    cand[0] = T.out[s] & ~used;
    while (depth >= 0) {
        if (depth == n - 1) { if (f(Q)) return true; used &= (uint16_t)~(1u << Q[depth]); depth--; continue; }
        if (!cand[depth]) { used &= (uint16_t)~(1u << Q[depth]); depth--; continue; }
        int v = __builtin_ctz(cand[depth]); cand[depth] &= cand[depth] - 1;
        depth++; Q[depth] = v; used |= (uint16_t)(1u << v); cand[depth] = T.out[v] & ~used;
    }
    return false;
}

int main(int argc, char **argv) {
    n = argc > 1 ? atoi(argv[1]) : 7;
    bool check = argc > 2 && strcmp(argv[2], "check") == 0;
    if (n < 3 || n > 12) { fprintf(stderr, "n in 3..12\n"); return 1; }
    char line[256];
    long long classes = 0, lemmaPaths = 0, lemmaFail = 0;
    map<int, long long> byLambda, noBlackDesc, noDesc, pureBlue;
    vector<string> exceptions;
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] != '0' && line[0] != '1') continue;
        Tour T; memset(&T, 0, sizeof T);
        int k = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
            if (line[k] == '1') T.out[i] |= (uint16_t)(1u << j); else T.out[j] |= (uint16_t)(1u << i);
            k++;
        }
        classes++;
        vector<uint16_t> cuts;
        int lam = lambdaOf(T, &cuts);
        byLambda[lam]++;
        // sep[s][e]: some minimum out-cut contains s and misses e
        bool sep[16][16]; memset(sep, 0, sizeof sep);
        for (uint16_t X : cuts) for (int s = 0; s < n; s++) if (X >> s & 1) for (int e = 0; e < n; e++) if (!(X >> e & 1)) sep[s][e] = true;
        bool anyBlack = false, blackDesc = false, anyDesc = false;
        for (int s = 0; s < n; s++) {
            bool stop = pathsFrom(T, s, [&](const int *Q) {
                bool b = isBlack(T, Q);
                bool d = sep[Q[0]][Q[n - 1]];
                if (b) anyBlack = true;
                if (d) anyDesc = true;
                if (b && d) blackDesc = true;
                if (check) {
                    Tour R = kmove(T, Q);
                    int l2 = lambdaOf(R, nullptr);
                    lemmaPaths++;
                    if (l2 < lam - 1 || (l2 == lam - 1) != d) lemmaFail++;
                    return false;
                }
                return blackDesc && anyDesc;
            });
            if (stop) break;
        }
        if (!anyBlack) pureBlue[lam]++;
        if (!blackDesc) { noBlackDesc[lam]++; if (exceptions.size() < 200) exceptions.push_back(string(line, strlen(line) - 1) + (anyBlack ? " (has black paths)" : " (pure blue)")); }
        if (!anyDesc) noDesc[lam]++;
    }
    printf("n=%d strong classes=%lld\n", n, classes);
    for (auto &p : byLambda)
        printf("  lambda=%d: classes=%lld  no-black-descent=%lld  no-descent-at-all=%lld  pure-blue=%lld\n", p.first, p.second,
               noBlackDesc[p.first], noDesc[p.first], pureBlue[p.first]);
    if (check) printf("  lemma check: %lld legal paths, %lld failures\n", lemmaPaths, lemmaFail);
    printf("  exceptions (no black descent), first %zu:\n", exceptions.size());
    for (auto &e : exceptions) printf("    %s\n", e.c_str());
    return 0;
}
