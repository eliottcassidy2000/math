// tour.hpp -- shared tournament utilities for the KEY LEMMA checkers (n <= 16).
#pragma once
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
using namespace std;

struct Tour {
    int n = 0;
    uint32_t out[20];
    int s[20];
    void finish() { for (int v = 0; v < n; v++) s[v] = __builtin_popcount(out[v]); }
    bool arc(int u, int v) const { return (out[u] >> v) & 1; }
    bool regular() const { for (int v = 0; v < n; v++) if (2 * s[v] != n - 1) return false; return true; }
    int smin() const { int m = 99; for (int v = 0; v < n; v++) m = min(m, s[v]); return m; }
    int smax() const { int m = -1; for (int v = 0; v < n; v++) m = max(m, s[v]); return m; }
};

static inline bool parseLine(const char *line, int n, Tour &T) {
    if (line[0] != '0' && line[0] != '1') return false;
    T.n = n; memset(T.out, 0, sizeof T.out);
    int k = 0;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
        if (line[k] == '1') T.out[i] |= 1u << j; else if (line[k] == '0') T.out[j] |= 1u << i; else return false;
        k++;
    }
    T.finish();
    return true;
}

static inline string encode(const Tour &T) {
    string r;
    for (int i = 0; i < T.n; i++) for (int j = i + 1; j < T.n; j++) r += T.arc(i, j) ? '1' : '0';
    return r;
}

// f(X) = number of arcs leaving X
static inline int fOut(const Tour &T, uint32_t X) {
    uint32_t full = (T.n == 32) ? ~0u : ((1u << T.n) - 1);
    int f = 0; for (int x = 0; x < T.n; x++) if (X >> x & 1) f += __builtin_popcount(T.out[x] & ~X & full);
    return f;
}

static inline int minCuts(const Tour &T, vector<uint32_t> &cuts) {
    uint32_t full = (1u << T.n) - 1; int best = 1 << 20; cuts.clear();
    for (uint32_t X = 1; X < full; X++) {
        int f = fOut(T, X);
        if (f < best) { best = f; cuts.clear(); }
        if (f == best) cuts.push_back(X);
    }
    return best;
}

static inline bool strongMask(const Tour &T, uint32_t mask) {
    if (!mask) return false;
    int st = __builtin_ctz(mask);
    uint32_t r = 1u << st, fr = r;
    while (fr) { int x = __builtin_ctz(fr); fr &= fr - 1; uint32_t nb = T.out[x] & mask & ~r; r |= nb; fr |= nb; }
    if (r != mask) return false;
    r = 1u << st; fr = r;
    while (fr) { int x = __builtin_ctz(fr); fr &= fr - 1; uint32_t nb = 0; for (int y = 0; y < T.n; y++) if ((mask >> y & 1) && !(r >> y & 1) && T.arc(y, x)) nb |= 1u << y; r |= nb; fr |= nb; }
    return r == mask;
}

// HP[a] = bitmask of b such that a Hamiltonian path a ~> b exists.
static inline void hamPairs(const Tour &T, uint32_t HP[20]) {
    int n = T.n; uint32_t full = (1u << n) - 1;
    vector<uint32_t> dp((size_t)1 << n << 5, 0);
    auto D = [&](uint32_t mask, int v) -> uint32_t & { return dp[((size_t)mask << 5) | v]; };
    for (int v = 0; v < n; v++) D(1u << v, v) = 1u << v;
    for (uint32_t mask = 1; mask <= full; mask++) for (int v = 0; v < n; v++) {
        uint32_t st = D(mask, v); if (!st) continue;
        uint32_t o = T.out[v] & ~mask;
        while (o) { int w = __builtin_ctz(o); o &= o - 1; D(mask | (1u << w), w) |= st; }
    }
    for (int a = 0; a < n; a++) HP[a] = 0;
    for (int e = 0; e < n; e++) { uint32_t st = D(full, e); while (st) { int a = __builtin_ctz(st); st &= st - 1; HP[a] |= 1u << e; } }
}

// positional reversal of Q is an anti-automorphism?
static inline bool isBlue(const Tour &T, const int *Q) {
    int n = T.n; int r[20];
    for (int i = 0; i < n; i++) r[Q[i]] = Q[n - 1 - i];
    for (int u = 0; u < n; u++) { uint32_t o = T.out[u]; while (o) { int v = __builtin_ctz(o); o &= o - 1; if (!T.arc(r[v], r[u])) return false; } }
    return true;
}

// enumerate Hamiltonian paths from a to b; callback returns true to stop. returns true if stopped.
template <class F> static bool pathsAB(const Tour &T, int a, int b, F &&f) {
    int n = T.n; int Q[20]; uint32_t cand[20]; int depth = 0; Q[0] = a; uint32_t used = 1u << a;
    if (n == 1) return f(Q);
    cand[0] = T.out[a] & ~used;
    while (depth >= 0) {
        if (depth == n - 1) { if (Q[depth] == b && f(Q)) return true; used &= ~(1u << Q[depth]); depth--; continue; }
        if (!cand[depth]) { used &= ~(1u << Q[depth]); depth--; continue; }
        int v = __builtin_ctz(cand[depth]); cand[depth] &= cand[depth] - 1;
        if (v == b && depth + 1 != n - 1) continue; // b must be last
        depth++; Q[depth] = v; used |= 1u << v; cand[depth] = T.out[v] & ~used;
    }
    return false;
}

// anti-automorphism search (score-respecting backtracking).  Returns true and fills sig if found.
struct AntiSearch {
    const Tour *T; int sig[20]; bool used[20];
    bool rec(int k) {
        int n = T->n;
        if (k == n) return true;
        for (int c = 0; c < n; c++) {
            if (used[c] || T->s[c] != n - 1 - T->s[k]) continue;
            bool ok = true;
            for (int j = 0; j < k && ok; j++) if (T->arc(k, j) != T->arc(sig[j], c)) ok = false;
            if (!ok) continue;
            sig[k] = c; used[c] = true;
            if (rec(k + 1)) return true;
            used[c] = false;
        }
        return false;
    }
    bool run(const Tour &TT) { T = &TT; memset(used, 0, sizeof used); return rec(0); }
};
static inline bool selfConverse(const Tour &T) { AntiSearch A; return A.run(T); }

// count automorphisms (backtracking)
struct AutCount {
    const Tour *T; int sig[20]; bool used[20]; long long cnt; int fix0; // fix0 >= 0: require sig(fix0)=fix0
    void rec(int k) {
        int n = T->n;
        if (k == n) { cnt++; return; }
        for (int c = 0; c < n; c++) {
            if (used[c] || T->s[c] != T->s[k]) continue;
            if (k == fix0 && c != fix0) continue;
            if (c == fix0 && k != fix0) continue;
            bool ok = true;
            for (int j = 0; j < k && ok; j++) if (T->arc(k, j) != T->arc(c, sig[j])) ok = false;
            if (!ok) continue;
            sig[k] = c; used[c] = true; rec(k + 1); used[c] = false;
        }
    }
    long long run(const Tour &TT, int fixv = -1) { T = &TT; memset(used, 0, sizeof used); cnt = 0; fix0 = fixv; rec(0); return cnt; }
};

struct LemmaResult {
    int lambda = 0;
    int ncuts = 0;
    bool strongForm = false;      // closed score-certified descent candidate exists
    bool openScore = false;       // (possibly open) score-certified candidate exists
    bool keyLemma = false;        // black descent candidate exists
    long long blueChecked = 0;    // number of candidate paths found blue during the search
};

// sep[a] = bitmask of b such that some minimum out-cut contains a and misses b
static inline LemmaResult analyse(const Tour &T, bool needKey = true) {
    LemmaResult R; int n = T.n;
    vector<uint32_t> cuts; R.lambda = minCuts(T, cuts); R.ncuts = (int)cuts.size();
    uint32_t full = (1u << n) - 1;
    uint32_t sep[20]; for (int a = 0; a < n; a++) sep[a] = 0;
    for (uint32_t X : cuts) for (int a = 0; a < n; a++) if (X >> a & 1) sep[a] |= (~X) & full;
    uint32_t HP[20]; hamPairs(T, HP);
    for (int a = 0; a < n; a++) {
        uint32_t c = HP[a] & sep[a];
        while (c) {
            int b = __builtin_ctz(c); c &= c - 1;
            if (T.s[a] + T.s[b] != n - 1) { R.openScore = true; R.keyLemma = true; if (T.arc(b, a)) R.strongForm = true; }
        }
    }
    if (!needKey || R.keyLemma) return R;
    // otherwise search for a black path among candidate pairs
    for (int a = 0; a < n && !R.keyLemma; a++) {
        uint32_t c = HP[a] & sep[a];
        while (c && !R.keyLemma) {
            int b = __builtin_ctz(c); c &= c - 1;
            pathsAB(T, a, b, [&](const int *Q) { if (!isBlue(T, Q)) { R.keyLemma = true; return true; } R.blueChecked++; return false; });
        }
    }
    return R;
}
