// aud.cpp -- independent auditor code for THM-4811 (written from the theorem's definitions only).
// Tournament on n <= 16 vertices: out[u] = bitmask of out-neighbours.
// Legal word Q = directed Hamiltonian path q_0 -> ... -> q_{n-1}.
// Native move tau_Q: keep the n-1 arcs of Q, reverse every other pair.
// Q blue iff sigma: q_i -> q_{n-1-i} is an anti-automorphism (u->v  =>  sigma(v)->sigma(u)), black otherwise.
// f(X) = e(X -> V\X) computed as sum_{x in X} s(x) - C(|X|,2); lambda = min over nonempty proper X.
//
// Modes:
//   census n            : stdin = gentourng -c n ; lambda histogram, black-descent existence, exceptions,
//                         lambda < delta0 count; at n = 10 the lambda = 2 classes must have a black descent whose
//                         target is not E10.
//   lemma n             : stdin = gentourng n (all classes) ; every Hamiltonian path: dipole law (scores, cut law
//                         for every X, c3 law) and the lambda Lipschitz / descent criterion, lambda recomputed by
//                         brute force on the target.
//   lemmalab n          : all 2^C(n,2) labelled tournaments (n <= 6), same checks.
//   bfs n               : stdin = gentourng n ; class graph with nauty, BFS distances to NS / SING / CORE with any
//                         and black moves; black-core components (merged and unmerged).
//   cobasin n k         : stdin = gentourng n ; sufficient certificate of B_n connectivity with at most k black
//                         descent words per lambda = 1 class plus one black line per non-strong class.
//   family              : the E(M;Z) members through order 16 and the K1-first variant.
//   e10                 : E10 statistics.
#define MAXN WORDSIZE
#include <nauty/nauty.h>
#include <nauty/nautinv.h>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <string>
#include <unordered_map>
#include <map>
#include <algorithm>
#include <functional>
#include <cmath>
#include <unordered_set>

typedef uint32_t msk;
static int n;
struct Tn { msk o[16]; };

static inline int popc(msk x) { return __builtin_popcount(x); }
static inline msk fullm() { return (msk)((1u << n) - 1); }

static Tn parseLine(const char *line) {
    Tn T; memset(&T, 0, sizeof T);
    int k = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++) {
            if (line[k] == '1') T.o[i] |= 1u << j; else T.o[j] |= 1u << i;
            k++;
        }
    return T;
}
static std::string code(const Tn &T) {
    std::string s;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) s += ((T.o[i] >> j) & 1) ? '1' : '0';
    return s;
}
static inline bool arc(const Tn &T, int u, int v) { return (T.o[u] >> v) & 1; }
static Tn converse(const Tn &T) {
    Tn R; memset(&R, 0, sizeof R);
    for (int u = 0; u < n; u++) for (int v = 0; v < n; v++) if (u != v && arc(T, u, v)) R.o[v] |= 1u << u;
    return R;
}
static void scores(const Tn &T, int *s) { for (int u = 0; u < n; u++) s[u] = popc(T.o[u]); }

// lambda via f(X) = sum s - C(k,2); returns lambda, fills desc[a] = set of b with some min cut containing a, missing b
static int lambdaDesc(const Tn &T, msk *desc, int *ncuts) {
    static int sumS[1 << 16];
    int s[16]; scores(T, s);
    int N = 1 << n; sumS[0] = 0;
    int best = 1 << 30, cnt = 0;
    for (int X = 1; X < N; X++) {
        int low = __builtin_ctz(X);
        sumS[X] = sumS[X & (X - 1)] + s[low];
    }
    for (int X = 1; X < N - 1; X++) {
        int k = popc(X);
        int f = sumS[X] - k * (k - 1) / 2;
        if (f < best) { best = f; cnt = 0; }
        if (f == best) cnt++;
    }
    if (desc) {
        for (int a = 0; a < n; a++) desc[a] = 0;
        for (int X = 1; X < N - 1; X++) {
            int k = popc(X);
            int f = sumS[X] - k * (k - 1) / 2;
            if (f == best) {
                msk out = fullm() & ~(msk)X;
                for (int a = 0; a < n; a++) if ((X >> a) & 1) desc[a] |= out;
            }
        }
    }
    if (ncuts) *ncuts = cnt;
    return best;
}
// brute-force f by counting arcs, used to cross-check
static int fBrute(const Tn &T, msk X) {
    int f = 0;
    for (int u = 0; u < n; u++) if ((X >> u) & 1) for (int v = 0; v < n; v++) if (!((X >> v) & 1) && arc(T, u, v)) f++;
    return f;
}
static int lambdaBrute(const Tn &T) {
    int best = 1 << 30;
    for (msk X = 1; X < fullm(); X++) best = std::min(best, fBrute(T, X));
    return best;
}
static bool strongOn(const Tn &T, msk V) {
    if (!V) return true;
    int st = __builtin_ctz(V);
    msk seen = 1u << st, fr = seen;
    while (fr) { msk nx = 0; for (msk f = fr; f; f &= f - 1) nx |= T.o[__builtin_ctz(f)] & V; nx &= ~seen; seen |= nx; fr = nx; }
    if (seen != V) return false;
    seen = 1u << st; fr = seen;
    while (fr) {
        msk nx = 0;
        for (msk w = V & ~seen; w; w &= w - 1) { int u = __builtin_ctz(w); if (T.o[u] & fr) nx |= 1u << u; }
        seen |= nx; fr = nx;
    }
    return seen == V;
}
static bool isBlue(const Tn &T, const int *Q) {
    int sg[16];
    for (int i = 0; i < n; i++) sg[Q[i]] = Q[n - 1 - i];
    for (int u = 0; u < n; u++)
        for (int v = 0; v < n; v++)
            if (u != v && arc(T, u, v) && !arc(T, sg[v], sg[u])) return false;
    return true;
}
static Tn moveQ(const Tn &T, const int *Q) {
    Tn R; memset(&R, 0, sizeof R);
    // reverse all pairs
    for (int u = 0; u < n; u++) for (int v = 0; v < n; v++) if (u != v && arc(T, u, v)) R.o[v] |= 1u << u;
    // restore path arcs
    for (int i = 0; i + 1 < n; i++) { int u = Q[i], v = Q[i + 1]; R.o[u] |= 1u << v; R.o[v] &= ~(1u << u); }
    return R;
}
// enumerate Hamiltonian paths; callback returns true to stop. startMask restricts the first vertex.
static bool hamPaths(const Tn &T, msk startMask, const std::function<bool(const int *)> &cb) {
    int Q[16];
    std::function<bool(int, msk)> rec = [&](int depth, msk used) -> bool {
        if (depth == n) return cb(Q);
        msk cand = T.o[Q[depth - 1]] & ~used;
        while (cand) {
            int v = __builtin_ctz(cand); cand &= cand - 1;
            Q[depth] = v;
            if (rec(depth + 1, used | (1u << v))) return true;
        }
        return false;
    };
    for (msk sm = startMask; sm; sm &= sm - 1) {
        int s = __builtin_ctz(sm);
        Q[0] = s;
        if (n == 1) { if (cb(Q)) return true; continue; }
        if (rec(1, 1u << s)) return true;
    }
    return false;
}
// nauty canonical code (upper triangle bits of canonical form)
static uint64_t canonCode(const Tn &T, double *grpsize = nullptr) {
    graph g[MAXN], cg[MAXN];
    int lab[MAXN], ptn[MAXN], orbits[MAXN];
    static DEFAULTOPTIONS_DIGRAPH(options);
    statsblk stats;
    options.getcanon = TRUE;
    EMPTYGRAPH(g, 1, n);
    for (int u = 0; u < n; u++) for (int v = 0; v < n; v++) if (u != v && arc(T, u, v)) ADDONEARC(g, u, v, 1);
    densenauty(g, lab, ptn, orbits, &options, &stats, 1, n, cg);
    if (grpsize) *grpsize = stats.grpsize1 * pow(10.0, stats.grpsize2);
    uint64_t c = 0; int b = 0;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (ISELEMENT(GRAPHROW(cg, i, 1), j)) c |= 1ull << b; b++; }
    return c;
}
static uint64_t mergedCode(const Tn &T) { uint64_t a = canonCode(T), b = canonCode(converse(T)); return a < b ? a : b; }
static int delta0(const Tn &T) {
    int s[16]; scores(T, s); int d = 1 << 30;
    for (int u = 0; u < n; u++) d = std::min(d, std::min(s[u], n - 1 - s[u]));
    return d;
}
static bool inSING(const Tn &T) { int s[16]; scores(T, s); for (int u = 0; u < n; u++) if (s[u] == 0 || s[u] == n - 1) return true; return false; }
static bool inCORE(const Tn &T) {
    int s[16]; scores(T, s);
    for (int u = 0; u < n; u++) if (s[u] == 0 || s[u] == n - 1) { if (n - 1 >= 3 && strongOn(T, fullm() & ~(1u << u))) return true; }
    return false;
}
static long long c3(const Tn &T) { int s[16]; scores(T, s); long long c = (long long)n * (n - 1) * (n - 2) / 6; for (int u = 0; u < n; u++) c -= (long long)s[u] * (s[u] - 1) / 2; return c; }

// ---------------------------------------------------------------------------------------------------------------
static const char *E10CODE = "111001111101111111111111111111111111000111101";
static uint64_t e10canon = 0;

// does T have a black descent word? if avoidCanon != 0, require target canon != avoidCanon (target identified exactly)
static bool hasBlackDescent(const Tn &T, int lam, const msk *desc, uint64_t avoidCanon, int avoidLam, long long *blackDescCount = nullptr) {
    msk starts = 0;
    for (int a = 0; a < n; a++) if (desc[a]) starts |= 1u << a;
    bool found = false;
    hamPaths(T, starts, [&](const int *Q) {
        if (!((desc[Q[0]] >> Q[n - 1]) & 1)) return false;
        if (isBlue(T, Q)) return false;
        if (avoidCanon && lam - 1 == avoidLam) {
            Tn R = moveQ(T, Q);
            if (canonCode(R) == avoidCanon) return false;
        }
        found = true;
        if (blackDescCount) { (*blackDescCount)++; return false; }
        return true;
    });
    return found;
}

static int modeCensus() {
    char line[512];
    std::map<int, long long> byLam, noBlack, noDesc, pureBlue;
    std::vector<std::string> exc;
    long long cls = 0, lamLtDelta = 0, sumH = 0;
    std::map<std::pair<int,int>, long long> lamDelta;
    if (n == 10) { n = 10; e10canon = canonCode(parseLine(E10CODE)); }
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] != '0' && line[0] != '1') continue;
        Tn T = parseLine(line);
        cls++;
        msk desc[16]; int nc;
        int lam = lambdaDesc(T, desc, &nc);
        int d0 = delta0(T);
        lamDelta[{lam, d0}]++;
        if (lam < d0) lamLtDelta++;
        byLam[lam]++;
        bool anyDesc = false;
        for (int a = 0; a < n; a++) if (desc[a]) anyDesc = true;
        // any descent word at all (any colour)?
        bool descWord = false;
        { msk st = 0; for (int a = 0; a < n; a++) if (desc[a]) st |= 1u << a;
          hamPaths(T, st, [&](const int *Q) { if ((desc[Q[0]] >> Q[n - 1]) & 1) { descWord = true; return true; } return false; }); }
        if (!descWord) noDesc[lam]++;
        uint64_t avoid = (n == 10 && lam == 2) ? e10canon : 0;
        bool bd = hasBlackDescent(T, lam, desc, avoid, 1);
        if (!bd) {
            noBlack[lam]++;
            bool anyBlack = false;
            hamPaths(T, fullm(), [&](const int *Q) { if (!isBlue(T, Q)) { anyBlack = true; return true; } return false; });
            if (!anyBlack) pureBlue[lam]++;
            exc.push_back(std::string(line, strcspn(line, "\n")) + (anyBlack ? " (has black words)" : " (pure blue)"));
        }
        if (n <= 8) hamPaths(T, fullm(), [&](const int *) { sumH++; return false; });
        (void)anyDesc;
    }
    printf("census n=%d strong classes=%lld\n", n, cls);
    for (auto &p : byLam) printf("  lambda=%d: %lld classes, no black descent %lld (pure blue %lld), no descent word %lld\n",
                                  p.first, p.second, noBlack[p.first], pureBlue[p.first], noDesc[p.first]);
    printf("  (lambda,delta0) histogram:"); for (auto &p : lamDelta) printf(" (%d,%d):%lld", p.first.first, p.first.second, p.second); printf("\n");
    printf("  lambda < delta0: %lld\n", lamLtDelta);
    if (n <= 8) printf("  sum of H over strong classes: %lld\n", sumH);
    for (auto &e : exc) printf("  EXCEPTION %s\n", e.c_str());
    if (n == 10) printf("  (lambda = 2 classes were required to have a black descent whose target is not E10)\n");
    return 0;
}

static long long lemmaPaths = 0, lemmaFail = 0, dipFail = 0, cutFail = 0, c3Fail = 0, identFail = 0;
static bool bruteCheck = false;
static void fArray(const Tn &T, std::vector<int> &f) {
    int s[16]; scores(T, s); int N = 1 << n; f.assign(N, 0);
    std::vector<int> sumS(N, 0);
    for (int X = 1; X < N; X++) sumS[X] = sumS[X & (X - 1)] + s[__builtin_ctz(X)];
    for (int X = 0; X < N; X++) { int k = popc(X); f[X] = sumS[X] - k * (k - 1) / 2; }
}
static void lemmaOne(const Tn &T) {
    msk desc[16]; int nc;
    int lam = lambdaDesc(T, desc, &nc);
    std::vector<int> fT, fR; fArray(T, fT);
    if (bruteCheck) { for (msk X = 1; X < fullm(); X++) if (fT[X] != fBrute(T, X)) identFail++; if (lam != lambdaBrute(T)) identFail++; }
    int s[16]; scores(T, s);
    long long c3T = c3(T);
    hamPaths(T, fullm(), [&](const int *Q) {
        lemmaPaths++;
        Tn R = moveQ(T, Q);
        int a = Q[0], b = Q[n - 1];
        int s2[16]; scores(R, s2);
        for (int v = 0; v < n; v++) {
            int exp = n - 1 - s[v] + (v == a ? 1 : 0) - (v == b ? 1 : 0);
            if (s2[v] != exp) dipFail++;
        }
        Tn RR = moveQ(R, Q);
        if (memcmp(RR.o, T.o, sizeof(msk) * n) != 0) dipFail++;
        if (isBlue(T, Q) != isBlue(R, Q)) dipFail++;
        fArray(R, fR);
        int l2 = 1 << 30;
        for (msk X = 1; X < fullm(); X++) {
            int rhs = fT[fullm() & ~X] + (int)((X >> a) & 1) - (int)((X >> b) & 1);
            if (fR[X] != rhs) cutFail++;
            if (bruteCheck && fR[X] != fBrute(R, X)) identFail++;
            l2 = std::min(l2, fR[X]);
        }
        if (c3(R) - c3T != s[a] - s[b] - 1) c3Fail++;
        bool crit = (desc[a] >> b) & 1;
        if (std::abs(l2 - lam) > 1 || ((l2 == lam - 1) != crit)) lemmaFail++;
        return false;
    });
}
static int modeLemma() {
    char line[512]; long long cls = 0;
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] != '0' && line[0] != '1') continue;
        lemmaOne(parseLine(line)); cls++;
    }
    printf("lemma n=%d classes=%lld paths=%lld: lambda/descent failures %lld, dipole failures %lld, cut-law failures %lld, c3 failures %lld, identity failures %lld (brute=%d)\n",
           n, cls, lemmaPaths, lemmaFail, dipFail, cutFail, c3Fail, identFail, (int)bruteCheck);
    return 0;
}
static int modeLemmaLab() {
    int m = n * (n - 1) / 2; long long cnt = 0;
    for (long long bits = 0; bits < (1ll << m); bits++) {
        Tn T; memset(&T, 0, sizeof T); int k = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if ((bits >> k) & 1) T.o[i] |= 1u << j; else T.o[j] |= 1u << i; k++; }
        lemmaOne(T); cnt++;
    }
    printf("lemmalab n=%d labelled tournaments=%lld paths=%lld: lambda/descent failures %lld, dipole failures %lld, cut-law failures %lld, c3 failures %lld, identity failures %lld (brute=%d)\n",
           n, cnt, lemmaPaths, lemmaFail, dipFail, cutFail, c3Fail, identFail, (int)bruteCheck);
    return 0;
}

// ---------------------------------------------------------------------------------------------------------------
static int modeBFS() {
    char line[512];
    std::vector<Tn> cl; std::unordered_map<uint64_t, int> id;
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] != '0' && line[0] != '1') continue;
        Tn T = parseLine(line);
        uint64_t c = canonCode(T);
        if (id.count(c)) { printf("duplicate class?\n"); return 1; }
        id[c] = (int)cl.size(); cl.push_back(T);
    }
    int N = (int)cl.size();
    std::vector<std::vector<std::pair<int, char>>> adj(N);
    std::vector<int> lam(N), d0(N), strong(N), sing(N), core(N), conv(N);
    std::vector<char> pureBlue(N, 1);
    long long lines = 0, black = 0;
    for (int i = 0; i < N; i++) {
        const Tn &T = cl[i];
        lam[i] = lambdaDesc(T, nullptr, nullptr);
        d0[i] = delta0(T); strong[i] = lam[i] > 0; sing[i] = inSING(T); core[i] = inCORE(T);
        conv[i] = id.at(canonCode(converse(T)));
        hamPaths(T, fullm(), [&](const int *Q) {
            Tn R = moveQ(T, Q);
            int j = id.at(canonCode(R));
            char bl = isBlue(T, Q) ? 0 : 1;
            adj[i].push_back({j, bl});
            lines++; black += bl;
            if (bl) pureBlue[i] = 0;
            return false;
        });
    }
    auto bfs = [&](const std::vector<int> &target, bool blackOnly) {
        std::vector<int> dist(N, -1); std::vector<int> q;
        for (int i = 0; i < N; i++) if (target[i]) { dist[i] = 0; q.push_back(i); }
        for (size_t h = 0; h < q.size(); h++) {
            int u = q[h];
            for (auto &e : adj[u]) { if (blackOnly && !e.second) continue; if (dist[e.first] < 0) { dist[e.first] = dist[u] + 1; q.push_back(e.first); } }
        }
        return dist;
    };
    // edges are symmetric (inverse move along the same word, same colour): verify
    long long asym = 0;
    {
        std::unordered_set<uint64_t> E;
        for (int i = 0; i < N; i++) for (auto &e : adj[i]) E.insert(((uint64_t)i * (uint64_t)N + (uint64_t)e.first) * 2 + (uint64_t)e.second);
        for (int i = 0; i < N; i++) for (auto &e : adj[i]) if (!E.count(((uint64_t)e.first * (uint64_t)N + (uint64_t)i) * 2 + (uint64_t)e.second)) asym++;
    }
    std::vector<int> NS(N); for (int i = 0; i < N; i++) NS[i] = !strong[i];
    auto dNSa = bfs(NS, false), dNSb = bfs(NS, true), dSa = bfs(sing, false), dSb = bfs(sing, true), dCa = bfs(core, false), dCb = bfs(core, true);
    long long strongCnt = 0, badNSa = 0, badNSb = 0, badSa = 0, badSb = 0, badCa = 0, badCb = 0;
    std::vector<std::string> notesB;
    for (int i = 0; i < N; i++) {
        if (dNSa[i] != lam[i]) badNSa++;
        if (!strong[i]) continue;
        strongCnt++;
        if (dNSb[i] != lam[i]) { badNSb++; notesB.push_back("NS-black " + code(cl[i]) + " lam=" + std::to_string(lam[i]) + " d=" + std::to_string(dNSb[i])); }
        if (dSa[i] != d0[i]) badSa++;
        if (dSb[i] != d0[i]) { badSb++; notesB.push_back("SING-black " + code(cl[i]) + " d0=" + std::to_string(d0[i]) + " d=" + std::to_string(dSb[i])); }
        if (dCa[i] != d0[i]) badCa++;
        if (dCb[i] != d0[i]) { badCb++; notesB.push_back("CORE-black " + code(cl[i]) + " d0=" + std::to_string(d0[i]) + " d=" + std::to_string(dCb[i])); }
    }
    // black components: unmerged and merged, among non-pure-blue classes
    std::vector<int> par(N); for (int i = 0; i < N; i++) par[i] = i;
    std::function<int(int)> fd = [&](int x) { while (par[x] != x) { par[x] = par[par[x]]; x = par[x]; } return x; };
    for (int i = 0; i < N; i++) for (auto &e : adj[i]) if (e.second) { int a = fd(i), b = fd(e.first); if (a != b) par[a] = b; }
    std::map<int, int> compU; int npb = 0;
    for (int i = 0; i < N; i++) if (!pureBlue[i]) compU[fd(i)]++; else npb++;
    for (int i = 0; i < N; i++) { int a = fd(i), b = fd(conv[i]); if (a != b) par[a] = b; }
    std::map<int, int> compM; for (int i = 0; i < N; i++) if (!pureBlue[i]) compM[fd(i)]++;
    long long nsMerged = 0, nsBlackMerged = 0;
    for (int i = 0; i < N; i++) if (!strong[i] && i <= conv[i]) { nsMerged++; if (!pureBlue[i]) nsBlackMerged++; }
    printf("bfs n=%d classes=%d strong=%lld lines(labelled class-words)=%lld black=%lld asymmetric-edges=%lld pure-blue classes=%d\n",
           n, N, strongCnt, lines, black, asym, npb);
    printf("  d_any(NS) != lambda (all classes): %lld | strong: d_black(NS) != lambda: %lld\n", badNSa, badNSb);
    printf("  strong: d_any(SING) != delta0: %lld, d_black(SING) != delta0: %lld\n", badSa, badSb);
    printf("  strong: d_any(CORE) != delta0: %lld, d_black(CORE) != delta0: %lld\n", badCa, badCb);
    printf("  black components among non-pure-blue classes: unmerged %zu, merged %zu; non-strong merged classes %lld, with a black line %lld\n",
           compU.size(), compM.size(), nsMerged, nsBlackMerged);
    for (auto &s : notesB) printf("  %s\n", s.c_str());
    return 0;
}

// ---------------------------------------------------------------------------------------------------------------
static int modeCobasin(int k) {
    char line[512];
    std::unordered_map<uint64_t, int> id; std::vector<int> par; std::vector<char> kind; // kind 1 = lambda1 strong, 2 = non-strong
    auto getId = [&](uint64_t c, char kd) { auto it = id.find(c); if (it != id.end()) { if (kind[it->second] != kd) kind[it->second] = 3; return it->second; } int x = (int)par.size(); id[c] = x; par.push_back(x); kind.push_back(kd); return x; };
    std::function<int(int)> fd = [&](int x) { while (par[x] != x) { par[x] = par[par[x]]; x = par[x]; } return x; };
    auto un = [&](int a, int b) { a = fd(a); b = fd(b); if (a != b) par[a] = b; };
    long long lam1 = 0, lam1NoBD = 0, nsCls = 0, nsPure = 0, edgesUsed = 0;
    std::unordered_map<uint64_t, char> nsSeen; // merged non-strong classes and whether pure blue
    std::vector<std::string> noBD;
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] != '0' && line[0] != '1') continue;
        Tn T = parseLine(line);
        msk desc[16]; int nc;
        int lam = lambdaDesc(T, desc, &nc);
        if (lam == 0) {
            nsCls++;
            uint64_t mc = mergedCode(T);
            if (nsSeen.count(mc)) continue;
            // one black line (target is strong lambda = 1 by THM-4729 L3 + Lipschitz; checked)
            int me = -1; bool any = false;
            hamPaths(T, fullm(), [&](const int *Q) {
                if (isBlue(T, Q)) return false;
                Tn R = moveQ(T, Q);
                int l2 = lambdaDesc(R, nullptr, nullptr);
                if (l2 != 1) { printf("UNEXPECTED: black line from non-strong class to lambda=%d\n", l2); }
                me = getId(mc, 2);
                int t = getId(mergedCode(R), 1);
                un(me, t); edgesUsed++; any = true;
                return true;
            });
            nsSeen[mc] = any ? 1 : 0;
            if (!any) nsPure++;
            continue;
        }
        if (lam != 1) continue;
        lam1++;
        uint64_t mc = mergedCode(T);
        int me = -1, used = 0;
        msk starts = 0; for (int a = 0; a < n; a++) if (desc[a]) starts |= 1u << a;
        hamPaths(T, starts, [&](const int *Q) {
            if (!((desc[Q[0]] >> Q[n - 1]) & 1)) return false;
            if (isBlue(T, Q)) return false;
            Tn R = moveQ(T, Q);
            if (me < 0) me = getId(mc, 1);
            int t = getId(mergedCode(R), 2);
            un(me, t); edgesUsed++; used++;
            return used >= k;
        });
        if (me < 0) { lam1NoBD++; noBD.push_back(std::string(line, strcspn(line, "\n"))); }
    }
    std::map<int, int> comps; long long touchedNS = 0, touchedL1 = 0, mixed = 0;
    for (auto &p : id) { comps[fd(p.second)]++; if (kind[p.second] == 2) touchedNS++; else if (kind[p.second] == 1) touchedL1++; else mixed++; }
    long long nsMerged = (long long)nsSeen.size(), nsBlack = 0; for (auto &p : nsSeen) nsBlack += p.second;
    printf("cobasin n=%d k=%d: non-strong classes %lld (merged %lld, with a black line %lld, pure blue %lld); lambda=1 classes %lld, without black descent %lld;\n"
           "  certificate edges used %lld; touched merged nodes: non-strong %lld, lambda=1 %lld, conflicting kind %lld; components %zu\n",
           n, k, nsCls, nsMerged, nsBlack, nsPure, lam1, lam1NoBD, edgesUsed, touchedNS, touchedL1, mixed, comps.size());
    for (auto &s : noBD) printf("  lambda=1 without black descent: %s\n", s.c_str());
    return 0;
}

// ---------------------------------------------------------------------------------------------------------------
// stacks: sequence of blocks 1 (K1) or 3 (C3); ordered join
struct Blocks { std::vector<int> b; int size() const { int s = 0; for (int x : b) s += x; return s; } };
static std::vector<Blocks> stacksOfSize(int k) {
    std::vector<Blocks> res;
    if (k == 0) { res.push_back(Blocks()); return res; }
    if (k >= 1) for (auto &r : stacksOfSize(k - 1)) { Blocks B; B.b.push_back(1); for (int x : r.b) B.b.push_back(x); res.push_back(B); }
    if (k >= 3) for (auto &r : stacksOfSize(k - 3)) { Blocks B; B.b.push_back(3); for (int x : r.b) B.b.push_back(x); res.push_back(B); }
    return res;
}
static std::string stackName(const Blocks &B) { std::string s; for (size_t i = 0; i < B.b.size(); i++) { if (i) s += "=>"; s += B.b[i] == 1 ? "K1" : "C3"; } return s; }
// build stack on vertices base..base+size-1 into adjacency matrix A
static void buildStack(std::vector<std::vector<int>> &A, int base, const Blocks &B) {
    std::vector<std::pair<int,int>> blk; int p = base;
    for (int x : B.b) { blk.push_back({p, x}); p += x; }
    for (size_t i = 0; i < blk.size(); i++) {
        int s = blk[i].first, x = blk[i].second;
        if (x == 3) { A[s][s + 1] = 1; A[s + 1][s + 2] = 1; A[s + 2][s] = 1; }
        for (size_t j = i + 1; j < blk.size(); j++)
            for (int u = s; u < s + x; u++) for (int v = blk[j].first; v < blk[j].first + blk[j].second; v++) A[u][v] = 1;
    }
}
// E(M;Z): D = x0, M, Z on 0..N-1 (x0 = 0), D* = N..2N-1 converse copy, all D -> D* except y0 = N -> x0 = 0
static Tn buildE(const Blocks &M, const Blocks &Z, int *Nout) {
    int m = M.size(), z = Z.size(), N = 1 + m + z;
    std::vector<std::vector<int>> A(2 * N, std::vector<int>(2 * N, 0));
    for (int v = 1; v <= m; v++) A[0][v] = 1;           // x0 -> M
    for (int v = 1 + m; v < N; v++) A[v][0] = 1;        // Z -> x0
    for (int u = 1; u <= m; u++) for (int v = 1 + m; v < N; v++) A[u][v] = 1; // M -> Z
    buildStack(A, 1, M); buildStack(A, 1 + m, Z);
    for (int u = 0; u < N; u++) for (int v = 0; v < N; v++) if (A[u][v]) A[N + v][N + u] = 1; // converse copy
    for (int u = 0; u < N; u++) for (int v = N; v < 2 * N; v++) A[u][v] = 1;
    A[0][N] = 0; A[N][0] = 1;
    n = 2 * N; *Nout = N;
    Tn T; memset(&T, 0, sizeof T);
    for (int u = 0; u < n; u++) for (int v = 0; v < n; v++) if (u != v && A[u][v]) T.o[u] |= 1u << v;
    // sanity: tournament
    for (int u = 0; u < n; u++) for (int v = u + 1; v < n; v++) if (arc(T, u, v) == arc(T, v, u)) { printf("not a tournament!\n"); exit(1); }
    return T;
}
static bool blackDescentExists(const Tn &T) { msk desc[16]; int nc; int lam = lambdaDesc(T, desc, &nc); return hasBlackDescent(T, lam, desc, 0, -1); }

static int modeFamily() {
    printf("E(M;Z) members (M = C3 => M', Z nonempty stack), order 2(1+|M|+|Z|) <= 16\n");
    for (int order = 10; order <= 16; order += 2) {
        int s = order / 2 - 1; // |M| + |Z|
        for (int m = 3; m <= s - 1; m++) {
            int z = s - m;
            for (auto &Mr : stacksOfSize(m - 3)) {
                Blocks M; M.b.push_back(3); for (int x : Mr.b) M.b.push_back(x);
                for (auto &Z : stacksOfSize(z)) {
                    int N; Tn T = buildE(M, Z, &N);
                    msk desc[16]; int nc; int lam = lambdaDesc(T, desc, &nc);
                    msk Dstar = (msk)(((1u << N) - 1) << N);
                    // unique min cut is D*: then desc[a] = complement of D* for a in D*, 0 otherwise
                    bool uniq = (nc == 1);
                    for (int a = 0; a < n; a++) { msk expect = ((Dstar >> a) & 1) ? (fullm() & ~Dstar) : 0; if (desc[a] != expect) uniq = false; }
                    bool st = strongOn(T, fullm());
                    long long H = 0, nblackW = 0, ndesc = 0, ndescBlue = 0, ndescClosed = 0;
                    hamPaths(T, fullm(), [&](const int *Q) {
                        H++;
                        bool blue = isBlue(T, Q);
                        if (!blue) nblackW++;
                        if ((desc[Q[0]] >> Q[n - 1]) & 1) { ndesc++; if (blue) ndescBlue++; if (arc(T, Q[n - 1], Q[0])) ndescClosed++; }
                        return false;
                    });
                    // black distance 2 certificate: black move to lambda = 1 tournament with a black descent
                    bool cert = false;
                    if (order <= 16) {
                        hamPaths(T, fullm(), [&](const int *Q) {
                            if (isBlue(T, Q)) return false;
                            Tn R = moveQ(T, Q);
                            msk d2[16]; int nc2; int l2 = lambdaDesc(R, d2, &nc2);
                            if (l2 != 1) return false;
                            int sv = n; bool bd = hasBlackDescent(R, l2, d2, 0, -1); n = sv;
                            if (bd) { cert = true; return true; }
                            return false;
                        });
                    }
                    double gs; canonCode(T, &gs);
                    printf("  E(%s ; %s) n=%d strong=%d lambda=%d mincuts=%d uniqueD*=%d H=%lld black words=%lld descent words=%lld (blue %lld, closed %lld) |Aut|=%.0f  black-distance-2 certificate=%d\n",
                           stackName(M).c_str(), stackName(Z).c_str(), n, st, lam, nc, uniq, H, nblackW, ndesc, ndescBlue, ndescClosed, gs, cert);
                }
            }
        }
    }
    // K1-first variants at order 12: M = K1 => C3 and others with |M| >= 3 starting with K1
    printf("variants with M starting with K1 (not in the family):\n");
    for (int order = 10; order <= 12; order += 2) {
        int s = order / 2 - 1;
        for (int m = 3; m <= s - 1; m++) {
            int z = s - m;
            for (auto &Mr : stacksOfSize(m - 1)) {
                Blocks M; M.b.push_back(1); for (int x : Mr.b) M.b.push_back(x);
                for (auto &Z : stacksOfSize(z)) {
                    int N; Tn T = buildE(M, Z, &N);
                    msk desc[16]; int nc; int lam = lambdaDesc(T, desc, &nc);
                    long long ndesc = 0, nbd = 0;
                    hamPaths(T, fullm(), [&](const int *Q) { if ((desc[Q[0]] >> Q[n - 1]) & 1) { ndesc++; if (!isBlue(T, Q)) nbd++; } return false; });
                    printf("  E(%s ; %s) n=%d lambda=%d mincuts=%d descent words=%lld black descent words=%lld\n", stackName(M).c_str(), stackName(Z).c_str(), n, lam, nc, ndesc, nbd);
                }
            }
        }
    }
    // E(C3;K1) vs E10 code
    { Blocks M; M.b.push_back(3); Blocks Z; Z.b.push_back(1); int N; Tn T = buildE(M, Z, &N);
      uint64_t a = canonCode(T); Tn E = parseLine(E10CODE); uint64_t b = canonCode(E);
      printf("E(C3;K1) isomorphic to the E10 code: %s\n", a == b ? "YES" : "NO"); }
    return 0;
}

static int modeE10() {
    n = 10; Tn E = parseLine(E10CODE);
    msk desc[16]; int nc; int lam = lambdaDesc(E, desc, &nc);
    int s[16]; scores(E, s); std::vector<int> sv(s, s + n); std::sort(sv.begin(), sv.end());
    printf("E10: scores"); for (int x : sv) printf(" %d", x); printf("; lambda=%d mincuts=%d strong=%d\n", lam, nc, strongOn(E, fullm()));
    double gs; uint64_t ce = canonCode(E, &gs);
    // anti-automorphisms = isomorphisms E -> converse(E): count = |Aut| if E self-converse
    uint64_t cc = canonCode(converse(E));
    printf("  |Aut| = %.0f, self-converse = %d (so #anti-automorphisms = %s)\n", gs, ce == cc, ce == cc ? "|Aut|" : "0");
    long long H = 0, black = 0, nd = 0, ndBlue = 0, ndClosed = 0;
    std::map<std::pair<int,int>, long long> moves; // (lambda', hasBlackDescent) for black moves
    long long lam2nb = 0, lam2bad = 0;
    std::vector<Tn> lam2;
    hamPaths(E, fullm(), [&](const int *Q) {
        H++;
        bool blue = isBlue(E, Q);
        if (!blue) black++;
        if ((desc[Q[0]] >> Q[n - 1]) & 1) { nd++; if (blue) ndBlue++; if (arc(E, Q[n - 1], Q[0])) ndClosed++; }
        Tn R = moveQ(E, Q);
        msk d2[16]; int nc2; int l2 = lambdaDesc(R, d2, &nc2);
        if (!blue) { bool bd = (l2 >= 1) ? hasBlackDescent(R, l2, d2, 0, -1) : false; moves[{l2, bd ? 1 : 0}]++; }
        if (l2 == 2) lam2.push_back(R);
        return false;
    });
    printf("  H=%lld black=%lld descent words=%lld (blue %lld, closed %lld)\n", H, black, nd, ndBlue, ndClosed);
    printf("  black moves by (lambda', target has black descent):"); for (auto &p : moves) printf(" (%d,%d):%lld", p.first.first, p.first.second, p.second); printf("\n");
    for (auto &R : lam2) {
        msk d2[16]; int nc2; int l2 = lambdaDesc(R, d2, &nc2);
        lam2nb++;
        if (!hasBlackDescent(R, l2, d2, ce, 1)) lam2bad++;
    }
    printf("  lambda=2 neighbours (any colour): %lld; without a black descent avoiding E10 (exact canonical test): %lld\n", lam2nb, lam2bad);
    return 0;
}


// D2 by dynamic programming over delta0 levels: good[i] iff class i reaches the target locus in exactly delta0(i)
// moves of the allowed colour (by the proven 1-Lipschitz property each step must lower delta0 by one).
static int modeD2(bool blackOnly, bool coreTarget) {
    char line[512];
    std::vector<Tn> cl; std::unordered_map<uint64_t, int> id;
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] != '0' && line[0] != '1') continue;
        Tn T = parseLine(line); id[canonCode(T)] = (int)cl.size(); cl.push_back(T);
    }
    int N = (int)cl.size(), maxd = 0;
    std::vector<int> d0(N), good(N, 0), strong(N);
    for (int i = 0; i < N; i++) { d0[i] = delta0(cl[i]); maxd = std::max(maxd, d0[i]); strong[i] = strongOn(cl[i], fullm()); }
    for (int i = 0; i < N; i++) if (d0[i] == 0) good[i] = coreTarget ? inCORE(cl[i]) : inSING(cl[i]);
    long long nauty = 0;
    for (int k = 1; k <= maxd; k++) {
        for (int i = 0; i < N; i++) {
            if (d0[i] != k) continue;
            const Tn &T = cl[i];
            int s[16]; scores(T, s);
            hamPaths(T, fullm(), [&](const int *Q) {
                // target scores without building it: s'(v) = n-1-s(v) + [v = start] - [v = end]
                int dd = 1 << 30;
                for (int v = 0; v < n; v++) { int sv = n - 1 - s[v] + (v == Q[0]) - (v == Q[n - 1]); dd = std::min(dd, std::min(sv, n - 1 - sv)); }
                if (dd != k - 1) return false;
                if (blackOnly && isBlue(T, Q)) return false;
                Tn R = moveQ(T, Q); nauty++;
                int j = id.at(canonCode(R));
                if (good[j]) { good[i] = 1; return true; }
                return false;
            });
        }
    }
    long long st = 0, stGood = 0; std::map<int, std::pair<long long,long long>> byd;
    std::vector<std::string> bad;
    for (int i = 0; i < N; i++) if (strong[i]) { st++; stGood += good[i]; byd[d0[i]].first++; byd[d0[i]].second += good[i]; if (!good[i] && bad.size() < 20) bad.push_back(code(cl[i])); }
    printf("D2-dp n=%d %s moves, target %s: strong classes %lld, reaching the target in exactly delta0 moves %lld (nauty calls %lld); by delta0:",
           n, blackOnly ? "black" : "any", coreTarget ? "CORE" : "SING", st, stGood, nauty);
    for (auto &p : byd) printf(" %d:%lld/%lld", p.first, p.second.second, p.second.first);
    printf("\n");
    for (auto &b : bad) printf("  fails: %s\n", b.c_str());
    return 0;
}

// memory-light D2 dynamic programme (sorted canonical codes instead of a hash map), for n = 10
static uint64_t packT(const Tn &T) { uint64_t c = 0; int b = 0; for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (arc(T, i, j)) c |= 1ull << b; b++; } return c; }
static Tn unpackT(uint64_t c) { Tn T; memset(&T, 0, sizeof T); int b = 0; for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if ((c >> b) & 1) T.o[i] |= 1u << j; else T.o[j] |= 1u << i; b++; } return T; }
static int modeD2big() {
    char line[512];
    std::vector<std::pair<uint64_t, uint64_t>> v; v.reserve(10000000);
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] != '0' && line[0] != '1') continue;
        Tn T = parseLine(line); v.push_back({canonCode(T), packT(T)});
    }
    std::sort(v.begin(), v.end());
    size_t N = v.size();
    std::vector<uint8_t> d0(N), good(N, 0), strong(N);
    int maxd = 0;
    for (size_t i = 0; i < N; i++) { Tn T = unpackT(v[i].second); d0[i] = (uint8_t)delta0(T); maxd = std::max(maxd, (int)d0[i]); strong[i] = strongOn(T, fullm()); if (d0[i] == 0) good[i] = inCORE(T); }
    auto find = [&](uint64_t c) -> size_t { auto it = std::lower_bound(v.begin(), v.end(), std::make_pair(c, (uint64_t)0)); if (it == v.end() || it->first != c) { printf("target not found!\n"); exit(1); } return (size_t)(it - v.begin()); };
    long long nauty = 0;
    for (int k = 1; k <= maxd; k++) {
        for (size_t i = 0; i < N; i++) {
            if (d0[i] != k) continue;
            Tn T = unpackT(v[i].second);
            int s[16]; scores(T, s);
            // a delta0-lowering move must start at a vertex with s = k (out-degree attains delta0) or end at one with
            // in-degree k (proved: these are the only moves raising max |d|); enumerate those paths only
            msk lowS = 0, highS = 0;
            for (int w = 0; w < n; w++) { if (s[w] == k) lowS |= 1u << w; if (n - 1 - s[w] == k) highS |= 1u << w; }
            Tn Rv = converse(T);
            auto tryQ = [&](const int *Q) -> bool {
                int dd = 1 << 30;
                for (int w = 0; w < n; w++) { int sv = n - 1 - s[w] + (w == Q[0]) - (w == Q[n - 1]); dd = std::min(dd, std::min(sv, n - 1 - sv)); }
                if (dd != k - 1) return false;
                if (isBlue(T, Q)) return false;
                Tn R = moveQ(T, Q); nauty++;
                if (good[find(canonCode(R))]) { good[i] = 1; return true; }
                return false;
            };
            bool done = hamPaths(T, lowS, tryQ);
            if (!done && highS) hamPaths(Rv, highS, [&](const int *Qr) { int Q[16]; for (int w = 0; w < n; w++) Q[w] = Qr[n - 1 - w]; return tryQ(Q); });
        }
    }
    long long st = 0, ok = 0; std::map<int, std::pair<long long,long long>> byd; std::vector<std::string> bad;
    for (size_t i = 0; i < N; i++) if (strong[i]) { st++; ok += good[i]; byd[d0[i]].first++; byd[d0[i]].second += good[i]; if (!good[i] && bad.size() < 20) bad.push_back(code(unpackT(v[i].second))); }
    printf("D2-dp(big) n=%d black moves, target CORE: classes %zu, strong %lld, reaching CORE in exactly delta0 black moves %lld (nauty calls in search %lld); by delta0:", n, N, st, ok, nauty);
    for (auto &p : byd) printf(" %d:%lld/%lld", p.first, p.second.second, p.second.first);
    printf("\n");
    for (auto &b : bad) printf("  fails: %s\n", b.c_str());
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage\n"); return 1; }
    std::string mode = argv[1];
    if (argc > 2) n = atoi(argv[2]);
    if (mode == "census") return modeCensus();
    if (argc > 3 && std::string(argv[3]) == "brute") bruteCheck = true;
    if (mode == "lemma") return modeLemma();
    if (mode == "lemmalab") return modeLemmaLab();
    if (mode == "bfs") return modeBFS();
    if (mode == "cobasin") return modeCobasin(argc > 3 ? atoi(argv[3]) : 2);
    if (mode == "family") return modeFamily();
    if (mode == "e10") return modeE10();
    if (mode == "canon") { char line[512]; while (fgets(line, sizeof line, stdin)) { if (line[0] != '0' && line[0] != '1') continue; Tn T = parseLine(line); printf("%llu\n", (unsigned long long)canonCode(T)); } return 0; }
    if (mode == "d2big") return modeD2big();
    if (mode == "d2") return modeD2(std::string(argv[3]) == "black", std::string(argv[4]) == "core");
    fprintf(stderr, "unknown mode\n"); return 1;
}
