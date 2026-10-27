// dipole_exactness.cpp -- mac-mini-2026-10-27-dipole
// Exactness of every repository potential bound for native-move distances, on isomorphism classes of order n.
// Moves: tau_Q(T) = T^op xor E(Q) for a Hamiltonian path Q (K -> -K + 2A_Q); black unless the positional reversal of Q is an
// anti-automorphism.  The class graph is undirected (moves are involutions; colour is preserved).  For each target locus
// and potential we compute the exact distance (multi-source BFS, any colour / black only), the bound, the classes where
// the bound is not attained (failures), and the local minima (outside the target, no neighbour one potential step lower).
// Every failure and local minimum is printed with its symmetry data (|Aut|, self-converse), strong components, scores.
//   targets/potentials:  TT  (transitive)          L  = (1/2) sum |s_(i) - i|          (THM-4810 Theorem L)
//                        TTf (transitive)          ceil(fas/(n-1))                      (THM-4810 (C))
//                        NS  (non-strong)          lambda (arc-strong connectivity)     (THM-4811 (B))
//                        SING (source or sink)     delta0 = min_v min(s_v, n-1-s_v)     (THM-4811 (D))
//                        CORE (singleton + strong core)  delta0                          (THM-4811 (D2))
//                        REG (regular / almost regular)  R = (1/2) sum_v |s_v - (n-1)/2| (odd n), ceil of (1/2) sum_v (|2s_v-n+1|-1)/2 (even n)
// Build: g++ -O3 -march=native -std=c++17 -I/opt/homebrew/opt/nauty/include -o dipole_exactness dipole_exactness.cpp /opt/homebrew/opt/nauty/lib/libnauty1.a
// Run:   gentourng -q n | ./dipole_exactness n [details]
#define MAXN WORDSIZE
#include <nauty/nauty.h>
#include <nauty/nautinv.h>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <unordered_map>
#include <map>
#include <deque>
#include <string>
#include <algorithm>
static int n;
struct Tour { uint32_t out[16]; };
static uint64_t canonCode(const Tour &T, double *grp = nullptr) {
    graph g[MAXN], cg[MAXN]; int lab[MAXN], ptn[MAXN], orbits[MAXN];
    static DEFAULTOPTIONS_DIGRAPH(options); statsblk stats; options.getcanon = TRUE;
    EMPTYGRAPH(g, 1, n);
    for (int u = 0; u < n; u++) for (uint32_t o = T.out[u]; o; o &= o - 1) ADDONEARC(g, u, __builtin_ctz(o), 1);
    densenauty(g, lab, ptn, orbits, &options, &stats, 1, n, cg);
    uint64_t code = 0; int b = 0;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (ISELEMENT(GRAPHROW(cg, i, 1), j)) code |= 1ull << b; b++; }
    if (grp) { double s = stats.grpsize1; for (int e = 0; e < stats.grpsize2; e++) s *= 10; *grp = s; }
    return code;
}
static Tour converse(const Tour &T) { Tour R; uint32_t full = (1u << n) - 1; for (int u = 0; u < n; u++) R.out[u] = (~T.out[u]) & full & ~(1u << u); return R; }
static bool strongMask(const Tour &T, uint32_t vs) {
    if (!vs) return true; int start = __builtin_ctz(vs);
    uint32_t seen = 1u << start, fr = seen;
    while (fr) { uint32_t nx = 0; for (uint32_t f = fr; f; f &= f - 1) nx |= T.out[__builtin_ctz(f)] & vs; nx &= ~seen; seen |= nx; fr = nx; }
    if (seen != vs) return false;
    seen = 1u << start; fr = seen;
    while (fr) { uint32_t nx = 0; for (uint32_t f = fr; f; f &= f - 1) { int v = __builtin_ctz(f); for (uint32_t w = vs; w; w &= w - 1) { int u = __builtin_ctz(w); if (T.out[u] >> v & 1) nx |= 1u << u; } } nx &= ~seen; seen |= nx; fr = nx; }
    return seen == vs;
}
static inline bool isBlack(const Tour &T, const int *Q) {
    int r[16]; for (int i = 0; i < n; i++) r[Q[i]] = Q[n - 1 - i];
    for (int u = 0; u < n; u++) for (uint32_t o = T.out[u]; o; o &= o - 1) { int v = __builtin_ctz(o); if (!((T.out[r[v]] >> r[u]) & 1)) return true; }
    return false;
}
static Tour kmove(const Tour &T, const int *Q) {
    Tour R; uint32_t full = (1u << n) - 1;
    for (int u = 0; u < n; u++) R.out[u] = (~T.out[u]) & full & ~(1u << u);
    for (int i = 0; i + 1 < n; i++) { int u = Q[i], v = Q[i + 1]; R.out[u] |= 1u << v; R.out[v] &= ~(1u << u); }
    return R;
}
template <class F> static void paths(const Tour &T, F &&f) {
    for (int s = 0; s < n; s++) {
        int Q[16]; uint32_t cand[16]; int depth = 0; Q[0] = s; uint32_t used = 1u << s; cand[0] = T.out[s] & ~used;
        while (depth >= 0) {
            if (depth == n - 1) { f(Q); used &= ~(1u << Q[depth]); depth--; continue; }
            if (!cand[depth]) { used &= ~(1u << Q[depth]); depth--; continue; }
            int v = __builtin_ctz(cand[depth]); cand[depth] &= cand[depth] - 1;
            depth++; Q[depth] = v; used |= 1u << v; cand[depth] = T.out[v] & ~used;
        }
    }
}
static int fasDP(const Tour &T) {
    int full = (1 << n) - 1; static int f[1 << 13];
    f[0] = 0;
    for (int S = 1; S <= full; S++) { int best = 1 << 30; for (uint32_t s = S; s; s &= s - 1) { int v = __builtin_ctz(s); int c = f[S & ~(1 << v)] + __builtin_popcount(T.out[v] & (S & ~(1 << v))); if (c < best) best = c; } f[S] = best; }
    return f[full];
}
static std::string compStr(const Tour &T) {   // sizes of strong components, top first
    uint32_t full = (1u << n) - 1; int sc[16], v[16]; for (int i = 0; i < n; i++) { v[i] = i; sc[i] = __builtin_popcount(T.out[i]); }
    std::sort(v, v + n, [&](int a, int b) { return sc[a] < sc[b]; });
    std::vector<int> sizes; int cum = 0, last = 0;
    for (int i = 0; i < n; i++) { cum += sc[v[i]]; if (cum == (i + 1) * i / 2) { sizes.push_back(i + 1 - last); last = i + 1; } }
    std::reverse(sizes.begin(), sizes.end()); std::string s;
    for (size_t i = 0; i < sizes.size(); i++) { if (i) s += ">"; s += std::to_string(sizes[i]); }
    (void)full; return s;
}
int main(int argc, char **argv) {
    n = atoi(argv[1]); bool details = argc > 2 && !strcmp(argv[2], "details"); uint32_t full = (1u << n) - 1;
    FILE *dump = nullptr; for (int a = 2; a + 1 < argc; a++) if (!strcmp(argv[a], "--dump")) dump = fopen(argv[a + 1], "w");
    std::vector<Tour> reps; std::unordered_map<uint64_t, int> id;
    char line[256];
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] != '0' && line[0] != '1') continue;
        Tour T; memset(&T, 0, sizeof T); int k = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (line[k] == '1') T.out[i] |= 1u << j; else T.out[j] |= 1u << i; k++; }
        id[canonCode(T)] = (int)reps.size(); reps.push_back(T);
    }
    int C = (int)reps.size();
    std::vector<std::vector<int>> adjAny(C), adjBlack(C);
    std::vector<int> L2(C), lam(C), kap(C), d0(C), fas(C), R2(C), aut(C), H(C); std::vector<char> sc(C), strong(C);
    std::vector<char> tTT(C), tNS(C), tSING(C), tCORE(C), tREG(C);
    std::vector<std::string> scores(C), comps(C);
    for (int c = 0; c < C; c++) {
        const Tour &T = reps[c];
        double g; uint64_t code = canonCode(T, &g); aut[c] = (int)(g + 0.5);
        sc[c] = canonCode(converse(T)) == code;
        int best = 99; for (uint32_t X = 1; X < full; X++) { int f = 0; for (int x = 0; x < n; x++) if (X >> x & 1) f += __builtin_popcount(T.out[x] & ~X & full); best = std::min(best, f); }
        lam[c] = best; strong[c] = best > 0;
        // vertex connectivity: least |S| with T - S not strong (n - 1 if none)
        int kk = 0; if (best > 0) { kk = n - 1; for (uint32_t S = 1; S < full; S++) { int sz = __builtin_popcount(S); if (sz >= kk || sz > n - 2) continue; if (!strongMask(T, full & ~S)) kk = sz; } }
        kap[c] = kk;
        int s[16]; int dd = 99; bool sing = false, core = false; int r2 = 0;
        for (int v = 0; v < n; v++) { s[v] = __builtin_popcount(T.out[v]); dd = std::min(dd, std::min(s[v], n - 1 - s[v])); r2 += abs(2 * s[v] - (n - 1));
            if (s[v] == 0 || s[v] == n - 1) { sing = true; if (strongMask(T, full & ~(1u << v))) core = true; } }
        std::sort(s, s + n); int l2 = 0; for (int i = 0; i < n; i++) l2 += abs(s[i] - i);
        L2[c] = l2; d0[c] = dd; R2[c] = r2;
        std::string ss; for (int i = 0; i < n; i++) { ss += std::to_string(s[i]); }
        scores[c] = ss; comps[c] = compStr(T);
        tTT[c] = (l2 == 0); tNS[c] = best == 0; tSING[c] = sing; tCORE[c] = core;
        tREG[c] = (n % 2 == 1) ? (r2 == 0) : (r2 == n);
        fas[c] = fasDP(T);
        int h = 0;
        std::vector<std::pair<int,int>> nb;
        paths(T, [&](const int *Q) { h++;
            Tour Rr = kmove(T, Q); int t = id.at(canonCode(Rr));
            nb.push_back({t, isBlack(T, Q) ? 1 : 0});
        });
        H[c] = h;
        std::sort(nb.begin(), nb.end()); nb.erase(std::unique(nb.begin(), nb.end()), nb.end());
        for (auto &p : nb) { if (adjAny[c].empty() || adjAny[c].back() != p.first) adjAny[c].push_back(p.first); if (p.second) adjBlack[c].push_back(p.first); }
    }
    // symmetry of the move graph: check undirectedness
    long asym = 0; for (int c = 0; c < C; c++) for (int t : adjAny[c]) if (!std::binary_search(adjAny[t].begin(), adjAny[t].end(), c)) asym++;
    long asymB = 0; for (int c = 0; c < C; c++) for (int t : adjBlack[c]) if (!std::binary_search(adjBlack[t].begin(), adjBlack[t].end(), c)) asymB++;
    auto bfs = [&](const std::vector<char> &tgt, const std::vector<std::vector<int>> &adj) {
        std::vector<int> d(C, -1); std::deque<int> q; for (int c = 0; c < C; c++) if (tgt[c]) { d[c] = 0; q.push_back(c); }
        while (!q.empty()) { int c = q.front(); q.pop_front(); for (int t : adj[c]) if (d[t] < 0) { d[t] = d[c] + 1; q.push_back(t); } }
        return d;
    };
    printf("n=%d classes=%d asymmetric-edges any=%ld black=%ld\n", n, C, asym, asymB);
    if (dump) {   // one line per class: id code scores comps aut sc lambda delta0 L2 fas R2 H | any-neighbours | black-neighbours
        for (int c = 0; c < C; c++) {
            fprintf(dump, "%d %llu %s %s %d %d %d %d %d %d %d %d |", c, (unsigned long long)canonCode(reps[c]), scores[c].c_str(), comps[c].c_str(), aut[c], sc[c], lam[c], d0[c], L2[c], fas[c], R2[c], H[c]);
            for (int t : adjAny[c]) fprintf(dump, " %d", t);
            fprintf(dump, " |"); for (int t : adjBlack[c]) fprintf(dump, " %d", t);
            fprintf(dump, " | arcs ");
            for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) fputc((reps[c].out[i] >> j & 1) ? '1' : '0', dump);
            fprintf(dump, "\n");
        }
        fclose(dump);
    }
    struct Job { const char *name; std::vector<char> *tgt; const char *pname; std::vector<int> pot; int step; };
    // potentials as integers with unit step: L (=L2/2), ceil(fas/(n-1)), lambda, delta0, R (=R2/2 for odd n, (R2-n)/2 for even n)
    std::vector<int> Lp(C), Fp(C), Rp(C);
    for (int c = 0; c < C; c++) { Lp[c] = L2[c] / 2; Fp[c] = (fas[c] + n - 2) / (n - 1); Rp[c] = (n % 2) ? R2[c] / 4 : (R2[c] - n + 3) / 4; }
    std::vector<Job> jobs = { {"TT", &tTT, "L", Lp, 1}, {"TT", &tTT, "ceil(fas/(n-1))", Fp, 1}, {"NS", &tNS, "lambda", lam, 1},
                              {"SING", &tSING, "delta0", d0, 1}, {"CORE", &tCORE, "delta0", d0, 1}, {"REG", &tREG, "R", Rp, 1} };
    for (auto &J : jobs) {
        for (int col = 0; col < 2; col++) {
            auto &adj = col ? adjBlack : adjAny;
            auto d = bfs(*J.tgt, adj);
            int unreach = 0, below = 0, fail = 0, lmin = 0, lminSym = 0, failSym = 0, colour = 0, colourSC = 0;
            std::map<std::pair<int,int>, int> lminConn;
            std::map<std::pair<int,int>, int> hist;
            std::vector<int> lminList, failList;
            for (int c = 0; c < C; c++) {
                if (d[c] < 0) { unreach++; continue; }
                if (d[c] < J.pot[c]) below++;
                hist[{J.pot[c], d[c] - J.pot[c]}]++;
                if (d[c] > J.pot[c]) { fail++; failList.push_back(c); if (aut[c] > 1 || sc[c]) failSym++; }
                if (!(*J.tgt)[c]) {
                    bool desc = false; for (int t : adj[c]) if (J.pot[t] == J.pot[c] - 1) { desc = true; break; }
                    if (!desc) { lmin++; lminList.push_back(c); if (aut[c] > 1 || sc[c]) lminSym++; lminConn[{lam[c], kap[c]}]++;
                        if (col) { bool anyDesc = false; for (int t : adjAny[c]) if (J.pot[t] == J.pot[c] - 1) { anyDesc = true; break; }
                                   if (anyDesc) { colour++; colourSC += sc[c]; } } }
                }
            }
            printf("target %-4s %-5s potential %-16s unreachable %d  d<bound %d  failures %d (symmetric %d)  local-minima %d (symmetric %d)\n   (bound, d-bound): ",
                   J.name, col ? "black" : "any", J.pname, unreach, below, fail, failSym, lmin, lminSym);
            for (auto &p : hist) printf("(%d,+%d):%d ", p.first.first, p.first.second, p.second);
            printf("\n   local minima by (lambda, kappa): ");
            for (auto &p : lminConn) printf("(%d,%d):%d ", p.first.first, p.first.second, p.second);
            if (col) printf("\n   colour obstructions (descending moves exist, all blue): %d, self-converse %d", colour, colourSC);
            printf("\n");
            if (details) {
                for (int c : lminList) printf("   LOCALMIN class %d scores %s comps %s |Aut| %d SC %d lambda %d kappa %d delta0 %d L %d fas %d H %d pot %d dist %d\n", c, scores[c].c_str(), comps[c].c_str(), aut[c], sc[c], lam[c], kap[c], d0[c], Lp[c], fas[c], H[c], J.pot[c], d[c]);
                int shown = 0;
                for (int c : failList) { if (shown++ >= 60) break; printf("   FAIL class %d scores %s comps %s |Aut| %d SC %d lambda %d delta0 %d L %d fas %d H %d pot %d dist %d\n", c, scores[c].c_str(), comps[c].c_str(), aut[c], sc[c], lam[c], d0[c], Lp[c], fas[c], H[c], J.pot[c], d[c]); }
            }
        }
    }
    return 0;
}
