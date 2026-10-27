// family3.cpp -- sanity search in the regime 2 <= lambda < delta with t >= 2 at odd order n = 2a+1:
// bottom S (a vertices), sigma-fixed middle vertex c, top = converse copy, all inter-layer pairs down except the
// sigma-symmetric forward arcs:  alpha = 2: x1->c, x2->c, c->s(x1), c->s(x2);  alpha = 1: x1->c, c->s(x1), x->s(x);
// alpha = 0: two sigma-symmetric bottom->top arcs (c is a skipped layer).  Keeps T with lambda >= 2 and
// lambda < delta and reports KEY.  Usage: gentourng -cq a | ./family3 a
#include "tour.hpp"
#include <map>
int main(int argc, char **argv) {
    int a = atoi(argv[1]); int n = 2 * a + 1; char line[512]; Tour S;
    long long built = 0, regime = 0, keyFail = 0; map<string, long long> st;
    while (fgets(line, sizeof line, stdin)) {
        if (!parseLine(line, a, S)) continue;
        int low = 0; for (int x = 0; x < a; x++) if (S.s[x] < 3) low += 3 - S.s[x];
        if (low > 2) continue;
        auto top = [&](int i) { return n - 1 - i; }; int c = a;
        vector<vector<pair<int,int>>> configs;   // list of forward arcs (u,v)
        for (int x1 = 0; x1 < a; x1++) for (int x2 = x1 + 1; x2 < a; x2++) configs.push_back({{x1, c}, {x2, c}, {c, top(x1)}, {c, top(x2)}});
        for (int x1 = 0; x1 < a; x1++) for (int x = 0; x < a; x++) configs.push_back({{x1, c}, {c, top(x1)}, {x, top(x)}});
        for (int x = 0; x < a; x++) for (int z = x + 1; z < a; z++) { configs.push_back({{x, top(x)}, {z, top(z)}}); configs.push_back({{x, top(z)}, {z, top(x)}}); }
        for (auto &F : configs) {
            Tour T; T.n = n; memset(T.out, 0, sizeof T.out);
            auto addArc = [&](int u, int v) { T.out[u] |= 1u << v; T.out[v] &= ~(1u << u); };
            for (int u = 0; u < n; u++) for (int v = 0; v < u; v++) addArc(u, v);
            for (int i = 0; i < a; i++) for (int j = 0; j < a; j++) if (i != j && S.arc(i, j)) { addArc(i, j); addArc(top(j), top(i)); }
            for (auto &f : F) addArc(f.first, f.second);
            T.finish(); built++;
            int delta = min(T.smin(), n - 1 - T.smax());
            vector<uint32_t> cuts; int lam = minCuts(T, cuts);
            if (!(lam >= 2 && lam < delta) || !strongMask(T, (1u << n) - 1)) continue;
            regime++;
            LemmaResult R = analyse(T, true);
            char buf[120]; snprintf(buf, sizeof buf, "lambda=%d delta=%d ncuts=%d key=%d strongForm=%d", lam, delta, (int)cuts.size(), R.keyLemma, R.strongForm);
            st[buf]++; if (!R.keyLemma) keyFail++;
        }
    }
    printf("a=%d n=%d built=%lld in regime=%lld KEY fails=%lld\n", a, n, built, regime, keyFail);
    for (auto &p : st) printf("  %8lld %s\n", p.second, p.first.c_str());
    return 0;
}
