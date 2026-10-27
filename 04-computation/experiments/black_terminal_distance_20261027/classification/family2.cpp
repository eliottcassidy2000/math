// family2.cpp -- t = 1 configurations with lambda = 2: X = copy of S (a vertices), Y = sigma(X) carrying S^op,
// forward (X -> Y) arcs F = sigma-symmetric pair set {x->sigma(z), z->sigma(x)} (x != z) or {x->sigma(x), z->sigma(z)},
// all other X-Y pairs Y -> X.  Keeps only T with delta(T) > lambda(T) (the chain regime 2 <= lambda < delta) and
// reports KEY / strong form.  Usage: gentourng -cq a | ./family2 a
#include "tour.hpp"
#include <map>
int main(int argc, char **argv) {
    int a = argc > 1 ? atoi(argv[1]) : 7; int n = 2 * a;
    char line[512]; Tour S; long long built = 0, regime = 0, keyFail = 0, strongFail = 0;
    map<string, long long> st; vector<string> ex;
    while (fgets(line, sizeof line, stdin)) {
        if (!parseLine(line, a, S)) continue;
        // quick filter: every x needs s_S(x) + (#forward arcs at x) >= 3, forward arcs at most 2 per x total 2
        int low = 0; for (int x = 0; x < a; x++) if (S.s[x] < 3) low += 3 - S.s[x];
        if (low > 2) continue;
        for (int x = 0; x < a; x++) for (int z = x; z < a; z++) for (int type = 0; type < 2; type++) {
            vector<pair<int,int>> F;
            if (type == 0) { if (x == z) continue; F = {{x, x}, {z, z}}; }
            else { if (x == z) continue; F = {{x, z}, {z, x}}; }
            Tour T; T.n = n; memset(T.out, 0, sizeof T.out);
            auto top = [&](int i) { return n - 1 - i; };
            for (int i = 0; i < a; i++) for (int j = 0; j < a; j++) if (i != j && S.arc(i, j)) { T.out[i] |= 1u << j; T.out[top(j)] |= 1u << top(i); }
            for (int i = 0; i < a; i++) for (int j = 0; j < a; j++) {
                bool fwd = false; for (auto &f : F) if (f.first == i && f.second == j) fwd = true;
                if (fwd) T.out[i] |= 1u << top(j); else T.out[top(j)] |= 1u << i;
            }
            T.finish(); built++;
            int delta = min(T.smin(), n - 1 - T.smax());
            vector<uint32_t> cuts; int lam = minCuts(T, cuts);
            if (!(lam >= 2 && lam < delta)) continue;
            if (!strongMask(T, (1u << n) - 1)) continue;
            regime++;
            LemmaResult R = analyse(T, true);
            char buf[160]; snprintf(buf, sizeof buf, "lambda=%d delta=%d ncuts=%d strongForm=%d key=%d", lam, delta, (int)cuts.size(), R.strongForm, R.keyLemma);
            st[buf]++;
            if (!R.keyLemma) { keyFail++; if (ex.size() < 20) ex.push_back(encode(T) + " " + buf); }
            if (!R.strongForm) strongFail++;
        }
    }
    printf("a=%d n=%d built=%lld in-regime(2<=lambda<delta)=%lld strongForm fails=%lld KEY fails=%lld\n", a, n, built, regime, strongFail, keyFail);
    for (auto &p : st) printf("  %8lld %s\n", p.second, p.first.c_str());
    for (auto &e : ex) printf("  KEYFAIL %s\n", e.c_str());
    return 0;
}
