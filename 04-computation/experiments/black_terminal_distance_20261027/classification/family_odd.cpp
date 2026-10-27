// family_odd.cpp -- odd-order shape G(S,x0): bottom S (a vertices), a sigma-fixed middle vertex c, top = converse copy;
// x0 -> c -> sigma(x0) are the only upward arcs; every other pair between layers points down.  Every odd-order
// exception has this shape (S = T[X_{t/2}], c = c_{t/2}).  Counts KEY failures by number of minimum cuts.
// Usage: gentourng -cq a | ./family_odd a
#include "tour.hpp"
#include <map>
int main(int argc, char **argv) {
    int a = atoi(argv[1]); int n = 2 * a + 1; char line[512]; Tour S;
    long long built = 0, fails = 0; map<string, long long> st; vector<string> ex;
    while (fgets(line, sizeof line, stdin)) {
        if (!parseLine(line, a, S)) continue;
        for (int x0 = 0; x0 < a; x0++) {
            Tour T; T.n = n; memset(T.out, 0, sizeof T.out);
            int c = a; auto top = [&](int i) { return n - 1 - i; };
            auto addArc = [&](int u, int v) { T.out[u] |= 1u << v; T.out[v] &= ~(1u << u); };
            for (int u = 0; u < n; u++) for (int v = 0; v < u; v++) addArc(u, v);   // default: down
            for (int i = 0; i < a; i++) for (int j = 0; j < a; j++) if (i != j && S.arc(i, j)) { addArc(i, j); addArc(top(j), top(i)); }
            addArc(x0, c); addArc(c, top(x0));
            T.finish();
            if (!strongMask(T, (1u << n) - 1)) continue;
            built++;
            LemmaResult R = analyse(T, true);
            char buf[160]; snprintf(buf, sizeof buf, "lambda=%d delta=%d ncuts=%d key=%d", R.lambda, min(T.smin(), n - 1 - T.smax()), R.ncuts, R.keyLemma);
            st[buf]++;
            if (!R.keyLemma) { fails++; ex.push_back(encode(T) + " S=" + string(line, strlen(line) - 1) + " x0=" + to_string(x0) + " " + buf); }
        }
    }
    printf("a=%d n=%d built=%lld KEY fails=%lld\n", a, n, built, fails);
    for (auto &p : st) if (p.first.find("key=0") != string::npos) printf("  %8lld %s\n", p.second, p.first.c_str());
    for (auto &e : ex) printf("  EX %s\n", e.c_str());
    return 0;
}
