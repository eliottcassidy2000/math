// family1.cpp -- family F(S,x0): X = copy of S, Y = sigma(X) carrying S^op, single arc x0 -> sigma(x0), all other
// X-Y pairs Y -> X.  T is self-converse (sigma), strong if S strong, lambda = 1.
// Run: gentourng -cq a | ./family1 a
#include "tour.hpp"
#include <map>
int main(int argc, char **argv) {
    int a = argc > 1 ? atoi(argv[1]) : 5;
    int n = 2 * a;
    char line[512]; Tour S;
    long long tot = 0, keyFail = 0, strongFail = 0, openFail = 0, nonStrong = 0;
    map<string, long long> st;
    vector<string> ex;
    while (fgets(line, sizeof line, stdin)) {
        if (!parseLine(line, a, S)) continue;
        for (int x0 = 0; x0 < a; x0++) {
            Tour T; T.n = n; memset(T.out, 0, sizeof T.out);
            for (int i = 0; i < a; i++) for (int j = 0; j < a; j++) if (i != j && S.arc(i, j)) { T.out[i] |= 1u << j; T.out[a + j] |= 1u << (a + i); }
            for (int i = 0; i < a; i++) for (int j = 0; j < a; j++) { if (i == x0 && j == x0) T.out[i] |= 1u << (a + j); else T.out[a + j] |= 1u << i; }
            T.finish();
            if (!strongMask(T, (1u << n) - 1)) { nonStrong++; continue; }
            tot++;
            LemmaResult R = analyse(T, true);
            char buf[200];
            snprintf(buf, sizeof buf, "lambda=%d delta=%d ncuts=%d strongForm=%d openScore=%d key=%d", R.lambda, T.smin(), R.ncuts, R.strongForm, R.openScore, R.keyLemma);
            st[buf]++;
            if (!R.keyLemma) keyFail++;
            if (!R.strongForm) strongFail++;
            if (!R.openScore) openFail++;
            if ((!R.strongForm || !R.keyLemma) && ex.size() < 20) ex.push_back(encode(T) + "  S=" + string(line, strlen(line) - 1) + " x0=" + to_string(x0) + "  " + buf);
        }
    }
    printf("a=%d n=%d  strong T built=%lld (nonstrong skipped %lld)  strongForm fails=%lld  openScore fails=%lld  KEY fails=%lld\n", a, n, tot, nonStrong, strongFail, openFail, keyFail);
    for (auto &p : st) printf("  %8lld %s\n", p.second, p.first.c_str());
    for (auto &e : ex) printf("  EX %s\n", e.c_str());
    return 0;
}
