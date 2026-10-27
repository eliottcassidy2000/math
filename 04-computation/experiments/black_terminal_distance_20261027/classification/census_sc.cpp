// census_sc.cpp -- read gentourng -c output; keep the nonregular strong SELF-CONVERSE classes and test, exactly:
//   strongForm : a minimum out-cut X and a Ham-cycle arc u->w, u notin X, w in X, s(u)+s(w) != n-1
//   openScore  : a minimum out-cut X and a Ham path a~>b, a in X, b notin X, s(a)+s(b) != n-1
//   key        : a minimum out-cut X and a BLACK Ham path a~>b, a in X, b notin X        (the KEY LEMMA)
// Non-self-converse classes satisfy KEY trivially (every Ham path is black) and are only counted.
// Usage: gentourng -cq n | ./census_sc n [dump]   (dump: print every self-converse class with its flags)
#include "tour.hpp"
#include <map>
int main(int argc, char **argv) {
    int n = argc > 1 ? atoi(argv[1]) : 8;
    bool dump = argc > 2 && string(argv[2]) == "dump";
    char line[512]; Tour T;
    long long classes = 0, reg = 0, nsc = 0, sc = 0, keyFail = 0, strongFail = 0, openFail = 0;
    map<string, long long> st; vector<string> ex;
    while (fgets(line, sizeof line, stdin)) {
        if (!parseLine(line, n, T)) continue;
        classes++;
        if (T.regular()) { reg++; continue; }
        // quick score-sequence filter for self-converseness
        int sc_[20]; for (int v = 0; v < n; v++) sc_[v] = T.s[v]; sort(sc_, sc_ + n);
        bool symm = true; for (int i = 0; i < n; i++) if (sc_[i] + sc_[n - 1 - i] != n - 1) { symm = false; break; }
        if (!symm || !selfConverse(T)) { nsc++; continue; }
        sc++;
        LemmaResult R = analyse(T, true);
        int delta = T.smin();
        char buf[200];
        snprintf(buf, sizeof buf, "lambda=%d delta=%d lam<delta=%d ncuts=%s strongForm=%d openScore=%d key=%d", R.lambda, delta, R.lambda < delta,
                 R.ncuts == 1 ? "1" : (R.ncuts <= 3 ? "2-3" : ">3"), R.strongForm, R.openScore, R.keyLemma);
        st[buf]++;
        if (!R.keyLemma) keyFail++;
        if (!R.strongForm) strongFail++;
        if (!R.openScore) openFail++;
        string code(line, strlen(line) - 1);
        if (dump) printf("SC %s lambda=%d delta=%d ncuts=%d strongForm=%d openScore=%d key=%d\n", code.c_str(), R.lambda, delta, R.ncuts, R.strongForm, R.openScore, R.keyLemma);
        if ((!R.keyLemma || !R.strongForm) && ex.size() < 200) {
            string ss; for (int i = 0; i < n; i++) ss += to_string(sc_[i]) + (i + 1 < n ? "," : "");
            ex.push_back(string(R.keyLemma ? "strongFormOnlyFail " : "KEYFAIL ") + code + " scores " + ss + " " + buf);
        }
    }
    printf("n=%d strong classes=%lld regular=%lld nonregular-non-self-converse=%lld nonregular-self-converse=%lld\n", n, classes, reg, nsc, sc);
    printf("  among nonregular self-converse: strongForm fails=%lld  openScore fails=%lld  KEY LEMMA fails=%lld\n", strongFail, openFail, keyFail);
    for (auto &p : st) printf("  %8lld  %s\n", p.second, p.first.c_str());
    for (auto &e : ex) printf("  %s\n", e.c_str());
    return 0;
}
