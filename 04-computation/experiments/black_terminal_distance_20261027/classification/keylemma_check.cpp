// keylemma_check.cpp -- exhaustive checks of the key claims of PROOF.md on every strong class of order n
// (read from `gentourng -cq n`; intended for n <= 8).
//
// For every strong class T it verifies:
//  [C1] minimum out-cuts never cross (X,Z min cuts => nested, disjoint or co-disjoint);  when lambda < delta
//       (delta = min(min outdegree, min indegree)) they are nested (a chain).
//  [C2] every anti-automorphism tau of the form r_Q (Q a blue Ham path) maps each minimum cut X to a set with
//       V \ tau(X) again a minimum cut (fact F3).
//  [C3] KEY classification: T has a black descent candidate unless T is R5 (no exception of the family E exists
//       below order 10).
//  [C4] lambda = delta: there is a "singleton-type" black descent Q (q1 has outdegree lambda, or qn has indegree
//       lambda) such that K_Q(T) again satisfies lambda = delta  (propagation step of the distance theorem).
//  [C5] self-converse, nonregular, lambda = delta, and some min-score v with T - v nonstrong: every Hamiltonian
//       cycle cut starting at v is black with s(last) < Delta  (cycles.md Section 7, used in Part II).
//  [C6] self-converse, nonregular, lambda = delta, every min-score deletion strong (case C2 of Part II):
//       some Ham path starting at a min-score vertex or ending at a max-score vertex is black.
//  [C7] cut identity: f_{K_Q T}(V\X) = f_T(X) - [q1 in X] + [qn in X] for every Ham path Q and every X.
// Build: g++ -O2 -std=c++17 -o keylemma_check keylemma_check.cpp ; run: gentourng -cq n | ./keylemma_check n
#include "tour.hpp"
#include <map>
#include <set>
#include <array>
static Tour kmove(const Tour &T, const int *Q) {
    Tour R; R.n = T.n; uint32_t full = (1u << T.n) - 1;
    for (int u = 0; u < T.n; u++) R.out[u] = (~T.out[u]) & full & ~(1u << u);
    for (int i = 0; i + 1 < T.n; i++) { int u = Q[i], v = Q[i + 1]; R.out[u] |= 1u << v; R.out[v] &= ~(1u << u); }
    R.finish(); return R;
}
static int deltaOf(const Tour &T) { return min(T.smin(), T.n - 1 - T.smax()); }
int main(int argc, char **argv) {
    int n = argc > 1 ? atoi(argv[1]) : 7;
    bool cutid = !(argc > 2 && string(argv[2]) == "nocutid");
    char line[256]; Tour T;
    long long classes = 0, fail[8] = {0}, tested[8] = {0};
    map<string, long long> cases; vector<string> bad;
    uint32_t full = (1u << n) - 1;
    while (fgets(line, sizeof line, stdin)) {
        if (!parseLine(line, n, T)) continue;
        classes++;
        string code = encode(T);
        bool reg = T.regular(), sc = selfConverse(T);
        int delta = deltaOf(T), Dm = T.smax(), dm = T.smin();
        vector<uint32_t> cuts; int lam = minCuts(T, cuts);
        // [C1]
        tested[1]++;
        bool c1ok = true;
        for (size_t i = 0; i < cuts.size(); i++) for (size_t j = i + 1; j < cuts.size(); j++) {
            uint32_t X = cuts[i], Z = cuts[j];
            bool nested = ((X & Z) == X) || ((X & Z) == Z), disj = (X & Z) == 0, codisj = (X | Z) == full;
            if (!(nested || disj || codisj)) c1ok = false;
            if (lam < delta && !nested) c1ok = false;
        }
        if (!c1ok) { fail[1]++; if (bad.size() < 20) bad.push_back("C1 " + code); }
        // enumerate all Ham paths once
        vector<array<int, 20>> paths; vector<char> blue;
        for (int s0 = 0; s0 < n; s0++) for (int e = 0; e < n; e++) if (e != s0)
            pathsAB(T, s0, e, [&](const int *Q) { array<int, 20> a; for (int i = 0; i < n; i++) a[i] = Q[i]; paths.push_back(a); return false; });
        for (auto &P : paths) blue.push_back(isBlue(T, P.data()));
        auto isCand = [&](const int *Q) { for (uint32_t X : cuts) if ((X >> Q[0] & 1) && !(X >> Q[n - 1] & 1)) return true; return false; };
        // [C2]
        if (sc) {
            tested[2]++; bool ok = true;
            set<uint32_t> cutset(cuts.begin(), cuts.end());
            for (size_t k = 0; k < paths.size() && ok; k++) if (blue[k]) {
                int r[20]; for (int i = 0; i < n; i++) r[paths[k][i]] = paths[k][n - 1 - i];
                for (uint32_t X : cuts) { uint32_t img = 0; for (int v = 0; v < n; v++) if (X >> v & 1) img |= 1u << r[v]; if (!cutset.count(full & ~img)) { ok = false; break; } }
            }
            if (!ok) { fail[2]++; if (bad.size() < 20) bad.push_back("C2 " + code); }
        }
        // [C3]
        tested[3]++;
        bool key = false;
        for (size_t k = 0; k < paths.size() && !key; k++) if (!blue[k] && isCand(paths[k].data())) key = true;
        bool isR5 = (n == 5 && reg);
        if (key == isR5) { fail[3]++; if (bad.size() < 20) bad.push_back("C3 " + code); }
        // [C4]
        if (lam == delta && !isR5) {
            tested[4]++; bool found = false;
            for (size_t k = 0; k < paths.size() && !found; k++) {
                if (blue[k]) continue;
                const int *Q = paths[k].data();
                bool st = (T.s[Q[0]] == lam && (uint32_t)lam == (uint32_t)lam) && T.s[Q[0]] == lam; // q1 outdegree lambda -> {q1} min cut
                bool en = (n - 1 - T.s[Q[n - 1]]) == lam;                                      // qn indegree lambda -> V\{qn} min cut
                if (!st && !en) continue;
                Tour R = kmove(T, Q); vector<uint32_t> c2; int l2 = minCuts(R, c2);
                if (l2 == lam - 1 && l2 == deltaOf(R)) found = true;
            }
            if (!found) { fail[4]++; if (bad.size() < 20) bad.push_back("C4 " + code); }
        }
        // case label + [C5]/[C6]
        string label;
        if (!sc) label = "non-self-converse";
        else if (reg) label = "regular";
        else if (lam < delta) label = "self-converse lambda<delta";
        else {
            int vNon = -1; for (int v = 0; v < n; v++) if (T.s[v] == dm && !strongMask(T, full & ~(1u << v))) { vNon = v; break; }
            if (vNon >= 0) {
                label = "self-converse lambda=delta, a min deletion nonstrong";
                tested[5]++; bool ok = true; int cnt = 0;
                for (size_t k = 0; k < paths.size(); k++) {
                    const int *Q = paths[k].data();
                    if (Q[0] != vNon || !T.arc(Q[n - 1], Q[0])) continue;
                    cnt++;
                    if (blue[k] || T.s[Q[n - 1]] >= Dm) ok = false;
                }
                if (!ok || cnt == 0) { fail[5]++; if (bad.size() < 20) bad.push_back("C5 " + code); }
            } else {
                label = "self-converse lambda=delta, all min deletions strong";
                tested[6]++; bool ok = false;
                for (size_t k = 0; k < paths.size() && !ok; k++) {
                    const int *Q = paths[k].data();
                    if (!blue[k] && (T.s[Q[0]] == dm || T.s[Q[n - 1]] == Dm)) ok = true;
                }
                if (!ok) { fail[6]++; if (bad.size() < 20) bad.push_back("C6 " + code); }
            }
        }
        cases[label + (key ? " [black descent]" : " [NO black descent]")]++;
        // [C7]
        if (cutid) {
            tested[7]++; bool ok = true;
            for (size_t k = 0; k < paths.size() && ok; k++) {
                const int *Q = paths[k].data(); Tour R = kmove(T, Q);
                for (uint32_t X = 1; X < full && ok; X++) if (fOut(R, full & ~X) != fOut(T, X) - (int)(X >> Q[0] & 1) + (int)(X >> Q[n - 1] & 1)) ok = false;
            }
            if (!ok) { fail[7]++; if (bad.size() < 20) bad.push_back("C7 " + code); }
        }
    }
    printf("n=%d strong classes=%lld\n", n, classes);
    const char *nm[8] = {"", "C1 min cuts non-crossing / chain if lambda<delta", "C2 blue reversals permute min cuts (X -> V\\r(X))",
                         "C3 black descent exists iff T != R5", "C4 lambda=delta: singleton-type black descent keeping lambda=delta",
                         "C5 Section-7 case: all cycle cuts at v black, s(last)<Delta", "C6 case C2: black path from a min / into a max vertex",
                         "C7 cut identity on every (Q, X)"};
    for (int c = 1; c <= 7; c++) printf("  [%s] tested %lld classes, failures %lld\n", nm[c], tested[c], fail[c]);
    printf("  case distribution:\n");
    for (auto &p : cases) printf("    %8lld  %s\n", p.second, p.first.c_str());
    for (auto &b : bad) printf("  BAD %s\n", b.c_str());
    return 0;
}
