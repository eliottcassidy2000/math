// e_check.cpp -- verifies, for members E_t(S, x0) of the exceptional family (bottom layers read from
// bottom_layers output "S code x0=..."), the claims of PROOF.md Part VII/VIII:
//  (e1) strong, self-converse, lambda = 1, delta = 2;
//  (e2) the minimum out-cuts are exactly the chain X_k = L0 + {c_1..c_{k-1}}, k = 1..t;
//  (e3) no black descent candidate;
//  (e4) the explicit word Q = (r, r', r'', x0, c_1..c_{t-1}, P_top, P_rest) (r r' r'' the last C3 block of S - x0 in
//       cycle order, P_top a Ham path of the top layer from sigma(x0), P_rest a Ham path of the remaining bottom
//       blocks) is a black Ham path, K_Q(E) has lambda = 1 = min indegree, and K_Q(E) has a black descent: so the
//       black distance of E is exactly 2;
//  (e5) at least three vertices of outdegree 2 (used in the avoidance lemma).
// Usage: grep '^S ' bl_full_a.txt | ./e_check a tmax_n
#include "tour.hpp"
#include <functional>
static Tour kmove(const Tour &T, const int *Q) {
    Tour R; R.n = T.n; uint32_t full = (1u << T.n) - 1;
    for (int u = 0; u < T.n; u++) R.out[u] = (~T.out[u]) & full & ~(1u << u);
    for (int i = 0; i + 1 < T.n; i++) { int u = Q[i], v = Q[i + 1]; R.out[u] |= 1u << v; R.out[v] &= ~(1u << u); }
    R.finish(); return R;
}
// Ham path of T[mask] starting at s (s = -1: any start); returns false if none
static bool hamPathInMask(const Tour &T, uint32_t mask, int s, vector<int> &out) {
    int k = __builtin_popcount(mask); vector<int> P; uint32_t used = 0;
    function<bool(int)> dfs = [&](int v) -> bool {
        P.push_back(v); used |= 1u << v;
        if ((int)P.size() == k) return true;
        uint32_t c = T.out[v] & mask & ~used;
        while (c) { int w = __builtin_ctz(c); c &= c - 1; if (dfs(w)) return true; }
        P.pop_back(); used &= ~(1u << v); return false;
    };
    for (int st = 0; st < T.n; st++) if ((mask >> st & 1) && (s < 0 || st == s)) { P.clear(); used = 0; if (dfs(st)) { out = P; return true; } }
    return false;
}
int main(int argc, char **argv) {
    int a = atoi(argv[1]); int nmax = atoi(argv[2]);
    char line[512]; long long members = 0, okAll = 0;
    while (fgets(line, sizeof line, stdin)) {
        char code[256], x0s[64]; if (sscanf(line, "S %255s %63s", code, x0s) != 2) continue;
        int x0 = atoi(x0s + 3);
        Tour S; if (!parseLine(code, a, S)) continue;
        for (int t = 1; 2 * a + t - 1 <= nmax; t++) {
            int n = 2 * a + t - 1; members++;
            Tour T; T.n = n; memset(T.out, 0, sizeof T.out);
            auto top = [&](int i) { return n - 1 - i; }; auto mid = [&](int k) { return a + k - 1; };
            auto addArc = [&](int u, int v) { T.out[u] |= 1u << v; T.out[v] &= ~(1u << u); };
            for (int u = 0; u < n; u++) for (int v = 0; v < u; v++) addArc(u, v);
            for (int i = 0; i < a; i++) for (int j = 0; j < a; j++) if (i != j && S.arc(i, j)) { addArc(i, j); addArc(top(j), top(i)); }
            if (t == 1) addArc(x0, top(x0)); else { addArc(x0, mid(1)); for (int k = 1; k + 1 <= t - 1; k++) addArc(mid(k), mid(k + 1)); addArc(mid(t - 1), top(x0)); }
            T.finish();
            uint32_t full = (1u << n) - 1;
            bool e1 = strongMask(T, full) && selfConverse(T);
            vector<uint32_t> cuts; int lam = minCuts(T, cuts);
            int delta = min(T.smin(), n - 1 - T.smax());
            e1 = e1 && lam == 1 && delta == 2;
            // (e2)
            bool e2 = (int)cuts.size() == t;
            for (int k = 1; k <= t && e2; k++) {
                uint32_t X = (1u << a) - 1; for (int j = 1; j <= k - 1; j++) X |= 1u << mid(j);
                if (find(cuts.begin(), cuts.end(), X) == cuts.end()) e2 = false;
            }
            // (e3)
            LemmaResult R = analyse(T, true); bool e3 = !R.keyLemma;
            // (e4) last block of S - x0: the strong component with no out-arcs inside S - x0
            uint32_t rest = ((1u << a) - 1) & ~(1u << x0); uint32_t last = 0;
            for (int v = 0; v < a; v++) if (rest >> v & 1) {
                // component of v within rest: vertices mutually reachable
                uint32_t comp = 0;
                for (int w = 0; w < a; w++) if (rest >> w & 1) {
                    auto reach = [&](int s0, int g) { uint32_t r = 1u << s0, fr = r; while (fr) { int x = __builtin_ctz(fr); fr &= fr - 1; uint32_t nb = S.out[x] & rest & ~r; r |= nb; fr |= nb; } return (bool)(r >> g & 1); };
                    if (reach(v, w) && reach(w, v)) comp |= 1u << w;
                }
                bool sink = true; for (int w = 0; w < a; w++) if (comp >> w & 1) if (S.out[w] & rest & ~comp) sink = false;
                if (sink) { last = comp; break; }
            }
            bool e4 = false;
            if (__builtin_popcount(last) == 3) {
                int r0 = __builtin_ctz(last); vector<int> cyc{r0};
                for (int step = 0; step < 2; step++) { uint32_t nb = S.out[cyc.back()] & last; cyc.push_back(__builtin_ctz(nb)); }
                vector<int> Q = cyc; Q.push_back(x0);
                for (int k = 1; k <= t - 1; k++) Q.push_back(mid(k));
                uint32_t topMask = 0; for (int i = 0; i < a; i++) topMask |= 1u << top(i);
                vector<int> Pt, Pr; bool okp = hamPathInMask(T, topMask, top(x0), Pt);
                uint32_t restMask = rest & ~last; okp = okp && hamPathInMask(T, restMask, -1, Pr);
                if (okp) {
                    for (int v : Pt) Q.push_back(v); for (int v : Pr) Q.push_back(v);
                    bool valid = (int)Q.size() == n; for (int i = 0; i + 1 < n && valid; i++) if (!T.arc(Q[i], Q[i + 1])) valid = false;
                    if (valid && !isBlue(T, Q.data())) {
                        Tour K = kmove(T, Q.data()); vector<uint32_t> ck; int lk = minCuts(K, ck);
                        int minIn = n; for (int v = 0; v < n; v++) minIn = min(minIn, n - 1 - K.s[v]);
                        LemmaResult RK = analyse(K, true);
                        e4 = (lk == 1 && minIn == 1 && RK.keyLemma);
                    }
                }
            }
            int deg2 = 0; for (int v = 0; v < n; v++) if (T.s[v] == 2) deg2++;
            bool e5 = deg2 >= 3;
            bool all = e1 && e2 && e3 && e4 && e5; okAll += all;
            printf("S=%s x0=%d t=%d n=%d |Aut|=%lld  e1=%d e2=%d e3=%d e4=%d e5=%d(outdeg-2 vertices %d)%s\n", code, x0, t, n, AutCount().run(T), e1, e2, e3, e4, e5, deg2, all ? "" : "  <-- FAIL");
            fflush(stdout);
        }
    }
    printf("a=%d: members checked=%lld, all claims verified=%lld\n", a, members, okAll);
    return 0;
}
