// family_E.cpp -- build E_t(S, x0): bottom layer S (a vertices, code), t-1 middle singletons c_1..c_{t-1},
// top layer = converse copy of S; jump arcs x0 -> c_1 -> ... -> c_{t-1} -> sigma(x0)  (t = 1: x0 -> sigma(x0));
// every other pair between different layers points DOWN (top -> middles -> bottom, c_j -> c_i for j >= i+2).
// For each input line "Scode x0" and each t in [tmin, tmax], prints n, lambda, number of min cuts, delta,
// self-converse, strongForm/openScore/KEY, and (optionally, flag "dist") a two-move black-distance certificate:
// a black move to T' with lambda(T') = 1 that has a black descent candidate.
// Usage: printf "Scode x0\n" | ./family_E a tmin tmax [dist]
#include "tour.hpp"
static Tour buildE(const Tour &S, int x0, int t) {
    int a = S.n, n = 2 * a + (t - 1);
    Tour T; T.n = n; memset(T.out, 0, sizeof T.out);
    // layout: bottom 0..a-1, middles a..a+t-2, top a+t-1 .. a+t-1+a-1 with sigma(i) = n-1-i
    auto top = [&](int i) { return n - 1 - i; };
    auto mid = [&](int k) { return a + k - 1; };   // k = 1..t-1
    auto addArc = [&](int u, int v) { T.out[u] |= 1u << v; T.out[v] &= ~(1u << u); };
    // default: every pair oriented from higher index to lower index (down), then fix
    for (int u = 0; u < n; u++) for (int v = 0; v < u; v++) addArc(u, v);
    for (int i = 0; i < a; i++) for (int j = 0; j < a; j++) if (i != j && S.arc(i, j)) { addArc(i, j); addArc(top(j), top(i)); }
    // jumps
    if (t == 1) addArc(x0, top(x0));
    else {
        addArc(x0, mid(1));
        for (int k = 1; k + 1 <= t - 1; k++) addArc(mid(k), mid(k + 1));
        addArc(mid(t - 1), top(x0));
    }
    T.finish();
    return T;
}
static Tour kmove(const Tour &T, const int *Q) {
    Tour R; R.n = T.n; uint32_t full = (1u << T.n) - 1;
    for (int u = 0; u < T.n; u++) R.out[u] = (~T.out[u]) & full & ~(1u << u);
    for (int i = 0; i + 1 < T.n; i++) { int u = Q[i], v = Q[i + 1]; R.out[u] |= 1u << v; R.out[v] &= ~(1u << u); }
    R.finish(); return R;
}
int main(int argc, char **argv) {
    int a = atoi(argv[1]), tmin = atoi(argv[2]), tmax = atoi(argv[3]);
    bool dist = argc > 4 && string(argv[4]) == "dist";
    char line[512];
    while (fgets(line, sizeof line, stdin)) {
        char code[256]; int x0;
        if (sscanf(line, "%255s %d", code, &x0) != 2) continue;
        Tour S; if (!parseLine(code, a, S)) continue;
        for (int t = tmin; t <= tmax; t++) {
            Tour T = buildE(S, x0, t);
            int n = T.n;
            bool strong = strongMask(T, (1u << n) - 1);
            LemmaResult R = analyse(T, true);
            printf("S=%s x0=%d t=%d n=%d strong=%d lambda=%d ncuts=%d delta+=%d delta-=%d selfconv=%d |Aut|=%lld strongForm=%d openScore=%d KEY=%d\n",
                   code, x0, t, n, strong, R.lambda, R.ncuts, T.smin(), n - 1 - T.smax(), selfConverse(T), AutCount().run(T), R.strongForm, R.openScore, R.keyLemma);
            if (dist && !R.keyLemma) {
                // search a black move Q with lambda(K_Q T) = 1 and K_Q T having a black descent candidate
                bool found = false; long long tried = 0;
                for (int s0 = 0; s0 < n && !found; s0++) {
                    int Q[20]; uint32_t cand[20]; int depth = 0; Q[0] = s0; uint32_t used = 1u << s0; cand[0] = T.out[s0] & ~used;
                    while (depth >= 0 && !found) {
                        if (depth == n - 1) {
                            if (!isBlue(T, Q)) {
                                tried++;
                                Tour R2 = kmove(T, Q);
                                vector<uint32_t> c2; int l2 = minCuts(R2, c2);
                                if (l2 == 1) { LemmaResult RR = analyse(R2, true); if (RR.keyLemma) {
                                    found = true;
                                    printf("   two-move certificate: first black word Q=(");
                                    for (int i = 0; i < n; i++) printf("%d%s", Q[i], i + 1 < n ? "," : "");
                                    printf(") gives lambda=1 tournament %s with a black descent (after %lld black words)\n", encode(R2).c_str(), tried);
                                } }
                            }
                            used &= ~(1u << Q[depth]); depth--; continue;
                        }
                        if (!cand[depth]) { used &= ~(1u << Q[depth]); depth--; continue; }
                        int v = __builtin_ctz(cand[depth]); cand[depth] &= cand[depth] - 1;
                        depth++; Q[depth] = v; used |= 1u << v; cand[depth] = T.out[v] & ~used;
                    }
                }
                if (!found) printf("   NO two-move certificate found\n");
            }
            fflush(stdout);
        }
    }
    return 0;
}
