// bottom_layers.cpp -- classify the possible bottom layers (S, x0) of a lambda=1 exception:
//   S strong, Condition A (every nonempty Z subset of V(S)-x0 has at least 2 arcs leaving Z inside S),
//   and Stab_Aut(S)(x0) transitive on the Hamiltonian paths of S ending at x0  (h == |Stab|).
// Reads gentourng -c a output.  Prints every (S, x0) passing, with the strong-component structure of S - x0.
// Usage: gentourng -cq a | ./bottom_layers a
#include "tour.hpp"
#include <map>
static long long countPathsEndingAt(const Tour &S, int x0) {
    int n = S.n; uint32_t full = (1u << n) - 1;
    vector<long long> dp((size_t)1 << n << 5, 0);   // dp[mask][v] = # paths covering mask ending at v
    for (int v = 0; v < n; v++) dp[((size_t)(1u << v) << 5) | v] = 1;
    for (uint32_t mask = 1; mask <= full; mask++) for (int v = 0; v < n; v++) {
        long long c = dp[((size_t)mask << 5) | v]; if (!c) continue;
        uint32_t o = S.out[v] & ~mask;
        while (o) { int w = __builtin_ctz(o); o &= o - 1; dp[((size_t)(mask | (1u << w)) << 5) | w] += c; }
    }
    return dp[((size_t)full << 5) | x0];
}
int main(int argc, char **argv) {
    int a = argc > 1 ? atoi(argv[1]) : 5;
    char line[512]; Tour S; long long pairs = 0, condA = 0, trans = 0;
    map<long long, long long> hHist;
    while (fgets(line, sizeof line, stdin)) {
        if (!parseLine(line, a, S)) continue;
        uint32_t full = (1u << a) - 1;
        for (int x0 = 0; x0 < a; x0++) {
            pairs++;
            bool ok = true;
            uint32_t rest = full & ~(1u << x0);
            // iterate nonempty subsets Z of rest
            for (uint32_t Z = rest; Z && ok; Z = (Z - 1) & rest) if (fOut(S, Z) < 2) ok = false;
            if (!ok) continue;
            condA++;
            long long h = countPathsEndingAt(S, x0);
            long long st = AutCount().run(S, x0);
            if (h != st) continue;
            trans++;
            hHist[h]++;
            // describe S - x0: strong components in order (by reachability), x0 adjacency
            // compute condensation order via repeated extraction of a source component
            vector<vector<int>> comps; uint32_t left = rest;
            while (left) {
                // find a vertex of left with no in-arcs from other components remaining: take component of a vertex that reaches all of left
                for (int v = 0; v < a; v++) if (left >> v & 1) {
                    // reach set of v within left
                    uint32_t r = 1u << v, fr = r;
                    while (fr) { int x = __builtin_ctz(fr); fr &= fr - 1; uint32_t nb = S.out[x] & left & ~r; r |= nb; fr |= nb; }
                    if (r == left) {
                        // component of v = vertices of left reachable from v that also reach v
                        vector<int> comp;
                        for (int w = 0; w < a; w++) if (left >> w & 1) {
                            uint32_t r2 = 1u << w, f2 = r2;
                            while (f2) { int x = __builtin_ctz(f2); f2 &= f2 - 1; uint32_t nb = S.out[x] & left & ~r2; r2 |= nb; f2 |= nb; }
                            if (r2 >> v & 1) comp.push_back(w);
                        }
                        comps.push_back(comp);
                        for (int w : comp) left &= ~(1u << w);
                        break;
                    }
                }
            }
            printf("S %s x0=%d h=%lld |Stab|=%lld  S-x0 blocks:", encode(S).c_str(), x0, h, st);
            for (auto &c : comps) {
                int beats = 0, beaten = 0; for (int w : c) { if (S.arc(x0, w)) beats++; else beaten++; }
                printf(" [size %zu, x0->%d, ->x0 %d]", c.size(), beats, beaten);
            }
            printf("\n");
        }
    }
    printf("a=%d: (S,x0) pairs=%lld  conditionA=%lld  transitive=%lld  h-histogram:", a, pairs, condA, trans);
    for (auto &p : hHist) printf(" %lld:%lld", p.first, p.second);
    printf("\n");
    return 0;
}
