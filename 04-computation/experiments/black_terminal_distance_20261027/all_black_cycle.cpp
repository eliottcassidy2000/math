// all_black_cycle.cpp -- mac-mini-2026-10-27-potentials
// For every strong class of order n: does T have a Hamiltonian cycle ALL of whose n rotations are black?
// (Then every minimum out-cut X is entered by that cycle, and the rotation starting at the head of an entering
// arc is a black descent.)  Also counts Hamiltonian cycles and those with a blue rotation.
// Build: g++ -O2 -std=c++17 -o all_black_cycle all_black_cycle.cpp ; run: gentourng -cq n | ./all_black_cycle n
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
using namespace std;
static int n;
struct Tour { uint16_t out[16]; };
static inline bool isBlack(const Tour &T, const int *Q) {
    int r[16]; for (int i = 0; i < n; i++) r[Q[i]] = Q[n - 1 - i];
    for (int u = 0; u < n; u++) { uint16_t o = T.out[u]; while (o) { int v = __builtin_ctz(o); o &= o - 1; if (!((T.out[r[v]] >> r[u]) & 1)) return true; } }
    return false;
}
int main(int argc, char **argv) {
    n = argc > 1 ? atoi(argv[1]) : 7;
    char line[256]; long long classes = 0, withAllBlack = 0; vector<string> bad;
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] != '0' && line[0] != '1') continue;
        Tour T; memset(&T, 0, sizeof T); int k = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (line[k] == '1') T.out[i] |= (uint16_t)(1u << j); else T.out[j] |= (uint16_t)(1u << i); k++; }
        classes++;
        // enumerate Hamiltonian cycles through vertex 0 (each cycle once, as a path from 0 closing back to 0)
        int Q[16]; uint16_t cand[16]; int depth = 0; Q[0] = 0; uint16_t used = 1; cand[0] = T.out[0] & ~used;
        bool found = false; long long cyc = 0, cycBlue = 0;
        while (depth >= 0 && !found) {
            if (depth == n - 1) {
                if ((T.out[Q[n - 1]] >> 0) & 1) {
                    cyc++; bool allBlack = true; int R[16];
                    for (int s = 0; s < n && allBlack; s++) { for (int i = 0; i < n; i++) R[i] = Q[(s + i) % n]; if (!isBlack(T, R)) allBlack = false; }
                    if (allBlack) found = true; else cycBlue++;
                }
                used &= (uint16_t)~(1u << Q[depth]); depth--; continue;
            }
            if (!cand[depth]) { used &= (uint16_t)~(1u << Q[depth]); depth--; continue; }
            int v = __builtin_ctz(cand[depth]); cand[depth] &= cand[depth] - 1;
            depth++; Q[depth] = v; used |= (uint16_t)(1u << v); cand[depth] = T.out[v] & ~used;
        }
        if (found) withAllBlack++; else if (bad.size() < 50) bad.push_back(string(line, strlen(line) - 1) + "  HC=" + to_string(cyc));
    }
    printf("n=%d strong classes=%lld  with an all-black Hamiltonian cycle=%lld  without=%lld\n", n, classes, withAllBlack, classes - withAllBlack);
    for (auto &b : bad) printf("    %s\n", b.c_str());
    return 0;
}
