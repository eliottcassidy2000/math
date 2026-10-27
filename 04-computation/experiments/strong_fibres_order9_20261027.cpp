// strong_fibres_order9_20261027.cpp -- opus-2026-10-27-S63
// Strong fibres at order 9 for small H: reads the labelled strong 9-tournaments with H <= B produced by
// strong_low_spectrum_augment_20261027.cpp (sla_list9.txt, pair codes with bit k = 1 iff i -> j for the k-th pair (i, j),
// i < j, lexicographic), canonizes each by brute force over all 9! relabellings, and counts the isomorphism classes per
// Hamiltonian-path count.  With the complete orders <= 8 (h_spectrum_sweep) this gives the strong fibre F_s(v) exactly
// for every v < f(10) = 125, since every strong tournament with H < 125 has order <= 9.  Test of HYP-9420 (3).
// Build: g++ -O2 -std=c++17 -fopenmp -o sf9 strong_fibres_order9_20261027.cpp ; run: ./sf9 sla_list9.txt
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <numeric>
using namespace std;
typedef unsigned int u32; typedef unsigned long long u64;
static const int n = 9, M = 36;
struct Tour { u32 out[9]; };
static Tour fromCode(u64 c) { Tour T; for (int v = 0; v < n; v++) T.out[v] = 0; int k = 0; for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++, k++) { if (c >> k & 1) T.out[i] |= 1u << j; else T.out[j] |= 1u << i; } return T; }
static inline u64 codeOf(const Tour& T) { u64 c = 0; int k = 0; for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++, k++) if (T.out[i] >> j & 1) c |= 1ULL << k; return c; }
static u64 canon(const Tour& T) { int p[9]; iota(p, p + n, 0); u64 best = ~0ULL;
    do { Tour R; for (int v = 0; v < n; v++) { u32 m = 0, o = T.out[v]; while (o) { int w = __builtin_ctz(o); o &= o - 1; m |= 1u << p[w]; } R.out[p[v]] = m; } u64 c = codeOf(R); if (c < best) best = c; } while (next_permutation(p, p + n));
    return best; }
static u32 countH(const Tour& T) { static thread_local u32 dp[1 << 9][9]; const u32 full = (1u << n) - 1;
    for (u32 S = 1; S <= full; S++) for (int v = 0; v < n; v++) dp[S][v] = 0; for (int v = 0; v < n; v++) dp[1u << v][v] = 1;
    for (u32 S = 1; S <= full; S++) for (int v = 0; v < n; v++) { u32 c = dp[S][v]; if (!c) continue; u32 ext = T.out[v] & ~S; while (ext) { int w = __builtin_ctz(ext); ext &= ext - 1; dp[S | (1u << w)][w] += c; } }
    u32 h = 0; for (int v = 0; v < n; v++) h += dp[full][v]; return h; }
static bool strongT(const Tour& T) { const u32 full = (1u << n) - 1; u32 r = 1, r2; do { r2 = r; u32 t = r; while (t) { int v = __builtin_ctz(t); t &= t - 1; r |= T.out[v]; } } while (r != r2); if (r != full) return false;
    u32 in[9]; for (int v = 0; v < n; v++) { in[v] = 0; for (int w = 0; w < n; w++) if (T.out[w] >> v & 1) in[v] |= 1u << w; } r = 1; do { r2 = r; u32 t = r; while (t) { int v = __builtin_ctz(t); t &= t - 1; r |= in[v]; } } while (r != r2); return r == full; }
int main(int argc, char** argv) {
    vector<u64> codes; { FILE* f = fopen(argc > 1 ? argv[1] : "sla_list9.txt", "r"); if (!f) { perror("list"); return 1; } u64 c; while (fscanf(f, "%llu", &c) == 1) codes.push_back(c); fclose(f); }
    vector<u64> canons(codes.size());
    #pragma omp parallel for schedule(dynamic, 16)
    for (long i = 0; i < (long)codes.size(); i++) canons[i] = canon(fromCode(codes[i]));
    set<u64> cls(canons.begin(), canons.end()); map<u32, int> fib; int bad = 0;
    for (u64 k : cls) { Tour T = fromCode(k); if (!strongT(T)) bad++; fib[countH(T)]++; }
    printf("order 9: %zu labelled inputs, %zu isomorphism classes (%d non-strong, expected 0)\n", codes.size(), cls.size(), bad);
    printf("strong 9-classes per H:"); for (auto& kv : fib) printf(" %u:%d", kv.first, kv.second); printf("\n");
    return 0;
}
