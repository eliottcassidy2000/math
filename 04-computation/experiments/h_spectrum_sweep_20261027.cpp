// h_spectrum_sweep_20261027.cpp -- opus-2026-10-27-S63
// Exact raw-H spectrum of tournaments on n <= 8 vertices.  For every labelled tournament (2^(n(n-1)/2) codes; bit for the
// pair (i, j), i < j, is 1 iff i -> j) the program computes H = the number of Hamiltonian paths (subset DP) and whether
// the tournament is strong, and accumulates, per value H: the number of labelled tournaments, the number of isomorphism
// classes (Burnside: only permutations with all cycles odd fix a tournament; for each such cycle type the fixed tournaments
// are enumerated over their pair orbits and weighted by the number of permutations of that type), and the number of strong
// classes.  Printed: the spectrum table, the sets A_n (all values) and S_n (strong values), the minimum strong value f(n),
// the maximum, and the odd gaps of S_n between f(n) and the maximum.
// Build: g++ -O2 -std=c++17 -fopenmp -o hsw h_spectrum_sweep_20261027.cpp ; run: ./hsw n
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <map>
#include <algorithm>
#include <numeric>
#include <functional>
using namespace std;
typedef unsigned long long u64; typedef unsigned int u32;

static int n, M; static int pidx[8][8]; static vector<pair<int,int>> pairs;
static const int VMAX = 1 << 16;

static inline u32 countH(const u32* out) {
    u32 full = (1u << n) - 1; static thread_local u32 dp[256 * 8];
    memset(dp, 0, sizeof(u32) * ((size_t)(full + 1) * 8));
    for (int v = 0; v < n; v++) dp[(1u << v) * 8 + v] = 1;
    for (u32 S = 1; S <= full; S++) { u32* row = dp + S * 8;
        for (int v = 0; v < n; v++) { u32 c = row[v]; if (!c) continue;
            u32 ext = out[v] & ~S; while (ext) { int w = __builtin_ctz(ext); ext &= ext - 1; dp[(S | (1u << w)) * 8 + w] += c; } } }
    u32 h = 0; for (int v = 0; v < n; v++) h += dp[full * 8 + v]; return h; }

static inline bool strong(const u32* out, const u32* in) {
    u32 full = (1u << n) - 1, r = 1, r2;
    do { r2 = r; u32 t = r; while (t) { int v = __builtin_ctz(t); t &= t - 1; r |= out[v]; } } while (r != r2);
    if (r != full) return false;
    r = 1; do { r2 = r; u32 t = r; while (t) { int v = __builtin_ctz(t); t &= t - 1; r |= in[v]; } } while (r != r2);
    return r == full; }

static inline void decode(u64 code, u32* out, u32* in) {
    for (int v = 0; v < n; v++) out[v] = in[v] = 0;
    for (int k = 0; k < M; k++) { int i = pairs[k].first, j = pairs[k].second;
        if (code >> k & 1) { out[i] |= 1u << j; in[j] |= 1u << i; } else { out[j] |= 1u << i; in[i] |= 1u << j; } } }

int main(int argc, char** argv) {
    n = argc > 1 ? atoi(argv[1]) : 6; if (n < 2 || n > 8) { fprintf(stderr, "2 <= n <= 8\n"); return 1; }
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { pidx[i][j] = pidx[j][i] = (int)pairs.size(); pairs.push_back({i, j}); }
    M = (int)pairs.size(); u64 N = 1ULL << M;
    vector<u64> lab(VMAX, 0), labS(VMAX, 0);        // labelled counts, all / strong, per H value
    // ---- full labelled sweep (identity cycle type) ----
    #pragma omp parallel
    {   vector<u64> l(VMAX, 0), ls(VMAX, 0); u32 out[8], in[8];
        #pragma omp for schedule(dynamic, 8192)
        for (long long c = 0; c < (long long)N; c++) { decode((u64)c, out, in); u32 h = countH(out); bool s = strong(out, in); l[h]++; if (s) ls[h]++; }
        #pragma omp critical
        { for (int v = 0; v < VMAX; v++) { lab[v] += l[v]; labS[v] += ls[v]; } }
    }
    // ---- Burnside over odd-cycle types ----
    // enumerate partitions of n into odd parts (cycle lengths), excluding the identity (handled by the sweep)
    vector<vector<int>> types; vector<int> cur;
    function<void(int,int)> gen = [&](int rem, int maxp) { if (rem == 0) { types.push_back(cur); return; }
        int p0 = min(rem, maxp); if (p0 % 2 == 0) p0--;
        for (int p = p0; p >= 1; p -= 2) { cur.push_back(p); gen(rem - p, p); cur.pop_back(); } };
    gen(n, n % 2 ? n : n - 1);
    u64 nfact = 1; for (int i = 2; i <= n; i++) nfact *= i;
    vector<u64> cls(VMAX, 0), clsS(VMAX, 0);        // n! * (class counts), accumulated
    for (int v = 0; v < VMAX; v++) { cls[v] = lab[v]; clsS[v] = labS[v]; }     // identity contributes with multiplicity 1
    for (auto& t : types) { bool ident = true; for (int p : t) if (p != 1) ident = false; if (ident) continue;
        // representative permutation g with cycles of the given lengths on consecutive vertices
        int g[8]; int pos = 0; for (int p : t) { for (int k = 0; k < p; k++) g[pos + k] = pos + (k + 1) % p; pos += p; }
        // multiplicity: n! / prod(len_i) / prod(mult_j!)
        u64 mult = nfact; map<int,int> mm; for (int p : t) { mult /= p; mm[p]++; } for (auto& kv : mm) { for (int i = 2; i <= kv.second; i++) mult /= i; }
        // pair orbits under g
        vector<int> orb(M, -1); vector<vector<pair<int,int>>> orbits;      // each orbit: list of (pair index, sign) where sign 0 = same direction as the representative
        for (int k = 0; k < M; k++) { if (orb[k] >= 0) continue; int id = (int)orbits.size(); vector<pair<int,int>> o;
            int a = pairs[k].first, b = pairs[k].second;           // representative oriented a -> b ; image pairs (g^t a, g^t b)
            int x = a, y = b; do { int kk = pidx[x][y]; int sign = (x < y) ? 0 : 1; if (orb[kk] < 0) { orb[kk] = id; o.push_back({kk, sign}); } x = g[x]; y = g[y]; } while (x != a || y != b);
            orbits.push_back(o); }
        int K = (int)orbits.size(); u64 fixedCount = 0; u32 out[8], in[8];
        for (u64 asg = 0; asg < (1ULL << K); asg++) { u64 code = 0;
            for (int o = 0; o < K; o++) { int dir = asg >> o & 1;       // dir 0: representative a -> b
                for (auto& ps : orbits[o]) { int bit = dir ^ ps.second; if (bit == 0) code |= 1ULL << ps.first; } }   // bit 1 iff lower -> higher
            // explanation: representative a < b oriented a->b means bit 1 at pair k; image (x, y) with x > y oriented x -> y means bit 0 at pidx[x][y]
            decode(code, out, in); u32 h = countH(out); bool s = strong(out, in); cls[h] += mult; if (s) clsS[h] += mult; fixedCount++; }
        printf("# cycle type"); for (int p : t) printf(" %d", p); printf(": %llu permutations, %d pair orbits, %llu fixed tournaments\n", mult, K, fixedCount);
    }
    // ---- report ----
    u64 totalL = 0, totalC = 0, totalCS = 0; int fn = 0, mx = 0;
    printf("n = %d: %llu labelled tournaments\n", n, N);
    printf("H\tlabelled\tclasses\tstrong_classes\tstrong_labelled\n");
    for (int v = 0; v < VMAX; v++) if (lab[v]) { if (cls[v] % nfact || clsS[v] % nfact) { printf("BURNSIDE NON-INTEGRAL at H = %d\n", v); return 2; }
        printf("%d\t%llu\t%llu\t%llu\t%llu\n", v, lab[v], cls[v] / nfact, clsS[v] / nfact, labS[v]); totalL += lab[v]; totalC += cls[v] / nfact; totalCS += clsS[v] / nfact; mx = v; if (labS[v] && !fn) fn = v; }
    printf("totals: labelled %llu, classes %llu, strong classes %llu\n", totalL, totalC, totalCS);
    printf("A_%d (all values):", n); for (int v = 0; v < VMAX; v++) if (lab[v]) printf(" %d", v); printf("\n");
    printf("S_%d (strong values):", n); for (int v = 0; v < VMAX; v++) if (labS[v]) printf(" %d", v); printf("\n");
    printf("f(%d) = %d (least strong value), max = %d; odd gaps of S_%d in [f, max]:", n, fn, mx, n);
    int gaps = 0; for (int v = fn; v <= mx; v += 2) if (!labS[v]) { printf(" %d", v); gaps++; } printf(" (%d gaps)\n", gaps);
    return 0;
}
