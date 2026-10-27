// strong_low_spectrum_augment_20261027.cpp -- opus-2026-10-27-S63
// Exact low end of the strong H-spectrum at orders 9, 10 and 11 by augmentation.  Moon (1966, THM-1370's citation): every
// strong tournament on n >= 4 vertices has a vertex whose deletion leaves a strong tournament; the insertion lemma gives
// H(T - v) <= H(T).  Hence every strong n-tournament with H <= B is a one-vertex augmentation of a strong (n-1)-tournament
// with H <= B.  The program (1) canonizes all 2^21 labelled 7-tournaments by brute force over the 5040 relabellings (must
// give 456 classes, 353 strong), (2) augments every strong 7-class by one vertex in all 2^7 ways and canonizes the strong
// results over the 40320 relabellings (must give 6008 strong 8-classes), (3) augments every strong 8-class with H <= B in
// all 2^8 ways and keeps the strong results with H <= B (labelled, with repetitions), (4) augments those in all 2^9 ways to
// order 10 and keeps strong with H <= B, (5) with the flag 11, once more to order 11.  The attained value sets are exactly
// S_n ∩ [1, B] for n = 9, 10 (11); the least value is f(n) when it is <= B, and the odd gaps are exact.
// Build: g++ -O2 -std=c++17 -fopenmp -o sla strong_low_spectrum_augment_20261027.cpp ; run: ./sla B [11]
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <set>
#include <algorithm>
#include <numeric>
using namespace std;
typedef unsigned int u32; typedef unsigned long long u64;

struct Tour { int n; u32 out[11]; };
static inline bool strongT(const Tour& T) { u32 full = (1u << T.n) - 1, r = 1, r2;
    do { r2 = r; u32 t = r; while (t) { int v = __builtin_ctz(t); t &= t - 1; r |= T.out[v]; } } while (r != r2); if (r != full) return false;
    u32 in[11]; for (int v = 0; v < T.n; v++) { in[v] = 0; for (int w = 0; w < T.n; w++) if (T.out[w] >> v & 1) in[v] |= 1u << w; }
    r = 1; do { r2 = r; u32 t = r; while (t) { int v = __builtin_ctz(t); t &= t - 1; r |= in[v]; } } while (r != r2); return r == full; }
static inline u32 countH(const Tour& T) { static thread_local vector<u32> dp; int n = T.n; u32 full = (1u << n) - 1; dp.assign((size_t)(full + 1) * 11, 0);
    for (int v = 0; v < n; v++) dp[(1u << v) * 11 + v] = 1;
    for (u32 S = 1; S <= full; S++) { u32* row = &dp[(size_t)S * 11]; for (int v = 0; v < n; v++) { u32 c = row[v]; if (!c) continue; u32 ext = T.out[v] & ~S; while (ext) { int w = __builtin_ctz(ext); ext &= ext - 1; dp[(size_t)(S | (1u << w)) * 11 + w] += c; } } }
    u32 h = 0; for (int v = 0; v < n; v++) h += dp[(size_t)full * 11 + v]; return h; }
static inline u64 codeOf(const Tour& T) { u64 c = 0; int k = 0; for (int i = 0; i < T.n; i++) for (int j = i + 1; j < T.n; j++, k++) if (T.out[i] >> j & 1) c |= 1ULL << k; return c; }
static Tour fromCode(int n, u64 c) { Tour T; T.n = n; for (int v = 0; v < n; v++) T.out[v] = 0; int k = 0; for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++, k++) { if (c >> k & 1) T.out[i] |= 1u << j; else T.out[j] |= 1u << i; } return T; }
static u64 canon(const Tour& T) {                                 // brute force over all relabellings, n <= 8
    int n = T.n; int p[8]; iota(p, p + n, 0); u64 best = ~0ULL;
    do { Tour R; R.n = n; for (int v = 0; v < n; v++) { u32 m = 0, o = T.out[v]; while (o) { int w = __builtin_ctz(o); o &= o - 1; m |= 1u << p[w]; } R.out[p[v]] = m; }
        u64 c = codeOf(R); if (c < best) best = c; } while (next_permutation(p, p + n));
    return best; }
static inline Tour augment(const Tour& T, u32 m) {                 // new vertex n; bit i of m set: i -> n, else n -> i
    Tour R; R.n = T.n + 1; int v = T.n; R.out[v] = 0; for (int i = 0; i < T.n; i++) { R.out[i] = T.out[i]; if (m >> i & 1) R.out[i] |= 1u << v; else R.out[v] |= 1u << i; } return R; }

static void report(int n, const vector<char>& att, int B, size_t cand) {
    int fn = 0; for (int v = 1; v <= B; v++) if (att[v]) { fn = v; break; }
    printf("order %d: %zu labelled strong candidates with H <= %d; attained values:", n, cand, B); for (int v = 1; v <= B; v++) if (att[v]) printf(" %d", v); printf("\n");
    printf("order %d: f(%d) = %d%s; odd gaps in [f, %d]:", n, n, fn, fn ? "" : " (none <= B)", B); int g = 0; if (fn) for (int v = fn; v <= B; v += 2) if (!att[v]) { printf(" %d", v); g++; } printf(" (%d gaps)\n", g); fflush(stdout); }

int main(int argc, char** argv) {
    int B = argc > 1 ? atoi(argv[1]) : 200; bool do11 = argc > 2 && atoi(argv[2]) == 11;
    // (1) classes at order 7
    set<u64> all7, strong7;
    #pragma omp parallel
    {   set<u64> la, ls;
        #pragma omp for schedule(dynamic, 4096)
        for (long c = 0; c < (1L << 21); c++) { Tour T = fromCode(7, (u64)c); u64 k = canon(T); la.insert(k); if (strongT(T)) ls.insert(k); }
        #pragma omp critical
        { all7.insert(la.begin(), la.end()); strong7.insert(ls.begin(), ls.end()); }
    }
    printf("order 7: %zu classes, %zu strong (expected 456, 353)\n", all7.size(), strong7.size()); fflush(stdout);
    // (2) strong classes at order 8 by augmentation + canonization
    vector<u64> s7(strong7.begin(), strong7.end()); set<u64> strong8;
    #pragma omp parallel
    {   set<u64> ls;
        #pragma omp for schedule(dynamic, 1)
        for (long i = 0; i < (long)s7.size(); i++) { Tour T = fromCode(7, s7[i]); for (u32 m = 0; m < 128; m++) { Tour R = augment(T, m); if (strongT(R)) ls.insert(canon(R)); } }
        #pragma omp critical
        { strong8.insert(ls.begin(), ls.end()); }
    }
    printf("order 8: %zu strong classes (expected 6008)\n", strong8.size()); fflush(stdout);
    vector<u64> s8; for (u64 k : strong8) { Tour T = fromCode(8, k); if (countH(T) <= (u32)B) s8.push_back(k); }
    printf("order 8: %zu strong classes with H <= %d\n", s8.size(), B); fflush(stdout);
    // (3) order 9
    vector<u64> list9; vector<char> att9(B + 1, 0);
    #pragma omp parallel
    {   vector<u64> loc; vector<char> la(B + 1, 0);
        #pragma omp for schedule(dynamic, 1)
        for (long i = 0; i < (long)s8.size(); i++) { Tour T = fromCode(8, s8[i]); for (u32 m = 0; m < 256; m++) { Tour R = augment(T, m); if (!strongT(R)) continue; u32 h = countH(R); if (h <= (u32)B) { la[h] = 1; loc.push_back(codeOf(R)); } } }
        #pragma omp critical
        { list9.insert(list9.end(), loc.begin(), loc.end()); for (int v = 0; v <= B; v++) if (la[v]) att9[v] = 1; }
    }
    report(9, att9, B, list9.size());
    // (4) order 10
    vector<u64> list10; vector<char> att10(B + 1, 0);
    #pragma omp parallel
    {   vector<u64> loc; vector<char> la(B + 1, 0);
        #pragma omp for schedule(dynamic, 16)
        for (long i = 0; i < (long)list9.size(); i++) { Tour T = fromCode(9, list9[i]); for (u32 m = 0; m < 512; m++) { Tour R = augment(T, m); if (!strongT(R)) continue; u32 h = countH(R); if (h <= (u32)B) { la[h] = 1; if (do11) loc.push_back(codeOf(R)); } } }
        #pragma omp critical
        { list10.insert(list10.end(), loc.begin(), loc.end()); for (int v = 0; v <= B; v++) if (la[v]) att10[v] = 1; }
    }
    report(10, att10, B, list10.size());
    if (do11) { vector<char> att11(B + 1, 0); u64 cnt = 0;
        #pragma omp parallel
        {   vector<char> la(B + 1, 0); u64 lc = 0;
            #pragma omp for schedule(dynamic, 16)
            for (long i = 0; i < (long)list10.size(); i++) { Tour T = fromCode(10, list10[i]); for (u32 m = 0; m < 1024; m++) { Tour R = augment(T, m); if (!strongT(R)) continue; u32 h = countH(R); if (h <= (u32)B) { la[h] = 1; lc++; } } }
            #pragma omp critical
            { cnt += lc; for (int v = 0; v <= B; v++) if (la[v]) att11[v] = 1; }
        }
        report(11, att11, B, (size_t)cnt); }
    return 0;
}
