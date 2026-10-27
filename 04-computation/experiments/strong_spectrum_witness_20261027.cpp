// strong_spectrum_witness_20261027.cpp -- opus-2026-10-27-S63
// Witness search for the strong H-spectrum at orders 9 and 10, where exhaustion over all labelled tournaments is out of
// reach.  Phase 1: random tournaments at all densities (every pair of the transitive tournament reversed independently with
// a probability p drawn uniformly from [0, 1/2] per sample).  Phase 2: for every odd value v <= vmax still without a
// witness, hill-climbing on |H - v| over single-pair reversals from random strong starts.  Every witness is an exact
// certificate: the program prints "v code" with the pair code (bit k = 1 iff i -> j for the k-th pair (i, j), i < j,
// lexicographic), so an independent program can recount the Hamiltonian paths and check strong connectivity.
// Only membership "v is in S_n" is certified; a value without a witness is not shown to be absent.
// Build: g++ -O2 -std=c++17 -fopenmp -o ssw strong_spectrum_witness_20261027.cpp ; run: ./ssw n samples vmax certificate [seed]
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;
typedef unsigned int u32; typedef unsigned long long u64;

static int n, M; static vector<pair<int,int>> pairs; static u32 fullV;
static inline u64 rnd(u64& s) { s ^= s << 13; s ^= s >> 7; s ^= s << 17; return s * 0x2545F4914F6CDD1DULL; }
static inline void decode(u64 code, u32* out) { for (int v = 0; v < n; v++) out[v] = 0;
    for (int k = 0; k < M; k++) { int i = pairs[k].first, j = pairs[k].second; if (code >> k & 1) out[i] |= 1u << j; else out[j] |= 1u << i; } }
static inline bool strongT(const u32* out) {
    u32 r = 1, r2; do { r2 = r; u32 t = r; while (t) { int v = __builtin_ctz(t); t &= t - 1; r |= out[v]; } } while (r != r2); if (r != fullV) return false;
    u32 in[10]; for (int v = 0; v < n; v++) { in[v] = 0; for (int w = 0; w < n; w++) if (out[w] >> v & 1) in[v] |= 1u << w; }
    r = 1; do { r2 = r; u32 t = r; while (t) { int v = __builtin_ctz(t); t &= t - 1; r |= in[v]; } } while (r != r2); return r == fullV; }
static inline u32 countH(const u32* out) { static thread_local u32 dp[1 << 10][10];
    for (u32 S = 1; S <= fullV; S++) for (int v = 0; v < n; v++) dp[S][v] = 0;
    for (int v = 0; v < n; v++) dp[1u << v][v] = 1;
    for (u32 S = 1; S <= fullV; S++) for (int v = 0; v < n; v++) { u32 c = dp[S][v]; if (!c) continue; u32 ext = out[v] & ~S; while (ext) { int w = __builtin_ctz(ext); ext &= ext - 1; dp[S | (1u << w)][w] += c; } }
    u32 h = 0; for (int v = 0; v < n; v++) h += dp[fullV][v]; return h; }

int main(int argc, char** argv) {
    n = atoi(argv[1]); long samples = atol(argv[2]); int vmax = atoi(argv[3]); const char* cert = argv[4]; u64 seed = argc > 5 ? strtoull(argv[5], 0, 10) : 0x9E3779B97F4A7C15ULL;
    if (n < 5 || n > 10) { fprintf(stderr, "5 <= n <= 10\n"); return 1; }
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) pairs.push_back({i, j}); M = (int)pairs.size(); fullV = (1u << n) - 1; u64 fullP = (M == 64) ? ~0ULL : ((1ULL << M) - 1);
    vector<u64> wit(vmax + 1, 0); vector<char> have(vmax + 1, 0);
    // phase 1
    #pragma omp parallel
    {   u64 s = seed ^ (0x9E3779B97F4A7C15ULL * (u64)(omp_get_thread_num() + 1)); for (int i = 0; i < 10; i++) rnd(s); u32 out[10]; vector<u64> lw(vmax + 1, 0); vector<char> lh(vmax + 1, 0);
        #pragma omp for schedule(static)
        for (long t = 0; t < samples; t++) { double p = (double)(rnd(s) >> 11) / 9007199254740992.0 * 0.5; u64 code = fullP;
            for (int k = 0; k < M; k++) if ((double)(rnd(s) >> 11) / 9007199254740992.0 < p) code ^= 1ULL << k;
            decode(code, out); if (!strongT(out)) continue; u32 h = countH(out); if (h <= (u32)vmax && !lh[h]) { lh[h] = 1; lw[h] = code; } }
        #pragma omp critical
        { for (int v = 0; v <= vmax; v++) if (lh[v] && !have[v]) { have[v] = 1; wit[v] = lw[v]; } }
    }
    int covered = 0, missing = 0; vector<int> miss; for (int v = 1; v <= vmax; v += 2) { if (have[v]) covered++; else { missing++; miss.push_back(v); } }
    printf("n = %d: phase 1, %ld samples: %d odd values <= %d have strong witnesses, %d missing\n", n, samples, covered, vmax, missing); fflush(stdout);
    // phase 2: local search for the missing values (parallel over values)
    int found2 = 0;
    #pragma omp parallel for schedule(dynamic, 1) reduction(+:found2)
    for (int idx = 0; idx < (int)miss.size(); idx++) { int v = miss[idx]; u64 s = seed ^ (0xD1B54A32D192ED03ULL * (u64)(v + 7)); for (int i = 0; i < 10; i++) rnd(s); u32 out[10]; bool ok = false; u64 best = 0;
        for (int restart = 0; restart < 40 && !ok; restart++) {
            double p = (double)(rnd(s) >> 11) / 9007199254740992.0 * 0.5; u64 code = fullP; for (int k = 0; k < M; k++) if ((double)(rnd(s) >> 11) / 9007199254740992.0 < p) code ^= 1ULL << k;
            decode(code, out); if (!strongT(out)) continue; long cur = labs((long)countH(out) - v);
            for (int step = 0; step < 4000 && cur != 0; step++) { int k = rnd(s) % M; u64 c2 = code ^ (1ULL << k); decode(c2, out); if (!strongT(out)) continue; long d = labs((long)countH(out) - v);
                if (d <= cur) { code = c2; cur = d; } }
            if (cur == 0) { ok = true; best = code; } }
        if (ok) { found2++;
            #pragma omp critical
            { have[v] = 1; wit[v] = best; } } }
    covered = 0; vector<int> still; for (int v = 1; v <= vmax; v += 2) { if (have[v]) covered++; else still.push_back(v); }
    printf("n = %d: phase 2 found %d more; total %d odd values <= %d witnessed; %zu without witness:", n, found2, covered, vmax, still.size()); for (int v : still) printf(" %d", v); printf("\n");
    FILE* f = fopen(cert, "w"); fprintf(f, "# strong_spectrum_witness_20261027: n = %d, vmax = %d, seed %llu; lines: v code (bit k = 1 iff i -> j for the k-th pair (i,j), i<j, lexicographic); every listed code is a strong tournament with exactly v Hamiltonian paths\n", n, vmax, seed);
    for (int v = 1; v <= vmax; v += 2) if (have[v]) fprintf(f, "%d %llu\n", v, wit[v]); fclose(f);
    // re-verify every certificate line from scratch
    int bad = 0; { u32 out[10]; for (int v = 1; v <= vmax; v += 2) if (have[v]) { decode(wit[v], out); if (!strongT(out) || countH(out) != (u32)v) bad++; } }
    printf("certificate %s written; re-verification: %d bad lines\n", cert, bad);
    return 0;
}
