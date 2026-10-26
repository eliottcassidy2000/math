// aliquot_sources_preregistered_20261026.cpp -- opus-2026-10-26-S60
// Pre-registered tests for the S60 note (predictions written in the note before this output was read):
//  (a) twin-collision values v = p + 1 (p - 2, p prime) in (10^4, 10^5]: share with exactly the two twin sources
//      {2(p-2), p^2}, mean number of sources (parity lemma: even sources m <= 2v - 2, odd sources are odd squares q^2, q < v);
//  (b) the fifth perfect number N = 33550336: its full source set (is it its own sole source, like 28?);
//  (c) abc hits (a + b = c coprime, rad(abc) < c) with c <= 10^5: count, share containing an odd perfect power,
//      share whose even term is a power of two, share with 1 as a term.
// Build: g++ -O2 -std=c++17 -o aliquot_sources_pre aliquot_sources_preregistered_20261026.cpp ; run: ./aliquot_sources_pre
#include <cstdio>
#include <cstdint>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
using namespace std;
typedef unsigned long long u64;
typedef unsigned int u32;

static vector<u32> spf_sieve(u32 N) {            // smallest prime factor
    vector<u32> spf(N + 1, 0);
    for (u32 i = 2; i <= N; i++) if (!spf[i]) { for (u64 j = i; j <= N; j += i) if (!spf[j]) spf[j] = i; }
    return spf;
}
static u64 sigma_of_square(u32 q, const vector<u32>& spf) {   // sigma(q^2) by factoring q
    u64 res = 1;
    while (q > 1) { u32 p = spf[q]; int e = 0; while (q % p == 0) { q /= p; e++; }
        u64 pw = 1, sum = 1; for (int i = 0; i < 2 * e; i++) { pw *= p; sum += pw; } res *= sum; }
    return res;
}
static bool is_perfect_power(u64 n) {
    if (n < 4) return false;
    for (int k = 2; k < 64; k++) {
        double r = pow((double)n, 1.0 / k); if (r < 1.5) break;
        for (long long c = (long long)r - 1; c <= (long long)r + 1; c++) {
            if (c < 2) continue; u64 v = 1; bool over = false;
            for (int i = 0; i < k; i++) { v *= (u64)c; if (v > n) { over = true; break; } }
            if (!over && v == n) return true;
        }
    }
    return false;
}

int main() {
    // ---- (a) twin values in (10^4, 10^5] ----
    {
        const u32 V = 100000, M = 2 * V;                  // even sources up to 2V - 2
        vector<u64> sig(M + 1, 0);
        for (u32 d = 1; d <= M; d++) for (u32 m = d; m <= M; m += d) sig[m] += d;
        vector<u32> spf = spf_sieve(V);
        vector<u32> cnt(V + 1, 0);                        // number of sources of each even value v <= V
        for (u32 m = 2; m <= M; m += 2) { u64 s = sig[m] - m; if (s <= V) cnt[s]++; }
        for (u32 q = 3; q < V; q += 2) { u64 s = sigma_of_square(q, spf) - (u64)q * q; if (s <= V) cnt[s]++; }
        int tw = 0, two = 0; double sum = 0; int mx = 0;
        for (u32 p = 10001; p <= V - 1; p += 2) {
            if (spf[p] != p || spf[p - 2] != p - 2) continue;   // p and p - 2 prime
            u32 v = p + 1; tw++; sum += cnt[v]; if (cnt[v] == 2) two++; if ((int)cnt[v] > mx) mx = cnt[v];
        }
        printf("(a) twin values p+1 in (10^4, 10^5]: %d values; exactly two sources: %d (share %.4f); mean sources %.3f; max %d\n", tw, two, (double)two / tw, sum / tw, mx);
        // also the full range <= 10^5 for reference
        int tw2 = 0, two2 = 0; double sum2 = 0;
        for (u32 p = 5; p <= V - 1; p += 2) { if (spf[p] != p || spf[p - 2] != p - 2) continue; tw2++; sum2 += cnt[p + 1]; if (cnt[p + 1] == 2) two2++; }
        printf("    twin values p+1 <= 10^5: %d; exactly two: %d (share %.4f); mean %.3f\n", tw2, two2, (double)two2 / tw2, sum2 / tw2);
        // untouchable check in the same table: even values <= 10^5 with no source at all; count and the first 12
        int unt = 0; printf("    even untouchables <= 10^5 (parity lemma, exact): ");
        for (u32 v = 2; v <= V; v += 2) if (cnt[v] == 0) { unt++; if (unt <= 12) printf("%u ", v); }
        printf("... count %d (A057978: 1 untouchable even below 10 is 2; the odd 5 is not in this even table)\n", unt);
    }
    // ---- (b) the fifth perfect number ----
    {
        const u64 N = 33550336ULL; const u32 M = (u32)(2 * N - 2);
        vector<u32> s(M + 1, 0);                           // s(m) = sigma(m) - m fits u32 here (sigma(m) < 5m < 2^32)
        for (u32 d = 1; d <= M; d++) for (u64 m = d; m <= M; m += d) s[(u32)m] += d;
        vector<u64> src;
        for (u32 m = 2; m <= M; m += 2) if ((u64)s[m] - m == N) src.push_back(m);
        vector<u32> spf = spf_sieve((u32)N);
        for (u32 q = 3; q < N; q += 2) { u64 sq = sigma_of_square(q, spf) - (u64)q * q; if (sq == N) src.push_back((u64)q * q); }
        sort(src.begin(), src.end());
        printf("(b) s^-1(33550336) = {");
        for (size_t i = 0; i < src.size(); i++) printf("%s%llu", i ? ", " : "", src[i]);
        printf("} (%zu sources; self-only = %s)\n", src.size(), (src.size() == 1 && src[0] == N) ? "YES" : "NO");
    }
    // ---- (c) abc hits with c <= 10^5 ----
    {
        const u32 C = 100000;
        vector<u32> rad(C + 1, 1); vector<char> comp(C + 1, 0);
        for (u32 p = 2; p <= C; p++) { if (comp[p]) continue; for (u32 m = p; m <= C; m += p) { comp[m] = 1; rad[m] *= p; } comp[p] = 0; }
        // comp[] is now 1 for all m >= 2 that were marked (every m is marked by its own prime factors) -- not used further
        long long hits = 0, oddpp = 0, pow2 = 0, one = 0;
        double bestq = 0; u32 bc = 0, ba = 0;
        for (u32 c = 3; c <= C; c++) {
            for (u32 a = 1; a <= c / 2; a++) {
                u32 b = c - a; if (std::gcd(a, b) != 1) continue;
                double R = (double)rad[a] * rad[b] * rad[c];
                if (R < (double)c) {
                    hits++;
                    if ((a % 2 == 1 && is_perfect_power(a)) || (b % 2 == 1 && is_perfect_power(b)) || (c % 2 == 1 && is_perfect_power(c))) oddpp++;
                    u32 ev = (a % 2 == 0) ? a : (b % 2 == 0 ? b : c); if ((ev & (ev - 1)) == 0) pow2++;
                    if (a == 1) one++;
                    double q = log((double)c) / log(R); if (q > bestq) { bestq = q; bc = c; ba = a; }
                }
            }
        }
        printf("(c) abc hits c <= 10^5: %lld; with an odd perfect power: %lld (share %.4f); even term a power of two: %lld (share %.4f); a = 1: %lld (share %.4f); best quality %.4f at (%u, %u, %u)\n",
               hits, oddpp, (double)oddpp / hits, pow2, (double)pow2 / hits, one, (double)one / hits, bestq, ba, bc - ba, bc);
    }
    return 0;
}
