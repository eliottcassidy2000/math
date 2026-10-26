// goldbach_dependent_odd_values_20261026.cpp -- opus-2026-10-26-S61
// Odd aliquot values v <= N and their sources other than semiprimes (the "(O) sources" of S57's
// odd_goldbach_dependence_20261026.py): odd non-square n = p*m with m odd composite and p prime <= spf(m)
// (then s(n) >= 1 + m, so m <= v - 1), and even n = u^2 (u even) or n = 2u^2 (u odd), n <= 2v - 2, plus n = 4.
// An odd v is GOLDBACH-DEPENDENT when it has no (O) source (then it is touched only through semiprimes pq,
// i.e. through a distinct-prime partition v - 1 = p + q, or it is untouchable: v = 5).
// This program records the dependent set and classifies v - 1 by its factorisation pattern, counts dependent
// values per decade, and prints the dependent values below 3000.
// Build: g++ -O2 -std=c++17 -o gdep goldbach_dependent_odd_values_20261026.cpp ; run: ./gdep N   (4N bytes)
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <chrono>
using namespace std;
typedef unsigned long long u64; typedef unsigned int u32;
static double now() { return chrono::duration<double>(chrono::steady_clock::now().time_since_epoch()).count(); }

int main(int argc, char** argv) {
    const u64 N = argc > 1 ? strtoull(argv[1], 0, 10) : 10000000ULL; double T0 = now();
    vector<u32> spf(N + 1, 0);                       // smallest prime factor (linear sieve)
    vector<u32> primes;
    for (u64 i = 2; i <= N; i++) {
        if (!spf[i]) { spf[i] = (u32)i; primes.push_back((u32)i); }
        for (u32 p : primes) { if (p > spf[i] || i * p > N) break; spf[i * p] = p; }
    }
    fprintf(stderr, "sieve done %.1fs (%zu primes)\n", now() - T0, primes.size());
    vector<unsigned char> cov((N + 1) / 2 + 1, 0);   // cov[v>>1] = 1 if odd v has an (O) source
    auto sigma = [&](u64 m) { u64 s = 1; while (m > 1) { u64 p = spf[m], pk = 1, t = 1; while (m % p == 0) { m /= p; pk *= p; t += pk; } s *= t; } return s; };
    // (O1) odd non-square n = p*m, m odd composite, p prime <= spf(m)
    for (u64 m = 9; m <= N; m += 2) {
        u64 q = spf[m]; if (q == m) continue;         // m prime
        u64 sm = sigma(m);
        u64 e = 0, x = m, sq = 1, pk = 1; while (x % q == 0) { x /= q; e++; pk *= q; sq += pk; }
        for (size_t idx = 0; idx < primes.size(); idx++) {
            u64 p = primes[idx]; if (p > q) break;
            u64 s;
            if (p < q) s = (p + 1) * sm - p * m;                                  // p does not divide m
            else       s = (sm / sq) * (sq + pk * q) - q * m;                     // p = q: sigma(q^(e+1)) = sigma(q^e) + q^(e+1)
            if (s <= N && (s & 1)) cov[s >> 1] = 1;
        }
    }
    fprintf(stderr, "O1 done %.1fs\n", now() - T0);
    // (O2) even squares u^2 and twice squares 2u^2 (u of either parity; sigma(n) is odd iff n is a square or twice a square), n <= 2N - 2, and n = 4
    auto sigma_sq = [&](u64 u) { u64 s = 1; while (u > 1) { u64 p = spf[u], e = 0; while (u % p == 0) { u /= p; e++; } u64 pk = 1, t = 1; for (u64 i = 0; i < 2 * e; i++) { pk *= p; t += pk; } s *= t; } return s; };
    for (u64 u = 2; u * u <= 2 * N; u++) {                                    // u <= sqrt(2N) <= N, so spf[u] is valid
        if (u % 2 == 0) { u64 n = u * u; u64 s = sigma_sq(u) - n; if (s <= N && (s & 1)) cov[s >> 1] = 1; }      // even squares
        { u64 n = 2 * u * u; if (n <= 2 * N) { u64 s = (u % 2 == 1 ? 3 * sigma_sq(u) : sigma_sq(u) * 2 * ((1ULL << 1) - 1) + sigma_sq(u)) - n;   // 2u^2, all u: sigma(2u^2) = sigma(2^(2a+1)) sigma(odd part^2)
            if (u % 2 == 0) { u64 a = 0, w = u; while (w % 2 == 0) { w /= 2; a++; } u64 so = 1; { u64 ww = w; while (ww > 1) { u64 pp = spf[ww], e = 0; while (ww % pp == 0) { ww /= pp; e++; } u64 pk = 1, t = 1; for (u64 i = 0; i < 2 * e; i++) { pk *= pp; t += pk; } so *= t; } }
                s = ((1ULL << (2 * a + 2)) - 1) * so - n; }
            if (s <= N && (s & 1)) cov[s >> 1] = 1; } }
    }
    cov[3 >> 1] = 1; // s(4) = 3
    fprintf(stderr, "O2 done %.1fs\n", now() - T0);
    // classification of the dependent set
    vector<u64> dep; map<string, u64> pat; vector<u64> perdecade(12, 0), odddecade(12, 0);
    auto pattern = [&](u64 w) { // factorisation pattern of w = v - 1 (even): "2^a * shape of odd part"
        u64 a = 0; while (w % 2 == 0) { w /= 2; a++; }
        vector<pair<u64,u64>> f; while (w > 1) { u64 p = spf[w], e = 0; while (w % p == 0) { w /= p; e++; } f.push_back({p, e}); }
        string s = "2^" + to_string(a);
        if (f.empty()) return s;
        if (f.size() == 1) return s + (f[0].second == 1 ? string(" * p") : " * p^" + to_string(f[0].second));
        if (f.size() == 2 && f[0].second == 1 && f[1].second == 1) return s + " * p q";
        u64 om = 0; for (auto& x : f) om += x.second;
        return s + " * (" + to_string(f.size()) + " primes, Omega " + to_string(om) + ")";
    };
    for (u64 v = 3; v <= N; v += 2) {
        int dec = 0; for (u64 t = v; t >= 10; t /= 10) dec++;
        odddecade[dec]++;
        if (!cov[v >> 1]) { dep.push_back(v); perdecade[dec]++; pat[pattern(v - 1)]++; }
    }
    fprintf(stderr, "classification done %.1fs\n", now() - T0);
    printf("N = %llu: odd values %llu, Goldbach-dependent (no non-semiprime source) %zu (share %.5f)\n", N, (N - 1) / 2, dep.size(), (double)dep.size() / ((N - 1) / 2));
    for (int d = 0; d < 12; d++) if (odddecade[d]) printf("  decade 10^%d: dependent %llu of %llu odd (share %.5f)\n", d, perdecade[d], odddecade[d], (double)perdecade[d] / odddecade[d]);
    printf("dependent values below 3000: ");
    for (u64 v : dep) { if (v >= 3000) break; printf("%llu ", v); } printf("\n");
    printf("factorisation patterns of v - 1 over the dependent set (count, share):\n");
    vector<pair<u64,string>> ps; for (auto& kv : pat) ps.push_back({kv.second, kv.first}); sort(ps.rbegin(), ps.rend());
    for (size_t i = 0; i < ps.size() && i < 16; i++) printf("  %-34s %llu  (%.4f)\n", ps[i].second.c_str(), ps[i].first, (double)ps[i].first / dep.size());
    map<string, u64> pall; u64 cnt = 0;
    for (u64 v = 3; v <= N; v += 2) { pall[pattern(v - 1)]++; cnt++; }
    printf("the same patterns over all odd v (share): ");
    for (size_t i = 0; i < ps.size() && i < 6; i++) printf("[%s: %.4f] ", ps[i].second.c_str(), (double)pall[ps[i].second] / cnt);
    printf("\n");
    u64 twop = 0, twopdep = 0;
    for (u64 v = 3; v <= N; v += 2) { u64 w = (v - 1) / 2; if (w >= 2 && spf[w] == w && (v - 1) % 4 == 2) { twop++; if (!cov[v >> 1]) twopdep++; } }
    printf("v - 1 = 2 * prime: %llu values, %llu dependent (share %.5f)\n", twop, twopdep, (double)twopdep / twop);
    fprintf(stderr, "all done %.1fs\n", now() - T0);
    return 0;
}
