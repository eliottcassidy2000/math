// perfect_companions_shapes_20261026.cpp -- opus-2026-10-26-S61
// Companions (other aliquot sources) of the even perfect numbers N = 2^(p-1)(2^p - 1), p = 2, 3, 5, 7, 13, 17, 19, 31.
// (1) Row theorem check: a source 2^a q (q an odd prime) exists iff (2^a - 1) | (N - 1), which holds iff a | (p - 1),
//     and q_a = (N + 1 - 2^(a+1))/(2^a - 1) is an odd prime; a = p - 1 gives q = 2^p - 1 (N itself).
// (2) Shape census of even sources m = 2^a M (M odd, Omega(M) <= W) with s(m) = N, enumerating all but the largest
//     prime of M and solving for it: for a stem S (= 2^a times the smaller primes), s(S q) = q s(S) + sigma(S) = N
//     gives q = (N - sigma(S))/s(S), which must be a prime larger than every prime of S.  Shapes covered: 2^a q,
//     2^a q r, 2^a q^2 r, 2^a q^3 r, 2^a q^2 r^2 (direct), 2^a q^k (direct), 2^a q r t, 2^a q^2 r t (W = 3), and
//     2^a q r t u, 2^a q^2 r t u (W = 4) when the bound allows.  Odd sources (odd squares) are not enumerated here
//     (exact for N <= 33550336 in THM-4805: none).
// Build: g++ -O2 -std=c++17 -o pcomp perfect_companions_shapes_20261026.cpp ; run: ./pcomp
#include <cstdio>
#include <cstdint>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace std;
typedef unsigned long long u64; typedef __uint128_t u128;

static u64 mulmod(u64 a, u64 b, u64 m) { return (u64)((u128)a * b % m); }
static u64 powmod(u64 a, u64 e, u64 m) { u64 r = 1; a %= m; while (e) { if (e & 1) r = mulmod(r, a, m); a = mulmod(a, a, m); e >>= 1; } return r; }
static bool isprime(u64 n) {
    if (n < 2) return false;
    for (u64 p : {2ULL,3ULL,5ULL,7ULL,11ULL,13ULL,17ULL,19ULL,23ULL,29ULL,31ULL,37ULL}) { if (n % p == 0) return n == p; }
    u64 d = n - 1; int s = 0; while (!(d & 1)) { d >>= 1; s++; }
    for (u64 a : {2ULL,3ULL,5ULL,7ULL,11ULL,13ULL,17ULL,19ULL,23ULL,29ULL,31ULL,37ULL}) {
        u64 x = powmod(a, d, n); if (x == 1 || x == n - 1) continue;
        bool comp = true; for (int i = 1; i < s; i++) { x = mulmod(x, x, n); if (x == n - 1) { comp = false; break; } }
        if (comp) return false;
    }
    return true;
}
static vector<u64> primes_upto(u64 L) { vector<char> c(L + 1, 0); vector<u64> P; for (u64 i = 2; i <= L; i++) { if (!c[i]) { P.push_back(i); for (u64 j = i * i; j <= L; j += i) c[j] = 1; } } return P; }

struct Hit { u64 m; const char* shape; };

int main() {
    int ps[] = {2, 3, 5, 7, 13, 17, 19, 31};
    vector<u64> primes = primes_upto(2000000);          // small primes for stems
    for (int p : ps) {
        u64 M = (1ULL << p) - 1, N = (1ULL << (p - 1)) * M;
        printf("==== p = %d, N = %llu ====\n", p, N);
        // (1) row theorem
        for (int a = 1; a <= 2 * p + 2 && a <= 62; a++) {
            u64 d = (1ULL << a) - 1; bool div = ((N - 1) % d == 0); bool pred = ((p - 1) % a == 0);
            if (div != pred) printf("  ROW THEOREM VIOLATED at a = %d\n", a);
            if (div) {
                u128 num = (u128)N + 1 - ((u128)1 << (a + 1)); u64 q = (u64)(num / d);
                printf("  a = %d (a | p-1): q_a = %llu -> %s%s\n", a, q, (q % 2 == 1 && isprime(q)) ? "odd prime: source 2^a q = " : "not an odd prime", (q % 2 == 1 && isprime(q)) ? "" : "");
                if (q % 2 == 1 && isprime(q)) printf("      m = %llu%s\n", (1ULL << a) * q, (q == M) ? " (N itself)" : " (COMPANION)");
            }
        }
        // (2) shape census
        vector<Hit> hits; u64 bound = 2 * N - 2;
        auto try_stem = [&](u64 S, u64 sigS, u64 maxp, const char* shape) {
            // s(S q) = q s(S) + sigma(S) = N with q prime > maxp, q not dividing S
            u64 sS = sigS - S; if (sS == 0) return; if (sigS > N) return;
            u64 r = N - sigS; if (r % sS) return; u64 q = r / sS;
            if (q <= maxp || q % 2 == 0 || !isprime(q)) return;
            if (S % q == 0) return;
            hits.push_back({S * q, shape});
        };
        for (int a = 1; (1ULL << a) <= bound; a++) {
            u64 A = 1ULL << a, sigA = (A << 1) - 1;
            // 2^a q
            try_stem(A, sigA, 2, "2^a q");
            // 2^a q^k direct (k >= 2)
            for (u64 q : primes) { if (q < 3) continue; u64 qk = q * q, sq = 1 + q + qk; if (A > bound / qk) break;
                for (int k = 2; ; k++) { u64 m = A * qk; if (m > bound) break; u64 s = sigA * sq - m; if (s == N) hits.push_back({m, "2^a q^k"});
                    if (qk > bound / q) break; qk *= q; sq += qk; } }
            // 2^a q r (q < r): stem 2^a q
            for (u64 q : primes) { if (q < 3) continue; if (A * q > bound / q) break; try_stem(A * q, sigA * (1 + q), q, "2^a q r"); }
            // 2^a q^2 r: stem 2^a q^2
            for (u64 q : primes) { if (q < 3) continue; if (A * q * q > bound / 3) break; try_stem(A * q * q, sigA * (1 + q + q * q), q, "2^a q^2 r"); }
            // 2^a q^3 r
            for (u64 q : primes) { if (q < 3) continue; u64 q3 = q * q * q; if (A > bound / q3 / 3) break; try_stem(A * q3, sigA * (1 + q + q * q + q3), q, "2^a q^3 r"); }
            // 2^a q^2 r^2 direct
            for (u64 q : primes) { if (q < 3) continue; u64 q2 = q * q; if (A * q2 > bound / 9) break;
                for (u64 r : primes) { if (r <= q) continue; u64 r2 = r * r; if (A * q2 > bound / r2) break; u64 m = A * q2 * r2; u64 s = sigA * (1 + q + q2) * (1 + r + r2) - m; if (s == N) hits.push_back({m, "2^a q^2 r^2"}); } }
            // 2^a q r t (q < r < t): stem 2^a q r, with q r <= (bound/A)^(2/3)
            double lim3 = pow((double)bound / A, 2.0 / 3.0);
            if (lim3 < 1e9) {
                for (u64 q : primes) { if (q < 3) continue; if ((double)q * q * q > (double)bound / A) break;
                    for (u64 r : primes) { if (r <= q) continue; if ((double)q * r > lim3) break; try_stem(A * q * r, sigA * (1 + q) * (1 + r), r, "2^a q r t"); } }
                // 2^a q^2 r t: stem 2^a q^2 r
                for (u64 q : primes) { if (q < 3) continue; if ((double)q * q * q * q > (double)bound / A) break;
                    for (u64 r : primes) { if (r <= q) continue; if ((double)q * q * r > lim3) break; try_stem(A * q * q * r, sigA * (1 + q + q * q) * (1 + r), r, "2^a q^2 r t"); } }
            }
            // 2^a q r t u (W = 4): stem 2^a q r t with q r t <= (bound/A)^(3/4), only when cheap
            double lim4 = pow((double)bound / A, 0.75);
            if (lim4 < 2e8) {
                for (u64 q : primes) { if (q < 3) continue; if ((double)q * q * q * q > (double)bound / A) break;
                    for (u64 r : primes) { if (r <= q) continue; if ((double)q * r * r * r > (double)bound / A) break;
                        for (u64 t : primes) { if (t <= r) continue; if ((double)q * r * t > lim4) break; try_stem(A * q * r * t, sigA * (1 + q) * (1 + r) * (1 + t), t, "2^a q r t u"); } } }
            }
        }
        sort(hits.begin(), hits.end(), [](const Hit& x, const Hit& y) { return x.m < y.m; });
        printf("  even sources found (shapes up to the stated bounds): %zu\n", hits.size());
        for (auto& h : hits) printf("    m = %llu  shape %s%s\n", h.m, h.shape, (h.m == N) ? "  (N itself)" : "");
    }
    return 0;
}
