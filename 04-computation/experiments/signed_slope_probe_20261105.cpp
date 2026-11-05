// signed_slope_probe_20261105.cpp -- opus-2026-11-05-S80
// Small-witness coverage of the signed one-apex response a n - b g(n) over odd n < 2^N, for a list of slopes (a, b).
// g(1) = 1, g(2n) = g(n), g(2n-1) = g(n) + n - (n mod 2): the Hamiltonian-path count of the one-apex tournament with
// ordered two-cover count n (Session198 (7); HYP-9468's g with n = v + 1).  Checked below against the one-apex automaton
// (H, A, B): the chain word with binary value v = n - 1 has H = g(n).
// Usage: probe N L a1 b1 a2 b2 ...   -> for each slope: parity of the image, count of attained targets in [-L, L] of that
//        parity, the first misses, the residue classes modulo 4 and 8 of the positive values, and the least n for small targets.
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <map>
#include <set>
#include <cstdint>
using namespace std;
typedef long long ll;
static ll oneApexH(ll v) {   // H of the chain word with binary value v (most significant bit first), via (H, A, B)
    if (v == 0) return 1;
    int bits = 0; while ((v >> bits) > 0) bits++;
    ll H = 1, A = 0, B = 0;
    for (int i = bits - 1; i >= 0; i--) { int bit = (v >> i) & 1; if (bit) { A = 2 * A + B; B = 1; } else { H += 2 * (A + B); A = 2 * A + B; B = 0; } }
    return H; }
int main(int argc, char** argv) {
    int N = atoi(argv[1]); ll L = atoll(argv[2]); ll top = 1LL << N;
    vector<ll> g(top + 1, 0); g[1] = 1;
    for (ll n = 2; n <= top; n++) g[n] = (n % 2 == 0) ? g[n / 2] : g[(n + 1) / 2] + (n + 1) / 2 - (((n + 1) / 2) % 2);
    // check the dictionary g(n) = H(word with value n - 1) for n <= 2^12
    ll bad = 0; for (ll n = 1; n <= 4096; n++) if (g[n] != oneApexH(n - 1)) bad++;
    printf("dictionary check g(n) = H(one-apex word of value n-1) for n <= 4096: %lld mismatches\n", bad);
    for (int i = 3; i + 1 < argc; i += 2) { ll a = atoll(argv[i]), b = atoll(argv[i + 1]);
        int par = (int)(((a - b) % 2 + 2) % 2);   // parity of a n - b g(n) for odd n (D odd)
        vector<char> att(2 * L + 1, 0); vector<ll> least(2 * L + 1, 0);
        map<int, ll> mod4pos, mod8pos; ll posCount = 0, negCount = 0;
        for (ll n = 3; n < top; n += 2) { ll val = a * n - b * g[n];
            if (val > 0) { posCount++; mod4pos[(int)(val % 4)]++; mod8pos[(int)(val % 8)]++; } else negCount++;
            if (val >= -L && val <= L) { if (!att[val + L]) { att[val + L] = 1; least[val + L] = n; } } }
        ll cnt = 0, tot = 0; vector<ll> miss; for (ll t = -L; t <= L; t++) { if (((t % 2) + 2) % 2 != par) continue; tot++; if (att[t + L]) cnt++; else miss.push_back(t); }
        printf("slope %lld/%lld (parity %d, %s): attained %lld of %lld targets in [-%lld, %lld]; positive values %lld, nonpositive %lld; positive residues mod 4:", a, b, par, (2 * a > b && a < b) ? "inside (1/2,1)" : "outside", cnt, tot, L, L, posCount, negCount);
        for (auto& kv : mod4pos) printf(" %d:%lld", kv.first, kv.second); printf("; mod 8:"); for (auto& kv : mod8pos) printf(" %d:%lld", kv.first, kv.second);
        printf("\n   first misses:"); for (size_t j = 0; j < miss.size() && j < 24; j++) printf(" %lld", miss[j]); if (miss.size() > 24) printf(" ... (%zu)", miss.size());
        printf("\n   least n for targets 1..15 and -1..-15:"); for (ll t = 1; t <= 15; t += 2) if (t <= L && ((t % 2) + 2) % 2 == par) printf(" %lld:%lld", t, att[t + L] ? least[t + L] : -1); for (ll t = -1; t >= -15; t -= 2) if (((t % 2) + 2) % 2 == par) printf(" %lld:%lld", t, att[t + L] ? least[t + L] : -1); printf("\n"); }
    return 0; }
