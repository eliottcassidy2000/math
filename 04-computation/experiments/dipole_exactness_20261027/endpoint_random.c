/* endpoint_random.c -- mac-mini-2026-10-27-dipole
 * Random tests of the endpoint statements beyond the exhaustive orders:
 *   (C7) some Hamiltonian path runs from a maximum-score vertex to a minimum-score vertex;
 *   (C9) (odd n) some Hamiltonian path runs from a vertex above (n-1)/2 to one below (= R-descent, regular target);
 *   (K2) kappa >= 2 implies an L-descending move (some Hamiltonian path joins a pair of D_L).
 * Samplers (mode): 0 uniform; 1 near-regular (regular circulant, many random 3-cycle reversals, then f random arc
 * flips, f uniform in 1..4: ties among the extreme scores are typical); 2 planted ties: like 1 but then forces every
 * maximum-score vertex to beat every minimum-score vertex (the tied configuration not covered by Lu's theorem).
 * Build: cc -O3 -march=native endpoint_random.c -o endpoint_random ; Run: ./endpoint_random n samples mode seed
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
static int n; static uint32_t out_[20]; static uint32_t pi_[20];
static uint32_t *dp;   /* dp[S*n + v] = bitmask of starts */
static uint64_t rng_s;
static uint64_t rnd(void) { rng_s ^= rng_s << 13; rng_s ^= rng_s >> 7; rng_s ^= rng_s << 17; return rng_s; }
static void setarc(int u, int v) { out_[u] |= 1u << v; out_[v] &= ~(1u << u); }
static void endpoints(void) {
    int full = (1 << n) - 1;
    memset(dp, 0, sizeof(uint32_t) * (size_t)(full + 1) * n);
    for (int v = 0; v < n; v++) dp[(size_t)(1 << v) * n + v] = 1u << v;
    for (int S = 1; S <= full; S++) for (int v = 0; v < n; v++) { uint32_t st = dp[(size_t)S * n + v]; if (!st) continue;
        uint32_t nx = out_[v] & ~S; while (nx) { int w = __builtin_ctz(nx); nx &= nx - 1; dp[(size_t)(S | 1 << w) * n + w] |= st; } }
    for (int a = 0; a < n; a++) pi_[a] = 0;
    for (int b = 0; b < n; b++) { uint32_t st = dp[(size_t)full * n + b]; while (st) { int a = __builtin_ctz(st); st &= st - 1; pi_[a] |= 1u << b; } }
}
static int strongMask(uint32_t vs) {
    if (!vs) return 1; int start = __builtin_ctz(vs);
    uint32_t seen = 1u << start, fr = seen;
    while (fr) { uint32_t nx = 0; for (uint32_t f = fr; f; f &= f - 1) nx |= out_[__builtin_ctz(f)] & vs; nx &= ~seen; seen |= nx; fr = nx; }
    if (seen != vs) return 0;
    seen = 1u << start; fr = seen;
    while (fr) { uint32_t nx = 0; for (uint32_t f = fr; f; f &= f - 1) { int v = __builtin_ctz(f); for (uint32_t w = vs; w; w &= w - 1) { int u = __builtin_ctz(w); if (out_[u] >> v & 1) nx |= 1u << u; } } nx &= ~seen; seen |= nx; fr = nx; }
    return seen == vs;
}
static int Lof(const int *s) { int t[20]; memcpy(t, s, n * sizeof(int));
    for (int i = 1; i < n; i++) { int v = t[i], j = i - 1; while (j >= 0 && t[j] > v) { t[j + 1] = t[j]; j--; } t[j + 1] = v; }
    int l = 0; for (int i = 0; i < n; i++) l += abs(t[i] - i); return l / 2; }
int main(int argc, char **argv) {
    n = atoi(argv[1]); long samples = atol(argv[2]); int mode = atoi(argv[3]); rng_s = 0x9E3779B97F4A7C15ULL ^ (uint64_t)atol(argv[4]);
    dp = malloc(sizeof(uint32_t) * ((size_t)1 << n) * n);
    uint32_t full = (1u << n) - 1;
    long c7_bad = 0, c9_bad = 0, c9_n = 0, k2_n = 0, k2_bad = 0, strong_n = 0, tied = 0;
    for (long it = 0; it < samples; it++) {
        for (int i = 0; i < n; i++) out_[i] = 0;
        if (mode == 0) { for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (rnd() & 1) setarc(i, j); else setarc(j, i); } }
        else {
            int m = n / 2;   /* rotational base: i -> i+1..i+m (n odd regular; n even: almost regular, ties broken) */
            for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { int d = j - i; if (d <= m && !(n % 2 == 0 && d == m && i >= m)) setarc(i, j); else setarc(j, i); }
            for (int r = 0; r < 20 * n; r++) {   /* reverse random directed 3-cycles (scores preserved) */
                int a = rnd() % n, b = rnd() % n, c = rnd() % n; if (a == b || b == c || a == c) continue;
                if ((out_[a] >> b & 1) && (out_[b] >> c & 1) && (out_[c] >> a & 1)) { setarc(b, a); setarc(c, b); setarc(a, c); }
            }
            int f = 1 + rnd() % 4; for (int r = 0; r < f; r++) { int a = rnd() % n, b = rnd() % n; if (a == b) continue; if (out_[a] >> b & 1) setarc(b, a); else setarc(a, b); }
            if (mode == 2) {
                int s[20], mn = 99, mx = -1; for (int v = 0; v < n; v++) { s[v] = __builtin_popcount(out_[v]); if (s[v] < mn) mn = s[v]; if (s[v] > mx) mx = s[v]; }
                for (int x = 0; x < n; x++) for (int y = 0; y < n; y++) if (s[x] == mx && s[y] == mn && x != y) setarc(x, y);
            }
        }
        int s[20], mn = 99, mx = -1; for (int v = 0; v < n; v++) { s[v] = __builtin_popcount(out_[v]); if (s[v] < mn) mn = s[v]; if (s[v] > mx) mx = s[v]; }
        if (mn == mx) continue;
        endpoints();
        int nmax = 0, nmin = 0; uint32_t MN = 0, LO = 0, HI = 0;
        for (int v = 0; v < n; v++) { nmax += s[v] == mx; nmin += s[v] == mn; if (s[v] == mn) MN |= 1u << v; if (2 * s[v] < n - 1) LO |= 1u << v; if (2 * s[v] > n - 1) HI |= 1u << v; }
        if (nmax > 1 || nmin > 1) tied++;
        int ok7 = 0; for (int a = 0; a < n; a++) if (s[a] == mx && (pi_[a] & MN)) ok7 = 1;
        if (!ok7) { c7_bad++; if (c7_bad <= 3) { printf("C7 COUNTEREXAMPLE n=%d:", n); for (int i = 0; i < n; i++) printf(" %08x", out_[i]); printf("\n"); } }
        if (n % 2 == 1 && LO && HI) { c9_n++; int ok9 = 0; for (int a = 0; a < n; a++) if ((HI >> a & 1) && (pi_[a] & LO)) ok9 = 1; if (!ok9) c9_bad++; }
        if (strongMask(full)) {
            strong_n++;
            int kap = n - 1; for (uint32_t S = 1; S < full && kap >= 2; S++) { int sz = __builtin_popcount(S); if (sz >= 2 || sz > n - 2) continue; if (!strongMask(full & ~S)) kap = sz; }
            if (kap >= 2) {   /* kappa >= 2 (only single-vertex cuts tested; kap here is 1 or >= 2) */
                k2_n++;
                int L0 = Lof(s), de = 0, t[20];
                for (int a = 0; a < n && !de; a++) for (uint32_t e = pi_[a]; e; e &= e - 1) { int b = __builtin_ctz(e);
                    for (int v = 0; v < n; v++) t[v] = n - 1 - s[v]; t[a]++; t[b]--; if (Lof(t) == L0 - 1) { de = 1; break; } }
                if (!de && L0 > 0) { k2_bad++; if (k2_bad <= 3) { printf("K2 COUNTEREXAMPLE n=%d:", n); for (int i = 0; i < n; i++) printf(" %08x", out_[i]); printf("\n"); } }
            }
        }
    }
    printf("n=%d mode=%d samples=%ld nonregular-tested; tied extremes %ld | C7 failures %ld | C9 tested %ld failures %ld | strong %ld, kappa>=2 %ld, no L-descent %ld\n",
           n, mode, samples, tied, c7_bad, c9_n, c9_bad, strong_n, k2_n, k2_bad);
    return 0;
}
