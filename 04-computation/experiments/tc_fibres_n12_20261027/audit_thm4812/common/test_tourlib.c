// test_tourlib.c -- self-tests for tourlib.h
#include "tourlib.h"
static uint64_t rng = 88172645463325252ULL;
static uint64_t xr(void) { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return rng; }
static long long ham_perm(int n, const uint32_t *out) {   // brute force over permutations (n <= 9)
    int p[16]; for (int i = 0; i < n; i++) p[i] = i;
    long long cnt = 0;
    // Heap's algorithm
    int c[16] = {0};
    for (;;) {
        int ok = 1; for (int i = 0; i + 1 < n && ok; i++) if (!((out[p[i]] >> p[i + 1]) & 1u)) ok = 0;
        cnt += ok;
        int i = 1;
        while (i < n && c[i] >= i) { c[i] = 0; i++; }
        if (i >= n) break;
        if (i % 2 == 0) { int t = p[0]; p[0] = p[i]; p[i] = t; } else { int t = p[c[i]]; p[c[i]] = p[i]; p[i] = t; }
        c[i]++;
    }
    return cnt;
}
int main(void) {
    uint64_t *f = malloc(sizeof(uint64_t) * (1u << 12) * 12);
    long long bad = 0, tests = 0;
    for (int n = 1; n <= 9; n++) {
        for (int t = 0; t < 400; t++) {
            uint32_t out[32] = {0};
            for (int u = 0; u < n; u++) for (int v = u + 1; v < n; v++) { if (xr() & 1) out[u] |= 1u << v; else out[v] |= 1u << u; }
            uint32_t in[32]; make_in(n, out, in);
            long long h1 = (long long)ham_pull(n, in, f), h2 = ham_perm(n, out);
            long long a1 = aut_order(n, out), a2 = aut_brute(n, out);
            int s1 = strong_landau(n, out), s2 = strong_bfs(n, out);
            tests++;
            if (h1 != h2 || a1 != a2 || s1 != s2 || h1 % a1) { bad++; printf("FAIL n=%d H %lld/%lld Aut %lld/%lld strong %d/%d\n", n, h1, h2, a1, a2, s1, s2); }
        }
    }
    // structured: circulants and lexicographic products with big groups
    for (int n = 3; n <= 11; n += 2) {
        int K = (n - 1) / 2;
        for (int m = 0; m < (1 << K); m++) {
            uint32_t out[32] = {0};
            for (int u = 0; u < n; u++) for (int d = 1; d <= K; d++) { int s = (m >> (d - 1) & 1) ? n - d : d; out[u] |= 1u << ((u + s) % n); }
            long long a1 = aut_order(n, out), a2 = aut_brute(n, out);
            tests++; if (a1 != a2) { bad++; printf("FAIL circulant n=%d m=%d %lld %lld\n", n, m, a1, a2); }
        }
    }
    // C3[C3] (|Aut| = 81) and C3[C3[C3]]-like sizes via brute force for 9
    {
        int n = 9; uint32_t out[32] = {0};
        for (int u = 0; u < 9; u++) for (int v = 0; v < 9; v++) if (u != v) {
            int bu = u / 3, bv = v / 3;
            if (bu == bv) { if ((v % 3) == (u % 3 + 1) % 3) out[u] |= 1u << v; }
            else if (bv == (bu + 1) % 3) out[u] |= 1u << v;
        }
        long long a1 = aut_order(n, out), a2 = aut_brute(n, out); uint32_t in[32]; make_in(n, out, in);
        printf("C3[C3]: aut_order %lld brute %lld H %llu\n", a1, a2, (unsigned long long)ham_pull(n, in, f));
        tests++; if (a1 != a2 || a1 != 81) bad++;
    }
    printf("tests %lld failures %lld\n", tests, bad);
    return bad != 0;
}
