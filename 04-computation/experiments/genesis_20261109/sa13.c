/* sa13.c (derived from sa12.c) - order-13 search, min semidegree >= 4, T proxy H <= max H(12) = 545835; deletions exact in S(11); hits certified afterwards by wit.c of THM-4942's five-operation grammar
 * (genesis lane, mac-mini-2026-11-09).  EXPLORATORY: finds examples, proves nothing about completeness.
 *
 * Escaper: 12-tournament T, min semidegree >= 3 (equivalent to lambda >= 3 at n = 12 <= 4*3+1, THM-4811 (B)),
 * prime, H(T) in S(11), H(T - v) in S(10) for every v.
 * Energy: sum over v of pen10(H(T-v)) + pen11(H(T)), pen = 0 inside the exact spectrum, else 1 + excess/1000.
 * Moves: single arc reversals preserving min semidegree >= 3.  Optional constraint (argv[3] = 1): every pair of
 * vertices is distinguished by at least 2 others (twin distance >= 2, i.e. T is not one reversal from a
 * substitution with a 2-vertex module).
 * Output: "E code H twin lam" for every escaper found (prime checked; lambda by cut enumeration).
 * usage: sa12 seed steps twinmin  (reads unattained_10.txt, unattained_11.txt from cwd)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#define N 13
static uint32_t out[N];
static unsigned char hole10[15746], hole11[95096];
#define MINSEMI 4
static uint64_t rng = 88172645463325252ULL;
static inline uint64_t xr(void) { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return rng; }

static uint64_t hsub(uint32_t verts) {
    int idx[N], k = 0;
    for (int v = 0; v < N; v++) if (verts >> v & 1) idx[k++] = v;
    static uint64_t dp[1 << N][N];
    uint32_t lo[N];
    for (int a = 0; a < k; a++) { lo[a] = 0; for (int b = 0; b < k; b++) if (out[idx[a]] >> idx[b] & 1) lo[a] |= 1u << b; }
    uint32_t S = 1u << k;
    for (uint32_t m = 0; m < S; m++) for (int a = 0; a < k; a++) dp[m][a] = 0;
    for (int a = 0; a < k; a++) dp[1u << a][a] = 1;
    for (uint32_t m = 1; m < S; m++) for (int a = 0; a < k; a++) {
        uint64_t c = dp[m][a]; if (!c) continue;
        uint32_t nx = lo[a] & ~m;
        while (nx) { int b = __builtin_ctz(nx); nx &= nx - 1; dp[m | (1u << b)][b] += c; }
    }
    uint64_t h = 0; for (int a = 0; a < k; a++) h += dp[S - 1][a];
    return h;
}
static int inS10(uint64_t h) { return h % 2 == 1 && h <= 15745 && h != 7 && h != 21 && !hole10[h]; }
static int inS11(uint64_t h) { return h % 2 == 1 && h <= 95095 && h != 7 && h != 21 && !hole11[h]; }
static double pen10(uint64_t h) { return inS11(h) ? 0.0 : 1.0 + (h > 80451 ? (h - 80451) / 1000.0 : 0.0); }
static double pen11(uint64_t h) { return h <= 545835 ? 0.0 : 1.0 + (h - 545835) / 10000.0; }
static uint64_t Hfull, Hdel[N];
static double energy(void) {
    uint32_t full = (1u << N) - 1;
    Hfull = hsub(full);
    double e = pen11(Hfull);
    for (int v = 0; v < N; v++) { Hdel[v] = hsub(full & ~(1u << v)); e += pen10(Hdel[v]); }
    return e;
}
static int minsemi(void) {
    int m = 99;
    for (int v = 0; v < N; v++) { int o = __builtin_popcount(out[v]); int i = N - 1 - o; if (o < m) m = o; if (i < m) m = i; }
    return m;
}
static int twin(void) {
    int best = 99;
    for (int u = 0; u < N; u++) for (int w = u + 1; w < N; w++) {
        int d = 0;
        for (int x = 0; x < N; x++) if (x != u && x != w && (((out[x] >> u) ^ (out[x] >> w)) & 1)) d++;
        if (d < best) best = d;
    }
    return best;
}
static int prime(void) {
    uint32_t full = (1u << N) - 1;
    for (int u = 0; u < N; u++) for (int w = u + 1; w < N; w++) {
        uint32_t M = (1u << u) | (1u << w); int ch = 1;
        while (ch) { ch = 0; for (int x = 0; x < N; x++) { if (M >> x & 1) continue; uint32_t o = out[x] & M; if (o && o != M) { M |= 1u << x; ch = 1; } } }
        if (M != full) return 0;
    }
    return 1;
}
static int lam(void) {
    uint32_t full = (1u << N) - 1; int best = 99;
    for (uint32_t S = 1; S < full; S++) { int c = 0; uint32_t r = S; while (r) { int v = __builtin_ctz(r); r &= r - 1; c += __builtin_popcount(out[v] & ~S & full); } if (c < best) best = c; }
    return best;
}
static void flip(int a, int b) { out[a] ^= 1u << b; out[b] ^= 1u << a; }
static void load(const char *fn, unsigned char *hole, int max) {
    FILE *f = fopen(fn, "r"); char line[128];
    if (!f) { fprintf(stderr, "cannot open %s\n", fn); exit(1); }
    while (fgets(line, sizeof line, f)) { if (line[0] == '#') continue; long h = atol(line); if (h > 0 && h <= max) hole[h] = 1; }
    fclose(f);
}
int main(int argc, char **argv) {
    rng ^= (uint64_t)atoll(argv[1]) * 0x9E3779B97F4A7C15ULL; for (int i = 0; i < 20; i++) xr();
    long steps = atol(argv[2]); int twinmin = argc > 3 ? atoi(argv[3]) : 0;
    load("unattained_10.txt", hole10, 15745); load("unattained_11.txt", hole11, 95095);
    long found = 0;
    for (long restart = 0; ; restart++) {
        /* random start with min semidegree >= 3 (and twin >= twinmin) */
        do {
            memset(out, 0, sizeof out);
            for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++) { if (xr() & 1) out[i] |= 1u << j; else out[j] |= 1u << i; }
        } while (minsemi() < MINSEMI || twin() < twinmin);
        double e = energy(), T0 = 2.0;
        long per = 20000;
        for (long t = 0; t < per && steps > 0; t++, steps--) {
            double temp = T0 * (1.0 - (double)t / per) + 0.02;
            int a = xr() % N, b = xr() % N; if (a == b) continue;
            flip(a, b);
            if (minsemi() < MINSEMI || (twinmin && twin() < twinmin)) { flip(a, b); continue; }
            uint64_t oh = Hfull, od[N]; memcpy(od, Hdel, sizeof od);
            double e2 = energy();
            if (e2 <= e || exp((e - e2) / temp) > (double)(xr() % 1000000) / 1e6) {
                e = e2;
                if (e == 0.0 && prime()) {
                    int L = lam();
                    if (L >= MINSEMI) {
                        printf("E ");
                        for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++) putchar(out[i] >> j & 1 ? '1' : '0');
                        printf(" %llu %d %d\n", (unsigned long long)Hfull, twin(), L);
                        fflush(stdout); found++;
                        break; /* restart after a find */
                    }
                }
            } else { flip(a, b); Hfull = oh; memcpy(Hdel, od, sizeof od); }
        }
        if (steps <= 0) break;
    }
    fprintf(stderr, "done, found %ld\n", found);
    return 0;
}
