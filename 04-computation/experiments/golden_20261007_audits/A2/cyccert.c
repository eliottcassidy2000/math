/* A2 audit: is the 2-adic class of a negative cycle point p (n = p mod 2^K) certified at depth <= K (J = 0)?
   Prefix word = parity word of p (periodic). Same backward search as msieve, in unsigned __int128 with
   division-free halving mod 3^m.  Reports first certification depth or "uncertified up to K", with node counts. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef unsigned __int128 u128;
static u128 P3[90];
static int S_, A_; static long long nodes, budget;
static int over;
static u128 half(u128 r, u128 m) { return (r & 1) ? (r + m) >> 1 : r >> 1; }  /* r * 2^-1 mod m, r < m, m odd */
static int ratio_lt(int ab, int e2) { /* 3^ab < 2^e2 ? */ if (e2 <= 0) return 0; if (e2 >= 127) return 1; return P3[ab] < ((u128)1 << e2); }
static int bw(int b, int E, u128 rr) {
    if (++nodes > budget) { over = 1; return 0; }
    if (E + 1 < S_ - A_) {
        u128 m = P3[A_ - b]; u128 r2 = rr * 2; if (r2 >= m) r2 -= m;
        if (ratio_lt(A_ - b, S_ - b - E - 1)) return 1;
        if (bw(b, E + 1, r2)) return 1;
        if (over) return 0;
    }
    if (b < A_ && rr % 3 == 2) {
        u128 m = P3[A_ - b]; u128 t = rr * 2; if (t >= m) t -= m; t = (t + m - 1); if (t >= m) t -= m;
        u128 r2 = (t / 3) % P3[A_ - b - 1];
        if (ratio_lt(A_ - b - 1, S_ - b - 1 - E)) return 1;
        if (bw(b + 1, E, r2)) return 1;
    }
    return 0;
}
int main(int argc, char **argv) {
    long p = atol(argv[1]); int KMAX = atoi(argv[2]); budget = atoll(argv[3]);
    P3[0] = 1; for (int i = 1; i < 81; i++) P3[i] = P3[i - 1] * 3;
    /* parity word of p under T (p negative) */
    int word[400]; long x = p;
    for (int i = 0; i < KMAX; i++) { word[i] = (int)(((x % 2) + 2) % 2); x = word[i] ? (3 * x + 1) / 2 : x / 2; }
    int a = 0; u128 r = 0; long long tot = 0;
    for (int s = 1; s <= KMAX; s++) {
        int bt = word[s - 1];
        if (bt == 0) { r = half(r, P3[a]); }
        else { if (a + 1 > 80) { printf("a too large at s=%d\n", s); return 0; } u128 m = P3[a + 1]; u128 t = 3 * r + 1; if (t >= m) t -= m; r = half(t, m); a++; }
        /* descent? */
        if (ratio_lt(a, s)) { printf("p=%ld: certified by DESCENT at depth s=%d (a=%d)\n", p, s, a); return 0; }
        S_ = s; A_ = a; nodes = 0; over = 0;
        int c = bw(0, 0, r); tot += nodes;
        if (c) { printf("p=%ld: certified by BRANCH at depth s=%d (a=%d), nodes %lld\n", p, s, a, nodes); return 0; }
        if (over) { printf("p=%ld: budget exceeded at s=%d (a=%d); uncertified (exhaustively) for all s<%d; total nodes %lld\n", p, s, a, s, tot); return 0; }
        if (s % 10 == 0) { fprintf(stderr, "  p=%ld s=%d a=%d nodes=%lld\n", p, s, a, nodes); }
    }
    printf("p=%ld: uncertified (exhaustive) for every depth s<=%d; total backward nodes %lld\n", p, KMAX, tot);
    return 0;
}
