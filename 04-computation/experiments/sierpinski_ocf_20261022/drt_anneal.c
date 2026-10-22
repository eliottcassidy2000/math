/* Simulated annealing for doubly regular tournaments of order n (n = 3 mod 4): minimise
   C(S) = sum_(u<v) ((S S^T)_uv + 1)^2 over skew +-1 matrices S, by single-arc flips.  Each solution (C = 0) is printed
   as a gentourng-style upper-triangle code ('1' at (i,j), i<j, iff i -> j), readable by hk_count.
   Usage: drt_anneal n seed nsol.  mac-mini-2026-10-22-sierp. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
static int n, S[64][64], P[64][64];
static unsigned long long rs;
static double urand(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return (rs >> 11) * (1.0 / 9007199254740992.0); }
static long cost(void) { long c = 0; for (int u = 0; u < n; u++) for (int v = u + 1; v < n; v++) c += (long)(P[u][v] + 1) * (P[u][v] + 1); return c; }
static void recompute(void) { for (int u = 0; u < n; u++) for (int v = 0; v < n; v++) { int s = 0; for (int w = 0; w < n; w++) s += S[u][w] * S[v][w]; P[u][v] = s; } }
/* delta of flipping arc (a,b): rows a and b change: S[a][b] -> -S[a][b], S[b][a] -> -S[b][a] */
static long delta(int a, int b) {
    long d = 0;
    int sab = S[a][b], sba = S[b][a];
    for (int v = 0; v < n; v++) {
        if (v == a || v == b) continue;
        /* P[a][v] = sum_w S[a][w] S[v][w]; term w = b changes by -2 S[a][b] S[v][b] */
        int na = P[a][v] - 2 * sab * S[v][b];
        int nb = P[b][v] - 2 * sba * S[v][a];
        d += (long)(na + 1) * (na + 1) - (long)(P[a][v] + 1) * (P[a][v] + 1);
        d += (long)(nb + 1) * (nb + 1) - (long)(P[b][v] + 1) * (P[b][v] + 1);
    }
    /* P[a][b] = sum_w S[a][w] S[b][w]; terms w = a (S[a][a]=0) and w = b (S[b][b]=0): unchanged */
    return d;
}
static void apply(int a, int b) {
    int sab = S[a][b], sba = S[b][a];
    for (int v = 0; v < n; v++) {
        if (v == a || v == b) continue;
        P[a][v] -= 2 * sab * S[v][b]; P[v][a] = P[a][v];
        P[b][v] -= 2 * sba * S[v][a]; P[v][b] = P[b][v];
    }
    S[a][b] = -sab; S[b][a] = -sba;
}
int main(int argc, char **argv) {
    n = atoi(argv[1]); rs = 88172645463325252ULL ^ (unsigned long long)atoll(argv[2]) * 2654435761ULL; int nsol = atoi(argv[3]);
    int found = 0;
    while (found < nsol) {
        for (int u = 0; u < n; u++) { S[u][u] = 0; for (int v = u + 1; v < n; v++) { int s = urand() < 0.5 ? 1 : -1; S[u][v] = s; S[v][u] = -s; } }
        recompute();
        long c = cost(); double T = 8.0;
        for (long it = 0; it < 4000000 && c > 0; it++) {
            int a = (int)(urand() * n), b = (int)(urand() * (n - 1)); if (b >= a) b++;
            long d = delta(a, b);
            if (d <= 0 || urand() < exp(-d / T)) { apply(a, b); c += d; }
            if (it % 1000 == 0) { T *= 0.995; if (T < 0.3) T = 0.3; }
        }
        if (c == 0) {
            found++;
            for (int u = 0; u < n; u++) for (int v = u + 1; v < n; v++) putchar(S[u][v] == 1 ? '1' : '0');
            putchar('\n'); fflush(stdout);
        }
    }
    return 0;
}
