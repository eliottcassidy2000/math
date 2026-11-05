/*
 * greedy_bank.c (opus-2026-11-04-S78, exploratory): greedy choice of a small sub-bank Y of rooted triples such that
 * F2(Y) = { g(v) y_h - e(v) y_e + (1+v) y_B : v >= 0, y in Y } (optionally together with F1) covers every odd N in
 * [LO, 2^LOGX]. Prints the chosen triples in order with the number of new targets each covers.
 * Usage: gb LOGX bankfile useF1(0/1) LO
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef long long ll;
typedef struct { ll h, e, B; } Tr;
int main(int argc, char **argv) {
    int LOGX = argc > 1 ? atoi(argv[1]) : 20; ll X = 1LL << LOGX;
    const char *fn = argc > 2 ? argv[2] : "bank7.txt"; int useF1 = argc > 3 ? atoi(argv[3]) : 1; ll LO = argc > 4 ? atoll(argv[4]) : 23;
    FILE *f = fopen(fn, "r"); int cap = 4096, n = 0; Tr *t = malloc(sizeof(Tr) * cap); ll a, b, c;
    while (fscanf(f, "%lld %lld %lld", &a, &b, &c) == 3) { if (n == cap) { cap *= 2; t = realloc(t, sizeof(Tr) * cap); } t[n].h = a; t[n].e = b; t[n].B = c; n++; }
    fclose(f);
    uint32_t *g = malloc(sizeof(uint32_t) * (size_t)(X + 1)); g[0] = 1;
    for (ll v = 1; v <= X; v++) { if (v & 1) g[v] = g[v >> 1]; else { ll u = v >> 1; g[v] = (uint32_t)((ll)g[u] + u + (u & 1)); } }
    char *cov = calloc((size_t)X + 1, 1);
    for (ll N = 0; N < LO; N++) cov[N] = 1;
    for (ll N = LO; N <= X; N++) if (!(N & 1)) cov[N] = 1;
    if (useF1) for (ll u = 0; (1 + u) * (1 + u) <= X + 1; u++) { ll gu = g[u], eu = (u == 0 || !(u & 1));
        for (ll v = u; (1 + u) * (1 + v) <= X; v++) { ll ev = (v == 0 || !(v & 1)); ll N = (1 + u) * (1 + v) + gu * (ll)g[v] - eu * ev; if (N >= 0 && N <= X) cov[N] = 1; } }
    ll left = 0; for (ll N = 1; N <= X; N += 2) if (!cov[N]) left++;
    printf("# %s n=%d X=2^%d useF1=%d LO=%lld: uncovered odd targets before greedy %lld\n", fn, n, LOGX, useF1, LO, left);
    char *used = calloc(n, 1);
    while (left > 0) {
        ll best = -1; int bi = -1;
        for (int i = 0; i < n; i++) { if (used[i]) continue; ll cnt = 0;
            for (ll v = 0; (1 + v) * t[i].B <= X; v++) { ll ev = (v == 0 || !(v & 1)); ll N = (ll)g[v] * t[i].h - ev * t[i].e + (1 + v) * t[i].B; if (N >= 1 && N <= X && !cov[N]) cnt++; }
            if (cnt > best) { best = cnt; bi = i; } }
        if (best <= 0) break;
        used[bi] = 1;
        for (ll v = 0; (1 + v) * t[bi].B <= X; v++) { ll ev = (v == 0 || !(v & 1)); ll N = (ll)g[v] * t[bi].h - ev * t[bi].e + (1 + v) * t[bi].B; if (N >= 1 && N <= X && !cov[N]) { cov[N] = 1; left--; } }
        printf("pick (h,e,B)=(%lld,%lld,%lld) new %lld, left %lld\n", t[bi].h, t[bi].e, t[bi].B, best, left);
        fflush(stdout);
    }
    printf("uncovered:"); int k = 0; for (ll N = 1; N <= X && k < 50; N += 2) if (!cov[N]) { printf(" %lld", N); k++; } printf("\n");
    return 0;
}
