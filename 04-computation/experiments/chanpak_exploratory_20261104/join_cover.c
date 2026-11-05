/*
 * join_cover.c (opus-2026-11-04-S78, exploratory/out-of-sample coverage driver)
 * Odd N <= 2^LOGX attained by one-arc joins (session29 / THM-4866 (B)): H = x_h y_h - x_e y_e + x_B y_B,
 * x, y rooted triples (h, e, B). Families:
 *   F0 plain values: h of every bank triple; one-apex values g(v);
 *   F1 one-apex x one-apex: (1+u)(1+v) + g(u) g(v) - e(u) e(v);
 *   F2 one-apex x bank: g(v) y_h - e(v) y_e + (1+v) y_B, y in the bank file;
 *   F3 bank x bank (bank file 1 x bank file 2).
 * one-apex: g(0)=1, g(2v+1)=g(v), g(2v)=g(v)+v+(v&1); triple (g(v), e(v), 1+v), e(v) = [v==0 or v even].
 * Usage: jc LOGX bankA bankB [listmax]
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
typedef long long ll;
typedef struct { ll h, e, B; } Tr;
static Tr *readbank(const char *fn, int *cnt) {
    FILE *f = fopen(fn, "r"); if (!f) { fprintf(stderr, "no %s\n", fn); exit(1); }
    int cap = 1024, n = 0; Tr *t = malloc(sizeof(Tr) * cap); ll a, b, c;
    while (fscanf(f, "%lld %lld %lld", &a, &b, &c) == 3) { if (n == cap) { cap *= 2; t = realloc(t, sizeof(Tr) * cap); } t[n].h = a; t[n].e = b; t[n].B = c; n++; }
    fclose(f); *cnt = n; return t;
}
static uint8_t *bits; static ll X;
static inline void mark(ll N) { if (N >= 1 && N <= X && (N & 1)) bits[N >> 4] |= (uint8_t)(1u << ((N >> 1) & 7)); }
static inline int got(ll N) { return (bits[N >> 4] >> ((N >> 1) & 7)) & 1; }
int main(int argc, char **argv) {
    int LOGX = argc > 1 ? atoi(argv[1]) : 22; X = 1LL << LOGX;
    int na, nb; Tr *A = readbank(argc > 2 ? argv[2] : "bank8.txt", &na), *Bk = readbank(argc > 3 ? argv[3] : "bank7.txt", &nb);
    int listmax = argc > 4 ? atoi(argv[4]) : 60;
    bits = calloc((size_t)(X >> 4) + 2, 1);
    uint32_t *g = malloc(sizeof(uint32_t) * (size_t)(X + 1)); g[0] = 1;
    for (ll v = 1; v <= X; v++) { if (v & 1) g[v] = g[v >> 1]; else { ll u = v >> 1; ll t = (ll)g[u] + u + (u & 1); g[v] = t > 0xFFFFFFFFLL ? 0xFFFFFFFFu : (uint32_t)t; } }
    /* F0 */
    for (int i = 0; i < na; i++) mark(A[i].h);
    for (ll v = 0; v <= X; v++) mark(g[v]);
    ll miss0 = 0; for (ll N = 1; N <= X; N += 2) if (!got(N)) miss0++;
    /* F1 */
    for (ll u = 0; (1 + u) * (1 + u) <= X + 1; u++) {
        ll gu = g[u], eu = (u == 0 || !(u & 1));
        for (ll v = u; (1 + u) * (1 + v) <= X; v++) { ll ev = (v == 0 || !(v & 1)); mark((1 + u) * (1 + v) + gu * (ll)g[v] - eu * ev); }
    }
    ll miss1 = 0; for (ll N = 1; N <= X; N += 2) if (!got(N)) miss1++;
    /* F2: one-apex x bank A */
    for (int i = 0; i < na; i++) {
        ll yh = A[i].h, ye = A[i].e, yB = A[i].B;
        for (ll v = 0; (1 + v) * yB <= X; v++) { ll ev = (v == 0 || !(v & 1)); mark((ll)g[v] * yh - ev * ye + (1 + v) * yB); }
    }
    ll miss2 = 0; for (ll N = 1; N <= X; N += 2) if (!got(N)) miss2++;
    /* F3: bank A x bank B */
    for (int i = 0; i < na; i++) for (int j = 0; j < nb; j++) mark(A[i].h * Bk[j].h - A[i].e * Bk[j].e + A[i].B * Bk[j].B);
    ll miss3 = 0, last = -1; for (ll N = 1; N <= X; N += 2) if (!got(N)) { miss3++; last = N; }
    printf("X=2^%d bankA=%d bankB=%d | odd misses after F0 %lld, +F1 %lld, +F2 %lld, +F3 %lld (largest miss %lld)\n", LOGX, na, nb, miss0, miss1, miss2, miss3, last);
    printf("misses:"); ll c = 0; for (ll N = 1; N <= X && c < listmax; N += 2) if (!got(N)) { printf(" %lld", N); c++; } printf("\n");
    for (int b = 10; b < LOGX; b++) { ll lo = 1LL << b, hi = 1LL << (b + 1), m = 0; for (ll N = lo + 1; N < hi; N += 2) if (!got(N)) m++; printf("[2^%d,2^%d): %lld\n", b, b + 1, m); }
    return 0;
}
