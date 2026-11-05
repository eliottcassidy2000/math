/* aud_f1res.c (S78 audit): residue profile of F1-only misses in (2^(L-2), 2^L], using the auditor's g (as in aud_cover.c).
 * Checks the session's exploratory claim (note lines 210-217): uniform modulo powers of 2; N = 1 mod 3 missed ~3.4x more.
 * Usage: aud_f1res L */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
typedef long long ll;
typedef unsigned long long ull;
static ll g_direct(ull v) {
    if (v == 0) return 1;
    int len = 64 - __builtin_clzll(v);
    ll s = 1 + ((v & 1) == 0);
    for (int i = 0; i + 1 < len; i++) if (!((v >> (i + 1)) & 1) && ((v >> i) & 1)) s++;
    for (int i = 0; i < len; i++) if (!((v >> i) & 1)) s += (ll)(v >> (i + 1));
    return s;
}
static uint32_t GK[1 << 16], ZT[1 << 16], QT[1 << 16]; static uint8_t B15Z[1 << 16];
static inline ll gfast(ull v) {
    ull k = v >> 16, r = v & 0xFFFF;
    if (k == 0) return GK[r];
    return (ll)GK[k] + (ll)k * ZT[r] + QT[r] - ((k & 1) == 0 ? B15Z[r] : 0);
}
int main(int argc, char **argv) {
    int L = argc > 1 ? atoi(argv[1]) : 24; ll X = 1LL << L;
    for (ull r = 0; r < (1u << 16); r++) {
        GK[r] = (uint32_t)g_direct(r);
        ll Z = 0, Q = ((r & 1) == 0), S = 0, c01 = 0;
        for (int i = 0; i < 16; i++) if (!((r >> i) & 1)) { Z += 1LL << (15 - i); S += (ll)(r >> (i + 1)); }
        for (int i = 0; i < 15; i++) if (!((r >> (i + 1)) & 1) && ((r >> i) & 1)) c01++;
        Q += c01 + S; ZT[r] = (uint32_t)Z; QT[r] = (uint32_t)Q; B15Z[r] = !((r >> 15) & 1);
    }
    char *hit = calloc((size_t)X + 1, 1);
    for (ll u = 0; (1 + u) * (1 + u) <= X; u++) {
        ll gu = gfast((ull)u), eu = ((u & 1) == 0);
        for (ll v = u; (1 + u) * (1 + v) <= X; v++) { ll N = (1 + u) * (1 + v) + gu * gfast((ull)v) - eu * ((v & 1) == 0); if (N <= X) hit[N] = 1; }
    }
    ll lo = X / 4, tot = 0;
    for (ll N = lo + 1; N <= X; N += 2) if (!hit[N]) tot++;
    printf("F1 misses in (2^%d, 2^%d]: %lld\n", L - 2, L, tot);
    int mods[] = {3, 4, 8, 9, 32};
    for (int i = 0; i < 5; i++) {
        int m = mods[i]; ll c[64] = {0}, b[64] = {0};
        for (ll N = lo + 1; N <= X; N += 2) { b[N % m]++; if (!hit[N]) c[N % m]++; }
        printf("mod %d:", m); for (int r = 0; r < m; r++) if (b[r]) printf(" %d:%.4f", r, (double)c[r] / b[r]); printf("\n");
    }
    return 0;
}
