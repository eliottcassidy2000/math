// f1_miss_residues.c (opus-2026-11-04-S78, exploratory): residue profile of the odd numbers missed by
// F1 = {(1+u)(1+v) + g(u) g(v) - e(u) e(v)} up to 2^LOGX, compared with the uniform profile.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef long long ll;
int main(int argc, char **argv) {
    int LOGX = argc > 1 ? atoi(argv[1]) : 24; ll X = 1LL << LOGX;
    uint32_t *g = malloc(sizeof(uint32_t) * (size_t)(X + 1)); g[0] = 1;
    for (ll v = 1; v <= X; v++) { if (v & 1) g[v] = g[v >> 1]; else { ll u = v >> 1; g[v] = (uint32_t)((ll)g[u] + u + (u & 1)); } }
    char *hit = calloc((size_t)X + 1, 1);
    for (ll u = 0; (1 + u) * (1 + u) <= X + 1; u++) { ll gu = g[u], eu = (u == 0 || !(u & 1));
        for (ll v = u; (1 + u) * (1 + v) <= X; v++) { ll ev = (v == 0 || !(v & 1)); ll N = (1 + u) * (1 + v) + gu * (ll)g[v] - eu * ev; if (N <= X) hit[N] = 1; } }
    int mods[] = {3, 4, 5, 7, 8, 9, 16, 32};
    ll lo = X / 4, tot = 0;
    for (ll N = lo + 1; N <= X; N += 2) if (!hit[N]) tot++;
    printf("F1 misses in (2^%d, 2^%d]: %lld\n", LOGX - 2, LOGX, tot);
    for (int i = 0; i < 8; i++) { int m = mods[i]; ll cnt[64] = {0}, base[64] = {0};
        for (ll N = lo + 1; N <= X; N += 2) { base[N % m]++; if (!hit[N]) cnt[N % m]++; }
        printf("mod %2d:", m); for (int r = 0; r < m; r++) if (base[r]) printf(" %d:%.4f", r, (double)cnt[r] / base[r]); printf("\n"); }
    return 0;
}
