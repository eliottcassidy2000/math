/*
 * aud82_expsum.c  (independent audit of opus-2026-11-06-S82; checks the uncommitted section-7 bullet)
 * |S(theta)|/V with S(theta) = sum_{v in [V, 2V)} e(theta g(v)) at theta = a/q for small q, and the residue
 * frequencies of g(v) mod 7 (times 7) on [2^k, 2^(k+1)).
 * Usage: aud82_expsum KMAX
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
typedef long long ll;
#define SB 16
#define SN (1 << SB)
static ll gs[SN], Ct[SN], Et[SN], Kt[SN];
static ll g_closed(ll v) {
    if (v == 0) return 1;
    int top = 63 - __builtin_clzll((unsigned long long)v);
    ll s = 1;
    for (int i = 0; i < top; i++) if (!((v >> i) & 1)) s += (v >> (i + 1)) + ((v >> (i + 1)) & 1);
    return s;
}
static inline ll G(ll v) {
    if (v < SN) return gs[v];
    ll k = v >> SB, r = v & (SN - 1);
    return G(k) + k * Ct[r] + (k & 1) * Et[r] + Kt[r];
}
int main(int argc, char **argv) {
    int KMAX = argc > 1 ? atoi(argv[1]) : 24;
    for (ll v = 0; v < SN; v++) gs[v] = g_closed(v);
    for (int r = 0; r < SN; r++) {
        ll C = 0, K = 0;
        for (int i = 0; i < SB; i++) if (!((r >> i) & 1)) { C += 1LL << (SB - 1 - i); K += r >> (i + 1); }
        for (int i = 0; i + 1 < SB; i++) if (((r >> (i + 1)) & 1) && !((r >> i) & 1)) K++;
        Ct[r] = C; Et[r] = !((r >> (SB - 1)) & 1); Kt[r] = K;
    }
    int qs[] = {7, 7, 7, 3, 5, 5, 9, 21, 49, 11, 13};
    int as[] = {1, 2, 3, 1, 1, 2, 1, 1, 1, 1, 1};
    int nq = 11;
    printf("# k | |S(a/q)|/V for a/q = 1/7 2/7 3/7 1/3 1/5 2/5 1/9 1/21 1/49 1/11 1/13 | 7 * freq(g mod 7 = 0..6)\n");
    for (int k = 10; k <= KMAX; k++) {
        ll V = 1LL << k;
        static ll cnt[49 * 21];
        ll c7[7] = {0};
        /* residues mod lcm-friendly moduli: compute g mod 9*49*5*11*13*... separately per q */
        double re[16] = {0}, im[16] = {0};
        for (ll v = V; v < 2 * V; v++) {
            ll gv = G(v);
            c7[gv % 7]++;
            for (int j = 0; j < nq; j++) {
                double ang = 2 * M_PI * (double)((gv % qs[j]) * as[j] % qs[j]) / qs[j];
                re[j] += cos(ang); im[j] += sin(ang);
            }
        }
        printf("k=%2d |", k);
        for (int j = 0; j < nq; j++) printf(" %.4f", sqrt(re[j] * re[j] + im[j] * im[j]) / V);
        printf(" |");
        for (int r = 0; r < 7; r++) printf(" %.3f", 7.0 * c7[r] / V);
        printf("\n"); fflush(stdout);
        (void)cnt;
    }
    return 0;
}
