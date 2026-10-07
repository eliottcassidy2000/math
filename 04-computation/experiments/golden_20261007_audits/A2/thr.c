/* thresholds of first-found certificates (descent at first descent; else first branch found) for prefixes up to KMAX */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
typedef unsigned __int128 u128; typedef __int128 i128;
static uint64_t P3[64]; static int KMAX, S_, A_;
static double maxd = 0, maxb = 0; static int argd[3], argb[4]; static long nb = 0, nd = 0;
static i128 Cfound; static int bfound, ifound;
static int bw(int b, int E, uint64_t rr, i128 C) {
    if (E + 1 < S_ - A_) {
        uint64_t m = P3[A_ - b]; uint64_t r2 = (uint64_t)(((u128)rr * 2) % m);
        if ((u128)P3[A_ - b] << (b + E + 1) < ((u128)1 << S_)) { Cfound = 2 * C; bfound = b; ifound = b + E + 1; return 1; }
        if (bw(b, E + 1, r2, 2 * C)) return 1;
    }
    if (b < A_ && rr % 3 == 2) {
        uint64_t m = P3[A_ - b]; uint64_t t = (uint64_t)(((u128)rr * 2 + m - 1) % m); uint64_t r2 = (t / 3) % P3[A_ - b - 1];
        i128 C2 = (2 * C - ((i128)1 << S_)); if (C2 % 3) { fprintf(stderr, "nonint C\n"); exit(1); } C2 /= 3;
        if ((u128)P3[A_ - b - 1] << (b + 1 + E) < ((u128)1 << S_)) { Cfound = C2; bfound = b + 1; ifound = b + 1 + E; return 1; }
        if (bw(b + 1, E, r2, C2)) return 1;
    }
    return 0;
}
static void fw(int s, int a, uint64_t r, u128 c) {
    if (s > 0 && (u128)P3[a] < ((u128)1 << s)) {   /* first descent */
        double th = (double)c / (double)(((u128)1 << s) - P3[a]); nd++;
        if (th > maxd) { maxd = th; argd[0] = s; argd[1] = a; }
        return;
    }
    if (s > 0) { S_ = s; A_ = a; if (bw(0, 0, r, (i128)c)) {
        nb++;
        i128 den = ((i128)1 << s) - (i128)P3[a - bfound] * ((i128)1 << ifound);
        double th = (Cfound > 0) ? (double)Cfound / (double)den : 0;
        if (th > maxb) { maxb = th; argb[0] = s; argb[1] = a; argb[2] = bfound; argb[3] = ifound; }
        return; } }
    if (s == KMAX) return;
    { uint64_t m = P3[a]; fw(s + 1, a, (uint64_t)(((u128)r * ((m + 1) / 2)) % m), c); }
    { uint64_t m = P3[a + 1]; fw(s + 1, a + 1, (uint64_t)(((u128)(3 * (u128)r + 1) * ((m + 1) / 2)) % m), 3 * c + ((u128)1 << s)); }
}
int main(int argc, char **argv) {
    KMAX = atoi(argv[1]); P3[0] = 1; for (int i = 1; i < 40; i++) P3[i] = P3[i - 1] * 3;
    fw(0, 0, 0, 0);
    printf("KMAX=%d first-descent certs %ld, max threshold %.3f (log2 %.3f) at s=%d a=%d\n", KMAX, nd, maxd, log2(maxd), argd[0], argd[1]);
    printf("         branch certs (first found) %ld, max threshold %.3f (log2 %.3f) at s=%d a=%d b=%d i=%d\n", nb, maxb, log2(maxb), argb[0], argb[1], argb[2], argb[3]);
    return 0;
}
