/* A2 audit: descent sieve and maximal (backward-branch) sieve on classes n mod 2^K (no 3-adic refinement).
   Prefix state: s, a = #odd steps, r = x_s mod 3^a where x_s = (3^a n + c_s)/2^s.
   Branch certificate at prefix s: backward word with b <= a odd steps, E even steps, valid residues
   (odd backward step needs value = 2 mod 3), and 3^(a-b) 2^(b+E) < 2^s.  Pruning: E < s - a is necessary.
   Usage: msieve KMAX  -> counts of uncertified classes for every K <= KMAX (descent-only and maximal). */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef unsigned __int128 u128;
static uint64_t P3[64];
static int KMAX;
static uint64_t desc_unc[64], max_unc[64];
static uint64_t nodes_bw;
static int S_, A_;
/* backward DFS: current value mod 3^(A_-b) is rr; returns 1 if some extension certifies */
static int bw(int b, int E, uint64_t rr) {
    nodes_bw++;
    /* even step */
    if (E + 1 < S_ - A_) {
        uint64_t m = P3[A_ - b];
        uint64_t r2 = (uint64_t)(((u128)rr * 2) % m);
        /* success check after even step */
        if ((u128)P3[A_ - b] << (b + E + 1) < ((u128)1 << S_)) return 1;
        if (bw(b, E + 1, r2)) return 1;
    }
    /* odd step */
    if (b < A_ && rr % 3 == 2) {
        uint64_t m = P3[A_ - b];
        uint64_t t = (uint64_t)(((u128)rr * 2 + m - 1) % m);   /* (2rr-1) mod 3^(a-b) */
        uint64_t r2 = (t / 3) % P3[A_ - b - 1];
        if ((u128)P3[A_ - b - 1] << (b + 1 + E) < ((u128)1 << S_)) return 1;
        if (bw(b + 1, E, r2)) return 1;
    }
    return 0;
}
static uint64_t inv2mod(uint64_t m) { return (m + 1) / 2; }   /* 2^-1 mod odd m */
/* forward DFS over prefixes */
static void fw(int s, int a, uint64_t r, int desc_dead, int max_dead) {
    /* desc_dead / max_dead: already certified by descent / by maximal sieve at an earlier prefix */
    int dcert = desc_dead, mcert = max_dead;
    if (!dcert && s > 0 && (u128)P3[a] < ((u128)1 << s)) dcert = 1;
    if (dcert) mcert = 1;
    if (!mcert && s > 0) { S_ = s; A_ = a; if (bw(0, 0, r)) mcert = 1; }
    if (!dcert) desc_unc[s]++;
    if (!mcert) max_unc[s]++;
    if (dcert && mcert) return;    /* subtree fully certified for both sieves */
    if (s == KMAX) return;
    /* even child */
    { uint64_t m = P3[a]; uint64_t r2 = (uint64_t)(((u128)r * inv2mod(m)) % m); fw(s + 1, a, r2, dcert, mcert); }
    /* odd child */
    { uint64_t m = P3[a + 1]; uint64_t r2 = (uint64_t)(((u128)(3 * (u128)r + 1) * inv2mod(m)) % m); fw(s + 1, a + 1, r2, dcert, mcert); }
}
int main(int argc, char **argv) {
    KMAX = atoi(argv[1]);
    P3[0] = 1; for (int i = 1; i < 40; i++) P3[i] = P3[i - 1] * 3;
    fw(0, 0, 0, 0, 0);
    printf(" K   descent_uncert   maximal_uncert   ratio   (backward nodes %llu)\n", (unsigned long long)nodes_bw);
    for (int K = 1; K <= KMAX; K++)
        printf("%2d  %14llu  %14llu   %.4f\n", K, (unsigned long long)desc_unc[K], (unsigned long long)max_unc[K], (double)max_unc[K] / desc_unc[K]);
    return 0;
}
