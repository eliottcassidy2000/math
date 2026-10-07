// Audit A, items 5/8: exact pair-chain simulator (GMP). State (k, N), N = 3^max(0,-k) e.
// usage: a5_chain_sim k0 N0 npaths CAP seed
// prints survival q(T) on a grid, E[J] (departures from k=0 before absorption) with truncation info,
// the J distribution tail, and the tail of the number of steps of single excursions.
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
static uint64_t s[4];
static inline uint64_t rotl(const uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static uint64_t next(void) {
    const uint64_t result = rotl(s[1] * 5, 7) * 9; const uint64_t t = s[1] << 17;
    s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3]; s[2] ^= t; s[3] = rotl(s[3], 45); return result; }
static uint64_t bitbuf; static int bitcnt = 0;
static inline int rbit(void) { if (!bitcnt) { bitbuf = next(); bitcnt = 64; } int b = bitbuf & 1; bitbuf >>= 1; bitcnt--; return b; }
#define MAXH 200000
static mpz_t P3[MAXH]; static int p3n = 0;
static void ensure(int h) { while (p3n <= h) { mpz_init(P3[p3n]); if (p3n == 0) mpz_set_ui(P3[0], 1); else mpz_mul_ui(P3[p3n], P3[p3n-1], 3); p3n++; } }
static mpz_t N, tmp;
static int k;
static inline void step(int beta) {
    int sig = mpz_odd_p(N);
    if (k >= 1) {
        if (!sig) { if (!beta) mpz_tdiv_q_2exp(N, N, 1); else { ensure(k); mpz_mul_ui(N, N, 3); mpz_add_ui(N, N, 1); mpz_sub(N, N, P3[k]); mpz_tdiv_q_2exp(N, N, 1); } }
        else { if (!beta) { mpz_mul_ui(N, N, 3); mpz_add_ui(N, N, 1); mpz_tdiv_q_2exp(N, N, 1); k++; }
               else { ensure(k); mpz_sub(N, N, P3[k-1]); mpz_tdiv_q_2exp(N, N, 1); k--; } }
    } else if (k == 0) {
        if (!sig) { if (!beta) mpz_tdiv_q_2exp(N, N, 1); else { mpz_mul_ui(N, N, 3); mpz_tdiv_q_2exp(N, N, 1); } }
        else { if (!beta) { mpz_mul_ui(N, N, 3); mpz_add_ui(N, N, 1); mpz_tdiv_q_2exp(N, N, 1); k = 1; }
               else { mpz_mul_ui(N, N, 3); mpz_sub_ui(N, N, 1); mpz_tdiv_q_2exp(N, N, 1); k = -1; } }
    } else {
        int h = -k; ensure(h);
        if (!sig) { if (!beta) mpz_tdiv_q_2exp(N, N, 1); else { mpz_mul_ui(N, N, 3); mpz_add(N, N, P3[h]); mpz_sub_ui(N, N, 1); mpz_tdiv_q_2exp(N, N, 1); } }
        else { if (!beta) { mpz_add(N, N, P3[h-1]); mpz_tdiv_q_2exp(N, N, 1); k++; }
               else { mpz_mul_ui(N, N, 3); mpz_sub_ui(N, N, 1); mpz_tdiv_q_2exp(N, N, 1); k--; } }
    }
}
int main(int argc, char **argv) {
    int k0 = atoi(argv[1]); long N0 = atol(argv[2]); long npaths = atol(argv[3]); long CAP = atol(argv[4]); uint64_t seed = strtoull(argv[5], 0, 10);
    s[0] = seed * 0x9E3779B97F4A7C15ULL + 1; s[1] = seed ^ 0xD1B54A32D192ED03ULL; s[2] = 0x8CB92BA72F3D8DD7ULL; s[3] = seed + 12345;
    for (int i = 0; i < 20; i++) next();
    mpz_init(N); mpz_init(tmp); ensure(64);
    long grid[] = {100, 400, 1600, 6400, 25600, 102400, 409600, 1638400};
    int ng = sizeof(grid) / sizeof(grid[0]);
    long surv[16] = {0};
    double sumJ = 0; double sumJg[16] = {0}; long trunc = 0; long Jhist[64] = {0};
    long exc_total = 0; long exc_tail[40] = {0};   // excursion step-count tail: X > 2^i
    long absorbed_by_cap = 0;
    for (long p = 0; p < npaths; p++) {
        k = k0; mpz_set_si(N, N0);
        long t = 0, J = 0, exc_start = -1;
        int absorbed = 0;
        while (t < CAP) {
            if (k == 0 && mpz_sgn(N) == 0) { absorbed = 1; break; }
            if (k == 0 && mpz_odd_p(N)) { J++; exc_start = t; }
            int kprev = k;
            for (int g = 0; g < ng; g++) if (t == grid[g]) sumJg[g] += J;
            step(rbit()); t++;
            if (k == 0 && kprev != 0 && exc_start >= 0) {   // excursion ended (arrival at 0)
                long X = t - exc_start; exc_total++;
                for (int i = 0; i < 40 && (1L << i) < X; i++) exc_tail[i]++;
            }
        }
        if (absorbed) absorbed_by_cap++; else trunc++;
        for (int g = 0; g < ng; g++) if (t <= grid[g] && absorbed) sumJg[g] += J;
        for (int g = 0; g < ng; g++) if (!absorbed || t > grid[g]) surv[g] += (grid[g] < t || !absorbed) ? 1 : 0;
        sumJ += J; if (J < 63) Jhist[J]++; else Jhist[63]++;
    }
    printf("start (k0,N0)=(%d,%ld): paths %ld, cap %ld steps, absorbed by cap %ld, truncated %ld\n", k0, N0, npaths, CAP, absorbed_by_cap, trunc);
    for (int g = 0; g < ng; g++) if (grid[g] <= CAP) {
        double q = (double)surv[g] / npaths; printf("  T=%8ld  q=%.5f  sqrtT*q=%.3f  (+-%.3f)\n", grid[g], q, sqrt((double)grid[g]) * q, sqrt((double)grid[g]) * sqrt(q * (1 - q) / npaths));
    }
    printf("  E[J] (departures counted up to cap) = %.4f\n", sumJ / npaths);
    for (int g = 0; g < ng; g++) if (grid[g] <= CAP) printf("    E[J counted by T=%ld] = %.4f\n", grid[g], sumJg[g] / npaths);
    printf("  P(J > r): ");
    long cum = npaths; for (int r = 0; r < 40; r++) { cum -= Jhist[r]; if (r % 4 == 0) printf("r=%d:%.4f ", r, (double)cum / npaths); } printf("\n");
    printf("  completed excursions %ld; P(X > 2^i) * 2^(i/2):", exc_total);
    for (int i = 4; i < 40 && (1L << i) < CAP; i += 2) printf(" i=%d:%.3f", i, (double)exc_tail[i] / exc_total * sqrt((double)(1L << i)));
    printf("\n");
    return 0;
}
