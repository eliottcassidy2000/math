/* audit_q2mc.c -- independent audit: fast sampling of the signs (-1)^(c_3), (-1)^(c_3+c_5), (-1)^(c_5) of uniformly random
 * labelled tournaments through the score formulas of Theorem P (1) (c_3 = C(q1,2) + floor(n/4), c_3 + c_5 = C(q2,2) + floor(n/4)),
 * which audit_small.c / audit_census.c / audit_rand.c verified against direct cycle counts.  Purpose: an independent
 * high-precision sample at the orders of Run C (11..18) to compare with the closed forms of Theorem P (5).
 * Usage: audit_q2mc n N seed     Output: n N sum(-1)^c3 sum(-1)^c5 sum(-1)^(c3+c5)
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
static uint64_t rs[4];
static uint64_t splitmix(uint64_t *x) { uint64_t z = (*x += 0x9E3779B97F4A7C15ull); z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull; z = (z ^ (z >> 27)) * 0x94D049BB133111EBull; return z ^ (z >> 31); }
static uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static uint64_t rnd(void) { uint64_t r = rotl(rs[1] * 5, 7) * 9, t = rs[1] << 17; rs[2] ^= rs[0]; rs[3] ^= rs[1]; rs[1] ^= rs[2]; rs[0] ^= rs[3]; rs[2] ^= t; rs[3] = rotl(rs[3], 45); return r; }
int main(int argc, char **argv) {
    int n = atoi(argv[1]); long N = atol(argv[2]); uint64_t sx = strtoull(argv[3], NULL, 10);
    for (int i = 0; i < 4; i++) rs[i] = splitmix(&sx);
    long s3 = 0, s5 = 0, s35 = 0; uint32_t out[32];
    for (long it = 0; it < N; it++) {
        for (int i = 0; i < n; i++) out[i] = 0;
        uint64_t bits = 0; int nb = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
            if (!nb) { bits = rnd(); nb = 64; }
            if (bits & 1) out[i] |= 1u << j; else out[j] |= 1u << i;
            bits >>= 1; nb--;
        }
        uint32_t s = 0, t = 0;
        for (int a = 0; a < n; a++) if (__builtin_popcount(out[a]) & 1) s |= 1u << a;
        for (int a = 0; a < n; a++) if (__builtin_popcount(out[a] & s) & 1) t |= 1u << a;
        int q1 = __builtin_popcount(s), q2 = __builtin_popcount(t);
        int p3 = (q1 * (q1 - 1) / 2 + n / 4) & 1, p35 = (q2 * (q2 - 1) / 2 + n / 4) & 1;
        s3 += p3 ? -1 : 1; s35 += p35 ? -1 : 1; s5 += (p3 ^ p35) ? -1 : 1;
    }
    printf("%d %ld %ld %ld %ld\n", n, N, s3, s5, s35);
    return 0;
}
