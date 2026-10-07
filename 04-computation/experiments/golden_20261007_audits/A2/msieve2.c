/* counts per K: descent, Angeltveit fixed-depth (descent + depth-1 path merging + odd-even-even), maximal */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef unsigned __int128 u128;
static uint64_t P3[64]; static int KMAX, S_, A_;
static uint64_t dunc[64], aunc[64], munc[64];
static int bw(int b, int E, uint64_t rr) {
    if (E + 1 < S_ - A_) { uint64_t m = P3[A_ - b]; uint64_t r2 = (uint64_t)(((u128)rr * 2) % m);
        if ((u128)P3[A_ - b] << (b + E + 1) < ((u128)1 << S_)) return 1; if (bw(b, E + 1, r2)) return 1; }
    if (b < A_ && rr % 3 == 2) { uint64_t m = P3[A_ - b]; uint64_t t = (uint64_t)(((u128)rr * 2 + m - 1) % m); uint64_t r2 = (t / 3) % P3[A_ - b - 1];
        if ((u128)P3[A_ - b - 1] << (b + 1 + E) < ((u128)1 << S_)) return 1; if (bw(b + 1, E, r2)) return 1; }
    return 0;
}
static void fw(int s, int a, uint64_t r, int dc, int ac, int mc, int prevbit, int k0, int a0, int zeros) {
    if (!dc && s > 0 && (u128)P3[a] < ((u128)1 << s)) dc = 1;
    if (dc) { ac = 1; mc = 1; }
    if (!ac && s > 0 && a >= 1 && r % 3 == 2 && ((u128)P3[a - 1] << 1) < ((u128)1 << s)) ac = 1;     /* depth-1 */
    if (!ac && zeros == 2 && k0 >= 0 && (u128)P3[a0] < ((u128)1 << (k0 + 1))) ac = 1;               /* odd-even-even */
    if (ac && !mc) { /* consistency: Angeltveit certs are branches; maximal must also certify */ }
    if (!mc && s > 0) { S_ = s; A_ = a; if (bw(0, 0, r)) mc = 1; }
    if (!dc) dunc[s]++; if (!ac) aunc[s]++; if (!mc) munc[s]++;
    if (ac && !mc) { fprintf(stderr, "inconsistency at s=%d\n", s); }
    if (dc && ac && mc) return;
    if (s == KMAX) return;
    { uint64_t m = P3[a]; int z2 = (prevbit == 1) ? 1 : (zeros > 0 ? zeros + 1 : 0);
      fw(s + 1, a, (uint64_t)(((u128)r * ((m + 1) / 2)) % m), dc, ac, mc, 0, k0, a0, z2); }
    { uint64_t m = P3[a + 1]; int nk0 = (prevbit == 1) ? k0 : s, na0 = (prevbit == 1) ? a0 : a;
      fw(s + 1, a + 1, (uint64_t)(((u128)(3 * (u128)r + 1) * ((m + 1) / 2)) % m), dc, ac, mc, 1, nk0, na0, 0); }
}
int main(int argc, char **argv) {
    KMAX = atoi(argv[1]); P3[0] = 1; for (int i = 1; i < 40; i++) P3[i] = P3[i - 1] * 3;
    fw(0, 0, 0, 0, 0, 0, 0, -1, 0, 0);
    printf(" K   descent   Angeltveit(d1+OEE)   maximal   maximal/Angeltveit  maximal/descent\n");
    for (int K = 10; K <= KMAX; K += 2) printf("%2d %10llu %12llu %12llu   %.4f   %.4f\n", K, (unsigned long long)dunc[K], (unsigned long long)aunc[K], (unsigned long long)munc[K], (double)munc[K]/aunc[K], (double)munc[K]/dunc[K]);
    return 0;
}
