/* words_f2.c - E1 / F2 of the S23 pre-registration.
 * One-apex words: chain c_1 -> ... -> c_L (transitive), apex v, w_i = 1 iff v -> c_i.
 * Inversion pairs (i < j, w_i = 1, w_j = 0), d = j - i:
 *   H    = 1 + 2 sum c(d),              c(1) = 1, c(d) = 2^(d-2)
 *   Pi_2 = 2^(L+1) + sum (3^(d-1) + 1) 2^(L-d)        (n = L + 1 vertices)
 * Irreducible words only (w_1 = 1, w_L = 0).  For every level s with H = 2^s + 1 print the odd parts of Pi_2 seen
 * and flag those dividing 9.  usage: words_f2 Lmax */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
#define MAXS 40
#define MAXODD 4096
static u64 oddseen[MAXS][MAXODD]; static int nodd[MAXS], firstL[MAXS][MAXODD]; static u64 firstw[MAXS][MAXODD]; static int firste[MAXS][MAXODD];
static u64 cnt[MAXS];
int main(int argc, char **argv) {
    int Lmax = atoi(argv[1]);
    u64 p3[64]; p3[0] = 1; for (int i = 1; i < 40; i++) p3[i] = p3[i-1] * 3;
    for (int L = 2; L <= Lmax; L++) {
        u64 lim = 1ull << (L - 2);
        for (u64 mid = 0; mid < lim; mid++) {
            u64 w = (1ull << (L - 1)) | (mid << 1);        /* bit (L-1-i) = letter i (0-based), first 1, last 0 */
            u64 m = 0, pi = 1ull << (L + 1);
            for (int j = 1; j < L; j++) {
                if ((w >> (L - 1 - j)) & 1) continue;
                for (int i = 0; i < j; i++) {
                    if (!((w >> (L - 1 - i)) & 1)) continue;
                    int d = j - i;
                    m += (d == 1) ? 1 : (1ull << (d - 2));
                    pi += (p3[d-1] + 1) << (L - d);
                }
            }
            u64 h = 1 + 2 * m, t = h - 1;
            if (t & (t - 1)) continue;                      /* h - 1 not a power of two */
            int s = __builtin_ctzll(t);
            int e = __builtin_ctzll(pi); u64 od = pi >> e;
            cnt[s]++;
            int k; for (k = 0; k < nodd[s]; k++) if (oddseen[s][k] == od) break;
            if (k == nodd[s] && nodd[s] < MAXODD) { oddseen[s][k] = od; firstL[s][k] = L; firstw[s][k] = w; firste[s][k] = e; nodd[s]++; }
        }
    }
    for (int s = 1; s < MAXS; s++) {
        if (!cnt[s]) continue;
        printf("s=%2d H=2^s+1: %llu words, %d odd parts of Pi_2; dividing 9:", s, cnt[s], nodd[s]);
        int any = 0;
        for (int k = 0; k < nodd[s]; k++) if (9 % oddseen[s][k] == 0) { any = 1;
            printf(" [odd %llu, Pi_2 = odd*2^%d, L=%d, w=", oddseen[s][k], firste[s][k], firstL[s][k]);
            for (int i = 0; i < firstL[s][k]; i++) putchar('0' + (int)((firstw[s][k] >> (firstL[s][k] - 1 - i)) & 1));
            printf("]"); }
        if (!any) printf(" none");
        u64 mn = ~0ull; for (int k = 0; k < nodd[s]; k++) if (oddseen[s][k] < mn) mn = oddseen[s][k];
        printf("   (least odd part %llu)\n", mn);
    }
    return 0;
}
