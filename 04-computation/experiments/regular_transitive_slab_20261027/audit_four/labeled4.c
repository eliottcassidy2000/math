/* labeled4.c -- gentourng-free cross-check of THM-4823 (B), (C) on ALL labelled tournaments of order n <= 7.
 * Counts labelled TT4-free, W4-free, L4-free, S4-free and locally transitive tournaments by brute force over 4-sets
 * and compares with closed forms implied by the structure theorems:
 *   LT_lab(n)  = 2^(n-1) (n-1)!                         (Knuth; THM-468 corollary)
 *   W_lab(n)   = n! + sum_{k=3..n} C(n,k) (n-k)! (2^(k-1)(k-1)! - k!)    (Theorem C: L => TT_j, L sinkless LT or K1)
 *   S4free_lab = n! a(n), a(n) = a(n-1) + a(n-3)/3, a(0)=1               (ordered K1/C3 stacks, 2 labelled C3 per block)
 * Usage: ./labeled4 nmax      Build: cc -O3 -o labeled4 labeled4.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
static int n; static uint32_t out_[16];
static long fact(int k) { long r = 1; for (int i = 2; i <= k; i++) r *= i; return r; }
static long binom(int a, int b) { if (b < 0 || b > a) return 0; long r = 1; for (int i = 0; i < b; i++) r = r * (a - i) / (i + 1); return r; }
int main(int argc, char **argv) {
    int nmax = atoi(argv[1]);
    for (n = 4; n <= nmax; n++) {
        int E = n * (n - 1) / 2; int pi[64], pj[64], e = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { pi[e] = i; pj[e] = j; e++; }
        long cTT = 0, cW = 0, cL = 0, cS = 0, cLT = 0;
        for (uint64_t mask = 0; mask < (1ULL << E); mask++) {
            for (int i = 0; i < n; i++) out_[i] = 0;
            for (int b = 0; b < E; b++) { if (mask >> b & 1) out_[pi[b]] |= 1u << pj[b]; else out_[pj[b]] |= 1u << pi[b]; }
            int has[4] = {0,0,0,0};
            for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) for (int c = b + 1; c < n; c++) for (int d = c + 1; d < n; d++) {
                uint32_t m = (1u << a) | (1u << b) | (1u << c) | (1u << d); int v[4] = {a,b,c,d}; int src = 0, snk = 0;
                for (int i = 0; i < 4; i++) { int s = __builtin_popcount(out_[v[i]] & m); if (s == 3) src = 1; if (s == 0) snk = 1; }
                has[(src && snk) ? 0 : src ? 1 : snk ? 2 : 3] = 1;
            }
            if (!has[0]) cTT++; if (!has[1]) cW++; if (!has[2]) cL++; if (!has[3]) cS++; if (!has[1] && !has[2]) cLT++;
        }
        long LTf = (1L << (n - 1)) * fact(n - 1);
        long Wf = fact(n); for (int k = 3; k <= n; k++) Wf += binom(n, k) * fact(n - k) * ((1L << (k - 1)) * fact(k - 1) - fact(k));
        double a[16]; a[0] = 1; for (int k = 1; k <= n; k++) a[k] = a[k - 1] + (k >= 3 ? a[k - 3] / 3.0 : 0);
        double Sf = fact(n) * a[n];
        printf("n=%d labelled=%llu | TT4-free %ld | W4-free %ld (formula %ld) L4-free %ld | LT %ld (2^(n-1)(n-1)! = %ld) | S4-free %ld (formula %.1f)\n",
               n, (unsigned long long)(1ULL << E), cTT, cW, Wf, cL, cLT, LTf, cS, Sf);
    }
    return 0;
}
