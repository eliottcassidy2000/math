/* a9_words.c (AUDIT, independent) - one-apex words, prediction F2.
 * Chain c_1 -> ... -> c_L (transitive), apex v, w_i = 1 iff v -> c_i.
 * Formulas derived for the audit (NOT the inversion-pair formulas of the packet):
 *   a Hamiltonian path is  S (increasing), v, S' (increasing)  with S + S' = [L], (S empty or w_max(S) = 0),
 *   (S' empty or w_min(S') = 1).  Hence
 *     H = [w_1 = 1] + [w_L = 0] + #{a : w_a = 0, w_(a+1) = 1} + sum_{b < a, w_b = 1, w_a = 0} 2^(a-b-1),
 *   and with a third set R (the second path, increasing) and a factor 2 for the order of the two paths
 *     Pi_2 / 2 = 1 + sum_{w_b = 1} 2^(L-b) + sum_{w_a = 0} 2^(a-1) + sum_{a < b, w_a = 0, w_b = 1} 2^(a-1+L-b)
 *                + sum_{b < a, w_b = 1, w_a = 0} 2^(b-1) 3^(a-b-1) 2^(L-a).
 * mode 'check Lmax': print "L word H Pi2" for ALL words of length <= Lmax (for the comparison with direct counts).
 * mode 'f2 Lmax'   : all words with w_1 = 1, w_L = 0, L <= Lmax; per level s with H = 2^s + 1: number of words,
 *                    least odd part of Pi_2, words whose odd part divides 9.  (128-bit arithmetic.) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
typedef unsigned __int128 u128;
static u128 p3[64];
static void eval(int L, u64 w, u128 *H, u128 *P2) { /* letter i (1-based) = bit (i-1) of w */
    u128 h = 0, q = 1;
#define W(i) ((int)((w >> ((i) - 1)) & 1))
    if (W(1) == 1) h += 1; if (W(L) == 0) h += 1;
    for (int a = 1; a < L; a++) if (W(a) == 0 && W(a + 1) == 1) h += 1;
    for (int b = 1; b <= L; b++) { if (W(b)) q += (u128)1 << (L - b); else q += (u128)1 << (b - 1); }
    for (int a = 1; a <= L; a++) { if (W(a)) continue;
        for (int b = 1; b <= L; b++) { if (!W(b)) continue;
            if (b < a) { h += (u128)1 << (a - b - 1); q += (((u128)1 << (b - 1)) * p3[a - b - 1]) << (L - a); }
            else q += (u128)1 << (a - 1 + L - b); } }
    *H = h; *P2 = 2 * q; }
static void pr(u128 x) { char buf[64]; int p = 63; buf[p] = 0; if (!x) buf[--p] = '0'; while (x) { buf[--p] = '0' + (int)(x % 10); x /= 10; } printf("%s", buf + p); }
int main(int argc, char **argv) {
    p3[0] = 1; for (int i = 1; i < 64; i++) p3[i] = p3[i - 1] * 3;
    int Lmax = atoi(argv[2]);
    if (!strcmp(argv[1], "check")) {
        for (int L = 1; L <= Lmax; L++) for (u64 w = 0; w < ((u64)1 << L); w++) { u128 H, P; eval(L, w, &H, &P); printf("%d %llu ", L, w); pr(H); printf(" "); pr(P); printf("\n"); }
        return 0; }
    static u64 cnt[64]; static u128 least[64]; static int hits[64];
    for (int s = 0; s < 64; s++) least[s] = ~(u128)0;
    for (int L = 2; L <= Lmax; L++) for (u64 mid = 0; mid < ((u64)1 << (L - 2)); mid++) {
        u64 w = 1 | (mid << 1);                 /* w_1 = 1, w_L = 0 */
        u128 H, P; eval(L, w, &H, &P); u128 t = H - 1; if (t == 0 || (t & (t - 1))) continue;
        int s = 0; while (!((t >> s) & 1)) s++;
        u128 od = P; while (!(od & 1)) od >>= 1;
        cnt[s]++; if (od < least[s]) least[s] = od;
        if (od == 1 || od == 3 || od == 9) { hits[s]++; printf("   level %d: word of length %d (bits, letter 1 first): ", s, L); for (int i = 1; i <= L; i++) putchar('0' + W(i)); printf(" H = "); pr(H); printf(" Pi_2 = "); pr(P); printf("\n"); } }
    for (int s = 1; s < 64; s++) if (cnt[s]) { printf("s = %2d: %llu irreducible words with H = 2^s + 1, least odd part of Pi_2 = ", s, cnt[s]); pr(least[s]); printf(", words with odd part dividing 9: %d\n", hits[s]); }
    return 0;
}
