/* Audit of HYP-9391: Legendre triple products and fourth moments, by direct O(q^2) summation.
 *   L(l) = sum_w phi(w(w-1)(w-l))            (built by direct summation, O(q^2))
 * Sums over s != t, s, t not in {0, 1}:
 *   Y02 = - sum phi(s) phi(t) phi(t-1) L(s)^2 L(t/s)
 *   Y05 = - sum phi(s-1) phi(t) phi(t-1) L(s)^2 L(t/s)
 *   Y09 =   sum phi(s) phi(s-1) phi(t) L(s) L(t) L(t/s)
 *   Y18 =   sum phi(s-1) phi(t-1) L(s) L(t) L(t/s)
 * and the auditor's own closed forms of the cone reductions of H04 and H06 (cone set {4,5,6}, base {0,1,2,3},
 * x0 = 0, x1 = 1, s = x2, t = x3; derived in the audit report, checked against a literal cone sum in Python):
 *   Y04 =   sum phi(s) phi(t) phi(t-1) phi(t-s) L(s)^2 L(t)
 *   Y06 =   sum phi(s) phi(s-1) phi(t) phi(t-1) phi(t-s) L(s)^2 L(t)
 * Moments over lambda != 0, 1 with a_l = -L(l):  M2t = sum phi(l) a^2, M4t = sum phi(l) a^4, M4u = sum a^4, M2u = sum a^2.
 * Output per prime: q Y02 Y04 Y05 Y06 Y09 Y18 M2t M4t M2u M4u
 * Usage: ysums q1 q2 ...
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

typedef __int128 i128;
static void p128(i128 v) {
    char buf[64]; int n = 0; int neg = v < 0; unsigned __int128 u = neg ? -(unsigned __int128)v : (unsigned __int128)v;
    if (u == 0) buf[n++] = '0';
    while (u) { buf[n++] = '0' + (int)(u % 10); u /= 10; }
    if (neg) putchar('-');
    while (n) putchar(buf[--n]);
}

int main(int argc, char **argv) {
    for (int a = 1; a < argc; a++) {
        long q = atol(argv[a]);
        int *phi = calloc(q, sizeof(int));
        for (long x = 1; x < q; x++) phi[x] = -1;
        for (long x = 1; x < q; x++) phi[(x * x) % q] = 1;
        long *L = calloc(q, sizeof(long));
        for (long l = 0; l < q; l++) {
            long s = 0;
            for (long w = 0; w < q; w++) {
                long v = (w * ((w - 1 + q) % q)) % q;
                v = (v * ((w - l + q) % q)) % q;
                s += phi[v];
            }
            L[l] = s;
        }
        long *inv = calloc(q, sizeof(long));
        for (long x = 1; x < q; x++) {      /* inverse by Fermat */
            long r = 1, b = x, e = q - 2;
            while (e) { if (e & 1) r = r * b % q; b = b * b % q; e >>= 1; }
            inv[x] = r;
        }
        i128 Y02 = 0, Y04 = 0, Y05 = 0, Y06 = 0, Y09 = 0, Y18 = 0;
        for (long s = 2; s < q; s++) {
            long Ls = L[s], is = inv[s];
            int ps = phi[s], ps1 = phi[s - 1];
            i128 a02 = 0, a05 = 0, a09 = 0, a18 = 0, a04 = 0, a06 = 0;
            for (long t = 2; t < q; t++) {
                if (t == s) continue;
                long r = t * is % q;
                long Lr = L[r], Lt = L[t];
                int pt = phi[t], pt1 = phi[t - 1], pts = phi[(t - s + q) % q];
                a02 += (long)(pt * pt1) * Lr;          /* times -phi(s) L(s)^2 */
                a05 += (long)(pt * pt1) * Lr;          /* times -phi(s-1) L(s)^2 */
                a09 += (long)pt * Lt * Lr;             /* times phi(s) phi(s-1) L(s) */
                a18 += (long)pt1 * Lt * Lr;            /* times phi(s-1) L(s) */
                a04 += (long)(pt * pt1 * pts) * Lt;    /* times phi(s) L(s)^2 */
                a06 += (long)(pt * pt1 * pts) * Lt;    /* times phi(s) phi(s-1) L(s)^2 */
            }
            Y02 += -(i128)ps * Ls * Ls * a02;
            Y05 += -(i128)ps1 * Ls * Ls * a05;
            Y09 += (i128)ps * ps1 * Ls * a09;
            Y18 += (i128)ps1 * Ls * a18;
            Y04 += (i128)ps * Ls * Ls * a04;
            Y06 += (i128)ps * ps1 * Ls * Ls * a06;
        }
        i128 M2t = 0, M4t = 0, M2u = 0, M4u = 0;
        for (long l = 2; l < q; l++) {
            i128 x = -L[l], x2 = x * x;
            M2t += phi[l] * x2; M4t += phi[l] * x2 * x2; M2u += x2; M4u += x2 * x2;
        }
        printf("%ld ", q);
        i128 vals[10] = {Y02, Y04, Y05, Y06, Y09, Y18, M2t, M4t, M2u, M4u};
        for (int k = 0; k < 10; k++) { p128(vals[k]); putchar(k < 9 ? ' ' : '\n'); }
        fflush(stdout);
        free(phi); free(L); free(inv);
    }
    return 0;
}
