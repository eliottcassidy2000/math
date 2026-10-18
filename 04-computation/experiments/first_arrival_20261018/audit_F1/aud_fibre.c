/* aud_fibre.c -- audit F1: brute-force checks of THM-4690 statement 2 (preimage description) and of the sibling claim of
   statement 4 (P(x) = P(2x+1) for x = 1 mod 8, step clocks 2+2r and 3+2r, opposite signs).

   (1) every x = 1 (mod 8), x <= X1: P(x) = P(2x+1); clock(x) = 2 + 2r, clock(2x+1) = 3 + 2r with r = #strips of (3x+1)/4;
       a(x) = 2, a(2x+1) = 1; 2x+1 = 3 (mod 16) is a seed.
   (2) every seed s <= S: compute y = P(s) and the step clock c; for y <= Y0 check that (s, c) is predicted by the formula of
       statement 2: s = (4 B^j(y) - 1)/3 with c = 2j + 2 when B^j(y) = 1 (mod 3), or s = (2 B^j(y) - 1)/3 with c = 2j + 1 when
       B^j(y) = 2 (mod 3); then count, for every seed y <= Y0, the predicted preimages <= S and compare with the brute-force
       counts (so the formula lists exactly the preimages). Also checks the pair / singleton structure of each fibre.
   Usage: ./aud_fibre X1 S Y0 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

static inline int is_seed(int64_t s) { return (s & 1) && ((s & 7) != 5); }
static int64_t Pstep(int64_t s, int *a, int *r) { int64_t t = 3 * s + 1; int aa = 0; while (!(t & 1)) { t >>= 1; aa++; } int rr = 0; while ((t & 7) == 5) { t = (t - 1) >> 2; rr++; } *a = aa; *r = rr; return t; }

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: %s X1 S Y0\n", argv[0]); return 1; }
    int64_t X1 = atoll(argv[1]), S = atoll(argv[2]), Y0 = atoll(argv[3]);
    int fails = 0;
    /* (1) siblings */
    int64_t npair = 0, bad1 = 0;
    for (int64_t x = 1; x <= X1; x += 8) {
        int a1, r1, a2, r2; int64_t y1 = Pstep(x, &a1, &r1), y2 = Pstep(2 * x + 1, &a2, &r2);
        int64_t z = (3 * x + 1) / 4; int rz = 0; { int64_t t = z; while ((t & 7) == 5) { t = (t - 1) >> 2; rz++; } }
        int ok = (y1 == y2) && a1 == 2 && a2 == 1 && (a1 + 2 * r1 == 2 + 2 * rz) && (a2 + 2 * r2 == 3 + 2 * rz) && is_seed(2 * x + 1) && ((2 * x + 1) % 16 == 3) && (z & 1);
        npair++; if (!ok) { bad1++; if (bad1 < 5) printf("sibling failure at x=%lld\n", (long long)x); }
    }
    printf("(1) sibling pairs x = 1 mod 8, x <= %lld: %lld checked, %lld failures (common image, clocks 2+2r / 3+2r, a = 2 / 1)\n", (long long)X1, (long long)npair, (long long)bad1);
    if (bad1) fails++;
    /* (2) fibres */
    int64_t *cnt = calloc(Y0 + 1, 8); int64_t hits = 0, bad2 = 0;
    for (int64_t s = 1; s <= S; s += 2) {
        if (!is_seed(s)) continue;
        int a, r; int64_t y = Pstep(s, &a, &r); int c = a + 2 * r;
        if (a < 1 || a > 2) { bad2++; continue; }
        if (y > Y0) continue;
        hits++; cnt[y]++;
        /* predicted: B^r(y) = (3s+1)/2^a, residue rule, clock */
        int64_t b = y; for (int j = 0; j < r; j++) b = 4 * b + 1;
        int res = (int)(b % 3), resrule = (int)((y + r) % 3);
        int okp = (res == resrule) && ((a == 2 && res == 1 && s == (4 * b - 1) / 3 && c == 2 * r + 2) || (a == 1 && res == 2 && s == (2 * b - 1) / 3 && c == 2 * r + 1));
        if (!okp) { bad2++; if (bad2 < 5) printf("fibre failure s=%lld y=%lld a=%d r=%d\n", (long long)s, (long long)y, a, r); }
    }
    /* completeness: count predicted preimages <= S for each seed y <= Y0, and the pair/singleton structure */
    int64_t predtot = 0, bad3 = 0, badstruct = 0, ny = 0;
    for (int64_t y = 1; y <= Y0; y += 2) {
        if (!is_seed(y)) { if (cnt[y]) bad3++; continue; }
        ny++;
        int64_t b = y, pc = 0; int64_t elems[80]; int clk[80], sgn[80], ne = 0;
        for (int j = 0; j < 40; j++) {
            if (b > (INT64_MAX / 8)) break;
            int res = (int)(b % 3);
            if (res == 1) { int64_t s = (4 * b - 1) / 3; if (s <= S) pc++; elems[ne] = s; clk[ne] = 2 * j + 2; sgn[ne] = +1; ne++; }
            if (res == 2) { int64_t s = (2 * b - 1) / 3; if (s <= S) pc++; elems[ne] = s; clk[ne] = 2 * j + 1; sgn[ne] = -1; ne++; }
            b = 4 * b + 1;
        }
        predtot += pc; if (pc != cnt[y]) { bad3++; if (bad3 < 5) printf("count mismatch y=%lld predicted %lld brute %lld\n", (long long)y, (long long)pc, (long long)cnt[y]); }
        /* structure: every a=2 element s (clock 2j+2) is followed by 2s+1 with clock 2j+3; the only unpaired element is (2y-1)/3 with clock 1, iff y = 2 mod 3 */
        int unp = 0;
        for (int e = 0; e < ne; e++) {
            if (sgn[e] == +1) { if (!(e + 1 < ne && elems[e + 1] == 2 * elems[e] + 1 && clk[e + 1] == clk[e] + 1 && sgn[e + 1] == -1) && e + 1 < ne) badstruct++; }
            else { int paired = (e > 0 && sgn[e - 1] == +1 && elems[e] == 2 * elems[e - 1] + 1); if (!paired) { unp++; if (!(clk[e] == 1 && elems[e] == (2 * y - 1) / 3 && y % 3 == 2)) badstruct++; } }
        }
        if ((y % 3 == 2) != (unp == 1) || unp > 1) badstruct++;
    }
    printf("(2) seeds s <= %lld with P(s) <= %lld: %lld brute-force preimages, %lld violate the formula/clock; predicted total %lld; per-y count mismatches %lld; structure violations %lld (over %lld seeds y)\n",
           (long long)S, (long long)Y0, (long long)hits, (long long)bad2, (long long)predtot, (long long)bad3, (long long)badstruct, (long long)ny);
    if (bad2 || bad3 || badstruct || predtot != hits) fails++;
    /* edge case y = 1: fibre elements */
    { printf("    fibre of y = 1 (first elements): "); int64_t b = 1; for (int j = 0; j < 12; j++) { int res = (int)(b % 3); if (res == 1) printf("%lld[c=%d,+] ", (long long)((4 * b - 1) / 3), 2 * j + 2); if (res == 2) printf("%lld[c=%d,-] ", (long long)((2 * b - 1) / 3), 2 * j + 1); b = 4 * b + 1; } printf("\n"); }
    printf("%s\n", fails ? "SOME CHECKS FAILED" : "ALL CHECKS PASSED");
    return fails ? 2 : 0;
}
