/* Independent minimal-element scan for integer cycles of the Terras map T(x) = x/2, (3x+1)/2 on Z.
 * For every start s with 1 <= |s| <= N (both signs, signed 128-bit arithmetic, no conjugation trick):
 * iterate T until it returns to s (cycle with min |.| = |s|), or |T^t s| < |s| (s is not the min-|.| element
 * of any cycle), or t exceeds CAP (reported), or |value| >= 2^100 (reported).
 * usage: b2_scan N CAP */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef __int128 i128;
static i128 iabs(i128 v) { return v < 0 ? -v : v; }
static i128 T(i128 v) {
    if (v % 2 == 0) return v / 2;          /* exact for even v of either sign */
    return (3 * v + 1) / 2;                /* 3v+1 even, exact division */
}
int main(int argc, char **argv) {
    int64_t N = atoll(argv[1]); int64_t CAP = atoll(argv[2]);
    int64_t capped = 0, big = 0, maxstop = 0;
    i128 LIM = ((i128)1) << 100;
    for (int sg = -1; sg <= 1; sg += 2) {
        for (int64_t m = 1; m <= N; m++) {
            i128 s = (i128)sg * m, v = s; int64_t t;
            for (t = 1; t <= CAP; t++) {
                v = T(v);
                if (v == s) {
                    /* print cycle: period, odd count, elements */
                    int64_t k = 0; i128 w = s; printf("CYCLE start=%lld period=%lld", (long long)s, (long long)t);
                    for (int64_t u = 0; u < t; u++) { if (w % 2 != 0) k++; w = T(w); }
                    printf(" odd=%lld elems:", (long long)k);
                    w = s; for (int64_t u = 0; u < t && u < 12; u++) { printf(" %lld", (long long)w); w = T(w); }
                    printf("\n");
                    break;
                }
                if (iabs(v) < m) break;
                if (iabs(v) >= LIM) { big++; break; }
            }
            if (t > CAP) capped++;
            else if (t > maxstop) maxstop = t;
        }
    }
    printf("done N=%lld CAP=%lld capped=%lld big=%lld max_steps_to_return_or_descent=%lld\n",
           (long long)N, (long long)CAP, (long long)capped, (long long)big, (long long)maxstop);
    return 0;
}
