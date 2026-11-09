/* AUDIT G1: decomposability defect by brute force over reversal sets, independent of the module-defect formula.
 *  mode "random": random n-tournaments (n = 5..8): brute-force fewest reversals reaching a non-strong or non-prime
 *                 tournament, compared with min(lambda, formula) and with the formula alone.
 *  mode "codes" : read "n code" lines; for each, report the brute-force decomposability distance d (searching reversal
 *                 sets of size 0,1,2,3 exhaustively) and one witnessing reversal set.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
#define MAXN 16
static int popc(unsigned x) { return __builtin_popcount(x); }
static int n;
static unsigned out[MAXN];
static int pu[200], pv[200], np;

static int decomposable(void) {   /* non-strong or non-prime */
    unsigned full = (1u << n) - 1;
    for (unsigned M = 1; M < full; M++) {
        int sz = popc(M);
        if (sz < 2) {
            /* singleton cut: check strong connectivity separately below */
            continue;
        }
        int ok = 1;
        for (int x = 0; x < n && ok; x++) {
            if (M >> x & 1) continue;
            unsigned o = out[x] & M;
            if (o != 0 && o != M) ok = 0;
        }
        if (ok) return 1;
    }
    /* non-strong check: some nonempty proper S with no arc leaving */
    for (unsigned S = 1; S < full; S++) {
        int leave = 0;
        for (int v = 0; v < n && !leave; v++) if ((S >> v & 1) && (out[v] & ~S & full)) leave = 1;
        if (!leave) return 1;
    }
    return 0;
}
static void flip(int e) {
    int a = pu[e], b = pv[e];
    out[a] ^= 1u << b;
    out[b] ^= 1u << a;
}
static int formula(void) {
    unsigned full = (1u << n) - 1;
    int best = 1 << 30;
    for (unsigned M = 1; M < full; M++) {
        int sz = popc(M);
        if (sz < 2) continue;
        int c = 0;
        for (int x = 0; x < n; x++) {
            if (M >> x & 1) continue;
            int k = popc(out[x] & M);
            c += k < sz - k ? k : sz - k;
        }
        if (c < best) best = c;
    }
    return best;
}
static int lambda_cut(void) {
    unsigned full = (1u << n) - 1;
    int best = 1 << 30;
    for (unsigned S = 1; S < full; S++) {
        int c = 0;
        for (int v = 0; v < n; v++) if (S >> v & 1) c += popc(out[v] & ~S & full);
        if (c < best) best = c;
    }
    return best;
}
/* brute force: smallest d <= dmax such that some d-set of reversals makes T decomposable; -1 if none */
static int wit[4];
static int brute(int dmax) {
    if (decomposable()) return 0;
    for (int a = 0; a < np; a++) {
        flip(a);
        if (decomposable()) { flip(a); wit[0] = a; return 1; }
        flip(a);
    }
    if (dmax < 2) return -1;
    for (int a = 0; a < np; a++) {
        flip(a);
        for (int b = a + 1; b < np; b++) {
            flip(b);
            if (decomposable()) { flip(b); flip(a); wit[0] = a; wit[1] = b; return 2; }
            flip(b);
        }
        flip(a);
    }
    if (dmax < 3) return -1;
    for (int a = 0; a < np; a++) {
        flip(a);
        for (int b = a + 1; b < np; b++) {
            flip(b);
            for (int c = b + 1; c < np; c++) {
                flip(c);
                if (decomposable()) { flip(c); flip(b); flip(a); wit[0] = a; wit[1] = b; wit[2] = c; return 3; }
                flip(c);
            }
            flip(b);
        }
        flip(a);
    }
    return -1;
}
static u64 rng = 0x9E3779B97F4A7C15ULL;
static u64 xr(void) { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return rng; }
static void setpairs(void) {
    np = 0;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { pu[np] = i; pv[np] = j; np++; }
}
int main(int argc, char **argv) {
    if (!strcmp(argv[1], "random")) {
        int bad = 0, tests = 0, hist[8] = {0};
        for (n = 5; n <= 8; n++) {
            setpairs();
            for (int t = 0; t < 400; t++) {
                for (int i = 0; i < n; i++) out[i] = 0;
                for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
                    if (xr() & 1) out[i] |= 1u << j; else out[j] |= 1u << i;
                }
                int f = formula(), l = lambda_cut();
                int fm = f < l ? f : l;
                int b = brute(fm <= 3 ? fm : 3);
                tests++;
                if (fm <= 3 && b != fm) { bad++; printf("MISMATCH n=%d formula=%d lambda=%d brute=%d\n", n, f, l, b); }
                if (f > l) { bad++; printf("moddef > lambda at n=%d\n", n); }
                if (fm < 8) hist[fm]++;
            }
        }
        printf("random: %d tests, %d failures; distribution of dec:", tests, bad);
        for (int i = 0; i < 8; i++) printf(" %d:%d", i, hist[i]);
        printf("\n");
        return 0;
    }
    if (!strcmp(argv[1], "codes")) {
        char line[4096];
        while (fgets(line, sizeof line, stdin)) {
            char code[512];
            if (sscanf(line, "%d %511s", &n, code) != 2) continue;
            setpairs();
            int p = 0;
            for (int i = 0; i < n; i++) out[i] = 0;
            for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
                if (code[p] == '1') out[i] |= 1u << j; else out[j] |= 1u << i;
                p++;
            }
            int f = formula(), l = lambda_cut();
            int b = brute(3);
            printf("n=%d formula=%d lambda=%d brute_dec=%d witness:", n, f, l, b);
            for (int k = 0; k < b; k++) printf(" (%d,%d)", pu[wit[k]], pv[wit[k]]);
            printf("\n");
            fflush(stdout);
        }
        return 0;
    }
    return 1;
}
