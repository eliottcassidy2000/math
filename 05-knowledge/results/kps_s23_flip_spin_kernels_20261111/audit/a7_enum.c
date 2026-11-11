/* a7_enum.c (AUDIT, independent) - exhaustive enumeration of violation patterns of sign tournaments on PG(3,2).
 * A problem is an affine subspace  w0 + span(b_1..b_k)  of F_2^210 (violation bits, 14 per plane, 2 per 4-set;
 * a 4-set is EVEN iff its 2-bit field is not 00).  For every point of the subspace:
 *   E = number of even 4-sets (all 15 planes);
 *   detail mode: per plane the defect; o = number of planes of defect 0, k2 = number of planes of defect exactly 2.
 * Problem file: one problem per line:  label k w0 b_1 ... b_k   (hex, 210 bits).
 * usage: a7_enum file mode modulus residue      (processes the lines with index = residue mod modulus)
 * output per problem (fast, mode 0):   label k minE  hist: E:count ...
 *        detail (mode 1): additionally  maxk2 (largest k2 over all points), table of (o, k2) -> count and least E. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
typedef struct { u64 w[4]; } V;
static V parse(const char *s) {
    V v; memset(&v, 0, sizeof v); int n = (int)strlen(s);
    for (int i = 0; i < n; i++) { char ch = s[n - 1 - i]; int d = (ch >= '0' && ch <= '9') ? ch - '0' : (ch >= 'a' && ch <= 'f') ? ch - 'a' + 10 : ch - 'A' + 10;
        v.w[i / 16] |= (u64)d << (4 * (i % 16)); }
    return v; }
#define M5 0x5555555555555555ULL
static inline int total(const V *v) {
    return __builtin_popcountll((v->w[0] | v->w[0] >> 1) & M5) + __builtin_popcountll((v->w[1] | v->w[1] >> 1) & M5)
         + __builtin_popcountll((v->w[2] | v->w[2] >> 1) & M5) + __builtin_popcountll((v->w[3] | v->w[3] >> 1) & M5); }
static inline unsigned field(const V *v, int p) { /* 14 bits of plane p */
    int pos = 14 * p, wd = pos >> 6, sh = pos & 63; u64 x = v->w[wd] >> sh; if (sh > 50 && wd < 3) x |= v->w[wd + 1] << (64 - sh); return (unsigned)(x & 0x3FFF); }
int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage\n"); return 1; }
    FILE *fp = fopen(argv[1], "r"); int mode = atoi(argv[2]), modulus = atoi(argv[3]), residue = atoi(argv[4]);
    if (!fp) { fprintf(stderr, "cannot open\n"); return 1; }
    static char line[1 << 16]; int idx = -1;
    while (fgets(line, sizeof line, fp)) {
        idx++; if (idx % modulus != residue) continue;
        char *tok = strtok(line, " \t\r\n"); if (!tok) continue; char label[256]; strncpy(label, tok, 255); label[255] = 0;
        tok = strtok(NULL, " \t\r\n"); int k = atoi(tok);
        tok = strtok(NULL, " \t\r\n"); V w = parse(tok); static V b[64];
        for (int i = 0; i < k; i++) { tok = strtok(NULL, " \t\r\n"); b[i] = parse(tok); }
        static u64 hist[128]; memset(hist, 0, sizeof hist);
        static u64 tab[16][16]; static int tmin[16][16]; memset(tab, 0, sizeof tab); for (int i = 0; i < 16; i++) for (int j = 0; j < 16; j++) tmin[i][j] = 999;
        int maxk2 = 0, minE = 999; u64 npts = (u64)1 << k;
        for (u64 g = 0; g < npts; g++) {
            if (g) { int bit = __builtin_ctzll(g); w.w[0] ^= b[bit].w[0]; w.w[1] ^= b[bit].w[1]; w.w[2] ^= b[bit].w[2]; w.w[3] ^= b[bit].w[3]; }
            int E = total(&w); hist[E]++; if (E < minE) minE = E;
            if (mode) { int o = 0, k2 = 0, E2 = 0;
                for (int p = 0; p < 15; p++) { unsigned f = field(&w, p); int d = __builtin_popcount((f | f >> 1) & 0x1555); E2 += d; if (d == 0) o++; else if (d == 2) k2++; }
                if (E2 != E) { fprintf(stderr, "internal mismatch\n"); return 3; }
                tab[o][k2]++; if (E < tmin[o][k2]) tmin[o][k2] = E; if (k2 > maxk2) maxk2 = k2; }
        }
        printf("%s %d %d hist:", label, k, minE);
        for (int e = 0; e < 128; e++) if (hist[e]) printf(" %d:%llu", e, hist[e]);
        if (mode) { printf(" | maxk2 %d | (o,k2):count:minE", maxk2);
            for (int o = 0; o < 16; o++) for (int j = 0; j < 16; j++) if (tab[o][j]) printf(" (%d,%d):%llu:%d", o, j, tab[o][j], tmin[o][j]); }
        printf("\n"); fflush(stdout);
    }
    return 0;
}
