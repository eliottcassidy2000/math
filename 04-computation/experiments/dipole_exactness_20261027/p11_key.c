/* p11_key.c -- canonical key of the Paley tournament P_q (q prime, 3 mod 4) in the format of regular_landing.c */
#define MAXN 32
#include <nauty/nauty.h>
#include <nauty/nautinv.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
int main(int argc, char **argv) {
    int n = atoi(argv[1]); uint32_t q[32] = {0}; int qr[64] = {0};
    for (int x = 1; x < n; x++) qr[(x * x) % n] = 1;
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (i != j && qr[((j - i) % n + n) % n]) q[i] |= 1u << j;
    graph g[MAXN], cg[MAXN]; int lab[MAXN], ptn[MAXN], orbits[MAXN];
    DEFAULTOPTIONS_DIGRAPH(options); statsblk stats; options.getcanon = TRUE;
    EMPTYGRAPH(g, 1, n);
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (q[i] >> j & 1) ADDONEARC(g, i, j, 1);
    densenauty(g, lab, ptn, orbits, &options, &stats, 1, n, cg);
    uint64_t lo = 0, hi = 0; int pos = 0;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (ISELEMENT(GRAPHROW(cg, i, 1), j)) { if (pos < 64) lo |= 1ULL << pos; else hi |= 1ULL << (pos - 64); } pos++; }
    printf("%016llx%016llx\n", (unsigned long long)hi, (unsigned long long)lo);
    return 0;
}
