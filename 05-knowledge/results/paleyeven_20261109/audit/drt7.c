/* drt7.c - audit P1: every labelled doubly regular tournament on 7 vertices (each pair has exactly one common
 * out-neighbour); count them and compute the minimum module defect over all of them (Lemma A1 at q = 7). */
#include <stdio.h>
#include <stdint.h>
int main(void) {
    int n = 7, pairs[21][2], k = 0;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { pairs[k][0] = i; pairs[k][1] = j; k++; }
    long cnt = 0; int mdmin = 99, mdmax = -1;
    for (uint32_t code = 0; code < (1u << 21); code++) {
        uint32_t out[7] = {0};
        for (int e = 0; e < 21; e++) { int i = pairs[e][0], j = pairs[e][1]; if (code >> e & 1) out[i] |= 1u << j; else out[j] |= 1u << i; }
        int ok = 1;
        for (int i = 0; i < n && ok; i++) for (int j = i + 1; j < n; j++) if (__builtin_popcount(out[i] & out[j]) != 1) { ok = 0; break; }
        if (!ok) continue;
        cnt++;
        int md = 99;
        for (uint32_t M = 1; M < 127; M++) {
            int s = __builtin_popcount(M); if (s < 2) continue;
            int d = 0;
            for (int z = 0; z < n; z++) if (!(M >> z & 1)) { int a = __builtin_popcount(out[z] & M); d += a < s - a ? a : s - a; }
            if (d < md) md = d;
        }
        if (md < mdmin) mdmin = md; if (md > mdmax) mdmax = md;
    }
    printf("labelled doubly regular 7-tournaments: %ld (5040/|Aut(QR7)| = 5040/21 = 240); md range [%d, %d]\n", cnt, mdmin, mdmax);
    return 0;
}
