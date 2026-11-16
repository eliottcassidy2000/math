/* Wave 44, exhaustive form: can ONE added vertex C, with ANY relations to the spine and to the apexes A, B,
   deposit a head-uniform digit while leaving the 7-21 notes/22 width-4 window intact?
   C beats exactly the spine set S (all 2^m subsets) and has free arcs to A and B (4 choices).
   Reference: C a sink (S empty, A -> C, B -> C), which leaves the two-apex count unchanged.
   intact  : H(head, c) - H(head, 0) = 2c for every head word and window setting c
   uniform : D = H(head, 0) - H_ref(head, 0) the same for every head word.
   usage: deposit3 h */
#include <stdio.h>
#include <stdlib.h>
#include "w44_kapex.h"
static const int et[4] = {45, 51, 21, 45}, es[4] = {68, 98, 8, 8}, eo[4] = {0, 1, 0, 1};
static int h;
static ull H(unsigned head, int c, unsigned S, int cA, int cB) {
    int m = h + 7; M_ = m; K_ = 3;
    for (int i = 0; i < m; i++) {
        unsigned u, w, z;
        if (i < h) { u = head >> i & 1; w = 1; }
        else { int k = i - h; u = et[c] >> k & 1; w = !(es[c] >> k & 1); }
        z = (S >> i & 1) ? 0 : 1;                 /* z = 1: spine vertex beats C */
        beat_[i] = u | (w << 1) | (z << 2);
    }
    arc_[0][1] = eo[c]; arc_[1][0] = !eo[c];
    arc_[2][0] = cA; arc_[0][2] = !cA; arc_[2][1] = cB; arc_[1][2] = !cB;
    return kapex_count();
}
int main(int argc, char **argv) {
    h = atoi(argv[1]); int m = h + 7, nh = 1 << h;
    kapex_init_perms(3);
    ull *ref = malloc(sizeof(ull) * nh);
    for (int hd = 0; hd < nh; hd++) ref[hd] = H(hd, 0, 0, 0, 0);
    long tot = 0, intact = 0, uni = 0, nontriv = 0;
    int order[8] = {0, nh - 1, 0x15 & (nh - 1), 0x2a & (nh - 1), 1, 2, 3, 4};
    for (unsigned S = 0; S < (1u << m); S++) for (int cA = 0; cA < 2; cA++) for (int cB = 0; cB < 2; cB++) {
        tot++; int ok = 1;
        for (int t = 0; t < 8 && ok; t++) { int hd = order[t]; ull b = H(hd, 0, S, cA, cB);
            for (int c = 1; c < 4; c++) if (H(hd, c, S, cA, cB) != b + 2 * c) { ok = 0; break; } }
        if (!ok) continue;
        int un = 1; long long D0 = 0;
        for (int hd = 0; hd < nh && ok; hd++) {
            ull b = H(hd, 0, S, cA, cB);
            for (int c = 1; c < 4; c++) if (H(hd, c, S, cA, cB) != b + 2 * c) { ok = 0; break; }
            long long D = (long long)b - (long long)ref[hd];
            if (hd == 0) D0 = D; else if (D != D0) un = 0;
        }
        if (!ok) continue;
        intact++;
        if (!un) continue;
        uni++;
        if (D0 != 0) { nontriv++; if (nontriv <= 20) printf("  NONTRIVIAL: C beats spine set %#x, C->A %d, C->B %d: deposit %lld\n", S, cA, cB, D0); }
        else if (uni <= 6) printf("  trivial: C beats spine set %#x, C->A %d, C->B %d\n", S, cA, cB);
    }
    printf("h=%d (spine %d): %ld one-vertex extensions; window intact for every head: %ld; head-uniform deposit: %ld; nonzero uniform deposits: %ld\n",
           h, m, tot, intact, uni, nontriv);
    return 0;
}
