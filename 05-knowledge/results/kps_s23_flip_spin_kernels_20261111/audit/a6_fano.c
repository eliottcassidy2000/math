/* a6_fano.c (AUDIT, independent) - item 6: all 2^21 tournaments on the 7 points of PG(2,2) (nonzero vectors of Z_2^3).
 * For every tournament: for each of the 7 line complements (the zero-sum 4-sets) count the cyclic triangles
 * directly; the complement is EVEN if the count is 0 or 2, ODD if it is 1.  defect = number of even complements.
 * Output: distribution of the defect; among defect 0: number with all 7 lines cyclic; orbit of the defect-0 set
 * under switchings and line coboundaries; the 128 all-cyclic orientations by (defect, H); the 21 pairs of lines
 * realised by defect-2 patterns and the number of tournaments per pair. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static int pid[8][8];            /* pair index for u < v */
static int arc(unsigned x, int u, int v) { /* 1 iff u -> v ; bit = 1 iff smaller -> larger */
    if (u < v) return (x >> pid[u][v]) & 1; return !((x >> pid[v][u]) & 1); }
static int cyc3(unsigned x, int a, int b, int c) {
    int ab = arc(x, a, b), bc = arc(x, b, c), ca = arc(x, c, a);
    return (ab && bc && ca) || (!ab && !bc && !ca); }
static long long hamp(unsigned x) {
    static long long f[128][7]; memset(f, 0, sizeof f);
    for (int v = 0; v < 7; v++) f[1 << v][v] = 1;
    for (int m = 1; m < 128; m++) for (int v = 0; v < 7; v++) { long long c = f[m][v]; if (!c) continue;
        for (int u = 0; u < 7; u++) if (!(m >> u & 1) && arc(x, v + 1, u + 1)) f[m | 1 << u][u] += c; }
    long long t = 0; for (int v = 0; v < 7; v++) t += f[127][v]; return t; }
int main(void) {
    int k = 0; for (int u = 1; u <= 7; u++) for (int v = u + 1; v <= 7; v++) pid[u][v] = k++;
    int lines[7][3], nl = 0, comp[7][4];
    for (int a = 1; a <= 7; a++) for (int b = a + 1; b <= 7; b++) { int c = a ^ b; if (c > b) { lines[nl][0] = a; lines[nl][1] = b; lines[nl][2] = c; nl++; } }
    for (int l = 0; l < 7; l++) { int n = 0; for (int v = 1; v <= 7; v++) if (v != lines[l][0] && v != lines[l][1] && v != lines[l][2]) comp[l][n++] = v;
        if ((comp[l][0] ^ comp[l][1] ^ comp[l][2] ^ comp[l][3]) != 0) { printf("complement not zero-sum?\n"); return 1; } }
    long long dist[8]; memset(dist, 0, sizeof dist);
    long long d0_allcyc = 0, d0_cyc_hist[8]; memset(d0_cyc_hist, 0, sizeof d0_cyc_hist);
    static unsigned char isd0[1 << 21];
    long long pairs2[7][7]; memset(pairs2, 0, sizeof pairs2);
    long long tricount_hist[7][3]; memset(tricount_hist, 0, sizeof tricount_hist);
    for (unsigned x = 0; x < (1u << 21); x++) {
        int defect = 0, ev[7];
        for (int l = 0; l < 7; l++) { int *q = comp[l];
            int t = cyc3(x, q[0], q[1], q[2]) + cyc3(x, q[0], q[1], q[3]) + cyc3(x, q[0], q[2], q[3]) + cyc3(x, q[1], q[2], q[3]);
            if (t > 2) { printf("impossible\n"); return 1; }
            ev[l] = (t != 1); defect += ev[l]; }
        dist[defect]++;
        if (defect == 0) { isd0[x] = 1; int c = 0; for (int l = 0; l < 7; l++) c += cyc3(x, lines[l][0], lines[l][1], lines[l][2]); d0_cyc_hist[c]++; if (c == 7) d0_allcyc++; }
        if (defect == 2) { int a = -1, b = -1; for (int l = 0; l < 7; l++) if (ev[l]) { if (a < 0) a = l; else b = l; } pairs2[a][b]++; }
    }
    printf("number of even line complements -> number of tournaments:");
    for (int d = 0; d <= 7; d++) printf(" %d:%lld", d, dist[d]); printf("\n");
    printf("defect 0: %lld tournaments; by number of cyclic lines:", dist[0]);
    for (int c = 0; c <= 7; c++) printf(" %d:%lld", c, d0_cyc_hist[c]); printf("\n");
    printf("defect 0 with all 7 lines cyclic: %lld\n", d0_allcyc);
    /* defect-2: every pair of lines, how many tournaments */
    long long mn = 1LL << 60, mx = 0; int npairs = 0;
    for (int a = 0; a < 7; a++) for (int b = a + 1; b < 7; b++) { npairs += pairs2[a][b] > 0; if (pairs2[a][b] < mn) mn = pairs2[a][b]; if (pairs2[a][b] > mx) mx = pairs2[a][b]; }
    printf("defect 2: pairs of lines whose complements are the two even ones: %d of 21 occur, tournaments per pair between %lld and %lld\n", npairs, mn, mx);
    /* orbit of one defect-0 tournament (the quadratic-residue tournament in suitable labels: take the first defect-0
       tournament with all lines cyclic) under the 7 switchings and the 7 line coboundaries */
    unsigned gens[14]; int ng = 0;
    for (int u = 1; u <= 7; u++) { unsigned g = 0; for (int v = 1; v <= 7; v++) if (v != u) g |= 1u << (u < v ? pid[u][v] : pid[v][u]); gens[ng++] = g; }
    for (int f = 1; f <= 7; f++) { unsigned g = 0; for (int v = 1; v <= 7; v++) { int w = v ^ f; if (w > v && w != 0) g |= 1u << pid[v][w]; } gens[ng++] = g; }
    unsigned start = 0; for (unsigned x = 0; x < (1u << 21); x++) if (isd0[x]) { int c = 0; for (int l = 0; l < 7; l++) c += cyc3(x, lines[l][0], lines[l][1], lines[l][2]); if (c == 7) { start = x; break; } }
    static unsigned char seen[1 << 21]; static unsigned stack[1 << 21]; int sp = 0; long long orb = 0, orb_in_d0 = 0;
    seen[start] = 1; stack[sp++] = start;
    while (sp) { unsigned x = stack[--sp]; orb++; orb_in_d0 += isd0[x];
        for (int g = 0; g < ng; g++) { unsigned y = x ^ gens[g]; if (!seen[y]) { seen[y] = 1; stack[sp++] = y; } } }
    printf("orbit of an all-cyclic defect-0 tournament (H = %lld) under 7 switchings + 7 line coboundaries: size %lld, inside the defect-0 set: %lld (defect-0 set has %lld)\n",
           hamp(start), orb, orb_in_d0, dist[0]);
    /* rank of the 14 generators over F_2 */
    { unsigned bas[21]; int r = 0; memset(bas, 0, sizeof bas);
      for (int g = 0; g < ng; g++) { unsigned v = gens[g]; for (int b = 20; b >= 0 && v; b--) { if (!(v >> b & 1)) continue; if (bas[b]) v ^= bas[b]; else { bas[b] = v; r++; v = 0; } } }
      printf("rank of the 14 moves over F_2: %d\n", r); }
    /* the 128 all-cyclic orientations */
    long long hist[8][2]; long long hval[8][2]; memset(hist, 0, sizeof hist); memset(hval, 0, sizeof hval);
    int nkeys = 0; long long keys[16][3];
    for (int bits = 0; bits < 128; bits++) {
        unsigned x = 0;
        for (int l = 0; l < 7; l++) { int a = lines[l][0], b = lines[l][1], c = lines[l][2];
            int t[3]; if (bits >> l & 1) { t[0] = a; t[1] = b; t[2] = c; } else { t[0] = a; t[1] = c; t[2] = b; }
            for (int i = 0; i < 3; i++) { int u = t[i], v = t[(i + 1) % 3]; /* u -> v */ if (u < v) x |= 1u << pid[u][v]; } }
        int defect = 0;
        for (int l = 0; l < 7; l++) { int *q = comp[l];
            int t = cyc3(x, q[0], q[1], q[2]) + cyc3(x, q[0], q[1], q[3]) + cyc3(x, q[0], q[2], q[3]) + cyc3(x, q[1], q[2], q[3]); defect += (t != 1); }
        long long h = hamp(x); int found = 0;
        for (int i = 0; i < nkeys; i++) if (keys[i][0] == defect && keys[i][1] == h) { keys[i][2]++; found = 1; }
        if (!found) { keys[nkeys][0] = defect; keys[nkeys][1] = h; keys[nkeys][2] = 1; nkeys++; }
    }
    printf("the 128 orientations with all 7 lines cyclic, (defect, H) -> count:");
    for (int i = 0; i < nkeys; i++) printf(" (%lld, %lld):%lld", keys[i][0], keys[i][1], keys[i][2]); printf("\n");
    return 0;
}
