/* Wave 43 (mac-mini-2026-11-16-odometer). Includes the exact two-apex count formula from the 7-21 repository (code/twoapex/apexlib.h, absolute path below). */
/* Wave 43: shared-seam carry odometer on two-apex tapes, using the exact two-apex count formula of 7-21
   (code/twoapex/apexlib.h):  H = 1 + 2(c_A + c_B + c_AB) + 4 d_AB.
   Tape: spine 0..m-1 transitive; u bit i = 1 iff spine vertex i beats A; w bit i = 1 iff vertex i beats B; o = 1 iff A -> B.
   Level 0 (steps 2): the head-uniform width-4 window of 7-21 notes/22 on the last 7 spine sites (A's end word, B's
   out-set, arc).  Head sites (0..h-1): A's word free, B loses to every head vertex (w = 1 there).
   Modes:
     window h        : verify the window is uniform over all 2^h heads (offsets 0, 2, 4, 6)
     levels h L      : head split as [free prefix of length h-L][level-1 region of length L]; for every prefix, list the
                       base counts of all 2^L region words; search for a level-1 digit set {0, 6, 12, 18} (offsets from a
                       common base) that is UNIFORM over all prefixes, then test a level-2 digit set {0,18,36,54} in the
                       prefix that is uniform over the level-1 choices (additivity).
     gaps h          : max gap between consecutive distinct base counts over all 2^h heads, in the central half. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll;
static int m;
#include "/Users/e/Documents/GitHub/7-21/code/twoapex/apexlib.h"
static const int et[4] = {45, 51, 21, 45}, es[4] = {68, 98, 8, 8}, eo[4] = {0, 1, 0, 1};
static ll H(unsigned head, int h, int c) {           /* count for head word (bits 0..h-1) and window config c */
    m = h + 7;
    unsigned u = head | ((unsigned)et[c] << h);
    unsigned wblk = (~(unsigned)es[c]) & 0x7Fu;      /* block vertices not beaten by B beat B */
    unsigned w = ((1u << h) - 1) | (wblk << h);
    return formula(u, w, eo[c]);
}
static int cmpll(const void *a, const void *b) { ll x = *(ll *)a, y = *(ll *)b; return x < y ? -1 : x > y; }
int main(int argc, char **argv) {
    if (argc < 3) return 1;
    int h = atoi(argv[2]);
    if (!strcmp(argv[1], "window")) {
        long bad = 0; ll mn = 1LL << 62, mx = 0;
        for (unsigned hd = 0; hd < (1u << h); hd++) { ll b = H(hd, h, 0);
            for (int c = 1; c < 4; c++) if (H(hd, h, c) != b + 2 * c) bad++;
            if (b < mn) mn = b; if (b > mx) mx = b; }
        printf("window uniform over all %u heads of length %d: %ld violations; base counts %lld..%lld\n", 1u << h, h, bad, mn, mx);
    } else if (!strcmp(argv[1], "gaps")) {
        unsigned N = 1u << h; ll *v = malloc(sizeof(ll) * N);
        for (unsigned hd = 0; hd < N; hd++) v[hd] = H(hd, h, 0);
        qsort(v, N, sizeof(ll), cmpll);
        ll lo = v[N / 4], hi = v[3 * N / 4], mg = 0, prev = -1; long distinct = 0;
        for (unsigned i = 0; i < N; i++) { if (i && v[i] == v[i - 1]) continue; distinct++;
            if (prev >= 0 && v[i] >= lo && v[i] <= hi && v[i] - prev > mg) mg = v[i] - prev; prev = v[i]; }
        printf("heads of length %d: %ld distinct base counts; central half [%lld, %lld]: max gap %lld (gap <= 8 would give every odd count there)\n", h, distinct, lo, hi, mg);
    } else if (!strcmp(argv[1], "levels")) {
        int L = atoi(argv[3]); int P = h - L;               /* prefix length P, region length L at the end of the head */
        int nR = 1 << L, nP = 1 << P;
        /* response of each region word r: offset R(p, r) = H(p,r) - H(p,0); uniform if independent of p */
        ll *off = malloc(sizeof(ll) * nR); char *uni = malloc(nR);
        for (int r = 0; r < nR; r++) { uni[r] = 1; off[r] = 0; }
        for (int p = 0; p < nP; p++) { ll b0 = H((unsigned)p | 0u << P, h, 0);
            for (int r = 0; r < nR; r++) { ll d = H((unsigned)p | ((unsigned)r << P), h, 0) - b0;
                if (p == 0) off[r] = d; else if (d != off[r]) uni[r] = 0; } }
        int nu = 0; for (int r = 0; r < nR; r++) nu += uni[r];
        printf("head %d = prefix %d + region %d: %d of %d region words are prefix-uniform (offsets independent of the prefix)\n", h, P, L, nu, nR);
        /* among uniform region words, find digit sets {x, x+6, x+12, x+18} */
        int found = 0;
        for (int a = 0; a < nR && found < 5; a++) if (uni[a]) for (int b = 0; b < nR; b++) if (uni[b] && off[b] == off[a] + 6)
            for (int c = 0; c < nR; c++) if (uni[c] && off[c] == off[a] + 12) for (int d = 0; d < nR; d++) if (uni[d] && off[d] == off[a] + 18) {
                if (found < 5) printf("  level-1 digits (region words) %d %d %d %d with offsets %lld +6 +12 +18\n", a, b, c, d, off[a]);
                found++; goto next; }
        next:;
        int nus = 0; ll mnoff = 1LL << 62, mxoff = -(1LL << 62);
        for (int r = 0; r < nR; r++) if (uni[r]) { nus++; if (off[r] < mnoff) mnoff = off[r]; if (off[r] > mxoff) mxoff = off[r]; }
        printf("  uniform offsets range %lld..%lld; level-1 digit sets with step 6 found: %s\n", mnoff, mxoff, found ? "yes" : "none");
        /* distinct uniform offsets (sorted) */
        ll *uo = malloc(sizeof(ll) * nR); int k = 0; for (int r = 0; r < nR; r++) if (uni[r]) uo[k++] = off[r];
        qsort(uo, k, sizeof(ll), cmpll); int dk = 0; for (int i = 0; i < k; i++) if (!i || uo[i] != uo[i - 1]) uo[dk++] = uo[i];
        printf("  distinct uniform offsets:"); for (int i = 0; i < dk && i < 40; i++) printf(" %lld", uo[i]); printf("%s\n", dk > 40 ? " ..." : "");
    }
    return 0;
}
