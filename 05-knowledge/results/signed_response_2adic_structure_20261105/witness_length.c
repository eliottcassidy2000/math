/* witness_length.c -- kind-pasteur-2026-11-05-S12, HYP-9476 (B), prediction P6.
 *
 * Exhaustive search over all O/Z words w of length <= DEPTH applied to the two-vertex seed t = (D,E,U,d) = (1,0,2,1):
 *   O: (D,E,U) -> (D, 0, 2U),   Z: (D,E,U) -> (D - E + U, 1, 2U - 1)   (d = 1 throughout; Session 198 eq. (1)).
 * For every word ending in Z (strong source, U = n odd >= 3) record the signed response r = aU - bD.
 * For each target |r| <= R of the right parity, keep the minimal word length and the first minimal word found
 * (depth-first, O before Z).  Report maxlen over |r| <= 2^j, the fraction found, and the number of complete copies
 * of the clock word W at the end of each minimal witness.
 *
 * usage: witness_length a b DEPTH R W      (W = clock word, chronological, e.g. ZOOZ for 3/5)
 * output is deterministic text on stdout.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll;
typedef unsigned long long ull;
static int A_, B_, DEPTH;
static ll RR;
static unsigned char *best;
static ull *bestw;

static void dfs(ll U, ll D, int E, int len, ull w) {
    if (E == 1) {
        ll F = (ll)A_ * U - (ll)B_ * D;
        if (F >= -RR && F <= RR) {
            ll i = F + RR;
            if (best[i] > len) { best[i] = (unsigned char)len; bestw[i] = w; }
        }
    }
    if (len == DEPTH) return;
    dfs(2 * U, D, 0, len + 1, w << 1);                 /* O */
    dfs(2 * U - 1, D - E + U, 1, len + 1, (w << 1) | 1); /* Z */
}

int main(int argc, char **argv) {
    if (argc < 6) { fprintf(stderr, "usage: a b DEPTH R W\n"); return 1; }
    A_ = atoi(argv[1]); B_ = atoi(argv[2]); DEPTH = atoi(argv[3]); RR = atoll(argv[4]);
    const char *W = argv[5]; int k = (int)strlen(W);
    if (DEPTH > 60) { fprintf(stderr, "depth too large for 64-bit words\n"); return 1; }
    ull wbits = 0; for (int i = 0; i < k; i++) wbits = (wbits << 1) | (ull)(W[i] == 'Z');
    best = malloc((size_t)(2 * RR + 1)); memset(best, 255, (size_t)(2 * RR + 1));
    bestw = calloc((size_t)(2 * RR + 1), sizeof(ull));
    dfs(2, 1, 0, 0, 0);
    int par = (((A_ - B_) % 2) + 2) % 2;
    printf("slope %d/%d clock %s depth %d R %lld\n", A_, B_, W, DEPTH, RR);
    for (ll w = 2; w <= RR; w *= 2) {
        int mx = 0; ll am = 0, ms = 0, tot = 0;
        for (ll r = -w; r <= w; r++) {
            if ((((r % 2) + 2) % 2) != par) continue;
            tot++;
            int L = best[r + RR];
            if (L == 255) { ms++; continue; }
            if (L > mx) { mx = L; am = r; }
        }
        printf("window %lld targets %lld missing %lld maxlen %d at r=%lld\n", w, tot, ms, mx, am);
    }
    ll hist[64]; memset(hist, 0, sizeof hist);
    for (ll r = -RR; r <= RR; r++) {
        if ((((r % 2) + 2) % 2) != par) continue;
        int L = best[r + RR];
        if (L == 255) continue;
        ull w = bestw[r + RR]; int t = 0, len = L;
        while (len >= k && ((w & ((1ULL << k) - 1)) == wbits)) { t++; w >>= k; len -= k; }
        hist[t < 63 ? t : 63]++;
    }
    printf("trailing-clock histogram over |r| <= %lld:", RR);
    for (int t = 0; t < 64; t++) if (hist[t]) printf(" t%d=%lld", t, hist[t]);
    printf("\n");
    return 0;
}
