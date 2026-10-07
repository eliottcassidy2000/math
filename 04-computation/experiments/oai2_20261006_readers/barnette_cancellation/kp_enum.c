/* Exhaustive certification that every k-subset D of candidate edges misses some pool cycle.
   Input file format (text):
     nE W k Nc M CAP
     Nc candidate edge ids
     M lines, each W hex words (edge bitmask of a pool cycle, prescribed-path edges included)
   Output: number of uncertified k-sets (D hit by every pool cycle), and up to CAP of them.
   Method: nested loops over k-subsets in index order; S_j = pool cycles avoiding the first j
   chosen edges (bitset over cycles); at the last level the uncertified last edges are exactly
   the candidate edges lying in EVERY cycle of S_{k-1} (AND of their masks). */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

static int nE, W, K, Nc, M, CAP;
static int *cand;
static uint64_t *mask;      /* M x W */
static uint64_t **avoid;    /* Nc x MW : cycles NOT containing cand[i] */
static int MW;
static long long uncert = 0, leaves = 0, localU = 0;
static int *eu = 0, *ev = 0; static int filt = 0;
static int chosen[16];
static uint64_t *S[17];

static void rec(int level, int start) {
    if (level == K - 1) {
        /* AND of masks of cycles in S[level] */
        uint64_t andm[8];
        for (int w = 0; w < W; w++) andm[w] = ~0ULL;
        int any = 0;
        for (int b = 0; b < MW; b++) {
            uint64_t x = S[level][b];
            while (x) {
                int t = __builtin_ctzll(x); x &= x - 1;
                int c = b * 64 + t;
                uint64_t *mk = mask + (size_t)c * W;
                int nz = 0;
                for (int w = 0; w < W; w++) { andm[w] &= mk[w]; nz |= (andm[w] != 0); }
                any = 1;
                if (!nz) goto done;
            }
        }
    done:
        for (int i = start; i < Nc; i++) {
            int e = cand[i];
            leaves++;
            if ((andm[e >> 6] >> (e & 63)) & 1ULL) {
                if (filt) {
                    /* local = all K edges share a common endpoint */
                    int c1 = eu[e], c2 = ev[e], l1 = 1, l2 = 1;
                    for (int j = 0; j < level; j++) {
                        int f2 = cand[chosen[j]];
                        if (eu[f2] != c1 && ev[f2] != c1) l1 = 0;
                        if (eu[f2] != c2 && ev[f2] != c2) l2 = 0;
                    }
                    if (l1 || l2) { localU++; continue; }
                }
                uncert++;
                if (uncert <= CAP) {
                    printf("U");
                    for (int j = 0; j < level; j++) printf(" %d", cand[chosen[j]]);
                    printf(" %d\n", e);
                }
            }
        }
        (void)any;
        return;
    }
    for (int i = start; i <= Nc - (K - level); i++) {
        chosen[level] = i;
        uint64_t *src = S[level], *dst = S[level + 1], *av = avoid[i];
        for (int b = 0; b < MW; b++) dst[b] = src[b] & av[b];
        rec(level + 1, i + 1);
        if (uncert > CAP) return;
    }
}

int main(int argc, char **argv) {
    FILE *f = fopen(argv[1], "r");
    if (!f) { perror("open"); return 1; }
    if (fscanf(f, "%d %d %d %d %d %d", &nE, &W, &K, &Nc, &M, &CAP) != 6) return 2;
    cand = malloc(sizeof(int) * Nc);
    for (int i = 0; i < Nc; i++) if (fscanf(f, "%d", &cand[i]) != 1) return 3;
    mask = calloc((size_t)M * W, sizeof(uint64_t));
    for (int c = 0; c < M; c++)
        for (int w = 0; w < W; w++) {
            unsigned long long v;
            if (fscanf(f, "%llx", &v) != 1) return 4;
            mask[(size_t)c * W + w] = v;
        }
    {
        char tag[32];
        if (fscanf(f, "%31s", tag) == 1 && strcmp(tag, "ENDPOINTS") == 0) {
            eu = malloc(sizeof(int) * nE); ev = malloc(sizeof(int) * nE);
            for (int e = 0; e < nE; e++) if (fscanf(f, "%d %d", &eu[e], &ev[e]) != 2) return 5;
            filt = 1;
        }
    }
    fclose(f);
    {   /* restrict cycle masks to candidate edges, so the AND can vanish early */
        uint64_t cm[8] = {0};
        for (int i = 0; i < Nc; i++) cm[cand[i] >> 6] |= 1ULL << (cand[i] & 63);
        for (int c = 0; c < M; c++) for (int w = 0; w < W; w++) mask[(size_t)c * W + w] &= cm[w];
    }
    MW = (M + 63) / 64;
    avoid = malloc(sizeof(uint64_t *) * Nc);
    for (int i = 0; i < Nc; i++) {
        avoid[i] = calloc(MW, sizeof(uint64_t));
        int e = cand[i];
        for (int c = 0; c < M; c++)
            if (!((mask[(size_t)c * W + (e >> 6)] >> (e & 63)) & 1ULL))
                avoid[i][c >> 6] |= 1ULL << (c & 63);
    }
    for (int j = 0; j <= K; j++) S[j] = calloc(MW, sizeof(uint64_t));
    for (int c = 0; c < M; c++) S[0][c >> 6] |= 1ULL << (c & 63);
    rec(0, 0);
    printf("DONE uncertified=%lld leaves=%lld local_uncertified=%lld\n", uncert, leaves, localU);
    return 0;
}
