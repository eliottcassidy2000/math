/* Audit of HYP-9391 (A): brute-force seven-vertex Paley graph sums.
 *   S_H(q) = sum over injective x : {0..6} -> F_q with x_fa = 0, x_fb = 1 of prod_{(i<j) in H} chi(x_j - x_i)
 * For every injective assignment the 21-bit mask of pairs (i<j) with chi(x_j - x_i) = -1 is histogrammed; then
 * S_H = sum_m hist[m] (-1)^{|m & H|} for each requested H.  Independent of the lane's contraction code; the fixed
 * pair (fa, fb) is a parameter (the lane fixes (0, 1)).  Threads split the value of the first free vertex.
 * Pair index of (i<j): position in the lexicographic list (0,1),(0,2),...,(5,6).
 * Usage: brute_sums q fa fb nthreads graphs.txt     graphs.txt: one line per graph, "k i1 j1 i2 j2 ... ik jk"
 * Output: one line "q S_1 S_2 ... S_n" (graphs in file order).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <pthread.h>

#define NB (1 << 21)
static int q, fa, fb, nth;
static int pidx[7][7];
static unsigned char *neg;           /* neg[d] = 1 iff chi(d) = -1 */
static int freev[5];
typedef struct { int tid; int64_t *hist; } targ;

static void *work(void *vp) {
    targ *A = (targ *)vp;
    int64_t *hist = A->hist;
    int x[7];
    unsigned char *used = calloc(q, 1);
    x[fa] = 0; x[fb] = 1; used[0] = used[1] = 1;
    /* base mask from the fixed pair */
    uint32_t m0 = 0;
    {
        int i = fa < fb ? fa : fb, j = fa < fb ? fb : fa;
        int d = ((x[j] - x[i]) % q + q) % q;
        if (neg[d]) m0 |= 1u << pidx[i][j];
    }
    int assigned[7]; int na = 2; assigned[0] = fa; assigned[1] = fb;
    /* recursive loops unrolled by hand: levels 0..4 */
    for (int t0 = 2; t0 < q; t0++) {
        if ((t0 - 2) % nth != A->tid) continue;
        int v0 = freev[0]; x[v0] = t0; used[t0] = 1;
        uint32_t m1 = m0;
        for (int a = 0; a < 2; a++) { int u = assigned[a]; int i = u < v0 ? u : v0, j = u < v0 ? v0 : u;
            if (neg[((x[j] - x[i]) % q + q) % q]) m1 |= 1u << pidx[i][j]; }
        for (int t1 = 2; t1 < q; t1++) {
            if (used[t1]) continue;
            int v1 = freev[1]; x[v1] = t1; used[t1] = 1;
            uint32_t m2 = m1;
            int pre1[3] = {fa, fb, v0};
            for (int a = 0; a < 3; a++) { int u = pre1[a]; int i = u < v1 ? u : v1, j = u < v1 ? v1 : u;
                if (neg[((x[j] - x[i]) % q + q) % q]) m2 |= 1u << pidx[i][j]; }
            for (int t2 = 2; t2 < q; t2++) {
                if (used[t2]) continue;
                int v2 = freev[2]; x[v2] = t2; used[t2] = 1;
                uint32_t m3 = m2;
                int pre2[4] = {fa, fb, v0, v1};
                for (int a = 0; a < 4; a++) { int u = pre2[a]; int i = u < v2 ? u : v2, j = u < v2 ? v2 : u;
                    if (neg[((x[j] - x[i]) % q + q) % q]) m3 |= 1u << pidx[i][j]; }
                for (int t3 = 2; t3 < q; t3++) {
                    if (used[t3]) continue;
                    int v3 = freev[3]; x[v3] = t3; used[t3] = 1;
                    uint32_t m4 = m3;
                    int pre3[5] = {fa, fb, v0, v1, v2};
                    for (int a = 0; a < 5; a++) { int u = pre3[a]; int i = u < v3 ? u : v3, j = u < v3 ? v3 : u;
                        if (neg[((x[j] - x[i]) % q + q) % q]) m4 |= 1u << pidx[i][j]; }
                    /* last level: precompute per earlier vertex the bit and the sign of (t - x_u) */
                    int v4 = freev[4];
                    int pre4[6] = {fa, fb, v0, v1, v2, v3};
                    int sgn[6], xv[6]; uint32_t bit[6];
                    for (int a = 0; a < 6; a++) { int u = pre4[a]; xv[a] = x[u]; sgn[a] = (u < v4) ? 1 : -1;
                        int i = u < v4 ? u : v4, j = u < v4 ? v4 : u; bit[a] = 1u << pidx[i][j]; }
                    for (int t4 = 2; t4 < q; t4++) {
                        if (used[t4]) continue;
                        uint32_t m5 = m4;
                        for (int a = 0; a < 6; a++) {
                            int d = sgn[a] > 0 ? t4 - xv[a] : xv[a] - t4;
                            d %= q; if (d < 0) d += q;
                            if (neg[d]) m5 |= bit[a];
                        }
                        hist[m5]++;
                    }
                    used[t3] = 0;
                }
                used[t2] = 0;
            }
            used[t1] = 0;
        }
        used[t0] = 0;
    }
    free(used);
    return NULL;
}

int main(int argc, char **argv) {
    if (argc < 6) { fprintf(stderr, "usage\n"); return 1; }
    q = atoi(argv[1]); fa = atoi(argv[2]); fb = atoi(argv[3]); nth = atoi(argv[4]);
    int k = 0;
    for (int i = 0; i < 7; i++) for (int j = i + 1; j < 7; j++) { pidx[i][j] = k; pidx[j][i] = k; k++; }
    neg = calloc(q, 1);
    for (int d = 1; d < q; d++) neg[d] = 1;
    for (long y = 1; y < q; y++) neg[(y * y) % q] = 0;
    int nf = 0; for (int v = 0; v < 7; v++) if (v != fa && v != fb) freev[nf++] = v;
    /* graphs */
    FILE *fh = fopen(argv[5], "r");
    uint32_t gm[256]; int ng = 0, ne;
    while (fscanf(fh, "%d", &ne) == 1) {
        uint32_t m = 0;
        for (int e = 0; e < ne; e++) { int i, j; if (fscanf(fh, "%d %d", &i, &j) != 2) return 2; m |= 1u << pidx[i][j]; }
        gm[ng++] = m;
    }
    fclose(fh);
    pthread_t th[64]; targ ta[64];
    for (int t = 0; t < nth; t++) { ta[t].tid = t; ta[t].hist = calloc(NB, sizeof(int64_t)); pthread_create(&th[t], NULL, work, &ta[t]); }
    for (int t = 0; t < nth; t++) pthread_join(th[t], NULL);
    int64_t *H = ta[0].hist;
    for (int t = 1; t < nth; t++) { for (int m = 0; m < NB; m++) H[m] += ta[t].hist[m]; free(ta[t].hist); }
    printf("%d", q);
    for (int g = 0; g < ng; g++) {
        int64_t s = 0;
        for (int m = 0; m < NB; m++) if (H[m]) s += (__builtin_popcount((uint32_t)m & gm[g]) & 1) ? -H[m] : H[m];
        printf(" %lld", (long long)s);
    }
    printf("\n");
    return 0;
}
