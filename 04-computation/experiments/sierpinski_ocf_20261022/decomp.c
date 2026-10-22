/* Run decompositions of a core of the linear-forest loop expansion: a line-by-line C port of
 * loop_engine.decompositions_fast (THM-4757 / HYP-9367), for excess 3 (cores with up to 6 vertices, 9 edges).
 *
 * A core edge e has a tail half-edge 2e (at ea[e]) and a head half-edge 2e+1 (at eb[e]).  Summed layers:
 *   - at each core vertex a partial matching of its half-edges (pass-throughs), Moebius weight (-1)^(s-1)(s-1)!,
 *     s = #half-edges - #pairs;
 *   - each edge of type m1 (one segment) or m2 (>= 2 segments);
 *   - runs = components of the link graph (m1 path links + pass-through links); closed runs are excluded; runs that
 *     touch no pass-through are free (orientation pre-summed, descriptor '*'), the others get both orientations.
 * Descriptor indices as in loop_engine_fast.DESCS: m1 sign +1, -1, * -> 0, 1, 2; m2 (s_t, s_h) -> 3 + 3*i(s_t) + i(s_h),
 * i(+1, -1, *) = 0, 1, 2.  Output modes:
 *   decomp -k          : one text line "w d_0 ... d_(E-1)" per descriptor tuple with nonzero total weight (validation);
 *   decomp OUTFILE     : the weight tensor in the 5-matrix basis of loop_engine_fast (B0 = R_+, B1 = R_-, B2 = Q_++,
 *                        B3 = Q_+- = Q_-+, B4 = Q_--), dense int64, C order over edges 0..E-1 (= CoreB.Wb), written
 *                        as raw binary; each (configuration, orientation) is expanded through CMAT directly.
 * Input (stdin): V E then E lines "a b".  Build: cc -O2 -o decomp decomp.c.  mac-mini-2026-10-22-loop4. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXE 10
#define MAXH (2 * MAXE)
#define MAXV 8

static int V, E;
static int ea[MAXE], eb[MAXE];
static int nh_at[MAXV], h_at[MAXV][MAXH];

typedef struct { int np; signed char p[MAXH / 2][2]; } Matching;
static Matching *mat[MAXV];
static int nmat[MAXV], capmat[MAXV];

static void add_matching(int v, int np, signed char p[][2]) {
    if (nmat[v] == capmat[v]) { capmat[v] = capmat[v] ? 2 * capmat[v] : 64; mat[v] = realloc(mat[v], capmat[v] * sizeof(Matching)); }
    Matching *m = &mat[v][nmat[v]++];
    m->np = np; memcpy(m->p, p, np * 2);
}

/* all partial matchings of the list items[0..k-1] (same recursion as loop_engine.pairings) */
static void gen_match(int v, int *items, int k, int np, signed char p[][2]) {
    if (k == 0) { add_matching(v, np, p); return; }
    int h = items[0];
    int rest[MAXH];
    for (int i = 1; i < k; i++) rest[i - 1] = items[i];
    gen_match(v, rest, k - 1, np, p);                      /* h unpaired */
    for (int i = 0; i < k - 1; i++) {                       /* h paired with rest[i] */
        int rest2[MAXH], m = 0;
        for (int j = 0; j < k - 1; j++) if (j != i) rest2[m++] = rest[j];
        p[np][0] = h; p[np][1] = rest[i];
        gen_match(v, rest2, m, np + 1, p);
    }
}

static int64_t mob(int s) {
    int64_t f = 1;
    for (int i = 2; i <= s - 1; i++) f *= i;
    return ((s - 1) % 2 == 0) ? f : -f;
}

/* open-addressing hash map uint64 key -> int64 weight */
static uint64_t *hk; static int64_t *hv; static size_t hcap, hcnt;
static void hinit(size_t cap) { hcap = cap; hcnt = 0; hk = calloc(hcap, 8); hv = calloc(hcap, 8); }
static void hadd(uint64_t key, int64_t w);
static void hgrow(void) {
    uint64_t *ok = hk; int64_t *ov = hv; size_t oc = hcap;
    hinit(2 * oc);
    for (size_t i = 0; i < oc; i++) if (ok[i]) hadd(ok[i] - 1, ov[i]);
    free(ok); free(ov);
}
static void hadd(uint64_t key, int64_t w) {
    if (2 * (hcnt + 1) > hcap) hgrow();
    uint64_t k1 = key + 1;
    size_t i = (size_t)((k1 * 0x9E3779B97F4A7C15ULL) >> 20) & (hcap - 1);
    while (hk[i] && hk[i] != k1) i = (i + 1) & (hcap - 1);
    if (!hk[i]) { hk[i] = k1; hcnt++; }
    hv[i] += w;
}

/* CMAT rows (loop_engine_fast.desc_coeffs): descriptor -> list of (basis index, coefficient) */
static const int NB[12] = {1, 1, 2, 1, 1, 2, 1, 1, 2, 2, 2, 3};
static const int BI[12][3] = {{0}, {1}, {0, 1}, {2}, {3}, {2, 3}, {3}, {4}, {3, 4}, {2, 3}, {3, 4}, {2, 3, 4}};
static const int BC[12][3] = {{1}, {1}, {1, 1}, {1}, {1}, {1, 1}, {1}, {1}, {1, 1}, {1, 1}, {1, 1}, {1, 2, 1}};
static int64_t *dense; static int64_t p5[MAXE + 1];
static void expand(const int *d, int e, int64_t idx, int64_t w) {
    if (e == E) { dense[idx] += w; return; }
    for (int j = 0; j < NB[d[e]]; j++) expand(d, e + 1, idx + BI[d[e]][j] * p5[E - 1 - e], w * BC[d[e]][j]);
}

int main(int argc, char **argv) {
    int keymode = (argc > 1 && strcmp(argv[1], "-k") == 0);
    if (!keymode && argc < 2) { fprintf(stderr, "usage: decomp -k | decomp OUTFILE\n"); return 1; }
    if (scanf("%d %d", &V, &E) != 2) return 1;
    if (V > MAXV || E > MAXE) { fprintf(stderr, "too large\n"); return 1; }
    for (int e = 0; e < E; e++) {
        if (scanf("%d %d", &ea[e], &eb[e]) != 2) return 1;
        h_at[ea[e]][nh_at[ea[e]]++] = 2 * e;
        h_at[eb[e]][nh_at[eb[e]]++] = 2 * e + 1;
    }
    signed char p[MAXH / 2][2];
    for (int v = 0; v < V; v++) gen_match(v, h_at[v], nh_at[v], 0, p);
    p5[0] = 1;
    for (int e = 1; e <= E; e++) p5[e] = 5 * p5[e - 1];
    if (keymode) hinit(1 << 16);
    else dense = calloc(p5[E], sizeof(int64_t));
    int H = 2 * E;
    int idx[MAXV]; memset(idx, 0, sizeof idx);
    long long nconf = 0, nvalid = 0;
    for (;;) {
        int64_t mu = 1;
        int link[MAXH];
        for (int h = 0; h < H; h++) link[h] = -1;
        for (int v = 0; v < V; v++) {
            Matching *m = &mat[v][idx[v]];
            mu *= mob(nh_at[v] - m->np);
            for (int j = 0; j < m->np; j++) { link[m->p[j][0]] = m->p[j][1]; link[m->p[j][1]] = m->p[j][0]; }
        }
        for (int types = 0; types < (1 << E); types++) {      /* bit e set: edge e is m1 */
            nconf++;
            int comp[MAXH], ncomp = 0, bad = 0;
            int compfixed[MAXH], compstart[MAXH];
            for (int h = 0; h < H; h++) comp[h] = -1;
            for (int h0 = 0; h0 < H && !bad; h0++) {
                if (comp[h0] >= 0) continue;
                int c = ncomp++;
                int stack[MAXH], sp = 0, fixed = 0, start = -1;
                stack[sp++] = h0; comp[h0] = c;
                while (sp) {
                    int x = stack[--sp];
                    int nb[2], nn = 0;
                    if ((types >> (x >> 1)) & 1) nb[nn++] = x ^ 1;
                    if (link[x] >= 0) { nb[nn++] = link[x]; fixed = 1; }
                    if (nn < 2 && start < 0) start = x;
                    for (int j = 0; j < nn; j++) if (comp[nb[j]] < 0) { comp[nb[j]] = c; stack[sp++] = nb[j]; }
                }
                if (start < 0) { bad = 1; break; }              /* closed run */
                compfixed[c] = fixed; compstart[c] = start;
            }
            if (bad) continue;
            nvalid++;
            /* orient the fixed components: sg[h] = traversal sign of the segment at h for sigma = +1 */
            int sg[MAXH], fidx[MAXH], nf = 0;
            for (int c = 0; c < ncomp; c++) {
                if (!compfixed[c]) { fidx[c] = -1; continue; }
                fidx[c] = nf++;
                int cur = compstart[c];
                int entering = !(link[cur] >= 0 && !((types >> (cur >> 1)) & 1));
                for (;;) {
                    int k = cur & 1, e = cur >> 1;
                    int s = entering ? (k == 0 ? 1 : -1) : (k == 0 ? -1 : 1);
                    sg[cur] = s;
                    if (entering) {
                        if ((types >> e) & 1) {
                            int other = cur ^ 1; sg[other] = s;
                            if (link[other] >= 0) { cur = link[other]; entering = 1; continue; }
                        }
                        break;
                    } else {
                        if (link[cur] >= 0) { cur = link[cur]; entering = 1; continue; }
                        break;
                    }
                }
            }
            for (int sig = 0; sig < (1 << nf); sig++) {
                uint64_t key = 0, pw = 1; int dd[MAXE];
                for (int e = 0; e < E; e++) {
                    int t = 2 * e, hh = 2 * e + 1, d;
                    if ((types >> e) & 1) {
                        if (!compfixed[comp[t]]) d = 2;
                        else { int s = sg[t] * (((sig >> fidx[comp[t]]) & 1) ? -1 : 1); d = (s == 1) ? 0 : 1; }
                    } else {
                        int ia, ib;
                        if (!compfixed[comp[t]]) ia = 2;
                        else { int s = sg[t] * (((sig >> fidx[comp[t]]) & 1) ? -1 : 1); ia = (s == 1) ? 0 : 1; }
                        if (!compfixed[comp[hh]]) ib = 2;
                        else { int s = sg[hh] * (((sig >> fidx[comp[hh]]) & 1) ? -1 : 1); ib = (s == 1) ? 0 : 1; }
                        d = 3 + 3 * ia + ib;
                    }
                    key += (uint64_t)d * pw; pw *= 12; dd[e] = d;
                }
                if (keymode) hadd(key, mu);
                else expand(dd, 0, 0, mu);
            }
        }
        int v = 0;
        while (v < V && ++idx[v] == nmat[v]) { idx[v] = 0; v++; }
        if (v == V) break;
    }
    if (!keymode) {
        FILE *f = fopen(argv[1], "wb");
        if (!f || fwrite(dense, sizeof(int64_t), p5[E], f) != (size_t)p5[E]) { fprintf(stderr, "write failed\n"); return 1; }
        fclose(f);
        long long nzd = 0;
        for (int64_t i = 0; i < p5[E]; i++) nzd += dense[i] != 0;
        fprintf(stderr, "V=%d E=%d: %lld configurations, %lld without closed runs, 5-basis tensor %lld/%lld nonzero\n",
                V, E, nconf, nvalid, nzd, (long long)p5[E]);
        return 0;
    }
    size_t nz = 0;
    for (size_t i = 0; i < hcap; i++) {
        if (!hk[i] || hv[i] == 0) continue;
        nz++;
        uint64_t key = hk[i] - 1;
        printf("%lld", (long long)hv[i]);
        for (int e = 0; e < E; e++) { printf(" %d", (int)(key % 12)); key /= 12; }
        printf("\n");
    }
    fprintf(stderr, "V=%d E=%d: %lld configurations, %lld without closed runs, %zu keys (%zu nonzero)\n", V, E, nconf, nvalid, hcnt, nz);
    return 0;
}
