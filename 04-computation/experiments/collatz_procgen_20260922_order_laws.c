/*
 * collatz_procgen_20260922_order_laws.c
 *
 * Census engine for the procedural search for sign-specific ORDER laws of
 * odd-only Syracuse orbits U_b(x) = oddpart(3x + b), b in {+1,-1}, on the
 * positive odd integers.  Driven by collatz_procgen_20260922_order_laws.py.
 *
 * usage: order_laws b N outprefix [lagmax]
 *
 * For every odd n in [1, N] the STOPPED orbit x_0 = n, ..., x_T is followed
 * until the first element of the cycle set C_b:
 *     C_+ = {1},  C_- = {1, 5, 7, 17, 25, 37, 41, 55, 61, 91}.
 * (Values in a stopped orbit are pairwise distinct.)  The word is
 * k_i = v_2(3 x_i + b), K_j = k_0 + ... + k_{j-1}.
 *
 * Outputs (text unless noted):
 *   <prefix>_patterns.txt   ordinal-pattern counts, windows of w = 3..6
 *                           values, four scopes (all / initial / interior /
 *                           large = all window values >= 1000), with first
 *                           and last witness; plus actual-vs-generic counts.
 *   <prefix>_lags.txt       per lag l <= lagmax: 2x2 tables (generic sign of
 *                           3^l vs 2^K') x (actual sign of x_{i+l} - x_i),
 *                           from the start (i = 0) and anywhere (all i).
 *   <prefix>_devs.txt       every distinct deviation window (x_i, l) with
 *                           x_{i+l}, K', from-start flag, terminal flag, and
 *                           whether the start x_i violates sigma = tau.
 *   <prefix>_margins.txt    near-gate windows: positive-gate windows (plus:
 *                           decay words, minus: growth words) with margin
 *                           mu = x_i/gate - 1 < 1, deduplicated by (x_i, l).
 *   <prefix>_records.txt    running-max / running-min record histograms.
 *   <prefix>_cst.txt        sigma/tau disagreements at any orbit point.
 *   <prefix>_pern.bin       per-n summary (binary, struct pern_t below).
 *
 * All comparisons 2^K vs 3^l are decided by K > l*log2(3) in double
 * precision; for l <= 4096 the smallest |K - l log2 3| is > 6e-5, far above
 * rounding error.  Overflow of 3x+b aborts the run.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

typedef unsigned long long u64;
#define MAXT 8192
#define EXTRA 16            /* cycle continuation steps for sigma / tau */
#define LOG23 1.5849625007211561815

static int B;               /* +1 or -1 */
static int LAGMAX = 64;

static int is_cycle(u64 x) {
    if (x == 1) return 1;
    if (B > 0) return 0;
    switch (x) {
        case 5: case 7: case 17: case 25: case 37: case 41:
        case 55: case 61: case 91: return 1;
    }
    return 0;
}
static int cycle_id(u64 x) {      /* 0: {1}; 1: {5,7}; 2: seven-cycle */
    if (x == 1) return 0;
    if (x == 5 || x == 7) return 1;
    return 2;
}
static inline u64 ustep(u64 x, int *k) {
    if (x > (0xFFFFFFFFFFFFFFFFULL - 2ULL) / 3ULL) {
        fprintf(stderr, "overflow at x=%llu\n", x); exit(3);
    }
    u64 y = 3ULL * x + (B > 0 ? 1ULL : 0ULL) - (B > 0 ? 0ULL : 1ULL);
    int kk = __builtin_ctzll(y);
    *k = kk;
    return y >> kk;
}
static inline int sres(u64 x, int M) {    /* signed residue b*x mod M */
    int r = (int)(x % (u64)M);
    if (B > 0) return r;
    return (M - r) % M;
}

/* ---------------- patterns ---------------- */
static const int FACT[8] = {1, 1, 2, 6, 24, 120, 720, 5040};
#define NSCOPE 4
typedef struct { u64 cnt; u64 first_n; int first_t; u64 last_n; } pstat_t;
static pstat_t PS[7][NSCOPE][720];        /* [w][scope][idx] */
static u64 PDEV[7][NSCOPE];               /* windows with actual != generic */
static u64 PDEV_first[7][NSCOPE];
static u64 PWIN[7][NSCOPE];

static int lehmer_u(const u64 *v, int w) {
    int idx = 0;
    for (int i = 0; i < w; i++) {
        int c = 0;
        for (int j = i + 1; j < w; j++) if (v[j] < v[i]) c++;
        idx += c * FACT[w - 1 - i];
    }
    return idx;
}
static int lehmer_d(const double *v, int w) {
    int idx = 0;
    for (int i = 0; i < w; i++) {
        int c = 0;
        for (int j = i + 1; j < w; j++) if (v[j] < v[i]) c++;
        idx += c * FACT[w - 1 - i];
    }
    return idx;
}

/* ---------------- lag tables ---------------- */
/* [lag][from_start?1:0][generic growth?1:0][actual up?1:0] */
static u64 LT[4097][2][2][2];
static u64 LT_first[4097][2][2][2];

/* ---------------- hash of windows (x_i, l) ---------------- */
typedef struct { u64 x; int l; int used; u64 xj; int K; int flags; u64 first_n; double mu; } hent_t;
typedef struct { hent_t *t; size_t cap; size_t n; } hash_t;
static void h_init(hash_t *h, size_t cap) { h->cap = cap; h->n = 0; h->t = calloc(cap, sizeof(hent_t)); if (!h->t) { fprintf(stderr, "oom\n"); exit(4);} }
static inline size_t hkey(u64 x, int l, size_t cap) {
    u64 z = x * 0x9E3779B97F4A7C15ULL ^ ((u64)l * 0xC2B2AE3D27D4EB4FULL);
    z ^= z >> 29; z *= 0xBF58476D1CE4E5B9ULL; z ^= z >> 32;
    return (size_t)(z & (cap - 1));
}
static hent_t *h_get(hash_t *h, u64 x, int l, int *isnew);
static void h_grow(hash_t *h) {
    hash_t nh; h_init(&nh, h->cap * 2);
    for (size_t i = 0; i < h->cap; i++) if (h->t[i].used) {
        int isn; hent_t *e = h_get(&nh, h->t[i].x, h->t[i].l, &isn); *e = h->t[i];
    }
    free(h->t); *h = nh;
}
static hent_t *h_get(hash_t *h, u64 x, int l, int *isnew) {
    if (h->n * 2 > h->cap) h_grow(h);
    size_t i = hkey(x, l, h->cap);
    while (h->t[i].used) {
        if (h->t[i].x == x && h->t[i].l == l) { *isnew = 0; return &h->t[i]; }
        i = (i + 1) & (h->cap - 1);
    }
    h->t[i].used = 1; h->t[i].x = x; h->t[i].l = l; h->n++; *isnew = 1;
    return &h->t[i];
}
static hash_t HDEV, HMAR, HCST;

/* ---------------- records ---------------- */
/* feature ids */
enum { RF_VS16, RF_VS27, RF_PS16, RF_PS27, RF_VR16, RF_VR27, RF_PR16, RF_PR27, RF_K, RF_GAP, RF_IDX, RF_NF };
static const int RF_SIZE[RF_NF] = {16, 27, 16, 27, 16, 27, 16, 27, 64, 256, 1024};
static const char *RF_NAME[RF_NF] = {"val_s16", "val_s27", "pred_s16", "pred_s27", "val_r16", "val_r27", "pred_r16", "pred_r27", "kprev", "gap", "idx"};
/* [type 0=max 1=min][scope 0=all 1=interior][feature][value] */
static u64 *RC[2][2][RF_NF];
static u64 *RCF[2][2][RF_NF];    /* first witness n */

/* ---------------- per-n summary ---------------- */
typedef struct {
    int16_t sigma, tau, T, argmax, rho, rmax, rmin, Ksig;
    int16_t fmr, nD, nU, KT;
    int8_t entry, cst;
    int8_t pad[6];
    u64 M;
} pern_t;

static u64 xs[MAXT + EXTRA + 2];
static int ks[MAXT + EXTRA + 2];
static int Kp[MAXT + EXTRA + 3];

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: %s b N prefix [lagmax]\n", argv[0]); return 1; }
    B = atoi(argv[1]);
    u64 N = strtoull(argv[2], 0, 10);
    const char *pre = argv[3];
    if (argc > 4) LAGMAX = atoi(argv[4]);
    if (LAGMAX > 4096) LAGMAX = 4096;
    if (B != 1 && B != -1) { fprintf(stderr, "b must be +-1\n"); return 1; }

    memset(PS, 0, sizeof(PS)); memset(LT, 0, sizeof(LT)); memset(LT_first, 0, sizeof(LT_first));
    h_init(&HDEV, 1 << 16); h_init(&HMAR, 1 << 16); h_init(&HCST, 1 << 10);
    for (int t = 0; t < 2; t++) for (int s = 0; s < 2; s++) for (int f = 0; f < RF_NF; f++) {
        RC[t][s][f] = calloc(RF_SIZE[f], sizeof(u64)); RCF[t][s][f] = calloc(RF_SIZE[f], sizeof(u64));
    }
    char fn[1024];
    snprintf(fn, sizeof fn, "%s_pern.bin", pre);
    FILE *fb = fopen(fn, "wb");
    if (!fb) { perror(fn); return 2; }

    u64 norbits = 0, maxval = 0; int maxT = 0; u64 maxT_n = 0;
    u64 cst_start_fail = 0;

    for (u64 n = 1; n <= N; n += 2) {
        /* stopped orbit */
        int T = 0; u64 x = n; xs[0] = n;
        while (!is_cycle(x)) {
            int k; u64 y = ustep(x, &k); ks[T] = k; T++;
            if (T >= MAXT) { fprintf(stderr, "orbit too long n=%llu\n", n); return 3; }
            xs[T] = y; x = y;
        }
        /* continuation into the cycle (for sigma, tau, k_T) */
        int E = 0;
        {
            u64 z = xs[T];
            for (E = 0; E < EXTRA; E++) { int k; u64 y = ustep(z, &k); ks[T + E] = k; xs[T + E + 1] = y; z = y; }
        }
        Kp[0] = 0;
        for (int j = 0; j < T + EXTRA; j++) Kp[j + 1] = Kp[j] + ks[j];
        norbits++;
        if (T > maxT) { maxT = T; maxT_n = n; }

        pern_t P; memset(&P, 0, sizeof P);
        P.T = (int16_t)T; P.entry = (int8_t)cycle_id(xs[T]); P.KT = (int16_t)Kp[T];
        /* sigma, tau over stopped orbit + continuation */
        int sig = -1, tau = -1;
        for (int j = 1; j <= T + EXTRA; j++) { if (xs[j] < n) { sig = j; break; } }
        for (int j = 1; j <= T + EXTRA; j++) { if ((double)Kp[j] > j * LOG23) { tau = j; break; } }
        P.sigma = (int16_t)sig; P.tau = (int16_t)tau; P.Ksig = (int16_t)(sig > 0 ? Kp[sig] : -1);
        P.cst = (int8_t)(sig == tau);
        if (!is_cycle(n) && sig != tau) cst_start_fail++;
        /* max, argmax, rho, records */
        u64 M = xs[0]; int am = 0; int rho = 0;
        for (int j = 1; j <= T; j++) { if (xs[j] > M) { M = xs[j]; am = j; } if (xs[j] > n) rho++; if (xs[j] > maxval) maxval = xs[j]; }
        if (xs[0] > maxval) maxval = xs[0];
        P.M = M; P.argmax = (int16_t)am; P.rho = (int16_t)rho;
        {
            u64 cmax = xs[0], cmin = xs[0]; int lastmax = 0, lastmin = 0; int rmax = 0, rmin = 0; int fmr = -1;
            for (int j = 1; j <= T; j++) {
                for (int typ = 0; typ < 2; typ++) {
                    int isrec = (typ == 0) ? (xs[j] > cmax) : (xs[j] < cmin);
                    if (!isrec) continue;
                    int last = (typ == 0) ? lastmax : lastmin;
                    int fv[RF_NF];
                    fv[RF_VS16] = sres(xs[j], 16); fv[RF_VS27] = sres(xs[j], 27);
                    fv[RF_PS16] = sres(xs[j - 1], 16); fv[RF_PS27] = sres(xs[j - 1], 27);
                    fv[RF_VR16] = (int)(xs[j] % 16); fv[RF_VR27] = (int)(xs[j] % 27);
                    fv[RF_PR16] = (int)(xs[j - 1] % 16); fv[RF_PR27] = (int)(xs[j - 1] % 27);
                    fv[RF_K] = ks[j - 1] < 63 ? ks[j - 1] : 63;
                    fv[RF_GAP] = (j - last) < 255 ? (j - last) : 255;
                    fv[RF_IDX] = j < 1023 ? j : 1023;
                    for (int s = 0; s < 2; s++) {
                        if (s == 1 && j == T) continue;
                        for (int f = 0; f < RF_NF; f++) {
                            if (RC[typ][s][f][fv[f]] == 0) RCF[typ][s][f][fv[f]] = n;
                            RC[typ][s][f][fv[f]]++;
                        }
                    }
                    if (typ == 0) { cmax = xs[j]; lastmax = j; rmax++; if (fmr < 0) fmr = j; }
                    else { cmin = xs[j]; lastmin = j; rmin++; }
                }
            }
            P.rmax = (int16_t)rmax; P.rmin = (int16_t)rmin; P.fmr = (int16_t)fmr;
        }
        /* patterns */
        for (int w = 3; w <= 6; w++) {
            for (int t = 0; t + w - 1 <= T; t++) {
                int ia = lehmer_u(&xs[t], w);
                double L[6];
                u64 mn = xs[t];
                for (int i = 0; i < w; i++) { L[i] = i * LOG23 - (double)(Kp[t + i] - Kp[t]); if (xs[t + i] < mn) mn = xs[t + i]; }
                int ig = lehmer_d(L, w);
                int sc[NSCOPE];
                sc[0] = 1; sc[1] = (t == 0); sc[2] = (t + w - 1 < T); sc[3] = (mn >= 1000);
                for (int s = 0; s < NSCOPE; s++) {
                    if (!sc[s]) continue;
                    pstat_t *p = &PS[w][s][ia];
                    if (p->cnt == 0) { p->first_n = n; p->first_t = t; }
                    p->cnt++; p->last_n = n;
                    PWIN[w][s]++;
                    if (ia != ig) { if (PDEV[w][s] == 0) PDEV_first[w][s] = n; PDEV[w][s]++; }
                }
            }
        }
        /* all pairs: lag tables, deviations, margins, CST along the orbit */
        int nD = 0, nU = 0;
        for (int i = 0; i < T; i++) {
            u64 xi = xs[i]; int K0 = Kp[i];
            int sig_i = -1, tau_i = -1;
            for (int j = i + 1; j <= T; j++) {
                int l = j - i; int kk = Kp[j] - K0;
                int growth = ((double)kk < l * LOG23);
                int up = xs[j] > xi;
                if (sig_i < 0 && xs[j] < xi) sig_i = j;
                if (tau_i < 0 && !growth) tau_i = j;
                if (l <= LAGMAX) {
                    int fs = (i == 0);
                    if (LT[l][0][growth][up] == 0) LT_first[l][0][growth][up] = n;
                    LT[l][0][growth][up]++;
                    if (fs) { if (LT[l][1][growth][up] == 0) LT_first[l][1][growth][up] = n; LT[l][1][growth][up]++; }
                }
                if (growth != up) {
                    if (growth) nD++; else nU++;
                    int isn; hent_t *e = h_get(&HDEV, xi, l, &isn);
                    if (isn) { e->xj = xs[j]; e->K = kk; e->flags = (j == T ? 1 : 0); e->first_n = n; e->mu = 0; }
                    if (i == 0) e->flags |= 2;
                }
                /* positive-gate windows: plus decay, minus growth */
                int posgate = (B > 0) ? !growth : growth;
                if (posgate) {
                    double r = exp2(l * LOG23 - (double)kk);
                    double num = (double)xi - (double)xs[j];
                    double den = (double)xs[j] - r * (double)xi;
                    double mu = num / den;
                    if (mu < 1.0) {
                        int isn; hent_t *e = h_get(&HMAR, xi, l, &isn);
                        if (isn) { e->xj = xs[j]; e->K = kk; e->flags = (j == T ? 1 : 0); e->first_n = n; e->mu = mu; }
                    }
                }
            }
            /* CST at orbit point x_i, decided only if both exist inside the stopped orbit */
            if (sig_i >= 0 && tau_i >= 0 && sig_i != tau_i) {
                int isn; hent_t *e = h_get(&HCST, xi, sig_i - i, &isn);
                if (isn) { e->xj = xs[sig_i]; e->K = tau_i - i; e->first_n = n; e->flags = (i == 0 ? 2 : 0); }
            }
        }
        P.nD = (int16_t)(nD < 32767 ? nD : 32767); P.nU = (int16_t)(nU < 32767 ? nU : 32767);
        fwrite(&P, sizeof P, 1, fb);
    }
    fclose(fb);

    /* ------------- write outputs ------------- */
    FILE *f;
    snprintf(fn, sizeof fn, "%s_patterns.txt", pre); f = fopen(fn, "w");
    fprintf(f, "# b=%d N=%llu orbits=%llu maxT=%d (n=%llu) maxval=%llu cst_start_fail=%llu\n", B, N, norbits, maxT, maxT_n, maxval, cst_start_fail);
    for (int w = 3; w <= 6; w++) for (int s = 0; s < NSCOPE; s++) {
        fprintf(f, "W %d %d windows=%llu dev=%llu devfirst=%llu\n", w, s, PWIN[w][s], PDEV[w][s], PDEV_first[w][s]);
        for (int i = 0; i < FACT[w]; i++) {
            pstat_t *p = &PS[w][s][i];
            fprintf(f, "P %d %d %d %llu %llu %d %llu\n", w, s, i, p->cnt, p->first_n, p->first_t, p->last_n);
        }
    }
    fclose(f);

    snprintf(fn, sizeof fn, "%s_lags.txt", pre); f = fopen(fn, "w");
    fprintf(f, "# lag fs growth up count first_n\n");
    for (int l = 1; l <= LAGMAX; l++) for (int fs = 0; fs < 2; fs++) for (int g = 0; g < 2; g++) for (int u = 0; u < 2; u++)
        fprintf(f, "L %d %d %d %d %llu %llu\n", l, fs, g, u, LT[l][fs][g][u], LT_first[l][fs][g][u]);
    fclose(f);

    snprintf(fn, sizeof fn, "%s_devs.txt", pre); f = fopen(fn, "w");
    fprintf(f, "# x_i lag x_j K flags(1=terminal,2=from-start-seen) first_n\n");
    for (size_t i = 0; i < HDEV.cap; i++) if (HDEV.t[i].used) {
        hent_t *e = &HDEV.t[i];
        fprintf(f, "D %llu %d %llu %d %d %llu\n", e->x, e->l, e->xj, e->K, e->flags, e->first_n);
    }
    fclose(f);

    snprintf(fn, sizeof fn, "%s_margins.txt", pre); f = fopen(fn, "w");
    fprintf(f, "# x_i lag x_j K mu flags first_n\n");
    for (size_t i = 0; i < HMAR.cap; i++) if (HMAR.t[i].used) {
        hent_t *e = &HMAR.t[i];
        fprintf(f, "M %llu %d %llu %d %.12g %d %llu\n", e->x, e->l, e->xj, e->K, e->mu, e->flags, e->first_n);
    }
    fclose(f);

    snprintf(fn, sizeof fn, "%s_cst.txt", pre); f = fopen(fn, "w");
    fprintf(f, "# x sigma_lag x_sigma tau_lag flags first_n\n");
    for (size_t i = 0; i < HCST.cap; i++) if (HCST.t[i].used) {
        hent_t *e = &HCST.t[i];
        fprintf(f, "C %llu %d %llu %d %d %llu\n", e->x, e->l, e->xj, e->K, e->flags, e->first_n);
    }
    fclose(f);

    snprintf(fn, sizeof fn, "%s_records.txt", pre); f = fopen(fn, "w");
    fprintf(f, "# type(0max,1min) scope(0all,1interior) feature value count first_n\n");
    for (int t = 0; t < 2; t++) for (int s = 0; s < 2; s++) for (int ff = 0; ff < RF_NF; ff++)
        for (int v = 0; v < RF_SIZE[ff]; v++) if (RC[t][s][ff][v])
            fprintf(f, "R %d %d %s %d %llu %llu\n", t, s, RF_NAME[ff], v, RC[t][s][ff][v], RCF[t][s][ff][v]);
    fclose(f);
    fprintf(stderr, "done b=%d N=%llu orbits=%llu maxT=%d devs=%zu margins=%zu cst=%zu\n", B, N, norbits, maxT, HDEV.n, HMAR.n, HCST.n);
    return 0;
}
