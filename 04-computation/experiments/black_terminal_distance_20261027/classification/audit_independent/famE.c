/* famE.c -- auditor's own construction of the family E from the DEFINITION in PROOF.md Section 3, and checks of
 * Part VII / VIII claims for every member with n <= nmax.
 *
 * Bottom layer (S, x0): S - x0 = R_1 => ... => R_q (q >= 2), blocks K1 or C3, R_q = C3, x0 uniform on every block,
 * x0 => R_1, R_q => x0, relation of x0 to R_2..R_(q-1) free.  Enumerated as (composition of |S|-4 into parts 1,3 for
 * R_1..R_(q-1)) x (2^(q-2) bits).
 * E_t(S, x0): L_0 = S, middles c_1..c_(t-1), top L_t = sigma(S) with the converse; only upward arcs
 * x0 -> c_1 -> ... -> c_(t-1) -> sigma(x0) (t = 1: x0 -> sigma(x0)); every other inter-layer pair points down.
 *
 * Checks: valid, strong; sigma (with c_i <-> c_(t-i)) involutive anti-automorphism; lambda = 1; delta = 2;
 * minimum cuts exactly X_k = L_0 + {c_1..c_(k-1)}, k = 1..t; candidate count = h^2 with h = 3^(#C3 blocks);
 * all candidates blue (no black descent); >= 3 vertices of outdegree 2; Part VIII.4 certificate: the word
 * Q = (r, r', r'', x0, c_1..c_(t-1), P_top, P_rest) is a black Ham path, K = K_Q T has lambda = delta = 1, and K has a
 * black descent Q' with K_Q'(K) non-strong (so d_black <= 2; d_black >= 2 from 'no black descent' + F1).
 * Canonical forms (nauty) give the number of distinct members per order (claimed e(n) = 1,1,3,3,8,8,20,20,...).
 * Usage: ./famE nmax [printcodes]
 */
#include "at_nauty.h"
#include "cand.h"

typedef struct { int a; int q; int sz[16]; int bit[16]; } blayer;   /* bit[i] = 1 if x0 -> R_i */

static int nBL = 0; static blayer BL[4096];
static void gen_comp(int rem, int *parts, int k, int a) {
    if (rem == 0) {
        if (k == 0) return;                                  /* need q - 1 >= 1 blocks before R_q */
        int q = k + 1; int nb = q >= 3 ? q - 2 : 0;
        for (int bits = 0; bits < (1 << nb); bits++) {
            blayer B; memset(&B, 0, sizeof B); B.a = a; B.q = q;
            for (int i = 0; i < k; i++) B.sz[i] = parts[i];
            B.sz[k] = 3;
            B.bit[0] = 1;                                     /* x0 => R_1 */
            for (int i = 1; i <= q - 2; i++) B.bit[i] = (bits >> (i - 1)) & 1;
            B.bit[q - 1] = 0;                                 /* R_q => x0 */
            BL[nBL++] = B;
        }
        return;
    }
    if (rem >= 1) { parts[k] = 1; gen_comp(rem - 1, parts, k + 1, a); }
    if (rem >= 3) { parts[k] = 3; gen_comp(rem - 3, parts, k + 1, a); }
}

/* S on a vertices: x0 = 0, blocks consecutive; returns block start indices */
static void build_S(const blayer *B, tour *S, int *bstart) {
    int a = B->a; t_clear(S, a);
    int idx = 1;
    for (int i = 0; i < B->q; i++) { bstart[i] = idx; idx += B->sz[i]; }
    bstart[B->q] = idx;
    if (idx != a) { fprintf(stderr, "size mismatch\n"); exit(1); }
    for (int i = 0; i < B->q; i++) {
        int b0 = bstart[i];
        if (B->sz[i] == 3) { t_arc(S, b0, b0 + 1); t_arc(S, b0 + 1, b0 + 2); t_arc(S, b0 + 2, b0); }
        for (int j = i + 1; j < B->q; j++) for (int u = bstart[i]; u < bstart[i + 1]; u++) for (int v = bstart[j]; v < bstart[j + 1]; v++) t_arc(S, u, v);
        for (int u = bstart[i]; u < bstart[i + 1]; u++) { if (B->bit[i]) t_arc(S, 0, u); else t_arc(S, u, 0); }
    }
    t_fin(S);
}

/* E_t(S, x0=0): bottom 0..a-1, middles c_k = a+k-1 (k=1..t-1), top sigma(u) = a+t-1+u */
static void build_E(const tour *S, int t, tour *T, int *sigma) {
    int a = S->n, n = 2 * a + t - 1; t_clear(T, n);
    int top0 = a + t - 1;
    /* default: higher layer -> lower layer; order of layers: bottom (0), middles by k, top */
    int layer[MAXV];
    for (int u = 0; u < a; u++) layer[u] = 0;
    for (int k = 1; k <= t - 1; k++) layer[a + k - 1] = k;
    for (int u = 0; u < a; u++) layer[top0 + u] = t;
    for (int u = 0; u < n; u++) for (int v = u + 1; v < n; v++) {
        if (layer[u] != layer[v]) { if (layer[u] > layer[v]) t_arc(T, u, v); else t_arc(T, v, u); }
    }
    for (int u = 0; u < a; u++) for (int v = 0; v < a; v++) if (u != v && t_has(S, u, v)) { t_arc(T, u, v); t_arc(T, top0 + v, top0 + u); }
    /* jumps */
    if (t == 1) t_arc(T, 0, top0 + 0);
    else { t_arc(T, 0, a); for (int k = 1; k + 1 <= t - 1; k++) t_arc(T, a + k - 1, a + k); t_arc(T, a + t - 2, top0 + 0); }
    t_fin(T);
    for (int u = 0; u < a; u++) { sigma[u] = top0 + u; sigma[top0 + u] = u; }
    for (int k = 1; k <= t - 1; k++) sigma[a + k - 1] = a + (t - k) - 1;
}

typedef struct { const tour *T; long long cnt, blue; } cnt_ctx;
static int cnt_cb(const int *Q, void *vc) { cnt_ctx *c = (cnt_ctx *)vc; c->cnt++; if (t_blue(c->T, Q)) c->blue++; return 0; }

int main(int argc, char **argv) {
    int nmax = argc > 1 ? atoi(argv[1]) : 16; int printcodes = argc > 2;
    nauty_check(WORDSIZE, 1, nmax, NAUTYVERSIONID);
    int parts[16];
    int amax = nmax / 2;
    long long members = 0, allok = 0;
    /* canonical forms per n for distinctness */
    static char canon[64][2000][200]; static int ncanon[64];
    for (int a = 5; a <= amax; a++) {
        int before = nBL; gen_comp(a - 4, parts, 0, a);
        printf("a=%d: bottom layers from the definition b(a)=%d\n", a, nBL - before);
    }
    for (int bi = 0; bi < nBL; bi++) {
        blayer *B = &BL[bi]; tour S; int bstart[17]; build_S(B, &S, bstart);
        int nC3 = 0; for (int i = 0; i < B->q; i++) if (B->sz[i] == 3) nC3++;
        long long h = 1; for (int i = 0; i < nC3; i++) h *= 3;
        /* check S is strong and the definition's features */
        if (!t_strong_mask(&S, fullmask(S.n))) { printf("S not strong?!\n"); return 1; }
        for (int t = 1; 2 * B->a + t - 1 <= nmax; t++) {
            int n = 2 * B->a + t - 1; tour T; int sg[MAXV]; build_E(&S, t, &T, sg); members++;
            int ok = 1; char why[400] = "";
            if (!t_valid(&T)) { ok = 0; strcat(why, " invalid"); }
            if (!t_strong_mask(&T, fullmask(n))) { ok = 0; strcat(why, " notstrong"); }
            int sgsg = 1; for (int v = 0; v < n; v++) if (sg[sg[v]] != v) sgsg = 0;
            if (!t_is_anti(&T, sg) || !sgsg) { ok = 0; strcat(why, " sigma"); }
            static uint32_t cuts[1 << 16]; int nc; int lam = t_mincuts(&T, cuts, 1 << 16, &nc);
            int del = t_delta(&T);
            if (lam != 1) { ok = 0; strcat(why, " lambda"); }
            if (del != 2) { ok = 0; strcat(why, " delta"); }
            /* min cuts exactly the chain */
            int chainok = (nc == t);
            for (int k = 1; k <= t && chainok; k++) {
                uint32_t X = fullmask(B->a); for (int j = 1; j <= k - 1; j++) X |= 1u << (B->a + j - 1);
                int fnd = 0; for (int i = 0; i < nc; i++) if (cuts[i] == X) fnd = 1;
                if (!fnd) chainok = 0;
            }
            if (!chainok) { ok = 0; strcat(why, " cuts"); }
            uint32_t sep[MAXV]; make_sep(&T, cuts, nc, sep);
            cnt_ctx cc = { &T, 0, 0 }; cand_enum(&T, sep, cnt_cb, &cc);
            if (cc.cnt != h * h) { ok = 0; strcat(why, " candcount"); }
            if (cc.blue != cc.cnt) { ok = 0; strcat(why, " BLACKCAND"); }
            int deg2 = 0; for (int v = 0; v < n; v++) if (T.s[v] == 2) deg2++;
            if (deg2 < 3) { ok = 0; strcat(why, " deg2"); }
            /* Part VIII.4 certificate */
            int rq = bstart[B->q - 1]; int r = rq, r1 = rq + 1, r2 = rq + 2;   /* r -> r1 -> r2 -> r in build_S */
            int Q[MAXV], k = 0; Q[k++] = r; Q[k++] = r1; Q[k++] = r2; Q[k++] = 0;
            for (int j = 1; j <= t - 1; j++) Q[k++] = B->a + j - 1;
            int top0 = B->a + t - 1; uint32_t topmask = 0; for (int u = 0; u < B->a; u++) topmask |= 1u << (top0 + u);
            int P[MAXV];
            int lt = t_hp_in(&T, topmask, top0, topmask, P); for (int i = 0; i < lt; i++) Q[k++] = P[i];
            uint32_t restmask = 0; for (int u = 1; u < rq; u++) restmask |= 1u << u;
            int lr = t_hp_in(&T, restmask, -1, restmask, P); for (int i = 0; i < lr; i++) Q[k++] = P[i];
            int certok = (k == n) && t_is_hpath(&T, Q) && !t_blue(&T, Q);
            tour K; int lamK = -1, delK = -1, lamK2 = -1;
            if (certok) {
                t_move(&T, Q, &K);
                int nck; static uint32_t ck[1 << 16]; lamK = t_mincuts(&K, ck, 1 << 16, &nck); delK = t_delta(&K);
                int inr = n - 1 - K.s[r];
                certok = (lamK == 1 && delK == 1 && inr == 1);
                int Q2[MAXV]; int l2;
                if (certok && black_descent(&K, 0, Q2, &l2, NULL, NULL)) {
                    tour K2; t_move(&K, Q2, &K2); int nc2; static uint32_t c2[1 << 16]; lamK2 = t_mincuts(&K2, c2, 1 << 16, &nc2);
                    certok = (lamK2 == 0) && !t_blue(&K, Q2) && !t_strong_mask(&K2, fullmask(n));
                } else certok = 0;
            }
            if (!certok) { ok = 0; strcat(why, " cert"); }
            double aut = nt_autsize(&T);
            char cs[200]; nt_canon_str(&T, cs);
            int dup = 0; for (int i = 0; i < ncanon[n]; i++) if (!strcmp(canon[n][i], cs)) dup = 1;
            if (!dup) strcpy(canon[n][ncanon[n]++], cs);
            allok += ok;
            char code[512]; t_code(&T, code);
            printf("n=%2d a=%d t=%d q=%d blocks=", n, B->a, t, B->q);
            for (int i = 0; i < B->q; i++) printf("%s%c", B->sz[i] == 3 ? "C3" : "K1", B->bit[i] ? '+' : '-');
            printf(" ('+': x0 beats the block) h=%lld cand=%lld blue=%lld |Aut|=%g deg2=%d lamK=%d delK=%d lamK2=%d %s%s\n",
                   h, cc.cnt, cc.blue, aut, deg2, lamK, delK, lamK2, ok ? "OK" : "FAIL:", why);
            if (printcodes) printf("   code %s\n", code);
        }
    }
    printf("members built=%lld all-checks-ok=%lld\n", members, allok);
    printf("distinct members (nauty canonical form) by order:");
    for (int n = 10; n <= nmax; n++) printf(" n=%d:%d", n, ncanon[n]);
    printf("\n");
    return 0;
}
