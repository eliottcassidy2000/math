/* lrc_lift_pyramid_S53.c -- mac-mini-2026-07-05-S53, HYP-4109.
 *
 * The l=3,4(,5) FLOOR-LEVEL lift sweeps at beta = 2/25, C speed.
 * Same mathematics as lrc_l3_floor_macmini_S53.py (see its header for the
 * chain-domain derivation); this adds levels for l=4/5 and runs ~100x faster.
 *
 * Chain domain (sorted lifted values w1 <= ... <= wl), all citation-derived:
 *   l=3: 13*w1 <= 3600;              51*w2 <= 1100*w1; w3 <= 24*w2
 *   l=4:  7*w1 <= 2400; 13*w2 <= 300*w1; 51*w3 <= 1100*w2; w4 <= 24*w3
 *   l=5:  3*w1 <= 1600;  7*w2 <= 200*w1; 13*w3 <= 300*w2; 51*w4 <= 1100*w3;
 *         w5 <= 24*w4
 * (anchors A_l = l*beta*12/((1/(13-l)-beta)(1-2*l*beta)), ratios
 *  R_j = j*beta/((1/(13-j)-beta)(1-2*j*beta)), R_1 = 24, R_2 = 1100/51,
 *  R_3 = 300/13, R_4 = 200/7.)
 *
 * Pyramid: at a k-prefix node, try q = 13u (u = 2..5) witnesses with all
 * fixed values clear at full q precision and every free coordinate's residue
 * orbit inside the mod-13 window [u+1, 12-u]; success kills the whole
 * subtree.  Last axis cleared by residue classes with general-q witnesses.
 * Leftover cells: sieve, primitivity, adaptive scan; survivors printed for
 * the python exact-M post-pass (expected: none or a handful).
 *
 * Output: per-l stats + "CELL ..." lines for unresolved sets.
 * gcc -O2 -o lrc_lift_pyramid_S53 lrc_lift_pyramid_S53.c && ./lrc_lift_pyramid_S53 3 4
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef long long ll;

static int gcd_i(int a, int b){ while(b){ int t=a%b; a=b; b=t; } return a; }

/* margin >= 2/25 at residue x mod q  <=>  25*dist >= 2q */
static inline int okq(int x, int q){ int r = x % q; if (r < 0) r += q; int d = r < q - r ? r : q - r; return 25*d >= 2*q; }
static inline int distq(int x, int q){ int r = x % q; if (r < 0) r += q; return r < q - r ? r : q - r; }
static inline int dist13(int x){ int r = x % 13; if (r < 0) r += 13; return r < 13 - r ? r : 13 - r; }

/* ---------------- global config ---------------- */
static int QGEN[128]; static int nQGEN = 0;   /* general witness moduli */
static void init_qgen(void){
    for (int q = 8; q <= 41; q++) if (q % 13) QGEN[nQGEN++] = q;
    for (int u = 2; u <= 21; u++) QGEN[nQGEN++] = 13*u;
    QGEN[nQGEN++] = 25; QGEN[nQGEN++] = 50;
}

/* chain caps: value_cap[l][pos]: w_{pos} bounds (pos = 0 smallest) as
 * (num, den) pairs vs previous value (pos>0) or absolute (pos=0). */
static const ll ANCH_N[6] = {0,0,0, 3600, 2400, 1600};
static const ll ANCH_D[6] = {1,1,1,   13,    7,    3};
static const ll RAT_N[6]  = {0, 24, 1100, 300, 200, 0};  /* R_j numerators, j=1..4 */
static const ll RAT_D[6]  = {1,  1,   51,  13,   7, 1};

/* stats */
static ll st_nodes[8], st_killed[8], st_cells, st_class, st_scan, st_exact_todo, st_sieved, st_rows, st_rowkill;

/* per-C state */
static int base_[12], nbase, Cset[8], L;
static int cwq[4096][2]; static int ncw;    /* general (q,a) with base clear */
static int pw[6][64]; static int npw[6];    /* 13u-family a's with base clear */

static int sieve_ok_set(int *W, int n){
    for (int m = 2; m <= 12; m++){
        int hit = 0;
        for (int i = 0; i < n; i++) if (W[i] % m == 0){ hit = 1; break; }
        if (!hit) return 0;
    }
    return 1;
}
static int gcd_set(int *W, int n){ int g = 0; for (int i = 0; i < n; i++) g = gcd_i(g, W[i]); return g; }

/* adaptive per-cell scan; W sorted ascending, n = 12 */
static int cell_scan(int *W, int n){
    int qs[96]; int nq = 0;
    int wl = W[n-1], wm = W[n-2], ws = W[n-3];
    qs[nq++] = wl+1; qs[nq++] = wl-1; qs[nq++] = wm+wl; qs[nq++] = wl-wm;
    qs[nq++] = ws+wm; qs[nq++] = wm+1; qs[nq++] = ws+wl; qs[nq++] = wl-ws;
    qs[nq++] = wm-ws; qs[nq++] = ws+1;
    for (int u = 2; u <= 39; u++) qs[nq++] = 13*u;
    qs[nq++] = 25; qs[nq++] = 50;
    for (int qi = 0; qi < nq; qi++){
        int q = qs[qi];
        if (q < 8) continue;
        int bad = 0;
        for (int i = 0; i < n; i++) if (W[i] % q == 0){ bad = 1; break; }
        if (bad) continue;
        for (int a = 1; a <= q/2; a++){
            if (gcd_i(a, q) != 1) continue;
            int all = 1;
            for (int i = 0; i < n; i++){
                ll x = (ll)a * W[i] % q;
                int d = (int)(x < q - x ? x : q - x);
                if (25*d < 2*q){ all = 0; break; }
            }
            if (all) return 1;
        }
    }
    return 0;
}

/* recursive enumeration over sorted-value positions.
 * perm[pos] = coordinate residue at value-rank pos (0 = smallest).
 * kv[pos] = height; wv[pos] = value.  At each node try pyramid kill. */
static int perm_[6], kv_[6], wv_[6];
static FILE *cellout;

static void do_last_axis(void){
    /* positions 0..L-2 fixed; last coord = perm_[L-1], w in (wv_[L-2], 24*wv_[L-2]] */
    int c = perm_[L-1];
    int wprev = wv_[L-2];
    int kclo = (wprev - c) / 13 + 1; if (kclo < 1) kclo = 1;
    int kchi = (24*wprev - c) / 13;
    if (kchi < kclo) return;
    st_rows++;
    /* row kill: q = 13u, all fixed values clear full-q, last-coord orbit window */
    for (int u = 2; u <= 5; u++){
        int q = 13*u;
        for (int ai = 0; ai < npw[u]; ai++){
            int a = pw[u][ai];
            int all = 1;
            for (int p = 0; p < L-1; p++)
                if (wv_[p] % q == 0 || !okq((int)((ll)a*wv_[p] % q), q)){ all = 0; break; }
            if (!all) continue;
            if (dist13(a*c) >= u+1){ st_rowkill++; return; }
        }
    }
    int ncell = kchi - kclo + 1;
    st_cells += ncell;
    /* class clearing on k axis: shrinking uncleared list (compaction) */
    int *unc = malloc(sizeof(int) * ncell);
    int nunc = ncell;
    for (int i = 0; i < ncell; i++) unc[i] = kclo + i;
    unsigned char goodbuf[1024];
    for (int wi = 0; wi < ncw && nunc; wi++){
        int q = cwq[wi][0], a = cwq[wi][1];
        int all = 1;
        for (int p = 0; p < L-1; p++)
            if (wv_[p] % q == 0 || !okq((int)((ll)a*wv_[p] % q), q)){ all = 0; break; }
        if (!all) continue;
        int step = (int)((ll)13*a % q);
        int b0 = (int)((ll)c*a % q);
        if (step == 0){ if (okq(b0, q)) nunc = 0; continue; }
        int g = gcd_i(step, q), period = q / g;
        if (period > 1024) continue;
        int any = 0;
        for (int t = 0; t < period; t++){
            goodbuf[t] = okq(b0 + step*t, q) ? 1 : 0;
            any |= goodbuf[t];
        }
        if (!any) continue;
        int w2 = 0;
        for (int i = 0; i < nunc; i++)
            if (!goodbuf[unc[i] % period]) unc[w2++] = unc[i];
        nunc = w2;
    }
    st_class += ncell - nunc;
    if (nunc){
        for (int i = 0; i < nunc; i++){
            int kc = unc[i];
            int wc = c + 13*kc;
            int W[12]; int n = 0;
            for (int b = 0; b < nbase; b++) W[n++] = base_[b];
            for (int p = 0; p < L-1; p++) W[n++] = wv_[p];
            W[n++] = wc;
            /* insertion sort */
            for (int x = 1; x < n; x++){ int v = W[x], y = x-1; while (y >= 0 && W[y] > v){ W[y+1] = W[y]; y--; } W[y+1] = v; }
            if (!sieve_ok_set(W, n)){ st_sieved++; continue; }
            if (gcd_set(W, n) != 1) continue;
            if (cell_scan(W, n)){ st_scan++; continue; }
            st_exact_todo++;
            fprintf(cellout, "CELL");
            for (int x = 0; x < n; x++) fprintf(cellout, " %d", W[x]);
            fprintf(cellout, "\n");
            fflush(cellout);
        }
    }
    free(unc);
}

static void recurse(int pos){
    /* pyramid kill at this node: fixed values wv_[0..pos-1]; free coords perm_[pos..L-1] */
    st_nodes[pos]++;
    for (int u = 2; u <= 5; u++){
        int q = 13*u;
        for (int ai = 0; ai < npw[u]; ai++){
            int a = pw[u][ai];
            int all = 1;
            for (int p = 0; p < pos; p++)
                if (wv_[p] % q == 0 || !okq((int)((ll)a*wv_[p] % q), q)){ all = 0; break; }
            if (!all) continue;
            for (int p = pos; p < L; p++)
                if (dist13(a*perm_[p]) < u+1){ all = 0; break; }
            if (all){ st_killed[pos]++; return; }
        }
    }
    if (pos == L-1){ do_last_axis(); return; }
    int c = perm_[pos];
    if (pos == 0){
        ll an = ANCH_N[L], ad = ANCH_D[L];
        for (int k = 1;; k++){
            ll w = c + 13LL*k;
            if (ad*w > an) break;
            kv_[0] = k; wv_[0] = (int)w;
            recurse(1);
        }
    } else {
        ll rn = RAT_N[L-pos], rd = RAT_D[L-pos];   /* w_{pos} <= R_{L-pos} * w_{pos-1} */
        ll wp = wv_[pos-1];
        for (int k = 1;; k++){
            ll w = c + 13LL*k;
            if (w < wp) continue;           /* enforce sorted order */
            if (rd*w > rn*wp) break;
            kv_[pos] = k; wv_[pos] = (int)w;
            recurse(pos+1);
        }
    }
}

static int nextperm(int *a, int n){
    int i = n-2; while (i >= 0 && a[i] >= a[i+1]) i--;
    if (i < 0) return 0;
    int j = n-1; while (a[j] <= a[i]) j--;
    int t = a[i]; a[i] = a[j]; a[j] = t;
    for (int x = i+1, y = n-1; x < y; x++, y--){ t = a[x]; a[x] = a[y]; a[y] = t; }
    return 1;
}

static void run_level_range(int l, int lo, int hi){
    L = l;
    int cidx = -1;
    memset(st_nodes, 0, sizeof st_nodes); memset(st_killed, 0, sizeof st_killed);
    st_cells = st_class = st_scan = st_exact_todo = st_sieved = st_rows = st_rowkill = 0;
    int idx[8];
    for (int i = 0; i < l; i++) idx[i] = i+1;
    while (1){
        cidx++;
        if (cidx >= hi) break;
        if (cidx < lo) goto nextcombo;
        /* pattern C = idx[0..l-1] */
        for (int i = 0; i < l; i++) Cset[i] = idx[i];
        nbase = 0;
        for (int v = 1; v <= 12; v++){
            int inC = 0; for (int i = 0; i < l; i++) if (idx[i] == v) inC = 1;
            if (!inC) base_[nbase++] = v;
        }
        /* per-C witness tables */
        ncw = 0;
        for (int qi = 0; qi < nQGEN; qi++){
            int q = QGEN[qi];
            int bad = 0;
            for (int b = 0; b < nbase; b++) if (base_[b] % q == 0){ bad = 1; break; }
            if (bad) continue;
            for (int a = 1; a <= q/2; a++){
                if (gcd_i(a, q) != 1) continue;
                int all = 1;
                for (int b = 0; b < nbase; b++)
                    if (!okq((int)((ll)a*base_[b] % q), q)){ all = 0; break; }
                if (all && ncw < 4096){ cwq[ncw][0] = q; cwq[ncw][1] = a; ncw++; }
            }
        }
        for (int u = 2; u <= 5; u++){
            int q = 13*u; npw[u] = 0;
            for (int a = 1; a <= q/2; a++){
                if (gcd_i(a, q) != 1) continue;
                int all = 1;
                for (int b = 0; b < nbase; b++)
                    if (!okq((int)((ll)a*base_[b] % q), q)){ all = 0; break; }
                if (all && npw[u] < 64) pw[u][npw[u]++] = a;
            }
        }
        /* permutations of C as value-rank assignment */
        int pa[8]; for (int i = 0; i < l; i++) pa[i] = idx[i];
        /* sort pa for permutation start */
        for (int x = 1; x < l; x++){ int v = pa[x], y = x-1; while (y >= 0 && pa[y] > v){ pa[y+1] = pa[y]; y--; } pa[y+1] = v; }
        do {
            for (int i = 0; i < l; i++) perm_[i] = pa[i];
            recurse(0);
        } while (nextperm(pa, l));
        /* next combination */
        nextcombo: ;
        int i = l-1;
        while (i >= 0 && idx[i] == 12 - (l-1-i)) i--;
        if (i < 0) break;
        idx[i]++;
        for (int j = i+1; j < l; j++) idx[j] = idx[j-1] + 1;
    }
    printf("l=%d [%d,%d) done: rows=%lld rowkill=%lld cells=%lld class=%lld scan=%lld sieved=%lld EXACT_TODO=%lld\n",
           l, lo, hi, st_rows, st_rowkill, st_cells, st_class, st_scan, st_sieved, st_exact_todo);
    printf("  node/kill by level:");
    for (int p = 0; p < l; p++) printf("  L%d %lld/%lld", p, st_killed[p], st_nodes[p]);
    printf("\n"); fflush(stdout);
}

int main(int argc, char **argv){
    init_qgen();
    cellout = stdout;
    for (int i = 1; i < argc; i++){
        int l, lo = 0, hi = 1000000;
        if (sscanf(argv[i], "%d:%d:%d", &l, &lo, &hi) < 1) continue;
        if (l >= 3 && l <= 5) run_level_range(l, lo, hi);
    }
    return 0;
}
