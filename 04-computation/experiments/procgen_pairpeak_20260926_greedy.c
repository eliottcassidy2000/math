/* procgen_pairpeak_20260926_greedy -- consistent (single-member) barrier construction in the Collatz pairing family.
 *
 * Pairing family (THM-4470/4475): pairs {2i-1, 2i}, one bit per pair; flipped pair: 2i-1 -> i-1, 2i -> 3i.
 * P_L: every n >= 3 falls below itself within L steps.  Sequential construction (THM-4475 style): process n = 3..N;
 * the default path uses frozen bits and reads free bits as 0 (Collatz); if it descends within L its pairs are frozen,
 * otherwise n is rescued.  Rescue order:
 *   (1) barrier rescues at scales W = WMAX, WMAX-1, ..., 2: follow frozen bits; at a FREE odd point x that is an
 *       up-image (previous step an unflipped odd step, so x = 2 mod 3 and its partner x+1 = 0 mod 3 is isolated),
 *       with x = 3 mod 4 (partner x+1 -> 3(x+1)/2 -> 3(x+1)/4 descends in 2 steps) and x >= n 3^(W-1), flip it
 *       (at most JMAX flips), protecting pair(T(x)) at 0; accept if n descends within L; freeze the certificate.
 *   (2) fallback = THM-4475's options A/F/B (then P, then a depth-first search with <= 3 flips).
 * MODE 0: every certificate is frozen (a valid section of a member of P_L; the runner verifies all n <= N).
 * MODE 1: default paths of Collatz-good n are NOT frozen (interference diagnostic; may break earlier n).
 * MODE 2: no barrier rescues (THM-4475's G_L; simulator check against the audited densities).
 * Instrumentation: private ("ideal") best scale (pure Collatz + own barrier rule), first deviation type of the actual
 * barrier path from the private one when the actual scale is lower (degraded), harmed rescues (Collatz-good n that
 * need rescue), riders hitting a high flip first (early / late with respect to the flip owner's time).
 * usage: greedy L N WMAX MODE [JMAX]
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
typedef unsigned long long u64;
static int L, WMAX, MODE, JMAX = 1000;
static u64 N, DENSE;
static signed char *dense; static u64 *hkey; static signed char *hval; static u64 HCAP, HCNT;
static inline u64 hsh(u64 k) { k ^= k >> 33; k *= 0xff51afd7ed558ccdULL; k ^= k >> 33; k *= 0xc4ceb9fe1a85ec53ULL; k ^= k >> 33; return k; }
static int getb(u64 i) { if (i < DENSE) return dense[i]; u64 h = hsh(i) & (HCAP - 1); while (hkey[h]) { if (hkey[h] == i) return hval[h]; h = (h + 1) & (HCAP - 1); } return -1; }
static void setb(u64 i, int b) {
    if (i < DENSE) { dense[i] = (signed char)b; return; }
    u64 h = hsh(i) & (HCAP - 1);
    while (hkey[h]) { if (hkey[h] == i) { hval[h] = (signed char)b; return; } h = (h + 1) & (HCAP - 1); }
    if (HCNT * 10 > HCAP * 7) { fprintf(stderr, "hash full\n"); exit(2); }
    hkey[h] = i; hval[h] = (signed char)b; HCNT++;
}
static inline u64 pairof(u64 v) { return (v + 1) >> 1; }
static inline u64 Tc(u64 v) { return (v & 1) ? (3 * v + 1) / 2 : v / 2; }
static inline u64 stepb(u64 v, int b) { if (b == 1) return (v & 1) ? (v - 1) / 2 : 3 * (v / 2); return Tc(v); }
typedef struct { u64 i[200]; int b[200]; int k; } Ext;
static int extget(const Ext *e, u64 i) { if (!e) return -1; for (int j = 0; j < e->k; j++) if (e->i[j] == i) return e->b[j]; return -1; }
static void extset(Ext *e, u64 i, int b) { e->i[e->k] = i; e->b[e->k] = b; e->k++; }
static int bitat(const Ext *e, u64 i) { int b = getb(i); if (b >= 0) return b; b = extget(e, i); return b >= 0 ? b : 0; }
static int isfree(const Ext *e, u64 i) { return getb(i) < 0 && extget(e, i) < 0; }
static u64 usedI[256]; static int usedB[256]; static int usedK;
static int run(u64 n, const Ext *e) { u64 v = n; usedK = 0; for (int s = 1; s <= L; s++) { u64 i = pairof(v); int b = bitat(e, i); usedI[usedK] = i; usedB[usedK] = b; usedK++; v = stepb(v, b); if (v < n) return s; } return 0; }
static void freeze_used(void) { for (int j = 0; j < usedK; j++) if (getb(usedI[j]) < 0) setb(usedI[j], usedB[j]); }
static void commit(const Ext *e) { for (int j = 0; j < e->k; j++) if (getb(e->i[j]) < 0) setb(e->i[j], e->b[j]); }
static int collatz_bad(u64 n) { u64 v = n; for (int s = 1; s <= L; s++) { v = Tc(v); if (v < n) return 0; } return 1; }

/* barrier rescue at scale W with frozen bits; returns 1 and fills e on success */
static int nfl_last;
static int barrier(u64 n, int W, Ext *e) {
    long double thr = (long double)n * powl(3.0L, (long double)(W - 1));
    u64 x = n; int prevup = 0, nfl = 0; e->k = 0;
    for (int s = 1; s <= L; s++) {
        u64 i = pairof(x); int b = getb(i); if (b < 0) b = extget(e, i);
        u64 xn;
        if (b >= 0) { xn = stepb(x, b); prevup = ((x & 1) && b == 0); }
        else {
            int flip = 0;
            if ((x & 1) && prevup && (x % 4 == 3) && (long double)x >= thr && nfl < JMAX) {
                u64 tx = (3 * x + 1) / 2, it = pairof(tx); int bt = getb(it); if (bt < 0) bt = extget(e, it);
                if (bt != 1 && it != i) flip = 1;
            }
            if (flip) { extset(e, i, 1); u64 tx = (3 * x + 1) / 2, it = pairof(tx); if (isfree(e, it)) extset(e, it, 0); nfl++; xn = (x - 1) / 2; prevup = 0; }
            else { extset(e, i, 0); xn = Tc(x); prevup = (x & 1); }
        }
        x = xn;
        if (x < n) { nfl_last = nfl; return 1; }
        if (e->k > 190) return 0;
    }
    return 0;
}
/* private (ideal) barrier: pure Collatz + own barrier rule, no frozen bits; ideal_cost = sum over flips of n/(pair index) */
static double ideal_cost;
static int ideal_barrier(u64 n, int W) {
    long double thr = (long double)n * powl(3.0L, (long double)(W - 1));
    u64 x = n; int prevup = 0, nfl = 0; ideal_cost = 0.0;
    for (int s = 1; s <= L; s++) {
        int flip = ((x & 1) && prevup && (x % 4 == 3) && (long double)x >= thr && nfl < JMAX);
        if (flip) { ideal_cost += (double)n / (double)((x + 1) / 2); x = (x - 1) / 2; prevup = 0; nfl++; } else { prevup = (x & 1); x = Tc(x); }
        if (x < n) return 1;
    }
    return 0;
}
static int dg_type;   /* first deviation of the actual barrier path from the private one: 1 blocked, 2 high flip, 3 low flip/partner */
static void diagnose(u64 n, int W) {
    long double thr = (long double)n * powl(3.0L, (long double)(W - 1));
    u64 x = n; int prevup = 0, nfl = 0; dg_type = 0;
    for (int s = 1; s <= L; s++) {
        int flip = ((x & 1) && prevup && (x % 4 == 3) && (long double)x >= thr && nfl < JMAX);
        int b = getb(pairof(x));
        if (b >= 0 && b != (flip ? 1 : 0)) { dg_type = (flip ? 1 : (((x & 1) && (x - 1) / 2 >= n) ? 2 : 3)); return; }
        if (flip) { x = (x - 1) / 2; prevup = 0; nfl++; } else { prevup = (x & 1); x = Tc(x); }
        if (x < n) return;
    }
}
/* flip owner info (pair index -> owner, owner's time at the flip point) for rider statistics */
static u64 *fkey, *fown; static int *ftim; static u64 FCAP = 1ULL << 22, FCNT = 0;
static void finfo_set(u64 i, u64 own, int t) { if (FCNT * 10 > FCAP * 7) return; u64 h = hsh(i) & (FCAP - 1); while (fkey[h]) { if (fkey[h] == i) { fown[h] = own; ftim[h] = t; return; } h = (h + 1) & (FCAP - 1); } fkey[h] = i; fown[h] = own; ftim[h] = t; FCNT++; }
static int finfo_get(u64 i, u64 *own, int *t) { u64 h = hsh(i) & (FCAP - 1); while (fkey[h]) { if (fkey[h] == i) { *own = fown[h]; *t = ftim[h]; return 1; } h = (h + 1) & (FCAP - 1); } return 0; }
static void record_flips(u64 n) { u64 v = n; for (int s = 0; s < L; s++) { u64 i = pairof(v); int b = getb(i); if (b == 1 && (v & 1)) finfo_set(i, n, s); v = stepb(v, b < 0 ? 0 : b); if (v < n) break; } }
static long long rid_hi = 0, rid_early = 0, rid_late = 0, rid_late_fail = 0;
static void rider_scan(u64 n) {
    u64 v = n;
    for (int s = 0; s < L; s++) {
        u64 i = pairof(v); int b = getb(i); if (b < 0) b = 0;
        if (b == 1 && (v & 1) && (v - 1) / 2 >= n) { u64 own; int ot; if (finfo_get(i, &own, &ot)) { rid_hi++; if (s <= ot) rid_early++; else { rid_late++; if (!run(n, NULL)) rid_late_fail++; } } return; }
        v = stepb(v, b); if (v < n) return;
    }
}
/* THM-4475 fallback options + DFS (as in procgen_price_20260925_greedy.c, offset 0, q = 3) */
static Ext best; static int bestf; static int MAXF = 3;
static void dfs(u64 n, u64 v, int depth, Ext *e, int nf) {
    if (bestf >= 0 && nf >= bestf) return; if (depth == L) return;
    u64 i = pairof(v); int fb = getb(i), eb = extget(e, i); int opts[2], no = 0;
    if (fb >= 0) opts[no++] = fb; else if (eb >= 0) opts[no++] = eb; else { opts[no++] = 0; opts[no++] = 1; }
    for (int t = 0; t < no; t++) {
        int b = opts[t]; int fresh = (fb < 0 && eb < 0); int nf2 = nf + (fresh && b == 1); if (nf2 > MAXF) continue;
        if (fresh) extset(e, i, b);
        u64 w = stepb(v, b);
        if (w < n) { if (bestf < 0 || nf2 < bestf) { bestf = nf2; best = *e; } } else if (e->k < 40) dfs(n, w, depth + 1, e, nf2);
        if (fresh) e->k--;
    }
}
int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage: greedy L N WMAX MODE [JMAX]\n"); return 1; }
    L = atoi(argv[1]); N = strtoull(argv[2], 0, 10); WMAX = atoi(argv[3]); MODE = atoi(argv[4]); if (argc > 5) JMAX = atoi(argv[5]);
    DENSE = 16 * N + 64; dense = malloc(DENSE); memset(dense, -1, DENSE);
    HCAP = 1ULL << 23; hkey = calloc(HCAP, sizeof(u64)); hval = calloc(HCAP, 1);
    fkey = calloc(FCAP, 8); fown = calloc(FCAP, 8); ftim = calloc(FCAP, sizeof(int));
    long long nres = 0, nbar = 0, flipsbar = 0, fb[5] = {0}, stuck = 0, harmed = 0;
    long long byW[64] = {0}, ideal_hist[64] = {0}, degraded = 0, dev[4] = {0}, ideal_low = 0; double ideal_sum = 0.0;
    for (u64 n = 3; n <= N; n++) {
        rider_scan(n);
        if (run(n, NULL)) { if (MODE != 1 || collatz_bad(n)) freeze_used(); continue; }
        nres++; if (!collatz_bad(n)) harmed++;
        int wi = 0; for (int W = WMAX; W >= 2; W--) if (ideal_barrier(n, W)) { wi = W; break; }
        ideal_hist[wi]++; if (wi < 2) { ideal_low++; ideal_sum += 2.0; } else ideal_sum += ideal_cost;
        int done = 0, wa = 0; Ext e;
        if (MODE != 2) {
            if (wi >= 2) diagnose(n, wi);
            for (int W = WMAX; W >= 2 && !done; W--) if (barrier(n, W, &e)) { commit(&e); record_flips(n); done = 1; wa = W; nbar++; byW[W]++; flipsbar += nfl_last; }
            if (wa < wi) { degraded++; if (dg_type >= 1 && dg_type <= 3) dev[dg_type]++; else dev[0]++; }
        }
        if (done) continue;
        Ext cand[4]; int nc = 0, tag[4];
        Ext eA; eA.k = 0; int okA = 0; if (n & 1) { u64 i = pairof(n); if (isfree(NULL, i)) { extset(&eA, i, 1); okA = 1; } }
        Ext eB; eB.k = 0; int okB = 0;
        { int nb = bitat(NULL, pairof(n)); int goesup = ((n & 1) ^ nb) & 1; u64 r = goesup ? ((n & 1) ? (3 * n + 1) / 2 : 3 * (n / 2)) : 0;
          if (goesup && (r & 1)) { u64 i = pairof(r); if (isfree(NULL, i)) { extset(&eB, i, 1); okB = 1; } } }
        Ext eF; eF.k = 0; int okF = 0;
        { u64 v = n; for (int s = 0; s < L; s++) { u64 i = pairof(v); int b = bitat(NULL, i); u64 w = stepb(v, b); int upmove = ((v & 1) ^ b) & 1;
            if (upmove && (v & 1) && (w & 1) && w < 2 * n) { u64 pw = w + 1; u64 pu = 3 * (pw / 2); u64 iw = pairof(w); if (!(pu & 1) && isfree(NULL, iw)) { extset(&eF, iw, 1); okF = 1; break; } }
            v = w; } }
        int preferA = (n & 1) && (n % 8 == 3);
        if (n & 1) { if (preferA && okA) { cand[nc] = eA; tag[nc++] = 0; } if (okF) { cand[nc] = eF; tag[nc++] = 2; } if (okB) { cand[nc] = eB; tag[nc++] = 1; } if (!preferA && okA) { cand[nc] = eA; tag[nc++] = 0; } }
        else { if (okB) { cand[nc] = eB; tag[nc++] = 3; } if (okF) { cand[nc] = eF; tag[nc++] = 2; } }
        for (int c = 0; c < nc && !done; c++) if (run(n, &cand[c])) { freeze_used(); commit(&cand[c]); fb[tag[c]]++; done = 1; }
        if (done) continue;
        Ext e2; e2.k = 0; bestf = -1; dfs(n, n, 0, &e2, 0);
        if (bestf >= 0) { if (!run(n, &best)) { fprintf(stderr, "internal\n"); return 3; } freeze_used(); commit(&best); fb[4]++; continue; }
        stuck++;
    }
    u64 H = N / 2, fl = 0; for (u64 i = 1; i <= H; i++) if (getb(i) == 1) fl++;
    long long fail = 0; u64 firstfail = 0;
    for (u64 n = 3; n <= N; n++) { u64 v = n; int ok = 0; for (int t = 1; t <= L; t++) { v = stepb(v, bitat(NULL, pairof(v))); if (v < n) { ok = 1; break; } } if (!ok) { fail++; if (!firstfail) firstfail = n; } }
    printf("RESULT L=%d N=%llu WMAX=%d MODE=%d JMAX=%d density=%.8f flipped=%llu rescues=%lld barrier=%lld barrierflips=%lld A=%lld B=%lld F=%lld P=%lld dfs=%lld stuck=%lld harmed=%lld verifyfail=%lld firstfail=%llu hash=%llu\n",
           L, N, WMAX, MODE, JMAX, (double)fl / H, fl, nres, nbar, flipsbar, fb[0], fb[1], fb[2], fb[3], fb[4], stuck, harmed, fail, firstfail, HCNT);
    printf("SCALES"); for (int W = WMAX; W >= 2; W--) printf(" %d:%lld", W, byW[W]); printf("\n");
    printf("IDEAL"); for (int W = WMAX; W >= 0; W--) if (ideal_hist[W]) printf(" %d:%lld", W, ideal_hist[W]); printf("\n");
    printf("PRIVATE ideal_cost_density=%.8f\n", ideal_sum / (double)N);
    printf("INTERF degraded=%lld blocked=%lld highflip=%lld lowflip=%lld other=%lld ideal_needs_low=%lld riders=%lld early=%lld late=%lld late_fail=%lld\n",
           degraded, dev[1], dev[2], dev[3], dev[0], ideal_low, rid_hi, rid_early, rid_late, rid_late_fail);
    return 0;
}
