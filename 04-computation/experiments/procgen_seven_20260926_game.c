/*
 * procgen_seven_20260926_game.c -- memory-lean exact solver for the min-max cycle density rho*(q,k)
 * of the q n +- 1 strategy cube (THM-4486), lane "seven" of session collatz-procgen-20260922
 * (2026-09-26).  Standalone program; nothing it prints is used as a proof step: every value it
 * reports is backed by two certificate files that the independent checker
 * procgen_seven_20260926_verify.c re-checks edge by edge with exact integer arithmetic.
 *
 * PAIR FORM OF THE GAME.  Level k, odd q, N = 2^k, H = 2^(k-1).  Pairs P in Z/H; the pair P has the
 * two lifts x_b = P + b H (b = 0, 1), both of the parity of P.  The options of a node x are the pairs
 *     x even : x/2 mod H,        x odd : (q x + s)/2 mod H,  s = +1 or -1  (Min's sign).
 * Max (the adversary) picks the lift b at every pair, Min picks the sign at every odd node.
 * Threshold F = fn/fd; weight e(P) = fd - fn (P odd), -fn (P even).
 *
 *   MAX energy (lower bound rho* >= F):  g(P) = max(0, -e(P) + min_b max_{Q option of x_b} g(Q)).
 *     If W = {g < TOP} is nonempty then, with tau(P) = argmin_b, every option Q of x_tau(P) lies in W
 *     and g(Q) <= g(P) + e(P): every sign strategy has a cycle of density >= F  (lower certificate).
 *   MIN energy (upper bound rho* <= F):  G(P) = max(0, e(P) + max_b min_{Q option of x_b} G(Q)).
 *     If G < TOP everywhere then, with sigma(x) = argmin sign, G(T_sigma(x) mod H) + e(x) <= G(x mod H)
 *     for every node x: every cycle of G_sigma has density <= F  (upper certificate).
 *
 * SYMMETRY.  Negation nu : x -> -x maps the option set of x onto minus the option set of -x and the
 * lift b of P onto the lift 1-b of -P (P != 0, H/2 are paired with H-P).  Both operators commute with
 * nu, so their least fixed points from 0 are nu-invariant and are computed on the orbit
 * representatives r = min(P, H-P), 0 <= r <= H/2 (half the memory and work).  Certificates are
 * written on all H pairs (value of P = value of rep(P)); the checker does not use the symmetry.
 *
 * SEARCH.  Stern-Brocot search on F between Farey neighbours; mediants above a known upper bound are
 * skipped.  At each mediant the two energy tests are dovetailed with doubling work budgets (the
 * winning side usually converges fast; the losing side climbs to the cap).  A completed MAX success
 * whose tight graph has no cycle proves rho* > F (no Min response reaches density F against that tau);
 * a completed MIN success whose tight graph has no cycle proves rho* < F.  Warm starts: MAX from the
 * last completed MAX fixed point at a smaller threshold, MIN from the last completed MIN fixed point
 * at a larger threshold, scaled (both least fixed points are monotone and positively homogeneous in
 * the weights).  Caps only turn finite values into TOP (never the reverse), so a small cap can only
 * make the search fail -- detected (mediant denominator > N) and retried with larger caps.
 *
 * usage:  procgen_seven_20260926_game value q k lon lod hin hid kun kud [prefix]
 *   (initial Farey bracket lon/lod < hin/hid; mediants above kun/kud (> 0) are skipped: pass 0 0 for none)
 *   prints "VALUE fn fd", statistics, and writes (if prefix given)
 *     prefix.lo_tau  uint8[H]   lift bit tau(P) (lower certificate)
 *     prefix.lo_g    int32[H]   potential g(P), -1 = not in W
 *     prefix.up_sig  uint8[H]   1 iff the sign of the odd node 2i+1 is '-' (index i)
 *     prefix.up_psi  int32[H]   potential G(P)
 *         procgen_seven_20260926_game tauval q k taufile   (value of a frozen Max lift strategy;
 *   exploratory: prints the least density Min can force against tau from its worst start, i.e.
 *   max over starts of the min reachable cycle density; Stern-Brocot with one-player energy tests)
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

typedef int32_t i32;
typedef int64_t i64;
#define TOP ((i32)0x7fffffff)

static int K;
static i64 N, H, HH, QQ, QINV;
static int USEQ = 1;        /* 1: work on nu-orbit representatives 0..H/2; 0: on all pairs 0..H-1 */
static i64 NR;              /* largest representative: H/2 (USEQ) or H-1 */

static uint8_t *TAUFIX = 0; /* optional frozen Max strategy (tauval mode), per pair */

static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }

static i64 inv_mod_pow2(i64 q) {
    uint64_t x = (uint64_t)q;             /* Newton: x <- x (2 - q x), correct mod 2^64 after 6 steps */
    for (int i = 0; i < 6; i++) x = x * (2ULL - (uint64_t)q * x);
    return (i64)(x & (uint64_t)(N - 1));
}

static inline i64 rep(i64 P) { return (!USEQ || P <= HH) ? P : H - P; }

/* the options of the lift b of pair P: returns the number of options (1 or 2) in Q[] */
static inline int options(i64 P, int b, i64 *Qo) {
    i64 x = P + (b ? H : 0);
    if (!(P & 1)) { Qo[0] = (x >> 1) & (H - 1); return 1; }
    i64 t = QQ * x;
    Qo[0] = ((t + 1) >> 1) & (H - 1);    /* sign + */
    Qo[1] = ((t - 1) >> 1) & (H - 1);    /* sign - */
    return 2;
}

typedef struct {
    int mode;            /* 1 = MAX energy (lower bound), 0 = MIN energy (upper bound) */
    i32 *g;              /* values on reps 0..HH */
    i32 *qu;             /* circular queue of reps */
    uint8_t *inq;
    i64 qh, qt, qn;      /* queue head, tail, count */
    i64 cap;
    i64 work;
    int done;
    i64 fn, fd;
} State;

static void st_alloc(State *S, int mode) {
    S->mode = mode;
    S->g = (i32 *)malloc(sizeof(i32) * (NR + 1));
    S->qu = (i32 *)malloc(sizeof(i32) * (NR + 2));
    S->inq = (uint8_t *)malloc(NR + 1);
    if (!S->g || !S->qu || !S->inq) { fprintf(stderr, "out of memory\n"); exit(2); }
}

/* start a fixed-point computation at threshold fn/fd with cap; if warm != NULL, start from
   ceil(warm * fd / wfd) (TOP stays TOP), else from 0 */
static void st_start(State *S, i64 fn, i64 fd, i64 cap, const i32 *warm, i64 wfd) {
    S->fn = fn; S->fd = fd; S->cap = cap; S->work = 0; S->done = 0;
    for (i64 r = 0; r <= NR; r++) {
        i32 v = 0;
        if (warm) {
            if (warm[r] == TOP) v = TOP;
            else {
                i64 w = ((i64)warm[r] * fd + wfd - 1) / wfd;
                v = w > cap ? TOP : (i32)w;
            }
        }
        S->g[r] = v;
        S->qu[r] = (i32)r;
        S->inq[r] = 1;
    }
    S->qh = 0; S->qt = (NR + 1) % (NR + 2); S->qn = NR + 1;
}

static inline i32 eval_rep(const State *S, i64 P) {
    const i32 *g = S->g;
    i64 e = (P & 1) ? (S->fd - S->fn) : -S->fn;
    i64 Qo[2];
    if (S->mode == 1) {                              /* MAX energy: min over b of max over options */
        i64 best = TOP;
        for (int b = 0; b < 2; b++) {
            if (TAUFIX && TAUFIX[P] != (uint8_t)b && TAUFIX[P] != 2) continue;
            int no = options(P, b, Qo);
            i64 v = -1;
            for (int j = 0; j < no; j++) { i32 x = g[rep(Qo[j])]; if (x > v) v = x; }
            if (v < best) best = v;
        }
        if (best >= TOP) return TOP;
        i64 v = best - e;
        if (v < 0) v = 0;
        return v > S->cap ? TOP : (i32)v;
    } else {                                         /* MIN energy: max over b of min over options */
        i64 best = -1;
        for (int b = 0; b < 2; b++) {
            int no = options(P, b, Qo);
            i64 v = TOP;
            for (int j = 0; j < no; j++) { i32 x = g[rep(Qo[j])]; if (x < v) v = x; }
            if (v > best) best = v;
        }
        if (best >= TOP) return TOP;
        i64 v = best + e;
        if (v < 0) v = 0;
        return v > S->cap ? TOP : (i32)v;
    }
}

/* run the worklist until the queue is empty (returns 1) or the work budget is used (returns 0) */
static int st_run(State *S, i64 budget) {
    i64 stop = S->work + budget;
    while (S->qn > 0) {
        if (S->work >= stop) return 0;
        i64 r = S->qu[S->qh]; S->qh++; if (S->qh == NR + 2) S->qh = 0; S->qn--; S->inq[r] = 0;
        S->work++;
        if (S->g[r] == TOP) continue;
        i32 v = eval_rep(S, r);
        if (v > S->g[r]) {
            S->g[r] = v;
            i64 pr[3];
            pr[0] = rep((2 * r) & (H - 1));
            pr[1] = rep((QINV * ((2 * r - 1) & (H - 1))) & (H - 1));
            pr[2] = rep((QINV * ((2 * r + 1) & (H - 1))) & (H - 1));
            for (int j = 0; j < 3; j++) {
                i64 w = pr[j];
                if (!S->inq[w] && S->g[w] != TOP) {
                    S->inq[w] = 1; S->qu[S->qt] = (i32)w; S->qt++; if (S->qt == NR + 2) S->qt = 0; S->qn++;
                }
            }
        }
    }
    S->done = 1;
    return 1;
}

static i64 count_top(const State *S) {
    i64 c = 0;
    for (i64 r = 0; r <= NR; r++) if (S->g[r] == TOP) c++;
    return c;
}

/* value of a pair under the state's fixed point */
static inline i32 gval(const State *S, i64 P) { return S->g[rep(P)]; }

/* MAX certificate lift choice at pair P (argmin over b of the max over options; ties -> 0) */
static int max_tau(const State *S, i64 P) {
    i64 Qo[2]; i64 bv[2];
    for (int b = 0; b < 2; b++) {
        int no = options(P, b, Qo);
        i64 v = -1;
        for (int j = 0; j < no; j++) { i32 x = gval(S, Qo[j]); if (x > v) v = x; }
        bv[b] = v;
    }
    if (TAUFIX && TAUFIX[P] != 2) return TAUFIX[P];
    return bv[1] < bv[0] ? 1 : 0;
}

/* MIN certificate sign at the odd node x (1 = sign '-'; ties -> '+') */
static int min_sig(const State *S, i64 x) {
    i64 t = QQ * x;
    i64 Qp = ((t + 1) >> 1) & (H - 1), Qm = ((t - 1) >> 1) & (H - 1);
    return gval(S, Qm) < gval(S, Qp) ? 1 : 0;
}

/* Does the tight graph of a completed certificate contain a cycle?
   MAX (lower) certificate: nodes = pairs in W, edges P -> Q for the options Q of x_tau(P) with
   g(Q) == g(P) + e(P).  A cycle of density exactly F against tau exists iff this graph has a cycle.
   MIN (upper) certificate: nodes = pairs, edges P -> Q (Q the sigma-option of a lift x_b of P) with
   G(Q) + e(P) == G(P).  A cycle of G_sigma of density exactly F exists iff this graph has a cycle.
   Peeling (Kahn) on out-degrees; the predecessors of Q are among 2Q, q^-1 (2Q -+ 1) mod H. */
static int tight_edges(const State *S, i64 P, i64 *out) {
    i64 e = (P & 1) ? (S->fd - S->fn) : -S->fn;
    i32 gp = gval(S, P);
    if (gp == TOP) return 0;
    int n = 0;
    if (S->mode == 1) {
        int b = max_tau(S, P);
        i64 Qo[2];
        int no = options(P, b, Qo);
        for (int j = 0; j < no; j++) { i32 gq = gval(S, Qo[j]); if (gq != TOP && (i64)gq == (i64)gp + e) out[n++] = Qo[j]; }
    } else {
        for (int b = 0; b < 2; b++) {
            i64 x = P + (b ? H : 0);
            i64 Q;
            if (!(P & 1)) Q = (x >> 1) & (H - 1);
            else {
                i64 t = QQ * x;
                Q = min_sig(S, x) ? (((t - 1) >> 1) & (H - 1)) : (((t + 1) >> 1) & (H - 1));
            }
            i32 gq = gval(S, Q);
            if ((i64)gq + e == (i64)gp) {
                int dup = 0;
                for (int j = 0; j < n; j++) if (out[j] == Q) dup = 1;
                if (!dup) out[n++] = Q;
            }
        }
    }
    return n;
}

static int tight_has_cycle(const State *S, i64 *remaining) {
    uint8_t *deg = (uint8_t *)malloc(H);
    i32 *stk = (i32 *)malloc(sizeof(i32) * H);
    if (!deg || !stk) { fprintf(stderr, "out of memory (tight)\n"); exit(2); }
    i64 sp = 0, alive = 0;
    i64 out[4];
    for (i64 P = 0; P < H; P++) {
        if (gval(S, P) == TOP) { deg[P] = 255; continue; }
        alive++;
        deg[P] = (uint8_t)tight_edges(S, P, out);
        if (deg[P] == 0) stk[sp++] = (i32)P;
    }
    while (sp > 0) {
        i64 Q = stk[--sp];
        alive--;
        i64 pr[3];
        pr[0] = (2 * Q) & (H - 1);
        pr[1] = (QINV * ((2 * Q - 1) & (H - 1))) & (H - 1);
        pr[2] = (QINV * ((2 * Q + 1) & (H - 1))) & (H - 1);
        for (int j = 0; j < 3; j++) {
            i64 P = pr[j];
            if (deg[P] == 255 || deg[P] == 0) continue;
            int n = tight_edges(S, P, out);
            for (int i = 0; i < n; i++) if (out[i] == Q) { deg[P]--; if (deg[P] == 0) stk[sp++] = (i32)P; }
        }
    }
    free(deg); free(stk);
    if (remaining) *remaining = alive;
    return alive > 0;
}

static void write_certs(const char *prefix, const State *SX, const State *SN) {
    char fn[1024];
    FILE *f;
    const i64 CH = 1 << 20;
    uint8_t *b8 = (uint8_t *)malloc(CH);
    i32 *b32 = (i32 *)malloc(sizeof(i32) * CH);
    snprintf(fn, sizeof fn, "%s.lo_tau", prefix); f = fopen(fn, "wb");
    for (i64 lo = 0; lo < H; lo += CH) { i64 n = H - lo < CH ? H - lo : CH;
        for (i64 i = 0; i < n; i++) b8[i] = gval(SX, lo + i) == TOP ? 0 : (uint8_t)max_tau(SX, lo + i);
        fwrite(b8, 1, n, f); }
    fclose(f);
    snprintf(fn, sizeof fn, "%s.lo_g", prefix); f = fopen(fn, "wb");
    for (i64 lo = 0; lo < H; lo += CH) { i64 n = H - lo < CH ? H - lo : CH;
        for (i64 i = 0; i < n; i++) { i32 v = gval(SX, lo + i); b32[i] = v == TOP ? -1 : v; }
        fwrite(b32, 4, n, f); }
    fclose(f);
    snprintf(fn, sizeof fn, "%s.up_sig", prefix); f = fopen(fn, "wb");
    for (i64 lo = 0; lo < H; lo += CH) { i64 n = H - lo < CH ? H - lo : CH;
        for (i64 i = 0; i < n; i++) b8[i] = (uint8_t)min_sig(SN, 2 * (lo + i) + 1);
        fwrite(b8, 1, n, f); }
    fclose(f);
    snprintf(fn, sizeof fn, "%s.up_psi", prefix); f = fopen(fn, "wb");
    for (i64 lo = 0; lo < H; lo += CH) { i64 n = H - lo < CH ? H - lo : CH;
        for (i64 i = 0; i < n; i++) b32[i] = gval(SN, lo + i);
        fwrite(b32, 4, n, f); }
    fclose(f);
    free(b8); free(b32);
}

static void init_level(int k, i64 q) {
    K = k; N = (i64)1 << k; H = N >> 1; HH = H >> 1; QQ = q; QINV = inv_mod_pow2(q);
    NR = USEQ ? HH : H - 1;
    if (((QQ * QINV) & (N - 1)) != 1) { fprintf(stderr, "inverse failed\n"); exit(2); }
}

/* ------------------------------------------------------------------------------------------------ */
static int value_mode(i64 q, int k, i64 lon0, i64 lod0, i64 hin0, i64 hid0, i64 kun, i64 kud, const char *prefix) {
    init_level(k, q);
    double t0 = now();
    State SX, SN;
    st_alloc(&SX, 1); st_alloc(&SN, 0);
    i32 *saveX = (i32 *)malloc(sizeof(i32) * (NR + 1));   /* completed MAX fixed point at lo */
    i32 *saveN = (i32 *)malloc(sizeof(i32) * (NR + 1));   /* completed MIN fixed point at hi */
    if (!saveX || !saveN) { fprintf(stderr, "out of memory\n"); return 2; }
    i64 capmul = 1;
    for (int attempt = 0; attempt < 4; attempt++, capmul *= 8) {
        i64 lon = lon0, lod = lod0, hin = hin0, hid = hid0;
        if (hin * lod - lon * hid != 1) { printf("initial bracket is not a pair of Farey neighbours\n"); return 2; }
        int haveX = 0, haveN = 0; i64 sXd = 1, sNd = 1;
        int tests = 0;
        while (1) {
            i64 mn = lon + hin, md = lod + hid;
            if (md > N) break;
            if (kud > 0 && mn * kud > kun * md) {          /* mediant above the known upper bound */
                hin = mn; hid = md;
                printf("  mediant %lld/%lld: above the known upper bound\n", (long long)mn, (long long)md);
                continue;
            }
            i64 capX = capmul * (12 * md + 64), capN = capmul * (6 * md + 64);
            st_start(&SX, mn, md, capX, haveX ? saveX : 0, sXd);
            st_start(&SN, mn, md, capN, haveN ? saveN : 0, sNd);
            i64 budget = 4 * (NR + 1);
            int verdict = 0;   /* +1: rho* > m ; -1: rho* < m ; 2: rho* = m */
            const char *why = "";
            while (!verdict) {
                if (!SX.done) st_run(&SX, budget);
                if (SX.done && !verdict) {
                    i64 top = count_top(&SX);
                    if (top == NR + 1) { verdict = -1; why = "MAX energy lost everywhere"; break; }
                    i64 rem = 0;
                    if (!tight_has_cycle(&SX, &rem)) { verdict = 1; why = "MAX won, tight graph acyclic"; break; }
                    /* need the MIN test to completion */
                    st_run(&SN, (i64)1 << 62);
                    if (count_top(&SN) == 0) { verdict = 2; why = "both certificates"; }
                    else { verdict = 1; why = "MAX won (tight cycle), MIN lost"; }
                    break;
                }
                if (!SN.done) st_run(&SN, budget);
                if (SN.done && !verdict) {
                    if (count_top(&SN) > 0) { verdict = 1; why = "MIN energy lost somewhere"; break; }
                    i64 rem = 0;
                    if (!tight_has_cycle(&SN, &rem)) { verdict = -1; why = "MIN won, tight graph acyclic"; break; }
                    st_run(&SX, (i64)1 << 62);
                    if (count_top(&SX) < NR + 1) { verdict = 2; why = "both certificates"; }
                    else { verdict = -1; why = "MIN won (tight cycle), MAX lost"; }
                    break;
                }
                budget *= 2;
            }
            tests++;
            printf("  mediant %lld/%lld: %s (%s; work MAX %lld, MIN %lld; %.1f s)\n", (long long)mn, (long long)md,
                   verdict == 2 ? "VALUE" : (verdict > 0 ? "rho* >" : "rho* <"), why,
                   (long long)SX.work, (long long)SN.work, now() - t0);
            fflush(stdout);
            if (verdict == 2) {
                printf("VALUE %lld %lld\n", (long long)mn, (long long)md);
                printf("STATS k=%d q=%lld tests=%d time=%.2f capmul=%lld\n", k, (long long)q, tests, now() - t0, (long long)capmul);
                if (prefix) write_certs(prefix, &SX, &SN);
                return 0;
            }
            if (verdict > 0) {
                lon = mn; lod = md;
                if (SX.done) { memcpy(saveX, SX.g, sizeof(i32) * (NR + 1)); haveX = 1; sXd = md; }
            } else {
                hin = mn; hid = md;
                if (SN.done) { memcpy(saveN, SN.g, sizeof(i32) * (NR + 1)); haveN = 1; sNd = md; }
            }
        }
        printf("  search failed with cap multiplier %lld; retrying\n", (long long)capmul);
    }
    printf("FAILED\n");
    return 1;
}


/* ------------------------------------------------------------------------------------------------
   LEAN mode (large k): one threshold F, one fixed point at a time, uint16 values on the nu-orbit
   representatives and a dirty bitset processed by in-order sweeps (2.125 bytes per representative).
   Writes the lower certificate (tau bits, g as uint8/uint16) and then the upper certificate
   (sigma bits, psi as uint8/uint16); the checker reads these compact formats.                      */
typedef uint16_t u16;
#define TOP16 ((u16)0xffff)
static u16 *LV;            /* values on reps */
static uint64_t *DIRTY;    /* dirty bits on reps */
static int LMODE;          /* 1 = MAX energy, 0 = MIN energy */
static i64 LFN, LFD, LCAP;

static inline u16 lval(i64 P) { return LV[rep(P)]; }

static inline u16 leval(i64 P) {
    i64 e = (P & 1) ? (LFD - LFN) : -LFN;
    i64 Qo[2];
    if (LMODE == 1) {
        i64 best = TOP16;
        for (int b = 0; b < 2; b++) {
            int no = options(P, b, Qo);
            i64 v = -1;
            for (int j = 0; j < no; j++) { u16 x = LV[rep(Qo[j])]; if ((i64)x > v) v = x; }
            if (v < best) best = v;
        }
        if (best >= TOP16) return TOP16;
        i64 v = best - e; if (v < 0) v = 0;
        return v > LCAP ? TOP16 : (u16)v;
    } else {
        i64 best = -1;
        for (int b = 0; b < 2; b++) {
            int no = options(P, b, Qo);
            i64 v = TOP16;
            for (int j = 0; j < no; j++) { u16 x = LV[rep(Qo[j])]; if ((i64)x < v) v = x; }
            if (v > best) best = v;
        }
        if (best >= TOP16) return TOP16;
        i64 v = best + e; if (v < 0) v = 0;
        return v > LCAP ? TOP16 : (u16)v;
    }
}

static i64 lsolve(int mode, i64 fn, i64 fd, i64 cap, int *sweeps_out) {
    LMODE = mode; LFN = fn; LFD = fd; LCAP = cap;
    i64 nw = (NR + 1 + 63) / 64;
    for (i64 r = 0; r <= NR; r++) LV[r] = 0;
    for (i64 w = 0; w < nw; w++) DIRTY[w] = ~0ULL;
    if ((NR + 1) % 64) DIRTY[nw - 1] = (1ULL << ((NR + 1) % 64)) - 1;
    i64 work = 0; int sweeps = 0, any = 1;
    while (any) {
        any = 0; sweeps++;
        for (i64 w = 0; w < nw; w++) {
            while (DIRTY[w]) {
                int bit = __builtin_ctzll(DIRTY[w]);
                DIRTY[w] &= DIRTY[w] - 1;
                i64 r = w * 64 + bit;
                work++;
                if (LV[r] == TOP16) continue;
                u16 v = leval(r);
                if (v > LV[r]) {
                    LV[r] = v; any = 1;
                    i64 pr[3];
                    pr[0] = rep((2 * r) & (H - 1));
                    pr[1] = rep((QINV * ((2 * r - 1) & (H - 1))) & (H - 1));
                    pr[2] = rep((QINV * ((2 * r + 1) & (H - 1))) & (H - 1));
                    for (int j = 0; j < 3; j++) if (LV[pr[j]] != TOP16) DIRTY[pr[j] >> 6] |= 1ULL << (pr[j] & 63);
                }
            }
        }
    }
    if (sweeps_out) *sweeps_out = sweeps;
    return work;
}

static int lean_mode(i64 q, int k, i64 fn, i64 fd, const char *prefix) {
    init_level(k, q);
    double t0 = now();
    LV = (u16 *)malloc(sizeof(u16) * (NR + 1));
    DIRTY = (uint64_t *)malloc(sizeof(uint64_t) * ((NR + 1 + 63) / 64));
    if (!LV || !DIRTY) { fprintf(stderr, "out of memory\n"); return 2; }
    char fnm[1024]; FILE *f;
    const i64 CH = 1 << 20;
    uint8_t *b8 = (uint8_t *)malloc(CH);
    u16 *b16 = (u16 *)malloc(2 * CH);
    int ok = 1;
    /* ---- lower (MAX) ---- */
    i64 cap = 12 * fd + 64; int sw = 0;
    if (cap > 60000) cap = 60000;
    i64 work = lsolve(1, fn, fd, cap, &sw);
    i64 ntop = 0, vmax = 0;
    for (i64 r = 0; r <= NR; r++) { if (LV[r] == TOP16) ntop++; else if (LV[r] > vmax) vmax = LV[r]; }
    printf("  LEAN MAX F=%lld/%lld cap %lld: work %lld, sweeps %d, TOP reps %lld of %lld, max g %lld, %.1f s\n",
           (long long)fn, (long long)fd, (long long)cap, (long long)work, sw, (long long)ntop, (long long)(NR + 1), (long long)vmax, now() - t0);
    fflush(stdout);
    if (ntop == NR + 1) { printf("LOWER-FAILED\n"); ok = 0; }
    else {
        int wide = vmax > 254;
        snprintf(fnm, sizeof fnm, "%s.lo_taub", prefix); f = fopen(fnm, "wb");
        for (i64 lo = 0; lo < H; lo += 8 * CH) {
            i64 n = H - lo < 8 * CH ? H - lo : 8 * CH;
            memset(b8, 0, (n + 7) / 8);
            for (i64 i = 0; i < n; i++) {
                i64 P = lo + i;
                if (lval(P) == TOP16) continue;
                i64 Qo[2]; i64 bv[2];
                for (int b = 0; b < 2; b++) {
                    int no = options(P, b, Qo); i64 v = -1;
                    for (int j = 0; j < no; j++) { u16 x = lval(Qo[j]); if ((i64)x > v) v = x; }
                    bv[b] = v;
                }
                if (bv[1] < bv[0]) b8[i >> 3] |= (uint8_t)(1u << (i & 7));
            }
            fwrite(b8, 1, (n + 7) / 8, f);
        }
        fclose(f);
        snprintf(fnm, sizeof fnm, "%s.%s", prefix, wide ? "lo_g16" : "lo_g8"); f = fopen(fnm, "wb");
        for (i64 lo = 0; lo < H; lo += CH) {
            i64 n = H - lo < CH ? H - lo : CH;
            if (wide) { for (i64 i = 0; i < n; i++) b16[i] = lval(lo + i); fwrite(b16, 2, n, f); }
            else { for (i64 i = 0; i < n; i++) { u16 v = lval(lo + i); b8[i] = v == TOP16 ? 255 : (uint8_t)v; } fwrite(b8, 1, n, f); }
        }
        fclose(f);
    }
    /* ---- upper (MIN) ---- */
    cap = 6 * fd + 64;
    if (cap > 60000) cap = 60000;
    work = lsolve(0, fn, fd, cap, &sw);
    ntop = 0; vmax = 0;
    for (i64 r = 0; r <= NR; r++) { if (LV[r] == TOP16) ntop++; else if (LV[r] > vmax) vmax = LV[r]; }
    printf("  LEAN MIN F=%lld/%lld cap %lld: work %lld, sweeps %d, TOP reps %lld, max psi %lld, %.1f s\n",
           (long long)fn, (long long)fd, (long long)cap, (long long)work, sw, (long long)ntop, (long long)vmax, now() - t0);
    fflush(stdout);
    if (ntop > 0) { printf("UPPER-FAILED\n"); ok = 0; }
    else {
        int wide = vmax > 254;
        snprintf(fnm, sizeof fnm, "%s.up_sigb", prefix); f = fopen(fnm, "wb");
        for (i64 lo = 0; lo < H; lo += 8 * CH) {
            i64 n = H - lo < 8 * CH ? H - lo : 8 * CH;
            memset(b8, 0, (n + 7) / 8);
            for (i64 i = 0; i < n; i++) {
                i64 x = 2 * (lo + i) + 1, t = QQ * x;
                i64 Qp = ((t + 1) >> 1) & (H - 1), Qm = ((t - 1) >> 1) & (H - 1);
                if (lval(Qm) < lval(Qp)) b8[i >> 3] |= (uint8_t)(1u << (i & 7));
            }
            fwrite(b8, 1, (n + 7) / 8, f);
        }
        fclose(f);
        snprintf(fnm, sizeof fnm, "%s.%s", prefix, wide ? "up_psi16" : "up_psi8"); f = fopen(fnm, "wb");
        for (i64 lo = 0; lo < H; lo += CH) {
            i64 n = H - lo < CH ? H - lo : CH;
            if (wide) { for (i64 i = 0; i < n; i++) b16[i] = lval(lo + i); fwrite(b16, 2, n, f); }
            else { for (i64 i = 0; i < n; i++) b8[i] = (uint8_t)lval(lo + i); fwrite(b8, 1, n, f); }
        }
        fclose(f);
    }
    printf("LEAN %s %lld %lld time=%.2f\n", ok ? "BOTH" : "INCOMPLETE", (long long)fn, (long long)fd, now() - t0);
    return ok ? 0 : 1;
}


/* ------------------------------------------------------------------------------------------------
   LEAN SEARCH (large k): the Stern-Brocot search of value_mode with the compact engine: two
   resumable sweep states (MAX and MIN, uint16 values + dirty bits) dovetailed with doubling budgets;
   no tight-cycle shortcut (so when one side wins with the answer "rho* >= m" or "rho* <= m" the other
   side is run to completion).  Certificates are written in the compact format at the value.       */
typedef struct {
    int mode; i64 fn, fd, cap;
    u16 *v; uint64_t *dirty;
    i64 w;          /* current word of the current sweep */
    int anysweep;   /* a value changed in the current sweep */
    int done; i64 work; int sweeps;
    int stopfirst;  /* MIN: stop at the first representative that reaches TOP (verdict "lost") */
    int lostfirst;  /* set when stopped that way */
} LState;

static void ls_alloc(LState *S, int mode) {
    S->mode = mode;
    S->v = (u16 *)malloc(sizeof(u16) * (NR + 1));
    S->dirty = (uint64_t *)malloc(sizeof(uint64_t) * ((NR + 1 + 63) / 64));
    if (!S->v || !S->dirty) { fprintf(stderr, "out of memory\n"); exit(2); }
}

static void ls_start(LState *S, i64 fn, i64 fd, i64 cap) {
    S->fn = fn; S->fd = fd; S->cap = cap; S->w = 0; S->anysweep = 0; S->done = 0; S->work = 0; S->sweeps = 1;
    S->stopfirst = 0; S->lostfirst = 0;
    i64 nw = (NR + 1 + 63) / 64;
    for (i64 r = 0; r <= NR; r++) S->v[r] = 0;
    for (i64 w = 0; w < nw; w++) S->dirty[w] = ~0ULL;
    if ((NR + 1) % 64) S->dirty[nw - 1] = (1ULL << ((NR + 1) % 64)) - 1;
}

static int ls_run(LState *S, i64 budget) {
    i64 nw = (NR + 1 + 63) / 64;
    i64 stop = S->work + budget;
    LV = S->v; LMODE = S->mode; LFN = S->fn; LFD = S->fd; LCAP = S->cap;
    while (1) {
        for (; S->w < nw; S->w++) {
            while (S->dirty[S->w]) {
                if (S->work >= stop) return 0;
                int bit = __builtin_ctzll(S->dirty[S->w]);
                S->dirty[S->w] &= S->dirty[S->w] - 1;
                i64 r = S->w * 64 + bit;
                S->work++;
                if (LV[r] == TOP16) continue;
                u16 v = leval(r);
                if (v > LV[r]) {
                    LV[r] = v; S->anysweep = 1;
                    if (v == TOP16 && S->stopfirst) { S->done = 1; S->lostfirst = 1; return 1; }
                    i64 pr[3];
                    pr[0] = rep((2 * r) & (H - 1));
                    pr[1] = rep((QINV * ((2 * r - 1) & (H - 1))) & (H - 1));
                    pr[2] = rep((QINV * ((2 * r + 1) & (H - 1))) & (H - 1));
                    for (int j = 0; j < 3; j++) if (LV[pr[j]] != TOP16) S->dirty[pr[j] >> 6] |= 1ULL << (pr[j] & 63);
                }
            }
        }
        if (!S->anysweep) { S->done = 1; return 1; }
        S->w = 0; S->anysweep = 0; S->sweeps++;
    }
}

static i64 ls_top(const LState *S) { i64 c = 0; for (i64 r = 0; r <= NR; r++) if (S->v[r] == TOP16) c++; return c; }

static void lean_write(const char *prefix, LState *SX, LState *SN) {
    char fnm[1024]; FILE *f;
    const i64 CH = 1 << 20;
    uint8_t *b8 = (uint8_t *)malloc(CH);
    u16 *b16 = (u16 *)malloc(2 * CH);
    i64 vmax = 0;
    LV = SX->v;
    for (i64 r = 0; r <= NR; r++) if (LV[r] != TOP16 && LV[r] > vmax) vmax = LV[r];
    int wide = vmax > 254;
    snprintf(fnm, sizeof fnm, "%s.lo_taub", prefix); f = fopen(fnm, "wb");
    for (i64 lo = 0; lo < H; lo += 8 * CH) {
        i64 n = H - lo < 8 * CH ? H - lo : 8 * CH;
        memset(b8, 0, (n + 7) / 8);
        for (i64 i = 0; i < n; i++) {
            i64 P = lo + i;
            if (lval(P) == TOP16) continue;
            i64 Qo[2]; i64 bv[2];
            for (int b = 0; b < 2; b++) {
                int no = options(P, b, Qo); i64 v = -1;
                for (int j = 0; j < no; j++) { u16 x = lval(Qo[j]); if ((i64)x > v) v = x; }
                bv[b] = v;
            }
            if (bv[1] < bv[0]) b8[i >> 3] |= (uint8_t)(1u << (i & 7));
        }
        fwrite(b8, 1, (n + 7) / 8, f);
    }
    fclose(f);
    snprintf(fnm, sizeof fnm, "%s.%s", prefix, wide ? "lo_g16" : "lo_g8"); f = fopen(fnm, "wb");
    for (i64 lo = 0; lo < H; lo += CH) {
        i64 n = H - lo < CH ? H - lo : CH;
        if (wide) { for (i64 i = 0; i < n; i++) b16[i] = lval(lo + i); fwrite(b16, 2, n, f); }
        else { for (i64 i = 0; i < n; i++) { u16 v = lval(lo + i); b8[i] = v == TOP16 ? 255 : (uint8_t)v; } fwrite(b8, 1, n, f); }
    }
    fclose(f);
    LV = SN->v; vmax = 0;
    for (i64 r = 0; r <= NR; r++) if (LV[r] > vmax) vmax = LV[r];
    wide = vmax > 254;
    snprintf(fnm, sizeof fnm, "%s.up_sigb", prefix); f = fopen(fnm, "wb");
    for (i64 lo = 0; lo < H; lo += 8 * CH) {
        i64 n = H - lo < 8 * CH ? H - lo : 8 * CH;
        memset(b8, 0, (n + 7) / 8);
        for (i64 i = 0; i < n; i++) {
            i64 x = 2 * (lo + i) + 1, t = QQ * x;
            i64 Qp = ((t + 1) >> 1) & (H - 1), Qm = ((t - 1) >> 1) & (H - 1);
            if (lval(Qm) < lval(Qp)) b8[i >> 3] |= (uint8_t)(1u << (i & 7));
        }
        fwrite(b8, 1, (n + 7) / 8, f);
    }
    fclose(f);
    snprintf(fnm, sizeof fnm, "%s.%s", prefix, wide ? "up_psi16" : "up_psi8"); f = fopen(fnm, "wb");
    for (i64 lo = 0; lo < H; lo += CH) {
        i64 n = H - lo < CH ? H - lo : CH;
        if (wide) { for (i64 i = 0; i < n; i++) b16[i] = lval(lo + i); fwrite(b16, 2, n, f); }
        else { for (i64 i = 0; i < n; i++) b8[i] = (uint8_t)lval(lo + i); fwrite(b8, 1, n, f); }
    }
    fclose(f);
    free(b8); free(b16);
}

static int lsearch_mode(i64 q, int k, i64 lon0, i64 lod0, i64 hin0, i64 hid0, i64 kun, i64 kud, const char *prefix) {
    init_level(k, q);
    double t0 = now();
    LState SX, SN; ls_alloc(&SX, 1); ls_alloc(&SN, 0);
    i64 capmul = 1;
    for (int attempt = 0; attempt < 3; attempt++, capmul *= 4) {
        i64 lon = lon0, lod = lod0, hin = hin0, hid = hid0;
        if (hin * lod - lon * hid != 1) { printf("initial bracket is not a pair of Farey neighbours\n"); return 2; }
        int tests = 0;
        while (1) {
            i64 mn = lon + hin, md = lod + hid;
            if (md > N) break;
            if (kud > 0 && mn * kud > kun * md) {
                hin = mn; hid = md;
                printf("  mediant %lld/%lld: above the known upper bound\n", (long long)mn, (long long)md);
                continue;
            }
            i64 capX = capmul * (mn * (k + 4) + 4 * md + 64), capN = capmul * (6 * md + 64);
            if (capX > 60000) capX = 60000;
            if (capN > 60000) capN = 60000;
            ls_start(&SX, mn, md, capX); ls_start(&SN, mn, md, capN);
            i64 budget = 4 * (NR + 1);
            int verdict = 0; const char *why = "";
            while (!verdict) {
                if (!SX.done) ls_run(&SX, budget);
                if (SX.done) {
                    if (ls_top(&SX) == NR + 1) { verdict = -1; why = "MAX energy lost everywhere"; break; }
                    ls_run(&SN, (i64)1 << 62);
                    if (ls_top(&SN) == 0) { verdict = 2; why = "both certificates"; } else { verdict = 1; why = "MAX won, MIN lost"; }
                    break;
                }
                if (!SN.done) ls_run(&SN, budget);
                if (SN.done) {
                    if (ls_top(&SN) > 0) { verdict = 1; why = "MIN energy lost somewhere"; break; }
                    ls_run(&SX, (i64)1 << 62);
                    if (ls_top(&SX) < NR + 1) { verdict = 2; why = "both certificates"; } else { verdict = -1; why = "MIN won, MAX lost"; }
                    break;
                }
                budget *= 2;
            }
            tests++;
            printf("  mediant %lld/%lld: %s (%s; work MAX %lld (%d sweeps), MIN %lld (%d sweeps); %.1f s)\n",
                   (long long)mn, (long long)md, verdict == 2 ? "VALUE" : (verdict > 0 ? "rho* >" : "rho* <"), why,
                   (long long)SX.work, SX.sweeps, (long long)SN.work, SN.sweeps, now() - t0);
            fflush(stdout);
            if (verdict == 2) {
                printf("VALUE %lld %lld\n", (long long)mn, (long long)md);
                printf("STATS k=%d q=%lld tests=%d time=%.2f capmul=%lld\n", k, (long long)q, tests, now() - t0, (long long)capmul);
                if (prefix) lean_write(prefix, &SX, &SN);
                return 0;
            }
            if (verdict > 0) { lon = mn; lod = md; } else { hin = mn; hid = md; }
        }
        printf("  search failed with cap multiplier %lld; retrying\n", (long long)capmul);
    }
    printf("FAILED\n");
    return 1;
}


/* LEAN SEARCH 2: at each mediant m the two tests are dovetailed until one of them finishes; its
   one-sided verdict is taken (MAX finite somewhere: rho* >= m; MAX TOP everywhere: rho* < m; MIN finite
   everywhere: rho* <= m; MIN reaches TOP somewhere (stopped at the first TOP): rho* > m, a verdict that a
   too-small cap could make wrong, which can only make the search fail, never certify a wrong value).
   After a one-sided "rho* >= m" or "rho* <= m", the other test gets a confirmation budget of CONF times
   the work of the finished one; if it finishes with a win, m is the value (both certificates).  If the
   endpoint that is not known to differ from the value (an initial endpoint, or one set by a non-strict
   verdict) stays fixed while the other endpoint moves towards it three times in a row, that endpoint is
   re-tested to completion (it is then usually the value).                                            */
static int lsearch2_mode(i64 q, int k, i64 lon0, i64 lod0, i64 hin0, i64 hid0, i64 kun, i64 kud, const char *prefix, i64 conf) {
    init_level(k, q);
    double t0 = now();
    LState SX, SN; ls_alloc(&SX, 1); ls_alloc(&SN, 0);
    for (int attempt = 0; attempt < 3; attempt++, conf *= 4) {
        i64 lon = lon0, lod = lod0, hin = hin0, hid = hid0;
        if (hin * lod - lon * hid != 1) { printf("initial bracket is not a pair of Farey neighbours\n"); return 2; }
        int tests = 0;
        int lo_open = 1, hi_open = 1;       /* endpoint not yet decided (initial endpoints may be hints) */
        int run_hi_moves = 0, run_lo_moves = 0;
        while (1) {
            i64 mn, md;
            int forced = 0;                  /* 1: re-test lo to completion, 2: re-test hi */
            if (lo_open && run_hi_moves >= 3) { mn = lon; md = lod; forced = 1; }
            else if (hi_open && run_lo_moves >= 3) { mn = hin; md = hid; forced = 2; }
            else { mn = lon + hin; md = lod + hid; }
            if (md > N || md > 100000) break;
            if (!forced && kud > 0 && mn * kud > kun * md) {
                hin = mn; hid = md; hi_open = 0; run_lo_moves = 0; run_hi_moves++;
                printf("  mediant %lld/%lld: above the known upper bound\n", (long long)mn, (long long)md);
                continue;
            }
            i64 capX = mn * (k + 4) + 4 * md + 64, capN = 6 * md + 64;
            if (capX > 60000) capX = 60000;
            if (capN > 60000) capN = 60000;
            ls_start(&SX, mn, md, capX); ls_start(&SN, mn, md, capN); SN.stopfirst = 1;
            int verdict = 0; const char *why = "";
            if (forced) {
                ls_run(&SX, (i64)1 << 62);
                int xw = ls_top(&SX) < NR + 1;
                int nw = 0;
                if (xw) { ls_run(&SN, (i64)1 << 62); nw = !SN.lostfirst && ls_top(&SN) == 0; }
                else nw = 0;
                if (xw && nw) { verdict = 2; why = "endpoint re-tested: both certificates"; }
                else if (forced == 1) { verdict = 3; why = xw ? "endpoint re-tested: MIN lost (strict)" : "endpoint re-tested: MAX lost?!"; }
                else { verdict = -1; why = xw ? "endpoint re-tested: MIN lost?!" : "endpoint re-tested: MAX lost (strict)"; }
            } else {
                i64 budget = 4 * (NR + 1);
                while (!verdict) {
                    if (!SX.done) ls_run(&SX, budget);
                    if (SX.done) {
                        if (ls_top(&SX) == NR + 1) { verdict = -1; why = "MAX lost everywhere"; break; }
                        if (!SN.done) ls_run(&SN, conf * SX.work + 4 * (NR + 1));
                        if (SN.done && !SN.lostfirst && ls_top(&SN) == 0) { verdict = 2; why = "both certificates"; }
                        else if (SN.done) { verdict = 3; why = "MAX won, MIN lost"; }
                        else { verdict = 1; why = "MAX won, MIN confirmation budget exhausted"; }
                        break;
                    }
                    if (!SN.done) ls_run(&SN, budget);
                    if (SN.done) {
                        if (SN.lostfirst || ls_top(&SN) > 0) { verdict = 3; why = "MIN reached TOP"; break; }
                        if (!SX.done) ls_run(&SX, conf * SN.work + 4 * (NR + 1));
                        if (SX.done && ls_top(&SX) < NR + 1) { verdict = 2; why = "both certificates"; }
                        else if (SX.done) { verdict = -1; why = "MIN won, MAX lost"; }
                        else { verdict = -2; why = "MIN won, MAX confirmation budget exhausted"; }
                        break;
                    }
                    budget *= 2;
                }
            }
            tests++;
            const char *vs = verdict == 2 ? "VALUE" : verdict == 3 ? "rho* >" : verdict == 1 ? "rho* >=" : verdict == -1 ? "rho* <" : "rho* <=";
            printf("  %s %lld/%lld: %s (%s; work MAX %lld (%d sweeps), MIN %lld (%d sweeps); %.1f s)\n",
                   forced ? "endpoint" : "mediant", (long long)mn, (long long)md, vs, why,
                   (long long)SX.work, SX.sweeps, (long long)SN.work, SN.sweeps, now() - t0);
            fflush(stdout);
            if (verdict == 2) {
                printf("VALUE %lld %lld\n", (long long)mn, (long long)md);
                printf("STATS k=%d q=%lld tests=%d time=%.2f conf=%lld\n", k, (long long)q, tests, now() - t0, (long long)conf);
                if (prefix) lean_write(prefix, &SX, &SN);
                return 0;
            }
            if (forced == 1) { lo_open = 0; run_hi_moves = 0; if (verdict < 0) break; continue; }
            if (forced == 2) { hi_open = 0; run_lo_moves = 0; if (verdict > 0) break; continue; }
            if (verdict > 0) {
                lon = mn; lod = md; lo_open = (verdict == 1); run_lo_moves++; run_hi_moves = 0;
            } else {
                hin = mn; hid = md; hi_open = (verdict == -2); run_hi_moves++; run_lo_moves = 0;
            }
        }
        printf("  search failed with confirmation factor %lld; retrying\n", (long long)conf);
    }
    printf("FAILED\n");
    return 1;
}

/* LOWER-ONLY mode: the MAX energy least fixed point at one threshold F, and the compact lower certificate
   (proves rho*(q,k) >= F if the checker accepts it).  One state: about 2.1 bytes per representative. */
static int lower_mode(i64 q, int k, i64 fn, i64 fd, const char *prefix) {
    init_level(k, q);
    double t0 = now();
    LState SX; ls_alloc(&SX, 1);
    i64 capX = fn * (k + 4) + 4 * fd + 64;
    if (capX > 60000) capX = 60000;
    ls_start(&SX, fn, fd, capX);
    ls_run(&SX, (i64)1 << 62);
    i64 top = ls_top(&SX);
    printf("  LOWER F=%lld/%lld cap %lld: work %lld (%d sweeps), TOP reps %lld of %lld, %.1f s\n", (long long)fn,
           (long long)fd, (long long)capX, (long long)SX.work, SX.sweeps, (long long)top, (long long)(NR + 1), now() - t0);
    if (top == NR + 1) { printf("LOWER-FAILED\n"); return 1; }
    /* write only the lower half of the compact format */
    char fnm[1024]; FILE *f;
    const i64 CH = 1 << 20;
    uint8_t *b8 = (uint8_t *)malloc(CH);
    u16 *b16 = (u16 *)malloc(2 * CH);
    LV = SX.v;
    i64 vmax = 0;
    for (i64 r = 0; r <= NR; r++) if (LV[r] != TOP16 && LV[r] > vmax) vmax = LV[r];
    int wide = vmax > 254;
    snprintf(fnm, sizeof fnm, "%s.lo_taub", prefix); f = fopen(fnm, "wb");
    for (i64 lo = 0; lo < H; lo += 8 * CH) {
        i64 n = H - lo < 8 * CH ? H - lo : 8 * CH;
        memset(b8, 0, (n + 7) / 8);
        for (i64 i = 0; i < n; i++) {
            i64 P = lo + i;
            if (lval(P) == TOP16) continue;
            i64 Qo[2]; i64 bv[2];
            for (int b = 0; b < 2; b++) {
                int no = options(P, b, Qo); i64 v = -1;
                for (int j = 0; j < no; j++) { u16 x = lval(Qo[j]); if ((i64)x > v) v = x; }
                bv[b] = v;
            }
            if (bv[1] < bv[0]) b8[i >> 3] |= (uint8_t)(1u << (i & 7));
        }
        fwrite(b8, 1, (n + 7) / 8, f);
    }
    fclose(f);
    snprintf(fnm, sizeof fnm, "%s.%s", prefix, wide ? "lo_g16" : "lo_g8"); f = fopen(fnm, "wb");
    for (i64 lo = 0; lo < H; lo += CH) {
        i64 n = H - lo < CH ? H - lo : CH;
        if (wide) { for (i64 i = 0; i < n; i++) b16[i] = lval(lo + i); fwrite(b16, 2, n, f); }
        else { for (i64 i = 0; i < n; i++) { u16 v = lval(lo + i); b8[i] = v == TOP16 ? 255 : (uint8_t)v; } fwrite(b8, 1, n, f); }
    }
    fclose(f);
    printf("LOWER-WRITTEN %lld %lld time=%.2f\n", (long long)fn, (long long)fd, now() - t0);
    return 0;
}

/* UPPER-ONLY mode: the MIN energy least fixed point at one threshold F, and the compact upper certificate
   (proves rho*(q,k) <= F if the checker accepts it). */
static int upper_mode(i64 q, int k, i64 fn, i64 fd, const char *prefix) {
    init_level(k, q);
    double t0 = now();
    LState SN; ls_alloc(&SN, 0);
    i64 capN = 6 * fd + 64;
    if (capN > 60000) capN = 60000;
    ls_start(&SN, fn, fd, capN);
    SN.stopfirst = 1;
    ls_run(&SN, (i64)1 << 62);
    i64 top = ls_top(&SN);
    printf("  UPPER F=%lld/%lld cap %lld: work %lld (%d sweeps), TOP reps %lld, %.1f s\n", (long long)fn, (long long)fd,
           (long long)capN, (long long)SN.work, SN.sweeps, (long long)top, now() - t0);
    if (top > 0 || SN.lostfirst) { printf("UPPER-FAILED\n"); return 1; }
    char fnm[1024]; FILE *f;
    const i64 CH = 1 << 20;
    uint8_t *b8 = (uint8_t *)malloc(CH);
    u16 *b16 = (u16 *)malloc(2 * CH);
    LV = SN.v;
    i64 vmax = 0;
    for (i64 r = 0; r <= NR; r++) if (LV[r] > vmax) vmax = LV[r];
    int wide = vmax > 254;
    snprintf(fnm, sizeof fnm, "%s.up_sigb", prefix); f = fopen(fnm, "wb");
    for (i64 lo = 0; lo < H; lo += 8 * CH) {
        i64 n = H - lo < 8 * CH ? H - lo : 8 * CH;
        memset(b8, 0, (n + 7) / 8);
        for (i64 i = 0; i < n; i++) {
            i64 x = 2 * (lo + i) + 1, t = QQ * x;
            i64 Qp = ((t + 1) >> 1) & (H - 1), Qm = ((t - 1) >> 1) & (H - 1);
            if (lval(Qm) < lval(Qp)) b8[i >> 3] |= (uint8_t)(1u << (i & 7));
        }
        fwrite(b8, 1, (n + 7) / 8, f);
    }
    fclose(f);
    snprintf(fnm, sizeof fnm, "%s.%s", prefix, wide ? "up_psi16" : "up_psi8"); f = fopen(fnm, "wb");
    for (i64 lo = 0; lo < H; lo += CH) {
        i64 n = H - lo < CH ? H - lo : CH;
        if (wide) { for (i64 i = 0; i < n; i++) b16[i] = lval(lo + i); fwrite(b16, 2, n, f); }
        else { for (i64 i = 0; i < n; i++) b8[i] = (uint8_t)lval(lo + i); fwrite(b8, 1, n, f); }
    }
    fclose(f);
    printf("UPPER-WRITTEN %lld %lld time=%.2f\n", (long long)fn, (long long)fd, now() - t0);
    return 0;
}

/* exploratory: value of a frozen Max strategy tau (one-player for Min), on all pairs (no quotient):
   the largest F (Stern-Brocot, denominators <= 4096) at which the MAX energy fixed point with tau
   frozen is finite somewhere, i.e. max over starts of the least cycle density Min can reach */
static int tauval_mode(i64 q, int k, const char *taufile) {
    USEQ = 0;
    init_level(k, q);
    TAUFIX = (uint8_t *)malloc(H);
    FILE *f = fopen(taufile, "rb");
    if (!f || fread(TAUFIX, 1, H, f) != (size_t)H) { fprintf(stderr, "cannot read tau\n"); return 2; }
    fclose(f);
    State SX; st_alloc(&SX, 1);
    i64 lon = 0, lod = 1, hin = 1, hid = 1;
    while (1) {
        i64 mn = lon + hin, md = lod + hid;
        if (md > N || md > 4096) break;
        st_start(&SX, mn, md, 64 * md + 64, 0, 1);
        st_run(&SX, (i64)1 << 62);
        if (count_top(&SX) < NR + 1) { lon = mn; lod = md; } else { hin = mn; hid = md; }
    }
    printf("TAUVAL %lld %lld %lld %lld\n", (long long)lon, (long long)lod, (long long)hin, (long long)hid);
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr, "usage: %s value q k lon lod hin hid kun kud [prefix] | lean q k fn fd prefix | tauval q k taufile\n", argv[0]);
        return 2;
    }
    i64 q = atoll(argv[2]);
    int k = atoi(argv[3]);
    if (k < 3 || k > 30 || q < 1 || !(q & 1)) { fprintf(stderr, "bad q/k\n"); return 2; }
    if (!strcmp(argv[1], "value")) {
        if (argc < 10) { fprintf(stderr, "value q k lon lod hin hid kun kud [prefix]\n"); return 2; }
        return value_mode(q, k, atoll(argv[4]), atoll(argv[5]), atoll(argv[6]), atoll(argv[7]),
                          atoll(argv[8]), atoll(argv[9]), argc > 10 ? argv[10] : 0);
    }
    if (!strcmp(argv[1], "lean")) {
        if (argc < 7) { fprintf(stderr, "lean q k fn fd prefix\n"); return 2; }
        return lean_mode(q, k, atoll(argv[4]), atoll(argv[5]), argv[6]);
    }
    if (!strcmp(argv[1], "lsearch")) {
        if (argc < 10) { fprintf(stderr, "lsearch q k lon lod hin hid kun kud [prefix]\n"); return 2; }
        return lsearch_mode(q, k, atoll(argv[4]), atoll(argv[5]), atoll(argv[6]), atoll(argv[7]),
                            atoll(argv[8]), atoll(argv[9]), argc > 10 ? argv[10] : 0);
    }
    if (!strcmp(argv[1], "lsearch2")) {
        if (argc < 10) { fprintf(stderr, "lsearch2 q k lon lod hin hid kun kud [prefix [conf]]\n"); return 2; }
        return lsearch2_mode(q, k, atoll(argv[4]), atoll(argv[5]), atoll(argv[6]), atoll(argv[7]),
                             atoll(argv[8]), atoll(argv[9]), argc > 10 ? argv[10] : 0, argc > 11 ? atoll(argv[11]) : 8);
    }
    if (!strcmp(argv[1], "lower")) {
        if (argc < 7) { fprintf(stderr, "lower q k fn fd prefix\n"); return 2; }
        return lower_mode(q, k, atoll(argv[4]), atoll(argv[5]), argv[6]);
    }
    if (!strcmp(argv[1], "upper")) {
        if (argc < 7) { fprintf(stderr, "upper q k fn fd prefix\n"); return 2; }
        return upper_mode(q, k, atoll(argv[4]), atoll(argv[5]), argv[6]);
    }
    if (!strcmp(argv[1], "tauval")) return tauval_mode(q, k, argv[4]);
    fprintf(stderr, "unknown mode\n");
    return 2;
}
