/*
 * procgen_drift_20260926_engine.c -- exact routines for the strategy cube of the maps q n +- 1 (q odd),
 * drift lane of session collatz-procgen-20260922 (2026-09-26).  Compiled as a shared library and driven
 * from Python by ctypes (procgen_drift_20260926_lib.py).  Every certificate produced here is re-checked
 * edge by edge in Python before it is used as a proof step.
 *
 * A level-k sign strategy sigma : odd residues mod 2^k -> {+1,-1};
 *   T(n) = n/2 (n even),  (q n + sigma(n mod 2^k))/2 (n odd).
 * Parity graph G_sigma: nodes Z/2^k; s -> the two lifts mod 2^k of T(s) mod 2^(k-1).
 * flip[s] = 1 means sigma(s) = -1 (s odd).  A cycle with a odd nodes and length p is expanding iff q^a > 2^p.
 *
 * Threshold weights (Lemma P): for F = fn/fd, w(odd) = fd - fn, w(even) = -fn; a cycle has positive weight iff
 * its odd density exceeds F.  The least potential psi(s) = sup over walks from s of the partial weight sums
 * satisfies psi(t) <= psi(s) - w(s) on every edge, and is finite iff no cycle has density > F.
 *
 * State (one instance per process): K, N = 2^K, H = N/2, Q, flip[], psi[], best[] (argmax successor pointer).
 * Positive cycles are detected by cycles of the pointer graph s -> best[s]; any such cycle has positive weight
 * (Bellman-Ford predecessor lemma; valid for any initial psi because psi only increases and best[s] is set at
 * every increase of psi[s]).
 */
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef long long ll;

static int K = 0, N = 0, H = 0;
static ll Q = 0, QINV = 0;
static ll WO = 0, WE = 0;
static uint8_t *FL = 0;
static ll *PSI = 0;
static int *BEST = 0;
/* work arrays */
static int *QU = 0;
static uint8_t *INQ = 0;
static int *STAMP = 0;
static int STAMPV = 0;
static int *TOUCH = 0;
static int NTOUCH = 0;
/* undo log */
static int *LOG_S = 0;
static ll *LOG_PSI = 0;
static int *LOG_BEST = 0;
static int NLOG = 0;
static int *LOGSTAMP = 0;
static int LOGV = 0;
/* pointer-cycle scratch */
static ll *PMARK = 0;
static ll WALK = 0;
static int *CYC = 0;
static int CLEN = 0;
static uint8_t *ALIVE = 0;   /* harvesting only: dead nodes are removed from the graph */
static int TOG = -1;         /* node toggled by eng_toggle (-1 otherwise) */


static ll mod_inverse_pow2(ll q, int k) {
    /* Newton iteration for q^{-1} mod 2^k (q odd) */
    unsigned long long x = (unsigned long long)q;  /* correct mod 8 */
    for (int i = 0; i < 6; i++) x = x * (2ULL - (unsigned long long)q * x);
    unsigned long long m = (k >= 64) ? ~0ULL : ((1ULL << k) - 1ULL);
    return (ll)(x & m);
}

static void freeall(void) {
    free(FL); free(PSI); free(BEST); free(QU); free(INQ); free(STAMP); free(TOUCH);
    free(LOG_S); free(LOG_PSI); free(LOG_BEST); free(LOGSTAMP); free(PMARK); free(CYC);
    FL = 0; PSI = 0; BEST = 0; QU = 0; INQ = 0; STAMP = 0; TOUCH = 0; LOG_S = 0; LOG_PSI = 0; LOG_BEST = 0;
    LOGSTAMP = 0; PMARK = 0; CYC = 0;
    free(ALIVE); ALIVE = 0;
}

/* initialise the instance; flip is an array of N bytes (only odd entries used) */
int eng_init(int k, ll q, const uint8_t *flip, ll fn, ll fd) {
    freeall();
    K = k; N = 1 << k; H = N >> 1; Q = q; QINV = mod_inverse_pow2(q, k);
    WO = fd - fn; WE = -fn;
    FL = (uint8_t *)calloc(N, 1);
    for (int s = 1; s < N; s += 2) FL[s] = flip[s] ? 1 : 0;
    PSI = (ll *)calloc(N, sizeof(ll));
    BEST = (int *)malloc(sizeof(int) * N);
    QU = (int *)malloc(sizeof(int) * (N + 1));
    INQ = (uint8_t *)calloc(N, 1);
    STAMP = (int *)calloc(N, sizeof(int));
    TOUCH = (int *)malloc(sizeof(int) * N);
    LOG_S = (int *)malloc(sizeof(int) * N);
    LOG_PSI = (ll *)malloc(sizeof(ll) * N);
    LOG_BEST = (int *)malloc(sizeof(int) * N);
    LOGSTAMP = (int *)calloc(N, sizeof(int));
    PMARK = (ll *)calloc(N, sizeof(ll));
    CYC = (int *)malloc(sizeof(int) * (N + 1));
    ALIVE = (uint8_t *)malloc(N); memset(ALIVE, 1, N);
    for (int s = 0; s < N; s++) BEST[s] = -1;
    STAMPV = 0; LOGV = 0; WALK = 0; NLOG = 0; CLEN = 0;
    return (int)QINV;
}

static inline int target(int s) {
    /* T(s) mod H */
    if (!(s & 1)) return (s >> 1) & (H - 1);
    ll v = Q * (ll)s + (FL[s] ? -1 : 1);
    return (int)((v >> 1) & (ll)(H - 1));
}

static inline ll wt(int s) { return (s & 1) ? WO : WE; }

/* predecessors of node t: fills out[] and returns count (<= 3) */
static inline int preds(int t, int *out) {
    int n = 0;
    ll tp = t & (H - 1);
    out[n++] = (int)((2 * tp) & (N - 1));                          /* even preimage */
    int sp = (int)((QINV * ((2 * tp - 1) & (ll)(N - 1))) & (ll)(N - 1));   /* q s + 1 = 2 t' mod N */
    if (!FL[sp]) out[n++] = sp;
    int sm = (int)((QINV * ((2 * tp + 1) & (ll)(N - 1))) & (ll)(N - 1));   /* q s - 1 = 2 t' mod N */
    if (FL[sm]) out[n++] = sm;
    return n;
}

static void logsave(int s) {
    if (LOGSTAMP[s] != LOGV) {
        LOGSTAMP[s] = LOGV;
        LOG_S[NLOG] = s; LOG_PSI[NLOG] = PSI[s]; LOG_BEST[NLOG] = BEST[s]; NLOG++;
    }
}

/* look for a cycle of the pointer graph s -> best[s] reachable from the touched nodes; fills CYC/CLEN.
   Walk ids are global and increasing; a node marked by an earlier walk of the same check lies on a chain
   already known to be acyclic, so the walk stops there.  Each node is visited at most once per check. */
static int pointer_cycle_from_touched(void) {
    ll base = WALK;
    for (int i = 0; i < NTOUCH; i++) {
        ll id = ++WALK;
        int v = TOUCH[i];
        while (v >= 0) {
            if (PMARK[v] == id) break;              /* closed a cycle in this walk */
            if (PMARK[v] > base) { v = -1; break; } /* seen in an earlier walk of this check */
            PMARK[v] = id;
            v = BEST[v];
        }
        if (v >= 0) {
            int u = v, L = 0;
            do { CYC[L++] = u; u = BEST[u]; } while (u != v && L <= N);
            CLEN = L;
            return L;
        }
    }
    return 0;
}

/* SPFA propagation from the current queue; returns 1 if converged, 0 if a positive cycle was found (CYC) */
static int propagate(int *qh, int *qt, int *qn, int uselog) {
    ll relax = 0, nextcheck = N / 4 + 1024;
    int pr[3];
    while (*qn > 0) {
        int t = QU[*qh]; *qh = (*qh + 1) % (N + 1); (*qn)--; INQ[t] = 0;
        int np = preds(t, pr);
        for (int i = 0; i < np; i++) {
            int s = pr[i];
            if (!ALIVE[s]) continue;
            ll cand = wt(s) + PSI[t];
            if (cand > PSI[s]) {
                if (s == TOG) {
                    /* second raise of the toggled node: every new cycle passes through TOG, and a cycle through TOG of
                       non-positive weight can never raise it again, so a positive cycle through TOG exists */
                    while (*qn > 0) { int u = QU[*qh]; *qh = (*qh + 1) % (N + 1); (*qn)--; INQ[u] = 0; }
                    CLEN = 0;
                    return 0;
                }
                if (uselog) logsave(s);
                PSI[s] = cand; BEST[s] = t; relax++;
                if (STAMP[s] != STAMPV) { STAMP[s] = STAMPV; TOUCH[NTOUCH++] = s; }
                if (!INQ[s]) { QU[*qt] = s; *qt = (*qt + 1) % (N + 1); (*qn)++; INQ[s] = 1; }
            }
        }
        if (relax >= nextcheck) {
            nextcheck = relax + (ll)NTOUCH + N / 4 + 1024;
            if (pointer_cycle_from_touched()) {
                /* drain the queue */
                while (*qn > 0) { int u = QU[*qh]; *qh = (*qh + 1) % (N + 1); (*qn)--; INQ[u] = 0; }
                return 0;
            }
        }
    }
    /* the queue is empty: every edge satisfies psi(t) <= psi(s) - w(s), so there is no positive cycle */
    return 1;
}

/* least potential from scratch (psi = 0, best = -1); returns 1 (OK) or 0 (positive cycle in CYC) */
int eng_solve(void) {
    TOG = -1;
    for (int s = 0; s < N; s++) { PSI[s] = 0; BEST[s] = -1; INQ[s] = 0; }
    STAMPV++; NTOUCH = 0;
    int qh = 0, qt = 0, qn = 0;
    for (int t = 0; t < N; t++) if (ALIVE[t]) { QU[qt] = t; qt = (qt + 1) % (N + 1); qn++; INQ[t] = 1; }
    return propagate(&qh, &qt, &qn, 0);
}

/* harvesting: remove nodes from the graph (dead nodes are ignored by eng_solve; do not mix with eng_toggle) */
void eng_kill(const int *nodes, int n) { for (int i = 0; i < n; i++) ALIVE[nodes[i]] = 0; }
void eng_revive_all(void) { memset(ALIVE, 1, N); }

/* toggle the sign at odd residue r and restore a valid potential incrementally.
   Returns 1 (toggle kept, potential valid) or 0 (the toggle creates a cycle above the threshold; the toggle and all
   potential changes are undone).  The instance must hold a valid potential (after eng_solve() == 1).
   Detection: every new cycle passes through r; the propagation declares a positive cycle when r would be raised a
   second time (TOG rule) or when the pointer graph closes a cycle (checked periodically). */
int eng_toggle(int r) {
    LOGV++; NLOG = 0; STAMPV++; NTOUCH = 0;
    logsave(r);
    FL[r] ^= 1;
    BEST[r] = -1;
    int t0 = target(r);
    int qh = 0, qt = 0, qn = 0;
    ll need = PSI[t0] > PSI[t0 + H] ? PSI[t0] : PSI[t0 + H];
    int bt = PSI[t0] >= PSI[t0 + H] ? t0 : t0 + H;
    if (wt(r) + need > PSI[r]) {
        PSI[r] = wt(r) + need; BEST[r] = bt;
        STAMP[r] = STAMPV; TOUCH[NTOUCH++] = r;
        QU[qt] = r; qt = (qt + 1) % (N + 1); qn++; INQ[r] = 1;
    }
    TOG = r;
    int ok = propagate(&qh, &qt, &qn, 1);
    TOG = -1;
    if (!ok) {
        for (int i = NLOG - 1; i >= 0; i--) { PSI[LOG_S[i]] = LOG_PSI[i]; BEST[LOG_S[i]] = LOG_BEST[i]; }
        FL[r] ^= 1;
        return 0;
    }
    return 1;
}

/* ---------------------------------------------------------------- short expanding cycles (DFS) */
/* every simple cycle of length <= P whose odd count a satisfies q^a > 2^p (tested as a/p > fn/fd with the exact
   threshold, valid for p <= fd's range: the caller passes the best lower approximation with denominator >= P),
   rooted at its minimum node (each cycle once); needflip = 1: only cycles through a flipped node.
   Cycles are written consecutively into out (node lists), lengths into lens; returns the count (<= cap). */
static int SP_P, SP_cap, SP_count, SP_needflip, SP_outpos, SP_outcap;
static int *SP_path, *SP_out, *SP_lens;
static uint8_t *SP_on;
static void sp_dfs(int root, int v, int depth, ll w) {
    if (SP_count >= SP_cap) return;
    ll w2 = w + wt(v);
    int t0 = target(v);
    int ts[2] = { t0, t0 + H };
    for (int i = 0; i < 2; i++) {
        int t = ts[i];
        if (t == root) {
            if (w2 > 0) {
                int hasflip = 0;
                for (int j = 0; j <= depth; j++) if ((SP_path[j] & 1) && FL[SP_path[j]]) hasflip = 1;
                if ((!SP_needflip || hasflip) && SP_outpos + depth + 1 <= SP_outcap) {
                    for (int j = 0; j <= depth; j++) SP_out[SP_outpos++] = SP_path[j];
                    SP_lens[SP_count++] = depth + 1;
                    if (SP_count >= SP_cap) return;
                }
            }
            continue;
        }
        if (t < root || SP_on[t] || depth + 1 >= SP_P) continue;
        if (w2 + (ll)(SP_P - depth - 1) * WO <= 0) continue;
        SP_on[t] = 1; SP_path[depth + 1] = t;
        sp_dfs(root, t, depth + 1, w2);
        SP_on[t] = 0;
    }
}
int eng_short_cycles(int P, int cap, int needflip, int *out, int outcap, int *lens) {
    SP_P = P; SP_cap = cap; SP_needflip = needflip; SP_count = 0; SP_outpos = 0; SP_outcap = outcap;
    SP_out = out; SP_lens = lens;
    SP_path = (int *)malloc(sizeof(int) * (P + 2));
    SP_on = (uint8_t *)calloc(N, 1);
    for (int root = 0; root < N && SP_count < SP_cap; root++) {
        SP_path[0] = root; SP_on[root] = 1;
        sp_dfs(root, root, 0, 0);
        SP_on[root] = 0;
    }
    free(SP_path); free(SP_on);
    return SP_count;
}

/* read-outs */
void eng_get_psi(ll *out) { memcpy(out, PSI, sizeof(ll) * N); }
void eng_get_flip(uint8_t *out) { memcpy(out, FL, N); }
int eng_get_cycle(int *out) { memcpy(out, CYC, sizeof(int) * CLEN); return CLEN; }
int eng_n(void) { return N; }

/* ---------------------------------------------------------------- exact Karp (max odd density) for small k */
/* returns numerator/denominator of the maximum cycle mean of the 0/1 odd indicator (unreduced) */
void eng_karp(ll *num, ll *den) {
    const int NEG = -1000000000;
    int *Dp = (int *)malloc(sizeof(int) * N), *Dc = (int *)malloc(sizeof(int) * N), *Dn = (int *)malloc(sizeof(int) * N);
    int *tg = (int *)malloc(sizeof(int) * N);
    for (int s = 0; s < N; s++) tg[s] = target(s);
    for (int v = 0; v < N; v++) Dp[v] = 0;
    for (int j = 1; j <= N; j++) {
        for (int v = 0; v < N; v++) Dc[v] = NEG;
        for (int u = 0; u < N; u++) {
            if (Dp[u] == NEG) continue;
            int val = Dp[u] + (u & 1);
            int a = tg[u], b = tg[u] + H;
            if (val > Dc[a]) Dc[a] = val;
            if (val > Dc[b]) Dc[b] = val;
        }
        int *t = Dp; Dp = Dc; Dc = t;
    }
    memcpy(Dn, Dp, sizeof(int) * N);
    ll *bn = (ll *)malloc(sizeof(ll) * N), *bd = (ll *)malloc(sizeof(ll) * N);
    for (int v = 0; v < N; v++) { bn[v] = 1; bd[v] = 0; }
    for (int v = 0; v < N; v++) Dp[v] = 0;
    for (int j = 0; j < N; j++) {
        if (j > 0) {
            for (int v = 0; v < N; v++) Dc[v] = NEG;
            for (int u = 0; u < N; u++) {
                if (Dp[u] == NEG) continue;
                int val = Dp[u] + (u & 1);
                int a = tg[u], b = tg[u] + H;
                if (val > Dc[a]) Dc[a] = val;
                if (val > Dc[b]) Dc[b] = val;
            }
            int *t = Dp; Dp = Dc; Dc = t;
        }
        for (int v = 0; v < N; v++) {
            if (Dn[v] == NEG || Dp[v] == NEG) continue;
            ll nu = (ll)Dn[v] - Dp[v], de = N - j;
            if (bd[v] == 0 || nu * bd[v] < bn[v] * de) { bn[v] = nu; bd[v] = de; }
        }
    }
    ll An = -1, Ad = 1;
    for (int v = 0; v < N; v++) {
        if (Dn[v] == NEG || bd[v] == 0) continue;
        if (An < 0 || bn[v] * Ad > An * bd[v]) { An = bn[v]; Ad = bd[v]; }
    }
    *num = An; *den = Ad;
    free(Dp); free(Dc); free(Dn); free(bn); free(bd); free(tg);
}
