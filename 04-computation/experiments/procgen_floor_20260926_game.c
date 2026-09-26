/*
 * procgen_floor_20260926_game.c -- the parity arena of the q n +- 1 strategy cube as a two-player
 * mean-payoff game; floor lane of session collatz-procgen-20260922 (2026-09-26).
 *
 * Compiled as a shared library and driven from Python by ctypes (procgen_floor_20260926_lib.py).
 * Nothing computed here is used as a proof step directly: every certificate it returns (a potential
 * for Min, a lift strategy plus a potential for Max) is re-checked edge by edge in Python with exact
 * integer arithmetic.
 *
 * Arena (level k, odd q; N = 2^k, H = 2^(k-1)).  Nodes s in Z/N; pairs P in Z/H (pair P = {P, P+H}).
 *   even s : one option, the pair s/2;
 *   odd s  : two options, pair (q s + 1)/2 mod H  (sign +)  and  (q s - 1)/2 mod H  (sign -).
 * Min (the strategy maker) picks the option at every odd node, once and for all (a sign strategy);
 * Max (the adversary) picks the lift P or P+H at every pair.  Payoff: the odd density of the play.
 *
 * Threshold F = fn/fd; node weight e(s) = fd - fn (s odd), -fn (s even); a cycle has e-weight > 0,
 * = 0, < 0 iff its odd density is > F, = F, < F.
 *
 * MODE_MIN (Min is the energy player): least fixed point of
 *     f(s) = max(0, e(s) + min_{P option of s} max(f(P), f(P+H)))      (values > cap become TOP)
 *   f(s) < TOP everywhere  ==>  the extracted sign strategy (argmin option) has a potential f with
 *   f(t) + e(s) <= f(s) on every edge s -> t of G_sigma, so every cycle of G_sigma has density <= F.
 * MODE_MAX (Max is the energy player): least fixed point of
 *     f(s) = max(0, -e(s) + max_{P option of s} min(f(P), f(P+H)))
 *   W = {s : f(s) < TOP} nonempty  ==>  with tau(P) = argmin(f(P), f(P+H)), every option P of every
 *   s in W has tau(P) in W and f(tau(P)) <= f(s) + e(s); so for EVERY sign strategy the (sigma,tau)
 *   play from W stays in W and closes a cycle of G_sigma with density >= F.
 * Caps only shrink the non-TOP region (they never create false certificates).
 */
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

typedef long long ll;
#define TOPV (LLONG_MAX / 4)

static int K = 0, N = 0, H = 0;
static ll Q = 0, QINV = 0;
static ll *FV = 0;
static int *QU = 0;
static uint8_t *INQ = 0;
static uint8_t *ALLOW = 0;   /* odd s: bit0 = option + allowed, bit1 = option - allowed */
static uint8_t *TAUF = 0;    /* optional: Max frozen to the lift P + TAUF[P] H (TAUF[P] in {0,1}); 2 = free */

static ll inv_pow2(ll q, int k) {
    unsigned long long x = (unsigned long long)q;
    for (int i = 0; i < 7; i++) x = x * (2ULL - (unsigned long long)q * x);
    return (ll)(x & ((1ULL << k) - 1ULL));
}

int g_init(int k, ll q) {
    free(FV); free(QU); free(INQ); free(ALLOW); free(TAUF);
    K = k; N = 1 << k; H = N >> 1; Q = q; QINV = inv_pow2(q, k);
    FV = (ll *)calloc(N, sizeof(ll));
    QU = (int *)malloc(sizeof(int) * (N + 1));
    INQ = (uint8_t *)calloc(N, 1);
    ALLOW = (uint8_t *)malloc(N);
    memset(ALLOW, 3, N);
    TAUF = (uint8_t *)malloc(H);
    memset(TAUF, 2, H);
    return (int)QINV;
}

void g_set_allow(const uint8_t *a) { memcpy(ALLOW, a, N); }
void g_set_tau(const uint8_t *t) { memcpy(TAUF, t, H); }

static inline int pplus(int s) { return (int)(((Q * (ll)s + 1) >> 1) & (ll)(H - 1)); }
static inline int pminus(int s) { return (int)(((Q * (ll)s - 1) >> 1) & (ll)(H - 1)); }
static inline ll mx(ll a, ll b) { return a > b ? a : b; }
static inline ll mn(ll a, ll b) { return a < b ? a : b; }
/* value of pair P for the MIN-energy player (Max maximizes) and the MAX-energy player (Max minimizes) */
static inline ll pv_min(int P) { return TAUF[P] == 2 ? mx(FV[P], FV[P + H]) : FV[P + (TAUF[P] ? H : 0)]; }
static inline ll pv_max(int P) { return TAUF[P] == 2 ? mn(FV[P], FV[P + H]) : FV[P + (TAUF[P] ? H : 0)]; }

static inline ll eval_node(int s, int mode, ll fn, ll fd, ll cap) {
    ll e = (s & 1) ? (fd - fn) : -fn;
    ll nxt;
    if (mode == 0) {                     /* MIN energy: min over options of max over lifts */
        if (!(s & 1)) {
            int P = (s >> 1) & (H - 1);
            nxt = pv_min(P);
        } else {
            nxt = TOPV;
            if (ALLOW[s] & 1) { int P = pplus(s); nxt = mn(nxt, pv_min(P)); }
            if (ALLOW[s] & 2) { int P = pminus(s); nxt = mn(nxt, pv_min(P)); }
        }
        if (nxt >= TOPV) return TOPV;
        ll v = nxt + e;
        if (v < 0) v = 0;
        return v > cap ? TOPV : v;
    } else {                             /* MAX energy: max over options of min over lifts */
        if (!(s & 1)) {
            int P = (s >> 1) & (H - 1);
            nxt = pv_max(P);
        } else {
            nxt = -1;
            if (ALLOW[s] & 1) { int P = pplus(s); nxt = mx(nxt, pv_max(P)); }
            if (ALLOW[s] & 2) { int P = pminus(s); nxt = mx(nxt, pv_max(P)); }
        }
        if (nxt >= TOPV) return TOPV;
        ll v = nxt - e;
        if (v < 0) v = 0;
        return v > cap ? TOPV : v;
    }
}

/* least fixed point from f = 0 (FIFO chaotic iteration); returns the number of TOP nodes, or -1 if
   the work limit (node evaluations) was exceeded */
ll g_solve(int mode, ll fn, ll fd, ll cap, ll maxwork) {
    for (int s = 0; s < N; s++) { FV[s] = 0; INQ[s] = 1; QU[s] = s; }
    int qh = 0, qt = N % (N + 1), qn = N;
    ll work = 0;
    while (qn > 0) {
        int s = QU[qh]; qh = (qh + 1) % (N + 1); qn--; INQ[s] = 0;
        if (++work > maxwork) return -1;
        if (FV[s] >= TOPV) continue;
        ll v = eval_node(s, mode, fn, fd, cap);
        if (v > FV[s]) {
            FV[s] = v;
            int P = s & (H - 1);
            int u[3];
            u[0] = (int)((2 * (ll)P) & (ll)(N - 1));
            u[1] = (int)((QINV * ((2 * (ll)P - 1) & (ll)(N - 1))) & (ll)(N - 1));
            u[2] = (int)((QINV * ((2 * (ll)P + 1) & (ll)(N - 1))) & (ll)(N - 1));
            for (int i = 0; i < 3; i++) {
                int w = u[i];
                if (!INQ[w] && FV[w] < TOPV) { INQ[w] = 1; QU[qt] = w; qt = (qt + 1) % (N + 1); qn++; }
            }
        }
    }
    ll ntop = 0;
    for (int s = 0; s < N; s++) if (FV[s] >= TOPV) ntop++;
    return ntop;
}

void g_get_f(ll *out) {
    for (int s = 0; s < N; s++) out[s] = FV[s] >= TOPV ? -1 : FV[s];
}

/* Min's sign strategy from the current f: flip[s] = 1 iff the '-' option is chosen (ties: '+') */
void g_extract_sigma(uint8_t *flip) {
    for (int s = 0; s < N; s++) {
        flip[s] = 0;
        if (!(s & 1)) continue;
        ll vp = TOPV + 1, vm = TOPV + 1;
        if (ALLOW[s] & 1) { int P = pplus(s); vp = pv_min(P); }
        if (ALLOW[s] & 2) { int P = pminus(s); vm = pv_min(P); }
        flip[s] = (vm < vp) ? 1 : 0;
    }
}

/* Max's lift strategy from the current f: tau[P] = 1 iff the lift P+H is chosen (ties: P) */
void g_extract_tau(uint8_t *tau) {
    for (int P = 0; P < H; P++) tau[P] = TAUF[P] != 2 ? TAUF[P] : ((FV[P + H] < FV[P]) ? 1 : 0);
}

/* ------------------------------------------------------------------ one-player: min cycle mean
   of the Min graph G^tau restricted to the set reachable from a start node (Min controls the
   options, Max is frozen to tau); used only as a cross-check of the Max certificates.
   Karp's algorithm on the reachable set, exact rational result (num/den, unreduced). */
void g_min_mean_tau(const uint8_t *tau, int start, ll *num, ll *den) {
    int *idx = (int *)malloc(sizeof(int) * N);
    int *lst = (int *)malloc(sizeof(int) * N);
    for (int s = 0; s < N; s++) idx[s] = -1;
    int n = 0;
    lst[n] = start; idx[start] = n; n++;
    for (int i = 0; i < n; i++) {
        int s = lst[i];
        int opts[2], no = 0;
        if (!(s & 1)) opts[no++] = (s >> 1) & (H - 1);
        else {
            if (ALLOW[s] & 1) opts[no++] = pplus(s);
            if (ALLOW[s] & 2) opts[no++] = pminus(s);
        }
        for (int j = 0; j < no; j++) {
            int t = opts[j] + (tau[opts[j]] ? H : 0);
            if (idx[t] < 0) { idx[t] = n; lst[n++] = t; }
        }
    }
    /* Karp: D_j(v) = min weight of a j-edge walk from start to v (weight = odd indicator of the source) */
    const ll INF = (ll)1 << 60;
    ll *Dall = (ll *)malloc(sizeof(ll) * (size_t)(n + 1) * (size_t)n);
    if (!Dall) { *num = -1; *den = 0; free(idx); free(lst); return; }
    for (int v = 0; v < n; v++) Dall[v] = INF;
    Dall[0] = 0;
    for (int j = 1; j <= n; j++) {
        ll *Dp = Dall + (size_t)(j - 1) * n, *Dc = Dall + (size_t)j * n;
        for (int v = 0; v < n; v++) Dc[v] = INF;
        for (int a = 0; a < n; a++) {
            if (Dp[a] >= INF) continue;
            int s = lst[a];
            ll val = Dp[a] + (s & 1);
            int opts[2], no = 0;
            if (!(s & 1)) opts[no++] = (s >> 1) & (H - 1);
            else {
                if (ALLOW[s] & 1) opts[no++] = pplus(s);
                if (ALLOW[s] & 2) opts[no++] = pminus(s);
            }
            for (int jj = 0; jj < no; jj++) {
                int t = opts[jj] + (tau[opts[jj]] ? H : 0);
                int b = idx[t];
                if (val < Dc[b]) Dc[b] = val;
            }
        }
    }
    ll bn = -1, bd = 1;
    ll *Dn = Dall + (size_t)n * n;
    for (int v = 0; v < n; v++) {
        if (Dn[v] >= INF) continue;
        ll wn = -1, wd = 1;   /* max over j of (Dn - Dj)/(n - j) */
        for (int j = 0; j < n; j++) {
            ll dj = Dall[(size_t)j * n + v];
            if (dj >= INF) continue;
            ll nu = Dn[v] - dj, de = n - j;
            if (wn < 0 || nu * wd > wn * de) { wn = nu; wd = de; }
        }
        if (wn < 0) continue;
        if (bn < 0 || wn * bd < bn * wd) { bn = wn; bd = wd; }
    }
    *num = bn; *den = bd;
    free(Dall); free(idx); free(lst);
}

/* ------------------------------------------------------------------ exhaustive cross-check (small k)
   For every flip set R (all 2^H of them), the maximum odd density of a cycle of G_sigma by Karp's
   algorithm (all nodes as sources via a virtual source); returns the minimum over R as num/den and
   the number of flip sets attaining it.  Independent of the game solver.  k <= 5 only. */
static void karp_max(const int *tg, ll *num, ll *den) {
    const int NEG = -1000000;
    int D[33][32];
    for (int v = 0; v < N; v++) D[0][v] = 0;           /* virtual source to every node */
    for (int j = 1; j <= N; j++) {
        for (int v = 0; v < N; v++) D[j][v] = NEG;
        for (int u = 0; u < N; u++) {
            if (D[j - 1][u] == NEG) continue;
            int val = D[j - 1][u] + (u & 1);
            int a = tg[u], b = tg[u] + H;
            if (val > D[j][a]) D[j][a] = val;
            if (val > D[j][b]) D[j][b] = val;
        }
    }
    ll bn = -1, bd = 1;
    for (int v = 0; v < N; v++) {
        if (D[N][v] == NEG) continue;
        ll wn = -1, wd = 1;                              /* min over j of (D_N - D_j)/(N - j) */
        for (int j = 0; j < N; j++) {
            if (D[j][v] == NEG) continue;
            ll nu = D[N][v] - D[j][v], de = N - j;
            if (wn < 0 || nu * wd < wn * de) { wn = nu; wd = de; }
        }
        if (wn < 0) continue;
        if (bn < 0 || wn * bd > bn * wd) { bn = wn; bd = wd; }
    }
    *num = bn; *den = bd;
}

int g_exhaustive_minmax(ll *num, ll *den, ll *count) {
    if (K > 5) return -1;
    int odd[16], no = 0;
    for (int s = 1; s < N; s += 2) odd[no++] = s;
    int tg[32];
    ll bn = -1, bd = 1, cnt = 0;
    for (long m = 0; m < (1L << no); m++) {
        for (int s = 0; s < N; s++) tg[s] = (s >> 1) & (H - 1);
        for (int i = 0; i < no; i++) {
            int s = odd[i];
            tg[s] = ((m >> i) & 1) ? pminus(s) : pplus(s);
        }
        ll a, b;
        karp_max(tg, &a, &b);
        if (bn < 0 || a * bd < bn * b) { bn = a; bd = b; cnt = 1; }
        else if (a * bd == bn * b) cnt++;
    }
    *num = bn; *den = bd; *count = cnt;
    return 0;
}
