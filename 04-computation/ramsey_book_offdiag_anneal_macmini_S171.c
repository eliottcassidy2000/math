/*
 * ramsey_anneal.c — fast two-level annealer for R(B_{n-1}, B_n) >= 4n-1
 * witnesses on N = 4n-2 = 2m vertices (m = 2n-1 odd).
 *
 * Model: vertices Z_m x {0,1};
 *   (x,0)~(y,0) iff y-x in U0 (symmetric), (x,1)~(y,1) iff y-x in U1 (symmetric),
 *   (x,0)~(y,1) iff y-x in D (arbitrary).
 * Exact codegree conditions via correlations:
 *   within level L, diff d!=0: W_L(d) = A_{UL}(d) + A_D(d)
 *     limit: d in UL ? n-2 : n+1+2*kL-N     (kL = |UL| + |D|)
 *   cross diff t: X(t) = conv(U0,D)(t) + corr(D,U1)(t)
 *     limit: t in D ? n-2 : n+1+k0+k1-N
 * Energy = sum of positive excesses (class-level).  Zero => witness; final
 * brute check done by the Python verifier on the emitted sets.
 *
 * Incremental O(m) updates per toggle; simulated annealing with reheats.
 * Usage: ./ramsey_anneal n iters seed restarts [emit_prefix]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static int m, n, N;
static int U0[512], U1[512], D[512];
static long AU0[512], AU1[512], AD[512], XC1[512], XC2[512];
/* XC1(t) = sum_d D(d) U0(t-d);  XC2(t) = sum_u U1(u) D(t+u) */

static unsigned long long rng_s;
static inline unsigned long long rnd(void) {
    rng_s ^= rng_s << 13; rng_s ^= rng_s >> 7; rng_s ^= rng_s << 17;
    return rng_s;
}
static inline double rndu(void) { return (rnd() >> 11) * (1.0 / 9007199254740992.0); }

static void recompute_all(void) {
    for (int d = 0; d < m; d++) {
        long s0 = 0, s1 = 0, sd = 0, x1 = 0, x2 = 0;
        for (int e = 0; e < m; e++) {
            int ed = e - d; if (ed < 0) ed += m;
            int eD = e + d; if (eD >= m) eD -= m;
            s0 += U0[e] & U0[ed];
            s1 += U1[e] & U1[ed];
            sd += D[e] & D[ed];
            x1 += D[e] & U0[ed];        /* U0(d-e)? careful below */
            x2 += U1[e] & D[eD];
        }
        AU0[d] = s0; AU1[d] = s1; AD[d] = sd; XC2[d] = x2;
        (void)x1;
    }
    /* XC1(t) = sum_d D(d) U0(t-d) */
    for (int t = 0; t < m; t++) {
        long s = 0;
        for (int d = 0; d < m; d++) {
            int td = t - d; if (td < 0) td += m;
            s += D[d] & U0[td];
        }
        XC1[t] = s;
    }
}

static long energy(void) {
    int su0 = 0, su1 = 0, sd = 0;
    for (int i = 0; i < m; i++) { su0 += U0[i]; su1 += U1[i]; sd += D[i]; }
    int k0 = su0 + sd, k1 = su1 + sd;
    int lim_on = n - 2;
    int off0 = n + 1 + 2 * k0 - N, off1 = n + 1 + 2 * k1 - N, offX = n + 1 + k0 + k1 - N;
    long e = 0;
    for (int d = 1; d < m; d++) {
        long w0 = AU0[d] + AD[d];
        long l0 = U0[d] ? lim_on : off0;
        if (w0 > l0) e += (w0 - l0);
        long w1 = AU1[d] + AD[d];
        long l1 = U1[d] ? lim_on : off1;
        if (w1 > l1) e += (w1 - l1);
    }
    for (int t = 0; t < m; t++) {
        long x = XC1[t] + XC2[t];
        long lx = D[t] ? lim_on : offX;
        if (x > lx) e += (x - lx);
    }
    return e;
}

/* toggle a single position in U0 (delta = +-1), updating arrays O(m) */
static void toggleU0(int p) {
    int delta = U0[p] ? -1 : 1;
    /* AU0(d) = sum_e U0(e)U0(e-d): changing U0(p): += delta*(U0(p-d)+U0(p+d)) for d!=0,
       and AU0(0) += delta (self, since U0(p)U0(p) term) -- with delta sign handling:
       do updates using OLD values then flip. */
    for (int d = 1; d < m; d++) {
        int a = p - d; if (a < 0) a += m;
        int b = p + d; if (b >= m) b -= m;
        AU0[d] += (long)delta * (U0[a] + U0[b]);
    }
    AU0[0] += delta; /* U0(p)^2 term */
    for (int t = 0; t < m; t++) {
        int td = t - p; if (td < 0) td += m; /* XC1(t) has term D(t-p)*U0(p)?? careful:
            XC1(t) = sum_d D(d)U0(t-d): term with t-d=p: d = t-p: D(t-p)*U0(p) */
        XC1[t] += (long)delta * D[td];
    }
    U0[p] ^= 1;
}
static void toggleU1(int p) {
    int delta = U1[p] ? -1 : 1;
    for (int d = 1; d < m; d++) {
        int a = p - d; if (a < 0) a += m;
        int b = p + d; if (b >= m) b -= m;
        AU1[d] += (long)delta * (U1[a] + U1[b]);
    }
    AU1[0] += delta;
    /* XC2(t) = sum_u U1(u)D(t+u): term u=p: D(t+p) */
    for (int t = 0; t < m; t++) {
        int tp = t + p; if (tp >= m) tp -= m;
        XC2[t] += (long)delta * D[tp];
    }
    U1[p] ^= 1;
}
static void toggleD(int p) {
    int delta = D[p] ? -1 : 1;
    for (int d = 1; d < m; d++) {
        int a = p - d; if (a < 0) a += m;
        int b = p + d; if (b >= m) b -= m;
        AD[d] += (long)delta * (D[a] + D[b]);
    }
    AD[0] += delta;
    /* XC1(t): term d=p: U0(t-p) */
    for (int t = 0; t < m; t++) {
        int tp = t - p; if (tp < 0) tp += m;
        XC1[t] += (long)delta * U0[tp];
    }
    /* XC2(t) = sum_u U1(u) D(t+u): term t+u=p: u = p-t: U1(p-t) */
    for (int t = 0; t < m; t++) {
        int pt = p - t; if (pt < 0) pt += m;
        XC2[t] += (long)delta * U1[pt];
    }
    D[p] ^= 1;
}

static void emit(const char *prefix, int seed, int r) {
    char fn[256];
    snprintf(fn, sizeof fn, "%s_n%d_seed%d_r%d.sets", prefix, n, seed, r);
    FILE *f = fopen(fn, "w");
    fprintf(f, "n %d m %d\nU0 ", n, m);
    for (int i = 0; i < m; i++) fprintf(f, "%d", U0[i]);
    fprintf(f, "\nU1 ");
    for (int i = 0; i < m; i++) fprintf(f, "%d", U1[i]);
    fprintf(f, "\nD ");
    for (int i = 0; i < m; i++) fprintf(f, "%d", D[i]);
    fprintf(f, "\n");
    fclose(f);
    printf("EMITTED %s\n", fn);
    fflush(stdout);
}

int main(int argc, char **argv) {
    n = atoi(argv[1]);
    long long iters = atoll(argv[2]);
    int seed = argc > 3 ? atoi(argv[3]) : 1;
    int restarts = argc > 4 ? atoi(argv[4]) : 8;
    const char *prefix = argc > 5 ? argv[5] : "ramsey_c_witness";
    N = 4 * n - 2; m = 2 * n - 1;
    long bestglobal = 1L << 40;
    for (int r = 0; r < restarts; r++) {
        rng_s = 0x9E3779B97F4A7C15ULL * (seed + 1000 * r + 1);
        memset(U0, 0, sizeof U0); memset(U1, 0, sizeof U1); memset(D, 0, sizeof D);
        /* size-targeted random init: |U| = utarget (even), |D| = dtarget
           (defaults chosen by budget-slack scan; override via env) */
        int utarget = 94, dtarget = 104;
        {
            const char *eu = getenv("RB_U"), *ed = getenv("RB_D");
            if (eu) utarget = atoi(eu);
            if (ed) dtarget = atoi(ed);
        }
        int placed = 0;
        while (placed < utarget / 2) {
            int a = 1 + (int)(rnd() % (m - 1));
            if (a > m / 2) a = m - a;
            if (!U0[a]) { U0[a] = U0[m - a] = 1; placed++; }
        }
        placed = 0;
        while (placed < utarget / 2) {
            int a = 1 + (int)(rnd() % (m - 1));
            if (a > m / 2) a = m - a;
            if (!U1[a]) { U1[a] = U1[m - a] = 1; placed++; }
        }
        placed = 0;
        while (placed < dtarget) {
            int i = (int)(rnd() % m);
            if (!D[i]) { D[i] = 1; placed++; }
        }
        recompute_all();
        long e = energy(), best = e;
        double T0 = 6.0;
        long long since = 0;
        for (long long it = 0; it < iters; it++) {
            double frac = (double)(it % 40000000LL) / 40000000.0; /* reheat cycles */
            double T = T0 * (1.0 - frac) + 0.02;
            int mv = (int)(rnd() % 100);
            int p;
            long olde = e;
            if (0) {
                /* guided move: find worst violating class, toggle a contributor */
                int su0 = 0, su1 = 0, sd = 0;
                for (int i = 0; i < m; i++) { su0 += U0[i]; su1 += U1[i]; sd += D[i]; }
                int k0 = su0 + sd, k1 = su1 + sd;
                int lim_on = n - 2;
                int off0 = n + 1 + 2 * k0 - N, off1 = n + 1 + 2 * k1 - N, offX = n + 1 + k0 + k1 - N;
                long worst = 0; int wd = -1, wt = 0; /* wt: 0=within0,1=within1,2=cross */
                for (int d = 1; d < m; d++) {
                    long w0 = AU0[d] + AD[d] - (U0[d] ? lim_on : off0);
                    if (w0 > worst) { worst = w0; wd = d; wt = 0; }
                    long w1 = AU1[d] + AD[d] - (U1[d] ? lim_on : off1);
                    if (w1 > worst) { worst = w1; wd = d; wt = 1; }
                }
                for (int t = 0; t < m; t++) {
                    long x = XC1[t] + XC2[t] - (D[t] ? lim_on : offX);
                    if (x > worst) { worst = x; wd = t; wt = 2; }
                }
                if (wd < 0) { continue; }
                /* pick a random contributor to that class and toggle it */
                int tries = 0, done = 0;
                while (tries++ < 40 && !done) {
                    int eidx = (int)(rnd() % m);
                    int ed = eidx - wd; if (ed < 0) ed += m;
                    if (wt == 0) {
                        /* contributor to AU0(wd): U0(e)&U0(e-wd); or to AD via D */
                        if (U0[eidx] && U0[ed]) {
                            p = (rnd() & 1) ? eidx : ed;
                            toggleU0(p); if (p != (m - p) % m) toggleU0((m - p) % m);
                            done = 1;
                        } else if (D[eidx] && D[ed]) {
                            p = (rnd() & 1) ? eidx : ed;
                            toggleD(p); done = 1;
                        }
                    } else if (wt == 1) {
                        if (U1[eidx] && U1[ed]) {
                            p = (rnd() & 1) ? eidx : ed;
                            toggleU1(p); if (p != (m - p) % m) toggleU1((m - p) % m);
                            done = 1;
                        } else if (D[eidx] && D[ed]) {
                            p = (rnd() & 1) ? eidx : ed;
                            toggleD(p); done = 1;
                        }
                    } else {
                        /* cross class wd: XC1 contributors: D(dd)&U0(wd-dd);
                           XC2: U1(u)&D(wd+u) */
                        int u = eidx;
                        int a = wd - u; if (a < 0) a += m;
                        int b2 = wd + u; if (b2 >= m) b2 -= m;
                        if (D[u] && U0[a]) {
                            if (rnd() & 1) { toggleD(u); }
                            else { toggleU0(a); if (a != (m - a) % m) toggleU0((m - a) % m); }
                            done = 1;
                        } else if (U1[u] && D[b2]) {
                            if (rnd() & 1) { toggleD(b2); }
                            else { toggleU1(u); if (u != (m - u) % m) toggleU1((m - u) % m); }
                            done = 1;
                        }
                    }
                }
                if (!done) continue;
                e = energy();
                if (!(e <= olde || rndu() < exp((olde - e) / T))) {
                    /* revert is complex; recompute via re-toggle is not tracked —
                       accept-only-if-better policy for guided moves: revert by
                       full state restore is avoided; instead we allow uphill with
                       the same Metropolis rule but must revert: we cannot — so
                       guided moves are ACCEPT-ALWAYS (they remove a violator). */
                    /* keep it: guided moves always accepted */
                }
                if (e < best) { best = e; since = 0; } else since++;
                if (e == 0) {
                    printf("[r%d] ZERO at it=%lld (guided)\n", r, it);
                    emit(prefix, seed, r);
                    return 0;
                }
                continue;
            }
            if (mv < 27) {
                /* U0 pair-swap: move one symmetric pair out, another in */
                int a = 0, b = 0, tries = 0;
                while (tries++ < 60) { a = 1 + (int)(rnd() % (m - 1)); if (U0[a]) break; }
                while (tries++ < 120) { b = 1 + (int)(rnd() % (m - 1)); if (!U0[b] && b != a && b != m - a) break; }
                if (!U0[a] || U0[b]) continue;
                toggleU0(a); if (a != m - a) toggleU0(m - a);
                toggleU0(b); if (b != m - b) toggleU0(m - b);
                e = energy();
                if (!(e <= olde || rndu() < exp((olde - e) / T))) {
                    toggleU0(b); if (b != m - b) toggleU0(m - b);
                    toggleU0(a); if (a != m - a) toggleU0(m - a); e = olde;
                }
            } else if (mv < 54) {
                int a = 0, b = 0, tries = 0;
                while (tries++ < 60) { a = 1 + (int)(rnd() % (m - 1)); if (U1[a]) break; }
                while (tries++ < 120) { b = 1 + (int)(rnd() % (m - 1)); if (!U1[b] && b != a && b != m - a) break; }
                if (!U1[a] || U1[b]) continue;
                toggleU1(a); if (a != m - a) toggleU1(m - a);
                toggleU1(b); if (b != m - b) toggleU1(m - b);
                e = energy();
                if (!(e <= olde || rndu() < exp((olde - e) / T))) {
                    toggleU1(b); if (b != m - b) toggleU1(m - b);
                    toggleU1(a); if (a != m - a) toggleU1(m - a); e = olde;
                }
            } else if (mv < 90) {
                /* D swap: out one, in one */
                int a = 0, b = 0, tries = 0;
                while (tries++ < 60) { a = (int)(rnd() % m); if (D[a]) break; }
                while (tries++ < 120) { b = (int)(rnd() % m); if (!D[b]) break; }
                if (!D[a] || D[b]) continue;
                toggleD(a); toggleD(b);
                e = energy();
                if (!(e <= olde || rndu() < exp((olde - e) / T))) {
                    toggleD(b); toggleD(a); e = olde;
                }
            } else if (mv < 93) {
                p = 1 + (int)(rnd() % (m - 1));
                toggleU0(p); if (p != m - p) toggleU0(m - p);
                e = energy();
                if (!(e <= olde || rndu() < exp((olde - e) / T))) {
                    toggleU0(p); if (p != m - p) toggleU0(m - p); e = olde;
                }
            } else if (mv < 96) {
                p = 1 + (int)(rnd() % (m - 1));
                toggleU1(p); if (p != m - p) toggleU1(m - p);
                e = energy();
                if (!(e <= olde || rndu() < exp((olde - e) / T))) {
                    toggleU1(p); if (p != m - p) toggleU1(m - p); e = olde;
                }
            } else {
                p = (int)(rnd() % m);
                toggleD(p);
                e = energy();
                if (!(e <= olde || rndu() < exp((olde - e) / T))) {
                    toggleD(p); e = olde;
                }
            }
            if (e < best) { best = e; since = 0; } else since++;
            if (e == 0) {
                printf("[r%d] ZERO at it=%lld\n", r, it);
                emit(prefix, seed, r);
                return 0;
            }
            if ((it & ((1 << 24) - 1)) == 0) {
                long echk = e;
                recompute_all();
                long efull = energy();
                if (efull != echk) {
                    printf("[r%d] INCREMENTAL DRIFT at it=%lld: %ld vs %ld — rebuilt\n",
                           r, it, echk, efull);
                    e = efull;
                }
                printf("[r%d] it=%lld e=%ld best=%ld T=%.2f\n", r, it, e, best, T);
                fflush(stdout);
            }
        }
        printf("[r%d] done best=%ld\n", r, best);
        fflush(stdout);
        if (best < bestglobal) bestglobal = best;
    }
    printf("NO WITNESS best=%ld\n", bestglobal);
    return 0;
}
