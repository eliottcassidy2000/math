/*
 * aud_unitstart.c  (independent audit of THM-4872 (A) (vi)-(vii), opus-2026-11-04-S76 audit)
 *
 * Written from scratch; it does not use the FUNNELS F1 component decomposition.
 *
 * Normal form.  If (T, r) has exactly one Hamiltonian path starting at r, relabel the vertices 0..n-1 along that path.
 * A rooted isomorphism maps the unique path from the root onto the unique path from the root, so rooted classes with
 * s = 1 are in bijection with labelled tournaments on {0..n-1} containing the arcs i -> i+1 whose only Hamiltonian path
 * from 0 is 0, 1, ..., n-1 ("path-order tournaments").
 *
 * Generator (root prepending).  Every Hamiltonian path from 0 is 0 followed by a Hamiltonian path of T - 0 starting at an
 * out-neighbour of 0, so  s_T(0) = sum over j in N+(0) of s_(T-0)(j).  Since 0 -> 1 and s_(T-0)(1) >= 1, s_T(0) = 1 iff
 * s_(T-0)(1) = 1 and every other out-neighbour j of 0 has s_(T-0)(j) = 0.  Hence the path-order tournaments of order n
 * are exactly: a path-order tournament P of order n - 1 (relabelled 1..n-1), a new root 0 -> 1, the root beating an
 * arbitrary subset of Z(P) = {j : s_P(j) = 0}, and every other vertex beating the root.  Each child determines its parent
 * and its subset, so there is no duplication.  This is complete and exact at every order; strongness is tested
 * afterwards by reachability.
 *
 * Cross-check: at orders <= 8 the generated sets are compared, element by element, with a brute-force scan of all
 * 2^(C(n,2) - (n-1)) labelled tournaments containing the path 0 -> 1 -> ... -> n-1.
 *
 * Output per order: all path-order classes, strong ones, the profile (h, s, q, F, e) of each strong one, and all
 * +-2 port pairs (F = F', h' = h + 2) with their profiles.
 *
 * Build: gcc -O3 -march=native -fopenmp -o aud_unitstart aud_unitstart.c ; run: ./aud_unitstart 14
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>

#define MAXN 16
typedef long long ll;
typedef struct { uint16_t out[MAXN]; uint16_t zmask; uint8_t strong; } Obj;
typedef struct { ll h, s, q, F, e; size_t idx; } Prof;

static Obj *LV[MAXN + 1]; static size_t LN[MAXN + 1];

static int strong_n(int n, const uint16_t *out) {
    unsigned full = (1u << n) - 1, f = 1, ch = 1;
    while (ch) { ch = 0; for (int x = 0; x < n; x++) if ((f >> x) & 1u) { unsigned nf = f | out[x]; if (nf != f) { f = nf; ch = 1; } } }
    if (f != full) return 0;
    unsigned g = 1; ch = 1;
    while (ch) { ch = 0; for (int y = 0; y < n; y++) if (!((g >> y) & 1u) && (out[y] & g)) { g |= 1u << y; ch = 1; } }
    return g == full;
}

/* start counts s_T(x) for all x, by the end-count DP on the converse */
static void start_counts(int n, const uint16_t *out, ll *dp, ll *sc) {
    unsigned full = (1u << n) - 1;
    uint16_t conv[MAXN];
    for (int x = 0; x < n; x++) { conv[x] = 0; for (int y = 0; y < n; y++) if (y != x && !((out[x] >> y) & 1u)) conv[x] |= (uint16_t)(1u << y); }
    memset(dp, 0, sizeof(ll) * ((size_t)1 << n) * n);
    for (int x = 0; x < n; x++) dp[(size_t)(1u << x) * n + x] = 1;
    for (unsigned S = 1; S <= full; S++) for (int x = 0; x < n; x++) {
        ll c = dp[(size_t)S * n + x]; if (!c) continue;
        unsigned nb = conv[x] & ~S;
        while (nb) { int y = __builtin_ctz(nb); nb &= nb - 1; dp[(size_t)(S | (1u << y)) * n + y] += c; }
    }
    for (int x = 0; x < n; x++) sc[x] = dp[(size_t)full * n + x];
}

static Prof profile(int n, const uint16_t *out, ll *dp, ll *HS) {
    unsigned full = (1u << n) - 1;
    memset(dp, 0, sizeof(ll) * ((size_t)1 << n) * n);
    for (int x = 0; x < n; x++) dp[(size_t)(1u << x) * n + x] = 1;
    for (unsigned S = 1; S <= full; S++) for (int x = 0; x < n; x++) {
        ll c = dp[(size_t)S * n + x]; if (!c) continue;
        unsigned nb = out[x] & ~S;
        while (nb) { int y = __builtin_ctz(nb); nb &= nb - 1; dp[(size_t)(S | (1u << y)) * n + y] += c; }
    }
    HS[0] = 1;
    for (unsigned S = 1; S <= full; S++) { ll t = 0; for (int x = 0; x < n; x++) t += dp[(size_t)S * n + x]; HS[S] = t; }
    Prof p; p.h = HS[full]; p.q = HS[full & ~1u]; p.e = dp[(size_t)full * n + 0];
    ll F = 0; for (unsigned S = 1; S <= full; S += 2) F += dp[(size_t)S * n + 0] * HS[full ^ S];
    p.F = F;
    return p;
}

static ll fib(int k) { ll a = 0, b = 1; for (int i = 0; i < k; i++) { ll t = a + b; a = b; b = t; } return a; }

static int cmpprof(const void *a, const void *b) {
    const Prof *x = a, *y = b;
    if (x->F != y->F) return x->F < y->F ? -1 : 1;
    if (x->h != y->h) return x->h < y->h ? -1 : 1;
    if (x->q != y->q) return x->q < y->q ? -1 : 1;
    return 0;
}
static int cmpu64(const void *a, const void *b) { uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b; return x < y ? -1 : x > y; }

/* encode the off-path arcs of a path-order tournament (pairs i < j, j >= i + 2; bit = 1 iff i -> j) */
static uint64_t offcode(int n, const uint16_t *out) {
    uint64_t c = 0; int k = 0;
    for (int i = 0; i < n; i++) for (int j = i + 2; j < n; j++) { if ((out[i] >> j) & 1u) c |= 1ull << k; k++; }
    return c;
}

int main(int argc, char **argv) {
    int NMAX = argc > 1 ? atoi(argv[1]) : 12;
    if (NMAX > 15) NMAX = 15;
    printf("# aud_unitstart NMAX=%d (root-prepending generator of path-order unit-start rooted tournaments)\n", NMAX);
    /* order 1 */
    LV[1] = calloc(1, sizeof(Obj)); LN[1] = 1; LV[1][0].zmask = 0; LV[1][0].strong = 1;
    for (int n = 2; n <= NMAX; n++) {
        /* children counts */
        size_t total = 0;
        for (size_t i = 0; i < LN[n - 1]; i++) total += (size_t)1 << __builtin_popcount(LV[n - 1][i].zmask);
        LV[n] = malloc(sizeof(Obj) * total); LN[n] = total;
        size_t pos = 0;
        for (size_t i = 0; i < LN[n - 1]; i++) {
            const Obj *P = &LV[n - 1][i];
            unsigned Z = P->zmask;   /* labels of P (0..n-2) */
            unsigned sub = 0;
            do {
                Obj *T = &LV[n][pos++]; memset(T, 0, sizeof *T);
                for (int x = 0; x < n - 1; x++) T->out[x + 1] = (uint16_t)(P->out[x] << 1);
                T->out[0] = (uint16_t)(1u << 1);
                for (int j = 1; j < n - 1; j++) {
                    if ((sub >> j) & 1u) T->out[0] |= (uint16_t)(1u << (j + 1));
                    else T->out[j + 1] |= 1u;          /* j+1 -> 0 */
                }
                sub = (sub - Z) & Z;                    /* next subset of Z */
            } while (sub != 0);
        }
        if (pos != total) { fprintf(stderr, "count mismatch\n"); return 1; }
        /* start counts (for Z), strongness, sanity s_T(0) = 1 */
        long bad = 0, nstrong = 0;
        #pragma omp parallel
        {
            ll *dp = malloc(sizeof(ll) * ((size_t)1 << n) * n); ll sc[MAXN];
            #pragma omp for schedule(dynamic, 64) reduction(+:bad,nstrong)
            for (size_t i = 0; i < total; i++) {
                Obj *T = &LV[n][i];
                for (int x = 0; x + 1 < n; x++) if (!((T->out[x] >> (x + 1)) & 1u)) bad++;
                for (int x = 0; x < n; x++) for (int y = x + 1; y < n; y++) if ((((T->out[x] >> y) & 1u) + ((T->out[y] >> x) & 1u)) != 1) bad++;
                start_counts(n, T->out, dp, sc);
                if (sc[0] != 1) bad++;
                unsigned z = 0; for (int x = 0; x < n; x++) if (sc[x] == 0) z |= 1u << x;
                T->zmask = (uint16_t)z;
                T->strong = (uint8_t)strong_n(n, T->out);
                nstrong += T->strong;
            }
            free(dp);
        }
        /* brute-force cross-check of the whole set at small orders */
        char bf[200] = "";
        if (n >= 2 && n <= 8) {
            int npairs = 0; for (int i = 0; i < n; i++) for (int j = i + 2; j < n; j++) npairs++;
            uint64_t NB = 1ull << npairs; uint64_t *found = malloc(sizeof(uint64_t) * (total + 16)); size_t nf = 0;
            #pragma omp parallel
            {
                ll *dp = malloc(sizeof(ll) * ((size_t)1 << n) * n);
                #pragma omp for schedule(dynamic, 1024)
                for (uint64_t c = 0; c < NB; c++) {
                    uint16_t out[MAXN]; memset(out, 0, sizeof out); int k = 0;
                    for (int i = 0; i + 1 < n; i++) out[i] |= (uint16_t)(1u << (i + 1));
                    for (int i = 0; i < n; i++) for (int j = i + 2; j < n; j++) { if ((c >> k) & 1ull) out[i] |= (uint16_t)(1u << j); else out[j] |= (uint16_t)(1u << i); k++; }
                    /* paths from 0 */
                    unsigned full = (1u << n) - 1;
                    memset(dp, 0, sizeof(ll) * ((size_t)1 << n) * n);
                    dp[1 * n + 0] = 1;
                    for (unsigned S = 1; S <= full; S += 2) for (int x = 0; x < n; x++) {
                        ll v = dp[(size_t)S * n + x]; if (!v) continue;
                        unsigned nb = out[x] & ~S;
                        while (nb) { int y = __builtin_ctz(nb); nb &= nb - 1; dp[(size_t)(S | (1u << y)) * n + y] += v; }
                    }
                    ll s = 0; for (int x = 0; x < n; x++) s += dp[(size_t)full * n + x];
                    if (s == 1) {
                        size_t slot;
                        #pragma omp atomic capture
                        slot = nf++;
                        if (slot < total + 16) found[slot] = c;
                    }
                }
                free(dp);
            }
            uint64_t *gen = malloc(sizeof(uint64_t) * total);
            for (size_t i = 0; i < total; i++) gen[i] = offcode(n, LV[n][i].out);
            qsort(gen, total, sizeof(uint64_t), cmpu64);
            int same = (nf == total);
            if (same) { qsort(found, nf, sizeof(uint64_t), cmpu64); for (size_t i = 0; i < nf; i++) if (found[i] != gen[i]) { same = 0; break; } }
            int dup = 0; for (size_t i = 1; i < total; i++) if (gen[i] == gen[i - 1]) dup++;
            snprintf(bf, sizeof bf, " | brute force over 2^%d: %zu unit-start, generated set identical: %s, duplicates %d", npairs, nf, same ? "yes" : "NO", dup);
            free(found); free(gen);
        }
        printf("n=%2d unit-start rooted classes %zu (F(2n-3) = %lld), strong %ld (F(2n-4) = %lld), sanity failures %ld%s\n",
               n, total, fib(2 * n - 3), nstrong, fib(2 * n - 4), bad, bf);
        fflush(stdout);
        if (n < 3) continue;
        /* profiles of the strong ones */
        Prof *PR = malloc(sizeof(Prof) * nstrong); size_t np = 0;
        size_t *ids = malloc(sizeof(size_t) * nstrong);
        for (size_t i = 0; i < total; i++) if (LV[n][i].strong) ids[np++] = i;
        #pragma omp parallel
        {
            ll *dp = malloc(sizeof(ll) * ((size_t)1 << n) * n), *HS = malloc(sizeof(ll) * ((size_t)1 << n));
            #pragma omp for schedule(dynamic, 16)
            for (size_t k = 0; k < np; k++) { PR[k] = profile(n, LV[n][ids[k]].out, dp, HS); PR[k].s = 1; PR[k].idx = ids[k]; }
            free(dp); free(HS);
        }
        qsort(PR, np, sizeof(Prof), cmpprof);
        long pairs = 0, pairs_q = 0; long distinctF = 0;
        for (size_t a = 0; a < np; a++) {
            if (a == 0 || PR[a].F != PR[a - 1].F) distinctF++;
            for (size_t b = a + 1; b < np && PR[b].F == PR[a].F; b++) {
                if (PR[b].h == PR[a].h + 2) {
                    pairs++; if (PR[b].q == PR[a].q) pairs_q++;
                    const uint16_t *oa = LV[n][PR[a].idx].out, *ob = LV[n][PR[b].idx].out;
                    printf("   pair F=%lld: (h,q,e) = (%lld,%lld,%lld) / (%lld,%lld,%lld)%s  rows A:", PR[a].F, PR[a].h, PR[a].q, PR[a].e, PR[b].h, PR[b].q, PR[b].e,
                           PR[a].q == PR[b].q ? " [q equal]" : "");
                    for (int x = 0; x < n; x++) printf(" %u", oa[x]);
                    printf("  B:"); for (int x = 0; x < n; x++) printf(" %u", ob[x]); printf("\n");
                }
            }
        }
        printf("   n=%2d: distinct F %ld, +-2 port pairs %ld (with q = q' also: %ld)\n", n, distinctF, pairs, pairs_q);
        fflush(stdout);
        free(PR); free(ids);
    }
    return 0;
}
