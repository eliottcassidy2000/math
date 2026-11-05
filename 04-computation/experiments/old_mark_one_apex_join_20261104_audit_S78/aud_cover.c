/*
 * aud_cover.c  (S78 audit, independent re-implementation; does not use the session's g routine or suffix tables)
 *
 * Coverage of the odd numbers N <= X = 2^LOGX by
 *   F0: h of every bank triple, and H(T(w)) = g(v) for v <= X;
 *   F1: one-arc joins of two one-apex tournaments, H = (1+u)(1+v) + g(u) g(v) - e(u) e(v);
 *   F2: one-arc joins of a one-apex tournament with a bank rooted tournament, H = g(v) y_h - e(v) y_e + (1+v) y_R.
 * Here e(v) = [v even] (e(0) = 1).
 *
 * g is computed from the auditor's own path-count formula for the one-apex tournament T(w) (chain 1->...->n, apex x,
 * x -> i iff w_i = 1): a Hamiltonian path is Q1 x Q2 with Q1, Q2 increasing, so
 *     H(T(w)) = [w_1 = 1] + [w_n = 0] + #{a : w_a = 0, w_(a+1) = 1} + sum_{b < a, w_b = 1, w_a = 0} 2^(a-b-1),
 * i.e. for v >= 1 (standard binary, MSB first)  g(v) = 1 + [v even] + #01(v) + sum_{zero bits i of v} (v >> (i+1)).
 * Fast evaluation (derived here from that formula): v = k 2^16 + r, k >= 1:
 *     g(v) = g(k) + k Z(r) + Q(r) - [k even][bit 15 of r = 0],
 *     Z(r) = sum_{zero bits i < 16 of r} 2^(15-i),  Q(r) = [r even] + #01(16-bit r) + sum_{zero bits i < 16} (r >> (i+1)).
 * The fast form is checked against the direct formula, and the direct formula against the plain recursion
 * g(2v+1) = g(v), g(2v) = g(v) + v + (v mod 2), before the sieve runs.
 *
 * Coverage marks bytes (one byte per odd N; relaxed atomic stores of 1, so no lost updates).
 * Work is split into (row, v-chunk) items so the long rows (u = 0, small y_R) are parallel.
 * Usage: aud_cover LOGX bankfile
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <omp.h>
typedef long long ll;
typedef unsigned long long ull;

static ll g_direct(ull v) {                    /* auditor's closed formula */
    if (v == 0) return 1;
    int len = 64 - __builtin_clzll(v);
    ll s = 1 + ((v & 1) == 0);
    for (int i = 0; i + 1 < len; i++) {        /* "01" pattern: bit i+1 = 0, bit i = 1 (MSB-first reading) */
        if (!((v >> (i + 1)) & 1) && ((v >> i) & 1)) s++;
    }
    for (int i = 0; i < len; i++) if (!((v >> i) & 1)) s += (ll)(v >> (i + 1));
    return s;
}
static ll g_rec_def(ull v) { if (v == 0) return 1; if (v & 1) return g_rec_def(v >> 1); ull u = v >> 1; return g_rec_def(u) + (ll)u + (ll)(u & 1); }

static uint32_t GK[1 << 16], ZT[1 << 16], QT[1 << 16]; static uint8_t B15Z[1 << 16];
static inline ll gfast(ull v) {
    ull k = v >> 16, r = v & 0xFFFF;
    if (k == 0) return GK[r];
    return (ll)GK[k] + (ll)k * ZT[r] + QT[r] - ((k & 1) == 0 ? B15Z[r] : 0);
}

static uint8_t *covA, *covF1; static ll X;
static inline void mark(uint8_t *c, ll N) { if (N >= 1 && N <= X && (N & 1)) __atomic_store_n(&c[N >> 1], (uint8_t)1, __ATOMIC_RELAXED); }

typedef struct { ll row; ll v0, v1; } Item;   /* row: u (F1) or triple index (F2); v in [v0, v1) */

int main(int argc, char **argv) {
    int LOGX = argc > 1 ? atoi(argv[1]) : 24; X = 1LL << LOGX;
    const char *fn = argc > 2 ? argv[2] : "aud_bank7.txt";
    if (LOGX > 32) { fprintf(stderr, "LOGX <= 32 (16+16 split)\n"); return 1; }
    FILE *f = fopen(fn, "r"); if (!f) { fprintf(stderr, "no bank %s\n", fn); return 1; }
    int cap = 4096, nb = 0; ll *yh = malloc(sizeof(ll) * cap), *ye = malloc(sizeof(ll) * cap), *yR = malloc(sizeof(ll) * cap); ll a, b, c;
    while (fscanf(f, "%lld %lld %lld", &a, &b, &c) == 3) { if (nb == cap) { cap *= 2; yh = realloc(yh, sizeof(ll) * cap); ye = realloc(ye, sizeof(ll) * cap); yR = realloc(yR, sizeof(ll) * cap); } yh[nb] = a; ye[nb] = b; yR[nb] = c; nb++; }
    fclose(f);
    double t0 = omp_get_wtime();
    printf("# aud_cover LOGX=%d bank=%s triples=%d threads=%d\n", LOGX, fn, nb, omp_get_max_threads());
    /* tables */
    for (ull r = 0; r < (1u << 16); r++) {
        GK[r] = (uint32_t)g_direct(r);
        ll Z = 0, Q = ((r & 1) == 0), S = 0, c01 = 0;
        for (int i = 0; i < 16; i++) if (!((r >> i) & 1)) { Z += 1LL << (15 - i); S += (ll)(r >> (i + 1)); }
        for (int i = 0; i < 15; i++) if (!((r >> (i + 1)) & 1) && ((r >> i) & 1)) c01++;
        Q += c01 + S; ZT[r] = (uint32_t)Z; QT[r] = (uint32_t)Q; B15Z[r] = !((r >> 15) & 1);
    }
    /* self-checks: direct formula vs recursion (exhaustive to 2^22, plus random), fast vs direct (exhaustive to 2^24, random to X) */
    {
        ll bad1 = 0, bad2 = 0, n1 = 0, n2 = 0;
        for (ull v = 0; v < (1u << 22); v++) { n1++; if (g_direct(v) != g_rec_def(v)) bad1++; }
        #pragma omp parallel for reduction(+:bad2,n2) schedule(static)
        for (ll v = 0; v < (1LL << 24); v++) { n2++; if (gfast((ull)v) != g_direct((ull)v)) bad2++; }
        ull s = 0x243F6A8885A308D3ULL; ll bad3 = 0, bad4 = 0;
        for (int t = 0; t < 1000000; t++) {
            s ^= s << 13; s ^= s >> 7; s ^= s << 17;
            ull v = s % (ull)(X + 1);
            if (gfast(v) != g_direct(v)) bad3++;
            if (g_direct(v) != g_rec_def(v)) bad4++;
        }
        printf("self-check: direct vs recursion v<2^22: %lld/%lld bad; fast vs direct v<2^24: %lld/%lld bad; random v<=X: fast-direct %lld bad, direct-recursion %lld bad (1e6 samples)\n",
               bad1, n1, bad2, n2, bad3, bad4);
        if (bad1 || bad2 || bad3 || bad4) return 2;
    }
    size_t nbytes = (size_t)(X >> 1) + 2;
    covA = calloc(nbytes, 1); covF1 = calloc(nbytes, 1);
    if (!covA || !covF1) { fprintf(stderr, "alloc\n"); return 1; }
    /* F0 (compile with -DNOF0 to test whether F1 + F2 alone suffice) */
#ifndef NOF0
    for (int i = 0; i < nb; i++) mark(covA, yh[i]);
    #pragma omp parallel for schedule(static)
    for (ll v = 0; v <= X; v++) { ll gv = gfast((ull)v); if (gv <= X) mark(covA, gv); }
    printf("F0 done %.1fs\n", omp_get_wtime() - t0); fflush(stdout);
#else
    printf("F0 SKIPPED (NOF0 build)\n");
#endif
    const ll CH = 1LL << 22;
    /* F1 work list */
    {
        size_t ni = 0, capi = 1 << 16; Item *it = malloc(sizeof(Item) * capi);
        for (ll u = 0; (1 + u) * (1 + u) <= X; u++) {
            ll vmax = X / (1 + u) - 1;            /* (1+u)(1+v) <= X */
            for (ll v0 = u; v0 <= vmax; v0 += CH) {
                if (ni == capi) { capi *= 2; it = realloc(it, sizeof(Item) * capi); }
                it[ni].row = u; it[ni].v0 = v0; it[ni].v1 = (v0 + CH - 1 < vmax ? v0 + CH - 1 : vmax) + 1; ni++;
            }
        }
        ll cnt = 0;
        #pragma omp parallel for schedule(dynamic, 1) reduction(+:cnt)
        for (size_t j = 0; j < ni; j++) {
            ll u = it[j].row, gu = gfast((ull)u), eu = ((u & 1) == 0);
            for (ll v = it[j].v0; v < it[j].v1; v++) {
                ll N = (1 + u) * (1 + v) + gu * gfast((ull)v) - eu * ((v & 1) == 0);
                cnt++;
                mark(covF1, N);
            }
        }
        printf("F1 done %.1fs (%zu work items, %lld pairs)\n", omp_get_wtime() - t0, ni, cnt); fflush(stdout);
        free(it);
    }
    /* F2 work list */
    {
        size_t ni = 0, capi = 1 << 16; Item *it = malloc(sizeof(Item) * capi);
        for (int i = 0; i < nb; i++) {
            ll vmax = X / yR[i] - 1;              /* (1+v) y_R <= X */
            for (ll v0 = 0; v0 <= vmax; v0 += CH) {
                if (ni == capi) { capi *= 2; it = realloc(it, sizeof(Item) * capi); }
                it[ni].row = i; it[ni].v0 = v0; it[ni].v1 = (v0 + CH - 1 < vmax ? v0 + CH - 1 : vmax) + 1; ni++;
            }
        }
        ll cnt = 0;
        #pragma omp parallel for schedule(dynamic, 1) reduction(+:cnt)
        for (size_t j = 0; j < ni; j++) {
            int i = (int)it[j].row; ll h = yh[i], e = ye[i], R = yR[i];
            for (ll v = it[j].v0; v < it[j].v1; v++) {
                ll N = gfast((ull)v) * h - (((v & 1) == 0) ? e : 0) + (1 + v) * R;
                cnt++;
                mark(covA, N);
            }
        }
        printf("F2 done %.1fs (%zu work items, %lld evaluations)\n", omp_get_wtime() - t0, ni, cnt); fflush(stdout);
        free(it);
    }
    /* misses */
    ll miss = 0;
    printf("odd N in [1, 2^%d] covered by none of F0, F1, F2:", LOGX);
    for (ll N = 1; N <= X; N += 2) if (!covA[N >> 1] && !covF1[N >> 1]) { miss++; if (miss <= 100) printf(" %lld", N); }
    printf("  [total %lld]\n", miss);
    ll missF0F1 = 0; for (ll N = 1; N <= X; N += 2) if (!covF1[N >> 1]) missF0F1++;
    printf("F1-only misses per dyadic block (2^k, 2^(k+1)): k count ratio density\n");
    ll prev = -1;
    for (int k = 10; k < LOGX; k++) {
        ll lo = 1LL << k, hi = 1LL << (k + 1), m = 0;
        for (ll N = lo + 1; N < hi; N += 2) if (!covF1[N >> 1]) m++;
        printf("  k=%d %lld %.4f %.6f\n", k, m, prev > 0 ? (double)m / prev : 0.0, (double)m / (double)(lo / 2));
        prev = m;
    }
    printf("total F1-only odd misses in [1, X]: %lld; elapsed %.1fs\n", missF0F1, omp_get_wtime() - t0);
    return 0;
}
