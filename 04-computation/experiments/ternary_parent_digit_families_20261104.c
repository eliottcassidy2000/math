/*
 * ternary_parent_digit_families_20261104.c  (opus-2026-11-04-S77)
 *
 * Part 1 (any-mark C2 law).  For every labelled tournament B of order n <= NLAW and every vertex a, build the C2 child
 *   X = B + x (x -> a only; every other old vertex beats x),  T' = X + y (every old vertex beats y; y -> x),
 * and check H(T') = 3 h(B) + 2 sigma_a(B), sigma_a(B) = [P A_a](B) - h(B), [P A_a] = sum_{S subset V-a} h(B[S]) s_a(B[V-S]).
 *
 * Part 2 (an explicit digit-1 family).  For t >= 1 and bits b in {0,1}^t, the unit-start rooted tournament C(b):
 *   vertices r, w, z_1..z_t;  r -> w;  w -> every z_i;  z_i -> z_j (i < j);  r -> z_i iff b_i = 1.
 * Its converse B = C^op with the sole exit a -> r appended (all other vertices beat a) is digit-1 ready (sigma_a = e_r(B) =
 * s_r(C) = 1) with H = h(C) + F(C, r) - 1 (THM-4862 (a) applied to the converse; F = rooted two-slot count).  Each such H
 * gives the C2 child 3H + 2. This program checks only the digit; the identity for H is confirmed by the S77 audit's
 * closed form. The bits include non-strong members (b_t = 1), which FUNNELS F1 excludes; they add only H = 1.  Reports which odd numbers are attained as H, up to a bound, as t grows.
 *
 * Usage: ternfam NLAW TMAX BOUND
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>

typedef long long ll;
#define MAXV 24

static ll hpaths(int n, const uint32_t *out, ll *dp) {     /* number of Hamilton paths */
    int full = (1 << n) - 1;
    memset(dp, 0, sizeof(ll) * (size_t)(full + 1) * n);
    for (int x = 0; x < n; x++) dp[(size_t)(1 << x) * n + x] = 1;
    for (int S = 1; S <= full; S++) for (int x = 0; x < n; x++) { ll c = dp[(size_t)S * n + x]; if (!c) continue;
        uint32_t nb = out[x] & ~(uint32_t)S; while (nb) { int y = __builtin_ctz(nb); nb &= nb - 1; dp[(size_t)(S | (1 << y)) * n + y] += c; } }
    ll t = 0; for (int x = 0; x < n; x++) t += dp[(size_t)full * n + x]; return t;
}

/* subset tables: Hs[S] = h(B[S]); St[S] = paths of B[S] starting at a (a in S) */
static void sigma_parts(int n, const uint32_t *out, int a, ll *dp, ll *Hs, ll *Sa, ll *h_out, ll *pa_out) {
    int full = (1 << n) - 1;
    memset(dp, 0, sizeof(ll) * (size_t)(full + 1) * n);
    for (int x = 0; x < n; x++) dp[(size_t)(1 << x) * n + x] = 1;
    for (int S = 1; S <= full; S++) for (int x = 0; x < n; x++) { ll c = dp[(size_t)S * n + x]; if (!c) continue;
        uint32_t nb = out[x] & ~(uint32_t)S; while (nb) { int y = __builtin_ctz(nb); nb &= nb - 1; dp[(size_t)(S | (1 << y)) * n + y] += c; } }
    Hs[0] = 1; for (int S = 1; S <= full; S++) { ll t = 0; for (int x = 0; x < n; x++) t += dp[(size_t)S * n + x]; Hs[S] = t; }
    /* starts at a: dp from a */
    memset(dp, 0, sizeof(ll) * (size_t)(full + 1) * n);
    dp[(size_t)(1 << a) * n + a] = 1;
    for (int S = 1; S <= full; S++) { if (!((S >> a) & 1)) continue; for (int x = 0; x < n; x++) { ll c = dp[(size_t)S * n + x]; if (!c) continue;
        uint32_t nb = out[x] & ~(uint32_t)S; while (nb) { int y = __builtin_ctz(nb); nb &= nb - 1; dp[(size_t)(S | (1 << y)) * n + y] += c; } } }
    for (int S = 0; S <= full; S++) { ll t = 0; if ((S >> a) & 1) for (int x = 0; x < n; x++) t += dp[(size_t)S * n + x]; Sa[S] = t; }
    ll pa = 0; for (int S = 0; S <= full; S++) if (!((S >> a) & 1)) pa += Hs[S] * Sa[full & ~S];
    *h_out = Hs[full]; *pa_out = pa;
}

int main(int argc, char **argv) {
    int NLAW = argc > 1 ? atoi(argv[1]) : 6, TMAX = argc > 2 ? atoi(argv[2]) : 14;
    ll BOUND = argc > 3 ? atoll(argv[3]) : 20000;
    printf("# ternary_parent_digit_families_20261104 NLAW=%d TMAX=%d BOUND=%lld\n", NLAW, TMAX, BOUND);
    /* ---------- Part 1 ---------- */
    for (int n = 2; n <= NLAW; n++) {
        int narc = n * (n - 1) / 2; long bad = 0, cases = 0;
        #pragma omp parallel reduction(+:bad,cases)
        {
            ll *dp = malloc(sizeof(ll) * ((size_t)1 << (n + 2)) * (n + 2)), *Hs = malloc(sizeof(ll) << n), *Sa = malloc(sizeof(ll) << n);
            #pragma omp for schedule(dynamic, 64)
            for (long code = 0; code < (1L << narc); code++) {
                uint32_t out[MAXV] = {0}; int k = 0;
                for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++, k++) { if ((code >> k) & 1) out[i] |= 1u << j; else out[j] |= 1u << i; }
                for (int a = 0; a < n; a++) {
                    ll h, pa; sigma_parts(n, out, a, dp, Hs, Sa, &h, &pa);
                    uint32_t o2[MAXV]; for (int v = 0; v < n; v++) o2[v] = out[v];
                    int x = n, y = n + 1;
                    o2[x] = 1u << a;                                     /* x -> a only */
                    for (int v = 0; v < n; v++) if (v != a) o2[v] |= 1u << x;
                    o2[y] = 1u << x;                                     /* y -> x */
                    for (int v = 0; v < n; v++) o2[v] |= 1u << y;        /* every old vertex beats y */
                    ll H2 = hpaths(n + 2, o2, dp);
                    cases++; if (H2 != 3 * h + 2 * (pa - h)) bad++;
                }
            }
            free(dp); free(Hs); free(Sa);
        }
        printf("any-mark C2 law, order %d: %ld (tournament, mark) cases, %ld failures\n", n, cases, bad);
    }
    /* ---------- Part 2 ---------- */
    char *hit = calloc((size_t)BOUND + 1, 1);
    for (int t = 1; t <= TMAX; t++) {
        int n = t + 2;                    /* r = 0, w = 1, z_i = 1 + i */
        long newhits = 0, checked = 0, bad = 0;
        ll *dp = malloc(sizeof(ll) * ((size_t)1 << (n + 1)) * (n + 1)), *Hs = malloc(sizeof(ll) << (n + 1)), *Sa = malloc(sizeof(ll) << (n + 1));
        for (uint32_t b = 0; b < (1u << t); b++) {
            uint32_t C[MAXV] = {0};
            C[0] |= 1u << 1;
            for (int i = 1; i <= t; i++) { C[1] |= 1u << (1 + i); for (int j = i + 1; j <= t; j++) C[1 + i] |= 1u << (1 + j);
                if ((b >> (i - 1)) & 1) C[0] |= 1u << (1 + i); else C[1 + i] |= 1u << 0; }
            /* B = converse of C, plus sole exit a = n with a -> r (= 0), every other vertex beats a */
            uint32_t B[MAXV] = {0};
            for (int u = 0; u < n; u++) for (int v = 0; v < n; v++) if (u != v && !((C[u] >> v) & 1)) B[u] |= 1u << v;
            B[n] = 1u << 0; for (int v = 1; v < n; v++) B[v] |= 1u << n;
            ll H = hpaths(n + 1, B, dp);
            if (t <= 9) {   /* direct check: the digit of a is 1 */
                ll h, pa; sigma_parts(n + 1, B, n, dp, Hs, Sa, &h, &pa);
                if (pa - h != 1) bad++;
                checked++;
            }
            if (H <= BOUND && !hit[H]) { hit[H] = 1; newhits++; }
        }
        long odd_missing = 0; ll first_miss = -1, last_hit = 0;
        for (ll v = 1; v <= BOUND; v += 2) { if (hit[v]) last_hit = v; }
        for (ll v = 1; v <= BOUND; v += 2) if (!hit[v]) { odd_missing++; if (first_miss < 0) first_miss = v; }
        printf("t=%2d (order %d, digit-1-ready order %d): new values %ld; digit check failures %ld/%ld | odd values <= %lld not yet attained: %ld (first %lld); largest attained %lld\n",
               t, n, n + 1, newhits, bad, checked, BOUND, odd_missing, first_miss, last_hit);
        fflush(stdout);
        free(dp); free(Hs); free(Sa);
    }
    /* list the missing odd values below a small window */
    printf("odd values <= 400 never attained by the family up to t = %d:", TMAX);
    for (ll v = 1; v <= 400 && v <= BOUND; v += 2) if (!hit[v]) printf(" %lld", v);
    printf("\n");
    return 0;
}
