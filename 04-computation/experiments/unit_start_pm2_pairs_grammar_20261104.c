/*
 * unit_start_pm2_pairs_grammar_20261104.c  (opus-2026-11-04-S76)
 *
 * Generates every unique-start rooted strong tournament (s_T(r) = 1) by the unique decomposition of FUNNELS F1/F2
 * (05-knowledge/results/polarized_path_calculus_20261103/FUNNELS.md):
 *   T - r = U => v_1 -> ... -> v_t (transitive tail dominated by U), r beats exactly one vertex w of the first strong
 *   component U, s_U(w) = 1 recursively; if t > 0 then v_t -> r and the arcs between r and v_1..v_(t-1) are free;
 *   U a singleton needs t = n - 2 >= 1; U a smaller core of order m >= 3 allows t >= 0.
 * F2 proves each rooted class arises exactly once, giving u_n = F(2n - 4).
 *
 * For each class: the profile (h, s, q, F) of adaptive_metagraph_reformulation_20261103/ADAPTIVE.md (s is re-checked = 1),
 * then the unit-start +-2 port pairs: F = F', h' = h + 2.  Session189's dyadic wrappers H(D_m) = h + 2^m F then differ by
 * exactly 2 at every clock m.
 *
 * Usage: unitgram NMAX   (NMAX <= 15)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>

#define MAXN 16
typedef long long ll;
typedef struct { int n; uint32_t out[MAXN]; } RT;       /* root is vertex 0, its unique entry is vertex 1 */
typedef struct { ll h, s, q, F; } Prof;

static RT *L[MAXN + 1]; static size_t Lcount[MAXN + 1], Lcap[MAXN + 1];
static void push(int n, const RT *t) {
    if (Lcount[n] == Lcap[n]) { Lcap[n] = Lcap[n] ? 2 * Lcap[n] : 1024; L[n] = realloc(L[n], sizeof(RT) * Lcap[n]); }
    L[n][Lcount[n]++] = *t;
}

/* assemble: root 0, core U (rooted at its own root, which becomes vertex 1 = w), tail of length t with free bits */
static void assemble(int n, const RT *U, int t, uint32_t freebits, RT *res) {
    int m = U ? U->n : 1;                 /* U occupies vertices 1..m, its root (unique entry w) at vertex 1 */
    memset(res, 0, sizeof *res); res->n = n;
    if (U) for (int a = 0; a < m; a++) for (int b = 0; b < m; b++) if ((U->out[a] >> b) & 1) res->out[1 + a] |= 1u << (1 + b);
    /* r = 0 beats exactly w = 1 inside U; the rest of U beats r */
    res->out[0] |= 1u << 1;
    for (int a = 2; a <= m; a++) res->out[a] |= 1u << 0;
    /* tail vertices m+1 .. m+t, transitive, dominated by U */
    for (int i = 0; i < t; i++) {
        int vi = m + 1 + i;
        for (int a = 1; a <= m; a++) res->out[a] |= 1u << vi;
        for (int j = i + 1; j < t; j++) res->out[vi] |= 1u << (m + 1 + j);
        if (i == t - 1) res->out[vi] |= 1u << 0;                 /* v_t -> r */
        else if ((freebits >> i) & 1) res->out[0] |= 1u << vi;  /* r -> v_i */
        else res->out[vi] |= 1u << 0;                           /* v_i -> r */
    }
}

static Prof profile(const RT *T, ll *dpe, ll *Hs) {
    int n = T->n, full = (1 << n) - 1, v = 0;
    for (size_t i = 0; i < (size_t)(full + 1) * n; i++) dpe[i] = 0;
    for (int x = 0; x < n; x++) dpe[(size_t)(1 << x) * n + x] = 1;
    for (int S = 1; S <= full; S++) for (int x = 0; x < n; x++) { ll c = dpe[(size_t)S * n + x]; if (!c) continue;
        uint32_t nb = T->out[x] & ~(uint32_t)S; while (nb) { int y = __builtin_ctz(nb); nb &= nb - 1; dpe[(size_t)(S | (1 << y)) * n + y] += c; } }
    Hs[0] = 1; for (int S = 1; S <= full; S++) { ll t = 0; for (int x = 0; x < n; x++) t += dpe[(size_t)S * n + x]; Hs[S] = t; }
    Prof p; p.h = Hs[full]; p.q = Hs[full & ~1];
    ll F = 0; for (int S = 1; S <= full; S++) if (S & 1) F += dpe[(size_t)S * n + v] * Hs[full & ~S];
    p.F = F;
    /* s via the converse */
    uint32_t rev[MAXN]; for (int a = 0; a < n; a++) { rev[a] = 0; for (int b = 0; b < n; b++) if (b != a && !((T->out[a] >> b) & 1)) rev[a] |= 1u << b; }
    for (size_t i = 0; i < (size_t)(full + 1) * n; i++) dpe[i] = 0;
    for (int x = 0; x < n; x++) dpe[(size_t)(1 << x) * n + x] = 1;
    for (int S = 1; S <= full; S++) for (int x = 0; x < n; x++) { ll c = dpe[(size_t)S * n + x]; if (!c) continue;
        uint32_t nb = rev[x] & ~(uint32_t)S; while (nb) { int y = __builtin_ctz(nb); nb &= nb - 1; dpe[(size_t)(S | (1 << y)) * n + y] += c; } }
    p.s = dpe[(size_t)full * n + v];
    return p;
}

typedef struct { ll F, h, q; size_t idx; } Key;
static int cmpk(const void *a, const void *b) { const Key *x = a, *y = b;
    if (x->F != y->F) return x->F < y->F ? -1 : 1; if (x->h != y->h) return x->h < y->h ? -1 : 1; return 0; }

static ll fib(int k) { ll a = 0, b = 1; for (int i = 0; i < k; i++) { ll t = a + b; a = b; b = t; } return a; }

int main(int argc, char **argv) {
    int NMAX = argc > 1 ? atoi(argv[1]) : 12;
    printf("# unit_start_pm2_pairs_grammar_20261104 NMAX=%d\n", NMAX);
    for (int n = 3; n <= NMAX; n++) {
        RT t;
        /* U singleton (vertex 1 alone), t = n - 2 */
        { int tl = n - 2; for (uint32_t fb = 0; fb < (1u << (tl - 1)); fb++) { assemble(n, NULL, tl, fb, &t); push(n, &t); } }
        /* U a smaller core of order m >= 3 */
        for (int m = 3; m <= n - 1; m++) {
            int tl = n - 1 - m;
            for (size_t i = 0; i < Lcount[m]; i++) {
                if (tl == 0) { assemble(n, &L[m][i], 0, 0, &t); push(n, &t); }
                else for (uint32_t fb = 0; fb < (1u << (tl - 1)); fb++) { assemble(n, &L[m][i], tl, fb, &t); push(n, &t); }
            }
        }
        size_t cnt = Lcount[n];
        Prof *P = malloc(sizeof(Prof) * cnt);
        long bad_s = 0;
        #pragma omp parallel
        {
            ll *dpe = malloc(sizeof(ll) * ((size_t)1 << n) * n), *Hs = malloc(sizeof(ll) * ((size_t)1 << n));
            #pragma omp for schedule(dynamic, 16) reduction(+:bad_s)
            for (size_t i = 0; i < cnt; i++) { P[i] = profile(&L[n][i], dpe, Hs); if (P[i].s != 1) bad_s++; }
            free(dpe); free(Hs);
        }
        Key *K = malloc(sizeof(Key) * cnt);
        for (size_t i = 0; i < cnt; i++) { K[i].F = P[i].F; K[i].h = P[i].h; K[i].q = P[i].q; K[i].idx = i; }
        qsort(K, cnt, sizeof(Key), cmpk);
        long pairs = 0, pairs_q = 0; ll minh = -1, minF = -1; long distinctF = 0;
        for (size_t a = 0; a < cnt; a++) {
            if (a == 0 || K[a].F != K[a - 1].F) distinctF++;
            for (size_t b = a + 1; b < cnt && K[b].F == K[a].F; b++) if (K[b].h - K[a].h == 2) {
                pairs++; if (K[b].q == K[a].q) pairs_q++;
                if (minh < 0 || K[a].h < minh) { minh = K[a].h; minF = K[a].F; }
            }
        }
        printf("n=%2d unit-start rooted strong classes %zu (F(2n-4) = %lld) s!=1: %ld | distinct F values %ld | +-2 pairs (F=F', h'=h+2) %ld (with q=q' also: %ld)",
               n, cnt, fib(2 * n - 4), bad_s, distinctF, pairs, pairs_q);
        if (pairs) printf("  smallest h pair: (h, h+2) = (%lld, %lld), F = %lld", minh, minh + 2, minF);
        printf("\n"); fflush(stdout);
        free(P); free(K);
    }
    return 0;
}
