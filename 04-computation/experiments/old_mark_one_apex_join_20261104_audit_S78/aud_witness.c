/*
 * aud_witness.c  (S78 audit, independent code)
 * (1) One-apex triple: for every binary word w of length 0..LMAX (leading zeros included), build T(w) (chain 1->...->n,
 *     apex x, x -> i iff w_i = 1, w_1 most significant) and compare (h, e_x, R_x), counted directly by subset DP, with
 *     (g(v), [v even], 1 + v) where g is the auditor's closed formula (aud_cover.c, g_direct).
 * (2) Explicit witnesses: for target odd N, find a representation in F1 (two one-apex tournaments) or F2 (one-apex with
 *     an order-<=7 rooted tournament, the rooted tournament found by enumeration), BUILD the actual join tournament
 *     (J = A => B with b -> a reversed, B = T(w2)^op or S^op) and count its Hamiltonian paths by direct subset DP.
 * Usage: aud_witness LMAX NTARGETS SEED MAXORDER [LGMIN LGMAX]   (targets N in [2^lg, 2^(lg+1)), lg uniform in [LGMIN, LGMAX])
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>
typedef long long ll;
typedef unsigned long long ull;

static ll g_direct(ull v) {
    if (v == 0) return 1;
    int len = 64 - __builtin_clzll(v);
    ll s = 1 + ((v & 1) == 0);
    for (int i = 0; i + 1 < len; i++) if (!((v >> (i + 1)) & 1) && ((v >> i) & 1)) s++;
    for (int i = 0; i < len; i++) if (!((v >> i) & 1)) s += (ll)(v >> (i + 1));
    return s;
}
static int bitlen(ull v) { return v ? 64 - __builtin_clzll(v) : 0; }

/* direct Hamiltonian path count by subset DP (pull form), counts fit in uint32 for H < 2^32 by monotonicity of h */
static ll hpaths(int n, const uint32_t *out) {
    uint32_t inn[32]; for (int v = 0; v < n; v++) { inn[v] = 0; for (int u = 0; u < n; u++) if ((out[u] >> v) & 1) inn[v] |= 1u << u; }
    size_t full = ((size_t)1 << n) - 1;
    uint32_t *dp = calloc((full + 1) * (size_t)n, sizeof(uint32_t));
    if (!dp) { fprintf(stderr, "alloc %d\n", n); exit(1); }
    for (int v = 0; v < n; v++) dp[((size_t)1 << v) * n + v] = 1;
    /* process by popcount layers so each layer can be parallel */
    for (int k = 2; k <= n; k++) {
        #pragma omp parallel for schedule(dynamic, 4096)
        for (ll S = 1; S <= (ll)full; S++) {
            if (__builtin_popcountll(S) != k) continue;
            for (int v = 0; v < n; v++) if ((S >> v) & 1) {
                ull P = (ull)S & ~(1ULL << v); uint32_t m = inn[v] & (uint32_t)P; ull acc = 0;
                while (m) { int u = __builtin_ctz(m); m &= m - 1; acc += dp[P * n + u]; }
                dp[(size_t)S * n + v] = (uint32_t)acc;
            }
        }
    }
    ull t = 0; for (int v = 0; v < n; v++) t += dp[full * n + v];
    free(dp); return (ll)t;
}
/* small-n full triple at root r (h, e_r, R_r) */
static void triple(int n, const uint32_t *out, int r, ll *h, ll *e, ll *R) {
    size_t full = ((size_t)1 << n) - 1;
    ll *end = calloc((full + 1) * n, sizeof(ll)), *st = calloc((full + 1) * n, sizeof(ll)), *hs = calloc(full + 1, sizeof(ll));
    uint32_t inn[32]; for (int v = 0; v < n; v++) { inn[v] = 0; for (int u = 0; u < n; u++) if ((out[u] >> v) & 1) inn[v] |= 1u << u; }
    hs[0] = 1;
    for (size_t S = 1; S <= full; S++) {
        ll t = 0;
        for (int v = 0; v < n; v++) if ((S >> v) & 1) {
            size_t P = S & ~((size_t)1 << v);
            if (!P) { end[S * n + v] = 1; st[S * n + v] = 1; t++; continue; }
            ll a = 0, b = 0;
            for (int u = 0; u < n; u++) if ((P >> u) & 1) { if ((out[u] >> v) & 1) a += end[P * n + u]; if ((out[v] >> u) & 1) b += st[P * n + u]; }
            end[S * n + v] = a; st[S * n + v] = b; t += a;
        }
        hs[S] = t;
    }
    *h = hs[full]; *e = end[full * n + r]; ll s = 0;
    for (size_t X = 0; X <= full; X++) if (!((X >> r) & 1)) s += hs[X] * st[(full & ~X) * n + r];
    *R = s; free(end); free(st); free(hs);
}
/* one-apex T(w): chain vertices base..base+L-1 (i -> j for i < j), apex base+L; x -> chain i iff bit (L-1-i) of w */
static void oneapex(uint32_t *out, int base, int L, ull w) {
    for (int i = 0; i < L; i++) for (int j = i + 1; j < L; j++) out[base + i] |= 1u << (base + j);
    int x = base + L;
    for (int i = 0; i < L; i++) { if ((w >> (L - 1 - i)) & 1) out[x] |= 1u << (base + i); else out[base + i] |= 1u << x; }
}
/* J = A => B^op with b -> a reversed; A on vertices 0..na-1 (root a), B given (root b) on nb vertices; result n = na+nb */
static int join_op(const uint32_t *A, int na, int a, const uint32_t *B, int nb, int b, uint32_t *J) {
    int n = na + nb; memset(J, 0, sizeof(uint32_t) * 32);
    for (int i = 0; i < na; i++) J[i] = A[i];
    for (int i = 0; i < nb; i++) for (int j = 0; j < nb; j++) if ((B[i] >> j) & 1) J[na + j] |= 1u << (na + i);   /* B^op */
    for (int i = 0; i < na; i++) for (int j = 0; j < nb; j++) { if (i == a && j == b) J[na + j] |= 1u << i; else J[i] |= 1u << (na + j); }
    return n;
}
/* bank with witnesses: triple -> (order, code) of a labelled tournament rooted at 0 */
typedef struct { ll h, e, R; int n; long code; } BT;
static BT *bank; static int nbank;
static void tour_from_code(int n, long code, uint32_t *out) {
    memset(out, 0, sizeof(uint32_t) * 32); int k = 0;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++, k++) { if ((code >> k) & 1) out[i] |= 1u << j; else out[j] |= 1u << i; }
}
static void build_bank(int NMAX) {
    int cap = 4096; bank = malloc(sizeof(BT) * cap); nbank = 0;
    for (int n = 1; n <= NMAX; n++) {
        int na = n * (n - 1) / 2;
        for (long code = 0; code < (1L << na); code++) {
            uint32_t out[32]; tour_from_code(n, code, out); ll h, e, R; triple(n, out, 0, &h, &e, &R);
            int found = 0; for (int i = 0; i < nbank; i++) if (bank[i].h == h && bank[i].e == e && bank[i].R == R) { found = 1; break; }
            if (!found) { if (nbank == cap) { cap *= 2; bank = realloc(bank, sizeof(BT) * cap); } bank[nbank].h = h; bank[nbank].e = e; bank[nbank].R = R; bank[nbank].n = n; bank[nbank].code = code; nbank++; }
        }
    }
}
static ull rs = 0x9E3779B97F4A7C15ULL;
static ull rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }

int main(int argc, char **argv) {
    int LMAX = argc > 1 ? atoi(argv[1]) : 11, NT = argc > 2 ? atoi(argv[2]) : 40; if (argc > 3) rs ^= (ull)atoll(argv[3]) * 0x2545F4914F6CDD1DULL;
    int MAXORD = argc > 4 ? atoi(argv[4]) : 22, LGMIN = argc > 5 ? atoi(argv[5]) : 12, LGMAX = argc > 6 ? atoi(argv[6]) : 22;
    /* (1) one-apex triples */
    long words = 0, bad = 0;
    for (int L = 0; L <= LMAX; L++) for (ull w = 0; w < (1ULL << L); w++) {
        uint32_t out[32]; memset(out, 0, sizeof out); oneapex(out, 0, L, w);
        ll h, e, R; triple(L + 1, out, L, &h, &e, &R); words++;
        if (h != g_direct(w) || e != ((w & 1) == 0) || R != (ll)(1 + w)) { bad++; if (bad < 5) printf("BAD word L=%d w=%llu: (%lld,%lld,%lld)\n", L, w, h, e, R); }
    }
    printf("(1) one-apex triple (g(v), [v even], 1+v) vs direct counts: %ld words of length 0..%d, %ld failures\n", words, LMAX, bad);
    fflush(stdout);
    /* (2) witnesses */
    build_bank(6);   /* witnesses of order <= 6 suffice for most targets; order-7 triples are searched in aud_cover only */
    printf("(2) bank of order <= 6 with explicit witnesses: %d triples\n", nbank);
    long ok = 0, fail = 0, nf1 = 0, nf2 = 0;
    for (int t = 0; t < NT; t++) {
        /* target: odd N with log2 in [12, 22] */
        int lg = LGMIN + (int)(rnd() % (LGMAX - LGMIN + 1)); ll N = (ll)((1ULL << lg) + (rnd() % (1ULL << lg))) | 1;
        /* F1: minimise bitlen(u) + bitlen(v) */
        ll bu = -1, bv = -1; int best = 1 << 30;
        for (ll u = 0; (1 + u) * (1 + u) <= N; u++) {
            ll gu = g_direct(u), eu = ((u & 1) == 0);
            for (ll v = u; (1 + u) * (1 + v) <= N; v++) {
                ll H = (1 + u) * (1 + v) + gu * g_direct(v) - eu * ((v & 1) == 0);
                if (H == N) { int c = bitlen(u) + bitlen(v); if (c < best) { best = c; bu = u; bv = v; } }
            }
        }
        int fam = 0; uint32_t J[32]; int n = 0; char desc[200];
        if (bu >= 0 && best + 2 <= MAXORD) {
            uint32_t A[32], B[32]; memset(A, 0, sizeof A); memset(B, 0, sizeof B);
            int L1 = bitlen(bv), L2 = bitlen(bu);
            oneapex(A, 0, L1, bv); oneapex(B, 0, L2, bu);
            n = join_op(A, L1 + 1, L1, B, L2 + 1, L2, J); fam = 1;
            snprintf(desc, sizeof desc, "F1 v=%lld u=%lld", bv, bu);
        } else {
            /* F2 with the order-<=6 witness bank */
            int bi = -1; ll bvv = -1; int bo = 1 << 30;
            for (int i = 0; i < nbank; i++) for (ll v = 0; (1 + v) * bank[i].R <= N; v++) {
                ll H = g_direct(v) * bank[i].h - ((v & 1) == 0 ? bank[i].e : 0) + (1 + v) * bank[i].R;
                if (H == N) { int c = bitlen(v) + 1 + bank[i].n; if (c < bo) { bo = c; bi = i; bvv = v; } }
            }
            if (bi >= 0 && bo <= MAXORD) {
                uint32_t A[32], S[32]; memset(A, 0, sizeof A);
                int L = bitlen(bvv); oneapex(A, 0, L, bvv); tour_from_code(bank[bi].n, bank[bi].code, S);
                n = join_op(A, L + 1, L, S, bank[bi].n, 0, J); fam = 2;
                snprintf(desc, sizeof desc, "F2 v=%lld y=(%lld,%lld,%lld) order(S)=%d", bvv, bank[bi].h, bank[bi].e, bank[bi].R, bank[bi].n);
            }
        }
        if (!fam) { printf("  N=%lld: no witness of order <= %d in F1 or F2(order<=6 bank) [skipped]\n", N, MAXORD); continue; }
        double t0 = omp_get_wtime();
        ll H = hpaths(n, J);
        printf("  N=%lld: %s, order %d: direct H = %lld %s (%.1fs)\n", N, desc, n, H, H == N ? "OK" : "MISMATCH", omp_get_wtime() - t0);
        fflush(stdout);
        if (H == N) ok++; else fail++;
        if (fam == 1) nf1++; else nf2++;
    }
    printf("(2) witnesses checked %ld (F1 %ld, F2 %ld): %ld OK, %ld mismatches\n", ok + fail, nf1, nf2, ok, fail);
    return 0;
}
