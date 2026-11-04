/*
 * audit_joins.c -- independent audit of the one-arc join law (THM-4866 (B)) and the join coverage (D).
 * Uses the auditor's own census code (audit_census.c, compiled with NO_MAIN).
 * Rooted triples are computed DIRECTLY from their definitions (no converse trick):
 *   x(A,a) = (h(A), R_a(A), e_a(A)),  R_a(A) = sum_{X subset V(A)-a} p(X) s_a(A - X)
 *   y(B,b) = (h(B), R'_b(B), s_b(B)), R'_b(B) = sum_{Y subset V(B)-b} e_b(B - Y) p(Y)
 * J(A,a;B,b) = A => B with a->b reversed. Claim: H(J) = hA hB + R_a R'_b - e_a s_b, hc(J) = s_a(A) e_b(B) (|A|+|B|>=3).
 */
#define NO_MAIN
#include "audit_census.c"

typedef struct { u64 h, R, e; } Tr;

/* all rooted data of tournament t (order n <= 9): p(S) for all S; for each vertex v: s_v(T), e_v(T), R_v, R'_v */
static void rooted_data(int n, const uint16_t *t, u64 *h, u64 *s, u64 *e, u64 *R, u64 *Rp) {
    int full = (1 << n) - 1;
    static __thread u64 dpE[(1 << 9) * 9], dpS[(1 << 9) * 9], p[1 << 9];
    /* dpE[S][v]: paths covering S ending at v (any start) */
    memset(dpE, 0, sizeof(u64) * (size_t)(full + 1) * n);
    for (int v = 0; v < n; v++) dpE[(1 << v) * n + v] = 1;
    for (int S = 1; S <= full; S++) for (int v = 0; v < n; v++) {
        u64 c = dpE[S * n + v]; if (!c) continue;
        for (int w = 0; w < n; w++) if (!((S >> w) & 1) && ((t[v] >> w) & 1)) dpE[(S | (1 << w)) * n + w] += c;
    }
    p[0] = 1;
    for (int S = 1; S <= full; S++) { u64 x = 0; for (int v = 0; v < n; v++) x += dpE[S * n + v]; p[S] = x; }
    *h = p[full];
    for (int a = 0; a < n; a++) {
        /* dpS[S][v]: paths starting at a covering S ending at v */
        memset(dpS, 0, sizeof(u64) * (size_t)(full + 1) * n);
        dpS[(1 << a) * n + a] = 1;
        for (int S = 1; S <= full; S++) {
            if (!((S >> a) & 1)) continue;
            for (int v = 0; v < n; v++) {
                u64 c = dpS[S * n + v]; if (!c) continue;
                for (int w = 0; w < n; w++) if (!((S >> w) & 1) && ((t[v] >> w) & 1)) dpS[(S | (1 << w)) * n + w] += c;
            }
        }
        u64 rr = 0, rp = 0;
        int rest = full & ~(1 << a);
        /* iterate X over subsets of rest (including empty) */
        for (int X = rest;; X = (X - 1) & rest) {
            int Y = full & ~X;       /* contains a */
            u64 sa = 0; for (int v = 0; v < n; v++) sa += dpS[Y * n + v];
            rr += p[X] * sa;                 /* p(X) s_a(T - X) */
            rp += dpE[Y * n + a] * p[X];     /* e_a(T - X) p(X) */
            if (X == 0) break;
        }
        u64 sfull = 0; for (int v = 0; v < n; v++) sfull += dpS[full * n + v];
        s[a] = sfull; e[a] = dpE[full * n + a]; R[a] = rr; Rp[a] = rp;
    }
}

static void build_join(int na, const uint16_t *A, int a, int nb, const uint16_t *B, int b, uint16_t *J) {
    int n = na + nb;
    for (int i = 0; i < n; i++) J[i] = 0;
    for (int i = 0; i < na; i++) { J[i] = A[i]; for (int j = 0; j < nb; j++) J[i] |= 1u << (na + j); }
    for (int j = 0; j < nb; j++) J[na + j] = (uint16_t)(B[j] << na);
    J[a] &= ~(1u << (na + b)); J[na + b] |= 1u << a;
}

static int test_join(int na, const uint16_t *A, int a, int nb, const uint16_t *B, int b, int *hcbad) {
    u64 hA, sA[9], eA[9], RA[9], RpA[9], hB, sB[9], eB[9], RB[9], RpB[9];
    rooted_data(na, A, &hA, sA, eA, RA, RpA);
    rooted_data(nb, B, &hB, sB, eB, RB, RpB);
    uint16_t J[MAXN]; build_join(na, A, a, nb, B, b, J);
    int n = na + nb;
    u64 H = count_paths(n, J), hc = count_cycles(n, J);
    u64 pred = hA * hB + RA[a] * RpB[b] - eA[a] * sB[b];
    u64 predc = (na + nb >= 3) ? sA[a] * eB[b] : 0;
    *hcbad = (hc != predc);
    return H != pred;
}

static int cmp_tr(const void *x, const void *y) {
    const Tr *a = x, *b = y;
    if (a->h != b->h) return a->h < b->h ? -1 : 1;
    if (a->R != b->R) return a->R < b->R ? -1 : 1;
    if (a->e != b->e) return a->e < b->e ? -1 : 1;
    return 0;
}
static size_t dedupe(Tr *v, size_t m) { qsort(v, m, sizeof(Tr), cmp_tr); size_t k = 0; for (size_t i = 0; i < m; i++) if (!k || cmp_tr(&v[i], &v[k - 1])) v[k++] = v[i]; return k; }

static void coverage(int K, long VMAX, Tr *X, size_t nx) {
    long *r = calloc(VMAX + 1, sizeof(long));
    for (size_t i = 0; i < nx; i++) for (size_t j = 0; j < nx; j++) {
        /* y-set equals x-set (checked); y = (h, R', s) is stored as a Tr with fields (h, R, e) = (h, R', s) */
        long long N = (long long)(X[i].h * X[j].h + X[i].R * X[j].R) - (long long)(X[i].e * X[j].e);
        if (N >= 1 && N <= VMAX) r[N]++;
    }
    printf("K=%d VMAX=%ld: %zu distinct triples; odd N with r(N)=0:", K, VMAX, nx);
    int shown = 0; long first_above21 = -1;
    for (long v = 1; v <= VMAX; v += 2) if (!r[v]) { if (shown < 12) { printf(" %ld", v); shown++; } if (v > 21 && first_above21 < 0) first_above21 = v; }
    printf(" ...; first gap above 21: %ld\n", first_above21);
    long ev = 0; for (long v = 2; v <= VMAX; v += 2) ev += r[v];
    printf("   even N hit: %ld\n", ev);
    for (long lo = 128; lo * 2 <= VMAX; lo *= 2) {
        long mn = -1; for (long v = lo + 1; v < 2 * lo; v += 2) if (mn < 0 || r[v] < mn) mn = r[v];
        printf("   window [%ld,%ld) min r = %ld\n", lo, 2 * lo, mn);
    }
    long lo = VMAX / 8, hi = VMAX / 4;
    int mods[3] = {3, 5, 7};
    for (int q = 0; q < 3; q++) {
        int md = mods[q]; double s[8] = {0}; long c[8] = {0};
        for (long v = lo | 1; v < hi; v += 2) { s[v % md] += r[v]; c[v % md]++; }
        double tot = 0; long tc = 0; for (int k = 0; k < md; k++) { tot += s[k]; tc += c[k]; }
        printf("   [%ld,%ld) mean r by N mod %d (relative):", lo, hi, md);
        for (int k = 0; k < md; k++) printf(" %d:%.3f", k, (s[k] / c[k]) / (tot / tc));
        if (md == 3) {
            double m0 = s[0] / c[0], m12 = (s[1] + s[2]) / (c[1] + c[2]);
            printf("  ratio res0/others = %.3f", m0 / m12);
        }
        printf("\n");
    }
    /* triple residues mod 3 */
    long hist[3][3][3] = {{{0}}};
    for (size_t i = 0; i < nx; i++) hist[X[i].h % 3][X[i].R % 3][X[i].e % 3]++;
    long h0 = 0; for (int b = 0; b < 3; b++) for (int c2 = 0; c2 < 3; c2++) h0 += hist[0][b][c2];
    long R0 = 0; for (int a2 = 0; a2 < 3; a2++) for (int c2 = 0; c2 < 3; c2++) R0 += hist[a2][0][c2];
    printf("   triples with h = 0 mod 3: %ld/%zu; R = 0 mod 3: %ld/%zu\n", h0, nx, R0, nx);
    /* local factors: over ALL pairs (no window), share of Q = 0 mod p versus Q = c != 0, against the generic
       bilinear-form value 1 + p/(p^3 - 1) for triples equidistributed mod p */
    int ps[3] = {3, 5, 7};
    for (int q = 0; q < 3; q++) {
        int p = ps[q]; double cnt[8] = {0}; long z = 0;
        for (size_t i = 0; i < nx; i++) if (X[i].h % p == 0 && X[i].R % p == 0 && X[i].e % p == 0) z++;
        for (size_t i = 0; i < nx; i++) for (size_t j = 0; j < nx; j++) {
            long long N = (long long)(X[i].h * X[j].h + X[i].R * X[j].R) - (long long)(X[i].e * X[j].e);
            cnt[N % p] += 1;
        }
        double nz = 0; for (int c = 1; c < p; c++) nz += cnt[c]; nz /= (p - 1);
        printf("   mod %d over all pairs: P(Q=0)/P(Q=c!=0) = %.3f (generic 1 + p/(p^3-1) = %.3f); triples = 0 mod %d: %ld/%zu\n",
               p, cnt[0] / nz, 1.0 + (double)p / (p * p * p - 1), p, z, nx);
    }
    free(r);
}

int main(int argc, char **argv) {
    int K = 7;
    for (int n = 1; n <= 9; n++) build(n);
    printf("classes per order 1..9:"); for (int n = 1; n <= 9; n++) printf(" %d", ncls[n]); printf("\n");
    /* 1. exhaustive join test, all rooted class pairs with |A| + |B| <= 9 */
    {
        long tests = 0, bad = 0, hcbad = 0, tiny = 0;
        for (int na = 1; na <= 8; na++) for (int nb = 1; na + nb <= 9; nb++) {
            long lt = 0, lb = 0, lhb = 0;
#pragma omp parallel for schedule(dynamic, 4) reduction(+:lt, lb, lhb)
            for (int ca = 0; ca < ncls[na]; ca++) {
                uint16_t A[MAXN]; decode(na, cls[na][ca], A);
                for (int cb = 0; cb < ncls[nb]; cb++) {
                    uint16_t B[MAXN]; decode(nb, cls[nb][cb], B);
                    for (int a = 0; a < na; a++) for (int b = 0; b < nb; b++) {
                        int hb; lb += test_join(na, A, a, nb, B, b, &hb); lhb += hb; lt++;
                    }
                }
            }
            tests += lt; bad += lb; hcbad += lhb;
            if (na == 1 && nb == 1) tiny = lhb;
        }
        printf("exhaustive joins |A|+|B| <= 9 over rooted class pairs: %ld tests, H mismatches %ld, hc mismatches %ld (|A|=|B|=1 case uses hc=0 by convention: %ld)\n", tests, bad, hcbad, tiny);
        /* the |A|=|B|=1 exception, explicitly */
        uint16_t A1[1] = {0}, B1[1] = {0}, J[2]; build_join(1, A1, 0, 1, B1, 0, J);
        printf("   |A|=|B|=1: H(J) = %llu, hc(J) = %llu, s_a e_b = 1 (so hc = s_a e_b fails here, as the exception states)\n", count_paths(2, J), count_cycles(2, J));
    }
    /* 2. random labelled joins, sides of order <= 6 (J of order <= 12 = MAXN) */
    {
        long tests = 0, bad = 0, hcbad = 0;
#pragma omp parallel for schedule(dynamic, 16) reduction(+:tests, bad, hcbad)
        for (int it = 0; it < 4000; it++) {
            unsigned long long s = 0x9E3779B97F4A7C15ULL * (it + 12345);
            #define RND() (s ^= s << 13, s ^= s >> 7, s ^= s << 17, s)
            int na = 1 + (int)(RND() % 6), nb = 1 + (int)(RND() % 6);
            uint16_t A[MAXN] = {0}, B[MAXN] = {0};
            for (int i = 0; i < na; i++) for (int j = i + 1; j < na; j++) { if (RND() & 1) A[i] |= 1u << j; else A[j] |= 1u << i; }
            for (int i = 0; i < nb; i++) for (int j = i + 1; j < nb; j++) { if (RND() & 1) B[i] |= 1u << j; else B[j] |= 1u << i; }
            int a = (int)(RND() % na), b = (int)(RND() % nb);
            int hb; bad += test_join(na, A, a, nb, B, b, &hb); hcbad += hb; tests++;
        }
        printf("random labelled joins (sides 1..6, independent seeds): %ld tests, H mismatches %ld, hc mismatches %ld\n", tests, bad, hcbad);
    }
    /* 3. triple sets and coverage */
    for (int KK = 5; KK <= K; KK++) {
        size_t cap = 0; for (int n = 1; n <= KK; n++) cap += (size_t)ncls[n] * n;
        Tr *X = malloc(sizeof(Tr) * cap), *Y = malloc(sizeof(Tr) * cap); size_t m = 0;
        for (int n = 1; n <= KK; n++) for (int c = 0; c < ncls[n]; c++) {
            uint16_t t[MAXN]; decode(n, cls[n][c], t);
            u64 h, s[9], e[9], R[9], Rp[9]; rooted_data(n, t, &h, s, e, R, Rp);
            for (int v = 0; v < n; v++) { X[m].h = h; X[m].R = R[v]; X[m].e = e[v]; Y[m].h = h; Y[m].R = Rp[v]; Y[m].e = s[v]; m++; }
        }
        size_t nx = dedupe(X, m), ny = dedupe(Y, m);
        int same = (nx == ny) && !memcmp(X, Y, sizeof(Tr) * nx);
        printf("orders <= %d: %zu x-triples, %zu y-triples, sets equal: %s\n", KK, nx, ny, same ? "yes" : "NO");
        long VMAX = KK == 5 ? 3000 : KK == 6 ? 40000 : 200000;
        coverage(KK, VMAX, X, nx);
        free(X); free(Y);
    }
    return 0;
}
