/* Independent audit search for HYP-9429: rational n = alpha^4 + beta^4 + gamma^4 + delta^4, alpha.. in Z[w].
 * Written from scratch for the S66 audit; shares no code with the author's program.
 *
 * Elements: (x, y) = x + y w, w^2 + w + 1 = 0.  Box: |a|, |b| <= B for a base a + b w.
 *
 * mode win B lo hi
 *   All targets n in [lo, hi] at once.  Every pair sum P = f_i + f_j (i <= j, f = fourth powers of the box up to sign,
 *   0 included) is stored as (y, x); the array is sorted in the lexicographic order of Z^2 (y first), which is a
 *   translation-invariant total order.  For n in [lo, hi] we need P + Q = (0, n) [as (y, x)], and the set
 *   {(0, n) : lo <= n <= hi} is an interval [(0, lo), (0, hi)] of that order.  A single sweep with i increasing and two
 *   pointers L_i >= R_i ... (both non-increasing in i) lists every (i, j), j >= i, with S_i + S_j in the interval.
 *   (The author's program buckets by y and binary-searches x per element; this is a different algorithm.)
 *
 * mode deep n B [maxwit]
 *   One target n.  Classes c(alpha) = alpha^4 mod 16 (7 values).  All 4-multisets of classes with sum = n mod 16
 *   (as an element of Z[w]/16) are enumerated by machine; for each, a split into two pairs is chosen; the pair sums are
 *   partitioned by their residue mod RM (default 105; in Z[w]/RM) and matched class by class (rho on the left with n - rho on the
 *   right) by hashing + radix sort + merge.  Exhaustive for the box; every hit is verified by recomputation.
 *
 * mode classes n : print the admissible class multisets for n (the residue analysis done by machine).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <omp.h>

#ifndef RM
#define RM 105
#endif
#define RM2 (RM * RM)

typedef long long ll;
typedef unsigned long long ull;

static void p4(ll a, ll b, ll *x, ll *y) {
    /* (a + b w)^2 = (a^2 - b^2) + (2ab - b^2) w ; square again */
    ll s = a * a - b * b, t = 2 * a * b - b * b;
    *x = s * s - t * t;
    *y = 2 * s * t - t * t;
}

static ll mod(ll v, ll m) { ll r = v % m; return r < 0 ? r + m : r; }

static inline ull mix(ull z) {
    z += 0x9e3779b97f4a7c15ULL; z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL; z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}
static inline ull hkey(ll x, ll y) { return mix((ull)x ^ mix((ull)y + 0x1234567ULL)); }

/* ------------------------------------------------------------------ base list */
typedef struct { ll x, y; int a, b; } fp;

static fp *bases(int B, int *cnt) {
    long cap = (long)(2 * B + 1) * (2 * B + 1) / 2 + 2, n = 0;
    fp *F = malloc(sizeof(fp) * cap);
    for (int a = 0; a <= B; a++)
        for (int b = -B; b <= B; b++) {
            if (a == 0 && b < 0) continue;            /* canonical sign: a > 0, or a = 0 and b >= 0 */
            ll x, y; p4(a, b, &x, &y);
            F[n].x = x; F[n].y = y; F[n].a = a; F[n].b = b; n++;
        }
    *cnt = (int)n;
    return F;
}

static int verify4(ll n, const int *a, const int *b) {
    __int128 sx = 0, sy = 0;
    for (int k = 0; k < 4; k++) {
        __int128 A = a[k], Bb = b[k];
        __int128 s = A * A - Bb * Bb, t = 2 * A * Bb - Bb * Bb;
        sx += s * s - t * t; sy += 2 * s * t - t * t;
    }
    return sx == n && sy == 0;
}

/* ------------------------------------------------------------------ mode win */
typedef struct { ll y, x; } yx;

static int cmp_yx(const void *p, const void *q) {
    const yx *u = p, *v = q;
    if (u->y != v->y) return u->y < v->y ? -1 : 1;
    if (u->x != v->x) return u->x < v->x ? -1 : 1;
    return 0;
}

static int lexcmp(ll y, ll x, ll Y, ll X) {      /* compare (y, x) with (Y, X) */
    if (y != Y) return y < Y ? -1 : 1;
    if (x != X) return x < X ? -1 : 1;
    return 0;
}

static int mode_win(int B, ll lo, ll hi) {
    clock_t t0 = clock();
    int nf; fp *F = bases(B, &nf);
    size_t N = (size_t)nf * (nf + 1) / 2;
    ll ymax = 0;
    for (int i = 0; i < nf; i++) { ll v = F[i].y < 0 ? -F[i].y : F[i].y; if (v > ymax) ymax = v; }
    ymax *= 2;
    const int NB = 1 << 16;
    size_t *cnt = calloc(NB + 1, sizeof(size_t));
    /* bucket index monotone in y */
    #define BUCKET(yy) ((int)(((__int128)((yy) + ymax) * NB) / (2 * ymax + 1)))
    for (int i = 0; i < nf; i++) for (int j = i; j < nf; j++) cnt[BUCKET(F[i].y + F[j].y) + 1]++;
    for (int k = 0; k < NB; k++) cnt[k + 1] += cnt[k];
    yx *S = malloc(sizeof(yx) * N);
    if (!S) { fprintf(stderr, "alloc %zu\n", N); return 2; }
    size_t *pos = malloc(sizeof(size_t) * NB);
    memcpy(pos, cnt, sizeof(size_t) * NB);
    for (int i = 0; i < nf; i++) for (int j = i; j < nf; j++) {
        ll y = F[i].y + F[j].y; size_t p = pos[BUCKET(y)]++;
        S[p].y = y; S[p].x = F[i].x + F[j].x;
    }
    #pragma omp parallel for schedule(dynamic, 64)
    for (int k = 0; k < NB; k++) qsort(S + cnt[k], cnt[k + 1] - cnt[k], sizeof(yx), cmp_yx);
    for (size_t p = 1; p < N; p++) if (cmp_yx(S + p - 1, S + p) > 0) { fprintf(stderr, "sort error\n"); return 3; }
    double tsort = (double)(clock() - t0) / CLOCKS_PER_SEC;
    int W = (int)(hi - lo + 1);
    long long *hits = calloc(W, sizeof(long long));
    yx *w1 = calloc(W, sizeof(yx)), *w2 = calloc(W, sizeof(yx));
    /* sweep: for each i, L = first j with S_i + S_j >= (0, lo) (N if none), R = last j with S_i + S_j <= (0, hi) */
    long long L = (long long)N, R = (long long)N - 1;
    for (size_t i = 0; i < N; i++) {
        while (R >= 0 && lexcmp(S[i].y + S[R].y, S[i].x + S[R].x, 0, hi) > 0) R--;
        if (L > R + 1) L = R + 1;
        while (L > 0 && lexcmp(S[i].y + S[L - 1].y, S[i].x + S[L - 1].x, 0, lo) >= 0) L--;
        if (R < (long long)i) break;
        long long js = L > (long long)i ? L : (long long)i;
        for (long long j = js; j <= R; j++) {
            ll sy = S[i].y + S[j].y, sx = S[i].x + S[j].x;
            if (sy != 0 || sx < lo || sx > hi) { fprintf(stderr, "window logic error\n"); return 4; }
            int t = (int)(sx - lo);
            if (hits[t] == 0) { w1[t] = S[i]; w2[t] = S[j]; }
            hits[t]++;
        }
    }
    /* witnesses: decompose each first-hit pair sum with an open-addressing table of the fourth powers */
    size_t cap = 1; while (cap < 4 * (size_t)nf) cap <<= 1;
    int *tab = malloc(sizeof(int) * cap);
    for (size_t k = 0; k < cap; k++) tab[k] = -1;
    for (int i = 0; i < nf; i++) {
        size_t h = (size_t)(hkey(F[i].x, F[i].y) & (cap - 1));
        while (tab[h] >= 0) h = (h + 1) & (cap - 1);
        tab[h] = i;
    }
    int none = 0;
    for (int t = 0; t < W; t++) {
        ll n = lo + t;
        if (!hits[t]) { printf("n = %lld: none\n", n); none++; continue; }
        int A[4], Bv[4], ok = 0;
        yx P[2] = {w1[t], w2[t]};
        for (int s = 0; s < 2; s++) {
            int found = 0;
            for (int i = 0; i < nf && !found; i++) {
                ll ry = P[s].y - F[i].y, rx = P[s].x - F[i].x;
                size_t h = (size_t)(hkey(rx, ry) & (cap - 1));
                while (tab[h] >= 0) {
                    int j = tab[h];
                    if (F[j].x == rx && F[j].y == ry) {
                        A[2 * s] = F[i].a; Bv[2 * s] = F[i].b; A[2 * s + 1] = F[j].a; Bv[2 * s + 1] = F[j].b; found = 1; break;
                    }
                    h = (h + 1) & (cap - 1);
                }
            }
            ok += found;
        }
        if (ok != 2 || !verify4(n, A, Bv)) { printf("n = %lld: WITNESS RECOVERY/VERIFY FAILED\n", n); return 5; }
        printf("n = %lld: %lld pair-pair hits; e.g. (%d,%d) (%d,%d) (%d,%d) (%d,%d)\n", n, hits[t], A[0], Bv[0], A[1], Bv[1],
               A[2], Bv[2], A[3], Bv[3]);
    }
    printf("win B = %d: %d fourth powers (up to sign), %zu pair sums; %d of %d targets in [%lld, %lld] without a "
           "representation [sort %.0f s, total %.0f s cpu]\n", B, nf, N, none, W, lo, hi, tsort,
           (double)(clock() - t0) / CLOCKS_PER_SEC);
    return 0;
}

/* ------------------------------------------------------------------ mode deep */
static const ll CLS[7][2] = {{0, 0}, {1, 0}, {9, 0}, {0, 1}, {0, 9}, {15, 15}, {7, 7}};   /* alpha^4 mod 16 (x, y) */
static const char *CLSNAME[7] = {"0", "1", "9", "w", "9w", "w2", "9w2"};

static int cls16(ll x, ll y) {
    ll X = mod(x, 16), Y = mod(y, 16);
    for (int c = 0; c < 7; c++) if (CLS[c][0] == X && CLS[c][1] == Y) return c;
    return -1;
}

typedef struct { int c[4]; } ms;

static int admissible(ll n, ms *out) {
    int k = 0;
    for (int a = 0; a < 7; a++) for (int b = a; b < 7; b++) for (int c = b; c < 7; c++) for (int d = c; d < 7; d++) {
        ll sx = CLS[a][0] + CLS[b][0] + CLS[c][0] + CLS[d][0], sy = CLS[a][1] + CLS[b][1] + CLS[c][1] + CLS[d][1];
        if (mod(sx - n, 16) == 0 && mod(sy, 16) == 0) { out[k].c[0] = a; out[k].c[1] = b; out[k].c[2] = c; out[k].c[3] = d; k++; }
    }
    return k;
}

typedef struct { ull h; int i, j; } ent;

static void rsort(ent *a, ent *tmp, size_t n) {          /* LSD radix sort on h, 16-bit digits */
    size_t *c = malloc(sizeof(size_t) * 65537);
    for (int pass = 0; pass < 4; pass++) {
        int sh = 16 * pass;
        memset(c, 0, sizeof(size_t) * 65537);
        for (size_t k = 0; k < n; k++) c[((a[k].h >> sh) & 0xffff) + 1]++;
        for (int k = 0; k < 65536; k++) c[k + 1] += c[k];
        for (size_t k = 0; k < n; k++) tmp[c[(a[k].h >> sh) & 0xffff]++] = a[k];
        ent *s = a; a = tmp; tmp = s;
    }
    free(c);   /* 4 passes: result back in the original array */
}

/* per class (c16) and residue r15 (0..224) bucket of fourth powers */
typedef struct { fp *v; int n; int *off; } cb;

static cb CB[7];

static void build_classes(fp *F, int nf) {
    for (int c = 0; c < 7; c++) { CB[c].n = 0; }
    int *cc = calloc(7, sizeof(int));
    int *cls = malloc(sizeof(int) * nf), *r15 = malloc(sizeof(int) * nf);
    for (int i = 0; i < nf; i++) {
        cls[i] = cls16(F[i].x, F[i].y);
        if (cls[i] < 0) { fprintf(stderr, "fourth power outside the 7 classes\n"); exit(7); }
        r15[i] = (int)(mod(F[i].x, RM) * RM + mod(F[i].y, RM));
        cc[cls[i]]++;
    }
    for (int c = 0; c < 7; c++) {
        CB[c].v = malloc(sizeof(fp) * (cc[c] + 1));
        int *cnt = calloc(RM2 + 1, sizeof(int));
        for (int i = 0; i < nf; i++) if (cls[i] == c) cnt[r15[i] + 1]++;
        for (int r = 0; r < RM2; r++) cnt[r + 1] += cnt[r];
        CB[c].off = malloc(sizeof(int) * (RM2 + 1));
        memcpy(CB[c].off, cnt, sizeof(int) * (RM2 + 1));
        int *pos = malloc(sizeof(int) * (RM2 + 1)); memcpy(pos, cnt, sizeof(int) * (RM2 + 1));
        for (int i = 0; i < nf; i++) if (cls[i] == c) CB[c].v[pos[r15[i]]++] = F[i];
        free(cnt); free(pos);
        CB[c].n = cc[c];
    }
    free(cc); free(cls); free(r15);
}

static inline int radd(int r1, int r2) { return (int)(((r1 / RM + r2 / RM) % RM) * RM + (r1 % RM + r2 % RM) % RM); }
static inline int rsub(int r1, int r2) { return (int)(((r1 / RM - r2 / RM + RM) % RM) * RM + (r1 % RM - r2 % RM + RM) % RM); }

/* enumerate pair sums of classes (c1, c2) with total residue rho; if neg, store hash of (n,0) - P; returns count.
   i, j are stored as global indices into CB[c1].v and CB[c2].v. */
static size_t gen(int c1, int c2, int rho, int neg, ll n, ent *out, int count_only) {
    size_t k = 0;
    for (int r1 = 0; r1 < RM2; r1++) {
        int s1 = CB[c1].off[r1], e1 = CB[c1].off[r1 + 1];
        if (s1 == e1) continue;
        int r2 = rsub(rho, r1);
        if (c1 == c2 && r2 < r1) continue;              /* unordered residue pair */
        int s2 = CB[c2].off[r2], e2 = CB[c2].off[r2 + 1];
        if (s2 == e2) continue;
        for (int i = s1; i < e1; i++) {
            int jstart = (c1 == c2 && r1 == r2) ? i : s2;
            if (count_only) { k += (size_t)(e2 - jstart); continue; }
            const fp *A = &CB[c1].v[i];
            for (int j = jstart; j < e2; j++) {
                const fp *Bq = &CB[c2].v[j];
                ll x = A->x + Bq->x, y = A->y + Bq->y;
                if (neg) { x = n - x; y = -y; }
                out[k].h = hkey(x, y); out[k].i = i; out[k].j = j; k++;
            }
        }
    }
    return k;
}

static ll PATTERN_N = -1;   /* negative control: if >= 0, use the class multisets admissible for PATTERN_N */
static int mode_deep(ll n, int B, int maxwit) {
    double t0 = omp_get_wtime();
    int nf; fp *F = bases(B, &nf);
    build_classes(F, nf);
    ms M[256]; int nm = admissible(PATTERN_N >= 0 ? PATTERN_N : n, M);
    printf("deep n = %lld, B = %d: %d fourth powers up to sign; class sizes:", n, B, nf);
    for (int c = 0; c < 7; c++) printf(" %s:%d", CLSNAME[c], CB[c].n);
    printf("\n  admissible class multisets mod 16 (machine-enumerated): %d\n", nm);
    long long tothits = 0; size_t totpairs = 0; int shown = 0;
    for (int m = 0; m < nm; m++) {
        int *c = M[m].c;
        /* choose the split with least work */
        int sp[3][4] = {{c[0], c[1], c[2], c[3]}, {c[0], c[2], c[1], c[3]}, {c[0], c[3], c[1], c[2]}};
        double best = 1e300; int bi = 0;
        for (int s = 0; s < 3; s++) {
            double w1 = (sp[s][0] == sp[s][1]) ? 0.5 * CB[sp[s][0]].n * (CB[sp[s][0]].n + 1.0) : (double)CB[sp[s][0]].n * CB[sp[s][1]].n;
            double w2 = (sp[s][2] == sp[s][3]) ? 0.5 * CB[sp[s][2]].n * (CB[sp[s][2]].n + 1.0) : (double)CB[sp[s][2]].n * CB[sp[s][3]].n;
            if (w1 + w2 < best) { best = w1 + w2; bi = s; }
        }
        int c1 = sp[bi][0], c2 = sp[bi][1], c3 = sp[bi][2], c4 = sp[bi][3];
        printf("  multiset {%s,%s,%s,%s}: split (%s,%s | %s,%s), %.3g pair sums\n", CLSNAME[c[0]], CLSNAME[c[1]],
               CLSNAME[c[2]], CLSNAME[c[3]], CLSNAME[c1], CLSNAME[c2], CLSNAME[c3], CLSNAME[c4], best);
        int nrho = RM2;
        long long hits_m = 0; size_t pairs_m = 0;
        #pragma omp parallel for schedule(dynamic, 1) reduction(+:hits_m, pairs_m)
        for (int rho = 0; rho < nrho; rho++) {
            int nr = (int)((mod(n, RM) * RM) % RM2);              /* residue of (n, 0) */
            int sig = rsub(nr, rho);
            size_t nl = gen(c1, c2, rho, 1, n, NULL, 1), nrr = gen(c3, c4, sig, 0, n, NULL, 1);
            pairs_m += nl + nrr;
            if (!nl || !nrr) continue;
            ent *Lh = malloc(sizeof(ent) * nl), *Rh = malloc(sizeof(ent) * nrr);
            size_t mx = nl > nrr ? nl : nrr;
            ent *tmp = malloc(sizeof(ent) * mx);
            if (!Lh || !Rh || !tmp) { fprintf(stderr, "alloc\n"); exit(2); }
            gen(c1, c2, rho, 1, n, Lh, 0); gen(c3, c4, sig, 0, n, Rh, 0);
            rsort(Lh, tmp, nl); rsort(Rh, tmp, nrr);
            size_t p = 0, q = 0;
            while (p < nl && q < nrr) {
                if (Lh[p].h < Rh[q].h) p++;
                else if (Lh[p].h > Rh[q].h) q++;
                else {
                    ull h = Lh[p].h; size_t p1 = p, q1 = q;
                    while (p1 < nl && Lh[p1].h == h) p1++;
                    while (q1 < nrr && Rh[q1].h == h) q1++;
                    for (size_t u = p; u < p1; u++) for (size_t v = q; v < q1; v++) {
                        const fp *A = &CB[c1].v[Lh[u].i], *Bq = &CB[c2].v[Lh[u].j];
                        const fp *C = &CB[c3].v[Rh[v].i], *D = &CB[c4].v[Rh[v].j];
                        ll lx = n - (A->x + Bq->x), ly = -(A->y + Bq->y), rx = C->x + D->x, ry = C->y + D->y;
                        if (lx != rx || ly != ry) continue;          /* hash collision */
                        int aa[4] = {A->a, Bq->a, C->a, D->a}, bb[4] = {A->b, Bq->b, C->b, D->b};
                        if (!verify4(n, aa, bb)) { fprintf(stderr, "VERIFY FAILED\n"); exit(6); }
                        hits_m++;
                        #pragma omp critical
                        {
                            if (shown < maxwit) {
                                printf("    hit: (%d,%d) (%d,%d) (%d,%d) (%d,%d)\n", aa[0], bb[0], aa[1], bb[1], aa[2], bb[2], aa[3], bb[3]);
                                shown++;
                            }
                        }
                    }
                    p = p1; q = q1;
                }
            }
            free(Lh); free(Rh); free(tmp);
        }
        tothits += hits_m; totpairs += pairs_m;
    }
    printf("deep n = %lld, B = %d: %zu pair sums processed; %lld hits (multiset-split solutions, verified) [%.0f s wall]\n",
           n, B, totpairs, tothits, omp_get_wtime() - t0);
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: win B lo hi | deep n B [maxwit] | classes n\n"); return 1; }
    if (!strcmp(argv[1], "win")) return mode_win(atoi(argv[2]), atoll(argv[3]), atoll(argv[4]));
    if (!strcmp(argv[1], "deep")) return mode_deep(atoll(argv[2]), atoi(argv[3]), argc > 4 ? atoi(argv[4]) : 3);
    if (!strcmp(argv[1], "deepas")) { PATTERN_N = atoll(argv[3]); return mode_deep(atoll(argv[2]), atoi(argv[4]), 3); }
    if (!strcmp(argv[1], "classes")) {
        ms M[256]; ll n = atoll(argv[2]); int k = admissible(n, M);
        printf("n = %lld (mod 16 = %lld): %d admissible multisets:", n, mod(n, 16), k);
        for (int m = 0; m < k; m++) printf(" {%s,%s,%s,%s}", CLSNAME[M[m].c[0]], CLSNAME[M[m].c[1]], CLSNAME[M[m].c[2]], CLSNAME[M[m].c[3]]);
        printf("\n");
        return 0;
    }
    return 1;
}
