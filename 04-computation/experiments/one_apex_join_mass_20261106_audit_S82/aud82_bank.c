/*
 * aud82_bank.c  (independent audit of opus-2026-11-06-S82; shares no code with the session or with S78)
 *
 * Rooted triples (h, e_t, R_t) of tournaments:
 *   h      = number of Hamiltonian paths of T,
 *   e_t    = number of Hamiltonian paths of T ending at t,
 *   R_t    = sum over X subset V - t of h(T[X]) * s_t(T[V - X])   (h(empty) = 1; s_t = paths starting at t).
 * Derivation of R_t from the one-arc join (A => B except b -> a reversed): the paths through b -> a are
 *   A1 B1 a A2 B2  with  A1 = path on X subset A - a,  a A2 = path on A - X from a,  B1 ending at b,  B2 the rest.
 *
 * Mode 1 (labelled):  all 2^(n(n-1)/2) labelled tournaments of order n <= 7, rooted at vertex 0.
 * Mode 2 (classes):   isomorphism-class representatives built order by order (extend every rep of order n-1 by a new
 *                     vertex in all 2^(n-1) ways; canonical form = lexicographically least arc code over all n!
 *                     relabellings), then the triple at the new vertex t for every (rep of order n-1, arc pattern of t):
 *                     this is every rooted tournament (T, t) of order n up to isomorphism.  Used for orders 1..8.
 * Output: per order m the number of distinct triples of order <= m and M_m = sum (2/y_h) ln(1 + y_h/y_R);
 *         files aud82_bank_le<m>.txt (sorted "h e R").
 * Usage: aud82_bank NMAX_LABELLED NMAX_CLASSES
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <omp.h>
typedef long long ll;
typedef struct { ll h, e, R; } Tri;

/* adjacency: out[v] bitmask of out-neighbours */
static void triple_at(int n, const unsigned *out, int t, Tri *res) {
    static __thread ll endp[1 << 9][9], stp[1 << 9][9];
    int full = (1 << n) - 1;
    /* end-DP: endp[S][v] = paths of T[S] ending at v */
    for (int S = 1; S <= full; S++) for (int v = 0; v < n; v++) endp[S][v] = 0;
    for (int v = 0; v < n; v++) endp[1 << v][v] = 1;
    for (int S = 1; S <= full; S++)
        for (int v = 0; v < n; v++) {
            ll c = endp[S][v]; if (!c || !((S >> v) & 1)) continue;
            unsigned nb = out[v] & ~(unsigned)S;
            while (nb) { int w = __builtin_ctz(nb); nb &= nb - 1; endp[S | (1 << w)][w] += c; }
        }
    /* start-DP from t: stp[S][v] = paths of T[S] starting at t, ending at v (S contains t) */
    for (int S = 0; S <= full; S++) for (int v = 0; v < n; v++) stp[S][v] = 0;
    stp[1 << t][t] = 1;
    for (int S = 1; S <= full; S++) {
        if (!((S >> t) & 1)) continue;
        for (int v = 0; v < n; v++) {
            ll c = stp[S][v]; if (!c) continue;
            unsigned nb = out[v] & ~(unsigned)S;
            while (nb) { int w = __builtin_ctz(nb); nb &= nb - 1; stp[S | (1 << w)][w] += c; }
        }
    }
    ll h = 0; for (int v = 0; v < n; v++) h += endp[full][v];
    ll e = endp[full][t];
    ll R = 0;
    for (int Y = 1; Y <= full; Y++) {
        if (!((Y >> t) & 1)) continue;
        ll sY = 0; for (int v = 0; v < n; v++) sY += stp[Y][v];
        if (!sY) continue;
        int X = full & ~Y; ll hX = 1;
        if (X) { hX = 0; for (int v = 0; v < n; v++) hX += endp[X][v]; }
        R += hX * sY;
    }
    res->h = h; res->e = e; res->R = R;
}

static int cmp_tri(const void *a, const void *b) {
    const Tri *x = a, *y = b;
    if (x->h != y->h) return x->h < y->h ? -1 : 1;
    if (x->e != y->e) return x->e < y->e ? -1 : 1;
    if (x->R != y->R) return x->R < y->R ? -1 : 1;
    return 0;
}
static ll dedup(Tri *a, ll n) {
    qsort(a, n, sizeof(Tri), cmp_tri);
    ll m = 0;
    for (ll i = 0; i < n; i++) if (m == 0 || cmp_tri(&a[m - 1], &a[i])) a[m++] = a[i];
    return m;
}

/* ---- canonical form by brute force over all permutations (n <= 8) ---- */
static unsigned long long code_of(int n, const unsigned *out, const int *perm) {
    /* perm[i] = old vertex placed at new position i; code bit for pair (i<j): 1 iff new i -> new j */
    unsigned long long c = 0;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) c = (c << 1) | ((out[perm[i]] >> perm[j]) & 1u);
    return c;
}
static unsigned long long canon(int n, const unsigned *out) {
    int perm[9]; for (int i = 0; i < n; i++) perm[i] = i;
    unsigned long long best = ~0ULL;
    /* Heap's algorithm over all n! permutations */
    int c[9] = {0};
    unsigned long long cc = code_of(n, out, perm); if (cc < best) best = cc;
    int i = 0;
    while (i < n) {
        if (c[i] < i) {
            if (i % 2 == 0) { int t = perm[0]; perm[0] = perm[i]; perm[i] = t; }
            else { int t = perm[c[i]]; perm[c[i]] = perm[i]; perm[i] = t; }
            cc = code_of(n, out, perm); if (cc < best) best = cc;
            c[i]++; i = 0;
        } else { c[i] = 0; i++; }
    }
    return best;
}
static void decode(int n, unsigned long long code, unsigned *out) {
    for (int v = 0; v < n; v++) out[v] = 0;
    int npairs = n * (n - 1) / 2, k = npairs - 1;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
        int bit = (code >> k) & 1; k--;
        if (bit) out[i] |= 1u << j; else out[j] |= 1u << i;
    }
}
static int cmp_u64(const void *a, const void *b) {
    unsigned long long x = *(const unsigned long long *)a, y = *(const unsigned long long *)b;
    return x < y ? -1 : x > y;
}

int main(int argc, char **argv) {
    int NL = argc > 1 ? atoi(argv[1]) : 7, NC = argc > 2 ? atoi(argv[2]) : 8;
    /* ---------- mode 1: labelled ---------- */
    Tri *lab_all = NULL; ll lab_n = 0;
    for (int n = 1; n <= NL; n++) {
        int np = n * (n - 1) / 2; ll cnt = 1LL << np;
        int pi[64], pj[64], k = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { pi[k] = i; pj[k] = j; k++; }
        int nth = omp_get_max_threads();
        Tri **loc = malloc(sizeof(Tri *) * nth); ll *locn = calloc(nth, sizeof(ll)); ll *cap = malloc(sizeof(ll) * nth);
        for (int q = 0; q < nth; q++) { cap[q] = 1 << 16; loc[q] = malloc(sizeof(Tri) * cap[q]); }
        #pragma omp parallel
        {
            int q = omp_get_thread_num();
            #pragma omp for schedule(dynamic, 4096)
            for (ll code = 0; code < cnt; code++) {
                unsigned out[9] = {0};
                for (int a = 0; a < np; a++) { if ((code >> a) & 1) out[pi[a]] |= 1u << pj[a]; else out[pj[a]] |= 1u << pi[a]; }
                Tri t; triple_at(n, out, 0, &t);
                if (locn[q] == cap[q]) { locn[q] = dedup(loc[q], locn[q]); if (locn[q] > cap[q] / 2) { cap[q] *= 2; loc[q] = realloc(loc[q], sizeof(Tri) * cap[q]); } }
                loc[q][locn[q]++] = t;
            }
        }
        for (int q = 0; q < nth; q++) {
            lab_all = realloc(lab_all, sizeof(Tri) * (lab_n + locn[q]));
            memcpy(lab_all + lab_n, loc[q], sizeof(Tri) * locn[q]); lab_n += locn[q]; free(loc[q]);
        }
        free(loc); free(locn); free(cap);
        lab_n = dedup(lab_all, lab_n);
        long double M = 0; for (ll i = 0; i < lab_n; i++) M += (2.0L / lab_all[i].h) * logl(1.0L + (long double)lab_all[i].h / lab_all[i].R);
        printf("LABELLED order<=%d: %lld labelled tournaments at order %d rooted at 0; distinct triples of order <= %d: %lld; M = %.9Lf\n",
               n, cnt, n, n, lab_n, M);
        fflush(stdout);
    }
    /* ---------- mode 2: classes ---------- */
    unsigned long long *reps = malloc(sizeof(unsigned long long)); reps[0] = 0; ll nreps = 1;   /* order 1 */
    Tri *cls_all = malloc(sizeof(Tri)); cls_all[0].h = 1; cls_all[0].e = 1; cls_all[0].R = 1; ll cls_n = 1;
    printf("CLASSES order 1: 1 class; distinct triples of order <= 1: 1; M = %.9Lf\n", 2.0L * logl(2.0L));
    for (int n = 2; n <= NC; n++) {
        /* rooted tournaments of order n: (rep U of order n-1) + new vertex t = n-1 with pattern p (bit i: t -> i) */
        ll ncand = nreps << (n - 1);
        Tri *tr = malloc(sizeof(Tri) * ncand);
        unsigned long long *cand = (n <= 7) ? malloc(sizeof(unsigned long long) * ncand) : NULL;
        #pragma omp parallel for schedule(dynamic, 64)
        for (ll idx = 0; idx < ncand; idx++) {
            ll r = idx >> (n - 1); unsigned p = (unsigned)(idx & ((1LL << (n - 1)) - 1));
            unsigned out[9] = {0};
            decode(n - 1, reps[r], out);
            int t = n - 1;
            for (int i = 0; i < n - 1; i++) { if ((p >> i) & 1) out[t] |= 1u << i; else out[i] |= 1u << t; }
            triple_at(n, out, t, &tr[idx]);
            if (cand) cand[idx] = canon(n, out);
        }
        cls_all = realloc(cls_all, sizeof(Tri) * (cls_n + ncand));
        memcpy(cls_all + cls_n, tr, sizeof(Tri) * ncand); cls_n += ncand; free(tr);
        cls_n = dedup(cls_all, cls_n);
        long double M = 0; for (ll i = 0; i < cls_n; i++) M += (2.0L / cls_all[i].h) * logl(1.0L + (long double)cls_all[i].h / cls_all[i].R);
        if (cand) {
            qsort(cand, ncand, sizeof(unsigned long long), cmp_u64);
            ll m = 0; for (ll i = 0; i < ncand; i++) if (m == 0 || cand[m - 1] != cand[i]) cand[m++] = cand[i];
            free(reps); reps = cand; nreps = m;
            printf("CLASSES order %d: %lld isomorphism classes; ", n, nreps);
        } else printf("CLASSES order %d: (classes not enumerated; %lld rooted candidates) ", n, ncand);
        printf("distinct triples of order <= %d: %lld; M = %.9Lf\n", n, cls_n, M);
        char fn[64]; sprintf(fn, "aud82_bank_le%d.txt", n);
        FILE *f = fopen(fn, "w"); for (ll i = 0; i < cls_n; i++) fprintf(f, "%lld %lld %lld\n", cls_all[i].h, cls_all[i].e, cls_all[i].R); fclose(f);
        if (n <= NL) {
            /* compare with the labelled set of order <= NL only at n == NL */
        }
        fflush(stdout);
        if (n == NL) {
            /* the labelled set is for order <= NL: compare */
            int same = (cls_n == lab_n);
            for (ll i = 0; same && i < cls_n; i++) if (cmp_tri(&cls_all[i], &lab_all[i])) same = 0;
            printf("  labelled-mode set (order <= %d) == class-mode set: %s\n", NL, same ? "YES" : "NO");
        }
        /* extra statistics */
        ll maxh = 0, maxR = 0, bad_e = 0; double sinv = 0;
        for (ll i = 0; i < cls_n; i++) { if (cls_all[i].h > maxh) maxh = cls_all[i].h; if (cls_all[i].R > maxR) maxR = cls_all[i].R; if (cls_all[i].e > cls_all[i].h || !(cls_all[i].h & 1)) bad_e++; sinv += 1.0 / cls_all[i].R; }
        printf("  max h %lld, max R %lld, #(e > h or h even) %lld, sum 1/R = %.4f\n", maxh, maxR, bad_e, sinv);
    }
    return 0;
}
