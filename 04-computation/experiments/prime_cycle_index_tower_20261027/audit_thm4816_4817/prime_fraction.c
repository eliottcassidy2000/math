/*
 * audit_census.c -- independent auditor's census of the NON-RIGID prime tournaments of order n
 * (audit of THM-4816; written from scratch, without reading nonrigid_primes.c).
 *
 * Method (THM-4816 A1/A2, re-derived by the auditor):
 *   every non-rigid prime Q has an automorphism of odd prime order p, cycle type 1^f p^c (n = pc + f),
 *   and by the cycle lemma c >= 2, or c = 1 and f = 0.  Relabel so that this automorphism is the fixed
 *   permutation g0: fixed points 0..f-1, cycle i = { f + i*p + s : s in Z_p },  s -> s+1.
 *   A g0-invariant tournament is one bit per <g0>-orbit of pairs, plus T[F]; T[F] is taken up to
 *   isomorphism (Sym(F) centralises g0).  T[F] representatives are generated here by brute force over
 *   all labelled tournaments on f vertices + nauty canonical form (independent of gentourng).
 *   Each invariant tournament is tested for primality (closure of every pair = V), canonicalised with
 *   nauty (digraph), deduplicated by its canonical upper triangle (exact, 128-bit key), and for every
 *   NEW class the automorphism group is closed from nauty's generators (checked against grpsize), and
 *   self-converseness / anti-automorphisms are computed from canon(T) == canon(T^op).
 *
 * Usage: ./audit_census n types_out.txt [classes_out.txt|-] [stats_out.txt]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <nauty/nauty.h>
#include <nauty/nautinv.h>

#define NV 16          /* max order handled */
typedef uint16_t row_t; /* out-neighbourhood bitmask: bit j of out[i] set iff i -> j */

static int M_WORDS = 1;
static graph ng[NV], ncan[NV];
static int nlab[NV], nptn[NV], norb[NV];
static int gens[256][NV];
static int ngens;
static double last_grpsize;

static void storegen(int count, int *perm, int *orbits, int numorbits, int stabvertex, int nn)
{
    (void)count; (void)orbits; (void)numorbits; (void)stabvertex;
    if (ngens < 256) { memcpy(gens[ngens], perm, nn * sizeof(int)); ngens++; }
    else { fprintf(stderr, "too many generators\n"); exit(2); }
}

/* canonical form; key = packed upper triangle (i<j: bit set iff i->j in canonical graph) */
static void canon(const row_t *out, int nn, uint64_t key[2], int *labout)
{
    static DEFAULTOPTIONS_DIGRAPH(options);
    statsblk stats;
    options.getcanon = TRUE;
    options.userautomproc = storegen;
    EMPTYGRAPH(ng, M_WORDS, nn);
    for (int i = 0; i < nn; i++)
        for (int j = 0; j < nn; j++)
            if ((out[i] >> j) & 1) ADDONEARC(ng, i, j, M_WORDS);
    ngens = 0;
    densenauty(ng, nlab, nptn, norb, &options, &stats, M_WORDS, nn, ncan);
    last_grpsize = stats.grpsize1;
    for (int k = 0; k < stats.grpsize2; k++) last_grpsize *= 10.0;
    key[0] = key[1] = 0;
    int b = 0;
    for (int i = 0; i < nn; i++) {
        set *gi = GRAPHROW(ncan, i, M_WORDS);
        for (int j = i + 1; j < nn; j++, b++)
            if (ISELEMENT(gi, j)) key[b >> 6] |= (1ULL << (b & 63));
    }
    key[1] |= (1ULL << 63); /* occupancy marker (n <= 14 uses < 92 bits) */
    if (labout) memcpy(labout, nlab, nn * sizeof(int));
}

/* -------- primality: the smallest module containing each pair must be V -------- */
static int is_prime(const row_t *out, int nn)
{
    row_t full = (row_t)((1u << nn) - 1);
    for (int u = 0; u < nn; u++)
        for (int v = u + 1; v < nn; v++) {
            row_t M = (row_t)((1u << u) | (1u << v));
            for (;;) {
                row_t add = 0, rest = (row_t)(full & ~M);
                while (rest) {
                    int w = __builtin_ctz(rest);
                    rest &= (row_t)(rest - 1);
                    row_t x = out[w] & M;
                    if (x && x != M) add |= (row_t)(1u << w);
                }
                if (!add) break;
                M |= add;
            }
            if (M != full) return 0;
        }
    return 1;
}

/* -------- hash set of canonical keys -------- */
static uint64_t *HT = NULL;
static size_t HTcap = 0, HTcount = 0;
static uint64_t hmix(const uint64_t k[2])
{
    uint64_t h = k[0] * 0x9E3779B97F4A7C15ULL ^ (k[1] + 0x632BE59BD9B4E019ULL);
    h ^= h >> 31; h *= 0xBF58476D1CE4E5B9ULL; h ^= h >> 29; h *= 0x94D049BB133111EBULL; h ^= h >> 32;
    return h;
}
static void ht_init(size_t cap)
{
    HTcap = cap; HTcount = 0;
    HT = calloc(2 * HTcap, sizeof(uint64_t));
    if (!HT) { fprintf(stderr, "alloc fail\n"); exit(3); }
}
static int ht_insert(const uint64_t k[2]); /* returns 1 if new */
static void ht_grow(void)
{
    uint64_t *old = HT; size_t oldcap = HTcap;
    ht_init(oldcap * 2);
    for (size_t i = 0; i < oldcap; i++)
        if (old[2 * i + 1]) { uint64_t k[2] = {old[2 * i], old[2 * i + 1]}; ht_insert(k); }
    free(old);
}
static int ht_insert(const uint64_t k[2])
{
    if ((HTcount + 1) * 10 > HTcap * 7) ht_grow();
    size_t i = hmix(k) & (HTcap - 1);
    for (;;) {
        if (!HT[2 * i + 1]) { HT[2 * i] = k[0]; HT[2 * i + 1] = k[1]; HTcount++; return 1; }
        if (HT[2 * i] == k[0] && HT[2 * i + 1] == k[1]) return 0;
        i = (i + 1) & (HTcap - 1);
    }
}

/* -------- permutation groups -------- */
#define GMAX 4096
static unsigned char grp[GMAX][NV];
static int gsize;
static int perm_eq(const unsigned char *a, const unsigned char *b, int nn) { return memcmp(a, b, nn) == 0; }
static void close_group(int nn)
{
    gsize = 1;
    for (int i = 0; i < nn; i++) grp[0][i] = (unsigned char)i;
    for (int idx = 0; idx < gsize; idx++)
        for (int g = 0; g < ngens; g++) {
            unsigned char prod[NV];
            for (int i = 0; i < nn; i++) prod[i] = (unsigned char)gens[g][grp[idx][i]]; /* gen o elem */
            int found = 0;
            for (int t = 0; t < gsize; t++) if (perm_eq(grp[t], prod, nn)) { found = 1; break; }
            if (!found) {
                if (gsize >= GMAX) { fprintf(stderr, "group too large\n"); exit(4); }
                memcpy(grp[gsize++], prod, nn);
            }
        }
}
static int is_aut(const row_t *out, const unsigned char *p, int nn)
{
    for (int x = 0; x < nn; x++)
        for (int y = 0; y < nn; y++)
            if (x != y && (((out[x] >> y) & 1) != ((out[p[x]] >> p[y]) & 1))) return 0;
    return 1;
}
static int is_anti(const row_t *out, const unsigned char *p, int nn)
{
    for (int x = 0; x < nn; x++)
        for (int y = 0; y < nn; y++)
            if (x != y && (((out[x] >> y) & 1) != ((out[p[y]] >> p[x]) & 1))) return 0;
    return 1;
}
/* cycle type as string "1^5.3^2" (lengths ascending) */
static void cycle_type(const unsigned char *p, int nn, char *buf)
{
    int seen[NV] = {0}, cnt[NV + 1] = {0};
    for (int i = 0; i < nn; i++) {
        if (seen[i]) continue;
        int len = 0, j = i;
        while (!seen[j]) { seen[j] = 1; j = p[j]; len++; }
        cnt[len]++;
    }
    buf[0] = 0;
    for (int l = 1; l <= nn; l++)
        if (cnt[l]) { char t[16]; sprintf(t, "%s%d^%d", buf[0] ? "." : "", l, cnt[l]); strcat(buf, t); }
}

/* multiset of cycle types -> canonical string "ct:k ct:k" sorted */
typedef struct { char ct[64]; int k; } ctcount;
static int cmp_ct(const void *a, const void *b) { return strcmp(((const ctcount *)a)->ct, ((const ctcount *)b)->ct); }
static void multiset_string(char cts[][64], int cnt, char *out)
{
    ctcount arr[GMAX]; int na = 0;
    for (int i = 0; i < cnt; i++) {
        int f = -1;
        for (int j = 0; j < na; j++) if (!strcmp(arr[j].ct, cts[i])) { f = j; break; }
        if (f < 0) { strcpy(arr[na].ct, cts[i]); arr[na].k = 1; na++; } else arr[f].k++;
    }
    qsort(arr, na, sizeof(ctcount), cmp_ct);
    out[0] = 0;
    for (int i = 0; i < na; i++) { char t[96]; sprintf(t, "%s%s:%d", i ? " " : "", arr[i].ct, arr[i].k); strcat(out, t); }
}

/* -------- type table -------- */
#define TMAX 512
static char tkey[TMAX][1024];
static long long tcount[TMAX];
static int ntypes;
static void add_type(const char *key)
{
    for (int i = 0; i < ntypes; i++) if (!strcmp(tkey[i], key)) { tcount[i]++; return; }
    if (ntypes >= TMAX) { fprintf(stderr, "too many types\n"); exit(5); }
    strcpy(tkey[ntypes], key); tcount[ntypes] = 1; ntypes++;
}

/* -------- tournament classes on f vertices (brute force) -------- */
#define FMAX 8
static row_t *reps[FMAX + 1];
static int nreps[FMAX + 1];
static void make_reps(int f)
{
    if (f <= 1) { reps[f] = calloc(NV, sizeof(row_t)); nreps[f] = 1; return; }
    int np = f * (f - 1) / 2;
    long total = 1L << np;
    /* small private hash on keys */
    size_t cap = 1; while (cap < 4096) cap <<= 1;
    uint64_t *h = calloc(2 * cap, sizeof(uint64_t));
    row_t *store = malloc(sizeof(row_t) * NV * 1024);
    int ns = 0;
    for (long x = 0; x < total; x++) {
        row_t out[NV] = {0};
        int b = 0;
        for (int i = 0; i < f; i++)
            for (int j = i + 1; j < f; j++, b++) {
                if ((x >> b) & 1) out[i] |= (row_t)(1u << j); else out[j] |= (row_t)(1u << i);
            }
        uint64_t key[2];
        canon(out, f, key, NULL);
        size_t i = hmix(key) & (cap - 1);
        int isnew = 0;
        for (;;) {
            if (!h[2 * i + 1]) { h[2 * i] = key[0]; h[2 * i + 1] = key[1]; isnew = 1; break; }
            if (h[2 * i] == key[0] && h[2 * i + 1] == key[1]) break;
            i = (i + 1) & (cap - 1);
        }
        if (isnew) {
            /* store the canonical graph itself as the representative */
            row_t *r = store + NV * ns;
            memset(r, 0, sizeof(row_t) * NV);
            for (int a = 0; a < f; a++) {
                set *ga = GRAPHROW(ncan, a, M_WORDS);
                for (int c = 0; c < f; c++) if (ISELEMENT(ga, c)) r[a] |= (row_t)(1u << c);
            }
            ns++;
            if (ns >= 1024) { fprintf(stderr, "too many reps\n"); exit(6); }
        }
    }
    free(h);
    reps[f] = store; nreps[f] = ns;
}


#define CV(i, s) (f + (i) * p + (((s) % p + p) % p))
#define SETARC(arr, a, b) do { (arr)[a] |= (row_t)(1u << (b)); (arr)[b] &= (row_t)~(1u << (a)); } while (0)
#define FLIP(arr, a, b) do { (arr)[a] ^= (row_t)(1u << (b)); (arr)[b] ^= (row_t)(1u << (a)); } while (0)

/* prime_fraction: for the type 1^f 3^2 (n = f + 6) count prime g0-invariant tournaments per T[F] class R and
   weight by the number f!/|Aut R| of labellings of R: labelled prime fraction of the invariant tournaments. */
int main(int argc, char **argv)
{
    int n = atoi(argv[1]); int p = 3, c = 2, f = n - 6;
    nauty_check(WORDSIZE, M_WORDS, NV, NAUTYVERSIONID);
    make_reps(f);
    int nb = 0; static row_t delta[64][NV]; row_t base[NV]; memset(base, 0, sizeof base); memset(delta, 0, sizeof delta);
    for (int v = 0; v < f; v++) for (int i = 0; i < c; i++) { for (int s = 0; s < p; s++) { SETARC(base, CV(i, s), v); FLIP(delta[nb], CV(i, s), v); } nb++; }
    for (int i = 0; i < c; i++) for (int d = 1; d <= (p - 1) / 2; d++) { for (int s = 0; s < p; s++) { SETARC(base, CV(i, s + d), CV(i, s)); FLIP(delta[nb], CV(i, s + d), CV(i, s)); } nb++; }
    for (int i = 0; i < c; i++) for (int j = i + 1; j < c; j++) for (int e = 0; e < p; e++) { for (int s = 0; s < p; s++) { SETARC(base, CV(j, s + e), CV(i, s)); FLIP(delta[nb], CV(j, s + e), CV(i, s)); } nb++; }
    double fact = 1; for (int i = 2; i <= f; i++) fact *= i;
    long double wprime = 0, wtot = 0; long long uprime = 0, utot = 0;
    for (int r = 0; r < nreps[f]; r++) {
        row_t out[NV]; memcpy(out, base, sizeof out);
        for (int a = 0; a < f; a++) out[a] = (row_t)((out[a] & ~((1u << f) - 1)) | (reps[f][NV * r + a] & ((1u << f) - 1)));
        double autR = 1;
        if (f >= 2) { uint64_t k[2]; canon(reps[f] + NV * r, f, k, NULL); autR = last_grpsize; }
        long long pr = 0, tot = 1LL << nb;
        for (long long x = 0; x < tot; x++) {
            if (x) { int bit = __builtin_ctzll(x); for (int a = 0; a < n; a++) out[a] ^= delta[bit][a]; }
            pr += is_prime(out, n);
        }
        wprime += (long double)(fact / autR) * pr; wtot += (long double)(fact / autR) * tot; uprime += pr; utot += tot;
    }
    printf("n=%d f=%d bits=%d reps=%d unweighted prime fraction %.6f (%lld/%lld); labelled prime fraction %.6Lf (labelled total %.0Lf = 2^%d)\n",
           n, f, nb, nreps[f], (double)uprime / utot, uprime, utot, wprime / wtot, wtot, f * (f - 1) / 2 + nb);
    return 0;
}
