/* nonrigid_primes.c -- every NON-RIGID prime tournament of order n, by invariant enumeration.
 *
 * A non-rigid tournament has an automorphism g of odd prime order p (Cauchy; automorphisms of
 * tournaments have odd order).  If T is prime and n > p, g has at least two p-cycles (cycle lemma:
 * a single p-cycle C with fixed points elsewhere is a module, since every fixed vertex v satisfies
 * v -> c iff v -> g(c)).  So g is conjugate to one of the representatives
 *     g = (c p-cycles, f = n - cp fixed points),  c >= 2, or c = 1 and f = 0 (n = p, circulants),
 * and every non-rigid prime class has a labelled representative invariant under one of them.
 * The g-invariant tournaments are parametrised by one bit per <g>-orbit on pairs:
 *     fixed x fixed  : the tournament T[F], taken up to Sym(F) (Sym(F) centralises g),
 *     fixed x cycle i: v -> C_i or C_i -> v                               (c f bits),
 *     inside cycle i : a circulant on Z_p, one bit per difference 1..(p-1)/2   (c (p-1)/2 bits),
 *     cycle i x cycle j: one bit per offset e in Z_p, c_{i,s} -> c_{j,s+e} or back  (C(c,2) p bits).
 * For each prime one the nauty canonical form is deduplicated; for each new class we enumerate Aut(T)
 * from nauty's generators, test self-converseness, and list the cycle types of Aut and of the
 * anti-automorphism coset.  Output lines use the TYPE format of klein's tourn_types.c (THM-4814), so the
 * species programs read them unchanged; the rigid primes are not enumerated (their number comes from
 * the ledger: THM-4816).
 *
 * Build: cc -O3 -I/opt/homebrew/opt/nauty/include nonrigid_primes.c /opt/homebrew/opt/nauty/lib/libnauty.a -o nonrigid_primes
 * Use:   ./nonrigid_primes n TYPESFILE [CLASSFILE] [only_p only_c]
 */
#define MAXN 32
#include <nauty/nauty.h>
#include <nauty/nautinv.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int N;
static uint32_t out_[MAXN];

/* ---------------- primality: no module of size 2..n-1 ---------------- */
static uint32_t closure(uint32_t V, uint32_t S) {
    int changed = 1;
    while (changed) {
        changed = 0;
        uint32_t rest = V & ~S;
        while (rest) {
            int z = __builtin_ctz(rest); rest &= rest - 1;
            uint32_t a = out_[z] & S;
            if (a != 0 && a != S) { S |= 1u << z; changed = 1; }
        }
    }
    return S;
}
static int is_prime(void) {
    uint32_t V = (N == 32) ? 0xffffffffu : ((1u << N) - 1);
    for (int x = 0; x < N; x++)
        for (int y = x + 1; y < N; y++)
            if (closure(V, (1u << x) | (1u << y)) != V) return 0;
    return 1;
}

/* ---------------- nauty ---------------- */
typedef struct { uint64_t lo, hi; } key128;
static int gen_count;
static int gens[512][MAXN];
static void autoproc(int count, int *perm, int *orbits, int numorbits, int stabvertex, int n) {
    (void)count; (void)orbits; (void)numorbits; (void)stabvertex;
    if (gen_count < 512) { memcpy(gens[gen_count], perm, n * sizeof(int)); gen_count++; }
}
static key128 canon(const uint32_t *q, int k, int *lab_out, int want_gens, double *grp) {
    graph g[MAXN], cg[MAXN];
    int lab[MAXN], ptn[MAXN], orbits[MAXN];
    DEFAULTOPTIONS_DIGRAPH(options);
    statsblk stats;
    options.getcanon = TRUE;
    gen_count = 0;
    options.userautomproc = want_gens ? autoproc : NULL;
    EMPTYGRAPH(g, 1, k);
    for (int i = 0; i < k; i++) for (int j = 0; j < k; j++) if (q[i] >> j & 1) ADDONEARC(g, i, j, 1);
    densenauty(g, lab, ptn, orbits, &options, &stats, 1, k, cg);
    key128 key = {0, 0}; int pos = 0;
    for (int i = 0; i < k; i++) for (int j = i + 1; j < k; j++) {
        if (ISELEMENT(GRAPHROW(cg, i, 1), j)) { if (pos < 64) key.lo |= 1ULL << pos; else key.hi |= 1ULL << (pos - 64); }
        pos++;
    }
    if (lab_out) memcpy(lab_out, lab, k * sizeof(int));
    if (grp) { double s = stats.grpsize1; for (int e = 0; e < stats.grpsize2; e++) s *= 10.0; *grp = s; }
    return key;
}

/* ---------------- class registry (canonical keys) ---------------- */
static key128 *hk; static unsigned char *hu; static uint64_t HCAP = 1ULL << 24;  /* raised to 2^26 for n >= 13 in main */ static long nclasses = 0;
static int reg_insert(key128 k) {
    uint64_t h = (k.lo * 0x9E3779B97F4A7C15ULL) ^ (k.hi * 0xC2B2AE3D27D4EB4FULL);
    h ^= h >> 31; h &= HCAP - 1;
    while (hu[h]) { if (hk[h].lo == k.lo && hk[h].hi == k.hi) return 0; h = (h + 1) & (HCAP - 1); }
    hu[h] = 1; hk[h] = k; nclasses++;
    if ((uint64_t)nclasses * 2 > HCAP) { fprintf(stderr, "registry full\n"); exit(5); }
    return 1;
}

/* ---------------- group enumeration (perms packed 5 bits per point; n <= 12 fits 64 bits) ---------------- */
/* use arrays instead of packing to allow n up to 16 */
#define GMAX 4096
static unsigned char grp[GMAX][MAXN]; static int grp_n;
static int perm_eq(const unsigned char *a, const unsigned char *b) { return memcmp(a, b, N) == 0; }
static int grp_find(const unsigned char *p) { for (int e = 0; e < grp_n; e++) if (perm_eq(grp[e], p)) return e; return -1; }
static void enumerate_group(void) {
    grp_n = 0;
    for (int i = 0; i < N; i++) grp[0][i] = i;
    grp_n = 1;
    for (int head = 0; head < grp_n; head++) {
        for (int gi = 0; gi < gen_count; gi++) {
            unsigned char r[MAXN];
            for (int i = 0; i < N; i++) r[i] = gens[gi][grp[head][i]];
            if (grp_find(r) < 0) { if (grp_n >= GMAX) { fprintf(stderr, "group too big\n"); exit(6); } memcpy(grp[grp_n++], r, N); }
        }
    }
}
static void cycle_type_str(const unsigned char *p, char *s) {
    int seen[MAXN] = {0}, ct[MAXN + 1] = {0};
    for (int i = 0; i < N; i++) if (!seen[i]) { int len = 0, j = i; while (!seen[j]) { seen[j] = 1; j = p[j]; len++; } ct[len]++; }
    int o = 0; s[0] = 0;
    for (int l = 1; l <= N; l++) if (ct[l]) o += sprintf(s + o, "%d^%d,", l, ct[l]);
}
static int cmpstr(const void *a, const void *b) { return strcmp(*(char * const *)a, *(char * const *)b); }
static char ctpool[GMAX][64]; static char *ctptr[GMAX];
static void ct_multiset(const int *psi, char *outbuf) {
    for (int e = 0; e < grp_n; e++) {
        unsigned char r[MAXN];
        for (int i = 0; i < N; i++) r[i] = psi ? (unsigned char)psi[grp[e][i]] : grp[e][i];
        cycle_type_str(r, ctpool[e]); ctptr[e] = ctpool[e];
    }
    qsort(ctptr, grp_n, sizeof(char *), cmpstr);
    int ob = 0; outbuf[0] = 0;
    for (int e = 0; e < grp_n; ) { int f = e; while (f < grp_n && !strcmp(ctptr[f], ctptr[e])) f++; ob += sprintf(outbuf + ob, " %s:%d", ctptr[e], f - e); e = f; }
}

/* ---------------- type aggregation ---------------- */
#define TCAP 65536
static char *tkeys[TCAP]; static long tcount[TCAP]; static int ntypes = 0;
static void type_add(const char *s) {
    uint64_t h = 1469598103934665603ULL; for (const char *p = s; *p; p++) { h ^= (unsigned char)*p; h *= 1099511628211ULL; }
    uint32_t i = (uint32_t)(h >> 40) & (TCAP - 1);
    for (;;) { if (!tkeys[i]) { tkeys[i] = strdup(s); tcount[i] = 1; ntypes++; return; } if (!strcmp(tkeys[i], s)) { tcount[i]++; return; } i = (i + 1) & (TCAP - 1); }
}

static FILE *fclass = NULL;
static long autohist[GMAX + 1], sc_count = 0;
static void new_class(key128 key) {
    int lab[MAXN], labc[MAXN]; double gs;
    key128 k2 = canon(out_, N, lab, 1, &gs);
    if (k2.lo != key.lo || k2.hi != key.hi) { fprintf(stderr, "canon mismatch\n"); exit(8); }
    enumerate_group();
    if ((int)(gs + 0.5) != grp_n) { fprintf(stderr, "group order mismatch %g vs %d\n", gs, grp_n); exit(7); }
    uint32_t conv[MAXN];
    for (int i = 0; i < N; i++) { conv[i] = 0; for (int j = 0; j < N; j++) if (i != j && !(out_[i] >> j & 1)) conv[i] |= 1u << j; }
    key128 kc = canon(conv, N, labc, 0, NULL);
    int isc = (kc.lo == key.lo && kc.hi == key.hi);
    static char sa[1 << 16], sb[1 << 16], sig[1 << 17];
    ct_multiset(NULL, sa);
    if (isc) {
        int psi[MAXN]; for (int i = 0; i < N; i++) psi[lab[i]] = labc[i];
        /* psi maps T -> T^op: check it */
        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) if (i != j && (out_[i] >> j & 1) && (out_[psi[i]] >> psi[j] & 1)) { fprintf(stderr, "anti check failed\n"); exit(9); }
        ct_multiset(psi, sb);
    } else sb[0] = 0;
    snprintf(sig, sizeof sig, "k=%d aut=%d AUT:%s SC=%d ANTI:%s", N, grp_n, sa, isc, sb);
    type_add(sig);
    autohist[grp_n < GMAX ? grp_n : GMAX]++; sc_count += isc;
    if (fclass) {
        int sc[MAXN]; for (int i = 0; i < N; i++) sc[i] = __builtin_popcount(out_[i]);
        /* sorted score sequence */
        for (int i = 1; i < N; i++) { int v = sc[i], j = i - 1; while (j >= 0 && sc[j] > v) { sc[j + 1] = sc[j]; j--; } sc[j + 1] = v; }
        fprintf(fclass, "%016llx%016llx aut=%d sc=%d scores=", (unsigned long long)key.hi, (unsigned long long)key.lo, grp_n, isc);
        for (int i = 0; i < N; i++) fprintf(fclass, "%d%s", sc[i], i + 1 < N ? "," : "");
        fprintf(fclass, " arcs=");
        for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++) fputc((out_[i] >> j & 1) ? '1' : '0', fclass);
        fprintf(fclass, "\n");
    }
}

/* ---------------- fixed-point tournament representatives up to Sym(F) ---------------- */
static uint32_t frep[2048][MAXN]; static int nfrep;
static void build_frep(int f) {
    nfrep = 0;
    if (f <= 1) { for (int i = 0; i < MAXN; i++) frep[0][i] = 0; nfrep = 1; return; }
    int np = f * (f - 1) / 2;
    /* registry local */
    key128 seen[2048]; int nseen = 0;
    for (uint64_t m = 0; m < (1ULL << np); m++) {
        uint32_t q[MAXN] = {0}; int b = 0;
        for (int i = 0; i < f; i++) for (int j = i + 1; j < f; j++) { if (m >> b & 1) q[i] |= 1u << j; else q[j] |= 1u << i; b++; }
        key128 k = canon(q, f, NULL, 0, NULL);
        int dup = 0; for (int s = 0; s < nseen; s++) if (seen[s].lo == k.lo && seen[s].hi == k.hi) { dup = 1; break; }
        if (dup) continue;
        seen[nseen++] = k; memcpy(frep[nfrep++], q, sizeof q);
    }
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s n TYPESFILE [CLASSFILE] [only_p only_c]\n", argv[0]); return 1; }
    N = atoi(argv[1]);
    if (N > 16) { fprintf(stderr, "n <= 16\n"); return 1; }
    FILE *ft = fopen(argv[2], "w");
    if (argc > 3 && strcmp(argv[3], "-")) fclass = fopen(argv[3], "w");
    int only_p = argc > 5 ? atoi(argv[4]) : 0, only_c = argc > 5 ? atoi(argv[5]) : 0;
    if (N >= 13) HCAP = 1ULL << 26;
    hk = calloc(HCAP, sizeof(key128)); hu = calloc(HCAP, 1);
    int primes_list[] = {3, 5, 7, 11, 13};
    long total_enum = 0, total_prime = 0;
    for (int pi = 0; pi < 5; pi++) {
        int p = primes_list[pi];
        for (int c = 1; c * p <= N; c++) {
            int f = N - c * p;
            if (!(c >= 2 || f == 0)) continue;
            if (only_p && (p != only_p || c != only_c)) continue;
            build_frep(f);
            /* vertex layout: fixed 0..f-1; cycle i, position t -> f + i*p + t */
            int nb = c * f + c * (p - 1) / 2 + (c * (c - 1) / 2) * p;
            long cnt_enum = 0, cnt_prime = 0, cnt_new = 0;
            for (int r = 0; r < nfrep; r++) {
                for (uint64_t m = 0; m < (1ULL << nb); m++) {
                    for (int i = 0; i < N; i++) out_[i] = 0;
                    for (int i = 0; i < f; i++) out_[i] = frep[r][i];
                    int b = 0;
                    for (int v = 0; v < f; v++) for (int i = 0; i < c; i++) {
                        int bit = m >> b & 1; b++;
                        for (int t = 0; t < p; t++) { int w = f + i * p + t; if (bit) out_[v] |= 1u << w; else out_[w] |= 1u << v; }
                    }
                    for (int i = 0; i < c; i++) for (int d = 1; d <= (p - 1) / 2; d++) {
                        int bit = m >> b & 1; b++;
                        for (int s = 0; s < p; s++) { int u = f + i * p + s, w = f + i * p + (s + d) % p; if (bit) out_[u] |= 1u << w; else out_[w] |= 1u << u; }
                    }
                    for (int i = 0; i < c; i++) for (int j = i + 1; j < c; j++) for (int e = 0; e < p; e++) {
                        int bit = m >> b & 1; b++;
                        for (int s = 0; s < p; s++) { int u = f + i * p + s, w = f + j * p + (s + e) % p; if (bit) out_[u] |= 1u << w; else out_[w] |= 1u << u; }
                    }
                    cnt_enum++;
                    if (!is_prime()) continue;
                    cnt_prime++;
                    key128 key = canon(out_, N, NULL, 0, NULL);
                    if (reg_insert(key)) { cnt_new++; new_class(key); }
                }
            }
            fprintf(stderr, "n=%d p=%d c=%d f=%d: F-reps=%d bits=%d enumerated=%ld prime=%ld new classes=%ld\n", N, p, c, f, nfrep, nb, cnt_enum, cnt_prime, cnt_new);
            printf("CASE n=%d p=%d c=%d f=%d Freps=%d bits=%d enumerated=%ld prime=%ld new=%ld\n", N, p, c, f, nfrep, nb, cnt_enum, cnt_prime, cnt_new);
            total_enum += cnt_enum; total_prime += cnt_prime;
        }
    }
    printf("n=%d nonrigid_prime_classes=%ld selfconverse=%ld enumerated=%ld prime_labelled=%ld distinct_types=%d\n", N, nclasses, sc_count, total_enum, total_prime, ntypes);
    printf("aut_hist:"); for (int g = 2; g <= GMAX; g++) if (autohist[g]) printf(" %d:%ld", g, autohist[g]); printf("\n");
    for (int i = 0; i < TCAP; i++) if (tkeys[i]) fprintf(ft, "TYPE count=%ld %s\n", tcount[i], tkeys[i]);
    fclose(ft); if (fclass) fclose(fclass);
    return 0;
}
