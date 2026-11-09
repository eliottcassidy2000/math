/* AUDIT G1 (independent code; shares nothing with the session's tprop.c / canon.c / sa12.c / wit.c).
 *
 * Tournament code: upper triangle, row-major over pairs (i,j), i<j; '1' means i -> j.
 * out[i] = bitmask of out-neighbours of i.
 *
 * Modes (read "n code [rest]" lines on stdin unless stated):
 *   props   : n H lambda prime moddef twin |  H(T-v) v=0..n-1  | scores        (exact, by enumeration)
 *   subs K  : H of every induced subtournament on n-K vertices (all C(n,K) deletions), one line per tournament
 *   ext PARENTCODE11 S10FILE S11FILE : all 2^11 one-vertex extensions of the given 11-tournament (S17 integer code);
 *             prints label-set statistics (FC, FCV, SUM, SUB, SUM2) and every escaper's 12-code.
 *   band s n : print band_s(n) code
 *   selftest : DP vs brute-force permutations on random tournaments n <= 8
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef unsigned long long u64;
#define MAXN 20

static int popc(unsigned x) { return __builtin_popcount(x); }

static void decode(const char *code, int n, unsigned *out) {
    int p = 0;
    for (int i = 0; i < n; i++) out[i] = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++) {
            if (code[p] == '1') out[i] |= 1u << j; else out[j] |= 1u << i;
            p++;
        }
}

static void encode(int n, const unsigned *out, char *code) {
    int p = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++) code[p++] = (out[i] >> j & 1) ? '1' : '0';
    code[p] = 0;
}

/* number of Hamiltonian paths of the subtournament induced on vertex set vs (bitmask over 0..n-1) */
static u64 *dpbuf = NULL;
static u64 ham_sub(int n, const unsigned *out, unsigned vs) {
    int idx[MAXN], k = 0;
    for (int v = 0; v < n; v++) if (vs >> v & 1) idx[k++] = v;
    if (k <= 1) return 1;
    unsigned loc[MAXN];
    for (int a = 0; a < k; a++) {
        loc[a] = 0;
        for (int b = 0; b < k; b++) if (out[idx[a]] >> idx[b] & 1) loc[a] |= 1u << b;
    }
    unsigned full = (1u << k) - 1;
    /* dp[mask*k + last] : paths covering exactly mask, ending at last.  Backward-free formulation:
       we extend at the END (append), iterating masks in increasing numeric order. */
    size_t sz = ((size_t)1 << k) * k;
    if (k > 16) { fprintf(stderr, "ham_sub: k=%d > 16 unsupported\n", k); exit(2); }
    if (!dpbuf) dpbuf = malloc(sizeof(u64) * ((size_t)1 << 16) * 16);
    u64 *dp = dpbuf;
    memset(dp, 0, sizeof(u64) * sz);
    for (int a = 0; a < k; a++) dp[((size_t)1 << a) * k + a] = 1;
    for (unsigned m = 1; m < full; m++) {
        u64 *row = dp + (size_t)m * k;
        for (int a = 0; a < k; a++) {
            u64 c = row[a];
            if (!c) continue;
            unsigned nx = loc[a] & ~m;
            while (nx) {
                int b = __builtin_ctz(nx);
                nx &= nx - 1;
                dp[(size_t)(m | (1u << b)) * k + b] += c;
            }
        }
    }
    u64 t = 0;
    for (int a = 0; a < k; a++) t += dp[(size_t)full * k + a];
    return t;
}

/* arc-strong connectivity: min over nonempty proper S of #arcs S -> V\S (all 2^n - 2 cuts) */
static int lambda_cut(int n, const unsigned *out) {
    unsigned full = (1u << n) - 1;
    int best = 1 << 30;
    for (unsigned S = 1; S < full; S++) {
        int c = 0;
        for (int v = 0; v < n; v++) if (S >> v & 1) c += popc(out[v] & ~S & full);
        if (c < best) best = c;
    }
    return best;
}

/* module defect: min over M (2 <= |M| <= n-1) of sum_{x notin M} min(|out(x)&M|, |M|-|out(x)&M|).
   moddef == 0 iff not prime.  twin = same restricted to |M| = 2.  Also returns one minimizing M. */
static int moddef(int n, const unsigned *out, int *twin, unsigned *argM) {
    unsigned full = (1u << n) - 1;
    int best = 1 << 30, bt = 1 << 30;
    unsigned bm = 0;
    for (unsigned M = 1; M < full; M++) {
        int sz = popc(M);
        if (sz < 2) continue;
        int c = 0;
        for (int x = 0; x < n && c < best + (sz == 2 ? 1000 : 0); x++) {
            if (M >> x & 1) continue;
            int k = popc(out[x] & M);
            c += (k < sz - k) ? k : sz - k;
        }
        if (sz == 2) {
            /* recompute exactly for twin (the loop above may have stopped early only when sz != 2) */
            if (c < bt) bt = c;
        }
        if (c < best) { best = c; bm = M; }
    }
    if (twin) *twin = bt;
    if (argM) *argM = bm;
    return best;
}

static int is_prime_enum(int n, const unsigned *out) {
    unsigned full = (1u << n) - 1;
    for (unsigned M = 1; M < full; M++) {
        int sz = popc(M);
        if (sz < 2) continue;
        int ok = 1;
        for (int x = 0; x < n && ok; x++) {
            if (M >> x & 1) continue;
            unsigned o = out[x] & M;
            if (o != 0 && o != M) ok = 0;
        }
        if (ok) return 0; /* M is a module */
    }
    return 1;
}

/* spectra bitmaps */
static unsigned char *Sbm[12];
static int Smax[12];
static void load_S(int n, const char *fn) {
    FILE *f = fopen(fn, "r");
    if (!f) { fprintf(stderr, "cannot open %s\n", fn); exit(1); }
    long h, mx = 0;
    long *vals = malloc(sizeof(long) * 200000);
    int cnt = 0;
    while (fscanf(f, "%ld", &h) == 1) { vals[cnt++] = h; if (h > mx) mx = h; }
    fclose(f);
    Smax[n] = (int)mx;
    Sbm[n] = calloc(mx + 2, 1);
    for (int i = 0; i < cnt; i++) Sbm[n][vals[i]] = 1;
    free(vals);
}
static int inS(int n, u64 h) { return h <= (u64)Smax[n] && Sbm[n][h]; }

static u64 rng = 88172645463325252ULL;
static u64 xr(void) { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return rng; }

static u64 brute_paths(int n, const unsigned *out) {
    int perm[MAXN];
    for (int i = 0; i < n; i++) perm[i] = i;
    u64 cnt = 0;
    /* Heap's algorithm over all n! permutations */
    int c[MAXN] = {0};
    int ok = 1;
    for (int i = 0; i + 1 < n; i++) if (!(out[perm[i]] >> perm[i + 1] & 1)) { ok = 0; break; }
    cnt += ok;
    int i = 0;
    while (i < n) {
        if (c[i] < i) {
            if (i % 2 == 0) { int t = perm[0]; perm[0] = perm[i]; perm[i] = t; }
            else { int t = perm[c[i]]; perm[c[i]] = perm[i]; perm[i] = t; }
            ok = 1;
            for (int j = 0; j + 1 < n; j++) if (!(out[perm[j]] >> perm[j + 1] & 1)) { ok = 0; break; }
            cnt += ok;
            c[i]++;
            i = 0;
        } else { c[i] = 0; i++; }
    }
    return cnt;
}

static void rand_tour(int n, unsigned *out) {
    for (int i = 0; i < n; i++) out[i] = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++) {
            if (xr() & 1) out[i] |= 1u << j; else out[j] |= 1u << i;
        }
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage\n"); return 1; }
    char line[4096];
    unsigned out[MAXN];
    if (!strcmp(argv[1], "selftest")) {
        int bad = 0, tests = 0;
        for (int n = 1; n <= 8; n++)
            for (int t = 0; t < 300; t++) {
                rand_tour(n, out);
                u64 a = ham_sub(n, out, (1u << n) - 1), b = brute_paths(n, out);
                tests++;
                if (a != b) { bad++; fprintf(stderr, "mismatch n=%d %llu %llu\n", n, a, b); }
                /* cross-check primality enumeration vs moddef==0 */
                int tw; unsigned M;
                int d = moddef(n, out, &tw, &M);
                if (n >= 3 && ((d == 0) != !is_prime_enum(n, out))) { bad++; fprintf(stderr, "prime mismatch\n"); }
            }
        printf("selftest: %d tests, %d failures\n", tests, bad);
        return bad != 0;
    }
    if (!strcmp(argv[1], "band")) {
        int s = atoi(argv[2]), n = atoi(argv[3]);
        for (int i = 0; i < n; i++) out[i] = 0;
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++) {
                if (j - i <= s) out[j] |= 1u << i; else out[i] |= 1u << j;
            }
        char code[512];
        encode(n, out, code);
        printf("%d %s\n", n, code);
        return 0;
    }
    if (!strcmp(argv[1], "props")) {
        while (fgets(line, sizeof line, stdin)) {
            int n; char code[512];
            if (sscanf(line, "%d %511s", &n, code) != 2) continue;
            if ((int)strlen(code) != n * (n - 1) / 2) { fprintf(stderr, "bad length\n"); return 1; }
            decode(code, n, out);
            unsigned full = (1u << n) - 1;
            u64 H = ham_sub(n, out, full);
            int lam = lambda_cut(n, out);
            int tw; unsigned M;
            int d = moddef(n, out, &tw, &M);
            int pr = is_prime_enum(n, out);
            printf("%d %llu lam=%d prime=%d moddef=%d twin=%d |", n, H, lam, pr, d, tw);
            for (int v = 0; v < n; v++) printf(" %llu", ham_sub(n, out, full & ~(1u << v)));
            printf(" | sc");
            for (int v = 0; v < n; v++) printf(" %d", popc(out[v]));
            printf(" | M %x\n", M);
            fflush(stdout);
        }
        return 0;
    }
    if (!strcmp(argv[1], "subs")) {
        int K = atoi(argv[2]);
        while (fgets(line, sizeof line, stdin)) {
            int n; char code[512];
            if (sscanf(line, "%d %511s", &n, code) != 2) continue;
            decode(code, n, out);
            unsigned full = (1u << n) - 1;
            printf("%d", n);
            for (unsigned D = 0; D <= full; D++) {
                if (popc(D) != K) continue;
                printf(" %llu", ham_sub(n, out, full & ~D));
            }
            printf("\n");
            fflush(stdout);
        }
        return 0;
    }
    if (!strcmp(argv[1], "ext")) {
        /* parent: S17 integer code, bit b = pair b (row-major i<j), 1 iff i -> j */
        u64 pc = strtoull(argv[2], NULL, 10);
        load_S(10, argv[3]);
        load_S(11, argv[4]);
        int n0 = 11, n = 12;
        unsigned par[MAXN];
        for (int i = 0; i < n0; i++) par[i] = 0;
        int b = 0;
        for (int i = 0; i < n0; i++)
            for (int j = i + 1; j < n0; j++) {
                if (pc >> b & 1) par[i] |= 1u << j; else par[j] |= 1u << i;
                b++;
            }
        u64 Hpar = ham_sub(n0, par, (1u << n0) - 1);
        printf("# parent %llu H=%llu lam=%d prime=%d\n", pc, Hpar, lambda_cut(n0, par), is_prime_enum(n0, par));
        long stat[32] = {0};
        int nesc = 0;
        for (unsigned x = 0; x < (1u << n0); x++) {
            for (int i = 0; i < n0; i++) out[i] = par[i];
            out[11] = 0;
            for (int a = 0; a < n0; a++) {
                if (x >> a & 1) out[a] |= 1u << 11; else out[11] |= 1u << a;
            }
            unsigned full = (1u << n) - 1;
            u64 H = ham_sub(n, out, full);
            int FC = !inS(11, H);
            int FCV = 0;
            for (int v = 0; v < n; v++) {
                u64 h = ham_sub(n, out, full & ~(1u << v));
                if (!inS(10, h)) FCV = 1;
            }
            int lam = lambda_cut(n, out);
            int SUM = (lam == 0);
            int SUB = !is_prime_enum(n, out);
            int SUM2 = (lam <= 2);
            int key = FC | FCV << 1 | SUM << 2 | SUB << 3 | SUM2 << 4;
            stat[key]++;
            if (!key) {
                char code[512];
                encode(n, out, code);
                printf("E %s H=%llu lam=%d x=%u\n", code, H, lam, x);
                nesc++;
            }
        }
        const char *nm[5] = {"FC", "FCV", "SUM", "SUB", "SUM2"};
        for (int k = 0; k < 32; k++) {
            if (!stat[k]) continue;
            printf("# labels");
            if (!k) printf(" ESCAPER");
            for (int t = 0; t < 5; t++) if (k >> t & 1) printf(" %s", nm[t]);
            printf(" : %ld\n", stat[k]);
        }
        printf("# escapers %d\n", nesc);
        return 0;
    }
    fprintf(stderr, "unknown mode\n");
    return 1;
}
