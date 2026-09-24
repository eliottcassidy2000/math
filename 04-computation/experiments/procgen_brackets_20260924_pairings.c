/* procgen_brackets_20260924_pairings.c -- the "pairing ladder" of the Collatz map.
 *
 * Pairing maps.  offset 0: pairs {2i-1, 2i}, length i = ceil(n/2);  offset 1: pairs {2i, 2i+1}, i = floor(n/2).
 * A choice bit eps_i per pair; n moves UP (n -> n+i) iff  (n odd) XOR eps_i,  otherwise DOWN (n -> n-i).
 * eps == 0: offset 0 is the Collatz shortcut T ((3n+1)/2, n/2); offset 1 is the 3n-1 map ((3n-1)/2, n/2).
 * Every member preserves the sum of each pair and uses every length i once up and once down.
 *
 * Modes (all output is plain text, parsed by procgen_brackets_20260924_pairings.py):
 *   periodic K OFF N CAP M     every c : Z/2^K -> {0,1} (eps_i = c(i mod 2^K)); orbits of n <= N;
 *                              cycles, escapes (> CAP), exceptional classes at 2-adic level M,
 *                              landing=>down flag, Haar-measure-preservation flag
 *   random SEED P OFF N CAP    eps_i = [hash(SEED,i) < P] (P = flip probability vs the eps=0 map)
 *   single OFF I N             flip exactly one pair i (1 <= i <= I) of the eps=0 map; report new cycles
 *   landing I N                the landing=>down member closest to Collatz (free bits 0), density and check
 *   landingvc X                min vertex cover of the {2k,3k} constraint forest on [1,X] (lower bound)
 *   pair2 N                    3n-1 map with two flips (one per nontrivial cycle): which make it a tree to N
 *   cert MASK K OFF M          exceptional-class DFS + threshold for one periodic map (proof mode)
 *   qmap Q SGN N CAP           cycle census of n -> (Qn+SGN)/2 (odd), n/2 (even) for comparison
 * Memory: O(N) bytes; N = 10^6 uses ~20 MB.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

typedef unsigned __int128 u128;
enum { K_PERIODIC, K_RANDOM, K_ARRAY, K_FLIPS, K_QMAP };
typedef struct {
    int offset, kind, K;
    uint64_t mask;           /* periodic: bit (i mod 2^K) */
    uint64_t seed, thresh;   /* random */
    uint8_t *arr; uint64_t arrlen; int arrdefault;
    int64_t flips[4]; int nflips;
    uint64_t q; int sgn;     /* K_QMAP: n odd -> (q n + sgn)/2, n even -> n/2 (not a pairing map; used for comparison) */
} Pairing;

static inline uint64_t splitmix(uint64_t x) {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}
static inline int eps_of(const Pairing *P, uint64_t i) {
    switch (P->kind) {
    case K_PERIODIC: return (int)((P->mask >> (i & ((1ULL << P->K) - 1))) & 1ULL);
    case K_RANDOM: return splitmix(P->seed * 0x100000001B3ULL ^ splitmix(i)) < P->thresh;
    case K_ARRAY: return i < P->arrlen ? P->arr[i] : P->arrdefault;
    case K_FLIPS: for (int t = 0; t < P->nflips; t++) if ((int64_t)i == P->flips[t]) return 1; return 0;
    }
    return 0;
}
static inline uint64_t F(const Pairing *P, uint64_t n) {
    if (P->kind == K_QMAP) { if (!(n & 1)) return n >> 1; if (n > (UINT64_MAX >> 1) / P->q) return UINT64_MAX; return (P->q * n + (uint64_t)(int64_t)P->sgn) >> 1; }
    uint64_t i = P->offset == 0 ? (n + 1) >> 1 : n >> 1;
    int up = (int)(n & 1) ^ eps_of(P, i);
    return up ? n + i : n - i;
}

/* ---------------------------------------------------------------- orbit / cycle census */
typedef struct { uint64_t minv, maxv; int64_t len; int64_t basin; } Cycle;
#define MAXCYC 4096
static Cycle cyc[MAXCYC]; static int ncyc;
static uint8_t *status;      /* 0 unknown, 1 done, 2 on current path */
static int32_t *cid;         /* cycle id (>=0), -1 escape */
static uint64_t *stk; static int64_t stkcap;
/* open-addressing hash for large values on the current walk */
#define HB 16
static uint64_t hkey[1 << HB]; static int64_t hpos[1 << HB]; static int hused[1 << HB]; static int hlist[1 << HB]; static int nh;
static void hclear(void) { for (int t = 0; t < nh; t++) hused[hlist[t]] = 0; nh = 0; }
static int64_t hget(uint64_t v) { uint64_t h = splitmix(v) & ((1 << HB) - 1); while (hused[h]) { if (hkey[h] == v) return hpos[h]; h = (h + 1) & ((1 << HB) - 1); } return -1; }
static int hput(uint64_t v, int64_t p) { if (nh > (1 << (HB - 1))) return 0; uint64_t h = splitmix(v) & ((1 << HB) - 1); while (hused[h]) h = (h + 1) & ((1 << HB) - 1); hused[h] = 1; hkey[h] = v; hpos[h] = p; hlist[nh++] = (int)h; return 1; }

static int64_t nesc;
static void census(const Pairing *P, uint64_t N, uint64_t CAP) {
    memset(status, 0, N + 1); ncyc = 0; nesc = 0;
    for (uint64_t n0 = 0; n0 <= N; n0++) {
        if (status[n0]) continue;
        int64_t sp = 0; uint64_t v = n0; int outcome = -2; int32_t id = -1;
        hclear();
        for (;;) {
            if (v <= N) {
                if (status[v] == 1) { id = cid[v]; outcome = 0; break; }
                if (status[v] == 2) { outcome = 1; break; }
                status[v] = 2;
            } else {
                if (v > CAP) { outcome = 2; break; }
                if (hget(v) >= 0) { outcome = 1; break; }
                if (!hput(v, sp)) { outcome = 2; break; }   /* walk too long above N: treat as escape */
            }
            if (sp >= stkcap) { outcome = 2; break; }
            stk[sp++] = v; v = F(P, v);
        }
        if (outcome == 1) {                 /* new cycle through v: find v on the stack */
            int64_t p = sp - 1; while (stk[p] != v) p--;
            Cycle c; c.minv = v; c.maxv = v; c.len = sp - p; c.basin = 0;
            for (int64_t q = p; q < sp; q++) { if (stk[q] < c.minv) c.minv = stk[q]; if (stk[q] > c.maxv) c.maxv = stk[q]; }
            if (ncyc < MAXCYC) { cyc[ncyc] = c; id = ncyc++; } else id = MAXCYC - 1;
        } else if (outcome == 2) { id = -1; nesc++; }
        for (int64_t q = 0; q < sp; q++) if (stk[q] <= N) { status[stk[q]] = 1; cid[stk[q]] = id; }
    }
    for (uint64_t n = 1; n <= N; n++) if (cid[n] >= 0 && cid[n] < ncyc) cyc[cid[n]].basin++;
}
static void alloc_census(uint64_t N) {
    status = malloc(N + 2); cid = malloc((N + 2) * sizeof(int32_t));
    stkcap = 1 << 22; stk = malloc(stkcap * sizeof(uint64_t));
    if (!status || !cid || !stk) { fprintf(stderr, "alloc\n"); exit(1); }
}
static void print_cycles(const char *tag) {
    /* sort by min */
    for (int a = 0; a < ncyc; a++) for (int b = a + 1; b < ncyc; b++) if (cyc[b].minv < cyc[a].minv) { Cycle t = cyc[a]; cyc[a] = cyc[b]; cyc[b] = t; }
    printf("%s ncyc %d nesc %lld cycles", tag, ncyc, (long long)nesc);
    for (int a = 0; a < ncyc && a < 40; a++) printf(" %llu:%lld:%llu:%lld", (unsigned long long)cyc[a].minv, (long long)cyc[a].len, (unsigned long long)cyc[a].maxv, (long long)cyc[a].basin);
    printf("\n");
}

/* ---------------------------------------------------------------- 2-adic exceptional classes (periodic maps)
 * class x = r (mod 2^j).  Simulate the integer representative r: F^s(x) = F^s(r) + 3^a 2^(j-s) t, so the branch of
 * step s+1 is known while j - s >= K+1 (it needs F^s(x) mod 2^(K+1)).  Certificate: 3^a < 2^s at some known step.
 * threshold for the certified class: x > B/(2^s - 3^a) with B = 2^s F^s(r) - 3^a r.
 */
static int64_t exc_count, exc_cap = 4000000; static double thr_max; static int KK;
static const Pairing *PP;
static void dfs(u128 r, int j, int M) {
    if (exc_count > exc_cap) return;
    u128 v = r; int s = 0, a = 0; double lg = 0.0;
    while (j - s >= KK + 1) {
        uint64_t i = PP->offset == 0 ? (uint64_t)((v + 1) >> 1) : (uint64_t)(v >> 1);
        int up = (int)(v & 1) ^ eps_of(PP, i);
        v = up ? v + i : v - i; s++; if (up) { a++; lg += log2(3.0); }
        if (lg < (double)s - 1e-9) {        /* 3^a < 2^s: certified */
            /* threshold x0 = (2^s v - 3^a r)/(2^s - 3^a) evaluated in long double */
            long double p3 = powl(3.0L, a), p2 = powl(2.0L, s);
            long double B = p2 * (long double)v - p3 * (long double)r;
            long double x0 = B / (p2 - p3);
            if (x0 > thr_max) thr_max = (double)x0;
            return;
        }
    }
    if (j >= M) { exc_count++; return; }
    dfs(r, j + 1, M); dfs(r + ((u128)1 << j), j + 1, M);
}
static int64_t exceptional(const Pairing *P, int K, int M, double *thr) {
    PP = P; KK = K; exc_count = 0; thr_max = 0;
    for (uint64_t r = 0; r < (1ULL << (K + 1)); r++) dfs(r, K + 1, M);
    *thr = thr_max; return exc_count;
}
/* landing => down for a periodic map: for every i mod 2^(K+1): eps(j(i)) == [i + eps(i) even] */
static int landing_ok(const Pairing *P, int K) {
    for (uint64_t i = 2; i < 2 + (1ULL << (K + 2)); i++) {
        int e = eps_of(P, i);
        uint64_t L = P->offset == 0 ? 3 * i - 1 + e : 3 * i + 1 - e;   /* up-mover's landing */
        uint64_t j = P->offset == 0 ? (L + 1) >> 1 : L >> 1;
        int down = !((int)(L & 1) ^ eps_of(P, j));
        if (!down) return 0;
    }
    return 1;
}
/* Haar measure preservation: every class mod 2^K is hit exactly twice by the images of the 2^(K+1) classes */
static int measure_preserving(const Pairing *P, int K) {
    uint64_t Mo = 1ULL << K; int *cnt = calloc(Mo, sizeof(int));
    for (uint64_t r = (1ULL << (K + 1)); r < (1ULL << (K + 2)); r++) cnt[F(P, r) & (Mo - 1)]++;   /* r = representative of class r mod 2^(K+1), shifted to avoid n = 0,1 */
    int ok = 1; for (uint64_t c = 0; c < Mo; c++) if (cnt[c] != 2) ok = 0;
    free(cnt); return ok;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "mode?\n"); return 1; }
    const char *mode = argv[1];
    if (!strcmp(mode, "periodic")) {
        int K = atoi(argv[2]), off = atoi(argv[3]); uint64_t N = strtoull(argv[4], 0, 10), CAP = strtoull(argv[5], 0, 10); int M = atoi(argv[6]);
        alloc_census(N);
        uint64_t nmaps = 1ULL << (1ULL << K);
        for (uint64_t mask = 0; mask < nmaps; mask++) {
            Pairing P = {0}; P.offset = off; P.kind = K_PERIODIC; P.K = K; P.mask = mask;
            census(&P, N, CAP);
            double thr; int64_t ex = M > 0 ? exceptional(&P, K, M, &thr) : -1;
            int ones = __builtin_popcountll(mask);
            printf("MAP K %d off %d mask %llu ones %d mp %d land %d exc %lld thr %.3g ", K, off, (unsigned long long)mask, ones,
                   measure_preserving(&P, K), landing_ok(&P, K), (long long)ex, M > 0 ? thr : -1.0);
            print_cycles("|");
        }
        return 0;
    }
    if (!strcmp(mode, "cert")) {
        uint64_t mask = strtoull(argv[2], 0, 10); int K = atoi(argv[3]), off = atoi(argv[4]), M = atoi(argv[5]);
        Pairing P = {0}; P.offset = off; P.kind = K_PERIODIC; P.K = K; P.mask = mask;
        for (int m = K + 2; m <= M; m += 2) { double thr; int64_t ex = exceptional(&P, K, m, &thr); printf("CERT mask %llu K %d off %d level %d exc %lld%s thr %.6g\n", (unsigned long long)mask, K, off, m, (long long)ex, ex > exc_cap ? "+" : "", thr); if (ex == 0 || ex > exc_cap) break; }
        return 0;
    }
    if (!strcmp(mode, "random")) {
        uint64_t seed = strtoull(argv[2], 0, 10); double p = atof(argv[3]); int off = atoi(argv[4]);
        uint64_t N = strtoull(argv[5], 0, 10), CAP = strtoull(argv[6], 0, 10);
        alloc_census(N);
        Pairing P = {0}; P.offset = off; P.kind = K_RANDOM; P.seed = seed;
        P.thresh = p >= 1.0 ? UINT64_MAX : (uint64_t)(p * 18446744073709551616.0);
        census(&P, N, CAP);
        printf("RANDOM seed %llu p %g off %d N %llu ", (unsigned long long)seed, p, off, (unsigned long long)N);
        print_cycles("|");
        return 0;
    }
    if (!strcmp(mode, "single")) {
        /* flip exactly one pair i of the eps=0 map.  A new cycle must pass through the flipped pair {a, a+1};
           walk from each element with a visited-set; a repeat at the start = a cycle through it. */
        int off = atoi(argv[2]); uint64_t I = strtoull(argv[3], 0, 10); int64_t STEPS = atoll(argv[4]);
        Pairing P = {0}; P.offset = off; P.kind = K_FLIPS; P.nflips = 1;
        int64_t nfrag = 0, ncap = 0; int printed = 0; uint64_t rep = 1;
        for (uint64_t i = 1; i <= I; i++) {
            P.flips[0] = (int64_t)i;
            uint64_t a = off == 0 ? 2 * i - 1 : 2 * i;
            int found = 0; uint64_t cmin = 0, cmax = 0; int64_t clen = 0;
            for (int w = 0; w < 2 && !found; w++) {
                uint64_t start = a + w, v = start; int64_t st = 0; hclear();
                for (;;) {
                    if (hget(v) >= 0) { if (v == start) found = 1; break; }
                    if (!hput(v, st) || st > STEPS) { ncap++; break; }
                    if (v == 0 || (off == 1 && v == 1)) break;           /* fixed points */
                    v = F(&P, v); st++;
                }
                if (found) { uint64_t x = start; cmin = cmax = start; clen = 0; do { x = F(&P, x); clen++; if (x < cmin) cmin = x; if (x > cmax) cmax = x; } while (x != start); }
            }
            if (found) { nfrag++; if (printed < 80) { printf("FRAG i %llu pair %llu,%llu cycle_min %llu len %lld max %llu\n", (unsigned long long)i, (unsigned long long)a, (unsigned long long)(a + 1), (unsigned long long)cmin, (long long)clen, (unsigned long long)cmax); printed++; } }
            if (i == rep || i == I) { printf("SINGLE off %d upto %llu fragile %lld capped %lld\n", off, (unsigned long long)i, (long long)nfrag, (long long)ncap); rep *= 4; }
        }
        return 0;
    }
    if (!strcmp(mode, "landing")) {
        /* landing => down, free bits 0 (closest to Collatz), eps_1 = 0; recursion eps_3k = 1 - eps_2k,
           eps_(2k+1) = 0 => eps_(3k+1) = 0, eps_(2k+1) = 1 => eps_(3k+2) = 1. */
        uint64_t I = strtoull(argv[2], 0, 10), N = strtoull(argv[3], 0, 10);
        uint8_t *e = calloc(2 * I + 8, 1); uint8_t *forced = calloc(2 * I + 8, 1);
        int64_t ones = 0, conflicts = 0; uint64_t nextrep = 16;
        int64_t cnt_by_v3[40] = {0}, tot_by_v3[40] = {0};
        for (uint64_t i = 1; i <= I; i++) {
            if (!forced[i]) e[i] = 0;
            uint64_t L = 3 * i - 1 + e[i], j = (L + 1) >> 1;
            int need = (int)(L & 1);          /* L goes down iff eps_j == [L odd] */
            if (j <= 2 * I + 4) {
                if (j > i) { if (forced[j] && e[j] != need) conflicts++; e[j] = (uint8_t)need; forced[j] = 1; }
                else if (e[j] != need) conflicts++;   /* j == i only for i = 1 */
            }
            ones += e[i];
            int v3 = 0; uint64_t t = i; while (t % 3 == 0) { t /= 3; v3++; }
            tot_by_v3[v3]++; cnt_by_v3[v3] += e[i];
            if (i == nextrep || i == I) { printf("LANDING upto %llu ones %lld density %.6f conflicts %lld\n", (unsigned long long)i, (long long)ones, (double)ones / i, (long long)conflicts); nextrep *= 4; }
        }
        printf("LANDING_V3");
        for (int v = 0; v < 8; v++) printf(" v3=%d:%.4f", v, tot_by_v3[v] ? (double)cnt_by_v3[v] / tot_by_v3[v] : -1.0);
        printf("\n");
        printf("LANDING_PREFIX");
        for (uint64_t i = 1; i <= 60; i++) printf("%d", e[i]);
        printf("\n");
        /* tree check: all n <= N reach 1 */
        Pairing P = {0}; P.offset = 0; P.kind = K_ARRAY; P.arr = e; P.arrlen = 2 * I + 4; P.arrdefault = 0;
        alloc_census(N); census(&P, N, (uint64_t)1 << 62);
        print_cycles("LANDING_TREE |");
        /* descent within 2 steps for every n in [2, N] */
        int64_t bad = 0; for (uint64_t n = 3; n <= N; n++) { uint64_t v = F(&P, n); if (v >= n) { v = F(&P, v); if (v >= n) bad++; } }
        printf("LANDING_2STEP_DESCENT_FAILS %lld\n", (long long)bad);
        return 0;
    }
    if (!strcmp(mode, "landingvc")) {
        /* min vertex cover of the forest with edges {2k, 3k}, 3k <= X.  Components: for u coprime to 6 and s >= 0
           the anti-diagonal u 2^a 3^(s-a), a = 0..s, a path; its vertices <= X form a subpath (values decrease in a). */
        uint64_t X = strtoull(argv[2], 0, 10); int64_t vc = 0, v3odd = 0, nv = 0;
        for (uint64_t u = 1; u <= X; u += 1) {
            if (u % 2 == 0 || u % 3 == 0) continue;
            for (int s = 0; ; s++) {
                u128 v = (u128)u << s; if (v > X) break;       /* a = s: smallest vertex */
                int64_t t = 0;
                for (int a2 = s; a2 >= 0; a2--) { if (v > X) break; t++; if (a2 > 0) v = v / 2 * 3; }
                vc += t / 2; nv += t;
            }
        }
        for (uint64_t v = 1; v <= X; v++) { uint64_t t = v; int c = 0; while (t % 3 == 0) { t /= 3; c++; } v3odd += c & 1; }
        printf("LANDINGVC X %llu vertices %lld min_vertex_cover %lld density %.6f v3odd_density %.6f\n", (unsigned long long)X, (long long)nv, (long long)vc, (double)vc / X, (double)v3odd / X);
        return 0;
    }
    if (!strcmp(mode, "pair2")) {
        /* 3n-1 map (offset 1): flip one pair touching {5,7,10} and one touching {17,...}; tree to N? */
        uint64_t N = strtoull(argv[2], 0, 10);
        uint64_t c1[] = {5, 7, 10}; uint64_t c2[] = {17, 25, 37, 55, 82, 41, 61, 91, 136, 68, 34};
        alloc_census(N);
        for (int x = 0; x < 3; x++) for (int y = 0; y < 11; y++) {
            Pairing P = {0}; P.offset = 1; P.kind = K_FLIPS; P.nflips = 2; P.flips[0] = (int64_t)(c1[x] >> 1); P.flips[1] = (int64_t)(c2[y] >> 1);
            census(&P, N, (uint64_t)1 << 62);
            printf("PAIR2 flips %lld %lld ", (long long)P.flips[0], (long long)P.flips[1]);
            print_cycles("|");
        }
        /* also: all single flips i <= 200 of the 3n-1 map */
        for (int64_t i = 1; i <= 200; i++) {
            Pairing P = {0}; P.offset = 1; P.kind = K_FLIPS; P.nflips = 1; P.flips[0] = i;
            census(&P, N, (uint64_t)1 << 62);
            printf("MINUS1 flip %lld ", (long long)i); print_cycles("|");
        }
        return 0;
    }
    if (!strcmp(mode, "qmap")) {
        uint64_t q = strtoull(argv[2], 0, 10); int sgn = atoi(argv[3]); uint64_t N = strtoull(argv[4], 0, 10), CAP = strtoull(argv[5], 0, 10);
        alloc_census(N);
        Pairing P = {0}; P.kind = K_QMAP; P.q = q; P.sgn = sgn;
        census(&P, N, CAP);
        printf("QMAP q %llu sgn %d N %llu ", (unsigned long long)q, sgn, (unsigned long long)N);
        print_cycles("|");
        return 0;
    }
    fprintf(stderr, "unknown mode\n"); return 1;
}
