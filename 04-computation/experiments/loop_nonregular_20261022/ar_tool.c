/* Almost-regular tournaments (even n, scores n/2 and n/2 - 1): MCMC sampling by 3-cycle reversals (score preserving,
 * connected on a score class), exact Hamiltonian path counts by Held-Karp, and isomorphism classes by a canonical form
 * (colour refinement on out/in-neighbour colour counts plus individualisation-refinement; the canonical code is the
 * minimum upper-triangle code over all leaves of the search tree, and |Aut| is the number of leaves attaining it).
 * Session opus-2026-10-22-S54 (replaces gentourng, which is not installed on this machine).
 *
 * Usage:  ar_tool classes n target seed     -> sample until `target` distinct classes are seen (or 2e7 samples);
 *                                              print "code H aut" per class
 *         ar_tool sample n count seed       -> print "code H" for `count` samples (labelled, uniform in the limit)
 * code = upper-triangle bits in row-major order, '1' if i -> j (gentourng's tournament convention for i < j).      */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXN 24
static int n;
static int A[MAXN][MAXN];

static uint64_t rng_s[2];
static uint64_t rnd(void) {             /* xorshift128+ */
    uint64_t s1 = rng_s[0], s0 = rng_s[1];
    rng_s[0] = s0; s1 ^= s1 << 23; rng_s[1] = s1 ^ s0 ^ (s1 >> 17) ^ (s0 >> 26);
    return rng_s[1] + s0;
}

static void start(void) {               /* delete vertex n of the rotational (n+1)-tournament */
    int m = n + 1, d = n / 2;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            int t = ((j - i) % m + m) % m;
            A[i][j] = (i != j && t >= 1 && t <= d) ? 1 : 0;
        }
}

static void mcmc(long steps) {
    for (long s = 0; s < steps; s++) {
        int u = rnd() % n, v = rnd() % n, w = rnd() % n;
        if (u == v || v == w || u == w) continue;
        if (A[u][v] && A[v][w] && A[w][u]) {
            A[u][v] = A[v][w] = A[w][u] = 0; A[v][u] = A[w][v] = A[u][w] = 1;
        }
    }
}

static uint64_t heldkarp(void) {
    int N = 1 << n;
    uint64_t *dp = calloc((size_t)N * n, sizeof(uint64_t));
    for (int v = 0; v < n; v++) dp[((size_t)1 << v) * n + v] = 1;
    for (int mask = 1; mask < N; mask++)
        for (int v = 0; v < n; v++) {
            uint64_t c = dp[(size_t)mask * n + v];
            if (!c || !(mask >> v & 1)) continue;
            for (int w = 0; w < n; w++)
                if (!(mask >> w & 1) && A[v][w]) dp[(size_t)(mask | (1 << w)) * n + w] += c;
        }
    uint64_t tot = 0;
    for (int v = 0; v < n; v++) tot += dp[(size_t)(N - 1) * n + v];
    free(dp);
    return tot;
}

/* ---------- canonical form ---------- */
static int best[MAXN * MAXN], have_best, aut_count;

/* refine: colours col[v] in 0..k-1 (ordered partition); split by (col, out-counts per colour, in-counts per colour)
 * until stable; colours are re-indexed by sorting signatures, which keeps the procedure isomorphism invariant. */
static int refine(int *col) {
    int k = 0;
    for (int v = 0; v < n; v++) if (col[v] + 1 > k) k = col[v] + 1;
    for (;;) {
        static int sig[MAXN][2 * MAXN + 1];
        for (int v = 0; v < n; v++) {
            for (int c = 0; c < 2 * k + 1; c++) sig[v][c] = 0;
            sig[v][0] = col[v];
            for (int w = 0; w < n; w++) if (w != v) {
                if (A[v][w]) sig[v][1 + col[w]]++; else sig[v][1 + k + col[w]]++;
            }
        }
        int order[MAXN];
        for (int v = 0; v < n; v++) order[v] = v;
        /* insertion sort by signature */
        for (int i = 1; i < n; i++) {
            int x = order[i], j = i - 1;
            while (j >= 0 && memcmp(sig[order[j]], sig[x], sizeof(int) * (2 * k + 1)) > 0) { order[j + 1] = order[j]; j--; }
            order[j + 1] = x;
        }
        int newcol[MAXN], nk = 0;
        for (int i = 0; i < n; i++) {
            if (i > 0 && memcmp(sig[order[i]], sig[order[i - 1]], sizeof(int) * (2 * k + 1)) != 0) nk++;
            newcol[order[i]] = nk;
        }
        nk++;
        memcpy(col, newcol, sizeof(int) * n);
        if (nk == k) return k;
        k = nk;
    }
}

static void search(int *col) {
    int k = refine(col);
    if (k == n) {                                   /* leaf: position of v is col[v] */
        int inv[MAXN], code[MAXN * MAXN], c = 0;
        for (int v = 0; v < n; v++) inv[col[v]] = v;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) code[c++] = A[inv[i]][inv[j]];
        if (!have_best) { memcpy(best, code, sizeof(int) * c); have_best = 1; aut_count = 1; return; }
        int cmp = memcmp(code, best, sizeof(int) * c);
        if (cmp < 0) { memcpy(best, code, sizeof(int) * c); aut_count = 1; }
        else if (cmp == 0) aut_count++;
        return;
    }
    /* first non-singleton cell (smallest colour with >= 2 members) */
    int cnt[MAXN] = {0};
    for (int v = 0; v < n; v++) cnt[col[v]]++;
    int target = -1;
    for (int c = 0; c < k; c++) if (cnt[c] >= 2) { target = c; break; }
    for (int v = 0; v < n; v++) {
        if (col[v] != target) continue;
        int nc[MAXN];
        for (int w = 0; w < n; w++) nc[w] = (col[w] > target || (col[w] == target && w != v)) ? col[w] + 1 : col[w];
        search(nc);                                 /* v individualised: it keeps colour target, the rest of the cell moves up */
    }
}

static void canon(char *out, int *aut) {
    int col[MAXN];
    for (int v = 0; v < n; v++) { int s = 0; for (int w = 0; w < n; w++) s += A[v][w]; col[v] = s; }
    /* compress scores to 0..k-1 preserving order */
    int used[MAXN] = {0};
    for (int v = 0; v < n; v++) used[col[v]] = 1;
    int map[MAXN], k = 0;
    for (int s = 0; s < MAXN; s++) if (used[s]) map[s] = k++;
    for (int v = 0; v < n; v++) col[v] = map[col[v]];
    have_best = 0; aut_count = 0;
    search(col);
    int c = 0;
    for (int i = 0; i < n * (n - 1) / 2; i++) out[c++] = best[i] ? '1' : '0';
    out[c] = 0;
    *aut = aut_count;
}

/* simple open-addressing hash set of codes */
typedef struct { char *code; uint64_t H; int aut; } Entry;
static Entry *tab; static size_t tabsz, tabcnt;
static uint64_t hstr(const char *s) { uint64_t h = 1469598103934665603ULL; while (*s) { h ^= (unsigned char)*s++; h *= 1099511628211ULL; } return h; }
static int insert(const char *code, int aut) {
    size_t i = hstr(code) & (tabsz - 1);
    while (tab[i].code) { if (!strcmp(tab[i].code, code)) return 0; i = (i + 1) & (tabsz - 1); }
    tab[i].code = strdup(code); tab[i].aut = aut; tab[i].H = 0; tabcnt++;
    return 1;
}

/* forest numbers: c_L = (n - L)! N_L from sum over orderings of prod (1 + z S_(s_i s_(i+1))), S = A - A^T */
static void forest(void) {
    int N = 1 << n;
    long long *dp = calloc((size_t)N * n * n, sizeof(long long));
    for (int v = 0; v < n; v++) dp[((size_t)(1 << v) * n + v) * n + 0] = 1;
    for (int mask = 1; mask < N; mask++)
        for (int v = 0; v < n; v++) {
            if (!(mask >> v & 1)) continue;
            long long *src = dp + ((size_t)mask * n + v) * n;
            for (int w = 0; w < n; w++) {
                if (mask >> w & 1) continue;
                int s = A[v][w] ? 1 : -1;
                long long *dst = dp + ((size_t)(mask | (1 << w)) * n + w) * n;
                for (int L = n - 2; L >= 0; L--) { if (src[L]) { dst[L] += src[L]; dst[L + 1] += s * src[L]; } }
            }
        }
    for (int L = 0; L < n; L++) {
        long long c = 0;
        for (int v = 0; v < n; v++) c += dp[((size_t)(N - 1) * n + v) * n + L];
        printf("%s%lld", L ? " " : "", c);
    }
    putchar(10);
    free(dp);
}

int main(int argc, char **argv) {
    if (argc >= 3 && !strcmp(argv[1], "forest")) {      /* ar_tool forest n < codes : prints c_0 .. c_(n-1) per code */
        n = atoi(argv[2]);
        if (n > 16) { fprintf(stderr, "forest: n <= 16 (memory n^2 2^n)\n"); return 1; }
        char buf[4096];
        while (fgets(buf, sizeof buf, stdin)) {
            char code[MAXN * MAXN]; if (sscanf(buf, "%s", code) != 1) continue;
            int c = 0;
            for (int a = 0; a < n; a++) { A[a][a] = 0; for (int b = a + 1; b < n; b++) { A[a][b] = code[c++] == '1'; A[b][a] = !A[a][b]; } }
            forest();
        }
        return 0;
    }
    if (argc < 5) { fprintf(stderr, "usage: ar_tool classes|sample n target|count seed\n"); return 1; }
    n = atoi(argv[2]); long target = atol(argv[3]); uint64_t seed = strtoull(argv[4], 0, 10);
    rng_s[0] = seed * 0x9E3779B97F4A7C15ULL + 1; rng_s[1] = (seed + 17) * 0xBF58476D1CE4E5B9ULL + 3;
    start(); mcmc(200L * n * n * n);
    long thin = 4L * n * n;
    if (!strcmp(argv[1], "random")) {                     /* uniformly random tournaments (every arc a fair coin) */
        char code[MAXN * MAXN];
        for (long s = 0; s < target; s++) {
            for (int i = 0; i < n; i++) { A[i][i] = 0; for (int j = i + 1; j < n; j++) { A[i][j] = rnd() >> 63; A[j][i] = !A[i][j]; } }
            int c = 0;
            for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) code[c++] = A[i][j] ? '1' : '0';
            code[c] = 0;
            printf("%s %llu\n", code, (unsigned long long)heldkarp());
        }
        return 0;
    }
    if (!strcmp(argv[1], "sample")) {
        char code[MAXN * MAXN];
        for (long s = 0; s < target; s++) {
            mcmc(thin);
            int c = 0;
            for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) code[c++] = A[i][j] ? '1' : '0';
            code[c] = 0;
            printf("%s %llu\n", code, (unsigned long long)heldkarp());
            fflush(stdout);
        }
        return 0;
    }
    tabsz = 1 << 22; tab = calloc(tabsz, sizeof(Entry));
    char code[MAXN * MAXN]; int aut;
    long samples = 0, lastnew = 0;
    while ((long)tabcnt < target && samples < 20000000L) {
        mcmc(thin); samples++;
        canon(code, &aut);
        if (insert(code, aut)) lastnew = samples;
    }
    fprintf(stderr, "samples %ld, classes %zu, last new class at sample %ld\n", samples, tabcnt, lastnew);
    for (size_t i = 0; i < tabsz; i++) if (tab[i].code) {
        /* rebuild A from the canonical code and count paths */
        int c = 0;
        for (int a = 0; a < n; a++) { A[a][a] = 0; for (int b = a + 1; b < n; b++) { A[a][b] = tab[i].code[c++] == '1'; A[b][a] = !A[a][b]; } }
        printf("%s %llu %d\n", tab[i].code, (unsigned long long)heldkarp(), tab[i].aut);
    }
    return 0;
}
