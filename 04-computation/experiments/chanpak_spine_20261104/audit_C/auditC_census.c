/* auditC_census.c -- independent audit (audit C, 2026-11-04) of THM-4876 / HYP-9469.
 *
 * Written from scratch for the audit; shares no code with the session's scripts.
 * Input: gentourng default output on stdin (upper triangle row-by-row, '1' at (i,j), i<j, means i -> j).
 *
 * Modes:
 *   del n          deletion census over all tournaments T of order n read from stdin: for every vertex v,
 *                  H(T), H(T-v), the increment k = (H(T)-H(T-v))/2; k_min per h = H(T-v); the h-sets of the
 *                  increment levels k <= KLEV; end counts e_a(T) and the bound e_a <= (H-1)/2 for non-sinks;
 *                  prints "P <tour> <v> <H> <H-v>" for every (T,v) with k in {1,2} (for structural checks in Python).
 *   ext n          extension spectra: for every X of order n on stdin and every cut A (bitmask; v -> A, V\A -> v),
 *                  H(X+v) by a direct subset DP on X+v (no path-pair formula).  Prints "E <tour> <h> <values...>".
 *   rand n r S seed   random search: S tournaments of order n, each TT_n with r random arcs reversed (r = 0 means a
 *                  uniformly random tournament); checks: increment 2 => H(T-v) in {1,3} and H(T) in {3,5};
 *                  increment 4 => H(T-v) in {1,5}; e_a <= (H-1)/2 for non-sinks.  Prints violations.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXN 13
#define KLEV 12
#define HMAX 70000

static int n;
static uint32_t adj[MAXN];           /* out-neighbour masks */
static uint64_t *dp;                 /* dp[S*n + u] = Hamiltonian paths of T[S] ending at u */

/* dp over all subsets of [m] for out-masks a[]; returns nothing; fills dp (size (1<<m)*m). */
static void run_dp(int m, const uint32_t *a) {
    uint32_t in[MAXN];
    int u, w;
    for (u = 0; u < m; u++) { in[u] = 0; for (w = 0; w < m; w++) if ((a[w] >> u) & 1) in[u] |= 1u << w; }
    uint32_t S, N = 1u << m;
    memset(dp, 0, sizeof(uint64_t) * (size_t)N * m);
    for (S = 1; S < N; S++) {
        for (u = 0; u < m; u++) {
            if (!((S >> u) & 1)) continue;
            uint32_t R = S & ~(1u << u);
            if (R == 0) { dp[(size_t)S * m + u] = 1; continue; }
            uint32_t pred = R & in[u];
            uint64_t acc = 0;
            while (pred) { w = __builtin_ctz(pred); pred &= pred - 1; acc += dp[(size_t)R * m + w]; }
            dp[(size_t)S * m + u] = acc;
        }
    }
}
static uint64_t hsub(int m, uint32_t S) { uint64_t t = 0; int u; if (S == 0) return 1; for (u = 0; u < m; u++) if ((S >> u) & 1) t += dp[(size_t)S * m + u]; return t; }

static int parse(const char *line, int m) {
    int i, j, p = 0, need = m * (m - 1) / 2;
    if ((int)strlen(line) < need) return 0;
    for (i = 0; i < m; i++) adj[i] = 0;
    for (i = 0; i < m; i++) for (j = i + 1; j < m; j++) {
        char c = line[p++];
        if (c == '1') adj[i] |= 1u << j; else if (c == '0') adj[j] |= 1u << i; else return 0;
    }
    return 1;
}
static int is_strong(int m, const uint32_t *a) {
    uint32_t full = (1u << m) - 1, reach, fr; int u;
    /* forward from 0 */
    reach = 1; fr = 1;
    while (fr) { u = __builtin_ctz(fr); fr &= fr - 1; uint32_t nw = a[u] & ~reach; reach |= nw; fr |= nw; }
    if (reach != full) return 0;
    /* backward to 0 */
    uint32_t in[MAXN]; int w;
    for (u = 0; u < m; u++) { in[u] = 0; for (w = 0; w < m; w++) if ((a[w] >> u) & 1) in[u] |= 1u << w; }
    reach = 1; fr = 1;
    while (fr) { u = __builtin_ctz(fr); fr &= fr - 1; uint32_t nw = in[u] & ~reach; reach |= nw; fr |= nw; }
    return reach == full;
}

static char line[1 << 16];
static int cmp_u64(const void *x, const void *y) { uint64_t p = *(const uint64_t *)x, q = *(const uint64_t *)y; return p < q ? -1 : p > q; }

static void mode_del(void) {
    static int kmin[HMAX]; static unsigned char lev[KLEV + 1][HMAX];
    static long cnt_pairs_k[KLEV + 2];
    long ntour = 0, ntrans = 0, nstrong = 0, viol_inc2 = 0, viol_inc4 = 0, viol_end = 0, npairs = 0;
    uint64_t bestE = 0, bestH = 1; long nonsink_checked = 0;
    int h;
    for (h = 0; h < HMAX; h++) kmin[h] = -1;
    uint32_t full = (1u << n) - 1;
    while (fgets(line, sizeof line, stdin)) {
        if (!parse(line, n)) continue;
        ntour++;
        run_dp(n, adj);
        uint64_t H = hsub(n, full);
        if (H == 1) ntrans++;
        if (is_strong(n, adj)) nstrong++;
        int v;
        for (v = 0; v < n; v++) {
            uint64_t Hv = hsub(n, full & ~(1u << v));
            if (H < Hv) { printf("X MONOTONICITY VIOLATION %s v=%d\n", line, v); }
            uint64_t k2 = H - Hv;
            if (k2 & 1) { printf("X ODD INCREMENT %s v=%d\n", line, v); }
            long k = (long)(k2 / 2);
            if (k >= 1) {
                npairs++;
                if (Hv < HMAX && (kmin[Hv] < 0 || k < kmin[Hv])) kmin[Hv] = (int)k;
                if (k <= KLEV && Hv < HMAX) { lev[k][Hv] = 1; cnt_pairs_k[k]++; }
                if (k <= 2) {
                    line[strcspn(line, "\r\n")] = 0;
                    printf("P %s %d %llu %llu\n", line, v, (unsigned long long)H, (unsigned long long)Hv);
                    if (k == 1 && !(Hv == 1 || Hv == 3)) viol_inc2++;
                    if (k == 1 && !(H == 3 || H == 5)) viol_inc2++;
                    if (k == 2 && !(Hv == 1 || Hv == 5)) viol_inc4++;
                }
            }
            /* end-count bound for non-sinks (THM-4873 (i) sanity + HYP-9469 (4)) */
            if (adj[v] != 0 && H >= 3) {
                uint64_t e = dp[(size_t)full * n + v];
                nonsink_checked++;
                if (2 * e + 1 > H) { viol_end++; line[strcspn(line, "\r\n")] = 0; printf("V ENDCOUNT %s a=%d H=%llu e=%llu\n", line, v, (unsigned long long)H, (unsigned long long)e); }
                /* track max e/H */
                if ((__uint128_t)e * bestH > (__uint128_t)bestE * H) { bestE = e; bestH = H; }
            }
        }
    }
    printf("SUMMARY del n=%d tournaments=%ld transitive=%ld strong=%ld pairs_with_positive_increment=%ld\n", n, ntour, ntrans, nstrong, npairs);
    printf("SUMMARY violations: inc2_not_{1,3}or_H_not_{3,5}=%ld inc4_not_{1,5}=%ld endcount(e_a>(H-1)/2, non-sink)=%ld (non-sink checks %ld)\n", viol_inc2, viol_inc4, viol_end, nonsink_checked);
    printf("SUMMARY endcount max e_a/H over non-sinks (H>=3): %llu/%llu\n", (unsigned long long)bestE, (unsigned long long)bestH);
    printf("KMIN"); for (h = 0; h < HMAX; h++) if (kmin[h] >= 0) printf(" %d:%d", h, kmin[h]); printf("\n");
    int k;
    for (k = 1; k <= KLEV; k++) { printf("LEVEL k=%d (increment %d) pairs=%ld h-set:", k, 2 * k, cnt_pairs_k[k]); for (h = 0; h < HMAX; h++) if (lev[k][h]) printf(" %d", h); printf("\n"); }
}

static void mode_ext(void) {
    long ntour = 0;
    uint32_t fullx = (1u << n) - 1;
    uint32_t a2[MAXN];
    static uint64_t vals[1 << MAXN];
    while (fgets(line, sizeof line, stdin)) {
        if (!parse(line, n)) continue;
        ntour++;
        uint32_t A; int u;
        for (A = 0; A <= fullx; A++) {
            for (u = 0; u < n; u++) a2[u] = adj[u] | (((A >> u) & 1) ? 0 : (1u << n));
            a2[n] = A;
            run_dp(n + 1, a2);
            vals[A] = hsub(n + 1, (1u << (n + 1)) - 1);
        }
        /* sort distinct */
        uint32_t m = fullx + 1, i, j;
        qsort(vals, m, sizeof(uint64_t), cmp_u64);
        line[strcspn(line, "\r\n")] = 0;
        printf("E %s %llu", line, (unsigned long long)vals[0]);
        /* h = H(X) must equal the value at cut 0 (v a sink) = min */
        for (i = 0, j = 0; i < m; i++) if (i == 0 || vals[i] != vals[i - 1]) { printf(" %llu", (unsigned long long)vals[i]); j++; }
        printf("\n");
    }
    fprintf(stderr, "ext n=%d: %ld tournaments\n", n, ntour);
}

static uint64_t rng_state;
static uint64_t rng(void) { rng_state ^= rng_state << 13; rng_state ^= rng_state >> 7; rng_state ^= rng_state << 17; return rng_state; }

static void mode_rand(int r, long S, uint64_t seed) {
    long s, viol = 0, inc2 = 0, inc4 = 0, endv = 0; uint32_t full = (1u << n) - 1;
    uint64_t bestE = 0, bestH = 1;
    rng_state = seed * 2654435761ULL + 88172645463325252ULL;
    for (s = 0; s < S; s++) {
        int i, j;
        if (r == 0) {
            for (i = 0; i < n; i++) adj[i] = 0;
            for (i = 0; i < n; i++) for (j = i + 1; j < n; j++) { if (rng() & 1) adj[i] |= 1u << j; else adj[j] |= 1u << i; }
        } else {
            for (i = 0; i < n; i++) { adj[i] = 0; for (j = i + 1; j < n; j++) adj[i] |= 1u << j; }
            int done = 0;
            while (done < r) {
                i = (int)(rng() % n); j = (int)(rng() % n);
                if (i >= j) continue;
                if ((adj[i] >> j) & 1) { adj[i] &= ~(1u << j); adj[j] |= 1u << i; done++; }
            }
        }
        run_dp(n, adj);
        uint64_t H = hsub(n, full);
        int v;
        for (v = 0; v < n; v++) {
            uint64_t Hv = hsub(n, full & ~(1u << v));
            uint64_t k2 = H - Hv;
            if (k2 == 2) { inc2++; if (!((Hv == 1 || Hv == 3) && (H == 3 || H == 5))) { viol++; printf("VIOL inc2 n=%d H=%llu Hv=%llu v=%d masks:", n, (unsigned long long)H, (unsigned long long)Hv, v); for (i = 0; i < n; i++) printf(" %u", adj[i]); printf("\n"); } }
            if (k2 == 4) { inc4++; if (!(Hv == 1 || Hv == 5)) { viol++; printf("VIOL inc4 n=%d H=%llu Hv=%llu v=%d masks:", n, (unsigned long long)H, (unsigned long long)Hv, v); for (i = 0; i < n; i++) printf(" %u", adj[i]); printf("\n"); } }
            if (adj[v] != 0 && H >= 3) {
                uint64_t e = dp[(size_t)full * n + v];
                if (2 * e + 1 > H) { endv++; printf("VIOL endcount n=%d H=%llu e=%llu a=%d masks:", n, (unsigned long long)H, (unsigned long long)e, v); for (i = 0; i < n; i++) printf(" %u", adj[i]); printf("\n"); }
                if ((__uint128_t)e * bestH > (__uint128_t)bestE * H) { bestE = e; bestH = H; }
            }
        }
    }
    printf("SUMMARY rand n=%d r=%d samples=%ld seed=%llu: increment-2 pairs=%ld increment-4 pairs=%ld structural violations=%ld endcount violations=%ld max e/H=%llu/%llu\n", n, r, S, (unsigned long long)seed, inc2, inc4, viol, endv, (unsigned long long)bestE, (unsigned long long)bestH);
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s del|ext|rand n [r S seed]\n", argv[0]); return 1; }
    n = atoi(argv[2]);
    if (n < 2 || n + 1 > MAXN) { fprintf(stderr, "bad n\n"); return 1; }
    dp = (uint64_t *)malloc(sizeof(uint64_t) * ((size_t)1 << (n + 1)) * (n + 1));
    if (!strcmp(argv[1], "del")) mode_del();
    else if (!strcmp(argv[1], "ext")) mode_ext();
    else if (!strcmp(argv[1], "rand")) { if (argc < 6) { fprintf(stderr, "rand needs r S seed\n"); return 1; } mode_rand(atoi(argv[3]), atol(argv[4]), (uint64_t)strtoull(argv[5], 0, 10)); }
    else { fprintf(stderr, "unknown mode\n"); return 1; }
    return 0;
}
