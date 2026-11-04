/* aud_fewcycle.c -- independent audit code (fermat_self_duality_20261104_audit_S10).
 * Generates every strong tournament class with 1 <= hc <= K (hc = number of directed Hamiltonian cycles) at orders
 * 3..NMAX by vertex extension, and tests HYP-9464 (B) (the uniquely Hamiltonian deletion law) on the sole-exit-free
 * classes (minimum out-degree >= 2).
 * Completeness: a strong T of order n >= 4 has a vertex v with T - v strong (Moon 1966), and hc(T - v) <= hc(T)
 * (each Hamiltonian cycle of T - v takes v at an in->out switch, injectively).  So every strong class with hc <= K is a
 * one-vertex extension of a strong class with hc <= K one order down.
 * hc of an extension: Hamiltonian cycles through the new vertex z <-> Hamiltonian paths of the parent from an
 * out-neighbour of z to an in-neighbour of z, so hc = sum_{w in N+(z), u in N-(z)} P[w][u] with P the start/end path
 * matrix of the parent.
 * Canonical form (valid because hc >= 1): minimum, over every Hamiltonian cycle C and every starting vertex on C, of the
 * adjacency rows written in the cycle order.  Shares no code with the session's programs.
 * Output per order: t_k (all strong classes with hc = k) and u_k (no vertex of out-degree 1), and for the sole-exit-free
 * classes with hc = 2..5 the distribution of #{v : hc(T - v) = 1}; for hc = 3 also the (sole exits, sole entries)
 * profile of the uniquely Hamiltonian deletions. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <omp.h>

#define NM 14
#define K 5
typedef unsigned long long u64;
typedef struct { uint16_t r[NM]; } Key;

static int cmpkey(const void *a, const void *b) { return memcmp(a, b, sizeof(Key)); }

/* number of Hamiltonian cycles (subset DP from vertex 0) */
static u64 hc_dp(int n, const uint16_t *out, u64 *dp) {
    if (n < 3) return 0;
    int full = (1 << n) - 1;
    memset(dp, 0, sizeof(u64) * (size_t)(full + 1) * n);
    dp[1 * n + 0] = 1;
    for (int S = 1; S <= full; S += 2) for (int v = 0; v < n; v++) {
        u64 c = dp[(size_t)S * n + v]; if (!c) continue;
        unsigned nb = out[v] & ~S & full;
        while (nb) { int w = __builtin_ctz(nb); nb &= nb - 1; dp[(size_t)(S | (1 << w)) * n + w] += c; }
    }
    u64 h = 0; for (int v = 1; v < n; v++) if (out[v] & 1u) h += dp[(size_t)full * n + v];
    return h;
}
/* P[w][u] = Hamiltonian paths from w to u */
static void path_matrix(int n, const uint16_t *out, u64 P[NM][NM], u64 *dp) {
    int full = (1 << n) - 1;
    for (int w = 0; w < n; w++) {
        memset(dp, 0, sizeof(u64) * (size_t)(full + 1) * n);
        dp[(size_t)(1 << w) * n + w] = 1;
        for (int S = 1; S <= full; S++) {
            if (!((S >> w) & 1)) continue;
            for (int v = 0; v < n; v++) {
                u64 c = dp[(size_t)S * n + v]; if (!c) continue;
                unsigned nb = out[v] & ~S & full;
                while (nb) { int x = __builtin_ctz(nb); nb &= nb - 1; dp[(size_t)(S | (1 << x)) * n + x] += c; }
            }
        }
        for (int u = 0; u < n; u++) P[w][u] = (u == w) ? 0 : dp[(size_t)full * n + u];
    }
}
/* enumerate Hamiltonian cycles from vertex 0 (at most cap) */
static __thread int cyc[8][NM], ncyc;
static void dfs_cyc(int n, const uint16_t *out, int *path, int len, unsigned used, int cap) {
    if (ncyc > cap) return;
    int v = path[len - 1];
    if (len == n) { if (out[v] & 1u) { if (ncyc < 8) memcpy(cyc[ncyc], path, sizeof(int) * n); ncyc++; } return; }
    unsigned nb = out[v] & ~used & ((1u << n) - 1);
    while (nb) { int w = __builtin_ctz(nb); nb &= nb - 1; path[len] = w; dfs_cyc(n, out, path, len + 1, used | (1u << w), cap); }
}
static Key canon(int n, const uint16_t *out) {
    int path[NM]; path[0] = 0; ncyc = 0;
    dfs_cyc(n, out, path, 1, 1u, 7);
    if (ncyc < 1 || ncyc > 7) { fprintf(stderr, "canon: bad cycle count %d\n", ncyc); exit(1); }
    Key best; int have = 0;
    for (int c = 0; c < ncyc; c++) for (int r = 0; r < n; r++) {
        int ord[NM], pos[NM]; for (int i = 0; i < n; i++) { ord[i] = cyc[c][(r + i) % n]; pos[ord[i]] = i; }
        Key k; memset(&k, 0, sizeof k);
        for (int i = 0; i < n; i++) { uint16_t row = 0; for (int j = 0; j < n; j++) if ((out[ord[i]] >> ord[j]) & 1u) row |= (uint16_t)(1u << j); k.r[i] = row; }
        if (!have || memcmp(&k, &best, sizeof k) < 0) { best = k; have = 1; }
    }
    return best;
}
/* |Aut| of a Hamiltonian tournament = number of (cycle, start) orderings attaining the canonical minimum
 * (Aut acts freely on these orderings, and two optimal orderings differ by an automorphism) */
static int aut_count(int n, const uint16_t *out) {
    Key best = canon(n, out);
    int path[NM]; path[0] = 0; ncyc = 0;
    dfs_cyc(n, out, path, 1, 1u, 7);
    int cnt = 0;
    for (int c = 0; c < ncyc; c++) for (int r = 0; r < n; r++) {
        int ord[NM]; for (int i = 0; i < n; i++) ord[i] = cyc[c][(r + i) % n];
        Key k; memset(&k, 0, sizeof k);
        for (int i = 0; i < n; i++) { uint16_t row = 0; for (int j = 0; j < n; j++) if ((out[ord[i]] >> ord[j]) & 1u) row |= (uint16_t)(1u << j); k.r[i] = row; }
        if (!memcmp(&k, &best, sizeof k)) cnt++;
    }
    return cnt;
}
static void delete_v(int n, const uint16_t *out, int x, uint16_t *res) {
    int k = 0;
    for (int i = 0; i < n; i++) {
        if (i == x) continue;
        uint16_t r = 0; int kk = 0;
        for (int j = 0; j < n; j++) { if (j == x) continue; if ((out[i] >> j) & 1u) r |= (uint16_t)(1u << kk); kk++; }
        res[k++] = r;
    }
}

int main(int argc, char **argv) {
    int NMAX = argc > 1 ? atoi(argv[1]) : 11;
    printf("aud_fewcycle (independent audit code): strong classes with 1 <= hc <= %d, orders 3..%d\n", K, NMAX);
    u64 *dp = malloc(sizeof(u64) * (size_t)(1 << NM) * NM);
    /* order 3: C3 */
    Key *cur = malloc(sizeof(Key)); memset(cur, 0, sizeof(Key));
    cur[0].r[0] = 2; cur[0].r[1] = 4; cur[0].r[2] = 1;   /* 0 -> 1 -> 2 -> 0 */
    long ncur = 1;
    for (int n = 3; n <= NMAX; n++) {
        double t0 = omp_get_wtime();
        if (n > 3) {
            /* extend every class of order n-1 */
            int m = n - 1;
            long cap = 1 << 20, cnt = 0; Key *cand = malloc(sizeof(Key) * cap);
            #pragma omp parallel
            {
                u64 *tdp = malloc(sizeof(u64) * (size_t)(1 << NM) * NM);
                long lcap = 1 << 16, lcnt = 0; Key *loc = malloc(sizeof(Key) * lcap);
                #pragma omp for schedule(dynamic, 8)
                for (long c = 0; c < ncur; c++) {
                    const uint16_t *pr = cur[c].r;
                    u64 P[NM][NM]; path_matrix(m, pr, P, tdp);
                    for (unsigned mask = 0; mask < (1u << m); mask++) {
                        u64 h = 0;
                        for (int w = 0; w < m; w++) if ((mask >> w) & 1u) for (int u = 0; u < m; u++) if (!((mask >> u) & 1u)) h += P[w][u];
                        if (h < 1 || h > K) continue;
                        uint16_t t[NM];
                        for (int i = 0; i < m; i++) t[i] = pr[i] | (((mask >> i) & 1u) ? 0 : (uint16_t)(1u << m));
                        t[m] = (uint16_t)mask;
                        if (lcnt == lcap) { lcap *= 2; loc = realloc(loc, sizeof(Key) * lcap); }
                        loc[lcnt++] = canon(n, t);
                    }
                }
                #pragma omp critical
                {
                    if (cnt + lcnt > cap) { while (cnt + lcnt > cap) cap *= 2; cand = realloc(cand, sizeof(Key) * cap); }
                    memcpy(cand + cnt, loc, sizeof(Key) * lcnt); cnt += lcnt;
                }
                free(loc); free(tdp);
            }
            qsort(cand, cnt, sizeof(Key), cmpkey);
            long u = 0; for (long i = 0; i < cnt; i++) if (i == 0 || memcmp(&cand[i], &cand[i - 1], sizeof(Key))) cand[u++] = cand[i];
            free(cur); cur = cand; ncur = u;
        }
        /* statistics */
        long t[K + 1] = {0}, us[K + 1] = {0}, uh_se[4] = {0}, aut1[K + 1] = {0}, aut3[K + 1] = {0}, autx[K + 1] = {0};
        long dh[K + 1][NM + 1]; memset(dh, 0, sizeof dh);
        long prof_n = 0; struct { int key[3][2]; long c; } prof[64]; int nprof = 0;
        for (long c = 0; c < ncur; c++) {
            const uint16_t *o = cur[c].r;
            u64 h = hc_dp(n, o, dp);
            if (h < 1 || h > K) { fprintf(stderr, "bad hc %llu\n", h); return 1; }
            t[h]++;
            if (h == 1) { int se = 0; for (int v = 0; v < n; v++) if (__builtin_popcount(o[v]) == 1) se++; if (se <= 3) uh_se[se]++; }
            int mo = 99, ok_entry[NM], ok_exit[NM];
            for (int v = 0; v < n; v++) { int d = __builtin_popcount(o[v]); if (d < mo) mo = d; (void)ok_entry; (void)ok_exit; }
            if (mo < 2) continue;
            us[h]++;
            { int a = aut_count(n, o); if (a == 1) aut1[h]++; else if (a == 3) aut3[h]++; else autx[h]++; }
            if (h < 2) continue;
            int k = 0, pe[NM][2], np = 0;
            for (int v = 0; v < n; v++) {
                uint16_t r[NM]; delete_v(n, o, v, r);
                if (hc_dp(n - 1, r, dp) == 1) {
                    k++;
                    int e1 = 0, i1 = 0;
                    for (int x = 0; x < n - 1; x++) { int d = __builtin_popcount(r[x]); if (d == 1) e1++; if (n - 2 - d == 1) i1++; }
                    pe[np][0] = e1; pe[np][1] = i1; np++;
                }
            }
            dh[h][k]++;
            if (h == 3 && np == 3) {
                /* sort the three (sole exits, sole entries) pairs */
                for (int a = 0; a < 3; a++) for (int b = a + 1; b < 3; b++)
                    if (pe[b][0] < pe[a][0] || (pe[b][0] == pe[a][0] && pe[b][1] < pe[a][1])) { int x0 = pe[a][0], x1 = pe[a][1]; pe[a][0] = pe[b][0]; pe[a][1] = pe[b][1]; pe[b][0] = x0; pe[b][1] = x1; }
                int f = -1;
                for (int q = 0; q < nprof; q++) if (!memcmp(prof[q].key, pe, sizeof(int) * 6)) { f = q; break; }
                if (f < 0) { f = nprof++; memcpy(prof[f].key, pe, sizeof(int) * 6); prof[f].c = 0; }
                prof[f].c++; prof_n++;
            }
        }
        printf("\norder %d: %ld strong classes with hc <= %d  [%.1f s]\n", n, ncur, K, omp_get_wtime() - t0);
        printf("  t_k (k=1..5):"); for (int k = 1; k <= K; k++) printf(" %ld", t[k]); printf("\n");
        printf("  u_k (k=1..5, no out-degree-1 vertex):"); for (int k = 1; k <= K; k++) printf(" %ld", us[k]); printf("\n");
        printf("  uniquely Hamiltonian classes by number of sole exits: 1: %ld, 2: %ld, 3: %ld\n", uh_se[1], uh_se[2], uh_se[3]);
        printf("  sole-exit-free classes by |Aut| (hc = 3, 4, 5): |Aut|=1: %ld %ld %ld; |Aut|=3: %ld %ld %ld; other: %ld %ld %ld\n",
               aut1[3], aut1[4], aut1[5], aut3[3], aut3[4], aut3[5], autx[3], autx[4], autx[5]);
        for (int h = 2; h <= K; h++) if (us[h]) {
            printf("  DEL n=%d hc=%d classes=%ld  #(v: hc(T-v)=1) ->", n, h, us[h]);
            for (int k = 0; k <= NM; k++) if (dh[h][k]) printf(" %d:%ld", k, dh[h][k]);
            printf("\n");
        }
        if (prof_n) {
            printf("  hc=3 profiles of the three UH deletions ((sole exits, sole entries) x3):");
            for (int q = 0; q < nprof; q++) printf(" ((%d,%d),(%d,%d),(%d,%d)):%ld", prof[q].key[0][0], prof[q].key[0][1], prof[q].key[1][0], prof[q].key[1][1], prof[q].key[2][0], prof[q].key[2][1], prof[q].c);
            printf("\n");
        }
        fflush(stdout);
    }
    return 0;
}
