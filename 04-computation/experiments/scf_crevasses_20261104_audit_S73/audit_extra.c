/*
 * audit_extra.c -- further independent checks for the THM-4866 audit (auditor's own code).
 *  (1) uniquely Hamiltonian classes: vertices of in-degree 1 / out-degree 1 (Shelor's maximum-outdegree insertion);
 *  (2) self-converse classes at even order: H mod 4 (THM-791), and H mod 8 when the condensation has no middle;
 *  (3) order 10 via ears with H <= 225: holes of the strong spectrum in [f(10), f(11)) = [125, 225);
 *  (4) interval chains IC(p) (parts >= 2) through order 17, and upset chains (parts >= 1) through order 12:
 *      H = prod(2^(p-1)+1) (parts >= 2) and exactly one Hamiltonian cycle.
 */
#define NO_MAIN
#include "audit_census.c"

static u64 big_count(int n, const uint32_t *adj, u64 *hc_out) {
    int full = (1 << n) - 1;
    u64 *dp = calloc((size_t)(full + 1) * n, sizeof(u64));
    for (int v = 0; v < n; v++) dp[(size_t)(1 << v) * n + v] = 1;
    for (int S = 1; S <= full; S++) for (int v = 0; v < n; v++) {
        u64 c = dp[(size_t)S * n + v]; if (!c) continue;
        uint32_t nb = adj[v] & ~(uint32_t)S;
        while (nb) { int w = __builtin_ctz(nb); nb &= nb - 1; dp[(size_t)(S | (1 << w)) * n + w] += c; }
    }
    u64 h = 0; for (int v = 0; v < n; v++) h += dp[(size_t)full * n + v];
    /* cycles: paths starting at vertex 0 */
    memset(dp, 0, sizeof(u64) * (size_t)(full + 1) * n);
    dp[(size_t)1 * n + 0] = 1;
    for (int S = 1; S <= full; S += 2) for (int v = 0; v < n; v++) {
        u64 c = dp[(size_t)S * n + v]; if (!c) continue;
        uint32_t nb = adj[v] & ~(uint32_t)S;
        while (nb) { int w = __builtin_ctz(nb); nb &= nb - 1; dp[(size_t)(S | (1 << w)) * n + w] += c; }
    }
    u64 hc = 0; for (int v = 1; v < n; v++) if (adj[v] & 1u) hc += dp[(size_t)full * n + v];
    *hc_out = hc; free(dp); return h;
}

static int parts[40], np, bad, cnt, minpart;
static void check_chain(int N) {
    int n = N + 1; uint32_t adj[32];
    for (int i = 0; i < n; i++) { adj[i] = 0; for (int j = i + 1; j < n; j++) adj[i] |= 1u << j; }
    int l = 0; u64 prod = 1;
    for (int k = 0; k < np; k++) { int l2 = l + parts[k]; adj[l] &= ~(1u << l2); adj[l2] |= 1u << l; prod *= ((u64)1 << (parts[k] - 1)) + 1; l = l2; }
    u64 hc, h = big_count(n, adj, &hc);
    cnt++;
    if (hc != 1 || (minpart >= 2 && h != prod)) bad++;
}
static void comps(int rem, int N) {
    if (!rem) { check_chain(N); return; }
    for (int p = minpart; p <= rem; p++) { parts[np++] = p; comps(rem - p, N); np--; }
}

int main(void) {
    for (int n = 1; n <= 9; n++) build(n);
    /* (1) and (2) */
    for (int n = 3; n <= 9; n++) {
        int m = ncls[n]; long hist[4][4] = {{0}}; long nu = 0, shelor = 0;
        long scmod4[4] = {0}, nomid_mod8[8] = {0}, nsc = 0;
#pragma omp parallel for schedule(dynamic, 64) reduction(+:nu, shelor, nsc)
        for (int i = 0; i < m; i++) {
            uint16_t t[MAXN], op[MAXN]; decode(n, cls[n][i], t);
            u64 hc = count_cycles(n, t);
            if (hc == 1) {
                nu++;
                int in1 = 0, out1 = 0, sh = 0;
                for (int v = 0; v < n; v++) {
                    int s = __builtin_popcount(t[v]);
                    if (s == 1) out1++;
                    if (n - 1 - s == 1) {
                        in1++;
                        /* Shelor: deleting a maximum-outdegree vertex leaves a uniquely Hamiltonian tournament */
                        if (n >= 5) {
                            uint16_t sub[MAXN]; int q = 0, idx[MAXN];
                            for (int u = 0; u < n; u++) if (u != v) idx[q++] = u;
                            for (int x = 0; x < q; x++) { sub[x] = 0; for (int y = 0; y < q; y++) if ((t[idx[x]] >> idx[y]) & 1) sub[x] |= 1u << y; }
                            if (count_cycles(q, sub) == 1) sh = 1;
                        } else sh = 1;
                    }
                }
                shelor += sh;
#pragma omp critical
                hist[in1 < 3 ? in1 : 3][out1 < 3 ? out1 : 3]++;
            }
            conv(n, t, op);
            if (!(n & 1) && canon(n, op) == cls[n][i]) {
                nsc++;
                u64 h = count_paths(n, t);
                int comp[MAXN]; int k = strong_components(n, t, comp);
#pragma omp critical
                { scmod4[h % 4]++; if (!(k & 1)) nomid_mod8[h % 8]++; }
            }
        }
        printf("n=%d: U(n) = %ld classes; (#in-degree-1 vertices, #out-degree-1 vertices) histogram:", n, nu);
        for (int a = 0; a < 4; a++) for (int b = 0; b < 4; b++) if (hist[a][b]) printf(" (%d,%d):%ld", a, b, hist[a][b]);
        printf("; classes with an in-degree-1 vertex whose deletion leaves a uniquely Hamiltonian tournament: %ld\n", shelor);
        if (!(n & 1)) {
            printf("   self-converse classes at even n: %ld; H mod 4 = 0,1,2,3: %ld %ld %ld %ld; no middle component: H mod 8 =",
                   nsc, scmod4[0], scmod4[1], scmod4[2], scmod4[3]);
            for (int r = 0; r < 8; r++) if (nomid_mod8[r]) printf(" %d:%ld", r, nomid_mod8[r]);
            printf("\n");
        }
        fflush(stdout);
    }
    /* (3) order 10 ears with H <= 225 */
    {
        int m = ncls[9]; const u64 LIM = 223;
        char *st = calloc(LIM + 2, 1);
        /* H of every class of order 9 */
        u64 *H9 = malloc(8 * (size_t)m), *C9 = malloc(8 * (size_t)m);
#pragma omp parallel for schedule(dynamic, 64)
        for (int i = 0; i < m; i++) { uint16_t t[MAXN]; decode(9, cls[9][i], t); H9[i] = count_paths(9, t); C9[i] = count_cycles(9, t); }
#pragma omp parallel for schedule(dynamic, 16)
        for (int c = 0; c < m; c++) {
            if (!C9[c] || H9[c] > LIM) continue;
            uint16_t a[MAXN]; decode(9, cls[9][c], a);
            for (int mask = 0; mask < 512; mask++) {
                uint16_t t[MAXN];
                for (int i = 0; i < 9; i++) t[i] = a[i] | (((mask >> i) & 1) ? 0 : (1u << 9));
                t[9] = (uint16_t)mask;
                u64 h = count_paths(10, t);
                if (h > LIM) continue;
                if (count_cycles(10, t)) st[h] = 1;
            }
        }
        printf("order 10, strong values in [125, 225): holes =");
        int k = 0; for (u64 v = 125; v < 225; v += 2) if (!st[v]) { printf(" %llu", v); k++; }
        printf(" (%d holes)\n", k);
        free(st); free(H9); free(C9);
    }
    /* (4) chains */
    {
        long tot = 0;
        minpart = 2;
        for (int n = 3; n <= 17; n++) { np = 0; bad = 0; cnt = 0; comps(n - 1, n - 1); tot += cnt; printf("IC chains order %d: %d (bad %d)\n", n, cnt, bad); fflush(stdout); }
        printf("total interval chains of order <= 17: %ld\n", tot);
        minpart = 1; long tot1 = 0, bad1 = 0;
        for (int n = 3; n <= 12; n++) { np = 0; bad = 0; cnt = 0; comps(n - 1, n - 1); tot1 += cnt; bad1 += bad; }
        printf("upset chains (parts >= 1, labelled compositions) of orders 3..12: %ld, with hc != 1: %ld\n", tot1, bad1);
    }
    return 0;
}
