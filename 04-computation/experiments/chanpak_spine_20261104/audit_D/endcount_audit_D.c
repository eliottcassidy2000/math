/* endcount_audit_D.c -- audit D of THM-4876 (7), written independently of the session's scripts.
 * Modes:
 *   endcount_audit_D classes n [soleexit]   : read gentourng's default ASCII output (upper triangle, row by row,
 *                                             '1' at (i,j) = arc i->j) for order n from stdin; one class per line.
 *   endcount_audit_D labelled n [soleexit]  : all 2^C(n,2) labelled tournaments of order n.
 *   endcount_audit_D family m               : TT_m with the arc (1,m) reversed, vertices 1..m -> indices 0..m-1.
 * Checks, per tournament X (h = H(X)) and vertex a:
 *   (B)  non-sink a: 2 e_a <= h - 1; non-source a: 2 s_a <= h - 1; also the sharper 2 e_a + e_O <= h (e_O = paths ending in N+(a)).
 *   (R)  max e_a/h over non-sinks (exact fraction), and all equality cases 2 e_a = h - 1.
 *   (C)  for h >= 5: no vertex has e_a in {h-1, h-2}; for h = 3, e_a = 1 forces X = TT_j => C3 with a in the C3.
 *   (S)  [soleexit] direct count of H(X + v), v -> a, all others -> v, compared with 2h - e_a + sigma,
 *        sigma = [P A_a](X) - h, [P A_a] = sum_{S subset V-a} h(X[S]) s_a(X[V-S]) (h(empty) = 1); sigma >= 0; sigma = 0 count;
 *        2 H(X+v) >= 3h + 1.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MAXN 12
typedef unsigned long long u64;

static void dp_end(int nn, int a[MAXN][MAXN], u64 *E) {   /* E[mask*nn+v] = # Hamiltonian paths of X[mask] ending at v */
    int full = 1 << nn;
    memset(E, 0, sizeof(u64) * (size_t)full * nn);
    for (int v = 0; v < nn; v++) E[(1u << v) * nn + v] = 1;
    for (int mask = 1; mask < full; mask++)
        for (int last = 0; last < nn; last++) {
            if (!((mask >> last) & 1)) continue;
            u64 c = E[mask * nn + last];
            if (!c) continue;
            for (int w = 0; w < nn; w++) {
                if (((mask >> w) & 1) || !a[last][w]) continue;
                E[(mask | (1 << w)) * nn + w] += c;
            }
        }
}
static u64 hsub(int nn, const u64 *E, int mask) { if (!mask) return 1; u64 s = 0; for (int v = 0; v < nn; v++) if ((mask >> v) & 1) s += E[mask * nn + v]; return s; }

static int n, adj[MAXN][MAXN], adjT[MAXN][MAXN], adj2[MAXN][MAXN];
static u64 *E, *S, *E2;
/* tallies */
static u64 ntour = 0, npairs_nonsink = 0, viol_B = 0, viol_dual = 0, viol_sharp = 0, viol_C = 0, viol_C3 = 0, n_eq = 0;
static u64 best_e = 0, best_h = 1; static char best_line[256] = ""; static int best_a = -1;
static u64 se_cases = 0, se_lawfail = 0, se_sigma_neg = 0, se_sigma0 = 0, se_floorfail = 0, se_maxsigma = 0;
static u64 se_min2k = ~0ull; static u64 se_min2k_h = 0;  /* minimal 2*(H(X+v)-h) relative: track min of (H2-h)/((h+1)/2) ... keep as pair */
static u64 se_tight = 0; /* cases with 2 H(X+v) == 3h + 1 */
static u64 eq_by_h[1 << 12];
static char line[256];

static void process(void) {
    int full = (1 << n) - 1;
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) adjT[i][j] = adj[j][i];
    dp_end(n, adj, E); dp_end(n, adjT, S);
    u64 h = hsub(n, E, full);
    if (!(h & 1)) { fprintf(stderr, "Redei violated?! h=%llu line=%s\n", h, line); exit(1); }
    ntour++;
    for (int a = 0; a < n; a++) {
        int outdeg = 0, indeg = 0; for (int j = 0; j < n; j++) { outdeg += adj[a][j]; indeg += adj[j][a]; }
        u64 ea = E[full * n + a], sa = S[full * n + a];
        u64 eO = 0; for (int j = 0; j < n; j++) if (adj[a][j]) eO += E[full * n + j];
        if (outdeg == 0) { if (ea != h) { fprintf(stderr, "sink with e_a != h ?!\n"); exit(1); } }
        else {
            npairs_nonsink++;
            if (2 * ea > h - 1) { viol_B++; printf("VIOLATION(B): n=%d line=%s a=%d h=%llu e_a=%llu\n", n, line, a, h, ea); }
            if (2 * ea + eO > h) { viol_sharp++; printf("VIOLATION(sharp): n=%d line=%s a=%d h=%llu e_a=%llu e_O=%llu\n", n, line, a, h, ea, eO); }
            if (ea * best_h > best_e * h) { best_e = ea; best_h = h; best_a = a; strncpy(best_line, line, 255); }
            if (2 * ea == h - 1) { n_eq++; if (h < (1u << 12)) eq_by_h[h]++;
                printf("EQUALITY: n=%d line=%s a=%d h=%llu e_a=%llu outdeg=%d\n", n, line, a, h, ea, outdeg); }
        }
        if (indeg != 0 && 2 * sa > h - 1) { viol_dual++; printf("VIOLATION(dual): n=%d line=%s a=%d h=%llu s_a=%llu\n", n, line, a, h, sa); }
        /* (C) */
        if (h >= 5 && (ea == h - 1 || ea == h - 2)) { viol_C++; printf("VIOLATION(C): n=%d line=%s a=%d h=%llu e_a=%llu\n", n, line, a, h, ea); }
        if (h == 3 && ea == 1) {
            /* expect: C = {v : e_v > 0} is a 3-cycle, last strong component, rest transitive */
            int C = 0, cnt = 0; for (int v = 0; v < n; v++) if (E[full * n + v] > 0) { C |= 1 << v; cnt++; }
            int ok = (cnt == 3) && hsub(n, E, C) == 3 && hsub(n, E, full ^ C) == 1;
            for (int u = 0; u < n && ok; u++) if (!((C >> u) & 1)) for (int v = 0; v < n; v++) if (((C >> v) & 1) && !adj[u][v]) ok = 0;
            if (!ok) { viol_C3++; printf("VIOLATION(C3 case): n=%d line=%s a=%d\n", n, line, a); }
        }
    }
}
static void soleexit(void) {
    int full = (1 << n) - 1, full2 = (1 << (n + 1)) - 1;
    u64 h = hsub(n, E, full);
    for (int a = 0; a < n; a++) {
        int outdeg = 0; for (int j = 0; j < n; j++) outdeg += adj[a][j];
        if (outdeg == 0) continue;
        u64 ea = E[full * n + a];
        /* [P A_a] from the definition */
        u64 PA = 0; int rest = full ^ (1 << a);
        for (int sub = rest;; sub = (sub - 1) & rest) { PA += hsub(n, E, sub) * S[(full ^ sub) * n + a]; if (!sub) break; }
        u64 sigma_signed_ok = (PA >= h); u64 sigma = PA - h;
        /* direct count of H(X + v) */
        for (int i = 0; i <= n; i++) for (int j = 0; j <= n; j++) adj2[i][j] = 0;
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) adj2[i][j] = adj[i][j];
        adj2[n][a] = 1; for (int x = 0; x < n; x++) if (x != a) adj2[x][n] = 1;
        dp_end(n + 1, adj2, E2);
        u64 H2 = hsub(n + 1, E2, full2);
        se_cases++;
        if (!sigma_signed_ok) { se_sigma_neg++; printf("SIGMA<0: n=%d line=%s a=%d\n", n, line, a); }
        if (H2 != 2 * h - ea + sigma) { se_lawfail++; printf("LAWFAIL: n=%d line=%s a=%d h=%llu e_a=%llu PA=%llu H2=%llu\n", n, line, a, h, ea, PA, H2); }
        if (sigma == 0) se_sigma0++;
        if (sigma > se_maxsigma) se_maxsigma = sigma;
        if (2 * H2 < 3 * h + 1) { se_floorfail++; printf("FLOORFAIL: n=%d line=%s a=%d h=%llu H2=%llu\n", n, line, a, h, H2); }
        if (2 * H2 == 3 * h + 1) { se_tight++; printf("TIGHT(3h+1)/2: n=%d line=%s a=%d h=%llu H2=%llu\n", n, line, a, h, H2); }
        /* minimal normalised rise (H2-h)/h: track min by cross-multiplication */
        if ((H2 - h) * se_min2k_h < se_min2k * h || se_min2k_h == 0) { se_min2k = H2 - h; se_min2k_h = h; }
    }
}
static void report(const char *mode, int do_se) {
    printf("MODE=%s n=%d tournaments=%llu nonsink_pairs=%llu\n", mode, n, ntour, npairs_nonsink);
    printf("  (B) violations of 2e_a<=h-1: %llu ; dual violations: %llu ; violations of 2e_a+e_O<=h: %llu\n", viol_B, viol_dual, viol_sharp);
    printf("  (R) max e_a/h over non-sinks = %llu/%llu (= %.6f) at a=%d, tournament %s ; equality cases (X,a): %llu\n", best_e, best_h, (double)best_e / best_h, best_a, best_line, n_eq);
    printf("      equality cases by h:"); for (int hh = 0; hh < (1 << 12); hh++) if (eq_by_h[hh]) printf(" h=%d:%llu", hh, eq_by_h[hh]); printf("\n");
    printf("  (C) h>=5 with e_a in {h-1,h-2}: %llu ; h=3,e_a=1 not TT_j=>C3-with-a-in-C3: %llu\n", viol_C, viol_C3);
    if (do_se) printf("  (S) sole-exit cases=%llu lawfail=%llu sigma<0=%llu sigma=0:%llu maxsigma=%llu floorfail(2H<3h+1)=%llu tight(2H=3h+1)=%llu min rise (H2-h)/h = %llu/%llu = %.6f\n",
        se_cases, se_lawfail, se_sigma_neg, se_sigma0, se_maxsigma, se_floorfail, se_tight, se_min2k, se_min2k_h, (double)se_min2k / se_min2k_h);
}
int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage\n"); return 2; }
    n = atoi(argv[2]); int do_se = (argc > 3 && strcmp(argv[3], "soleexit") == 0);
    E = malloc(sizeof(u64) * (1u << MAXN) * MAXN); S = malloc(sizeof(u64) * (1u << MAXN) * MAXN); E2 = malloc(sizeof(u64) * (1u << MAXN) * MAXN);
    memset(eq_by_h, 0, sizeof eq_by_h);
    if (strcmp(argv[1], "family") == 0) {
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) adj[i][j] = (i < j);
        adj[0][n - 1] = 0; adj[n - 1][0] = 1;
        snprintf(line, 255, "TT_%d with (1,%d) reversed", n, n);
        process();
        int full = (1 << n) - 1; u64 h = hsub(n, E, full);
        printf("FAMILY m=%d: h=%llu (2^(m-2)+1=%llu) e_v =", n, h, (1ull << (n - 2)) + 1);
        for (int v = 0; v < n; v++) printf(" %llu", E[full * n + v]);
        printf(" ; e_(m-1)=%llu = 2^(m-3)=%llu ; vertices (1-based) with 2e=h-1:", E[full * n + n - 2], 1ull << (n - 3));
        for (int v = 0; v < n; v++) if (2 * E[full * n + v] == h - 1) printf(" %d", v + 1);
        printf("\n");
        soleexit(); report("family", 1); return 0;
    }
    if (strcmp(argv[1], "labelled") == 0) {
        int m = n * (n - 1) / 2; u64 total = 1ull << m;
        for (u64 bits = 0; bits < total; bits++) {
            int k = 0; for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { int b = (bits >> k) & 1; adj[i][j] = b; adj[j][i] = !b; k++; }
            snprintf(line, 255, "bits=%llu", bits); process(); if (do_se) soleexit();
        }
        report("labelled", do_se); return 0;
    }
    /* classes from gentourng */
    while (fgets(line, sizeof line, stdin)) {
        size_t L = strlen(line); while (L && (line[L - 1] == '\n' || line[L - 1] == '\r')) line[--L] = 0;
        if (L != (size_t)(n * (n - 1) / 2)) { fprintf(stderr, "bad line length %zu for n=%d: %s\n", L, n, line); return 1; }
        int k = 0; for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { int b = line[k++] == '1'; adj[i][j] = b; adj[j][i] = !b; }
        process(); if (do_se) soleexit();
    }
    report("classes", do_se); return 0;
}
