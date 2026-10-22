/* Independent auditor code (audit of opus-2026-10-22-S54 drafts).
 * Reads tournaments, one per line: "code ..." where code = upper-triangle bits row-major, '1' iff i -> j (i < j).
 * For each tournament prints:
 *   H  c_0 .. c_(n-1)  [| N_0 .. N_(n-1)]
 * H   : number of Hamiltonian paths (Held-Karp on arcs).
 * c_L : coefficient of t^L in sum over orderings s of prod_j (1 + t S[s_j][s_(j+1)])   (method 1).
 * N_L : signed count of directed linear forests with L arcs, computed by a DIFFERENT algorithm (method 2):
 *       g[B] = signed Hamiltonian-path sum of the subtournament on B (|B| >= 2), then a set-partition DP
 *       f[U] = sum over the block B containing min(U) (B a path block or a singleton) t^(|B|-1) g[B] f[U \ B].
 * Usage: forest_audit n [m2]   (m2 = also run method 2; n <= 16 for method 1, n <= 14 for method 2)        */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int n;
static int A[32][32];
static int S[32][32];

static unsigned long long held_karp(void) {
    size_t N = (size_t)1 << n;
    unsigned long long *dp = calloc(N * n, sizeof *dp);
    for (int v = 0; v < n; v++) dp[((size_t)1 << v) * n + v] = 1;
    for (size_t m = 1; m < N; m++)
        for (int v = 0; v < n; v++) {
            unsigned long long c = dp[m * n + v];
            if (!c) continue;
            for (int w = 0; w < n; w++)
                if (!((m >> w) & 1) && A[v][w]) dp[(m | ((size_t)1 << w)) * n + w] += c;
        }
    unsigned long long t = 0;
    for (int v = 0; v < n; v++) t += dp[(N - 1) * n + v];
    free(dp);
    return t;
}

/* method 1: dp[mask][v][L] */
static void perm_product(long long *c) {
    size_t N = (size_t)1 << n;
    long long *dp = calloc(N * n * n, sizeof *dp);
    for (int v = 0; v < n; v++) dp[(((size_t)1 << v) * n + v) * n] = 1;
    for (size_t m = 1; m < N; m++) {
        int pc = __builtin_popcountll(m);
        for (int v = 0; v < n; v++) {
            if (!((m >> v) & 1)) continue;
            long long *src = dp + (m * n + v) * n;
            for (int w = 0; w < n; w++) {
                if ((m >> w) & 1) continue;
                long long *dst = dp + ((m | ((size_t)1 << w)) * n + w) * n;
                int s = S[v][w];
                for (int L = 0; L < pc; L++) {      /* degree <= pc - 1; unsigned wrap-around = exact mod 2^64 */
                    unsigned long long x = (unsigned long long)src[L];
                    if (!x) continue;
                    dst[L] = (long long)((unsigned long long)dst[L] + x);
                    dst[L + 1] = (long long)((unsigned long long)dst[L + 1] + (s > 0 ? x : 0ULL - x));
                }
            }
        }
    }
    for (int L = 0; L < n; L++) {
        unsigned long long t = 0;
        for (int v = 0; v < n; v++) t += (unsigned long long)dp[((N - 1) * n + v) * n + L];
        c[L] = (long long)t;
    }
    free(dp);
}

/* method 2 */
static void forest_dp(long long *Nl) {
    size_t N = (size_t)1 << n;
    long long *g = calloc(N * n, sizeof *g);      /* signed path sums ending at v */
    long long *P = calloc(N, sizeof *P);          /* P[B] = sum_v g[B][v], |B| >= 2 */
    for (int v = 0; v < n; v++) g[((size_t)1 << v) * n + v] = 1;
    for (size_t m = 1; m < N; m++)
        for (int v = 0; v < n; v++) {
            long long x = g[m * n + v];
            if (!x) continue;
            for (int w = 0; w < n; w++)
                if (!((m >> w) & 1)) g[(m | ((size_t)1 << w)) * n + w] += S[v][w] * x;
        }
    for (size_t m = 1; m < N; m++) {
        if (__builtin_popcountll(m) < 2) continue;
        long long t = 0;
        for (int v = 0; v < n; v++) t += g[m * n + v];
        P[m] = t;
    }
    free(g);
    /* f[U][L], L <= n-1 */
    long long *f = calloc(N * n, sizeof *f);
    f[0] = 1;
    for (size_t U = 1; U < N; U++) {
        int lo = __builtin_ctzll(U);
        size_t rest = U & ~((size_t)1 << lo);
        /* singleton block {lo} */
        for (int L = 0; L < n; L++) f[U * n + L] += f[rest * n + L];
        /* blocks B = {lo} u sub, sub nonempty subset of rest */
        for (size_t sub = rest; sub; sub = (sub - 1) & rest) {
            size_t B = sub | ((size_t)1 << lo);
            long long p = P[B];
            if (!p) continue;
            int k = __builtin_popcountll(B) - 1;    /* arcs in the path block */
            size_t R = rest & ~sub;
            for (int L = 0; L + k < n; L++) {
                long long x = f[R * n + L];
                if (x) f[U * n + L + k] += p * x;
            }
        }
    }
    for (int L = 0; L < n; L++) Nl[L] = f[(N - 1) * n + L];
    free(f); free(P);
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: forest_audit n [m2]\n"); return 1; }
    n = atoi(argv[1]);
    int m2 = argc >= 3 && !strcmp(argv[2], "m2");
    char buf[8192];
    while (fgets(buf, sizeof buf, stdin)) {
        char code[1024];
        if (sscanf(buf, "%1023s", code) != 1) continue;
        if ((int)strlen(code) != n * (n - 1) / 2) { fprintf(stderr, "bad code length\n"); return 2; }
        int c = 0;
        for (int i = 0; i < n; i++) {
            A[i][i] = 0; S[i][i] = 0;
            for (int j = i + 1; j < n; j++) {
                A[i][j] = code[c++] == '1'; A[j][i] = !A[i][j];
                S[i][j] = A[i][j] ? 1 : -1; S[j][i] = -S[i][j];
            }
        }
        long long cl[32];
        perm_product(cl);
        unsigned long long H = held_karp();
        printf("%llu", H);
        for (int L = 0; L < n; L++) printf(" %lld", cl[L]);
        if (m2) {
            long long Nl[32];
            forest_dp(Nl);
            printf(" |");
            for (int L = 0; L < n; L++) printf(" %lld", Nl[L]);
        }
        printf("\n");
        fflush(stdout);
    }
    return 0;
}
