/* audit_dp.c -- independent audit of THM-4820 (f over count vectors).  Written from scratch for the audit.
 *
 * f(u) = minimum over binary merge trees whose leaves are the unit charges of the multiset u = (a,b,c)
 * (a ones, b omegas, c omega^2's) of the sum over merges of max(0, beta), beta = Q(y+z) - Q(y) - Q(z).
 *
 * Differences from the producer's merge_debt.c (deliberate, for independence):
 *   - beta is computed from Eisenstein VALUES via a precomputed table Qv[u] = Q(value(u)), value(a,b,c) = (a-c) + (b-c) w,
 *     beta = Qv[u] - Qv[u1] - Qv[u - u1]; the count formula 3u.u' - |u||u'| is NOT used in the DP;
 *   - states are processed in lexicographic order (every proper sub-vector precedes its parent), not by total size;
 *   - all ORDERED splits are scanned (no symmetry pruning);
 *   - optional second pass ('verify'): push-style check of the Bellman inequality f(u1+u2) <= beta+ + f(u1) + f(u2)
 *     for EVERY pair of nonzero sub-vectors (this certifies f as a lower bound on every tree, by induction on size),
 *     and that every non-unit state attains its value by some split (so f is achieved by a tree);
 *   - structural checks: f(u+(1,1,1)) <= f(u), f(n,0,0) = n(n-1), f(u) >= floor(Q/3), f(u) >= Q - |u|, f(a,b,c) = f(a,c,b).
 *
 * Usage: ./audit_dp A B C [verify] [lastmin] [tree a b c]
 *   prints "F n k f" for every n >= 0, k >= 0 with n + k <= A, k <= min(B,C)  (f(n+k,k,k)), and check summaries.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int A, B, C;
static inline long id3(int a, int b, int c) { return ((long)a * (B + 1) + b) * (C + 1) + c; }
static inline int Qab(int a, int b) { return a * a - a * b + b * b; }

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage\n"); return 1; }
    A = atoi(argv[1]); B = atoi(argv[2]); C = atoi(argv[3]);
    int verify = 0, tree = 0, ta = 0, tb = 0, tc = 0, lastmin = 0;
    for (int i = 4; i < argc; i++) {
        if (!strcmp(argv[i], "verify")) verify = 1;
        if (!strcmp(argv[i], "lastmin")) lastmin = 1;   /* break ties by the LAST minimising split (a different certificate) */
        if (!strcmp(argv[i], "tree") && i + 3 < argc) { tree = 1; ta = atoi(argv[i + 1]); tb = atoi(argv[i + 2]); tc = atoi(argv[i + 3]); i += 3; }
    }
    long N = (long)(A + 1) * (B + 1) * (C + 1);
    int *f = malloc(sizeof(int) * N), *Qv = malloc(sizeof(int) * N);
    long *how = tree ? malloc(sizeof(long) * N) : NULL;
    for (int a = 0; a <= A; a++) for (int b = 0; b <= B; b++) for (int c = 0; c <= C; c++)
        Qv[id3(a, b, c)] = Qab(a - c, b - c);
    const int INF = 1 << 29;
    long iters = 0;
    for (int a = 0; a <= A; a++) for (int b = 0; b <= B; b++) for (int c = 0; c <= C; c++) {
        long u = id3(a, b, c); int s = a + b + c;
        if (s == 0) { f[u] = INF; continue; }
        if (s == 1) { f[u] = 0; if (how) how[u] = -1; continue; }
        int best = INF; long barg = -1; int qu = Qv[u];
        for (int x = 0; x <= a; x++) for (int y = 0; y <= b; y++) {
            long base1 = id3(x, y, 0), base2 = id3(a - x, b - y, 0);
            for (int z = 0; z <= c; z++) {
                int s1 = x + y + z; if (s1 == 0 || s1 == s) continue;
                long u1 = base1 + z, u2 = base2 + (c - z);
                int beta = qu - Qv[u1] - Qv[u2];
                int v = (beta > 0 ? beta : 0) + f[u1] + f[u2];
                if (v < best || (lastmin && v == best)) { best = v; barg = u1; }
            }
            iters += c + 1;
        }
        f[u] = best; if (how) how[u] = barg;
    }
    fprintf(stderr, "DP done: %ld states, %ld split evaluations\n", N, iters);

    /* structural checks */
    long bad_a2 = 0, bad_a3 = 0, bad_phi = 0, bad_sum = 0, bad_sym = 0, n_a2 = 0;
    for (int a = 0; a <= A; a++) for (int b = 0; b <= B; b++) for (int c = 0; c <= C; c++) {
        long u = id3(a, b, c); if (a + b + c == 0) continue;
        if (a < A && b < B && c < C) { n_a2++; if (f[id3(a + 1, b + 1, c + 1)] > f[u]) bad_a2++; }
        if (f[u] < Qv[u] / 3) bad_phi++;
        if (f[u] < Qv[u] - (a + b + c)) bad_sum++;
        if (c <= B && b <= C && f[u] != f[id3(a, c, b)]) bad_sym++;
    }
    /* the next line needs c <= B and b <= C to index (a,c,b); it is skipped otherwise */
    for (int n = 1; n <= A; n++) if (f[id3(n, 0, 0)] != n * (n - 1)) bad_a3++;
    printf("CHECK A2 f(u+(1,1,1)) <= f(u): %ld violations among %ld pairs\n", bad_a2, n_a2);
    printf("CHECK A3 f(n,0,0) = n(n-1) for n <= %d: %ld violations\n", A, bad_a3);
    printf("CHECK B  f(u) >= floor(Q(value)/3): %ld violations\n", bad_phi);
    printf("CHECK E  f(u) >= Q(value) - |u|: %ld violations\n", bad_sum);
    printf("CHECK    f(a,b,c) = f(a,c,b) (where both in box): %ld violations\n", bad_sym);

    if (verify) {   /* push-style Bellman inequality over every pair (u1, u2) of nonzero vectors with u1 + u2 in the box */
        long viol = 0, pairs = 0;
        char *attained = calloc(N, 1);
        for (int x = 0; x <= A; x++) for (int y = 0; y <= B; y++) for (int z = 0; z <= C; z++) {
            if (x + y + z == 0) continue; long u1 = id3(x, y, z);
            for (int p = 0; p + x <= A; p++) for (int q = 0; q + y <= B; q++) for (int r = 0; r + z <= C; r++) {
                if (p + q + r == 0) continue; long u2 = id3(p, q, r), u = id3(x + p, y + q, z + r);
                int beta = Qv[u] - Qv[u1] - Qv[u2];
                int v = (beta > 0 ? beta : 0) + f[u1] + f[u2];
                pairs++;
                if (f[u] > v) viol++;
                if (f[u] == v) attained[u] = 1;
            }
        }
        long unatt = 0;
        for (int a = 0; a <= A; a++) for (int b = 0; b <= B; b++) for (int c = 0; c <= C; c++)
            if (a + b + c >= 2 && !attained[id3(a, b, c)]) unatt++;
        printf("VERIFY push-style Bellman inequality: %ld violations over %ld ordered pairs; %ld non-unit states not attained\n", viol, pairs, unatt);
    }
    int K = B < C ? B : C;
    for (int n = 0; n <= A; n++) for (int k = 0; k <= K && n + k <= A; k++)
        printf("F %d %d %d\n", n, k, f[id3(n + k, k, k)]);
    if (tree) {   /* print an optimal tree in the producer's MERGE format (for the audit's own certificate) */
        FILE *fo = fopen("audit_tree_out.txt", "w");
        long *st = malloc(sizeof(long) * 100000); int sp = 0; st[sp++] = id3(ta, tb, tc);
        while (sp) {
            long u = st[--sp]; long u1 = how[u]; if (u1 < 0) continue;
            int c = u % (C + 1), b = (u / (C + 1)) % (B + 1), a = u / ((long)(C + 1) * (B + 1));
            int z = u1 % (C + 1), y = (u1 / (C + 1)) % (B + 1), x = u1 / ((long)(C + 1) * (B + 1));
            fprintf(fo, "MERGE %d %d %d = %d %d %d + %d %d %d\n", a, b, c, x, y, z, a - x, b - y, c - z);
            st[sp++] = u1; st[sp++] = id3(a - x, b - y, c - z);
        }
        fclose(fo);
        printf("TREE written for (%d,%d,%d), value %d\n", ta, tb, tc, f[id3(ta, tb, tc)]);
    }
    return 0;
}
