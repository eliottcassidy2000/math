/* ready_search.c (klein-2026-11-04) -- strong "ready" tournaments for the two-branch C2 grammar.
   A tournament T on n vertices with root r and auxiliary a is sigma-ready iff N^+(a) = {r} and e_r(T - a) = sigma,
   where e_r counts Hamiltonian paths of T - a ending at r (JOINT_CLASSIFICATION, nonzero_preparation_20261103).
   modes:
     search n sigma target seed steps  : simulated annealing over the free arcs (root = 0, auxiliary = 1 fixed:
                                         1 -> 0 and v -> 1 for all v != 0,1) for a STRONG sigma-ready T with H = target.
                                         Prints the witness as out-row bit encodings.
     floor n                           : exhaustive over gentourng output on stdin: for sigma = 1, 2 the minimal H of a
                                         strong tournament having some pair (a, r) with N^+(a) = {r}, e_r(T-a) = sigma;
                                         also counts such classes.
     check n rows...                   : H, strongness, and all sole-exit auxiliaries with their e values.         */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
typedef unsigned long long u64;
static int n; static u64 adj[64]; static u64 *f;
static u64 count_paths_sub(u64 verts, int end) {  /* Hamiltonian paths of T[verts]; if end >= 0 only those ending at end */
    int m = __builtin_popcountll(verts); int map[64], k = 0;
    for (int v = 0; v < n; v++) if (verts >> v & 1) map[k++] = v;
    u64 N = 1ULL << m; memset(f, 0, sizeof(u64) * N * m);
    for (int i = 0; i < m; i++) f[(1ULL << i) * m + i] = 1;
    for (u64 S = 1; S < N; S++) for (int i = 0; i < m; i++) {
        u64 c = f[S * m + i]; if (!c) continue;
        u64 out = adj[map[i]];
        for (int j = 0; j < m; j++) if (!(S >> j & 1) && (out >> map[j] & 1)) f[(S | (1ULL << j)) * m + j] += c;
    }
    u64 h = 0; for (int i = 0; i < m; i++) if (end < 0 || map[i] == end) h += f[(N - 1) * m + i]; return h;
}
static int is_strong(void) {
    u64 all = (1ULL << n) - 1, seen = 1, fr = 1;
    while (fr) { u64 nf = 0; for (int v = 0; v < n; v++) if (fr >> v & 1) nf |= adj[v]; nf &= ~seen; seen |= nf; fr = nf; }
    if (seen != all) return 0;
    u64 radj[64]; for (int v = 0; v < n; v++) radj[v] = 0;
    for (int u = 0; u < n; u++) for (int v = 0; v < n; v++) if (adj[u] >> v & 1) radj[v] |= 1ULL << u;
    seen = 1; fr = 1;
    while (fr) { u64 nf = 0; for (int v = 0; v < n; v++) if (fr >> v & 1) nf |= radj[v]; nf &= ~seen; seen |= nf; fr = nf; }
    return seen == all;
}
static u64 xs(u64 *s) { u64 x = *s; x ^= x << 13; x ^= x >> 7; x ^= x << 17; *s = x; return x; }
static void flip(int i, int j) { if (adj[i] >> j & 1) { adj[i] &= ~(1ULL << j); adj[j] |= 1ULL << i; } else { adj[j] &= ~(1ULL << i); adj[i] |= 1ULL << j; } }
int main(int argc, char **argv) { setvbuf(stdout, NULL, _IOLBF, 0);
    if (argc < 3) { fprintf(stderr, "usage\n"); return 1; }
    n = atoi(argv[2]); f = malloc(sizeof(u64) * (1ULL << n) * n);
    u64 all = (1ULL << n) - 1;
    if (!strcmp(argv[1], "check")) {
        for (int i = 0; i < n; i++) adj[i] = strtoull(argv[3 + i], 0, 10);
        printf("H = %llu, strong = %d\n", count_paths_sub(all, -1), is_strong());
        for (int a = 0; a < n; a++) if (__builtin_popcountll(adj[a]) == 1) { int r = __builtin_ctzll(adj[a]); printf("  aux %d -> root %d: e_r(T-a) = %llu\n", a, r, count_paths_sub(all & ~(1ULL << a), r)); }
        return 0;
    }
    if (!strcmp(argv[1], "floor")) {
        char line[4096]; u64 minH[3] = {0, 0, 0}, cnt[3] = {0, 0, 0}, total = 0; char best[3][4096];
        while (fgets(line, sizeof line, stdin)) {
            int L = strlen(line); while (L && (line[L-1] == '\n' || line[L-1] == '\r')) line[--L] = 0;
            if (L != n * (n - 1) / 2) continue;
            memset(adj, 0, sizeof adj); int idx = 0;
            for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (line[idx++] == '1') adj[i] |= 1ULL << j; else adj[j] |= 1ULL << i; }
            total++; if (!is_strong()) continue;
            int ready[3] = {0, 0, 0};
            for (int a = 0; a < n; a++) if (__builtin_popcountll(adj[a]) == 1) { int r = __builtin_ctzll(adj[a]); u64 e = count_paths_sub(all & ~(1ULL << a), r); if (e == 1 || e == 2) ready[e] = 1; }
            if (ready[1] || ready[2]) { u64 h = count_paths_sub(all, -1);
                for (int s = 1; s <= 2; s++) if (ready[s]) { cnt[s]++; if (!minH[s] || h < minH[s]) { minH[s] = h; strcpy(best[s], line); } } }
        }
        printf("order %d: %llu tournaments; sigma=1-ready strong classes %llu, min H %llu (%s); sigma=2-ready %llu, min H %llu (%s)\n", n, total, cnt[1], minH[1], best[1], cnt[2], minH[2], best[2]);
        return 0;
    }
    if (!strcmp(argv[1], "pairs")) {   /* every strong tournament on stdin: for each sole-exit vertex a (out-degree 1) print H sigma */
        char line[4096];
        while (fgets(line, sizeof line, stdin)) {
            int L = strlen(line); while (L && (line[L-1] == '\n' || line[L-1] == '\r')) line[--L] = 0;
            if (L != n * (n - 1) / 2) continue;
            memset(adj, 0, sizeof adj); int idx = 0;
            for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (line[idx++] == '1') adj[i] |= 1ULL << j; else adj[j] |= 1ULL << i; }
            if (!is_strong()) continue;
            u64 h = 0; int have = 0;
            for (int a = 0; a < n; a++) if (__builtin_popcountll(adj[a]) == 1) {
                if (!have) { h = count_paths_sub(all, -1); have = 1; }
                int r = __builtin_ctzll(adj[a]); u64 e = count_paths_sub(all & ~(1ULL << a), r);
                printf("%llu %llu\n", h, e);
            }
        }
        return 0;
    }
    if (!strcmp(argv[1], "csearch")) {   /* constrained walk: start from E_{n-1} + a (sigma = 1 ready), keep e = sigma and strongness as HARD constraints */
        int sigma = atoi(argv[3]); double target = atof(argv[4]); u64 seed = (u64)atoll(argv[5]) * 2654435761ULL + 7; long steps = atol(argv[6]);
        /* vertices: 0 = root r, 1 = auxiliary a; U = T - a on {0,2,...,n-1}: build E_{n-1} on U with r = 0 the "vertex 1" of E (beaten only by the last vertex) */
        memset(adj, 0, sizeof adj);
        int U[64], m = 0; U[m++] = 0; for (int v = 2; v < n; v++) U[m++] = v;
        for (int i = 0; i < m; i++) for (int j = i + 1; j < m; j++) adj[U[i]] |= 1ULL << U[j];   /* transitive on U in the order U[0..m-1] */
        adj[U[0]] &= ~(1ULL << U[m-1]); adj[U[m-1]] |= 1ULL << U[0];                        /* reverse the extreme arc: last -> first */
        adj[1] = 1ULL; for (int v = 2; v < n; v++) adj[v] |= 1ULL << 1;                         /* a -> r only; everyone else -> a */
        u64 H = count_paths_sub(all, -1), E = count_paths_sub(all & ~2ULL, 0); int st = is_strong();
        fprintf(stderr, "start: H=%llu e=%llu strong=%d\n", H, E, st);
        /* if sigma = 2 is wanted, first random-walk until e == sigma (accepting any strong state), then constrain */
        double cost = fabs(log((double)H / target)); long it; double T0 = 1.0, T1 = 0.01;
        for (it = 0; it < steps; it++) {
            int i = xs(&seed) % n, j = xs(&seed) % n; if (i == j || i == 1 || j == 1) continue;
            flip(i, j);
            u64 E2 = count_paths_sub(all & ~2ULL, 0); int st2 = is_strong();
            int feasible_now = (E == (u64)sigma);
            if (!st2 || (feasible_now && E2 != (u64)sigma) || (!feasible_now && (E2 > (u64)sigma + 3))) { flip(i, j); continue; }
            u64 H2 = count_paths_sub(all, -1);
            double c2 = fabs(log((double)H2 / target)) + (E2 == (u64)sigma ? 0 : 0.5);
            double temp = T0 * pow(T1 / T0, (double)it / steps);
            if (c2 <= cost || exp((cost - c2) / temp) > (double)(xs(&seed) % 1000000) / 1e6) { cost = c2; H = H2; E = E2; st = st2; }
            else flip(i, j);
            if (H == (u64)target && E == (u64)sigma && st) {
                printf("FOUND n=%d sigma=%d H=%llu after %ld steps; rows:", n, sigma, H, it);
                for (int v = 0; v < n; v++) printf(" %llu", adj[v]); printf("\n"); return 0;
            }
            if (it % 20000 == 0 && it) fprintf(stderr, "  it %ld: H=%llu e=%llu strong=%d cost=%.4f\n", it, H, E, st, cost);
        }
        printf("not found n=%d sigma=%d target=%.0f; last H=%llu e=%llu\n", n, sigma, target, H, E); return 1;
    }
    if (!strcmp(argv[1], "births")) {   /* births n a rows...: first birth x (x -> a, others -> x); then every pattern for y; report children's H and readiness */
        int a = atoi(argv[3]); for (int i = 0; i < n; i++) adj[i] = strtoull(argv[4 + i], 0, 10);
        int n0 = n; u64 base[64]; memcpy(base, adj, sizeof adj);
        int x = n0, y = n0 + 1; n = n0 + 2; free(f); f = malloc(sizeof(u64) * (1ULL << n) * n); all = (1ULL << n) - 1;
        u64 H0; { n = n0; u64 a0 = (1ULL << n0) - 1; H0 = count_paths_sub(a0, -1); n = n0 + 2; }
        printf("input order %d H=%llu; scanning %llu y-patterns\n", n0, H0, 1ULL << (n0 + 1));
        for (u64 pat = 0; pat < (1ULL << (n0 + 1)); pat++) {   /* bit v of pat: v beats y (v in 0..n0, including x = n0) */
            memcpy(adj, base, sizeof(u64) * n0); adj[x] = 1ULL << a; adj[y] = 0;
            for (int v = 0; v < n0; v++) if (v != a) adj[v] |= 1ULL << x;
            for (int v = 0; v <= n0; v++) { if (pat >> v & 1) adj[v] |= 1ULL << y; else adj[y] |= 1ULL << v; }
            if (!is_strong()) continue;
            u64 h = count_paths_sub(all, -1);
            /* readiness: sole-exit vertices with e in {1,2} */
            char info[512]; info[0] = 0; int ready = 0;
            for (int v = 0; v < n; v++) if (__builtin_popcountll(adj[v]) == 1) { int r = __builtin_ctzll(adj[v]); u64 e = count_paths_sub(all & ~(1ULL << v), r); if (e == 1 || e == 2) { ready = 1; char b[64]; sprintf(b, " (aux %d root %d e=%llu)", v, r, e); strcat(info, b); } }
            if (ready) printf("pattern %llu: H=%llu%s\n", pat, h, info);
        }
        return 0;
    }
    if (!strcmp(argv[1], "csearch2")) {   /* csearch2 n target seed steps rows... : start from the given rows (a sigma=1 witness, aux 1 -> root 0),
                                             keep strongness and e in {1,2} as hard constraints, aim at e = 2 with H = target */
        double target = atof(argv[3]); u64 seed = (u64)atoll(argv[4]) * 2654435761ULL + 11; long steps = atol(argv[5]);
        for (int i = 0; i < n; i++) adj[i] = strtoull(argv[6 + i], 0, 10);
        u64 H = count_paths_sub(all, -1), E = count_paths_sub(all & ~2ULL, 0); int st = is_strong();
        fprintf(stderr, "start: H=%llu e=%llu strong=%d\n", H, E, st);
        double cost = fabs(log((double)H / target)) + (E == 2 ? 0 : 0.3); long it; double T0 = 0.6, T1 = 0.005;
        for (it = 0; it < steps; it++) {
            int i = xs(&seed) % n, j = xs(&seed) % n; if (i == j || i == 1 || j == 1) continue;
            flip(i, j);
            u64 E2 = count_paths_sub(all & ~2ULL, 0); int st2 = is_strong();
            if (!st2 || E2 < 1 || E2 > 2) { flip(i, j); continue; }
            u64 H2 = count_paths_sub(all, -1);
            double c2 = fabs(log((double)H2 / target)) + (E2 == 2 ? 0 : 0.3);
            double temp = T0 * pow(T1 / T0, (double)it / steps);
            if (c2 <= cost || exp((cost - c2) / temp) > (double)(xs(&seed) % 1000000) / 1e6) { cost = c2; H = H2; E = E2; st = st2; }
            else flip(i, j);
            if (H == (u64)target && E == 2 && st) {
                printf("FOUND n=%d sigma=2 H=%llu after %ld steps; rows:", n, H, it);
                for (int v = 0; v < n; v++) printf(" %llu", adj[v]); printf("\n"); return 0;
            }
            if (it % 20000 == 0 && it) fprintf(stderr, "  it %ld: H=%llu e=%llu cost=%.4f\n", it, H, E, cost);
        }
        printf("not found sigma=2 target=%.0f; last H=%llu e=%llu\n", target, H, E); return 1;
    }
    if (!strcmp(argv[1], "joincyc")) {   /* per strong tournament on stdin: H, number of directed Hamiltonian cycles C, and the list of digits sigma over sole-exit vertices */
        char line[4096];
        while (fgets(line, sizeof line, stdin)) {
            int L = strlen(line); while (L && (line[L-1] == '\n' || line[L-1] == '\r')) line[--L] = 0;
            if (L != n * (n - 1) / 2) continue;
            memset(adj, 0, sizeof adj); int idx = 0;
            for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (line[idx++] == '1') adj[i] |= 1ULL << j; else adj[j] |= 1ULL << i; }
            if (!is_strong()) continue;
            u64 h = count_paths_sub(all, -1);
            /* cycles: paths of T - 0 starting at an out-neighbour of 0 and ending at an in-neighbour of 0: count = sum over v in N+(0) of paths of T-0 from v ending in N-(0) */
            u64 cyc = 0; u64 rest = all & ~1ULL;
            { int m = n - 1; int map[64], k = 0; for (int v = 1; v < n; v++) map[k++] = v; u64 N = 1ULL << m; memset(f, 0, sizeof(u64) * N * m);
              for (int i = 0; i < m; i++) if (adj[0] >> map[i] & 1) f[(1ULL << i) * m + i] = 1;
              for (u64 S = 1; S < N; S++) for (int i = 0; i < m; i++) { u64 c = f[S * m + i]; if (!c) continue; u64 out = adj[map[i]]; for (int j = 0; j < m; j++) if (!(S >> j & 1) && (out >> map[j] & 1)) f[(S | (1ULL << j)) * m + j] += c; }
              for (int i = 0; i < m; i++) if (adj[map[i]] >> 0 & 1) cyc += f[(N - 1) * m + i]; }
            printf("%llu %llu", h, cyc);
            for (int a = 0; a < n; a++) if (__builtin_popcountll(adj[a]) == 1) { int r = __builtin_ctzll(adj[a]); printf(" %llu", count_paths_sub(all & ~(1ULL << a), r)); }
            printf("\n");
        }
        return 0;
    }
    if (!strcmp(argv[1], "sample")) {   /* sample n count seed: random tournaments with vertex 1 a sole-exit into 0; keep strong; print H sigma */
        long count = atol(argv[3]); u64 seed = (u64)atoll(argv[4]) * 2654435761ULL + 3;
        for (long t = 0; t < count; t++) {
            memset(adj, 0, sizeof adj);
            for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (xs(&seed) & 1) adj[i] |= 1ULL << j; else adj[j] |= 1ULL << i; }
            adj[1] = 1ULL; for (int v = 2; v < n; v++) adj[v] |= 1ULL << 1; adj[0] &= ~(1ULL << 1);
            if (!is_strong()) continue;
            u64 e = count_paths_sub(all & ~2ULL, 0); u64 h = count_paths_sub(all, -1);
            printf("%llu %llu\n", h, e);
        }
        return 0;
    }
    if (!strcmp(argv[1], "search")) {
        int sigma = atoi(argv[3]); double target = atof(argv[4]); u64 seed = (u64)atoll(argv[5]) * 2654435761ULL + 7; long steps = atol(argv[6]);
        /* initial random tournament with the fixed auxiliary structure: vertex 1 beats only vertex 0 */
        memset(adj, 0, sizeof adj);
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (xs(&seed) & 1) adj[i] |= 1ULL << j; else adj[j] |= 1ULL << i; }
        adj[1] = 1ULL; for (int v = 2; v < n; v++) adj[v] |= 1ULL << 1; adj[0] &= ~(1ULL << 1);
        double T0 = 2.0, T1 = 0.02; double best = 1e18; long it;
        u64 H = count_paths_sub(all, -1), E = count_paths_sub(all & ~2ULL, 0); int st = is_strong();
        double cost = fabs(log((double)H / target)) * 10 + fabs((double)E - sigma) * 3 + (st ? 0 : 5);
        for (it = 0; it < steps; it++) {
            int i = xs(&seed) % n, j = xs(&seed) % n; if (i == j || i == 1 || j == 1) continue;   /* keep the auxiliary's arcs */
            flip(i, j);
            u64 H2 = count_paths_sub(all, -1), E2 = count_paths_sub(all & ~2ULL, 0); int st2 = is_strong();
            double c2 = fabs(log((double)H2 / target)) * 10 + fabs((double)E2 - sigma) * 3 + (st2 ? 0 : 5);
            double temp = T0 * pow(T1 / T0, (double)it / steps);
            if (c2 <= cost || exp((cost - c2) / temp) > (double)(xs(&seed) % 1000000) / 1e6) { cost = c2; H = H2; E = E2; st = st2; }
            else flip(i, j);
            if (cost < best) best = cost;
            if (H == (u64)target && E == (u64)sigma && st) {
                printf("FOUND n=%d sigma=%d H=%llu after %ld steps; rows:", n, sigma, H, it);
                for (int v = 0; v < n; v++) printf(" %llu", adj[v]); printf("\n"); return 0;
            }
            if (it % 20000 == 0 && it) fprintf(stderr, "  it %ld: H=%llu e=%llu strong=%d cost=%.3f best=%.3f\n", it, H, E, st, cost, best);
        }
        printf("not found n=%d sigma=%d target=%.0f after %ld steps; last H=%llu e=%llu strong=%d\n", n, sigma, target, steps, H, E, st);
        return 1;
    }
    return 1;
}
