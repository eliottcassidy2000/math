/* hampath_counts_circulant.c
 * mac-mini-2026-06-10 (research subagent)
 *
 * Exact counts of Hamiltonian paths P(T) and Hamiltonian cycles C(T) via
 * Held-Karp bitmask DP, for two circulant tournament families:
 *   - R_n  : "consecutive" rotational tournament on Z_n (odd n),
 *            arc i->j iff (j-i) mod n in {1,...,(n-1)/2}
 *            (the family Alon/Wormald conjecture maximizes #Hamilton cycles)
 *   - QR_p : Paley/quadratic-residue tournament on F_p (p prime, p=3 mod 4),
 *            arc i->j iff (j-i) is a nonzero square mod p
 * Also: brute force over ALL labeled tournaments for n<=7 to find max #ham-paths
 * (checks OEIS A038375 and identifies maximizers).
 *
 * Context: comparison against A038375 = 1,1,3,5,15,45,189,661,3357,15745,95095.
 * Compile: cc -O2 -o hampath_counts hampath_counts_circulant.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static uint64_t *dp; /* dp[S*n+v] */

/* count Hamiltonian paths in tournament given by out-neighbor bitmasks */
static uint64_t count_ham_paths(int n, const uint32_t *out) {
    size_t N = (size_t)1 << n;
    memset(dp, 0, N * n * sizeof(uint64_t));
    for (int v = 0; v < n; v++) dp[((size_t)1 << v) * n + v] = 1;
    for (size_t S = 1; S < N; S++) {
        for (int v = 0; v < n; v++) {
            uint64_t c = dp[S * n + v];
            if (!c || !((S >> v) & 1)) continue;
            uint32_t ext = out[v] & ~(uint32_t)S;
            while (ext) {
                int u = __builtin_ctz(ext);
                ext &= ext - 1;
                dp[(S | ((size_t)1 << u)) * n + u] += c;
            }
        }
    }
    uint64_t tot = 0;
    for (int v = 0; v < n; v++) tot += dp[(N - 1) * n + v];
    return tot;
}

/* count Hamiltonian cycles: paths starting at 0, closing arc v->0 */
static uint64_t count_ham_cycles(int n, const uint32_t *out) {
    size_t N = (size_t)1 << n;
    memset(dp, 0, N * n * sizeof(uint64_t));
    dp[1 * (size_t)n + 0] = 1;
    for (size_t S = 1; S < N; S++) {
        if (!(S & 1)) continue;
        for (int v = 0; v < n; v++) {
            uint64_t c = dp[S * n + v];
            if (!c || !((S >> v) & 1)) continue;
            uint32_t ext = out[v] & ~(uint32_t)S;
            while (ext) {
                int u = __builtin_ctz(ext);
                ext &= ext - 1;
                dp[(S | ((size_t)1 << u)) * n + u] += c;
            }
        }
    }
    uint64_t tot = 0;
    for (int v = 1; v < n; v++)
        if ((out[v] >> 0) & 1) tot += dp[(N - 1) * n + v];
    return tot;
}

static void build_consecutive(int n, uint32_t *out) {
    for (int i = 0; i < n; i++) {
        out[i] = 0;
        for (int r = 1; r <= (n - 1) / 2; r++) out[i] |= 1u << ((i + r) % n);
    }
}

static int build_paley(int p, uint32_t *out) {
    if (p % 4 != 3) return 0;
    int *isqr = calloc(p, sizeof(int));
    for (int i = 1; i < p; i++) isqr[(int)(((long long)i * i) % p)] = 1;
    for (int i = 0; i < p; i++) {
        out[i] = 0;
        for (int j = 0; j < p; j++)
            if (j != i && isqr[((j - i) % p + p) % p]) out[i] |= 1u << j;
    }
    free(isqr);
    return 1;
}

/* brute force all labeled tournaments on n vertices (n<=7), max ham paths */
static void brute_force_max(int n) {
    int m = n * (n - 1) / 2;
    int ei[32], ej[32], k = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++) { ei[k] = i; ej[k] = j; k++; }
    uint64_t best = 0, count_best = 0;
    uint32_t out[8];
    uint32_t best_example = 0;
    for (uint32_t mask = 0; mask < (1u << m); mask++) {
        for (int i = 0; i < n; i++) out[i] = 0;
        for (int e = 0; e < m; e++) {
            if ((mask >> e) & 1) out[ei[e]] |= 1u << ej[e];
            else out[ej[e]] |= 1u << ei[e];
        }
        uint64_t c = count_ham_paths(n, out);
        if (c > best) { best = c; count_best = 1; best_example = mask; }
        else if (c == best) count_best++;
    }
    printf("BRUTE n=%d : max #hampaths = %llu ; achieved by %llu labeled tournaments ; example mask=%u\n",
           n, (unsigned long long)best, (unsigned long long)count_best, best_example);
    /* second pass: max hamiltonian cycles */
    best = 0; count_best = 0; best_example = 0;
    for (uint32_t mask = 0; mask < (1u << m); mask++) {
        for (int i = 0; i < n; i++) out[i] = 0;
        for (int e = 0; e < m; e++) {
            if ((mask >> e) & 1) out[ei[e]] |= 1u << ej[e];
            else out[ej[e]] |= 1u << ei[e];
        }
        uint64_t c = count_ham_cycles(n, out);
        if (c > best) { best = c; count_best = 1; best_example = mask; }
        else if (c == best) count_best++;
    }
    printf("BRUTE n=%d : max #hamcycles = %llu ; achieved by %llu labeled tournaments ; example mask=%u\n",
           n, (unsigned long long)best, (unsigned long long)count_best, best_example);
}

int main(int argc, char **argv) {
    int NMAX = 23;
    if (argc > 1) NMAX = atoi(argv[1]);
    dp = malloc(((size_t)1 << NMAX) * NMAX * sizeof(uint64_t));
    if (!dp) { fprintf(stderr, "alloc failed\n"); return 1; }
    uint32_t out[32];

    printf("== consecutive rotational tournaments R_n (arc i->i+r, 1<=r<=(n-1)/2) ==\n");
    for (int n = 3; n <= NMAX; n += 2) {
        build_consecutive(n, out);
        uint64_t hp = count_ham_paths(n, out);
        uint64_t hc = count_ham_cycles(n, out);
        printf("R_%-2d : hampaths = %20llu   hamcycles = %20llu\n",
               n, (unsigned long long)hp, (unsigned long long)hc);
        fflush(stdout);
    }

    printf("== Paley (quadratic residue) tournaments QR_p, p = 3 mod 4 ==\n");
    int primes[] = {3, 7, 11, 19, 23};
    for (int t = 0; t < 5; t++) {
        int p = primes[t];
        if (p > NMAX) break;
        if (!build_paley(p, out)) continue;
        uint64_t hp = count_ham_paths(p, out);
        uint64_t hc = count_ham_cycles(p, out);
        printf("QR_%-2d: hampaths = %20llu   hamcycles = %20llu\n",
               p, (unsigned long long)hp, (unsigned long long)hc);
        fflush(stdout);
    }

    printf("== other 9-vertex circulants (connection sets) ==\n");
    int sets9[][4] = {{1,2,3,4},{1,2,3,5},{1,2,4,5},{1,3,4,5},{2,3,4,5},{1,2,4,8},{1,2,3,7}};
    for (int t = 0; t < 7; t++) {
        int n = 9; int ok = 1;
        char used[9]; memset(used,0,9);
        for (int i = 0; i < n; i++) out[i] = 0;
        for (int s = 0; s < 4; s++) {
            int r = sets9[t][s] % 9; int rr = (9 - r) % 9;
            if (used[r] || used[rr] || r == 0) { ok = 0; break; }
            used[r] = 1;
            for (int i = 0; i < n; i++) out[i] |= 1u << ((i + r) % n);
        }
        if (!ok) continue;
        uint64_t hp = count_ham_paths(n, out);
        uint64_t hc = count_ham_cycles(n, out);
        printf("Z9{%d,%d,%d,%d}: hampaths = %12llu hamcycles = %12llu\n",
               sets9[t][0],sets9[t][1],sets9[t][2],sets9[t][3],
               (unsigned long long)hp, (unsigned long long)hc);
    }

    if (NMAX >= 7) { brute_force_max(5); brute_force_max(6); brute_force_max(7); }
    free(dp);
    return 0;
}
