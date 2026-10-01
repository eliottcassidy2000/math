/* procgen_petersen_20261001_hp.c
   Directed Hamiltonian-path counts with per-arc counts, for orientations of a simple graph
   (or for a given digraph, via mode L with mask 0 and the arcs listed as edges).

   stdin:  n m
           m lines "u v"              (edge k joins u_k and v_k)
           mode:
             A                         enumerate all 2^m orientations; print aggregate statistics
                                       and every orientation whose arcs all lie on an odd number of HPs
             L k mask_1 ... mask_k     evaluate the listed orientations; print "mask H c_0 ... c_{m-1}"
   Orientation bit k = 0 means the arc u_k -> v_k, bit k = 1 means v_k -> u_k (only the first 64 edges can be
   flipped by a mask; m <= 1024 edges in mode L, m <= 40 in mode A).
   c_k = number of directed Hamiltonian paths using edge k in its chosen direction.
   Method: subset DP, f[S][v] = #directed paths with vertex set S ending at v, g[S][v] = starting at v,
   c(a->b) = sum over S (a in S, b not in S) of f[S][a] * g[V \ S][b]. n <= 20. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXM 1024
static int n, m, eu[MAXM], ev[MAXM];
static int64_t *f, *g;
static int outa[32][32];

static void build(uint64_t mask) {
    memset(outa, 0, sizeof(outa));
    for (int k = 0; k < m; k++) {
        int bit = (k < 64) ? (int)((mask >> k) & 1) : 0;   /* edges beyond 64 keep the listed direction */
        if (bit) outa[ev[k]][eu[k]] = 1; else outa[eu[k]][ev[k]] = 1;
    }
}

static int64_t count(int64_t *c) {
    int N = 1 << n;
    memset(f, 0, sizeof(int64_t) * (size_t)N * n);
    memset(g, 0, sizeof(int64_t) * (size_t)N * n);
    for (int v = 0; v < n; v++) { f[(size_t)(1 << v) * n + v] = 1; g[(size_t)(1 << v) * n + v] = 1; }
    for (int S = 1; S < N; S++) {
        for (int v = 0; v < n; v++) {
            if (!((S >> v) & 1)) continue;
            int64_t fv = f[(size_t)S * n + v], gv = g[(size_t)S * n + v];
            if (!fv && !gv) continue;
            for (int w = 0; w < n; w++) {
                if ((S >> w) & 1) continue;
                if (fv && outa[v][w]) f[(size_t)(S | (1 << w)) * n + w] += fv;
                if (gv && outa[w][v]) g[(size_t)(S | (1 << w)) * n + w] += gv;
            }
        }
    }
    int64_t H = 0;
    for (int v = 0; v < n; v++) H += f[(size_t)(N - 1) * n + v];
    for (int k = 0; k < m; k++) {
        int a = eu[k], b = ev[k];
        if (!outa[a][b]) { int t = a; a = b; b = t; }
        int64_t s = 0;
        for (int S = 1; S < N; S++) {
            if (!((S >> a) & 1) || ((S >> b) & 1)) continue;
            int64_t x = f[(size_t)S * n + a];
            if (!x) continue;
            s += x * g[(size_t)((N - 1) ^ S) * n + b];
        }
        c[k] = s;
    }
    return H;
}

int main(void) {
    if (scanf("%d %d", &n, &m) != 2 || n < 1 || n > 20 || m < 0 || m > MAXM) return 1;
    for (int k = 0; k < m; k++) if (scanf("%d %d", &eu[k], &ev[k]) != 2) return 1;
    char mode[8];
    if (scanf("%7s", mode) != 1) return 1;
    size_t N = (size_t)1 << n;
    f = malloc(sizeof(int64_t) * N * n);
    g = malloc(sizeof(int64_t) * N * n);
    if (!f || !g) return 2;
    static int64_t c[MAXM];
    if (mode[0] == 'A') {
        if (m > 40) return 3;
        uint64_t total = (uint64_t)1 << m;
        int64_t Hodd = 0, allodd = 0, alleven = 0, alleven_Hodd = 0, Hmax = 0, Hsum = 0;
        static int64_t nodd_hist[MAXM + 1]; memset(nodd_hist, 0, sizeof(nodd_hist));
        for (uint64_t mask = 0; mask < total; mask++) {
            build(mask);
            int64_t H = count(c);
            Hsum += H;
            if (H > Hmax) Hmax = H;
            int nodd = 0;
            for (int k = 0; k < m; k++) nodd += (int)(c[k] & 1);
            nodd_hist[nodd]++;
            if (H & 1) Hodd++;
            if (nodd == m) { allodd++; printf("ALLODD %llu %lld\n", (unsigned long long)mask, (long long)H); }
            if (nodd == 0) { alleven++; if (H & 1) alleven_Hodd++; }
        }
        printf("TOTAL %llu HSUM %lld HMAX %lld HODD %lld ALLODD %lld ALLEVEN %lld ALLEVEN_HODD %lld\n",
               (unsigned long long)total, (long long)Hsum, (long long)Hmax, (long long)Hodd, (long long)allodd,
               (long long)alleven, (long long)alleven_Hodd);
        printf("NODD_HIST");
        for (int i = 0; i <= m; i++) printf(" %lld", (long long)nodd_hist[i]);
        printf("\n");
    } else {
        int k0;
        if (scanf("%d", &k0) != 1) return 1;
        for (int i = 0; i < k0; i++) {
            unsigned long long mask;
            if (scanf("%llu", &mask) != 1) return 1;
            build(mask);
            int64_t H = count(c);
            printf("%llu %lld", mask, (long long)H);
            for (int k = 0; k < m; k++) printf(" %lld", (long long)c[k]);
            printf("\n");
        }
    }
    free(f); free(g);
    return 0;
}
