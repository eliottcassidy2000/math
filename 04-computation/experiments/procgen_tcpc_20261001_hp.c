/*
 * procgen_tcpc_20261001_hp.c  (TCPC lane, collatz-procgen-20260922, 2026-10-01)
 *
 * Hamiltonian path / cycle counts of a digraph with arc MULTIPLICITIES
 * (loops are ignored: a Hamiltonian path never uses a loop).
 *
 * Input on stdin:  n  followed by the n x n matrix A (A[u][v] = multiplicity
 * of the arc u -> v, a small non-negative integer).
 *
 * Mode "exact" (argv[1] = "exact", n <= 18): unsigned 128-bit arithmetic.
 *   prints  H <H>            (number of Hamiltonian paths, weighted)
 *           HC <HC>          (number of directed Hamiltonian cycles, each cycle once)
 *           S v <start(v)>   (HPs starting at v)          for every v
 *           E v <end(v)>     (HPs ending at v)            for every v
 *           C u v <c(u->v)>  (HPs using the arc u -> v)   for every u != v with A[u][v] > 0
 * Mode "rooted" (n <= 20): prints  R <number of Hamiltonian paths starting at vertex 0>.
 * Mode "batchpar": K digraphs in one stream; prints "<i> <H mod 2> <#odd arcs> <#arcs>" per digraph.
 * Mode "par" (argv[1] = "par", n <= 26): everything modulo 2 (bitset DP).
 *   prints  H <H mod 2>, S/E/C lines with values mod 2 (no HC).
 *
 * Method: F[S][v] = weighted number of directed paths with vertex set S ending at v,
 * B[S][v] = the same starting at v (F of the transposed digraph);
 * c(u->v) = sum over S containing u, not v, of F[S][u] * A[u][v] * B[V\S][v].
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef unsigned __int128 u128;

static void print_u128(u128 x) {
    char buf[64]; int i = 63; buf[i] = 0;
    if (x == 0) { putchar('0'); return; }
    while (x > 0) { buf[--i] = '0' + (int)(x % 10); x /= 10; }
    fputs(buf + i, stdout);
}

static int n;
static int A[32][32];

static void run_exact(void) {
    if (n > 18) { fprintf(stderr, "exact mode needs n <= 18\n"); exit(2); }
    size_t N = (size_t)1 << n;
    u128 *F = calloc(N * n, sizeof(u128));
    u128 *B = calloc(N * n, sizeof(u128));
    if (!F || !B) { fprintf(stderr, "alloc failed\n"); exit(3); }
    for (int v = 0; v < n; v++) { F[((size_t)1 << v) * n + v] = 1; B[((size_t)1 << v) * n + v] = 1; }
    for (size_t S = 1; S < N; S++) {
        for (int v = 0; v < n; v++) {
            if (!(S >> v & 1)) continue;
            u128 f = F[S * n + v], b = B[S * n + v];
            if (f == 0 && b == 0) continue;
            for (int w = 0; w < n; w++) {
                if (S >> w & 1) continue;
                size_t T = S | ((size_t)1 << w);
                if (f && A[v][w]) F[T * n + w] += f * (u128)A[v][w];
                if (b && A[w][v]) B[T * n + w] += b * (u128)A[w][v];
            }
        }
    }
    size_t FULL = N - 1;
    u128 H = 0;
    for (int v = 0; v < n; v++) H += F[FULL * n + v];
    printf("H "); print_u128(H); printf("\n");
    /* Hamiltonian cycles: paths starting at vertex 0 */
    {
        u128 *G = calloc(N * n, sizeof(u128));
        if (!G) { fprintf(stderr, "alloc failed\n"); exit(3); }
        G[1 * n + 0] = 1;
        for (size_t S = 1; S < N; S += 2) {
            for (int v = 0; v < n; v++) {
                if (!(S >> v & 1)) continue;
                u128 g = G[S * n + v];
                if (!g) continue;
                for (int w = 0; w < n; w++) {
                    if (S >> w & 1) continue;
                    if (A[v][w]) G[(S | ((size_t)1 << w)) * n + w] += g * (u128)A[v][w];
                }
            }
        }
        u128 HC = 0;
        if (n == 1) HC = 0;
        else for (int v = 1; v < n; v++) HC += G[FULL * n + v] * (u128)A[v][0];
        printf("HC "); print_u128(HC); printf("\n");
        free(G);
    }
    for (int v = 0; v < n; v++) { printf("S %d ", v); print_u128(B[FULL * n + v]); printf("\n"); }
    for (int v = 0; v < n; v++) { printf("E %d ", v); print_u128(F[FULL * n + v]); printf("\n"); }
    u128 *C = calloc((size_t)n * n, sizeof(u128));
    for (size_t S = 1; S < FULL; S++) {
        size_t R = FULL ^ S;
        for (int u = 0; u < n; u++) {
            if (!(S >> u & 1)) continue;
            u128 f = F[S * n + u];
            if (!f) continue;
            for (int v = 0; v < n; v++) {
                if (!(R >> v & 1) || !A[u][v]) continue;
                u128 b = B[R * n + v];
                if (b) C[u * n + v] += f * (u128)A[u][v] * b;
            }
        }
    }
    for (int u = 0; u < n; u++) for (int v = 0; v < n; v++) {
        if (u == v || !A[u][v]) continue;
        printf("C %d %d ", u, v); print_u128(C[u * n + v]); printf("\n");
    }
    free(C); free(F); free(B);
}

static void run_par(void) {
    if (n > 26) { fprintf(stderr, "par mode needs n <= 26\n"); exit(2); }
    size_t N = (size_t)1 << n;
    uint32_t *F = calloc(N, sizeof(uint32_t));
    uint32_t *B = calloc(N, sizeof(uint32_t));
    if (!F || !B) { fprintf(stderr, "alloc failed\n"); exit(3); }
    uint32_t inA[32], outA[32];
    for (int w = 0; w < n; w++) {
        inA[w] = 0; outA[w] = 0;
        for (int u = 0; u < n; u++) {
            if (u == w) continue;
            if (A[u][w] & 1) inA[w] |= 1u << u;
            if (A[w][u] & 1) outA[w] |= 1u << u;
        }
    }
    for (int v = 0; v < n; v++) { F[(size_t)1 << v] = 1u << v; B[(size_t)1 << v] = 1u << v; }
    for (size_t S = 1; S < N; S++) {
        uint32_t f = F[S], b = B[S];
        if (!f && !b) continue;
        for (int w = 0; w < n; w++) {
            if (S >> w & 1) continue;
            size_t T = S | ((size_t)1 << w);
            if (__builtin_popcount(f & inA[w]) & 1) F[T] ^= 1u << w;
            if (__builtin_popcount(b & outA[w]) & 1) B[T] ^= 1u << w;
        }
    }
    size_t FULL = N - 1;
    printf("H %d\n", __builtin_popcount(F[FULL]) & 1);
    for (int v = 0; v < n; v++) printf("S %d %u\n", v, (B[FULL] >> v) & 1u);
    for (int v = 0; v < n; v++) printf("E %d %u\n", v, (F[FULL] >> v) & 1u);
    uint32_t acc[32];
    memset(acc, 0, sizeof(acc));
    for (size_t S = 1; S < FULL; S++) {
        uint32_t f = F[S];
        if (!f) continue;
        uint32_t b = B[FULL ^ S];
        if (!b) continue;
        while (f) {
            int u = __builtin_ctz(f); f &= f - 1;
            acc[u] ^= b & outA[u];
        }
    }
    for (int u = 0; u < n; u++) for (int v = 0; v < n; v++) {
        if (u == v || !A[u][v]) continue;
        printf("C %d %d %u\n", u, v, (acc[u] >> v) & 1u);
    }
    free(F); free(B);
}


static void run_rooted(void) {
    /* number of Hamiltonian paths starting at vertex 0 (weighted), n <= 20 */
    if (n > 20) { fprintf(stderr, "rooted mode needs n <= 20\n"); exit(2); }
    size_t N = (size_t)1 << n;
    size_t half = N >> 1;            /* subsets containing vertex 0 are S = 2*T+1 */
    u128 *G = calloc(half * n, sizeof(u128));
    if (!G) { fprintf(stderr, "alloc failed\n"); exit(3); }
    G[0 * n + 0] = 1;                /* S = {0} -> index 0 */
    for (size_t t = 0; t < half; t++) {
        size_t S = 2 * t + 1;
        for (int v = 0; v < n; v++) {
            if (!(S >> v & 1)) continue;
            u128 g = G[t * n + v];
            if (!g) continue;
            for (int w = 1; w < n; w++) {
                if (S >> w & 1) continue;
                if (A[v][w]) G[(((S | ((size_t)1 << w)) - 1) >> 1) * n + w] += g * (u128)A[v][w];
            }
        }
    }
    u128 R = 0;
    for (int v = 0; v < n; v++) R += G[(half - 1) * n + v];
    printf("R "); print_u128(R); printf("\n");
    free(G);
}

/* batch parity mode: stdin = K, then K digraphs (n + matrix); output per digraph:
   "<index> <H mod 2> <#arcs with c odd> <#arcs>"  (simple digraphs, n <= 22) */
static int batch_one(int idx) {
    if (scanf("%d", &n) != 1 || n < 1 || n > 22) return 0;
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++)
        if (scanf("%d", &A[i][j]) != 1) return 0;
    size_t N = (size_t)1 << n;
    static uint32_t *F = NULL, *B = NULL;
    static size_t cap = 0;
    if (cap < N) { free(F); free(B); F = malloc(N * 4); B = malloc(N * 4); cap = N; }
    memset(F, 0, N * 4); memset(B, 0, N * 4);
    uint32_t inA[32], outA[32];
    for (int w = 0; w < n; w++) {
        inA[w] = 0; outA[w] = 0;
        for (int u = 0; u < n; u++) {
            if (u == w) continue;
            if (A[u][w] & 1) inA[w] |= 1u << u;
            if (A[w][u] & 1) outA[w] |= 1u << u;
        }
    }
    for (int v = 0; v < n; v++) { F[(size_t)1 << v] = 1u << v; B[(size_t)1 << v] = 1u << v; }
    for (size_t S = 1; S < N; S++) {
        uint32_t f = F[S], b = B[S];
        if (!f && !b) continue;
        for (int w = 0; w < n; w++) {
            if (S >> w & 1) continue;
            size_t T = S | ((size_t)1 << w);
            if (__builtin_popcount(f & inA[w]) & 1) F[T] ^= 1u << w;
            if (__builtin_popcount(b & outA[w]) & 1) B[T] ^= 1u << w;
        }
    }
    size_t FULL = N - 1;
    uint32_t acc[32]; memset(acc, 0, sizeof(acc));
    for (size_t S = 1; S < FULL; S++) {
        uint32_t f = F[S];
        if (!f) continue;
        uint32_t b = B[FULL ^ S];
        if (!b) continue;
        while (f) { int u = __builtin_ctz(f); f &= f - 1; acc[u] ^= b & outA[u]; }
    }
    int odd = 0, tot = 0;
    for (int u = 0; u < n; u++) for (int v = 0; v < n; v++) {
        if (u == v || !A[u][v]) continue;
        tot++; odd += (acc[u] >> v) & 1u;
    }
    printf("%d %d %d %d\n", idx, __builtin_popcount(F[FULL]) & 1, odd, tot);
    return 1;
}
static void run_batchpar(void) {
    int K;
    if (scanf("%d", &K) != 1) return;
    for (int t = 0; t < K; t++) if (!batch_one(t)) { fprintf(stderr, "bad batch item %d\n", t); exit(1); }
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s exact|par|rooted|batchpar < digraph\n", argv[0]); return 1; }
    if (!strcmp(argv[1], "batchpar")) { run_batchpar(); return 0; }
    if (scanf("%d", &n) != 1 || n < 1 || n > 26) { fprintf(stderr, "bad n\n"); return 1; }
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) {
        if (scanf("%d", &A[i][j]) != 1) { fprintf(stderr, "bad matrix\n"); return 1; }
        if (A[i][j] < 0) { fprintf(stderr, "negative multiplicity\n"); return 1; }
    }
    if (!strcmp(argv[1], "exact")) run_exact();
    else if (!strcmp(argv[1], "par")) run_par();
    else if (!strcmp(argv[1], "rooted")) run_rooted();
    else if (!strcmp(argv[1], "batchpar")) { run_batchpar(); return 0; }
    else { fprintf(stderr, "unknown mode\n"); return 1; }
    return 0;
}
