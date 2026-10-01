/* procgen_selfie_20261001_paley2.c  (procgen selfie lane, 2026-10-01)
 *
 * Parity of c(e) = #{Hamiltonian paths through arc e} for the Paley tournament QR_q minus the vertex 0,
 * q = 3 mod 4 a prime or q = 27 (GF(27) = GF(3)[x]/(x^3 - x - 1)).  Mod-2 bitmask DP:
 *   f[S] = bitmask over v of (#paths covering S ending at v) mod 2      (one uint32 per subset S)
 *   g(S, v) = f(-S, -v)   because x -> -x is an anti-automorphism of QR_q - 0  (q = 3 mod 4).
 * c(u->v) mod 2 = sum over S containing u, not v, of f(S,u) g(S^c, v).
 * Memory: 4 * 2^(q-1) bytes (256 MB at q = 27).  Prints one summary line to stdout.
 * usage: procgen_selfie_paley2 q
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

static int q, N;
static int elt[32];            /* vertex index -> field element code */
static int idx_of[32];         /* field element code -> vertex index (or -1 for 0) */
static uint32_t outm[32], inm[32];
static int negidx[32];

static int add3(int a, int b) { int r = 0, m = 1; for (int k = 0; k < 3; k++) { r += (((a % 3) + (b % 3)) % 3) * m; a /= 3; b /= 3; m *= 3; } return r; }
static int neg3(int a) { int r = 0, m = 1; for (int k = 0; k < 3; k++) { r += ((3 - a % 3) % 3) * m; a /= 3; m *= 3; } return r; }
static int mul27(int a, int b) {
    int A[3], B[3], C[5] = {0};
    for (int k = 0; k < 3; k++) { A[k] = a % 3; a /= 3; B[k] = b % 3; b /= 3; }
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) C[i + j] += A[i] * B[j];
    /* x^3 = x + 1, x^4 = x^2 + x */
    C[1] += C[3]; C[0] += C[3]; C[2] += C[4]; C[1] += C[4];
    return (C[0] % 3) + 3 * (C[1] % 3) + 9 * (C[2] % 3);
}

int main(int argc, char **argv) {
    q = atoi(argv[1]); N = q - 1;
    int sq[32] = {0}, issq[32] = {0};
    for (int a = 1; a < q; a++) { int s = (q == 27) ? mul27(a, a) : (a * a) % q; issq[s] = 1; }
    (void)sq;
    for (int i = 0; i < 32; i++) idx_of[i] = -1;
    for (int i = 0; i < N; i++) { elt[i] = i + 1; idx_of[i + 1] = i; }
    for (int i = 0; i < N; i++) {
        outm[i] = inm[i] = 0;
        int ni = (q == 27) ? neg3(elt[i]) : (q - elt[i]) % q;
        negidx[i] = idx_of[ni];
    }
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) if (i != j) {
        int d = (q == 27) ? add3(elt[j], neg3(elt[i])) : ((elt[j] - elt[i]) % q + q) % q;
        if (issq[d]) { outm[i] |= 1u << j; inm[j] |= 1u << i; }
    }
    /* sanity: tournament */
    for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++) {
        int a = (outm[i] >> j) & 1, b = (outm[j] >> i) & 1;
        if (a + b != 1) { printf("NOT A TOURNAMENT q=%d\n", q); return 1; }
    }
    size_t M = (size_t)1 << N;
    uint32_t *f = calloc(M, sizeof(uint32_t));
    if (!f) { printf("alloc failed\n"); return 1; }
    for (int v = 0; v < N; v++) f[(size_t)1 << v] = 1u << v;
    for (size_t S = 1; S < M; S++) {
        uint32_t F = f[S];
        if (!F) continue;
        uint32_t rest = (uint32_t)(~S) & (uint32_t)(M - 1);
        while (rest) {
            int w = __builtin_ctz(rest); rest &= rest - 1;
            if (__builtin_popcount(F & inm[w]) & 1) f[S | ((size_t)1 << w)] ^= 1u << w;
        }
    }
    size_t full = M - 1;
    int Hpar = __builtin_popcount(f[full]) & 1;
    /* c parity */
    static uint8_t cpar[32][32];
    memset(cpar, 0, sizeof cpar);
    for (size_t S = 1; S < full; S++) {
        uint32_t F = f[S];
        if (!F) continue;
        size_t Cm = full ^ S;
        /* G = starts of paths covering Cm : g(Cm, v) = f(-Cm, -v) */
        size_t negC = 0; size_t t = Cm;
        while (t) { int b = __builtin_ctzll(t); t &= t - 1; negC |= (size_t)1 << negidx[b]; }
        uint32_t Fn = f[negC], G = 0;
        while (Fn) { int b = __builtin_ctz(Fn); Fn &= Fn - 1; G |= 1u << negidx[b]; }
        if (!G) continue;
        uint32_t FF = F;
        while (FF) {
            int u = __builtin_ctz(FF); FF &= FF - 1;
            uint32_t tg = G & outm[u];
            while (tg) { int v = __builtin_ctz(tg); tg &= tg - 1; cpar[u][v] ^= 1; }
        }
    }
    int odd = 0, arcs = 0, oddanti = 0, anti = 0;
    for (int u = 0; u < N; u++) for (int v = 0; v < N; v++) if ((outm[u] >> v) & 1) {
        arcs++; odd += cpar[u][v];
        if (negidx[u] == v) { anti++; oddanti += cpar[u][v]; }
    }
    printf("PALEY_MINUS_VERTEX q=%d N=%d arcs=%d odd_arcs=%d antipodal_arcs=%d odd_antipodal=%d H_parity=%d all_odd=%d\n",
           q, N, arcs, odd, anti, oddanti, Hpar, odd == arcs);
    free(f);
    return 0;
}
