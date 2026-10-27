/* audit_thomassen.c -- cross-check of the endpoint DP against Thomassen's characterisation (Proc. London Math. Soc.
 * 45 (1982) 151-168; Theorem 5.5 of Bang-Jensen & Gutin, "Paths, trees and cycles in tournaments"): for x of maximum
 * and y of minimum out-degree, T has a Hamiltonian (x,y)-path unless (i) the terminal strong component of T-x is one
 * vertex z with z->x and the preceding component is {y} with x->y, or (ii) the initial component of T-y is one vertex z
 * with y->z and the second component is {x} with x->y, or (iii) T is R(1,2) or one of two small exceptional tournaments.
 * Every (max, min) pair x != y of every class is tested.  Build: cc -O3 audit_thomassen.c -o audit_thomassen
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int n;
static uint32_t out[32], full;
static uint32_t *E;     /* E[S*n+v] */
static uint32_t Pi[32]; /* Pi[a] = set of ends b */

static void endpoint_relation(void) {
    size_t N = (size_t)1 << n;
    memset(E, 0, sizeof(uint32_t) * N * n);
    for (int v = 0; v < n; v++) E[((size_t)1 << v) * n + v] = 1u << v;
    for (size_t S = 1; S < N; S++) {
        for (int v = 0; v < n; v++) {
            uint32_t st = E[S * n + v];
            if (!st) continue;
            uint32_t nx = out[v] & ~(uint32_t)S & full;
            while (nx) { int w = __builtin_ctz(nx); nx &= nx - 1; E[(S | ((size_t)1 << w)) * n + w] |= st; }
        }
    }
    for (int a = 0; a < n; a++) Pi[a] = 0;
    for (int b = 0; b < n; b++) { uint32_t st = E[(N - 1) * n + b]; while (st) { int a = __builtin_ctz(st); st &= st - 1; Pi[a] |= 1u << b; } }
}

/* reachability inside vertex set vs from vertex v */
static uint32_t reach_in(uint32_t vs, int v) {
    uint32_t seen = 1u << v, fr = seen;
    while (fr) { uint32_t nx = 0; for (uint32_t f = fr; f; f &= f - 1) nx |= out[__builtin_ctz(f)] & vs; nx &= ~seen; seen |= nx; fr = nx; }
    return seen;
}
static int strong_in(uint32_t vs) {   /* induced subtournament on vs strongly connected (|vs| = 1 counts as strong) */
    if (!vs) return 1;
    int pc = __builtin_popcount(vs); if (pc == 1) return 1; if (pc == 2) return 0;
    for (uint32_t w = vs; w; w &= w - 1) if (reach_in(vs, __builtin_ctz(w)) != vs) return 0;
    return 1;
}
static int kappa_of(void) {    /* least |S| with T - S not strong (|T-S| >= 2); T assumed strong */
    for (int k = 1; k <= n - 2; k++) {
        for (uint32_t S = 1; S < full; S++) if (__builtin_popcount(S) == k) if (!strong_in(full & ~S)) return k;
    }
    return n - 1;
}
static int lambda_of(void) {
    int best = 1 << 20;
    for (uint32_t X = 1; X < full; X++) { int f = 0; for (int x = 0; x < n; x++) if (X >> x & 1) f += __builtin_popcount(out[x] & ~X & full); if (f < best) best = f; }
    return best;
}

/* ordered strong components of T[vs]: peel initial components */
static int comps_of(uint32_t vs, uint32_t *comp) {
    int k = 0;
    while (vs) {
        uint32_t init = 0;
        for (uint32_t w = vs; w; w &= w - 1) { int v = __builtin_ctz(w); if (reach_in(vs, v) == vs) init |= 1u << v; }
        comp[k++] = init; vs &= ~init;
    }
    return k;
}
int main(int argc, char **argv) {
    n = atoi(argv[1]); full = (1u << n) - 1;
    E = malloc(sizeof(uint32_t) * ((size_t)1 << n) * n);
    char line[1024]; long classes = 0, pairs = 0, nohp = 0, cond = 0, mis_cond_hp = 0, mis_nohp_nocond = 0;
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] != '0' && line[0] != '1') continue;
        classes++;
        for (int i = 0; i < n; i++) out[i] = 0;
        int k = 0; for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (line[k] == '1') out[i] |= 1u << j; else out[j] |= 1u << i; k++; }
        int s[32], Dl = -1, dl = 99; for (int v = 0; v < n; v++) { s[v] = __builtin_popcount(out[v]); if (s[v] > Dl) Dl = s[v]; if (s[v] < dl) dl = s[v]; }
        endpoint_relation();
        for (int x = 0; x < n; x++) if (s[x] == Dl) for (int y = 0; y < n; y++) if (s[y] == dl && y != x) {
            pairs++;
            int hp = Pi[x] >> y & 1; if (!hp) nohp++;
            uint32_t c[32]; int kc; int ci = 0, cii = 0;
            kc = comps_of(full & ~(1u << x), c);
            if (kc >= 2 && __builtin_popcount(c[kc - 1]) == 1) { int z = __builtin_ctz(c[kc - 1]); if ((out[z] >> x & 1) && c[kc - 2] == (1u << y) && (out[x] >> y & 1)) ci = 1; }
            kc = comps_of(full & ~(1u << y), c);
            if (kc >= 2 && __builtin_popcount(c[0]) == 1) { int z = __builtin_ctz(c[0]); if ((out[y] >> z & 1) && c[1] == (1u << x) && (out[x] >> y & 1)) cii = 1; }
            int cnd = ci || cii; cond += cnd;
            if (cnd && hp) { mis_cond_hp++; printf("condition (i)/(ii) holds but HP exists: %s", line); }
            if (!cnd && !hp) { mis_nohp_nocond++; printf("no HP x..y, neither (i) nor (ii) [x=%d y=%d, scores max %d min %d, regular %d]: %s", x, y, Dl, dl, Dl == dl, line); }
        }
    }
    printf("n=%d classes %ld: (max,min) pairs %ld, without HP x..y %ld; Thomassen (i)/(ii) holds %ld; (i)/(ii) but HP exists %ld; no HP and neither (i) nor (ii) %ld\n", n, classes, pairs, nohp, cond, mis_cond_hp, mis_nohp_nocond);
    return 0;
}
