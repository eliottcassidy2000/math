/* P111 test.  For each odd count h (h != 7, 21) with p = 2h-1 prime (then p = 1 mod 4 automatically), look for a
   tournament T with H(T) = h and ONE added vertex x (arbitrary out-set) with H(T+x) = 2h-1 (THM-5272's doubling D).
   Bases T: one-parent trees (each new vertex beats only its parent) with up to NB vertices; every out-set of x tried.
   Also record whether a one-parent child (x beats exactly one vertex) suffices.  Usage: p111 NB PMAX */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long ull;
static ull *dp;
static ull Hc(int n, const unsigned *out) { int N = 1 << n; memset(dp, 0, sizeof(ull) * (size_t)N * n);
    for (int v = 0; v < n; v++) dp[(size_t)(1u << v) * n + v] = 1;
    for (int S = 1; S < N; S++) { ull *row = dp + (size_t)S * n; for (int v = 0; v < n; v++) { ull c = row[v]; if (!c) continue;
        unsigned t = out[v] & ~(unsigned)S; while (t) { int u = __builtin_ctz(t); t &= t - 1; dp[(size_t)(S | 1u << u) * n + u] += c; } } }
    ull s = 0; for (int v = 0; v < n; v++) s += dp[(size_t)(N - 1) * n + v]; return s; }
static int isprime(ull x) { if (x < 2) return 0; for (ull d = 2; d * d <= x; d++) if (x % d == 0) return 0; return 1; }
#define HMAX 20000
static signed char status[HMAX];   /* 0 unknown target, 1 found (general), 2 found one-parent, -1 not a target */
static int found_n[HMAX]; static int tries[HMAX];
static int p[32]; static unsigned out[32]; static int NB;
static void build(int n) { for (int i = 0; i < n; i++) out[i] = 0;
    for (int j = 1; j < n; j++) for (int i = 0; i < j; i++) { if (p[j] == i) out[j] |= 1u << i; else out[i] |= 1u << j; } }
static void try_tree(int n) {
    build(n); ull h = Hc(n, out);
    if (h >= HMAX || status[h] != 0 || ++tries[h] > 300) return;
    unsigned base[32]; memcpy(base, out, sizeof(unsigned) * n);
    for (unsigned S = 0; S < (1u << n); S++) {          /* x beats exactly the vertices in S */
        for (int i = 0; i < n; i++) out[i] = base[i] | ((S >> i & 1) ? 0 : (1u << n));
        out[n] = S;
        if (Hc(n + 1, out) == 2 * h - 1) {
            int onep = __builtin_popcount(S) == 1;
            if (onep) { status[h] = 2; found_n[h] = n; return; }
            status[h] = 1; found_n[h] = n; return;
            if (status[h] == 0) { status[h] = 1; found_n[h] = n; }
        }
    }
}
static void rec(int j, int n) { if (j == n) { try_tree(n); return; } for (int q = 0; q < j; q++) { p[j] = q; rec(j + 1, n); } }
int main(int argc, char **argv) {
    NB = atoi(argv[1]); ull PMAX = strtoull(argv[2], 0, 10);
    dp = malloc(sizeof(ull) * ((size_t)1 << (NB + 1)) * (NB + 1));
    for (int h = 0; h < HMAX; h++) status[h] = -1;
    int ntarget = 0;
    for (ull h = 3; 2 * h - 1 <= PMAX && h < HMAX; h += 2) if (h != 7 && h != 21 && isprime(2 * h - 1)) { status[h] = 0; ntarget++; }
    for (int n = 1; n <= NB; n++) rec(1, n);
    int gen = 0, onep = 0, miss = 0;
    printf("targets h (2h-1 prime <= %llu, h != 7, 21): %d\n", PMAX, ntarget);
    printf("not found (h, p):");
    for (ull h = 3; h < HMAX; h += 2) { if (status[h] < 0) continue;
        if (status[h] == 2) onep++; else if (status[h] == 1) gen++; else { miss++; if (miss <= 40) printf(" (%llu,%llu)", h, 2 * h - 1); } }
    printf("\nfound via a one-parent child: %d; only via a general extension: %d; not found with bases <= %d vertices: %d\n", onep, gen, NB, miss);
}
