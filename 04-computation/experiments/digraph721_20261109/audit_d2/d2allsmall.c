/* d2allsmall.c - audit D2 (THM-5056): the final statement checked directly for n <= 7.
 * Enumerates ALL labelled graphs (connected or not) on n vertices, computes d(G) and whether some automorphism
 * reverses an odd number of edges (complete backtracking enumeration of Aut), and compares
 *   {d(G) : G odd on n vertices}   with   {odd k <= 2^n - 1} minus {1, 5, 37}.
 * Independent of the census and of the construction. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#define NMAX 7
static int n; static uint32_t adj[NMAX];
static int img[NMAX], used[NMAX], dg[NMAX], anyodd;
static int edge(int a, int b) { return (adj[a] >> b) & 1; }
static void bt(int v) {
    if (anyodd) return;
    if (v == n) {
        int rev = 0;
        for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) if (edge(a, b) && img[a] > img[b]) rev++;
        if (rev & 1) anyodd = 1;
        return;
    }
    for (int c = 0; c < n && !anyodd; c++) {
        if (used[c] || dg[c] != dg[v]) continue;
        int ok = 1;
        for (int j = 0; j < v; j++) if (edge(v, j) != edge(c, img[j])) { ok = 0; break; }
        if (!ok) continue;
        img[v] = c; used[c] = 1; bt(v + 1); used[c] = 0;
    }
}
int main(void) {
    int allok = 1;
    for (n = 1; n <= NMAX; n++) {
        int E = n * (n - 1) / 2, ea[32], eb[32], t = 0;
        for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) { ea[t] = a; eb[t] = b; t++; }
        static char oddval[1 << NMAX], anyval[1 << NMAX];
        memset(oddval, 0, sizeof oddval); memset(anyval, 0, sizeof anyval);
        for (uint64_t mask = 0; mask < (1ull << E); mask++) {
            memset(adj, 0, sizeof adj);
            for (int e = 0; e < E; e++) if ((mask >> e) & 1) { adj[ea[e]] |= 1u << eb[e]; adj[eb[e]] |= 1u << ea[e]; }
            uint32_t full = (1u << n) - 1; int d = 0;
            for (uint32_t S = 1; S <= full; S++) {
                uint32_t dom = 0;
                for (int v = 0; v < n; v++) if ((S >> v) & 1) dom |= adj[v] | (1u << v);
                if (dom == full) d++;
            }
            anyval[d] = 1;
            for (int v = 0; v < n; v++) { dg[v] = __builtin_popcount(adj[v]); used[v] = 0; }
            anyodd = 0; bt(0);
            if (anyodd) oddval[d] = 1;
        }
        int ok = 1, cnt = 0;
        printf("n=%d odd-graph d-values:", n);
        for (int k = 0; k < (1 << n); k++) {
            int expect = (k % 2 == 1) && k != 1 && k != 5 && k != 37;
            if (oddval[k]) { printf(" %d", k); cnt++; }
            if (oddval[k] != expect) ok = 0;
            if (anyval[k] != (k % 2 == 1)) { printf(" [all-graph value set wrong at %d]", k); ok = 0; }
        }
        printf("\n  count=%d  matches {odd k <= 2^n-1} minus {1,5,37}: %s\n", cnt, ok ? "YES" : "NO");
        if (!ok) allok = 0;
        fflush(stdout);
    }
    printf(allok ? "ALL n <= 7 MATCH\n" : "MISMATCH\n");
    return !allok;
}
