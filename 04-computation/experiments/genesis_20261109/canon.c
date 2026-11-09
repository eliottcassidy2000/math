/* canon.c - canonical forms of tournaments of order n <= 16 (genesis lane, mac-mini-2026-11-09).
 * Individualization-refinement with ordered colour refinement (the scheme of kind-pasteur S16/S17's
 * classgen/canon11.c), but with 128-bit codes so that n = 12..16 fit (C(16,2) = 120 bits).
 * Input lines:  n code   (code = upper triangle, row-major, '1' iff i -> j)
 * Output line:  n canonical_code(as 0/1 string) |Aut| code_in
 * The canonical code is the lexicographically smallest... (we minimise the 128-bit integer whose bit b is pair b).
 * Compile: cc -O2 -o canon canon.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned __int128 u128;
#define MAXN 16
static int nn;
static uint32_t A[MAXN];
static u128 best; static int have; static long leaves_equal;

static u128 code_of(const int *pos) {
    int inv[MAXN];
    for (int v = 0; v < nn; v++) inv[pos[v]] = v;
    u128 c = 0; int b = 0;
    for (int i = 0; i < nn; i++) for (int j = i + 1; j < nn; j++) {
        if (A[inv[i]] >> inv[j] & 1) c |= ((u128)1) << b;
        b++;
    }
    return c;
}
static int ncells(const int *cell) { int m = -1; for (int v = 0; v < nn; v++) if (cell[v] > m) m = cell[v]; return m + 1; }
static void refine(int *cell) {
    for (;;) {
        int k = ncells(cell);
        int sig[MAXN][MAXN + 1];
        for (int v = 0; v < nn; v++) {
            sig[v][0] = cell[v];
            for (int c = 0; c < k; c++) sig[v][c + 1] = 0;
            for (int u = 0; u < nn; u++) if (A[v] >> u & 1) sig[v][cell[u] + 1]++;
        }
        int ord[MAXN];
        for (int v = 0; v < nn; v++) ord[v] = v;
        for (int i = 1; i < nn; i++) {
            int x = ord[i], j = i - 1;
            while (j >= 0 && memcmp(sig[ord[j]], sig[x], (k + 1) * sizeof(int)) > 0) { ord[j + 1] = ord[j]; j--; }
            ord[j + 1] = x;
        }
        int nc[MAXN]; int c = 0; nc[ord[0]] = 0;
        for (int i = 1; i < nn; i++) { if (memcmp(sig[ord[i]], sig[ord[i - 1]], (k + 1) * sizeof(int))) c++; nc[ord[i]] = c; }
        int changed = (c + 1 != k);
        memcpy(cell, nc, sizeof(int) * nn);
        if (!changed) return;
    }
}
static void search(int *cell) {
    refine(cell);
    int k = ncells(cell);
    if (k == nn) {
        u128 c = code_of(cell);
        if (!have || c < best) { best = c; have = 1; leaves_equal = 1; }
        else if (c == best) leaves_equal++;
        return;
    }
    int cnt[MAXN] = {0};
    for (int v = 0; v < nn; v++) cnt[cell[v]]++;
    int target = -1;
    for (int c = 0; c < k; c++) if (cnt[c] > 1) { target = c; break; }
    for (int v = 0; v < nn; v++) if (cell[v] == target) {
        int c2[MAXN];
        for (int u = 0; u < nn; u++) {
            c2[u] = cell[u] > target ? cell[u] + 1 : cell[u];
            if (cell[u] == target && u != v) c2[u] = target + 1;
        }
        search(c2);
    }
}
int main(void) {
    char code[256];
    int n;
    while (scanf("%d %255s", &n, code) == 2) {
        nn = n; memset(A, 0, sizeof A);
        int p = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
            if (code[p] == '1') A[i] |= 1u << j; else A[j] |= 1u << i;
            p++;
        }
        have = 0;
        int cell[MAXN];
        for (int v = 0; v < n; v++) cell[v] = 0;
        search(cell);
        char out[256]; int b = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { out[b] = (best >> b & 1) ? '1' : '0'; b++; }
        out[b] = 0;
        printf("%d %s %ld %s\n", n, out, leaves_equal, code);
    }
    return 0;
}
