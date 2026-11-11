/* AUDIT program (independent): the pentagon cube sum
 *     G_35(R) = sum over the 2^n one-vertex extensions R + v_X of (-1)^((c_3 + c_5)(R + v_X)).
 * Two methods:
 *   full (n <= 12): c_3 and c_5 of every extension counted directly (all 3-subsets, all 5-subsets, all
 *        cyclic orders);
 *   walk (n <= 31): (c_3 + c_5)(R + v_X) = (c_3 + c_5)(R) + #{a in X, b not in X : a -> b}
 *        + #{a in X, b not in X} (A^3)_{ab}, because in a tournament every walk of length 3 between two
 *        different vertices is a path; Gray code over all X. (c_3 + c_5)(R) is counted directly.
 * usage: aud_cube < input      (input: n, then n rows of 0/1)
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define MAXN 33
static int n;
static int A[MAXN][MAXN];

static long long count35(int N, int a[MAXN][MAXN], long long *c3out, long long *c5out) {
    long long c3 = 0, c5 = 0;
    int i, j, k, l, m2;
    for (i = 0; i < N; i++) for (j = i + 1; j < N; j++) for (k = j + 1; k < N; k++) {
        if (a[i][j] && a[j][k] && a[k][i]) c3++;
        if (a[i][k] && a[k][j] && a[j][i]) c3++;
    }
    int v[5], p[4];
    static const int perms[24][4] = {
        {1,2,3,4},{1,2,4,3},{1,3,2,4},{1,3,4,2},{1,4,2,3},{1,4,3,2},
        {2,1,3,4},{2,1,4,3},{2,3,1,4},{2,3,4,1},{2,4,1,3},{2,4,3,1},
        {3,1,2,4},{3,1,4,2},{3,2,1,4},{3,2,4,1},{3,4,1,2},{3,4,2,1},
        {4,1,2,3},{4,1,3,2},{4,2,1,3},{4,2,3,1},{4,3,1,2},{4,3,2,1}};
    for (i = 0; i < N; i++) for (j = i + 1; j < N; j++) for (k = j + 1; k < N; k++)
        for (l = k + 1; l < N; l++) for (m2 = l + 1; m2 < N; m2++) {
            v[0] = i; v[1] = j; v[2] = k; v[3] = l; v[4] = m2;
            int t;
            for (t = 0; t < 24; t++) {
                p[0] = v[perms[t][0]]; p[1] = v[perms[t][1]]; p[2] = v[perms[t][2]]; p[3] = v[perms[t][3]];
                if (a[v[0]][p[0]] && a[p[0]][p[1]] && a[p[1]][p[2]] && a[p[2]][p[3]] && a[p[3]][v[0]]) c5++;
            }
        }
    *c3out = c3; *c5out = c5;
    return c3 + c5;
}

int main(void) {
    int i, j, k;
    char row[64];
    if (scanf("%d", &n) != 1 || n < 1 || n > 31) { fprintf(stderr, "bad n\n"); return 1; }
    for (i = 0; i < n; i++) {
        if (scanf("%63s", row) != 1) return 1;
        for (j = 0; j < n; j++) A[i][j] = row[j] == '1';
    }
    long long c3, c5;
    count35(n, A, &c3, &c5);
    int base = (int)((c3 + c5) & 1);
    printf("order %d: c_3 = %lld, c_5 = %lld\n", n, c3, c5);
    /* B = A + A^3 mod 2 */
    static long long A2[MAXN][MAXN], A3[MAXN][MAXN];
    for (i = 0; i < n; i++) for (j = 0; j < n; j++) { long long s = 0; for (k = 0; k < n; k++) s += A[i][k] * A[k][j]; A2[i][j] = s; }
    for (i = 0; i < n; i++) for (j = 0; j < n; j++) { long long s = 0; for (k = 0; k < n; k++) s += A2[i][k] * A[k][j]; A3[i][j] = s; }
    uint32_t Brow[MAXN], Bcol[MAXN];
    for (i = 0; i < n; i++) { Brow[i] = 0; Bcol[i] = 0; }
    for (i = 0; i < n; i++) for (j = 0; j < n; j++) if (i != j && ((A[i][j] + A3[i][j]) & 1)) { Brow[i] |= 1u << j; Bcol[j] |= 1u << i; }
    /* walk method, Gray code */
    uint64_t total = 1ULL << n, g;
    uint32_t x = 0; int par = base; /* X empty: no cycle through v */
    long long even = 0, odd = 0;
    if (par) odd++; else even++;
    for (g = 1; g < total; g++) {
        int b = __builtin_ctzll(g);           /* bit that flips */
        uint32_t bit = 1u << b, rest = x & ~bit;
        uint32_t full = (n == 32) ? 0xffffffffu : ((1u << n) - 1);
        /* pairs (b, y) with y outside X u {b}, and pairs (a, b) with a in X \ {b}, toggle */
        int d = (__builtin_popcount(Brow[b] & ~rest & ~bit & full) + __builtin_popcount(Bcol[b] & rest)) & 1;
        par ^= d;
        x ^= bit;
        if (par) odd++; else even++;
    }
    printf("walk method: extensions with c_3 + c_5 even: %lld, odd: %lld, G_35 = %lld\n", even, odd, even - odd);
    if (n <= 12) {
        long long e2 = 0, o2 = 0;
        static int E[MAXN][MAXN];
        uint32_t X;
        for (X = 0; X < (1u << n); X++) {
            for (i = 0; i < n; i++) for (j = 0; j < n; j++) E[i][j] = A[i][j];
            for (i = 0; i < n; i++) { E[n][i] = (X >> i) & 1; E[i][n] = 1 - E[n][i]; }
            E[n][n] = 0;
            long long d3, d5;
            count35(n + 1, E, &d3, &d5);
            if ((d3 + d5) & 1) o2++; else e2++;
        }
        printf("full method: extensions with c_3 + c_5 even: %lld, odd: %lld, G_35 = %lld (%s)\n", e2, o2, e2 - o2,
               (e2 == even && o2 == odd) ? "agrees with the walk method" : "DISAGREES");
    }
    return 0;
}
