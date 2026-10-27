/* merge_debt.c -- mac-mini-2026-10-27-slab.  Reading A of the owner's "minimum merge cost for charges (11+k,k,k)":
 * a pile of p unit charges 1, q charges omega, r charges omega^2 (Eisenstein units) is merged into one pile by a binary
 * tree; a merge of piles with count vectors u, v costs beta(x,y)_+ where beta(x,y) = Q(x+y) - Q(x) - Q(y),
 * Q(a + b omega) = a^2 - ab + b^2.  In counts, beta = 3 u.v - |u||v| (|u| = number of unit charges).  The minimum total
 * cost f(u) obeys f(e_i) = 0, f(u) = min over u = u1 + u2 of beta(u1,u2)_+ + f(u1) + f(u2) (a subtree's cost depends only
 * on its leaf multiset).  Prints f(n0 + k, k, k) for k = 0..K.  Usage: ./merge_debt n0 K [tree]  (with a third argument: print an optimal tree for (n0+K, K, K)) */
#include <stdio.h>
#include <stdlib.h>
int main(int argc, char **argv) {
    int n0 = atoi(argv[1]), K = atoi(argv[2]); int P = n0 + K, Qm = K, R = K;
    #define IDX(a,b,c) (((a) * (Qm + 1) + (b)) * (R + 1) + (c))
    int *f = malloc(sizeof(int) * (P + 1) * (Qm + 1) * (R + 1)); int *arg = malloc(sizeof(int) * (P + 1) * (Qm + 1) * (R + 1));
    for (int s = 1; s <= P + Qm + R; s++)
    for (int a = 0; a <= P; a++) for (int b = 0; b <= Qm; b++) { int c = s - a - b; if (c < 0 || c > R) continue;
        if (s == 1) { f[IDX(a,b,c)] = 0; arg[IDX(a,b,c)] = -1; continue; }
        int best = 1 << 30, bi = -1;
        for (int x = 0; x <= a; x++) for (int y = 0; y <= b; y++) for (int z = 0; z <= c; z++) {
            int s1 = x + y + z; if (s1 == 0 || s1 == s) continue;
            int x2 = a - x, y2 = b - y, z2 = c - z;
            if (x > x2 || (x == x2 && (y > y2 || (y == y2 && z > z2)))) continue;   /* unordered split */
            int beta = 3 * (x * x2 + y * y2 + z * z2) - s1 * (s - s1);
            int cost = (beta > 0 ? beta : 0) + f[IDX(x,y,z)] + f[IDX(x2,y2,z2)];
            if (cost < best) { best = cost; bi = IDX(x,y,z); }
        }
        f[IDX(a,b,c)] = best; arg[IDX(a,b,c)] = bi;
    }
    if (argc > 3) {   /* print an optimal tree for (n0 + K, K, K) as merge lines "parent = left + right" */
        int stack[4096][3], sp = 0; stack[sp][0] = n0 + K; stack[sp][1] = K; stack[sp][2] = K; sp++;
        while (sp) { sp--; int a = stack[sp][0], b = stack[sp][1], c = stack[sp][2]; int i = arg[IDX(a,b,c)]; if (i < 0) continue;
            int z = i % (R + 1), y = (i / (R + 1)) % (Qm + 1), x = i / ((R + 1) * (Qm + 1));
            printf("MERGE %d %d %d = %d %d %d + %d %d %d\n", a, b, c, x, y, z, a - x, b - y, c - z);
            stack[sp][0] = x; stack[sp][1] = y; stack[sp][2] = z; sp++; stack[sp][0] = a - x; stack[sp][1] = b - y; stack[sp][2] = c - z; sp++; }
        return 0; }
    for (int k = 0; k <= K; k++) printf("k=%d f(%d,%d,%d)=%d  (telescoped bound (Q(net) - |leaves|)_+ = %d)\n", k, n0 + k, k, k, f[IDX(n0 + k, k, k)], (n0 * n0 - n0 - 3 * k) > 0 ? n0 * n0 - n0 - 3 * k : 0);
    return 0;
}
