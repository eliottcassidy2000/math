/* lattice_debt.c -- mac-mini-2026-10-27-slab.  The catalytic limit of merge_debt.c: with unlimited neutral triples the
 * leaves are effectively all six Eisenstein units (-1 = omega + omega^2 costs nothing), and g(x) = min over binary merge
 * trees with root value x of sum over merges of beta(y,z)_+, beta(y,z) = Q(y+z) - Q(y) - Q(z) = 2<y,z>.  g is the
 * largest function with g(unit) = g(0) = 0 and g(y+z) <= g(y) + g(z) + beta(y,z)_+.  Computed by repeated relaxation
 * over all pairs in the window Q(x) <= W (intermediate values outside the window are not used; stability in W is checked
 * by the caller).  Prints g(n) for n = 1..nmax and, optionally, the whole table.  Usage: ./lattice_debt W nmax [dump] */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static int W, A;   /* coordinates a, b in [-A, A] */
#define ID(a,b) (((a) + A) * (2 * A + 1) + ((b) + A))
static int Qf(int a, int b) { return a * a - a * b + b * b; }
int main(int argc, char **argv) {
    W = atoi(argv[1]); int nmax = atoi(argv[2]); int dump = argc > 3;
    A = 1; while (Qf(A, 0) <= W || Qf(A, A / 2) <= W) A++; A += 2;
    int side = 2 * A + 1; int *g = malloc(sizeof(int) * side * side);
    int *pa = malloc(sizeof(int) * side * side), *pb = malloc(sizeof(int) * side * side), np = 0;
    for (int a = -A; a <= A; a++) for (int b = -A; b <= A; b++) { g[ID(a,b)] = 1 << 28; if (Qf(a,b) <= W) { pa[np] = a; pb[np] = b; np++; } }
    g[ID(0,0)] = 0; g[ID(1,0)] = g[ID(-1,0)] = g[ID(0,1)] = g[ID(0,-1)] = g[ID(1,1)] = g[ID(-1,-1)] = 0;
    for (int round = 0; ; round++) {
        int changed = 0;
        for (int i = 0; i < np; i++) { int ya = pa[i], yb = pb[i], gy = g[ID(ya,yb)]; if (gy >= 1 << 28) continue;
            int qy = Qf(ya, yb);
            for (int j = i; j < np; j++) { int za = pa[j], zb = pb[j], gz = g[ID(za,zb)]; if (gz >= 1 << 28) continue;
                int xa = ya + za, xb = yb + zb, qx = Qf(xa, xb); if (qx > W) continue;
                int beta = qx - qy - Qf(za, zb); int c = gy + gz + (beta > 0 ? beta : 0);
                if (c < g[ID(xa,xb)]) { g[ID(xa,xb)] = c; changed = 1; } } }
        if (!changed) { fprintf(stderr, "converged after %d rounds, %d points\n", round + 1, np); break; }
    }
    for (int n = 1; n <= nmax; n++) printf("g(%d) = %d   (n^2 = %d)\n", n, g[ID(n,0)], n * n);
    if (dump) for (int i = 0; i < np; i++) printf("P %d %d %d %d\n", pa[i], pb[i], Qf(pa[i], pb[i]), g[ID(pa[i],pb[i])]);
    return 0;
}
