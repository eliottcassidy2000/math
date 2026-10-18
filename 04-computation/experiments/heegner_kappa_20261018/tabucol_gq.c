/* tabucol_gq.c -- tabu-search colouring of G_q = Cay(F_(q^2), mu_(q+1)), q an odd prime (mac-mini-2026-10-18).
   F_(q^2) = F_q[t]/(t^2 - r), r the least non-residue; vertex a + b t has index a*q + b.
   Usage: ./tabucol_gq q k maxiter seed  -> prints the best conflict count; on success writes the colouring and verifies it. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static unsigned long long rs = 88172645463325252ULL;
static inline unsigned long long xr(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
int main(int argc, char **argv) {
    int q = atoi(argv[1]), k = atoi(argv[2]); long maxit = atol(argv[3]); rs ^= (unsigned long long)atol(argv[4]) * 2654435761ULL;
    int r = 2; while (1) { long p = 1; for (int e = 0; e < (q - 1) / 2; e++) p = p * r % q; if (p == q - 1) break; r++; }
    int nmu = 0, mua[4096], mub[4096];
    for (int a = 0; a < q; a++) for (int b = 0; b < q; b++) if (((long)a * a - (long)r * b % q * b % q + 2L * q * q) % q == 1) { mua[nmu] = a; mub[nmu] = b; nmu++; }
    int n = q * q, deg = nmu;
    int *adj = malloc(sizeof(int) * n * deg);
    for (int a = 0; a < q; a++) for (int b = 0; b < q; b++) for (int i = 0; i < deg; i++) adj[(a * q + b) * deg + i] = ((a + mua[i]) % q) * q + (b + mub[i]) % q;
    int *col = malloc(sizeof(int) * n), *gam = calloc((size_t)n * k, sizeof(int)); long *tabu = calloc((size_t)n * k, sizeof(long));
    for (int v = 0; v < n; v++) col[v] = xr() % k;
    for (int v = 0; v < n; v++) for (int i = 0; i < deg; i++) gam[v * k + col[adj[v * deg + i]]]++;
    long conf = 0; for (int v = 0; v < n; v++) conf += gam[v * k + col[v]]; conf /= 2;
    long best = conf;
    for (long it = 0; it < maxit && conf > 0; it++) {
        int bv = -1, bc = -1, bd = 1 << 30, cnt = 0;
        for (int v = 0; v < n; v++) { int cv = col[v]; if (gam[v * k + cv] == 0) continue;
            for (int c = 0; c < k; c++) { if (c == cv) continue; int d = gam[v * k + c] - gam[v * k + cv];
                if (tabu[v * k + c] > it && !(conf + d < best)) continue;
                if (d < bd) { bd = d; bv = v; bc = c; cnt = 1; } else if (d == bd && (xr() % (++cnt)) == 0) { bv = v; bc = c; } } }
        if (bv < 0) { bv = xr() % n; bc = xr() % k; if (bc == col[bv]) bc = (bc + 1) % k; bd = gam[bv * k + bc] - gam[bv * k + col[bv]]; }
        int old = col[bv]; col[bv] = bc; conf += bd;
        for (int i = 0; i < deg; i++) { int u = adj[bv * deg + i]; gam[u * k + old]--; gam[u * k + bc]++; }
        tabu[bv * k + old] = it + 10 + (long)(0.6 * conf) + xr() % 10;
        if (conf < best) best = conf;
    }
    long check = 0; for (int v = 0; v < n; v++) for (int i = 0; i < deg; i++) if (col[v] == col[adj[v * deg + i]]) check++;
    printf("q=%d k=%d |mu|=%d best_conflicts=%ld final_conflicts(verified)=%ld\n", q, k, deg, best, check / 2);
    if (check == 0) { char fn[64]; sprintf(fn, "kappa_%d_colouring_k%d.txt", q, k); FILE *fo = fopen(fn, "w");
        fprintf(fo, "# q=%d r=%d (t^2=r); vertex a*q+b for a+bt; colours 0..%d; proper colouring of Cay(F_q^2, mu_(q+1))\n", q, r, k - 1);
        for (int v = 0; v < n; v++) fprintf(fo, "%d%c", col[v], v == n - 1 ? '\n' : ' '); fclose(fo); }
    return 0;
}
