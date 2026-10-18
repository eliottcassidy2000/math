/* f3_tabu.c -- audit F3 (independent referee code, written from scratch): TabuCol (Galinier-Hao 1999 parameters,
   incremental conflict lists) for G_q = Cay(F_(q^2), mu_(q+1)), q an odd prime; mu built from a primitive element.
   Usage: ./f3_tabu q k maxiter seed [outfile]
   Prints the best conflict count and the iteration at which 0 was reached; writes a colouring when one is found
   (header in the session's format: "# q=.. r=.. (t^2=r); vertex a*q+b for a+bt; colours 0..k-1"). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static unsigned long long rs = 0x9E3779B97F4A7C15ULL;
static inline unsigned long long rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }

static int q, r;
static void mul(int a, int b, int c, int d, int *x, int *y) { /* (a+bt)(c+dt), t^2 = r */
    *x = (int)(((long)a * c + (long)r * b % q * d) % q); *y = (int)(((long)a * d + (long)b * c) % q); }

int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage\n"); return 1; }
    q = atoi(argv[1]); int k = atoi(argv[2]); long maxit = atol(argv[3]); rs ^= (unsigned long long)atol(argv[4]) * 0xBF58476D1CE4E5B9ULL;
    for (int i = 0; i < 20; i++) rnd();
    for (r = 2; ; r++) { long p = 1; for (int e = 0; e < (q - 1) / 2; e++) p = p * r % q; if (p == q - 1) break; }
    /* primitive element g of F_(q^2)^*: order q^2 - 1 */
    long n1 = (long)q * q - 1; int pf[32], npf = 0; { long m = n1; for (long p = 2; p * p <= m; p++) if (m % p == 0) { pf[npf++] = (int)p; while (m % p == 0) m /= p; } if (m > 1) pf[npf++] = (int)m; }
    int ga = -1, gb = -1;
    for (int a = 0; a < q && ga < 0; a++) for (int b = 1; b < q && ga < 0; b++) {
        int ok = 1;
        for (int i = 0; i < npf && ok; i++) { long e = n1 / pf[i]; int xa = 1, xb = 0, ba = a, bb = b;
            while (e) { if (e & 1) mul(xa, xb, ba, bb, &xa, &xb); mul(ba, bb, ba, bb, &ba, &bb); e >>= 1; }
            if (xa == 1 && xb == 0) ok = 0; }
        if (ok) { ga = a; gb = b; } }
    /* h = g^(q-1) generates mu_(q+1) */
    int ha = 1, hb = 0; { long e = q - 1; int ba = ga, bb = gb; while (e) { if (e & 1) mul(ha, hb, ba, bb, &ha, &hb); mul(ba, bb, ba, bb, &ba, &bb); e >>= 1; } }
    int deg = q + 1, *mua = malloc(sizeof(int) * deg), *mub = malloc(sizeof(int) * deg);
    { int xa = 1, xb = 0; for (int i = 0; i < deg; i++) { mua[i] = xa; mub[i] = xb; mul(xa, xb, ha, hb, &xa, &xb); } if (!(xa == 1 && xb == 0)) { fprintf(stderr, "mu order\n"); return 1; } }
    int n = q * q;
    int *adj = malloc(sizeof(int) * (size_t)n * deg);
    for (int a = 0; a < q; a++) for (int b = 0; b < q; b++) for (int i = 0; i < deg; i++)
        adj[(a * q + b) * deg + i] = ((a + mua[i]) % q) * q + (b + mub[i]) % q;
    int *col = malloc(sizeof(int) * n), *gam = calloc((size_t)n * k, sizeof(int));
    long *tabu = calloc((size_t)n * k, sizeof(long));
    int *clist = malloc(sizeof(int) * n), *cpos = malloc(sizeof(int) * n), nc = 0;
    for (int v = 0; v < n; v++) col[v] = (int)(rnd() % k);
    for (int v = 0; v < n; v++) for (int i = 0; i < deg; i++) gam[v * k + col[adj[v * deg + i]]]++;
    long conf = 0;
    for (int v = 0; v < n; v++) { conf += gam[v * k + col[v]]; cpos[v] = -1; if (gam[v * k + col[v]] > 0) { cpos[v] = nc; clist[nc++] = v; } }
    conf /= 2;
    long best = conf, bestit = 0, it;
    for (it = 0; it < maxit && conf > 0; it++) {
        int bv = -1, bc = -1, bd = 1 << 30, cnt = 0;
        for (int ii = 0; ii < nc; ii++) { int v = clist[ii], cv = col[v], gv = gam[v * k + cv];
            for (int c = 0; c < k; c++) { if (c == cv) continue; int d = gam[v * k + c] - gv;
                if (tabu[v * k + c] > it && !(conf + d < best)) continue;
                if (d < bd) { bd = d; bv = v; bc = c; cnt = 1; } else if (d == bd && rnd() % (++cnt) == 0) { bv = v; bc = c; } } }
        if (bv < 0) { int ii = (int)(rnd() % nc); bv = clist[ii]; bc = (int)(rnd() % k); if (bc == col[bv]) bc = (bc + 1) % k; bd = gam[bv * k + bc] - gam[bv * k + col[bv]]; }
        int old = col[bv]; col[bv] = bc; conf += bd;
        for (int i = 0; i < deg; i++) { int u = adj[bv * deg + i];
            gam[u * k + old]--; gam[u * k + bc]++;
            int cu = gam[u * k + col[u]] > 0;
            if (cu && cpos[u] < 0) { cpos[u] = nc; clist[nc++] = u; }
            else if (!cu && cpos[u] >= 0) { int last = clist[--nc]; clist[cpos[u]] = last; cpos[last] = cpos[u]; cpos[u] = -1; } }
        { int cu = gam[bv * k + bc] > 0;
          if (cu && cpos[bv] < 0) { cpos[bv] = nc; clist[nc++] = bv; }
          else if (!cu && cpos[bv] >= 0) { int last = clist[--nc]; clist[cpos[bv]] = last; cpos[last] = cpos[bv]; cpos[bv] = -1; } }
        tabu[bv * k + old] = it + (long)(rnd() % 10) + (long)(0.6 * nc);
        if (conf < best) { best = conf; bestit = it; }
    }
    long check = 0; for (int v = 0; v < n; v++) for (int i = 0; i < deg; i++) if (col[v] == col[adj[v * deg + i]]) check++;
    printf("q=%d k=%d maxiter=%ld seed=%s: best=%ld (at it %ld), final verified conflicts=%ld, iterations used=%ld\n", q, k, maxit, argv[4], best, bestit, check / 2, it);
    if (check == 0 && argc > 5) { FILE *fo = fopen(argv[5], "w");
        fprintf(fo, "# q=%d r=%d (t^2=r); vertex a*q+b for a+bt; colours 0..%d; proper colouring of Cay(F_q^2, mu_(q+1)) [audit F3 f3_tabu.c]\n", q, r, k - 1);
        for (int v = 0; v < n; v++) fprintf(fo, "%d%c", col[v], v == n - 1 ? '\n' : ' '); fclose(fo); }
    return 0;
}
