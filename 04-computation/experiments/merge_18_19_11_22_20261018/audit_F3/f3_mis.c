/* f3_mis.c -- audit F3 (independent referee code): k-fixed tabu search for an independent set of size k in
   G_q = Cay(F_(q^2), mu_(q+1)) (q odd prime). Minimises the number of edges inside a k-subset by swaps
   (v in S out, u not in S in), tabu on reversing a swap, aspiration on new best.
   Usage: ./f3_mis q k maxiter seed  -> prints best (edges inside) and, when 0 is reached, verifies independence. */
#include <stdio.h>
#include <stdlib.h>

static unsigned long long rs = 0x2545F4914F6CDD1DULL;
static inline unsigned long long rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }

int main(int argc, char **argv) {
    int q = atoi(argv[1]), K = atoi(argv[2]); long maxit = atol(argv[3]); rs ^= (unsigned long long)atol(argv[4]) * 0x9E3779B97F4A7C15ULL;
    for (int i = 0; i < 20; i++) rnd();
    int r; for (r = 2; ; r++) { long p = 1; for (int e = 0; e < (q - 1) / 2; e++) p = p * r % q; if (p == q - 1) break; }
    int deg = 0, mua[256], mub[256];
    for (int a = 0; a < q; a++) for (int b = 0; b < q; b++) if ((((long)a * a - (long)r * b * b) % q + q) % q == 1) { mua[deg] = a; mub[deg] = b; deg++; }
    if (deg != q + 1) { fprintf(stderr, "mu\n"); return 1; }
    int n = q * q; int *adj = malloc(sizeof(int) * (size_t)n * deg);
    for (int a = 0; a < q; a++) for (int b = 0; b < q; b++) for (int i = 0; i < deg; i++) adj[(a * q + b) * deg + i] = ((a + mua[i]) % q) * q + (b + mub[i]) % q;
    char *in = calloc(n, 1); int *din = calloc(n, sizeof(int)); long *tabu = calloc(n, sizeof(long));
    int *S = malloc(sizeof(int) * K), *pos = malloc(sizeof(int) * n);
    for (int v = 0; v < n; v++) pos[v] = -1;
    /* random initial K-set */
    int cnt = 0; while (cnt < K) { int v = (int)(rnd() % n); if (!in[v]) { in[v] = 1; pos[v] = cnt; S[cnt++] = v; } }
    for (int v = 0; v < n; v++) if (in[v]) for (int i = 0; i < deg; i++) din[adj[v * deg + i]]++;
    long f = 0; for (int i = 0; i < K; i++) f += din[S[i]]; f /= 2;
    long best = f, bestit = 0, it;
    int *cand = malloc(sizeof(int) * n), *vc = malloc(sizeof(int) * n);
    char *isMu = calloc((size_t)q * q, 1); for (int i = 0; i < deg; i++) isMu[mua[i] * q + mub[i]] = 1;
    #define ADJ(u, v) (isMu[((((u) / q) - ((v) / q) + q) % q) * q + ((((u) % q) - ((v) % q) + q) % q)])
    for (it = 0; it < maxit && f > 0; it++) {
        int mind = 1 << 30, nc = 0, maxd = -1, nv = 0;
        for (int u = 0; u < n; u++) if (!in[u]) { int d = din[u];
            if (tabu[u] > it && !(f + d - 1 < best)) continue;
            if (d < mind) { mind = d; nc = 0; } if (d == mind) cand[nc++] = u; }
        for (int i = 0; i < K; i++) { int v = S[i]; if (tabu[v] > it) continue; int d = din[v];
            if (d > maxd) { maxd = d; nv = 0; } if (d == maxd) vc[nv++] = v; }
        int bu = -1, bv = -1; long bd = 1L << 40; int ties = 0;
        if (nc > 0 && nv > 0) {
            int su = nc < 48 ? nc : 48;
            for (int ci = 0; ci < su; ci++) { int u = (nc <= 48) ? cand[ci] : cand[rnd() % nc];
                for (int vi = 0; vi < nv; vi++) { int v = vc[vi];
                    long d = (long)din[u] - ADJ(u, v) - din[v];
                    if (d < bd) { bd = d; bu = u; bv = v; ties = 1; } else if (d == bd && rnd() % (++ties) == 0) { bu = u; bv = v; } } } }
        if (bu < 0) { bv = S[rnd() % K]; do { bu = (int)(rnd() % n); } while (in[bu]); bd = (long)din[bu] - ADJ(bu, bv) - din[bv]; }
        in[bv] = 0; for (int t = 0; t < deg; t++) din[adj[bv * deg + t]]--;
        in[bu] = 1; for (int t = 0; t < deg; t++) din[adj[bu * deg + t]]++;
        int p = pos[bv]; S[p] = bu; pos[bu] = p; pos[bv] = -1;
        f += bd;
        tabu[bv] = it + 7 + (long)(rnd() % 10);
        tabu[bu] = it + 3 + (long)(rnd() % 5);
        if (f < best) { best = f; bestit = it; }
    }
    long chk = 0; for (int i = 0; i < K; i++) for (int j = i + 1; j < K; j++) chk += ADJ(S[i], S[j]);   /* direct pair check */
    printf("q=%d k=%d maxiter=%ld seed=%s: best edges inside=%ld (at it %ld), final verified=%ld, iterations=%ld\n", q, K, maxit, argv[4], best, bestit, chk, it);
    return 0;
}
