/* AUDIT G1: numerical tests of the band lemmas (own code).
 * band_s(n): for i<j, j->i iff j-i <= s, else i->j.
 *  (B2) lambda(band_s(n)) = s for 2s+1 <= n <= 24, s <= 6: unit-capacity max-flow (Edmonds-Karp), and full cut
 *       enumeration for n <= 20 as a cross-check.
 *  (B1) prime iff n >= 2s+3 (for n >= 2s+2), by module closure of every pair, and by full subset enumeration for n <= 18;
 *       at n = 2s+2 check that {0, 2s+1} is a module.
 *  (W)  two-step walk lemma: a,b <= 25 coprime, a+b <= n <= 80.
 *  (D)  B2's construction: for each s <= 8 find a coprime bijection sigma:{1..s}->{s+1..2s} by search; check the D_r
 *       are inside band_s(n), pairwise arc-disjoint, spanning and strongly connected for 3s <= n <= 40.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MAXN 81
static int n;
static unsigned long long adjm[MAXN]; /* for n <= 64 */
static int A[MAXN][MAXN];             /* A[i][j]=1 iff i->j */

static void band(int s, int nn) {
    n = nn;
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) A[i][j] = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++) {
            if (j - i <= s) A[j][i] = 1; else A[i][j] = 1;
        }
}
static int maxflow(int src, int dst) {
    static int cap[MAXN][MAXN];
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) cap[i][j] = A[i][j];
    int flow = 0;
    int par[MAXN], q[MAXN];
    for (;;) {
        for (int i = 0; i < n; i++) par[i] = -1;
        par[src] = src;
        int h = 0, t = 0;
        q[t++] = src;
        while (h < t && par[dst] < 0) {
            int v = q[h++];
            for (int w = 0; w < n; w++)
                if (cap[v][w] > 0 && par[w] < 0) { par[w] = v; q[t++] = w; }
        }
        if (par[dst] < 0) return flow;
        for (int w = dst; w != src; w = par[w]) { cap[par[w]][w]--; cap[w][par[w]]++; }
        flow++;
    }
}
static int lambda_flow(void) {
    int best = 1 << 30;
    for (int v = 1; v < n; v++) {
        int f = maxflow(0, v); if (f < best) best = f;
        f = maxflow(v, 0); if (f < best) best = f;
    }
    return best;
}
static int lambda_enum(void) {
    unsigned full = (1u << n) - 1;
    unsigned outm[32];
    for (int i = 0; i < n; i++) { outm[i] = 0; for (int j = 0; j < n; j++) if (A[i][j]) outm[i] |= 1u << j; }
    int best = 1 << 30;
    for (unsigned S = 1; S < full; S++) {
        int c = 0;
        for (int v = 0; v < n; v++) if (S >> v & 1) c += __builtin_popcount(outm[v] & ~S & full);
        if (c < best) best = c;
    }
    return best;
}
static int prime_closure(void) {
    for (int u = 0; u < n; u++)
        for (int w = u + 1; w < n; w++) {
            int in[MAXN] = {0};
            in[u] = in[w] = 1;
            int sz = 2, ch = 1;
            while (ch) {
                ch = 0;
                for (int x = 0; x < n; x++) {
                    if (in[x]) continue;
                    int o = 0, i = 0;
                    for (int m = 0; m < n; m++) if (in[m]) { if (A[x][m]) o++; else i++; }
                    if (o && i) { in[x] = 1; sz++; ch = 1; }
                }
            }
            if (sz < n) return 0;
        }
    return 1;
}
static int prime_enum(void) {
    unsigned full = (1u << n) - 1;
    unsigned outm[32];
    for (int i = 0; i < n; i++) { outm[i] = 0; for (int j = 0; j < n; j++) if (A[i][j]) outm[i] |= 1u << j; }
    for (unsigned M = 1; M < full; M++) {
        if (__builtin_popcount(M) < 2) continue;
        int ok = 1;
        for (int x = 0; x < n && ok; x++) {
            if (M >> x & 1) continue;
            unsigned o = outm[x] & M;
            if (o && o != M) ok = 0;
        }
        if (ok) return 0;
    }
    return 1;
}
static int is_module_pair(int a, int b) {
    for (int x = 0; x < n; x++) if (x != a && x != b && A[x][a] != A[x][b]) return 0;
    return 1;
}
static int strongly_connected_arcs(int nn, int (*arc)(int, int, void *), void *ctx) {
    /* reachability from 0 forward and backward */
    int seen[MAXN], q[MAXN];
    for (int dir = 0; dir < 2; dir++) {
        for (int i = 0; i < nn; i++) seen[i] = 0;
        int h = 0, t = 0;
        seen[0] = 1; q[t++] = 0;
        while (h < t) {
            int v = q[h++];
            for (int w = 0; w < nn; w++) {
                int e = dir ? arc(w, v, ctx) : arc(v, w, ctx);
                if (e && !seen[w]) { seen[w] = 1; q[t++] = w; }
            }
        }
        for (int i = 0; i < nn; i++) if (!seen[i]) return 0;
    }
    return 1;
}
struct walk { int a, b; };
static int walk_arc(int v, int w, void *c) { struct walk *k = c; return (w == v + k->a) || (w == v - k->b); }
static int gcd(int a, int b) { while (b) { int t = a % b; a = b; b = t; } return a; }
static int sig[20], usedv[20], S_;
static int findsig(int r) {
    if (r > S_) return 1;
    for (int v = S_ + 1; v <= 2 * S_; v++)
        if (!usedv[v] && gcd(r, v) == 1) {
            usedv[v] = 1; sig[r] = v;
            if (findsig(r + 1)) return 1;
            usedv[v] = 0;
        }
    return 0;
}
int main(void) {
    int bad = 0;
    for (int s = 1; s <= 6; s++) {
        printf("s=%d:", s);
        for (int nn = 2 * s + 1; nn <= 24; nn++) {
            band(s, nn);
            int lf = lambda_flow();
            int le = (nn <= 20) ? lambda_enum() : -1;
            int pc = prime_closure();
            int pe = (nn <= 18) ? prime_enum() : -1;
            int expP = (nn >= 2 * s + 3);
            if (lf != s || (le >= 0 && le != s)) { bad++; printf(" [LAMBDA FAIL n=%d %d %d]", nn, lf, le); }
            if (nn >= 2 * s + 2) {
                if (pc != expP || (pe >= 0 && pe != expP)) { bad++; printf(" [PRIME FAIL n=%d %d %d]", nn, pc, pe); }
            }
            if (nn == 2 * s + 2 && !is_module_pair(0, 2 * s + 1)) { bad++; printf(" [MODULE {0,2s+1} FAIL]"); }
            if (nn == 2 * s + 1) printf(" n=%d(rot: lam %d prime %d)", nn, lf, pc);
        }
        printf(" checked n=%d..24\n", 2 * s + 1);
    }
    /* W */
    int wb = 0, wt = 0;
    for (int a = 1; a <= 25; a++) for (int b = 1; b <= 25; b++) {
        if (gcd(a, b) != 1) continue;
        for (int nn = a + b; nn <= 80; nn++) {
            struct walk k = {a, b};
            wt++;
            if (!strongly_connected_arcs(nn, walk_arc, &k)) wb++;
        }
        /* necessity check: n = a+b-1 should fail (sanity of the threshold), when a+b-1 >= 2 */
    }
    printf("Lemma W: %d instances, %d failures\n", wt, wb);
    bad += wb;
    /* threshold sharpness of W: for n = a + b - 1 the walk digraph is NOT strongly connected? */
    int sharp = 0, sharpt = 0;
    for (int a = 1; a <= 25; a++) for (int b = 1; b <= 25; b++) {
        if (gcd(a, b) != 1 || a + b - 1 < 2) continue;
        struct walk k = {a, b};
        sharpt++;
        if (!strongly_connected_arcs(a + b - 1, walk_arc, &k)) sharp++;
    }
    printf("Lemma W threshold: at n = a+b-1, %d of %d coprime pairs are not strongly connected\n", sharp, sharpt);
    /* D construction with a coprime bijection */
    for (int s = 1; s <= 8; s++) {
        S_ = s;
        memset(usedv, 0, sizeof usedv);
        if (!findsig(1)) { printf("no coprime bijection s=%d\n", s); bad++; continue; }
        int dbad = 0;
        for (int nn = 3 * s; nn <= 40; nn++) {
            band(s, nn);
            static int used[MAXN][MAXN];
            memset(used, 0, sizeof used);
            for (int r = 1; r <= s; r++) {
                struct walk k = {sig[r], r};
                /* arcs inside band and disjoint */
                for (int v = 0; v < nn; v++) {
                    int ws[2] = {v + sig[r], v - r};
                    for (int t = 0; t < 2; t++) {
                        int w = ws[t];
                        if (w < 0 || w >= nn) continue;
                        if (!A[v][w]) dbad++;
                        if (used[v][w]) dbad++;
                        used[v][w] = 1;
                    }
                }
                if (!strongly_connected_arcs(nn, walk_arc, &k)) dbad++;
            }
        }
        printf("D-construction s=%d sigma=", s);
        for (int r = 1; r <= s; r++) printf("%d%s", sig[r], r < s ? "," : "");
        printf(" n=%d..40: %d failures\n", 3 * s, dbad);
        bad += dbad;
    }
    printf("TOTAL FAILURES %d\n", bad);
    return 0;
}
