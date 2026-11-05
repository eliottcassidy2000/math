/*
 * aud_descent.c  (S78 audit, independent code)
 * Checks of the S78 descent-lemma application (note section 3):
 *   for every labelled gadget (G, u, z) of order <= MMAX (closed-formula matrix, validated by aud_gadget) and every
 *   realisable triple y = (h, e, R) of order <= 7 (bank file), h' = (M y)_1 satisfies h' >= R >= h, and
 *   h' = h  iff  e = h = R (the base port is a sink)  and  G is transitive with u its source.
 *   Also: the determinant distribution of M, and the number of unimodular (det +-1) gadgets.
 * Usage: aud_descent MMAX bankfile
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef long long ll;
typedef struct { int m; ll P[1 << 7][7][7]; ll hS[1 << 7], eS[1 << 7][7], sS[1 << 7][7]; } GT;
static void tour(int n, long code, uint32_t *out) {
    memset(out, 0, sizeof(uint32_t) * 8); int k = 0;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++, k++) { if ((code >> k) & 1) out[i] |= 1u << j; else out[j] |= 1u << i; }
}
static void gtables(GT *g, int m, const uint32_t *out) {
    g->m = m; int full = (1 << m) - 1;
    memset(g->P, 0, sizeof g->P);
    for (int a = 0; a < m; a++) g->P[1 << a][a][a] = 1;
    for (int S = 1; S <= full; S++) for (int a = 0; a < m; a++) for (int b = 0; b < m; b++) {
        if (a == b || !((S >> a) & 1) || !((S >> b) & 1)) continue;
        int P = S & ~(1 << b); ll t = 0;
        for (int c = 0; c < m; c++) if (((P >> c) & 1) && ((out[c] >> b) & 1)) t += g->P[P][a][c];
        g->P[S][a][b] = t;
    }
    g->hS[0] = 1; for (int a = 0; a < m; a++) { g->eS[0][a] = g->sS[0][a] = 0; }
    for (int S = 1; S <= full; S++) {
        ll t = 0; for (int a = 0; a < m; a++) { g->eS[S][a] = g->sS[S][a] = 0; }
        for (int a = 0; a < m; a++) for (int b = 0; b < m; b++) { ll c = g->P[S][a][b]; t += c; g->sS[S][a] += c; g->eS[S][b] += c; }
        g->hS[S] = t;
    }
}
static ll EU(const GT *g, int S, int u) { ll t = 0; for (int B = S; ; B = (B - 1) & S) { if ((B >> u) & 1) t += g->eS[B][u] * g->hS[S & ~B]; if (!B) break; } return t; }
static void closed(const GT *g, int u, int z, ll M[3][3]) {
    int m = g->m, full = (1 << m) - 1;
    ll Euz = 0, Rz = 0, J1 = 0, J2 = 0, K = 0;
    for (int B = 0; B <= full; B++) if (((B >> u) & 1) && !((B >> z) & 1)) { Euz += g->eS[B][u] * g->eS[full & ~B][z]; J1 += g->sS[B][u] * g->sS[full & ~B][z]; J2 += EU(g, B, u) * g->sS[full & ~B][z]; }
    for (int X = 0; X <= full; X++) if (!((X >> z) & 1)) Rz += g->hS[X] * g->sS[full & ~X][z];
    for (int Y = 1; Y <= full; Y++) if (((Y >> u) & 1) && ((Y >> z) & 1)) {
        ll p = g->P[Y][z][u]; if (!p) continue; int rest = full & ~Y;
        for (int S = rest; ; S = (S - 1) & rest) { K += g->hS[S] * p * g->hS[rest & ~S]; if (!S) break; }
    }
    M[0][0] = g->hS[full]; M[0][1] = -g->sS[full][u]; M[0][2] = EU(g, full, u);
    M[1][0] = g->eS[full][z]; M[1][1] = -g->P[full][u][z]; M[1][2] = Euz;
    M[2][0] = Rz; M[2][1] = -J1; M[2][2] = J2 + K;
}
static ll det3(ll M[3][3]) { return M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]); }
int main(int argc, char **argv) {
    int MMAX = argc > 1 ? atoi(argv[1]) : 5; const char *fn = argc > 2 ? argv[2] : "aud_bank7.txt";
    FILE *f = fopen(fn, "r"); int nb = 0, cap = 4096; ll (*Y)[3] = malloc(sizeof(ll[3]) * cap); ll a, b, c;
    while (fscanf(f, "%lld %lld %lld", &a, &b, &c) == 3) { Y[nb][0] = a; Y[nb][1] = b; Y[nb][2] = c; nb++; }
    fclose(f);
    printf("# aud_descent MMAX=%d bank=%s (%d triples)\n", MMAX, fn, nb);
    GT *g = malloc(sizeof(GT));
    for (int m = 1; m <= MMAX; m++) {
        long cases = 0, viol = 0, eqcases = 0, eqmis = 0, nuni = 0, nsing = 0, nother = 0, ndetpos = 0;
        long dets[41] = {0};
        for (long code = 0; code < (1L << (m * (m - 1) / 2)); code++) {
            uint32_t G[8]; tour(m, code, G); gtables(g, m, G);
            int trans = 1; uint32_t seen = 0; for (int v = 0; v < m; v++) { int d = __builtin_popcount(G[v]); if ((seen >> d) & 1) trans = 0; seen |= 1u << d; }
            for (int u = 0; u < m; u++) for (int z = 0; z < m; z++) {
                ll M[3][3]; closed(g, u, z, M); cases++;
                ll d = det3(M); if (d == 1 || d == -1) nuni++; else if (d == 0) nsing++; else nother++;
                if (d > 0) ndetpos++;
                if (d >= -20 && d <= 20) dets[d + 20]++;
                int usrc = (__builtin_popcount(G[u]) == m - 1);
                for (int i = 0; i < nb; i++) {
                    ll hp = M[0][0] * Y[i][0] + M[0][1] * Y[i][1] + M[0][2] * Y[i][2];
                    if (hp < Y[i][2] || Y[i][2] < Y[i][0]) viol++;
                    int eq = (hp == Y[i][0]);
                    int pred = (Y[i][1] == Y[i][0] && Y[i][2] == Y[i][0] && trans && usrc);
                    if (eq) eqcases++;
                    if (eq != pred) eqmis++;
                }
            }
        }
        printf("order %d: gadgets (G,u,z) %ld | h' >= R >= h violations %ld (of %ld pairs) | h' = h cases %ld, mismatches with [base port a sink, G transitive, u its source] %ld | det: unimodular %ld, singular %ld, |det|>=2 %ld\n",
               m, cases, viol, cases * nb, eqcases, eqmis, nuni, nsing, nother);
        printf("   det histogram (-20..20, nonzero entries):");
        for (int d = 0; d <= 40; d++) if (dets[d]) printf(" %d:%ld", d - 20, dets[d]);
        printf("\n");
        fflush(stdout);
    }
    return 0;
}
