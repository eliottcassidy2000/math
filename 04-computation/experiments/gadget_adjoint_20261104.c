/*
 * gadget_adjoint_20261104.c  (opus-2026-11-04-S78)
 *
 * One-arc gadget transfer. For a tournament G with marks u (entry) and z (exit), and a rooted base (T, v), let
 *     N = T => G with u -> v reversed, rooted at z.
 * Then (h, e, R)(N, z) = M(G, u, z) (h, e, R)(T, v), where (h, e, R) = (H, paths ending at the root, root response) and
 *     M = [[h_G,  -s_u,       E_u      ],
 *          [e_z,  -P_{u->z},  E_{u,z}  ],
 *          [R_z,  -J_1,       J_2 + K  ]]
 * with  E_u = sum_{B ∋ u} e_u(G[B]) h(G - B),  E_{u,z} = sum_{B ∋ u, B ∌ z} e_u(G[B]) e_z(G - B),
 *       J_1 = sum_{S ∋ u, S ∌ z} s_u(G[S]) s_z(G - S),  J_2 = sum_{S ∋ u, S ∌ z} E_u(G[S]) s_z(G - S),
 *       K = sum over partitions G = S + Y + Z with u, z in Y of h(G[S]) P_{z->u}(G[Y]) h(G[Z]).
 * Checks:
 *   (1) the closed formula against the direct count of N for three probe bases (single vertex; TT_2 rooted at its source;
 *       C_3) and two further bases, for every labelled G of order <= MCHK and all (u, z);
 *   (2) the Lorentzian adjoint M(G^op, z, u) = J M(G, u, z)^T J, J = diag(1, -1, 1), for every labelled G of order <= MMAX.
 * Usage: gadj MMAX MCHK
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef long long ll;
#define MX 12
typedef struct { int m; uint32_t out[MX]; ll hS[1 << 8], eS[1 << 8][8], sS[1 << 8][8], P[1 << 8][8][8]; } GT;
static void tables(GT *G) {
    int m = G->m, full = (1 << m) - 1;
    memset(G->P, 0, sizeof(G->P));
    for (int a = 0; a < m; a++) G->P[1 << a][a][a] = 1;
    for (int S = 1; S <= full; S++) for (int a = 0; a < m; a++) if (S >> a & 1) for (int b = 0; b < m; b++) { ll c = G->P[S][a][b]; if (!c) continue;
        uint32_t nb = G->out[b] & ~(uint32_t)S; while (nb) { int w = __builtin_ctz(nb); nb &= nb - 1; G->P[S | (1 << w)][a][w] += c; } }
    G->hS[0] = 1; for (int a = 0; a < m; a++) { G->eS[0][a] = 0; G->sS[0][a] = 0; }
    for (int S = 1; S <= full; S++) { ll t = 0; for (int a = 0; a < m; a++) { G->eS[S][a] = 0; G->sS[S][a] = 0; }
        for (int a = 0; a < m; a++) for (int b = 0; b < m; b++) { ll c = G->P[S][a][b]; t += c; G->sS[S][a] += c; G->eS[S][b] += c; }
        G->hS[S] = t; }
}
static ll Eu_of(GT *G, int S, int u) { ll t = 0; for (int B = S; ; B = (B - 1) & S) { if (B >> u & 1) t += G->eS[B][u] * G->hS[S & ~B]; if (B == 0) break; } return t; }
static void matrix(GT *G, int u, int z, ll M[3][3]) {
    int m = G->m, full = (1 << m) - 1;
    ll hG = G->hS[full], su = G->sS[full][u], Eu = Eu_of(G, full, u), ez = G->eS[full][z], Puz = G->P[full][u][z], Euz = 0, Rz = 0, J1 = 0, J2 = 0, K = 0;
    for (int B1 = 0; B1 <= full; B1++) if ((B1 >> u & 1) && !(B1 >> z & 1)) Euz += G->eS[B1][u] * G->eS[full & ~B1][z];
    for (int S = 0; S <= full; S++) if (!(S >> z & 1)) Rz += G->hS[S] * G->sS[full & ~S][z];
    for (int S = 0; S <= full; S++) if ((S >> u & 1) && !(S >> z & 1)) { J1 += G->sS[S][u] * G->sS[full & ~S][z]; J2 += Eu_of(G, S, u) * G->sS[full & ~S][z]; }
    int rest = full & ~(1 << z) & ~(1 << u);
    for (int S = rest; ; S = (S - 1) & rest) { int Rm = full & ~S, Rr = Rm & ~(1 << z) & ~(1 << u);
        for (int X = Rr; ; X = (X - 1) & Rr) { int Y = X | (1 << z) | (1 << u); ll p = G->P[Y][z][u]; if (p) K += G->hS[S] * p * G->hS[Rm & ~Y]; if (X == 0) break; }
        if (S == 0) break; }
    M[0][0] = hG; M[0][1] = -su; M[0][2] = Eu; M[1][0] = ez; M[1][1] = -Puz; M[1][2] = Euz; M[2][0] = Rz; M[2][1] = -J1; M[2][2] = J2 + K;
}
/* direct (h, e_r, R_r) of a tournament with up to MX vertices */
static void triple(int n, const uint32_t *out, int r, ll *h, ll *e, ll *R) {
    int full = (1 << n) - 1; size_t sz = (size_t)(full + 1) * n;
    ll *dp = calloc(sz, sizeof(ll)), *st = calloc(sz, sizeof(ll)), *Hs = malloc(sizeof(ll) * (full + 1));
    for (int x = 0; x < n; x++) { dp[(size_t)(1 << x) * n + x] = 1; st[(size_t)(1 << x) * n + x] = 1; }
    uint32_t inn[MX]; for (int v = 0; v < n; v++) { inn[v] = 0; for (int w = 0; w < n; w++) if (out[w] >> v & 1) inn[v] |= 1u << w; }
    for (int S = 1; S <= full; S++) for (int x = 0; x < n; x++) {
        ll c = dp[(size_t)S * n + x]; if (c) { uint32_t nb = out[x] & ~(uint32_t)S; while (nb) { int y = __builtin_ctz(nb); nb &= nb - 1; dp[(size_t)(S | (1 << y)) * n + y] += c; } }
        ll d = st[(size_t)S * n + x]; if (d) { uint32_t nb = inn[x] & ~(uint32_t)S; while (nb) { int y = __builtin_ctz(nb); nb &= nb - 1; st[(size_t)(S | (1 << y)) * n + y] += d; } } }
    Hs[0] = 1; for (int S = 1; S <= full; S++) { ll t = 0; for (int x = 0; x < n; x++) t += dp[(size_t)S * n + x]; Hs[S] = t; }
    *h = Hs[full]; *e = dp[(size_t)full * n + r]; ll t = 0;
    for (int S = 0; S <= full; S++) if (!(S >> r & 1)) t += Hs[S] * st[(size_t)(full & ~S) * n + r];
    *R = t; free(dp); free(st); free(Hs);
}
int main(int argc, char **argv) {
    int MMAX = argc > 1 ? atoi(argv[1]) : 5, MCHK = argc > 2 ? atoi(argv[2]) : 5;
    /* probe bases: (order, out, port) */
    uint32_t pb[5][4] = {{0}, {2, 0}, {2, 4, 1}, {6, 4, 0}, {0, 0, 0, 0}};
    int pn[5] = {1, 2, 3, 3, 4}, pv[5] = {0, 0, 0, 1, 2};
    { /* order-4 probe: 0->1, 0->2, 1->2, 3->0, 1->3, 2->3 */ pb[4][0] = 2 | 4; pb[4][1] = 4 | 8; pb[4][2] = 8; pb[4][3] = 1; }
    ll ptr[5][3]; for (int i = 0; i < 5; i++) triple(pn[i], pb[i], pv[i], &ptr[i][0], &ptr[i][1], &ptr[i][2]);
    printf("# gadget_adjoint_20261104 MMAX=%d MCHK=%d; probe triples:", MMAX, MCHK);
    for (int i = 0; i < 5; i++) printf(" (%lld,%lld,%lld)", ptr[i][0], ptr[i][1], ptr[i][2]); printf("\n");
    GT *G = malloc(sizeof(GT)), *Gop = malloc(sizeof(GT));
    for (int m = 1; m <= MMAX; m++) {
        int na = m * (m - 1) / 2; long cases = 0, fadj = 0, fform = 0, chk = 0;
        for (long code = 0; code < (1L << na); code++) {
            G->m = m; Gop->m = m; memset(G->out, 0, sizeof G->out); memset(Gop->out, 0, sizeof Gop->out); int k = 0;
            for (int i = 0; i < m; i++) for (int j = i + 1; j < m; j++, k++) { if (code >> k & 1) G->out[i] |= 1u << j; else G->out[j] |= 1u << i; }
            for (int i = 0; i < m; i++) for (int j = 0; j < m; j++) if (G->out[i] >> j & 1) Gop->out[j] |= 1u << i;
            tables(G); tables(Gop);
            for (int u = 0; u < m; u++) for (int z = 0; z < m; z++) {
                ll M[3][3], Mo[3][3]; matrix(G, u, z, M); matrix(Gop, z, u, Mo); cases++;
                int sg[3] = {1, -1, 1};
                for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) if (Mo[i][j] != sg[i] * M[j][i] * sg[j]) { fadj++; i = 3; break; }
                if (m <= MCHK) for (int p = 0; p < 5; p++) {
                    int n1 = pn[p], n = n1 + m; uint32_t o[MX]; memset(o, 0, sizeof o);
                    for (int a = 0; a < n1; a++) o[a] = pb[p][a];
                    for (int a = 0; a < m; a++) o[n1 + a] = G->out[a] << n1;
                    for (int a = 0; a < n1; a++) for (int b = 0; b < m; b++) { if (a == pv[p] && b == u) o[n1 + b] |= 1u << a; else o[a] |= 1u << (n1 + b); }
                    ll h, e, R; triple(n, o, n1 + z, &h, &e, &R); chk++;
                    ll w[3]; for (int i = 0; i < 3; i++) w[i] = M[i][0] * ptr[p][0] + M[i][1] * ptr[p][1] + M[i][2] * ptr[p][2];
                    if (w[0] != h || w[1] != e || w[2] != R) fform++;
                }
            }
        }
        printf("order %d: (G,u,z) cases %ld | adjoint failures %ld | formula-vs-direct checks %ld, failures %ld\n", m, cases, fadj, chk, fform);
        fflush(stdout);
    }
    return 0;
}
