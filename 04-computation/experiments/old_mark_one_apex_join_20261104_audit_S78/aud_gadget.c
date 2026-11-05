/*
 * aud_gadget.c  (S78 audit, independent code)
 * For every labelled tournament G of order m <= MFIT and every (u, z):
 *   - FIT M(G,u,z) from direct counts on three spanning probe bases (single vertex (1,1,1); TT_2 at its source (1,0,2);
 *     C_3 (3,1,3)), solving the 3x3 system exactly over Q (Cramer, det 2);
 *   - TEST linearity on further bases: every labelled rooted base of order <= 4 (75 bases) when m <= 4, and 14 fixed
 *     bases of orders 2..6 (distinct triples) when m = 5;
 *   - COMPARE the fitted M with the auditor's own implementation of the closed formula in the S78 note, section 2;
 *   - ADJOINT: fitted M(G^op, z, u) against J M(G,u,z)^T J, J = diag(1,-1,1);
 *   - determinant distribution of M.
 * For m = MFIT+1 (default 6): closed-formula M only, adjoint and determinants.
 * N = T => G with u -> v reversed; rooted at z.
 * Usage: aud_gadget MFIT
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>
typedef long long ll;
#define MXV 12

static void triple(int n, const uint32_t *out, int r, ll *h, ll *e, ll *R) {
    int full = (1 << n) - 1;
    static __thread ll *buf = NULL;
    if (!buf) buf = malloc(sizeof(ll) * ((size_t)(1 << MXV) * MXV * 2 + (1 << MXV)));
    ll (*end)[MXV] = (ll (*)[MXV])buf, (*st)[MXV] = (ll (*)[MXV])(buf + (size_t)(1 << MXV) * MXV), *hs = buf + (size_t)(1 << MXV) * MXV * 2;
    hs[0] = 1;
    for (int S = 1; S <= full; S++) {
        ll t = 0;
        for (int v = 0; v < n; v++) {
            if (!((S >> v) & 1)) { end[S][v] = st[S][v] = 0; continue; }
            int P = S & ~(1 << v);
            if (!P) { end[S][v] = st[S][v] = 1; t++; continue; }
            ll a = 0, b = 0;
            for (int u = 0; u < n; u++) if ((P >> u) & 1) { if ((out[u] >> v) & 1) a += end[P][u]; if ((out[v] >> u) & 1) b += st[P][u]; }
            end[S][v] = a; st[S][v] = b; t += a;
        }
        hs[S] = t;
    }
    *h = hs[full]; *e = end[full][r]; ll s = 0;
    for (int X = 0; X <= full; X++) if (!((X >> r) & 1)) s += hs[X] * st[full & ~X][r];
    *R = s;
}
static void tour(int n, long code, uint32_t *out) {
    memset(out, 0, sizeof(uint32_t) * MXV); int k = 0;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++, k++) { if ((code >> k) & 1) out[i] |= 1u << j; else out[j] |= 1u << i; }
}
/* N = T => G with u -> v reversed; T on 0..nt-1 (port v), G on nt..nt+m-1 */
static void transfer(int nt, const uint32_t *T, int v, int m, const uint32_t *G, int u, int z, ll *o3) {
    uint32_t o[MXV]; memset(o, 0, sizeof o);
    for (int i = 0; i < nt; i++) o[i] = T[i];
    for (int i = 0; i < m; i++) o[nt + i] = G[i] << nt;
    for (int i = 0; i < nt; i++) for (int j = 0; j < m; j++) { if (i == v && j == u) o[nt + j] |= 1u << i; else o[i] |= 1u << (nt + j); }
    triple(nt + m, o, nt + z, &o3[0], &o3[1], &o3[2]);
}
/* ---- closed formula (auditor's implementation of the note's definitions) ---- */
typedef struct { int m; ll P[1 << 7][7][7]; ll hS[1 << 7], eS[1 << 7][7], sS[1 << 7][7]; } GT;
static void gtables(GT *g, int m, const uint32_t *out) {
    g->m = m; int full = (1 << m) - 1;
    memset(g->P, 0, sizeof g->P);
    for (int a = 0; a < m; a++) g->P[1 << a][a][a] = 1;
    for (int S = 1; S <= full; S++) for (int a = 0; a < m; a++) for (int b = 0; b < m; b++) {
        if (!((S >> a) & 1) || !((S >> b) & 1) || (a == b && S != (1 << a))) continue;
        if (a == b) continue;
        /* paths a .. b covering S: last step c -> b with c in S - b, path a..c covering S - b */
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
    ll hG = g->hS[full], su = g->sS[full][u], Eu = EU(g, full, u), ez = g->eS[full][z], Puz = g->P[full][u][z];
    ll Euz = 0, Rz = 0, J1 = 0, J2 = 0, K = 0;
    for (int B = 0; B <= full; B++) if (((B >> u) & 1) && !((B >> z) & 1)) { Euz += g->eS[B][u] * g->eS[full & ~B][z]; J1 += g->sS[B][u] * g->sS[full & ~B][z]; J2 += EU(g, B, u) * g->sS[full & ~B][z]; }
    for (int X = 0; X <= full; X++) if (!((X >> z) & 1)) Rz += g->hS[X] * g->sS[full & ~X][z];
    /* K: ordered (S, Y, Z) partition, u, z in Y */
    for (int Y = 1; Y <= full; Y++) if (((Y >> u) & 1) && ((Y >> z) & 1)) {
        ll p = g->P[Y][z][u]; if (!p) continue;
        int rest = full & ~Y;
        for (int S = rest; ; S = (S - 1) & rest) { K += g->hS[S] * p * g->hS[rest & ~S]; if (!S) break; }
    }
    M[0][0] = hG; M[0][1] = -su; M[0][2] = Eu;
    M[1][0] = ez; M[1][1] = -Puz; M[1][2] = Euz;
    M[2][0] = Rz; M[2][1] = -J1; M[2][2] = J2 + K;
}
static ll det3(ll M[3][3]) { return M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]); }

int main(int argc, char **argv) {
    int MFIT = argc > 1 ? atoi(argv[1]) : 5;
    /* probes */
    uint32_t P1[MXV] = {0}, P2[MXV] = {0}, P3[MXV] = {0};
    P2[0] = 2;                       /* 0 -> 1, root 0 (source) */
    P3[0] = 2; P3[1] = 4; P3[2] = 1; /* 0->1->2->0 */
    ll t[3][3]; triple(1, P1, 0, &t[0][0], &t[0][1], &t[0][2]); triple(2, P2, 0, &t[1][0], &t[1][1], &t[1][2]); triple(3, P3, 0, &t[2][0], &t[2][1], &t[2][2]);
    printf("# aud_gadget MFIT=%d; probe triples (%lld,%lld,%lld) (%lld,%lld,%lld) (%lld,%lld,%lld)\n", MFIT, t[0][0], t[0][1], t[0][2], t[1][0], t[1][1], t[1][2], t[2][0], t[2][1], t[2][2]);
    /* A = matrix with columns t_k; M A = O (columns outputs); M = O A^{-1}; A^{-1} = adj(A)/det(A) */
    ll A[3][3]; for (int i = 0; i < 3; i++) for (int k = 0; k < 3; k++) A[i][k] = t[k][i];
    ll dA = det3(A), adj[3][3];
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) {
        int r0 = (j + 1) % 3, r1 = (j + 2) % 3, c0 = (i + 1) % 3, c1 = (i + 2) % 3;
        adj[i][j] = A[r0][c0] * A[r1][c1] - A[r0][c1] * A[r1][c0];
    }
    printf("probe determinant %lld\n", dA);
    /* extra test bases */
    typedef struct { int n; uint32_t out[MXV]; int r; ll tr[3]; } Base;
    Base *bases = malloc(sizeof(Base) * 200); int nbase = 0;
    for (int n = 1; n <= 4; n++) for (long code = 0; code < (1L << (n * (n - 1) / 2)); code++) {
        Base *b = &bases[nbase++]; b->n = n; tour(n, code, b->out); b->r = 0; triple(n, b->out, 0, &b->tr[0], &b->tr[1], &b->tr[2]);
    }
    int n4 = nbase;
    /* 14 fixed bases of orders 5 and 6 with distinct triples (deterministic pseudo-random codes) */
    { unsigned long long s = 12345; int got = 0;
      while (got < 14) { s ^= s << 13; s ^= s >> 7; s ^= s << 17; int n = 5 + (int)(s % 2); long code = (long)((s >> 8) % (1ULL << (n * (n - 1) / 2)));
        Base b; b.n = n; tour(n, code, b.out); b.r = (int)((s >> 40) % n); triple(n, b.out, b.r, &b.tr[0], &b.tr[1], &b.tr[2]);
        int dup = 0; for (int i = n4; i < nbase; i++) if (bases[i].tr[0] == b.tr[0] && bases[i].tr[1] == b.tr[1] && bases[i].tr[2] == b.tr[2]) dup = 1;
        if (!dup) { bases[nbase++] = b; got++; } } }
    printf("test bases: %d of order <= 4 (all labelled, root 0) + %d of order 5-6\n", n4, nbase - n4);
    for (int m = 1; m <= MFIT + 1; m++) {
        int na = m * (m - 1) / 2; long ncodes = 1L << na;
        long cases = 0, nonint = 0, linfail = 0, linchk = 0, cfmis = 0, adjfail = 0; long detcnt[7] = {0}; /* det -3..3 */ long detother = 0;
        #pragma omp parallel for schedule(dynamic, 16) reduction(+:cases,nonint,linfail,linchk,cfmis,adjfail,detother)
        for (long code = 0; code < ncodes; code++) {
            uint32_t G[MXV], Gop[MXV]; tour(m, code, G); memset(Gop, 0, sizeof Gop);
            for (int i = 0; i < m; i++) for (int j = 0; j < m; j++) if ((G[i] >> j) & 1) Gop[j] |= 1u << i;
            GT *g = malloc(sizeof(GT)), *gop = malloc(sizeof(GT)); gtables(g, m, G); gtables(gop, m, Gop);
            for (int u = 0; u < m; u++) for (int z = 0; z < m; z++) {
                cases++;
                ll Mc[3][3], Mco[3][3]; closed(g, u, z, Mc); closed(gop, z, u, Mco);
                ll Mf[3][3]; int fitted = 0;
                if (m <= MFIT) {
                    ll O[3][3];   /* O[i][k] = output coordinate i on probe k */
                    ll o3[3];
                    transfer(1, P1, 0, m, G, u, z, o3); for (int i = 0; i < 3; i++) O[i][0] = o3[i];
                    transfer(2, P2, 0, m, G, u, z, o3); for (int i = 0; i < 3; i++) O[i][1] = o3[i];
                    transfer(3, P3, 0, m, G, u, z, o3); for (int i = 0; i < 3; i++) O[i][2] = o3[i];
                    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) {
                        ll num = 0; for (int k = 0; k < 3; k++) num += O[i][k] * adj[k][j];
                        if (num % dA) nonint++;
                        Mf[i][j] = num / dA;
                    }
                    fitted = 1;
                    /* linearity tests */
                    int b0 = (m <= 4) ? 0 : n4, b1 = nbase;
                    for (int bi = b0; bi < b1; bi++) {
                        Base *b = &bases[bi]; transfer(b->n, b->out, b->r, m, G, u, z, o3); linchk++;
                        for (int i = 0; i < 3; i++) { ll p = Mf[i][0] * b->tr[0] + Mf[i][1] * b->tr[1] + Mf[i][2] * b->tr[2]; if (p != o3[i]) { linfail++; break; } }
                    }
                    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) if (Mf[i][j] != Mc[i][j]) { cfmis++; i = 3; break; }
                }
                /* adjoint: use fitted matrices when available (fit G^op too), else closed */
                ll Ma[3][3], Mb[3][3];
                if (fitted) {
                    ll O[3][3], o3[3];
                    transfer(1, P1, 0, m, Gop, z, u, o3); for (int i = 0; i < 3; i++) O[i][0] = o3[i];
                    transfer(2, P2, 0, m, Gop, z, u, o3); for (int i = 0; i < 3; i++) O[i][1] = o3[i];
                    transfer(3, P3, 0, m, Gop, z, u, o3); for (int i = 0; i < 3; i++) O[i][2] = o3[i];
                    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) { ll num = 0; for (int k = 0; k < 3; k++) num += O[i][k] * adj[k][j]; Mb[i][j] = num / dA; }
                    memcpy(Ma, Mf, sizeof Ma);
                } else { memcpy(Ma, Mc, sizeof Ma); memcpy(Mb, Mco, sizeof Mb); }
                int sg[3] = {1, -1, 1};
                for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) if (Mb[i][j] != sg[i] * Ma[j][i] * sg[j]) { adjfail++; i = 3; break; }
                ll d = det3(Ma);
                if (d >= -3 && d <= 3) { __atomic_fetch_add(&detcnt[d + 3], 1, __ATOMIC_RELAXED); } else detother++;
            }
            free(g); free(gop);
        }
        printf("order %d: (G,u,z) %ld | %s | fitted-nonintegral %ld | linearity checks %ld failures %ld | fitted vs closed formula mismatches %ld | adjoint failures %ld (%s) | det: -1:%ld 0:%ld 1:%ld other:%ld\n",
               m, cases, m <= MFIT ? "fitted" : "closed only", nonint, linchk, linfail, cfmis, adjfail, m <= MFIT ? "fitted matrices" : "closed formula",
               detcnt[2], detcnt[3], detcnt[4], detother + detcnt[0] + detcnt[1] + detcnt[5] + detcnt[6]);
        fflush(stdout);
    }
    return 0;
}
