/*
 * aud_oldmark.c  (S78 audit, independent code)
 * Tournament classes of order 1..NMAX by canonical augmentation with the auditor's own canonical form (ordered
 * invariant partition refined by out-neighbour cell counts, then the minimum code over all cell-respecting
 * permutations). Class counts are compared with OEIS A000568.
 * For every class B and every vertex a of its representative, T1 = C2(B, a) (x -> a only, B - a -> x; B -> y; y -> x).
 * Checks, with failure counts:
 *   D     sigma_v >= d(d-1) at every vertex of B and of T1;
 *   D+    sigma_v >= E_v + 2 sum_{pairs p,q in N+(v)} h(T - {v,p,q}) at every vertex of T1  (auditor's strengthening;
 *         E_v = Hamiltonian paths of T - v ending in N+(v));
 *   nat   sigma_x = 0, sigma_y = h - e_a, H(T1) = h + 2R_a;
 *   sink  sigma_z = 3h(B') + 2R_a(B') - e_a(B') for the sink z != a of B, B' = B - z;
 *   min4  every old w != a has sigma_w >= 4; sigma_w = 4 iff B transitive, a penultimate, w the sink (both directions);
 *   ns9   every old w != a with d_B(w) >= 1 has sigma_w >= 9 (auditor's upgrade; the session proves 6);
 *   mark  RESULT (12): sigma_a = 2R_a(B) - h(B) - h(B - a) for every mark; a sink -> 0; a -> w only -> e_w(U) + R_w(U).
 * Also: minimum digits, and the (B, a, w) attaining 9 at non-sink old vertices (counted, first few printed).
 * Usage: aud_oldmark NMAX
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>
typedef long long ll;
typedef unsigned long long ull;
#define MX 10

/* ---------- canonical form ---------- */
static ull code_of(int n, const uint32_t *out, const int *perm) {     /* bit for each pair (i<j) in new labels: perm[i] -> perm[j] */
    ull c = 0; for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) c = (c << 1) | ((out[perm[i]] >> perm[j]) & 1u);
    return c;
}
static ull best_code; static int cells[MX], ncell, cellstart[MX + 1], cellv[MX];
static void perm_rec(int n, const uint32_t *out, int *perm, int pos, uint32_t used) {
    if (pos == n) { ull c = code_of(n, out, perm); if (c < best_code) best_code = c; return; }
    int ci = 0; while (cellstart[ci + 1] <= pos) ci++;
    for (int k = cellstart[ci]; k < cellstart[ci + 1]; k++) { int v = cellv[k]; if ((used >> v) & 1) continue; perm[pos] = v; perm_rec(n, out, perm, pos + 1, used | (1u << v)); }
}
static ull canon(int n, const uint32_t *out) {
    /* initial colour: score; refine by counts of out-neighbours per colour (ordered) */
    int col[MX], ncol;
    ll key[MX];
    for (int v = 0; v < n; v++) key[v] = __builtin_popcount(out[v]);
    for (int iter = 0; iter < 2 * n + 2; iter++) {
        /* rank keys */
        ll uniq[MX]; int nu = 0;
        for (int v = 0; v < n; v++) { int f = 0; for (int i = 0; i < nu; i++) if (uniq[i] == key[v]) f = 1; if (!f) uniq[nu++] = key[v]; }
        for (int i = 1; i < nu; i++) { ll x = uniq[i]; int j = i - 1; while (j >= 0 && uniq[j] > x) { uniq[j + 1] = uniq[j]; j--; } uniq[j + 1] = x; }
        for (int v = 0; v < n; v++) for (int i = 0; i < nu; i++) if (uniq[i] == key[v]) col[v] = i;
        ncol = nu;
        ll nk[MX];
        for (int v = 0; v < n; v++) {
            ll k2 = col[v];
            for (int c = 0; c < ncol; c++) { int cnt = 0; for (int u = 0; u < n; u++) if (col[u] == c && ((out[v] >> u) & 1)) cnt++; k2 = k2 * 16 + cnt; }
            nk[v] = k2;
        }
        int nu2 = 0; ll uq2[MX];
        for (int v = 0; v < n; v++) { int f = 0; for (int i = 0; i < nu2; i++) if (uq2[i] == nk[v]) f = 1; if (!f) uq2[nu2++] = nk[v]; }
        for (int v = 0; v < n; v++) key[v] = nk[v];
        if (nu2 == ncol) break;
    }
    /* final ranking of keys -> ordered cells */
    ll uniq[MX]; int nu = 0;
    for (int v = 0; v < n; v++) { int f = 0; for (int i = 0; i < nu; i++) if (uniq[i] == key[v]) f = 1; if (!f) uniq[nu++] = key[v]; }
    for (int i = 1; i < nu; i++) { ll x = uniq[i]; int j = i - 1; while (j >= 0 && uniq[j] > x) { uniq[j + 1] = uniq[j]; j--; } uniq[j + 1] = x; }
    int p = 0; ncell = nu;
    for (int i = 0; i < nu; i++) { cellstart[i] = p; for (int v = 0; v < n; v++) if (key[v] == uniq[i]) cellv[p++] = v; }
    cellstart[nu] = n; (void)cells;
    best_code = ~0ULL; int perm[MX];
    perm_rec(n, out, perm, 0, 0);
    return best_code;
}
static void decode(int n, ull code, uint32_t *out) {
    for (int i = 0; i < n; i++) out[i] = 0;
    int tot = n * (n - 1) / 2, b = tot - 1;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if ((code >> b) & 1) out[i] |= 1u << j; else out[j] |= 1u << i; b--; }
}
/* hash set of codes */
typedef struct { ull *k; char *u; size_t cap; } HS;
static int hs_add(HS *h, ull key) { size_t i = (size_t)((key * 0x9E3779B97F4A7C15ULL) >> 20) % h->cap; while (h->u[i]) { if (h->k[i] == key) return 0; i = (i + 1) % h->cap; } h->u[i] = 1; h->k[i] = key; return 1; }

/* ---------- path tables ---------- */
typedef struct { int n; ll hS[1 << MX], end[1 << MX][MX], st[1 << MX][MX]; } PT;
static void tables(PT *T, int n, const uint32_t *out) {
    T->n = n; int full = (1 << n) - 1; T->hS[0] = 1;
    for (int v = 0; v < n; v++) { T->end[0][v] = 0; T->st[0][v] = 0; }
    for (int S = 1; S <= full; S++) {
        ll t = 0;
        for (int v = 0; v < n; v++) {
            if (!((S >> v) & 1)) { T->end[S][v] = 0; T->st[S][v] = 0; continue; }
            int P = S & ~(1 << v);
            if (!P) { T->end[S][v] = 1; T->st[S][v] = 1; t++; continue; }
            ll a = 0, b = 0;
            for (int u = 0; u < n; u++) if ((P >> u) & 1) { if ((out[u] >> v) & 1) a += T->end[P][u]; if ((out[v] >> u) & 1) b += T->st[P][u]; }
            T->end[S][v] = a; T->st[S][v] = b; t += a;
        }
        T->hS[S] = t;
    }
}
static ll Rsub(const PT *T, int S, int v) { ll s = 0; int rest = S & ~(1 << v); for (int X = rest; ; X = (X - 1) & rest) { s += T->hS[X] * T->st[S & ~X][v]; if (!X) break; } return s; }

int main(int argc, char **argv) {
    int NMAX = argc > 1 ? atoi(argv[1]) : 8;
    const long A000568[] = {1, 1, 1, 2, 4, 12, 56, 456, 6880, 191536};
    ull *cur = malloc(sizeof(ull)); size_t ncur = 1; { uint32_t o[1] = {0}; cur[0] = canon(1, o); }
    printf("# aud_oldmark NMAX=%d\n", NMAX);
    for (int n = 1; n <= NMAX; n++) {
        if (n > 1) {
            HS h; h.cap = ncur * (1u << (n - 1)) * 2 + 1024; h.k = malloc(sizeof(ull) * h.cap); h.u = calloc(h.cap, 1);
            size_t cap = 1024, nn = 0; ull *nx = malloc(sizeof(ull) * cap);
            for (size_t i = 0; i < ncur; i++) {
                uint32_t out[MX]; decode(n - 1, cur[i], out);
                for (uint32_t O = 0; O < (1u << (n - 1)); O++) {
                    uint32_t o2[MX]; for (int v = 0; v < n - 1; v++) o2[v] = out[v] | (((O >> v) & 1) ? 0u : (1u << (n - 1)));
                    o2[n - 1] = O;
                    ull c = canon(n, o2);
                    if (hs_add(&h, c)) { if (nn == cap) { cap *= 2; nx = realloc(nx, sizeof(ull) * cap); } nx[nn++] = c; }
                }
            }
            free(h.k); free(h.u); free(cur); cur = nx; ncur = nn;
        }
        long fD = 0, fDp = 0, fnat = 0, fsink = 0, fmin4 = 0, feq = 0, fns9 = 0, fR12 = 0, fa0 = 0, fa1 = 0, pairs = 0, sinkc = 0, n9 = 0, n4 = 0, n9other = 0;
        ll minold = -1, minns = -1, mina1 = -1;
        #pragma omp parallel reduction(+:fD,fDp,fnat,fsink,fmin4,feq,fns9,fR12,fa0,fa1,pairs,sinkc,n9,n4,n9other)
        {
            PT *TB = malloc(sizeof(PT)), *T1 = malloc(sizeof(PT));
            ll lminold = -1, lminns = -1, lmina1 = -1;
            #pragma omp for schedule(dynamic, 4)
            for (size_t i = 0; i < ncur; i++) {
                uint32_t B[MX]; decode(n, cur[i], B);
                tables(TB, n, B);
                int full = (1 << n) - 1; ll h = TB->hS[full];
                int deg[MX], trans = 1; uint32_t seen = 0;
                for (int v = 0; v < n; v++) { deg[v] = __builtin_popcount(B[v]); if ((seen >> deg[v]) & 1) trans = 0; seen |= 1u << deg[v]; }
                for (int v = 0; v < n; v++) { ll sg = Rsub(TB, full, v) - h; if (sg < (ll)deg[v] * (deg[v] - 1)) fD++; }
                for (int a = 0; a < n; a++) {
                    pairs++;
                    int N2 = n + 2, x = n, y = n + 1; uint32_t o[MX];
                    for (int v = 0; v < n; v++) o[v] = B[v] | (1u << y) | (v != a ? (1u << x) : 0u);
                    o[x] = 1u << a; o[y] = 1u << x;
                    tables(T1, N2, o);
                    int F = (1 << N2) - 1; ll H1 = T1->hS[F]; ll sg[MX];
                    for (int v = 0; v < N2; v++) sg[v] = Rsub(T1, F, v) - H1;
                    ll Ra = Rsub(TB, full, a);
                    if (sg[x] != 0 || sg[y] != h - TB->end[full][a] || H1 != h + 2 * Ra) fnat++;
                    for (int v = 0; v < N2; v++) {
                        int d = __builtin_popcount(o[v]); if (sg[v] < (ll)d * (d - 1)) fD++;
                        /* D+ */
                        int rest = F & ~(1 << v); ll Ev = 0; for (int p = 0; p < N2; p++) if ((o[v] >> p) & 1) Ev += T1->end[rest][p];
                        ll pr = 0; for (int p = 0; p < N2; p++) if ((o[v] >> p) & 1) for (int q = p + 1; q < N2; q++) if ((o[v] >> q) & 1) pr += T1->hS[rest & ~(1 << p) & ~(1 << q)];
                        if (sg[v] < Ev + 2 * pr) fDp++;
                    }
                    ll mo = -1;
                    for (int w = 0; w < n; w++) if (w != a) {
                        if (mo < 0 || sg[w] < mo) mo = sg[w];
                        if (deg[w] == 0) {
                            sinkc++; int Bp = full & ~(1 << w);
                            ll hp = TB->hS[Bp], ep = TB->end[Bp][a], Rp = Rsub(TB, Bp, a);
                            if (sg[w] != 3 * hp + 2 * Rp - ep) fsink++;
                        } else {
                            if (sg[w] < 9) fns9++;
                            if (sg[w] == 9) { n9++; if (!(trans && deg[a] == 0 && deg[w] == 1)) n9other++; }
                            if (trans && deg[a] == 0 && deg[w] == 1 && sg[w] != 9) n9other++;
                            if (lminns < 0 || sg[w] < lminns) lminns = sg[w];
                        }
                        int eqcase = trans && deg[w] == 0 && deg[a] == 1;
                        if (sg[w] == 4) n4++;
                        if ((sg[w] == 4) != eqcase) feq++;
                    }
                    if (n >= 2) { if (mo < 4) fmin4++; if (lminold < 0 || mo < lminold) lminold = mo; }
                    /* marks */
                    int Ua = full & ~(1 << a);
                    if (sg[a] != 2 * Ra - h - TB->hS[Ua]) fR12++;
                    if (deg[a] == 0 && sg[a] != 0) fa0++;
                    if (deg[a] == 1) { int w = __builtin_ctz(B[a]); ll ew = TB->end[Ua][w], Rw = Rsub(TB, Ua, w); if (sg[a] != ew + Rw) fa1++; if (lmina1 < 0 || sg[a] < lmina1) lmina1 = sg[a]; }
                }
            }
            #pragma omp critical
            {
                if (lminold >= 0 && (minold < 0 || lminold < minold)) minold = lminold;
                if (lminns >= 0 && (minns < 0 || lminns < minns)) minns = lminns;
                if (lmina1 >= 0 && (mina1 < 0 || lmina1 < mina1)) mina1 = lmina1;
            }
            free(TB); free(T1);
        }
        printf("order %d: classes %zu (A000568: %ld) %s | (B,a) %ld | failures: D %ld, D+ %ld, native %ld, sink formula %ld (of %ld), min4 %ld, digit-4 iff (TT,penult,sink) %ld, non-sink>=9 %ld, RESULT(12) %ld, mark-sink %ld, mark d=1 %ld | min old %lld, min non-sink old %lld (count of 9s %ld; 9 iff (TT, a sink, w penultimate) failures %ld), digit-4 cases %ld, min mark d=1 %lld\n",
               n, ncur, A000568[n], (long)ncur == A000568[n] ? "OK" : "MISMATCH", pairs, fD, fDp, fnat, fsink, sinkc, fmin4, feq, fns9, fR12, fa0, fa1, minold, minns, n9, n9other, n4, mina1);
        fflush(stdout);
    }
    return 0;
}
