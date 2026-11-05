/* aud_at_family.c -- independent audit code (uh_deletion_law_20261104_audit_S10). Shares no code with the session.
 *
 * Tests THM-4874 (5) on explicitly built one-arc joins T = Q[z <- (D, u)] with the block (AT_m, {m-1, m}), m = 4..MMAX:
 *   D = chain w_1 => ... => w_p => AT_m, plus u; u beaten by A = {m-1, m}, u beats the rest of AT_m, u -> w_1, the arcs
 *   u -- w_2..w_p in every orientation; p = 0..PMAX;
 *   prefixes: the one-vertex prefix {l}, and every rooted uniquely Hamiltonian (Q, z) with z a sole exit of order 3..QMAX,
 *   up to isomorphism (found by brute force; the counts must be R(3..5) = 1, 2, 5 = F(1), F(3), F(5)).
 * For each T it computes hc(T) and every hc(T - v) by subset DP and checks:
 *   hc(T) = 2^(m-3) + 1; T has no vertex of out-degree 1;
 *   exactly three v with hc(T - v) = 1, namely the block vertices 1, m-1, m;
 *   hc(T - j) = 2^(m-4) + 1 for the block vertices 2..m-2;
 *   hc(T - l) = h0 = 2^(m-3) when R = {l} and p = 0 (k when R = {l}, p >= 1; 0 when |R| >= 2);
 *   chain vertices k, u 0, prefix vertices multiples of k.
 * Usage: aud_at_family [MMAX [PMAX [QMAX]]] */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef unsigned long long u64;
#define MAXN 20

static int popc(unsigned x) { return __builtin_popcount(x); }
static u64 *DP;

/* Hamiltonian cycles of the tournament on vertex set 'mask' of out[] (anchor = lowest vertex of mask) */
static u64 hc_mask(int n, const unsigned *out, unsigned mask) {
    int k = popc(mask); if (k < 3) return 0;
    int idx[MAXN], vs[MAXN], c = 0; for (int v = 0; v < n; v++) if ((mask >> v) & 1) { idx[v] = c; vs[c++] = v; }
    unsigned o[MAXN]; for (int i = 0; i < k; i++) { o[i] = 0; for (int j = 0; j < k; j++) if ((out[vs[i]] >> vs[j]) & 1) o[i] |= 1u << j; }
    unsigned full = (1u << k) - 1;
    memset(DP, 0, sizeof(u64) * ((size_t)full + 1) * k);
    DP[1 * k + 0] = 1;
    for (unsigned S = 1; S <= full; S += 2) for (int v = 0; v < k; v++) {
        u64 x = DP[(size_t)S * k + v]; if (!x) continue;
        unsigned nb = o[v] & ~S & full;
        while (nb) { int w = __builtin_ctz(nb); nb &= nb - 1; DP[(size_t)(S | (1u << w)) * k + w] += x; }
    }
    u64 h = 0; for (int v = 1; v < k; v++) if (o[v] & 1u) h += DP[(size_t)full * k + v];
    (void)idx;
    return h;
}
/* Hamiltonian paths of the sub-tournament on 'mask' starting at s */
static u64 paths_from(int n, const unsigned *out, unsigned mask, int s) {
    int k = popc(mask); if (k == 1) return 1;
    int vs[MAXN], c = 0, si = -1; for (int v = 0; v < n; v++) if ((mask >> v) & 1) { if (v == s) si = c; vs[c++] = v; }
    unsigned o[MAXN]; for (int i = 0; i < k; i++) { o[i] = 0; for (int j = 0; j < k; j++) if ((out[vs[i]] >> vs[j]) & 1) o[i] |= 1u << j; }
    unsigned full = (1u << k) - 1;
    memset(DP, 0, sizeof(u64) * ((size_t)full + 1) * k);
    DP[(size_t)(1u << si) * k + si] = 1;
    for (unsigned S = 1; S <= full; S++) { if (!((S >> si) & 1)) continue; for (int v = 0; v < k; v++) {
        u64 x = DP[(size_t)S * k + v]; if (!x) continue;
        unsigned nb = o[v] & ~S & full;
        while (nb) { int w = __builtin_ctz(nb); nb &= nb - 1; DP[(size_t)(S | (1u << w)) * k + w] += x; } } }
    u64 h = 0; for (int v = 0; v < k; v++) h += DP[(size_t)full * k + v];
    return h;
}

/* rooted uniquely Hamiltonian prefixes (Q, z), z a sole exit, up to isomorphism, by brute force */
typedef struct { int q; unsigned out[8]; int z; } Pref;
static Pref prefs[64]; static int npref = 0;
static int perm_next(int *p, int n) {
    int i = n - 2; while (i >= 0 && p[i] > p[i + 1]) i--; if (i < 0) return 0;
    int j = n - 1; while (p[j] < p[i]) j--; int t = p[i]; p[i] = p[j]; p[j] = t;
    for (int a = i + 1, b = n - 1; a < b; a++, b--) { t = p[a]; p[a] = p[b]; p[b] = t; }
    return 1;
}
static u64 rooted_key(int q, const unsigned *out, int z) {
    u64 best = ~0ull; int p[8]; for (int i = 0; i < q; i++) p[i] = i;
    do {
        if (p[z] != 0) continue;   /* root gets label 0 */
        unsigned nr[8] = {0}; for (int i = 0; i < q; i++) for (int j = 0; j < q; j++) if ((out[i] >> j) & 1) nr[p[i]] |= 1u << p[j];
        u64 code = 0; for (int i = 0; i < q; i++) code = (code << q) | nr[i];
        if (code < best) best = code;
    } while (perm_next(p, q));
    return best;
}
static void find_prefixes(int qmax) {
    u64 keys[64]; int nk = 0;
    for (int q = 3; q <= qmax; q++) {
        int np = q * (q - 1) / 2, pi[16], pj[16], c = 0, before = npref;
        for (int i = 0; i < q; i++) for (int j = i + 1; j < q; j++) { pi[c] = i; pj[c] = j; c++; }
        for (unsigned bits = 0; bits < (1u << np); bits++) {
            unsigned out[8] = {0};
            for (int b = 0; b < np; b++) { if ((bits >> b) & 1) out[pi[b]] |= 1u << pj[b]; else out[pj[b]] |= 1u << pi[b]; }
            if (hc_mask(q, out, (1u << q) - 1) != 1) continue;
            for (int z = 0; z < q; z++) {
                if (popc(out[z]) != 1) continue;
                u64 key = rooted_key(q, out, z) * 16 + (u64)q;
                int f = 0; for (int i = 0; i < nk; i++) if (keys[i] == key) { f = 1; break; }
                if (f) continue;
                keys[nk++] = key;
                prefs[npref].q = q; memcpy(prefs[npref].out, out, sizeof out); prefs[npref].z = z; npref++;
            }
        }
        printf("rooted uniquely Hamiltonian prefixes (Q, z) of order %d: %d\n", q, npref - before);
    }
}

int main(int argc, char **argv) {
    int MMAX = argc > 1 ? atoi(argv[1]) : 9, PMAX = argc > 2 ? atoi(argv[2]) : 2, QMAX = argc > 3 ? atoi(argv[3]) : 5;
    DP = malloc(sizeof(u64) * ((size_t)1 << MAXN) * 1);   /* resized below */
    free(DP); DP = malloc(sizeof(u64) * ((size_t)1 << 17) * 17);
    find_prefixes(QMAX);
    long nT = 0, nbad = 0, ndel = 0;
    for (int m = 4; m <= MMAX; m++) {
        u64 k = (1ull << (m - 3)) + 1, inner = (1ull << (m - 4)) + 1, h0 = 1ull << (m - 3);
        long nTm = 0;
        for (int p = 0; p <= PMAX; p++) {
            for (unsigned fa = 0; fa < (p >= 2 ? (1u << (p - 1)) : 1u); fa++) {
                for (int pr = -1; pr < npref; pr++) {
                    /* vertices: R (l = 0, then the rest of Q - z), chain, W = AT_m (block vertex j = 1..m), u */
                    int r = (pr < 0) ? 1 : prefs[pr].q - 1;
                    int n = r + p + m + 1;
                    if (n > 17) continue;
                    unsigned out[MAXN] = {0};
                    int l = 0;
                    if (pr >= 0) {
                        const Pref *P = &prefs[pr]; int zz = P->z; int ll = __builtin_ctz(P->out[zz]);
                        int lab[8], c = 0; lab[ll] = c++; for (int x = 0; x < P->q; x++) if (x != zz && x != ll) lab[x] = c++;
                        for (int x = 0; x < P->q; x++) if (x != zz) for (int y = 0; y < P->q; y++) if (y != zz && ((P->out[x] >> y) & 1)) out[lab[x]] |= 1u << lab[y];
                    }
                    int ch0 = r, W0 = r + p, u = n - 1;       /* chain w_i = ch0 + i - 1; block vertex j = W0 + j - 1 */
                    /* R beats all of D except u -> l */
                    for (int i = 0; i < r; i++) for (int y = ch0; y < n; y++) if (!(i == l && y == u)) out[i] |= 1u << y;
                    out[u] |= 1u << l;
                    for (int a = 0; a < p; a++) { for (int b = a + 1; b < p; b++) out[ch0 + a] |= 1u << (ch0 + b); for (int j = 0; j < m; j++) out[ch0 + a] |= 1u << (W0 + j); }
                    for (int i = 1; i <= m; i++) for (int j = i + 1; j <= m; j++) {
                        if (i == 1 && j == m) out[W0 + m - 1] |= 1u << W0;          /* m -> 1 */
                        else out[W0 + i - 1] |= 1u << (W0 + j - 1);
                    }
                    for (int j = 1; j <= m; j++) { if (j >= m - 1) out[W0 + j - 1] |= 1u << u; else out[u] |= 1u << (W0 + j - 1); }
                    if (p >= 1) { out[u] |= 1u << ch0; for (int t = 1; t < p; t++) { if ((fa >> (t - 1)) & 1) out[u] |= 1u << (ch0 + t); else out[ch0 + t] |= 1u << u; } }
                    /* sanity: tournament */
                    for (int x = 0; x < n; x++) for (int y = x + 1; y < n; y++) if ((((out[x] >> y) & 1) + ((out[y] >> x) & 1)) != 1) { printf("not a tournament\n"); return 1; }
                    unsigned full = (1u << n) - 1;
                    int mo = 99; for (int x = 0; x < n; x++) if (popc(out[x]) < mo) mo = popc(out[x]);
                    u64 h = hc_mask(n, out, full);
                    int bad = (h != k) || (mo < 2);
                    int k1 = 0, ones_ok = 1;
                    char line[512]; int pos = 0;
                    for (int v = 0; v < n; v++) {
                        u64 hv = hc_mask(n, out, full & ~(1u << v)); ndel++;
                        u64 pred; int isblock = (v >= W0 && v < W0 + m);
                        if (v < r && v != l) { pred = k * paths_from(n, out, ((1u << r) - 1) & ~(1u << v), l); if (hv % k) bad = 1; }
                        else if (v == u) pred = 0;
                        else if (v == l) pred = (r >= 2) ? 0 : (p >= 1 ? k : h0);
                        else if (v >= ch0 && v < W0) pred = k;
                        else { int j = v - W0 + 1; pred = (j == 1 || j >= m - 1) ? 1 : inner; }
                        if (hv != pred) { bad = 1; pos += snprintf(line + pos, sizeof line - pos, " [v=%d got %llu pred %llu]", v, hv, pred); }
                        if (hv == 1) { k1++; if (!isblock) ones_ok = 0; else { int j = v - W0 + 1; if (!(j == 1 || j >= m - 1)) ones_ok = 0; } }
                    }
                    if (k1 != 3 || !ones_ok) bad = 1;
                    nT++; nTm++;
                    if (bad) { nbad++; printf("  FAIL m=%d p=%d fa=%u prefix=%d n=%d hc=%llu minout=%d k1=%d%s\n", m, p, fa, pr, n, h, mo, k1, line); }
                }
            }
        }
        printf("m=%d: type k = 2^(m-3)+1 = %llu, interior value 2^(m-4)+1 = %llu, h0 = 2^(m-3) = %llu; one-arc joins built and checked: %ld\n",
               m, k, inner, h0, nTm);
        fflush(stdout);
    }
    printf("total: %ld tournaments, %ld vertex deletions, %ld failures\n", nT, ndel, nbad);
    return 0;
}
