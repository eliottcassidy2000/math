/* audit_census.c  --  independent audit program (not derived from the author's code).
 *
 * Reads isomorphism classes of tournaments of order n from `gentourng -q n` (upper triangle, row by row, in ascii).
 * Convention used here: the character for the pair (i,j), i<j, is '1' iff i -> j.  (Every statistic printed is either
 * invariant under taking the converse or is checked for the tournament AND its converse, and the list of classes is
 * closed under converse, so the opposite reading of the file gives the same tables.)
 *
 * For each class:
 *   H            number of Hamiltonian paths, subset DP (exact)
 *   hc[S]        number of directed Hamiltonian cycles of T[S] for every vertex set S, DP from the least vertex (exact)
 *   c_k          = sum_{|S|=k} hc[S];  alpha_1 = sum over odd |S|;  alpha_2, alpha_3 = numbers of pairs / triples of
 *                vertex-disjoint odd cycles (from precomputed lists of disjoint odd sets)
 *   check        H == 1 + 2 alpha_1 + 4 alpha_2 + 8 alpha_3   (exact; valid for n <= 11)
 *   check        c_3 = C(q1,2) + floor(n/4), c_3 + c_5 = C(q2,2) + floor(n/4) (mod 2) for T and for its converse
 *   |Aut|        backtracking;  self-converse: backtracking search for an anti-isomorphism;  strong: reachability
 *   option -b    (n <= 8) also brute force: all n! permutations for |Aut| and self-converse, DFS for cycles and paths
 * Output: one line per non-empty bin
 *   sc sym strong Hmod8 p3 p5 p7 p9  classes  labelled_weight        (weight = n!/|Aut|)
 * and a final line with totals and failure counts.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXN 11
static int n;
static unsigned outm[MAXN], inm[MAXN];
static unsigned char adj[MAXN][MAXN];
static int score[MAXN];

static int popc(unsigned x) { return __builtin_popcount(x); }
static int C2(int k) { return k * (k - 1) / 2; }

/* ---------- H by subset DP ---------- */
static uint32_t dpH[1 << MAXN][MAXN];
static uint64_t ham_paths(void) {
    unsigned full = (1u << n) - 1;
    for (unsigned mask = 1; mask <= full; mask++) for (int v = 0; v < n; v++) dpH[mask][v] = 0;
    for (int v = 0; v < n; v++) dpH[1u << v][v] = 1;
    for (unsigned mask = 1; mask <= full; mask++) {
        unsigned rem = mask;
        while (rem) {
            int last = __builtin_ctz(rem); rem &= rem - 1;
            uint32_t c = dpH[mask][last]; if (!c) continue;
            unsigned nx = outm[last] & ~mask;
            while (nx) { int w = __builtin_ctz(nx); nx &= nx - 1; dpH[mask | (1u << w)][w] += c; }
        }
    }
    uint64_t H = 0; for (int v = 0; v < n; v++) H += dpH[full][v];
    return H;
}

/* ---------- Hamiltonian cycles of every induced subtournament ---------- */
static uint32_t hc[1 << MAXN];
static uint32_t dpC[1 << (MAXN - 1)][MAXN];
static void ham_cycles(void) {
    memset(hc, 0, sizeof(uint32_t) << n);
    for (int v = 0; v < n - 2; v++) {
        int hi = n - 1 - v;                       /* vertices v+1 .. n-1, re-indexed 0 .. hi-1 */
        unsigned nsub = 1u << hi;
        for (unsigned sub = 0; sub < nsub; sub++) for (int j = 0; j < hi; j++) dpC[sub][j] = 0;
        for (int j = 0; j < hi; j++) if (adj[v][v + 1 + j]) dpC[1u << j][j] = 1;
        for (unsigned sub = 1; sub < nsub; sub++) {
            unsigned rem = sub;
            while (rem) {
                int j = __builtin_ctz(rem); rem &= rem - 1;
                uint32_t c = dpC[sub][j]; if (!c) continue;
                int last = v + 1 + j;
                if (adj[last][v] && popc(sub) >= 2) hc[(sub << (v + 1)) | (1u << v)] += c;
                unsigned nx = (outm[last] >> (v + 1)) & ~sub;
                while (nx) { int k = __builtin_ctz(nx); nx &= nx - 1; dpC[sub | (1u << k)][k] += c; }
            }
        }
    }
}

/* ---------- lists of disjoint odd sets ---------- */
static uint32_t *pairA, *pairB; static long npair;
static uint32_t *triA, *triB, *triC; static long ntri;
static void build_lists(void) {
    unsigned full = (1u << n) - 1; long cap = 1 << 16;
    pairA = malloc(cap * 4); pairB = malloc(cap * 4); npair = 0;
    long capt = 1 << 16; triA = malloc(capt * 4); triB = malloc(capt * 4); triC = malloc(capt * 4); ntri = 0;
    for (unsigned S1 = 1; S1 <= full; S1++) {
        int k1 = popc(S1); if (k1 < 3 || !(k1 & 1)) continue;
        unsigned comp = full & ~S1;
        for (unsigned S2 = comp; S2; S2 = (S2 - 1) & comp) {
            int k2 = popc(S2); if (k2 < 3 || !(k2 & 1) || S2 < S1) continue;
            if (npair == cap) { cap *= 2; pairA = realloc(pairA, cap * 4); pairB = realloc(pairB, cap * 4); }
            pairA[npair] = S1; pairB[npair] = S2; npair++;
            unsigned comp2 = comp & ~S2;
            for (unsigned S3 = comp2; S3; S3 = (S3 - 1) & comp2) {
                int k3 = popc(S3); if (k3 < 3 || !(k3 & 1) || S3 < S2) continue;
                if (ntri == capt) { capt *= 2; triA = realloc(triA, capt * 4); triB = realloc(triB, capt * 4); triC = realloc(triC, capt * 4); }
                triA[ntri] = S1; triB[ntri] = S2; triC[ntri] = S3; ntri++;
            }
        }
    }
}

/* ---------- automorphisms and anti-isomorphisms by backtracking ---------- */
static int img[MAXN]; static unsigned used; static long autcount;
static void aut_bt(int v) {
    if (v == n) { autcount++; return; }
    for (int w = 0; w < n; w++) {
        if ((used >> w) & 1) continue;
        if (score[w] != score[v]) continue;
        int ok = 1;
        for (int u = 0; u < v; u++) if (adj[u][v] != adj[img[u]][w]) { ok = 0; break; }
        if (!ok) continue;
        img[v] = w; used |= 1u << w; aut_bt(v + 1); used &= ~(1u << w);
    }
}
static int anti_bt(int v) {
    if (v == n) return 1;
    for (int w = 0; w < n; w++) {
        if ((used >> w) & 1) continue;
        if (score[w] != n - 1 - score[v]) continue;
        int ok = 1;
        for (int u = 0; u < v; u++) if (adj[u][v] != adj[w][img[u]]) { ok = 0; break; }
        if (!ok) continue;
        img[v] = w; used |= 1u << w;
        if (anti_bt(v + 1)) { used &= ~(1u << w); return 1; }
        used &= ~(1u << w);
    }
    return 0;
}
static int is_strong(void) {
    unsigned full = (1u << n) - 1, r = 1, f = 1;
    while (f) { int v = __builtin_ctz(f); f &= f - 1; unsigned nw = outm[v] & ~r; r |= nw; f |= nw; }
    if (r != full) return 0;
    r = 1; f = 1;
    while (f) { int v = __builtin_ctz(f); f &= f - 1; unsigned nw = inm[v] & ~r; r |= nw; f |= nw; }
    return r == full;
}

/* ---------- brute force (option -b) ---------- */
static long bcyc[MAXN + 1], bH;
static void dfs_cycle(int start, int last, int len, unsigned vis) {
    if (len >= 3 && adj[last][start]) bcyc[len]++;
    for (int w = start + 1; w < n; w++) if (!((vis >> w) & 1) && adj[last][w]) dfs_cycle(start, w, len + 1, vis | (1u << w));
}
static void dfs_ham(int last, unsigned vis, int cnt) {
    if (cnt == n) { bH++; return; }
    for (int w = 0; w < n; w++) if (!((vis >> w) & 1) && adj[last][w]) dfs_ham(w, vis | (1u << w), cnt + 1);
}
static void brute_sym(long *naut, int *sc) {
    int p[MAXN]; for (int i = 0; i < n; i++) p[i] = i;
    *naut = 0; *sc = 0;
    for (;;) {
        int ok = 1;
        for (int j = 1; j < n && ok; j++) for (int i = 0; i < j; i++) if (adj[i][j] != adj[p[i]][p[j]]) { ok = 0; break; }
        if (ok) (*naut)++;
        if (!*sc) {
            ok = 1;
            for (int j = 1; j < n && ok; j++) for (int i = 0; i < j; i++) if (adj[i][j] != adj[p[j]][p[i]]) { ok = 0; break; }
            if (ok) *sc = 1;
        }
        int i = n - 2; while (i >= 0 && p[i] > p[i + 1]) i--;
        if (i < 0) break;
        int j = n - 1; while (p[j] < p[i]) j--;
        int t = p[i]; p[i] = p[j]; p[j] = t;
        for (int a = i + 1, b = n - 1; a < b; a++, b--) { t = p[a]; p[a] = p[b]; p[b] = t; }
    }
}

/* option -g: exact cube sums of the base over all 2^n out-sets X of a new vertex, by direct evaluation of
 *   inc3(X) = number of arcs from X to V\X                 (increment of c_3, exact)
 *   inc5(X) = number of walks of length 3 from X to V\X    (increment of c_5 mod 2; a 4-vertex walk between distinct ends is a path)
 * giving the exact labelled sums of (-1)^(c_3), (-1)^(c_5), (-1)^(c_3+c_5) at order n+1. */
static void cube35(int p3, int p5, long *g3, long *g5, long *g35) {
    unsigned A3[MAXN], A2[MAXN], full = (1u << n) - 1;
    for (int a = 0; a < n; a++) { unsigned r = 0, o = outm[a]; while (o) { int k = __builtin_ctz(o); o &= o - 1; r ^= outm[k]; } A2[a] = r; }
    for (int a = 0; a < n; a++) { unsigned r = 0, o = outm[a]; while (o) { int k = __builtin_ctz(o); o &= o - 1; r ^= A2[k]; } A3[a] = r; }
    long s3 = 0, s5 = 0, s35 = 0;
    for (unsigned X = 0; X <= full; X++) {
        unsigned nX = full & ~X, rem = X; int e3 = 0, e5 = 0;
        while (rem) { int a = __builtin_ctz(rem); rem &= rem - 1; e3 ^= popc(outm[a] & nX) & 1; e5 ^= popc(A3[a] & nX) & 1; }
        s3 += e3 ? -1 : 1; s5 += e5 ? -1 : 1; s35 += (e3 ^ e5) ? -1 : 1;
    }
    *g3 = p3 ? -s3 : s3; *g5 = p5 ? -s5 : s5; *g35 = (p3 ^ p5) ? -s35 : s35;
}

int main(int argc, char **argv) {
    int brute = 0, ai = 1, docube = 0;
    __int128 G3 = 0, G5 = 0, G35 = 0;
    if (argc > 1 && strcmp(argv[1], "-b") == 0) { brute = 1; ai = 2; }
    if (argc > 1 && strcmp(argv[1], "-g") == 0) { docube = 1; ai = 2; }
    if (argc <= ai) { fprintf(stderr, "usage: audit_census [-b] n < classes\n"); return 1; }
    n = atoi(argv[ai]);
    if (n < 1 || n > MAXN) return 1;
    int E = C2(n);
    uint64_t nfact = 1; for (int i = 2; i <= n; i++) nfact *= i;
    build_lists();
    static uint64_t binC[1024], binW[1024];
    uint64_t classes = 0, totW = 0;
    long f_ocf = 0, f_p1a = 0, f_p1b = 0, f_brute = 0, f_mod4 = 0;
    char line[256];
    while (fgets(line, sizeof line, stdin)) {
        int len = 0; while (line[len] == '0' || line[len] == '1') len++;
        if (len != E) { if (len == 0 && n == 1) { /* empty line for n = 1 */ } else continue; }
        int pos = 0;
        for (int i = 0; i < n; i++) { outm[i] = inm[i] = 0; for (int j = 0; j < n; j++) adj[i][j] = 0; }
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
            int b = line[pos++] == '1';
            if (b) { adj[i][j] = 1; outm[i] |= 1u << j; inm[j] |= 1u << i; }
            else   { adj[j][i] = 1; outm[j] |= 1u << i; inm[i] |= 1u << j; }
        }
        for (int i = 0; i < n; i++) score[i] = popc(outm[i]);
        uint64_t H = ham_paths();
        ham_cycles();
        uint64_t ck[MAXN + 1]; memset(ck, 0, sizeof ck);
        uint64_t a1 = 0, a2 = 0, a3 = 0;
        for (unsigned S = 0; S < (1u << n); S++) if (hc[S]) { int k = popc(S); ck[k] += hc[S]; if (k & 1) a1 += hc[S]; }
        for (long i = 0; i < npair; i++) a2 += (uint64_t)hc[pairA[i]] * hc[pairB[i]];
        for (long i = 0; i < ntri; i++) a3 += (uint64_t)hc[triA[i]] * hc[triB[i]] * hc[triC[i]];
        if (H != 1 + 2 * a1 + 4 * a2 + 8 * a3) f_ocf++;
        if ((H & 3) != ((1 + 2 * a1) & 3)) f_mod4++;
        /* score formulas for T and its converse */
        for (int conv = 0; conv < 2; conv++) {
            unsigned s = 0, t = 0;
            for (int a = 0; a < n; a++) if (popc(conv ? inm[a] : outm[a]) & 1) s |= 1u << a;
            for (int a = 0; a < n; a++) if (popc((conv ? inm[a] : outm[a]) & s) & 1) t |= 1u << a;
            int q1 = popc(s), q2 = popc(t);
            if (((ck[3] + C2(q1) + n / 4) & 1) != 0) f_p1a++;
            if (((ck[3] + ck[5] + C2(q2) + n / 4) & 1) != 0) f_p1b++;
        }
        used = 0; autcount = 0; aut_bt(0);
        used = 0; int sc = anti_bt(0);
        int strong = is_strong();
        if (brute) {
            long na; int bsc; brute_sym(&na, &bsc);
            memset(bcyc, 0, sizeof bcyc); for (int v = 0; v < n; v++) dfs_cycle(v, v, 1, 1u << v);
            bH = 0; for (int v = 0; v < n; v++) dfs_ham(v, 1u << v, 1);
            int bad = (na != autcount) || (bsc != sc) || ((uint64_t)bH != H);
            for (int k = 3; k <= n; k++) if ((uint64_t)bcyc[k] != ck[k]) bad = 1;
            if (bad) f_brute++;
        }
        uint64_t w = nfact / (uint64_t)autcount;
        if (nfact % (uint64_t)autcount) f_brute += 1000000;
        int key = (sc << 8) | ((autcount > 1) << 7) | (strong << 6) | ((int)((H & 7) >> 1) << 4)
                | ((int)(ck[3] & 1) << 3) | ((int)(ck[5] & 1) << 2) | ((int)(ck[7] & 1) << 1) | (int)(ck[9] & 1);
        binC[key]++; binW[key] += w; classes++; totW += w;
        if (docube) {
            long g3, g5, g35; cube35((int)(ck[3] & 1), (int)(ck[5] & 1), &g3, &g5, &g35);
            G3 += (__int128)w * g3; G5 += (__int128)w * g5; G35 += (__int128)w * g35;
        }
    }
    if (docube)
        fprintf(stderr, "order %d from the bases of order %d: exact labelled sums over 2^%d: sum(-1)^c3 = %lld ; sum(-1)^c5 = %lld ; sum(-1)^(c3+c5) = %lld\n",
                n + 1, n, E + n, (long long)G3, (long long)G5, (long long)G35);
    for (int key = 0; key < 512; key++) if (binC[key])
        printf("%d %d %d %d %d %d %d %d %llu %llu\n", (key >> 8) & 1, (key >> 7) & 1, (key >> 6) & 1, 2 * ((key >> 4) & 3) + 1,
               (key >> 3) & 1, (key >> 2) & 1, (key >> 1) & 1, key & 1, (unsigned long long)binC[key], (unsigned long long)binW[key]);
    fprintf(stderr, "n %d classes %llu total_weight %llu (2^%d = %llu) failures: OCF %ld mod4 %ld P1_c3 %ld P1_c3c5 %ld brute %ld\n",
            n, (unsigned long long)classes, (unsigned long long)totW, E, (unsigned long long)1 << E, f_ocf, f_mod4, f_p1a, f_p1b, f_brute);
    return 0;
}
