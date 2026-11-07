/* Independent audit code (auditor's own implementation) for THM-4890 (C), (C'), (E) and HYP-9482 statistics.

   g(0)=1, g(2v+1)=g(v), g(2v)=g(v)+v+(v mod 2), tabulated for v < 2^n (uint32).
   Y_n = even v in [2^(n-1), 2^n).  Values are grouped by a 64-bit LSD radix sort of keys (g << 32 | index).

   For every ordered collision pair (p, q), p != q, in Y_n (converse pairs included):
     t, s = highest / lowest differing bit; a = n-1-t; w = t-s+1; x, x' = window values;
     X = x - x', Xi = C_w(x) - C_w(x'); dX = w - bitlen|X|, dXi = w - bitlen|Xi|; A = a - dXi, S = s - dX.
   Checked: a, s >= 1; X, Xi odd; the bounds 2^(A+w-2) <= h|Xi| < 2^(A+w), 2^(S+w-2) <= tau|X| < 2^(S+w);
   |h Xi + tau X| <= 2^w + 1; same-sign case => max(A,S) <= 2;
   (C)  |A-S| <= 2 or max(A,S) <= 3;   sharp: |A-S| <= 1 or max(A,S) <= 3;
   P5 (non-converse pairs): |a-s| <= 1 + max(dX, dXi);
   (C') (p*, q*) is a collision pair in Y_n with (a,s,dX,dXi) -> (s,a,dXi,dX).
   Statistics (HYP-9482): E^Y/|Y|, converse / non-converse parts, values/|Y|, multiplicity histogram,
   smoothed slopes (E(K)-E(0))/(K|Y|) for K = 2,4,8,16, level law L_a, max multiplicity;
   X_n energy E_n computed directly on all of X_n (independent of the layer decomposition), the layer sum
   S_n = sum_(j=2..n) E^Y_j, the cross-layer part (E_n - S_n - 1)/|X_n|, (E) range and overlap checks.
   Usage: aud_pairs n_min n_max
*/
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

static uint32_t *G = NULL;
static int GN = 0;

static void build_g(int n) {
    if (GN >= n) return;
    free(G);
    G = (uint32_t *)malloc(sizeof(uint32_t) << n);
    if (!G) { fprintf(stderr, "oom g\n"); exit(1); }
    G[0] = 1;
    for (uint64_t v = 1; v < (1ull << n); v++) {
        uint64_t u = v >> 1;
        G[v] = (v & 1) ? G[u] : (uint32_t)(G[u] + u + (u & 1));
    }
    GN = n;
}

static inline int bitlen64(uint64_t x) { return x ? 64 - __builtin_clzll(x) : 0; }

static inline uint32_t rev_w(uint32_t x, int w) {   /* reverse the w low bits */
    uint32_t r = 0;
    for (int i = 0; i < w; i++) { r = (r << 1) | (x & 1); x >>= 1; }
    return r;
}
static inline uint32_t Cw(uint32_t x, int w) { return ((w >= 32) ? 0xffffffffu : ((1u << w) - 1)) - rev_w(x, w); }

static void radix64(uint64_t *a, uint64_t *tmp, size_t m) {
    size_t cnt[65536];
    for (int pass = 0; pass < 4; pass++) {
        int sh = 16 * pass;
        memset(cnt, 0, sizeof(cnt));
        for (size_t i = 0; i < m; i++) cnt[(a[i] >> sh) & 0xffff]++;
        size_t s = 0;
        for (int k = 0; k < 65536; k++) { size_t c = cnt[k]; cnt[k] = s; s += c; }
        for (size_t i = 0; i < m; i++) tmp[cnt[(a[i] >> sh) & 0xffff]++] = a[i];
        uint64_t *t = a; a = tmp; tmp = t;   /* after 4 passes the data is back in the original array */
    }
}

typedef struct { int a, s, w, dX, dXi; long long X, Xi; uint32_t x, xp; } pinfo;

static void pair_info(uint32_t p, uint32_t q, int n, pinfo *I) {
    uint32_t d = p ^ q;
    int t = 31 - __builtin_clz(d);
    int s = __builtin_ctz(d);
    int w = t - s + 1;
    uint32_t mask = (w >= 32) ? 0xffffffffu : ((1u << w) - 1);
    uint32_t x = (p >> s) & mask, xp = (q >> s) & mask;
    long long X = (long long)x - (long long)xp;
    long long Xi = (long long)Cw(x, w) - (long long)Cw(xp, w);
    I->a = n - 1 - t; I->s = s; I->w = w; I->x = x; I->xp = xp; I->X = X; I->Xi = Xi;
    I->dX = w - bitlen64((uint64_t)(X < 0 ? -X : X));
    I->dXi = w - bitlen64((uint64_t)(Xi < 0 ? -Xi : Xi));
}

/* energy of Y_j only (for the layer sum) */
static long long energyY(int j, uint64_t *keys, uint64_t *tmp) {
    size_t m = (size_t)1 << (j - 2);
    uint32_t lo = 1u << (j - 1);
    for (size_t i = 0; i < m; i++) keys[i] = ((uint64_t)G[lo + 2 * i] << 32) | i;
    radix64(keys, tmp, m);
    long long E = 0; size_t i = 0;
    while (i < m) { size_t k = i + 1; while (k < m && (keys[k] >> 32) == (keys[i] >> 32)) k++; E += (long long)(k - i) * (k - i); i = k; }
    return E;
}

int main(int argc, char **argv) {
    int n0 = atoi(argv[1]), n1 = atoi(argv[2]);
    build_g(n1);
    size_t mmax = (size_t)1 << (n1 - 1);
    uint64_t *keys = (uint64_t *)malloc(sizeof(uint64_t) * mmax);
    uint64_t *tmp = (uint64_t *)malloc(sizeof(uint64_t) * mmax);
    if (!keys || !tmp) { fprintf(stderr, "oom\n"); return 1; }
    long long EYj[64]; memset(EYj, 0, sizeof(EYj));
    for (int j = 2; j < n0 && j <= n1; j++) EYj[j] = energyY(j, keys, tmp);
    for (int n = n0; n <= n1; n++) {
        size_t m = (size_t)1 << (n - 2);
        uint32_t lo = 1u << (n - 1);
        for (size_t i = 0; i < m; i++) keys[i] = ((uint64_t)G[lo + 2 * i] << 32) | i;
        radix64(keys, tmp, m);
        long long hist[256]; memset(hist, 0, sizeof(hist));
        long long E = 0, conv = 0, nonconv = 0, nvals = 0;
        long long viol_C = 0, viol_sharp = 0, viol_P5 = 0, viol_Cp = 0, viol_misc = 0, same_sign = 0, same_sign_big = 0;
        long long viol_sessionproved = 0;
        long long lvl[64]; memset(lvl, 0, sizeof(lvl));
        long long dist[81][16]; memset(dist, 0, sizeof(dist));   /* (A-S+40, min(max(A,S)+4,15)) over all pairs */
        long long absAS2 = 0, absAS_ge3 = 0, absAS_ge3_maxle3 = 0;
        int maxAS_when_big_gap = -99; int ex_printed = 0;
        size_t i = 0;
        while (i < m) {
            size_t k = i + 1;
            while (k < m && (keys[k] >> 32) == (keys[i] >> 32)) k++;
            long long r = (long long)(k - i);
            nvals++;
            hist[r < 255 ? r : 255]++;
            E += r * r;
            for (size_t ii = i; ii < k; ii++) for (size_t jj = i; jj < k; jj++) {
                if (ii == jj) continue;
                uint32_t p = lo + 2 * (uint32_t)(keys[ii] & 0xffffffffu);
                uint32_t q = lo + 2 * (uint32_t)(keys[jj] & 0xffffffffu);
                uint32_t pstar = Cw(p, n), qstar = Cw(q, n);
                int isconv = (pstar == q);
                if (isconv) conv++; else nonconv++;
                pinfo I; pair_info(p, q, n, &I);
                int a = I.a, s = I.s, w = I.w;
                if (a < 1 || s < 1 || !(I.X & 1) || !(I.Xi & 1)) viol_misc++;
                int A = a - I.dXi, S = s - I.dX;
                /* proof inequalities */
                uint64_t h = p >> (I.s + w);
                uint32_t sig = p & ((1u << s) - 1);
                uint64_t tau = Cw(sig, s);
                unsigned long long hXi = h * (unsigned long long)(I.Xi < 0 ? -I.Xi : I.Xi);
                unsigned long long tX = tau * (unsigned long long)(I.X < 0 ? -I.X : I.X);
                /* 2^(A+w-2) <= h|Xi| < 2^(A+w): compare via exponents (A+w-2 may be negative only if ... check) */
                int eA = A + w, eS = S + w;
                if (eA - 2 < 0 || eS - 2 < 0) viol_misc++;
                else {
                    if (!(hXi >= (1ull << (eA - 2)) && hXi < (1ull << eA))) viol_misc++;
                    if (!(tX >= (1ull << (eS - 2)) && tX < (1ull << eS))) viol_misc++;
                }
                long long lhs = (long long)h * I.Xi + (long long)tau * I.X;
                if ((lhs < 0 ? -lhs : lhs) > (1ll << w) + 1) viol_misc++;
                if ((I.X > 0) == (I.Xi > 0)) { same_sign++; if (A > 2 || S > 2) same_sign_big++; }
                int dAS = A - S; int adAS = dAS < 0 ? -dAS : dAS; int mx = A > S ? A : S;
                if (!(adAS <= 2 || mx <= 3)) viol_C++;
                if (!(adAS <= 1 || mx <= 3)) viol_sharp++;
                if (adAS == 2) absAS2++;
                if (adAS >= 3) { absAS_ge3++; if (mx <= 3) absAS_ge3_maxle3++; if (mx > maxAS_when_big_gap) maxAS_when_big_gap = mx; }
                if (adAS >= 2 && ex_printed < 3) {
                    printf("   example n=%d |A-S|=%d: p=%u q=%u a=%d s=%d w=%d dX=%d dXi=%d A=%d S=%d conv=%d\n",
                           n, adAS, p, q, a, s, w, I.dX, I.dXi, A, S, isconv);
                    ex_printed++;
                }
                int c1 = dAS + 40; if (c1 < 0) c1 = 0; if (c1 > 80) c1 = 80;
                int c2 = mx + 4; if (c2 < 0) c2 = 0; if (c2 > 15) c2 = 15;
                dist[c1][c2]++;
                int dm = I.dX > I.dXi ? I.dX : I.dXi;
                int ads = a - s < 0 ? s - a : a - s;
                if (!isconv) {
                    if (ads > 1 + dm) viol_P5++;
                    /* the session's code test (s84_energy.py): |a-s| <= 2+dm or (a <= 4+dXi and s <= 4+dX) */
                    if (!(ads <= 2 + dm || (a <= 4 + I.dXi && s <= 4 + I.dX))) viol_sessionproved++;
                    if (a < 64) lvl[a]++;
                }
                /* (C') */
                {
                    int okc = 1;
                    if (G[pstar] != G[p] || G[qstar] != G[q]) okc = 0;
                    if ((pstar & 1) || !(pstar >> (n - 1))) okc = 0;
                    if (((qstar & 1) || !(qstar >> (n - 1)))) okc = 0;
                    pinfo J; pair_info(pstar, qstar, n, &J);
                    if (!(J.a == s && J.s == a && J.dX == I.dXi && J.dXi == I.dX)) okc = 0;
                    if (!okc) viol_Cp++;
                }
            }
            i = k;
        }
        EYj[n] = E;
        double Y = (double)m;
        printf("n=%d |Y|=%zu E/|Y|=%.6f conv/|Y|=%.6f nonconv/|Y|=%.6f values/|Y|=%.6f\n", n, m, E / Y, conv / Y, nonconv / Y, nvals / Y);
        int maxmult = 0; for (int r = 1; r < 256; r++) if (hist[r]) maxmult = r;
        printf("   max multiplicity %d; multiplicity histogram:", maxmult);
        for (int r = 1; r < 256; r++) if (hist[r]) printf(" %d:%lld", r, hist[r]);
        printf("\n   values with exactly 2q preimages per |Y| (q=1..5):");
        for (int q = 1; q <= 5; q++) printf(" %.6f", hist[2 * q] / Y);
        printf("\n");
        /* smoothed energies via two pointers on the sorted values */
        {
            int Ks[4] = {2, 4, 8, 16};
            for (int kk = 0; kk < 4; kk++) {
                long long K = Ks[kk];
                long long EK = 0; size_t lo_i = 0, hi_i = 0;
                for (size_t z = 0; z < m; z++) {
                    long long gz = (long long)(keys[z] >> 32);
                    while ((long long)(keys[lo_i] >> 32) < gz - K) lo_i++;
                    if (hi_i < z) hi_i = z;
                    while (hi_i + 1 < m && (long long)(keys[hi_i + 1] >> 32) <= gz + K) hi_i++;
                    EK += (long long)(hi_i - lo_i + 1);
                }
                printf("   slope K=%lld: (E(K)-E(0))/(K|Y|) = %.6f\n", K, (double)(EK - E) / (K * Y));
            }
        }
        printf("   level law L_a (non-converse, per |Y|):");
        for (int a = 1; a < 64; a++) if (lvl[a]) printf(" %d:%.6f", a, lvl[a] / Y);
        printf("\n   level counts:");
        for (int a = 1; a < 64; a++) if (lvl[a]) printf(" %d:%lld", a, lvl[a]);
        printf("\n   ratios L_(a+1)/L_a:");
        for (int a = 1; a < 63; a++) if (lvl[a] && lvl[a + 1]) printf(" a=%d:%.4f", a, (double)lvl[a + 1] / lvl[a]);
        printf("\n   checks: misc(a,s>=1, X,Xi odd, bounds, |hXi+tauX|<=2^w+1) violations %lld; same-sign pairs %lld (with max(A,S)>2: %lld)\n",
               viol_misc, same_sign, same_sign_big);
        printf("   (C) violations %lld; sharp (|A-S|<=1 or max<=3) violations %lld; P5 (non-conv) violations %lld; session code-test violations %lld; (C') violations %lld\n",
               viol_C, viol_sharp, viol_P5, viol_sessionproved, viol_Cp);
        printf("   |A-S| = 2: %lld pairs; |A-S| >= 3: %lld pairs (with max(A,S) <= 3: %lld; largest max(A,S) among them %d)\n",
               absAS2, absAS_ge3, absAS_ge3_maxle3, maxAS_when_big_gap);
        printf("   joint distribution (A-S : max(A,S) -> count), |A-S| >= 2 only:");
        for (int c1 = 0; c1 <= 80; c1++) for (int c2 = 0; c2 < 16; c2++)
            if (dist[c1][c2] && (c1 - 40 >= 2 || c1 - 40 <= -2)) printf(" (%d:%d)->%lld", c1 - 40, c2 - 4, dist[c1][c2]);
        printf("\n   distribution of A-S over all pairs:");
        for (int c1 = 0; c1 <= 80; c1++) { long long tot = 0; for (int c2 = 0; c2 < 16; c2++) tot += dist[c1][c2]; if (tot) printf(" %d:%lld", c1 - 40, tot); }
        printf("\n");
        /* X_n energy, computed directly on all of X_n */
        {
            size_t mx = (size_t)1 << (n - 1);
            uint32_t lox = 1u << (n - 1);
            for (size_t z = 0; z < mx; z++) keys[z] = ((uint64_t)G[lox + z] << 32) | z;
            radix64(keys, tmp, mx);
            long long En = 0; size_t z = 0; int max_layers = 0; long long range_viol = 0;
            while (z < mx) {
                size_t k = z + 1;
                while (k < mx && (keys[k] >> 32) == (keys[z] >> 32)) k++;
                En += (long long)(k - z) * (k - z);
                /* layers present in this value class */
                uint64_t layermask = 0;
                for (size_t y = z; y < k; y++) {
                    uint32_t v = lox + (uint32_t)(keys[y] & 0xffffffffu);
                    int kk = __builtin_ctz(~v);         /* trailing ones */
                    int j = n - kk;                      /* layer Y_j (j = 0 for the all-ones word) */
                    if (v == (uint32_t)((1ull << n) - 1)) j = 0;
                    layermask |= 1ull << j;
                    if (j >= 2) {
                        uint32_t gv = (uint32_t)(keys[y] >> 32);
                        if (!(gv >= (1u << (j - 2)) + 1 && gv <= (1u << j) - 1)) range_viol++;
                    }
                }
                int nl = __builtin_popcountll(layermask);
                if (nl > max_layers) max_layers = nl;
                z = k;
            }
            long long Ssum = 0; for (int j = 2; j <= n; j++) Ssum += EYj[j];
            double Xn = (double)mx;
            printf("   X_n: E_n/|X_n| = %.6f; same-layer sum_j E^Y_j / |X_n| = %.6f; cross-layer (E_n - sum - 1)/|X_n| = %.6f;\n"
                   "        (E) E_n <= 1 + 2 sum: %s (ratio E_n/(1+2 sum) = %.4f); max #layers per value = %d; layer-range violations %lld\n",
                   En / Xn, Ssum / Xn, (En - Ssum - 1) / Xn, (En <= 1 + 2 * Ssum) ? "yes" : "NO", En / (1.0 + 2.0 * Ssum), max_layers, range_viol);
        }
        fflush(stdout);
    }
    return 0;
}
