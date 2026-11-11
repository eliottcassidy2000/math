/* audit_cube.c  --  independent audit program (not derived from the author's code; no Gauss-sum or rank shortcut).
 *
 * For a base tournament R on m vertices and an out-set X of a new vertex v:
 *     H(R + v_X) = sum_{a in X} St(a) + sum_{b notin X} En(b) + sum_{b notin X, a in X} M(b,a),
 * where St(a) / En(b) = number of Hamiltonian paths of R starting at a / ending at b, and
 *     M(b,a) = sum over S with b in S, a notin S of  (#Ham. paths of R[S] ending at b) * (#Ham. paths of R[V\S] starting at a)
 * (v first, v last, or v in the middle of the path).  f[S][b], g[S][a] by subset DP (arithmetic mod 2^32, enough for H mod 8).
 * The program evaluates H(R + v_X) mod 8 for ALL 2^m out-sets X (Gray code, direct evaluation of the formula above),
 * so  g_R = sum_X chi(R + v_X),  chi = +1 if H = 1 mod 4, -1 if H = 3 mod 4,  is an exact integer.
 *
 * Modes:
 *   census m           classes of order m from `gentourng -q m` on stdin ('1' at (i,j), i<j, read as i -> j);
 *                      weight m!/|Aut| (backtracking).  Prints the exact labelled counts of (m+1)-tournaments by H mod 8,
 *                      the numerator of D_(m+1), sum w |g_R| and sum w g_R^2.
 *   mc m N seed        N uniformly random labelled bases (xoshiro256**, seeded by splitmix64(seed)).
 *                      Prints sum g_R, sum g_R^2, sum |g_R|, the mod-8 totals.
 *   selftest m N seed  compares, for N random bases, the formula with a direct subset-DP count of H(R + v_X) for every X.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXM 12
static int m;
static unsigned outm[MAXM + 1], inm[MAXM + 1];
static unsigned char adj[MAXM + 1][MAXM + 1];
static int score[MAXM + 1];
static uint32_t f[1 << MAXM][MAXM], g[1 << MAXM][MAXM];
static uint32_t St[MAXM], En[MAXM], M[MAXM][MAXM];
static int popc(unsigned x) { return __builtin_popcount(x); }

static void build_forms(void) {
    unsigned full = (1u << m) - 1;
    for (unsigned S = 1; S <= full; S++) {
        if ((S & (S - 1)) == 0) { int b = __builtin_ctz(S); for (int j = 0; j < m; j++) f[S][j] = g[S][j] = 0; f[S][b] = g[S][b] = 1; continue; }
        unsigned rem = S;
        for (int j = 0; j < m; j++) f[S][j] = g[S][j] = 0;
        while (rem) {
            int b = __builtin_ctz(rem); rem &= rem - 1;
            unsigned T = S & ~(1u << b);
            uint32_t sf = 0, sg = 0;
            unsigned pin = inm[b] & T;            /* c -> b */
            while (pin) { int c = __builtin_ctz(pin); pin &= pin - 1; sf += f[T][c]; }
            unsigned pout = outm[b] & T;          /* b -> c */
            while (pout) { int c = __builtin_ctz(pout); pout &= pout - 1; sg += g[T][c]; }
            f[S][b] = sf; g[S][b] = sg;
        }
    }
    for (int a = 0; a < m; a++) { St[a] = g[full][a]; En[a] = f[full][a]; for (int b = 0; b < m; b++) M[a][b] = 0; }
    for (unsigned S = 1; S < full; S++) {
        unsigned comp = full & ~S, rb = S;
        while (rb) {
            int b = __builtin_ctz(rb); rb &= rb - 1;
            uint32_t fb = f[S][b]; if (!fb) continue;
            unsigned ra = comp;
            while (ra) { int a = __builtin_ctz(ra); ra &= ra - 1; M[b][a] += fb * g[comp][a]; }
        }
    }
}
/* sum over all X of chi, and the mod 8 histogram; returns g_R */
static long cube(uint64_t h8[8]) {
    uint32_t H = 0; for (int b = 0; b < m; b++) H += En[b];        /* X empty */
    unsigned X = 0; long gs = 0;
    unsigned long total = 1ul << m;
    for (unsigned long k = 0;; k++) {
        unsigned r = H & 7; h8[r]++; gs += ((r & 3) == 1) ? 1 : -1;
        if (k + 1 == total) break;
        int j = __builtin_ctzl(k + 1);
        /* derivative of H in x_j with the other coordinates fixed */
        uint32_t d = St[j] - En[j];
        for (int b = 0; b < m; b++) if (b != j) { if (!((X >> b) & 1)) d += M[b][j]; else d -= M[j][b]; }
        if ((X >> j) & 1) { H -= d; } else { H += d; }
        X ^= 1u << j;
    }
    return gs;
}
/* direct H of a tournament given by out-masks on nn vertices */
static uint32_t dpH[1 << (MAXM + 1)][MAXM + 1];
static uint32_t ham_direct(const unsigned *om, int nn) {
    unsigned full = (1u << nn) - 1;
    for (unsigned S = 1; S <= full; S++) for (int v = 0; v < nn; v++) dpH[S][v] = 0;
    for (int v = 0; v < nn; v++) dpH[1u << v][v] = 1;
    for (unsigned S = 1; S <= full; S++) { unsigned rem = S; while (rem) { int l = __builtin_ctz(rem); rem &= rem - 1; uint32_t c = dpH[S][l]; if (!c) continue;
        unsigned nx = om[l] & ~S; while (nx) { int w = __builtin_ctz(nx); nx &= nx - 1; dpH[S | (1u << w)][w] += c; } } }
    uint32_t H = 0; for (int v = 0; v < nn; v++) H += dpH[full][v];
    return H;
}

static int img[MAXM]; static unsigned used; static long autcount;
static void aut_bt(int v) {
    if (v == m) { autcount++; return; }
    for (int w = 0; w < m; w++) {
        if ((used >> w) & 1) continue;
        if (score[w] != score[v]) continue;
        int ok = 1;
        for (int u = 0; u < v; u++) if (adj[u][v] != adj[img[u]][w]) { ok = 0; break; }
        if (!ok) continue;
        img[v] = w; used |= 1u << w; aut_bt(v + 1); used &= ~(1u << w);
    }
}

static uint64_t rs[4];
static uint64_t splitmix(uint64_t *x) { uint64_t z = (*x += 0x9E3779B97F4A7C15ull); z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull; z = (z ^ (z >> 27)) * 0x94D049BB133111EBull; return z ^ (z >> 31); }
static uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static uint64_t rnd(void) { uint64_t r = rotl(rs[1] * 5, 7) * 9, t = rs[1] << 17; rs[2] ^= rs[0]; rs[3] ^= rs[1]; rs[1] ^= rs[2]; rs[0] ^= rs[3]; rs[2] ^= t; rs[3] = rotl(rs[3], 45); return r; }
static void random_base(void) {
    for (int i = 0; i < m; i++) { outm[i] = inm[i] = 0; }
    uint64_t bits = 0; int nb = 0;
    for (int i = 0; i < m; i++) for (int j = i + 1; j < m; j++) {
        if (!nb) { bits = rnd(); nb = 64; }
        int b = (int)(bits & 1); bits >>= 1; nb--;
        if (b) { outm[i] |= 1u << j; inm[j] |= 1u << i; } else { outm[j] |= 1u << i; inm[i] |= 1u << j; }
    }
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: audit_cube census m | mc m N seed | selftest m N seed\n"); return 1; }
    m = atoi(argv[2]);
    if (m < 1 || m > MAXM) return 1;
    if (strcmp(argv[1], "census") == 0) {
        int E = m * (m - 1) / 2; uint64_t mfact = 1; for (int i = 2; i <= m; i++) mfact *= i;
        uint64_t H8[8]; memset(H8, 0, sizeof H8);
        unsigned __int128 Sabs = 0, Ssq = 0; uint64_t classes = 0, totW = 0, zeroW = 0; __int128 Ssum = 0;
        char line[256];
        while (fgets(line, sizeof line, stdin)) {
            int len = 0; while (line[len] == '0' || line[len] == '1') len++;
            if (len != E) continue;
            int pos = 0;
            for (int i = 0; i < m; i++) { outm[i] = inm[i] = 0; for (int j = 0; j < m; j++) adj[i][j] = 0; }
            for (int i = 0; i < m; i++) for (int j = i + 1; j < m; j++) {
                int b = line[pos++] == '1';
                if (b) { adj[i][j] = 1; outm[i] |= 1u << j; inm[j] |= 1u << i; } else { adj[j][i] = 1; outm[j] |= 1u << i; inm[i] |= 1u << j; }
            }
            for (int i = 0; i < m; i++) score[i] = popc(outm[i]);
            used = 0; autcount = 0; aut_bt(0);
            uint64_t w = mfact / (uint64_t)autcount;
            build_forms();
            uint64_t h8[8]; memset(h8, 0, sizeof h8);
            long gs = cube(h8);
            for (int r = 0; r < 8; r++) H8[r] += w * h8[r];
            Ssum += (__int128)w * gs; Sabs += (unsigned __int128)w * (uint64_t)(gs < 0 ? -gs : gs); Ssq += (unsigned __int128)w * (uint64_t)(gs * gs);
            if (gs == 0) zeroW += w;
            classes++; totW += w;
        }
        printf("base order m = %d: classes %llu, total weight %llu (2^%d = %llu)\n", m, (unsigned long long)classes, (unsigned long long)totW, E, 1ull << E);
        printf("labelled tournaments of order %d by H mod 8:  1: %llu  3: %llu  5: %llu  7: %llu  (even residues: %llu)\n", m + 1,
               (unsigned long long)H8[1], (unsigned long long)H8[3], (unsigned long long)H8[5], (unsigned long long)H8[7],
               (unsigned long long)(H8[0] + H8[2] + H8[4] + H8[6]));
        long long num = (long long)Ssum;
        printf("numerator of D_%d over 2^%d: %lld\n", m + 1, E + m, num);
        printf("sum_R w |g_R| = %llu%018llu (decimal split 10^18), as double %.12e ; sum_R w g_R^2 as double %.12e ; weight of bases with g_R = 0: %llu\n",
               (unsigned long long)(Sabs / 1000000000000000000ull), (unsigned long long)(Sabs % 1000000000000000000ull), (double)Sabs, (double)Ssq, (unsigned long long)zeroW);
        return 0;
    }
    long N = argc > 3 ? atol(argv[3]) : 1000; uint64_t seed = argc > 4 ? strtoull(argv[4], NULL, 10) : 1;
    uint64_t sx = seed; for (int i = 0; i < 4; i++) rs[i] = splitmix(&sx);
    if (strcmp(argv[1], "selftest") == 0) {
        long bad = 0, tested = 0;
        for (long it = 0; it < N; it++) {
            random_base(); build_forms();
            uint32_t H0 = 0; for (int b = 0; b < m; b++) H0 += En[b];
            for (unsigned X = 0; X < (1u << m); X++) {
                unsigned om[MAXM + 1];
                for (int a = 0; a < m; a++) om[a] = outm[a] | (((X >> a) & 1) ? 0 : (1u << m));
                om[m] = X;
                uint32_t Hd = ham_direct(om, m + 1);
                uint32_t Hf = 0;
                for (int a = 0; a < m; a++) { if ((X >> a) & 1) Hf += St[a]; else Hf += En[a]; }
                for (int b = 0; b < m; b++) if (!((X >> b) & 1)) for (int a = 0; a < m; a++) if ((X >> a) & 1) Hf += M[b][a];
                if (Hd != Hf) bad++;
                tested++;
            }
            /* Gray-code evaluation against the direct formula: totals must agree */
            uint64_t h8[8]; memset(h8, 0, sizeof h8); long gs = cube(h8), gd = 0;
            for (unsigned X = 0; X < (1u << m); X++) {
                uint32_t Hf = 0;
                for (int a = 0; a < m; a++) { if ((X >> a) & 1) Hf += St[a]; else Hf += En[a]; }
                for (int b = 0; b < m; b++) if (!((X >> b) & 1)) for (int a = 0; a < m; a++) if ((X >> a) & 1) Hf += M[b][a];
                gd += ((Hf & 3) == 1) ? 1 : -1;
            }
            if (gs != gd) bad += 1000000;
            (void)H0;
        }
        printf("selftest m=%d: %ld (base, X) pairs compared with a direct subset-DP count of H; mismatches %ld\n", m, tested, bad);
        return 0;
    }
    /* mc */
    double S1 = 0, S2 = 0, SA = 0; uint64_t H8[8]; memset(H8, 0, sizeof H8); long nzero = 0, npos = 0;
    double scale = 1.0 / (double)(1ul << m);
    for (long it = 0; it < N; it++) {
        random_base(); build_forms();
        long gs = cube(H8);
        double b = gs * scale; S1 += b; S2 += b * b; SA += (b < 0 ? -b : b);
        if (gs == 0) nzero++; if (gs > 0) npos++;
    }
    printf("m %d N %ld seed %llu S1 %.12e S2 %.12e SA %.12e nzero %ld npos %ld h8 %llu %llu %llu %llu\n", m, N, (unsigned long long)seed, S1, S2, SA, nzero, npos,
           (unsigned long long)H8[1], (unsigned long long)H8[3], (unsigned long long)H8[5], (unsigned long long)H8[7]);
    return 0;
}
