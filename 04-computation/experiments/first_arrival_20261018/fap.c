/* fap.c -- first-arrival port decomposition of the stopped ROOT transport on a certified forward-closed frontier V
   (clock_duality_oct15 R6; V_1000 = 627 seeds).  mac-mini-2026-10-18 (first-arrival cancellation brief).

   Definitions (THM-4667 / HYP-9290): seeds are odd s with s % 8 != 5; seed(x) strips B(x) = 4x+1;
   P(s) = seed((3s+1)/2^a), a = nu_2(3s+1) in {1,2}, sign (-1)^a; literal shortcut clock of the step = a + 2*(#strips).
   Sources: every odd m <= 4N contributes K_N(m) = w(m)[m <= 4N] + thr(m)[m < 2N] to the seed seed(m)
   (w(m) = 2 if m%3==2 else -1; thr(m) = 1, -2, 1 for m%3 = 0, 1, 2), i.e. K_N(m) = -3 chi_3(m) below 2N.
   Every source seed x is followed to its FIRST ENTRY u into V; c_V(u) = sum sigma(x->u) c_0(x).
   The entry PORT of x is the last outside seed on its path (x itself if P(x) in V); sources already in V are DIRECT.
   Since V is forward closed, R(N) = sum_u f(u) c_V(u) with f(u) = (-1)^clock(u), checked against H(4N) - H(N).

   Tracked at every N (incremental, three source entries change per N):
     E_V    = sum_u c_V(u)^2                               (frozen energy)
     E_port = sum_u D(u)^2 + sum_ports A(p)^2              (port-diagonal: no cross terms between different ports)
     E_pair = sum_u D(u)^2 + sum_pairs (A(s)+A(2s+1))^2    (sibling pairs merged first; unpaired ports alone)
     E_dep  = sum_j sum_u g_j(u)^2                         (two-time Gram diagonal over first-arrival depth j)
     E_clk  = sum_t sum_u h_t(u)^2                         (diagonal over the literal integer clock t of the arrival)
   Usage: ./fap Vfile KMAX [dumpN ...]     (N runs over 1 .. 2^(KMAX+1)-1)
   With dumpN arguments, the per-depth Gram matrix and the per-port/pair currents are printed at those N. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#define MAXV 4096
#define MAXDEP 512
#define MAXCLK 4096
static int nV; static int64_t Vs[MAXV];
static int8_t *inVflag; static int64_t maxV;

static inline int64_t seedf(int64_t u) { while ((u & 7) == 5) u = (u - 1) >> 2; return u; }
static inline int wfun(int64_t m) { return (m % 3 == 2) ? 2 : -1; }
static inline int thr(int64_t m) { int r = (int)(m % 3); return r == 0 ? 1 : (r == 1 ? -2 : 1); }
static int vidx(int64_t s) { int lo = 0, hi = nV - 1; while (lo <= hi) { int mid = (lo + hi) / 2; if (Vs[mid] == s) return mid; if (Vs[mid] < s) lo = mid + 1; else hi = mid - 1; } return -1; }
static inline int inV(int64_t s) { return s <= maxV && inVflag[s]; }

/* ports: sorted list of outside fibre elements of V seeds, with pair ids */
static int64_t *ports; static int nports; static int *pairid; static int npairs; static int *portu;
static int pidx(int64_t s) { int lo = 0, hi = nports - 1; while (lo <= hi) { int mid = (lo + hi) / 2; if (ports[mid] == s) return mid; if (ports[mid] < s) lo = mid + 1; else hi = mid - 1; } return -1; }
static int cmp64(const void *a, const void *b) { int64_t x = *(const int64_t *)a, y = *(const int64_t *)b; return x < y ? -1 : x > y; }

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s Vfile KMAX [dumpN...]\n", argv[0]); return 1; }
    FILE *fv = fopen(argv[1], "r"); char line[256];
    while (fgets(line, sizeof line, fv)) { if (line[0] == '#') continue; long long v, r; if (sscanf(line, "%lld %lld", &v, &r) >= 1) Vs[nV++] = v; }
    fclose(fv);
    qsort(Vs, nV, sizeof(int64_t), cmp64);
    maxV = Vs[nV - 1];
    inVflag = calloc(maxV + 1, 1); for (int i = 0; i < nV; i++) inVflag[Vs[i]] = 1;
    int KMAX = atoi(argv[2]);
    int64_t Nmax = (int64_t)1 << (KMAX + 1), M = 4 * Nmax;
    int ndump = argc - 3; int64_t *dumpN = malloc(sizeof(int64_t) * (ndump + 1));
    for (int i = 0; i < ndump; i++) dumpN[i] = atoll(argv[3 + i]);

    /* clock parities par[n], n <= M */
    int8_t *par = malloc(M + 1); par[0] = 0; par[1] = 0;
    for (int64_t n = 2; n <= M; n++) { int64_t x = n; int st = 0; while (x >= n) { if (x & 1) x = (3 * x + 1) >> 1; else x >>= 1; st++; } par[n] = (st + par[x]) & 1; }

    /* ports: for each u in V, fibre elements s_j = (2^a B^j(u) - 1)/3 <= M that are outside V */
    int cap = 1 << 20; ports = malloc(sizeof(int64_t) * cap); nports = 0;
    for (int i = 0; i < nV; i++) {
        int64_t u = Vs[i], b = u;
        for (int j = 0; j < 40; j++) {
            if (b > ((int64_t)1 << 60)) break;          /* paths may climb far above M before their first entry */
            int r = (int)(b % 3);
            if (r != 0) {
                int64_t s = (r == 1) ? (4 * b - 1) / 3 : (2 * b - 1) / 3;
                if (!(s == u && u == 1) && !inV(s)) { if (nports == cap) { cap *= 2; ports = realloc(ports, sizeof(int64_t) * cap); } ports[nports++] = s; }
            }
            b = 4 * b + 1;
        }
    }
    qsort(ports, nports, sizeof(int64_t), cmp64);
    pairid = malloc(sizeof(int) * nports); portu = malloc(sizeof(int) * nports);
    for (int i = 0; i < nports; i++) pairid[i] = -1;
    npairs = 0;
    for (int i = 0; i < nports; i++) {
        int64_t s = ports[i];
        /* entry seed of the port */
        int64_t t = 3 * s + 1; while (!(t & 1)) t >>= 1; t = seedf(t); portu[i] = vidx(t);
        if (pairid[i] >= 0) continue;
        int k = pidx(2 * s + 1);
        if (k >= 0 && (3 * s + 1) % 4 == 0) { pairid[i] = npairs; pairid[k] = npairs; npairs++; }   /* s has a = 2, 2s+1 has a = 1 */
        else pairid[i] = npairs++;
    }
    int npaired = 0; for (int i = 0; i < nports; i++) { int k = pidx(2 * ports[i] + 1); if (k >= 0 && pairid[k] == pairid[i]) npaired++; }
    fprintf(stderr, "V: %d seeds (max %lld); ports <= M=%lld: %d (in %d sibling pairs: %d ports paired); pair ids %d\n",
            nV, (long long)maxV, (long long)M, nports, npaired, 2 * npaired, npairs);

    /* per-source records for every odd seed x <= M */
    int64_t nsrc = (M + 1) / 2;
    int16_t *su = malloc(sizeof(int16_t) * nsrc); int8_t *ss = malloc(nsrc); int32_t *sp = malloc(sizeof(int32_t) * nsrc);
    uint16_t *sd = malloc(sizeof(uint16_t) * nsrc); int16_t *sc = malloc(sizeof(int16_t) * nsrc);
    int maxdep = 0, maxclk = 0;
    for (int64_t i = 0; i < nsrc; i++) {
        int64_t x = 2 * i + 1;
        if ((x & 7) == 5) { su[i] = -1; continue; }
        int sign = 1, dep = 0, clk = 0; int64_t prev = -1;
        while (!inV(x)) {
            int64_t t = 3 * x + 1; int a = 0; while (!(t & 1)) { t >>= 1; a++; }
            int strips = 0; while ((t & 7) == 5) { t = (t - 1) >> 2; strips++; }
            prev = x; x = t; dep++; clk += a + 2 * strips; if (a != 2) sign = -sign;
        }
        su[i] = (int16_t)vidx(x); ss[i] = (int8_t)sign; sd[i] = (uint16_t)dep; sc[i] = (int16_t)clk; if (clk > 32000) { fprintf(stderr, "clock overflow\n"); return 4; }
        sp[i] = (prev < 0) ? -1 : pidx(prev);
        if (prev >= 0 && sp[i] < 0) { fprintf(stderr, "port %lld missing\n", (long long)prev); return 2; }
        if (dep > maxdep) maxdep = dep; if (clk > maxclk) maxclk = clk;
    }
    fprintf(stderr, "sources done: max depth %d, max clock %d\n", maxdep, maxclk);
    if (maxdep >= MAXDEP || maxclk >= MAXCLK) { fprintf(stderr, "increase MAXDEP/MAXCLK\n"); return 3; }

    int64_t *cV = calloc(nV, sizeof(int64_t)), *D = calloc(nV, sizeof(int64_t));
    int64_t *A = calloc(nports, sizeof(int64_t)), *S = calloc(npairs, sizeof(int64_t));
    int64_t *g = calloc((size_t)MAXDEP * nV, sizeof(int64_t)), *h = calloc((size_t)MAXCLK * nV, sizeof(int64_t));
    double EV = 0, Eport = 0, Epair = 0, Edep = 0, Eclk = 0, G1 = 0; /* exact as long as < 2^53; G1 = lag-1 clock Gram */
    int64_t R = 0, H4 = 0, H1 = 0; int ok = 1;
    double mx[5][64], sm[5][64]; int64_t cnt[64], am[5][64]; double xp_max[64], xd_max[64], xp_sum[64], xd_sum[64], xp_hi[64], xp_lo[64], r1_sum[64], r1_lo[64], r1_hi[64];
    for (int k = 0; k < 64; k++) { xp_max[k] = xd_max[k] = xp_sum[k] = xd_sum[k] = 0; xp_hi[k] = -1e9; xp_lo[k] = 1e9; r1_sum[k] = 0; r1_lo[k] = 1e9; r1_hi[k] = -1e9; }
    for (int q = 0; q < 5; q++) for (int k = 0; k < 64; k++) { mx[q][k] = -1; sm[q][k] = 0; am[q][k] = 0; }
    memset(cnt, 0, sizeof cnt);
    int di = 0;
    for (int64_t N = 1; N < Nmax; N++) {
        for (int t = 0; t < 3; t++) {
            int64_t m; int amt;
            if (t == 0) { m = 4 * N - 3; amt = wfun(m); } else if (t == 1) { m = 4 * N - 1; amt = wfun(m); } else { m = 2 * N - 1; amt = thr(m); }
            int64_t x = seedf(m); int64_t i = (x - 1) / 2;
            int u = su[i]; int64_t d = (int64_t)ss[i] * amt;
            EV += (double)(2 * cV[u] * d + d * d); cV[u] += d;
            R += (par[Vs[u]] ? -1 : 1) * d;
            if (sp[i] < 0) { Eport += (double)(2 * D[u] * d + d * d); Epair += (double)(2 * D[u] * d + d * d); D[u] += d; }
            else { int p = sp[i], q = pairid[p];
                   Eport += (double)(2 * A[p] * d + d * d); A[p] += d;
                   Epair += (double)(2 * S[q] * d + d * d); S[q] += d; }
            int64_t *gg = g + (size_t)sd[i] * nV + u; Edep += (double)(2 * (*gg) * d + d * d); *gg += d;
            int64_t *hh = h + (size_t)sc[i] * nV + u;
            { int64_t lo = sc[i] > 0 ? h[(size_t)(sc[i] - 1) * nV + u] : 0, hi = h[(size_t)(sc[i] + 1) * nV + u]; G1 += (double)(d * (lo + hi)); }
            Eclk += (double)(2 * (*hh) * d + d * d); *hh += d;
        }
        for (int64_t k = 4 * N - 3; k <= 4 * N; k++) H4 += wfun(k) * (par[k] ? -1 : 1);
        H1 += wfun(N) * (par[N] ? -1 : 1);
        if (R != H4 - H1) ok = 0;
        int k = 0; { int64_t xx = N; while (xx > 1) { xx >>= 1; k++; } }
        double vals[5] = {EV / N, Eport / N, Epair / N, Edep / N, Eclk / N};
        for (int q = 0; q < 5; q++) { if (vals[q] > mx[q][k]) { mx[q][k] = vals[q]; am[q][k] = N; } sm[q][k] += vals[q]; }
        cnt[k]++;
        { double xp = (EV - Eport) / N, xd = (EV - Edep) / N; { double r1 = Eclk > 0 ? G1 / Eclk : 0; r1_sum[k] += r1; if (r1 < r1_lo[k]) r1_lo[k] = r1; if (r1 > r1_hi[k]) r1_hi[k] = r1; } if (xp > xp_hi[k]) xp_hi[k] = xp; if (xp < xp_lo[k]) xp_lo[k] = xp; if (fabs(xp) > xp_max[k]) xp_max[k] = fabs(xp); if (fabs(xd) > xd_max[k]) xd_max[k] = fabs(xd); xp_sum[k] += xp; xd_sum[k] += xd; }
        if (di < ndump && N == dumpN[di]) {
            di++;
            printf("== dump N=%lld: E_V/N %.4f  E_port/N %.4f  E_pair/N %.4f  E_dep/N %.4f  E_clk/N %.4f  R^2/N %.4f\n", (long long)N,
                   EV / N, Eport / N, Epair / N, Edep / N, Eclk / N, (double)R * R / N);
            /* depth Gram matrix (normalised by N), depths 0..min(maxdep, 60) */
            int J = maxdep < 60 ? maxdep : 60; double diag = 0, off = 0, offneg = 0;
            for (int a = 0; a <= J; a++) for (int b = 0; b <= J; b++) {
                double G = 0; for (int uu = 0; uu < nV; uu++) G += (double)g[(size_t)a * nV + uu] * g[(size_t)b * nV + uu];
                if (a == b) diag += G; else { off += G; if (G < 0) offneg += G; }
            }
            printf("   depth Gram: diag %.4f  off-diag %.4f (negative part %.4f)  [/N]\n", diag / N, off / N, offneg / N);
            /* adjacent-depth correlations */
            printf("   <g_j,g_(j+1)>/sqrt(|g_j|^2|g_(j+1)|^2) for j = 0..20:");
            for (int a = 0; a < 21 && a < maxdep; a++) { double G = 0, n1 = 0, n2 = 0;
                for (int uu = 0; uu < nV; uu++) { double x1 = g[(size_t)a * nV + uu], x2 = g[(size_t)(a + 1) * nV + uu]; G += x1 * x2; n1 += x1 * x1; n2 += x2 * x2; }
                printf(" %.2f", (n1 > 0 && n2 > 0) ? G / sqrt(n1 * n2) : 0.0); }
            printf("\n");
            /* integer-clock lag Gram: G(D) = sum_u sum_t h_t(u) h_(t+D)(u), D = 0..12 ([/N]); E_V = G(0) + 2 sum_(D>0) G(D) */
            printf("   clock-lag Gram G(D)/N, D = 0..12:");
            { double tot = 0; for (int D = 0; D <= 12; D++) { double G = 0;
                for (int t = 0; t + D <= maxclk; t++) for (int uu = 0; uu < nV; uu++) G += (double)h[(size_t)t * nV + uu] * h[(size_t)(t + D) * nV + uu];
                printf(" %+.3f", G / N); tot += (D ? 2 : 1) * G; }
              double rest = EV - tot; printf("  | lags > 12 contribute %+.3f\n", rest / N); }
            /* sibling pair anti-correlation: over pairs with both ports, corr(A(s), A(2s+1)) and energy of halves vs merged */
            double sxy = 0, sxx = 0, syy = 0, merged = 0; int np2 = 0;
            for (int p = 0; p < nports; p++) { int kq = pidx(2 * ports[p] + 1); if (kq >= 0 && pairid[kq] == pairid[p]) {
                double x1 = A[p], x2 = A[kq]; sxy += x1 * x2; sxx += x1 * x1; syy += x2 * x2; merged += (x1 + x2) * (x1 + x2); np2++; } }
            printf("   sibling pairs (%d): sum A(s)^2+A(2s+1)^2 = %.4f N, merged sum (A(s)+A(2s+1))^2 = %.4f N, corr %.3f\n",
                   np2, (sxx + syy) / N, merged / N, (sxx > 0 && syy > 0) ? sxy / sqrt(sxx * syy) : 0.0);
            fflush(stdout);
        }
    }
    printf("V=%s |V|=%d KMAX=%d invariant R(N)=H(4N)-H(N) at every N: %s\n", argv[1], nV, KMAX, ok ? "OK" : "FAILED");
    const char *nm[5] = {"E_V", "E_port", "E_pair", "E_dep", "E_clk"};
    for (int k = 8; k <= KMAX; k++) {
        printf("octave 2^%2d:", k);
        for (int q = 0; q < 5; q++) printf("  %s max %.3f mean %.3f", nm[q], mx[q][k], sm[q][k] / cnt[k]);
        printf("  | cross-port (E_V-E_port)/N mean %+.3f range [%+.3f, %+.3f] ; cross-depth (E_V-E_dep)/N mean %+.3f max|.| %.3f ; lag-1 clock ratio G1/G0 mean %+.4f range [%+.4f, %+.4f]\n",
               xp_sum[k] / cnt[k], xp_lo[k], xp_hi[k], xd_sum[k] / cnt[k], xd_max[k], r1_sum[k] / cnt[k], r1_lo[k], r1_hi[k]);
    }
    return 0;
}
