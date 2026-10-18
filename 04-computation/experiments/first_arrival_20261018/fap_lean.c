/* fap_lean.c -- memory-lean pair ledger of first arrivals into a stopping set V (mac-mini-2026-10-18-sixcolour).
   Same transport and sources as fap.c (THM-4667 / HYP-9290 / THM-4670); only the integer-clock quantities are kept:
     E_V  = sum_u c_V(u)^2,  G0 = sum_u sum_t h_t(u)^2,  G1 = sum_u sum_t h_t(u) h_(t+1)(u),
     E0   = sum_x c_0(x)^2 (source energy, -> 34/7 N),  SIB = sum over sibling SOURCE pairs (x, 2x+1), x = 1 mod 8, with a
            common entry seed, of c_0(x) c_0(2x+1) (-> 5/7 N),  Gam1 = G1 + SIB (lag-1 pairs other than sibling sources),
     HL   = E_V - G0 - 2 G1 (lags >= 2),  E+ = sum_u (sum_(x->u) c_0(x))^2 (unsigned arrivals).
   Ledger identity: E_V = E0 - 2 SIB + (G0 - E0) + 2 Gam1 + HL.
   Per source: entry index (int16) and literal clock (int16; the sign is (-1)^clock) and c_0 (int8): 5 bytes; the clock
   parities of all n <= 4N are bit-packed.  Usage: ./fap_lean Vfile KMAX  (N runs over 1 .. 2^(KMAX+1)-1). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXV 4096
#define MAXCLK 4096
static int nV; static int64_t Vs[MAXV]; static int8_t *inVflag; static int64_t maxV;
static inline int64_t seedf(int64_t u) { while ((u & 7) == 5) u = (u - 1) >> 2; return u; }
static inline int wfun(int64_t m) { return (m % 3 == 2) ? 2 : -1; }
static inline int thr(int64_t m) { int r = (int)(m % 3); return r == 0 ? 1 : (r == 1 ? -2 : 1); }
static int vidx(int64_t s) { int lo = 0, hi = nV - 1; while (lo <= hi) { int mid = (lo + hi) / 2; if (Vs[mid] == s) return mid; if (Vs[mid] < s) lo = mid + 1; else hi = mid - 1; } return -1; }
static inline int inV(int64_t s) { return s <= maxV && inVflag[s]; }
static int cmp64(const void *a, const void *b) { int64_t x = *(const int64_t *)a, y = *(const int64_t *)b; return x < y ? -1 : x > y; }
static uint8_t *pbits;
static inline int getpar(int64_t n) { return (pbits[n >> 3] >> (n & 7)) & 1; }
static inline void setpar(int64_t n, int v) { if (v) pbits[n >> 3] |= (uint8_t)(1 << (n & 7)); }

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s Vfile KMAX\n", argv[0]); return 1; }
    FILE *fv = fopen(argv[1], "r"); char line[256];
    while (fgets(line, sizeof line, fv)) { if (line[0] == '#') continue; long long v; if (sscanf(line, "%lld", &v) == 1) Vs[nV++] = v; }
    fclose(fv); qsort(Vs, nV, sizeof(int64_t), cmp64); maxV = Vs[nV - 1];
    inVflag = calloc(maxV + 1, 1); for (int i = 0; i < nV; i++) inVflag[Vs[i]] = 1;
    int KMAX = atoi(argv[2]);
    int64_t Nmax = (int64_t)1 << (KMAX + 1), M = 4 * Nmax;
    pbits = calloc(M / 8 + 2, 1);
    for (int64_t n = 2; n <= M; n++) { int64_t x = n; int st = 0; while (x >= n) { if (x & 1) x = (3 * x + 1) >> 1; else x >>= 1; st++; } setpar(n, (st + getpar(x)) & 1); }
    int64_t nsrc = (M + 1) / 2;
    int16_t *su = malloc(2 * nsrc), *sc = malloc(2 * nsrc); int8_t *c0s = calloc(nsrc, 1);
    if (!su || !sc || !c0s) { fprintf(stderr, "out of memory\n"); return 5; }
    int maxclk = 0;
    for (int64_t i = 0; i < nsrc; i++) {
        int64_t x = 2 * i + 1;
        if ((x & 7) == 5) { su[i] = -1; sc[i] = 0; continue; }
        int clk = 0;
        while (!inV(x)) {
            int64_t t = 3 * x + 1; int a = 0; while (!(t & 1)) { t >>= 1; a++; }
            int strips = 0; while ((t & 7) == 5) { t = (t - 1) >> 2; strips++; }
            x = t; clk += a + 2 * strips;
        }
        if (clk >= MAXCLK - 1) { fprintf(stderr, "increase MAXCLK\n"); return 3; }
        su[i] = (int16_t)vidx(x); sc[i] = (int16_t)clk; if (clk > maxclk) maxclk = clk;
    }
    fprintf(stderr, "V=%s |V|=%d: sources done, max clock %d\n", argv[1], nV, maxclk);
    int64_t *cV = calloc(nV, 8), *cplus = calloc(nV, 8), *h = calloc((size_t)MAXCLK * nV, 8);
    double EV = 0, G0 = 0, G1 = 0, E0 = 0, SIB = 0, Eplus = 0;
    int64_t R = 0, H4 = 0, H1 = 0; int ok = 1;
    enum { EVm, EVx, G0m, R1m, R1lo, R1hi, Z0m, Z0x, GAm, GAx, HLm, HLx, EPm, DIPm, E0lo, E0hi, SBlo, SBhi, NST };
    double st[64][NST]; int64_t cnt[64];
    for (int k = 0; k < 64; k++) { cnt[k] = 0; for (int q = 0; q < NST; q++) st[k][q] = 0;
        st[k][EVx] = st[k][R1hi] = st[k][Z0x] = st[k][GAx] = st[k][HLx] = st[k][E0hi] = st[k][SBhi] = -1e18; st[k][R1lo] = st[k][E0lo] = st[k][SBlo] = 1e18; }
    for (int64_t N = 1; N < Nmax; N++) {
        for (int tt = 0; tt < 3; tt++) {
            int64_t m; int amt;
            if (tt == 0) { m = 4 * N - 3; amt = wfun(m); } else if (tt == 1) { m = 4 * N - 1; amt = wfun(m); } else { m = 2 * N - 1; amt = thr(m); }
            int64_t x = seedf(m), i = (x - 1) / 2; int u = su[i], c = sc[i];
            int64_t d = (c & 1) ? -amt : amt;
            EV += (double)(2 * cV[u] * d + d * d); cV[u] += d;
            Eplus += (double)(2 * cplus[u] * amt + (int64_t)amt * amt); cplus[u] += amt;
            R += (getpar(Vs[u]) ? -1 : 1) * d;
            E0 += (double)(2 * c0s[i] * amt + amt * amt);
            { int64_t pi = -1; if ((x & 7) == 1) pi = x; else if ((x & 15) == 3) pi = (x - 1) / 4;   /* index of 2x+1, resp. of (x-1)/2 */
              if (pi >= 0 && pi < nsrc && su[pi] == u) SIB += (double)amt * c0s[pi]; }
            c0s[i] += amt;
            int64_t *hh = h + (size_t)c * nV + u;
            int64_t lo = c > 0 ? h[(size_t)(c - 1) * nV + u] : 0, hi = h[(size_t)(c + 1) * nV + u];
            G1 += (double)(d * (lo + hi)); G0 += (double)(2 * (*hh) * d + d * d); *hh += d;
        }
        for (int64_t k = 4 * N - 3; k <= 4 * N; k++) H4 += wfun(k) * (getpar(k) ? -1 : 1);
        H1 += wfun(N) * (getpar(N) ? -1 : 1);
        if (R != H4 - H1) ok = 0;
        int k = 0; { int64_t xx = N; while (xx > 1) { xx >>= 1; k++; } }
        double *s = st[k]; double n = (double)N;
        double ev = EV / n, r1 = G0 > 0 ? G1 / G0 : 0, z0 = (G0 - E0) / n, ga = (G1 + SIB) / n, hl = (EV - G0 - 2 * G1) / n, e0 = E0 / n, sb = SIB / n;
        s[EVm] += ev; if (ev > s[EVx]) s[EVx] = ev; s[G0m] += G0 / n;
        s[R1m] += r1; if (r1 < s[R1lo]) s[R1lo] = r1; if (r1 > s[R1hi]) s[R1hi] = r1;
        s[Z0m] += z0; if (z0 > s[Z0x]) s[Z0x] = z0; s[GAm] += ga; if (ga > s[GAx]) s[GAx] = ga;
        s[HLm] += hl; if (hl > s[HLx]) s[HLx] = hl; s[EPm] += Eplus / n; s[DIPm] += Eplus > 0 ? EV / Eplus : 0;
        if (e0 < s[E0lo]) s[E0lo] = e0; if (e0 > s[E0hi]) s[E0hi] = e0; if (sb < s[SBlo]) s[SBlo] = sb; if (sb > s[SBhi]) s[SBhi] = sb;
        cnt[k]++;
    }
    printf("V=%s |V|=%d KMAX=%d invariant R(N)=H(4N)-H(N) at every N: %s\n", argv[1], nV, KMAX, ok ? "OK" : "FAILED");
    printf(" oct | E_V/N mean max | G0/N mean | G1/G0 mean [lo, hi] | (G0-E0)/N mean max | Gam1/N mean max | HL/N mean max | E+/N mean | E_V/E+ mean | E0/N range | SIB/N range\n");
    for (int k = 8; k <= KMAX; k++) { double *s = st[k]; double c = (double)cnt[k];
        printf(" %2d | %.3f %.3f | %.3f | %+.4f [%+.4f, %+.4f] | %+.3f %+.3f | %+.3f %+.3f | %+.3f %+.3f | %.3f | %.4f | [%.5f, %.5f] | [%.5f, %.5f]\n", k,
               s[EVm] / c, s[EVx], s[G0m] / c, s[R1m] / c, s[R1lo], s[R1hi], s[Z0m] / c, s[Z0x], s[GAm] / c, s[GAx], s[HLm] / c, s[HLx],
               s[EPm] / c, s[DIPm] / c, s[E0lo], s[E0hi], s[SBlo], s[SBhi]); }
    return 0;
}
