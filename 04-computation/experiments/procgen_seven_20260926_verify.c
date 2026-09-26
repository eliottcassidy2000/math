/*
 * procgen_seven_20260926_verify.c -- independent exact checker of the two certificates written by
 * procgen_seven_20260926_game.c (lane "seven", session collatz-procgen-20260922, 2026-09-26).
 * It shares no code with the solver and does not use the negation symmetry: it walks all 2^k nodes.
 *
 * Level k, odd q, N = 2^k, H = 2^(k-1), F = fn/fd in lowest terms, weight e(x) = fd - fn (x odd),
 * -fn (x even).  For a node x (0 <= x < N) the target pair of the sign s is
 *     x even: (x/2) mod H;       x odd: ((q x + s)/2) mod H.
 *
 * LOWER certificate (tau(P) in {0,1}, g(P) >= 0 or "outside W"):  W = {P : g(P) defined} nonempty,
 *   and for every P in W, with x = P + tau(P) H, and every target pair Q of x (one if x even, both
 *   signs if x odd):  Q in W  and  g(Q) <= g(P) + e(x).
 *   ==> for every sign strategy, the play P_0 in W, P_(i+1) = sigma-target of P_i + tau(P_i) H stays in
 *   W, its nodes form a walk of the parity graph G_sigma, and its cycle has e-weight >= 0, i.e. odd
 *   density >= F.  So rho*(q,k) >= F.
 * UPPER certificate (sign of every odd node; psi(P) >= 0):  for every node x, with Q the target pair
 *   of x under its sign:  psi(Q) + e(x) <= psi(x mod H).
 *   ==> every cycle x_0 -> x_1 -> ... of G_sigma (x_(i+1) any lift of the target pair of x_i) has
 *   e-weight <= 0, i.e. odd density <= F.  So rho*(q,k) <= F.
 *
 * File formats (prefix.EXT):
 *   full:  lo_tau uint8[H], lo_g int32[H] (-1 = outside W), up_sig uint8[H] (index i = odd node 2i+1,
 *          1 = sign '-'), up_psi int32[H];
 *   lean:  lo_taub bits[H] (LSB first), lo_g8 uint8[H] (255 = outside W) or lo_g16 uint16[H]
 *          (65535 = outside W), up_sigb bits[H], up_psi8 uint8[H] or up_psi16 uint16[H].
 * The lower certificate is checked and freed before the upper one is loaded (peak memory ~ 1.1 H
 * bytes in the lean uint8 format).
 *
 * usage: procgen_seven_20260926_verify q k fn fd prefix
 * prints "LOWER OK|FAIL ...", "UPPER OK|FAIL ...", a 64-bit FNV-1a digest of the files, and
 * "RESULT CERTIFIED fn/fd" iff both hold (exit code 0).  If no upper-certificate files exist, only the lower
 * certificate is checked and "RESULT LOWER-CERTIFIED fn/fd" (rho*(q,k) >= fn/fd) is printed; symmetrically,
 * without lower-certificate files, "RESULT UPPER-CERTIFIED fn/fd" (rho*(q,k) <= fn/fd).
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

static long long gcdll(long long a, long long b) { while (b) { long long t = a % b; a = b; b = t; } return a < 0 ? -a : a; }
static uint64_t DIG = 1469598103934665603ULL;

static int exists(const char *prefix, const char *ext) {
    char fn[1024]; snprintf(fn, sizeof fn, "%s.%s", prefix, ext);
    FILE *f = fopen(fn, "rb"); if (!f) return 0; fclose(f); return 1;
}

static void *load(const char *prefix, const char *ext, size_t bytes) {
    char fn[1024];
    snprintf(fn, sizeof fn, "%s.%s", prefix, ext);
    FILE *f = fopen(fn, "rb");
    if (!f) { printf("cannot open %s\n", fn); exit(3); }
    void *p = malloc(bytes ? bytes : 1);
    if (!p) { printf("out of memory\n"); exit(3); }
    if (fread(p, 1, bytes, f) != bytes) { printf("short file %s\n", fn); exit(3); }
    if (fgetc(f) != EOF) { printf("long file %s\n", fn); exit(3); }
    fclose(f);
    const unsigned char *c = (const unsigned char *)p;
    for (size_t i = 0; i < bytes; i++) { DIG ^= c[i]; DIG *= 1099511628211ULL; }
    return p;
}

/* a potential array in one of three widths; val() returns -1 for "undefined" */
typedef struct { int w; const void *p; } Pot;
static inline long long val(Pot a, long long i) {
    if (a.w == 4) { int32_t v = ((const int32_t *)a.p)[i]; return v < 0 ? -1 : v; }
    if (a.w == 2) { uint16_t v = ((const uint16_t *)a.p)[i]; return v == 0xffff ? -1 : v; }
    uint8_t v = ((const uint8_t *)a.p)[i]; return v == 0xff ? -1 : v;
}
/* a 0/1 array as bytes (w = 8) or bits (w = 1) */
typedef struct { int w; const uint8_t *p; } Bits;
static inline int bit(Bits b, long long i) {
    if (b.w == 8) return b.p[i];
    return (b.p[i >> 3] >> (i & 7)) & 1;
}

int main(int argc, char **argv) {
    if (argc != 6) { fprintf(stderr, "usage: %s q k fn fd prefix\n", argv[0]); return 2; }
    long long q = atoll(argv[1]);
    int k = atoi(argv[2]);
    long long fn = atoll(argv[3]), fd = atoll(argv[4]);
    const char *prefix = argv[5];
    if (!(q & 1) || q < 1 || k < 3 || k > 30 || fd <= 0 || fn < 0 || fn > fd || gcdll(fn, fd) != 1) {
        printf("bad parameters\n"); return 2;
    }
    long long N = 1LL << k, H = N / 2;
    int lean = exists(prefix, "lo_taub") || exists(prefix, "up_sigb");
    int have_lower = exists(prefix, "lo_taub") || exists(prefix, "lo_tau");
    int lower_ok = 0;
    if (!have_lower) { printf("LOWER ABSENT  (upper-only certificate)\n"); goto upper; }

    /* ---- lower certificate ---- */
    Bits tau; Pot g;
    if (lean) {
        tau.w = 1; tau.p = (const uint8_t *)load(prefix, "lo_taub", (size_t)((H + 7) / 8));
        if (exists(prefix, "lo_g8")) { g.w = 1; g.p = load(prefix, "lo_g8", (size_t)H); }
        else { g.w = 2; g.p = load(prefix, "lo_g16", (size_t)H * 2); }
    } else {
        tau.w = 8; tau.p = (const uint8_t *)load(prefix, "lo_tau", (size_t)H);
        g.w = 4; g.p = load(prefix, "lo_g", (size_t)H * 4);
    }
    long long nW = 0, bad = 0, ntau1 = 0, gmax = -1;
    for (long long P = 0; P < H; P++) {
        long long gp = val(g, P);
        if (gp < 0) continue;
        nW++;
        int t = bit(tau, P);
        if (t > 1) { bad++; continue; }
        if (t) ntau1++;
        if (gp > gmax) gmax = gp;
        long long x = P + (t ? H : 0);
        long long e = (x % 2 == 1) ? fd - fn : -fn;
        long long tg[2]; int nt;
        if (x % 2 == 0) { tg[0] = (x / 2) % H; nt = 1; }
        else { tg[0] = ((q * x + 1) / 2) % H; tg[1] = ((q * x - 1) / 2) % H; nt = 2; }
        for (int j = 0; j < nt; j++) {
            long long gq = val(g, tg[j]);
            if (gq < 0 || gq > gp + e) { bad++; break; }
        }
    }
    lower_ok = (nW > 0 && bad == 0);
    printf("LOWER %s  (|W| = %lld of %lld pairs, violations %lld, max g %lld, tau=1 at %lld pairs)\n",
           lower_ok ? "OK" : "FAIL", nW, H, bad, gmax, ntau1);
    fflush(stdout);
    free((void *)tau.p); free((void *)g.p);

    /* ---- upper certificate (optional: a lower-only certificate has no upper files) ---- */
upper:
    if (!exists(prefix, "up_sigb") && !exists(prefix, "up_sig")) {
        printf("UPPER ABSENT  (lower-only certificate)\n");
        printf("DIGEST %016llx\n", (unsigned long long)DIG);
        printf("RESULT %s %lld/%lld\n", lower_ok ? "LOWER-CERTIFIED" : "NOT-CERTIFIED", fn, fd);
        return lower_ok ? 0 : 1;
    }
    Bits sig; Pot psi;
    if (lean) {
        sig.w = 1; sig.p = (const uint8_t *)load(prefix, "up_sigb", (size_t)((H + 7) / 8));
        if (exists(prefix, "up_psi8")) { psi.w = 1; psi.p = load(prefix, "up_psi8", (size_t)H); }
        else { psi.w = 2; psi.p = load(prefix, "up_psi16", (size_t)H * 2); }
    } else {
        sig.w = 8; sig.p = (const uint8_t *)load(prefix, "up_sig", (size_t)H);
        psi.w = 4; psi.p = load(prefix, "up_psi", (size_t)H * 4);
    }
    long long badu = 0, pmax = -1;
    for (long long P = 0; P < H; P++) {
        long long v = val(psi, P);
        if (v < 0) badu++;
        if (v > pmax) pmax = v;
        if (bit(sig, P) > 1) badu++;
    }
    for (long long x = 0; x < N && badu == 0; x++) {
        long long Q;
        long long e = (x % 2 == 1) ? fd - fn : -fn;
        if (x % 2 == 0) Q = (x / 2) % H;
        else {
            int minus = bit(sig, (x - 1) / 2);
            Q = ((q * x + (minus ? -1 : 1)) / 2) % H;
        }
        if (val(psi, Q) + e > val(psi, x % H)) badu++;
    }
    int upper_ok = (badu == 0);
    printf("UPPER %s  (violations %lld, max psi %lld)\n", upper_ok ? "OK" : "FAIL", badu, pmax);
    printf("DIGEST %016llx\n", (unsigned long long)DIG);
    if (!have_lower) {
        printf("RESULT %s %lld/%lld\n", upper_ok ? "UPPER-CERTIFIED" : "NOT-CERTIFIED", fn, fd);
        free((void *)sig.p); free((void *)psi.p);
        return upper_ok ? 0 : 1;
    }
    printf("RESULT %s %lld/%lld\n", (lower_ok && upper_ok) ? "CERTIFIED" : "NOT-CERTIFIED", fn, fd);
    free((void *)sig.p); free((void *)psi.p);
    return (lower_ok && upper_ok) ? 0 : 1;
}
