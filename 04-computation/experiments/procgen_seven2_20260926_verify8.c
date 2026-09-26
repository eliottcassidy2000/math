/*
 * procgen_seven2_20260926_verify8.c -- memory-lean exact checker of one-sided certificates in the seven lane's compact
 * uint8/uint16 format (lane "seven2", session collatz-procgen-20260922, 2026-09-26).  Written separately from the
 * solvers.  Potentials: prefix.lo_g8 / up_psi8 (uint8, 255 = undefined) or prefix.lo_g16 / up_psi16 (uint16,
 * 65535 = undefined), or prefix.lo_h8 (uint8 h = g - fn v(P), where v(P) is the 2-adic valuation of P and v(0) =
 * k-1; 255 = undefined).  Below, 'NONE' denotes the undefined value of the width in use.
 *
 * Level k, odd q, N = 2^k, H = 2^(k-1), F = fn/fd in lowest terms, e(x) = fd - fn (x odd), -fn (x even).
 * Target pairs of a node x: x/2 mod H (x even); (q x + s)/2 mod H, s = +1 or -1 (x odd).
 * LOWER (files prefix.lo_taub: H bits, prefix.lo_g8: H bytes, 255 = outside W):  W = {P : g(P) != 255} nonempty,
 *   and for every P in W, x = P + tau(P) H, and every target pair Q of x (both signs if x odd):
 *   g(Q) != 255 and g(Q) <= g(P) + e(x)            ==> rho*(q,k) >= F   (THM-4486 Lemma G2 / Prop. P).
 * UPPER (files prefix.up_sigb: bit i = sign '-' at the odd node 2i+1, prefix.up_psi8: H bytes):  every psi(P) != 255,
 *   and for every node x with target pair Q under its sign: psi(Q) + e(x) <= psi(x mod H)   ==> rho*(q,k) <= F.
 * Memory: the potential file is first checked to be invariant under P -> H - P (streamed: the half 0..H/2 is kept in
 * memory, the other half is compared entry by entry); then every inequality is checked for ALL pairs / ALL nodes,
 * reading potentials from the kept half (valid because of the verified invariance) and the bits in streaming order.
 * Peak memory about 2^(k-2) bytes per potential byte.  Prints "RESULT LOWER-CERTIFIED fn/fd" / "RESULT UPPER-CERTIFIED fn/fd" (exit 0)
 * or "RESULT NOT-CERTIFIED" (exit 1), and a 64-bit FNV-1a digest of the files.
 * usage: procgen_seven2_20260926_verify8 lower|upper q k fn fd prefix
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

typedef long long ll;
static uint64_t DIG = 1469598103934665603ULL;
static void dig(const uint8_t *p, size_t n) { for (size_t i = 0; i < n; i++) { DIG ^= p[i]; DIG *= 1099511628211ULL; } }
static ll gcdll(ll a, ll b) { while (b) { ll t = a % b; a = b; b = t; } return a < 0 ? -a : a; }

int main(int argc, char **argv) {
    if (argc != 7) { fprintf(stderr, "usage: %s lower|upper q k fn fd prefix\n", argv[0]); return 2; }
    int lower = !strcmp(argv[1], "lower");
    ll q = atoll(argv[2]); int k = atoi(argv[3]); ll fn = atoll(argv[4]), fd = atoll(argv[5]);
    const char *prefix = argv[6];
    if (!(q & 1) || k < 4 || k > 32 || fd <= 0 || fn < 0 || fn > fd || gcdll(fn, fd) != 1) { printf("bad parameters\n"); return 2; }
    ll N = 1LL << k, H = N / 2, HH = H / 2;
    char fb[1024], fp[1024];
    snprintf(fb, sizeof fb, "%s.%s", prefix, lower ? "lo_taub" : "up_sigb");
    snprintf(fp, sizeof fp, "%s.%s", prefix, lower ? "lo_g8" : "up_psi8");
    int W = 1, SH = 0;
    FILE *F1 = fopen(fp, "rb");
    if (!F1 && lower) {
        snprintf(fp, sizeof fp, "%s.lo_h8", prefix);
        F1 = fopen(fp, "rb"); SH = 1;
        if (!F1) { SH = 0; snprintf(fp, sizeof fp, "%s.%s", prefix, lower ? "lo_g8" : "up_psi8"); }
    }
    if (!F1) {
        snprintf(fp, sizeof fp, "%s.%s", prefix, lower ? "lo_g16" : "up_psi16");
        F1 = fopen(fp, "rb"); W = 2;
    }
    if (!F1) { printf("cannot open the potential file of %s\n", prefix); return 3; }
    const ll NONE = W == 1 ? 255 : 65535;
    fseek(F1, 0, SEEK_END);
    if (ftell(F1) != H * W) { printf("wrong size %s\n", fp); return 3; }
    fseek(F1, 0, SEEK_SET);
    /* the half 0..HH of the potential */
    uint8_t *g8 = 0; uint16_t *g16 = 0;
    if (W == 1) g8 = (uint8_t *)malloc((size_t)(HH + 1)); else g16 = (uint16_t *)malloc((size_t)(HH + 1) * 2);
    if (!g8 && !g16) { printf("out of memory\n"); return 3; }
    void *gp0 = W == 1 ? (void *)g8 : (void *)g16;
    if (fread(gp0, (size_t)W, (size_t)(HH + 1), F1) != (size_t)(HH + 1)) { printf("short read\n"); return 3; }
    dig((const uint8_t *)gp0, (size_t)(HH + 1) * W);
    /* stream the other half and compare with the mirror */
    const ll CH = 1 << 22;
    uint8_t *buf = (uint8_t *)malloc((size_t)CH * 2);
    ll asym = 0;
    for (ll lo = HH + 1; lo < H; lo += CH) {
        ll n = H - lo < CH ? H - lo : CH;
        if (fread(buf, (size_t)W, (size_t)n, F1) != (size_t)n) { printf("short read\n"); return 3; }
        dig(buf, (size_t)n * W);
        for (ll i = 0; i < n; i++) {
            ll a = W == 1 ? buf[i] : ((uint16_t *)buf)[i];
            ll b = W == 1 ? g8[H - (lo + i)] : g16[H - (lo + i)];
            if (a != b) asym++;
        }
    }
    fclose(F1);
    printf("potential file %s: mirror asymmetries %lld\n", fp, asym);
    if (asym) { printf("RESULT NOT-CERTIFIED\n"); return 1; }
#define RAW(P) ((ll)(W == 1 ? g8[(P) <= HH ? (P) : H - (P)] : g16[(P) <= HH ? (P) : H - (P)]))
#define VV(P) ((P) == 0 ? (ll)(k - 1) : (ll)__builtin_ctzll((unsigned long long)(P)))
#define GV(P) (SH ? (RAW(P) == NONE ? NONE : RAW(P) + fn * VV(P)) : RAW(P))
    FILE *F2 = fopen(fb, "rb");
    if (!F2) { printf("cannot open %s\n", fb); return 3; }
    fseek(F2, 0, SEEK_END);
    if (ftell(F2) != (H + 7) / 8) { printf("wrong size %s\n", fb); return 3; }
    fseek(F2, 0, SEEK_SET);
    ll bad = 0, nW = 0, maxv = -1;
    if (lower) {
        for (ll lo = 0; lo < H; lo += 8 * CH) {
            ll n = H - lo < 8 * CH ? H - lo : 8 * CH;
            if (fread(buf, 1, (size_t)((n + 7) / 8), F2) != (size_t)((n + 7) / 8)) { printf("short read\n"); return 3; }
            dig(buf, (size_t)((n + 7) / 8));
            for (ll i = 0; i < n; i++) {
                ll P = lo + i;
                ll gp = GV(P);
                if (gp == NONE) continue;
                nW++; if (gp > maxv) maxv = gp;
                int t = (buf[i >> 3] >> (i & 7)) & 1;
                ll x = P + (t ? H : 0);
                ll e = (x & 1) ? fd - fn : -fn;
                ll tg[2]; int nt;
                if (!(x & 1)) { tg[0] = (x / 2) % H; nt = 1; }
                else { tg[0] = ((q * x + 1) / 2) % H; tg[1] = ((q * x - 1) / 2) % H; nt = 2; }
                for (int j = 0; j < nt; j++) { ll gq = GV(tg[j]); if (gq == NONE || gq > gp + e) { bad++; break; } }
            }
        }
        int ok = nW > 0 && bad == 0;
        printf("LOWER %s  (|W| = %lld of %lld pairs, violations %lld, max g %lld)\n", ok ? "OK" : "FAIL", nW, H, bad, maxv);
        printf("DIGEST %016llx\n", (unsigned long long)DIG);
        printf("RESULT %s %lld/%lld\n", ok ? "LOWER-CERTIFIED" : "NOT-CERTIFIED", fn, fd);
        return ok ? 0 : 1;
    }
    for (ll P = 0; P <= HH; P++) { ll v = GV(P); if (v == NONE) bad++; if (v > maxv) maxv = v; }
    /* all nodes x in increasing order; the sign bit of odd x = 2i+1 is bit i of the stream */
    ll i0 = 0; ll have = 0;
    for (ll x = 0; x < N && bad == 0; x++) {
        ll Qt;
        if (!(x & 1)) Qt = (x / 2) % H;
        else {
            ll i = (x - 1) / 2;
            if (i >= i0 + have) {
                i0 = i; ll n = H - i0 < 8 * CH ? H - i0 : 8 * CH;
                if (fread(buf, 1, (size_t)((n + 7) / 8), F2) != (size_t)((n + 7) / 8)) { printf("short read\n"); return 3; }
                dig(buf, (size_t)((n + 7) / 8)); have = n;
            }
            ll j = i - i0;
            int minus = (buf[j >> 3] >> (j & 7)) & 1;
            Qt = ((q * x + (minus ? -1 : 1)) / 2) % H;
        }
        ll e = (x & 1) ? fd - fn : -fn;
        if (GV(Qt) + e > GV(x % H)) bad++;
    }
    int ok = bad == 0;
    printf("UPPER %s  (violations %lld, max psi %lld)\n", ok ? "OK" : "FAIL", bad, maxv);
    printf("DIGEST %016llx\n", (unsigned long long)DIG);
    printf("RESULT %s %lld/%lld\n", ok ? "UPPER-CERTIFIED" : "NOT-CERTIFIED", fn, fd);
    return ok ? 0 : 1;
}
