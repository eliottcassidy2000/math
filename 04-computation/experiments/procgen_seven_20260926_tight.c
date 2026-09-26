/*
 * procgen_seven_20260926_tight.c -- extract the tight core of a certificate written by
 * procgen_seven_20260926_game.c (lane "seven", session collatz-procgen-20260922, 2026-09-26).
 * Exploratory/structural tool: the certificates themselves are checked by ..._verify.c.
 *
 * mode lo: nodes = pairs P in W; tight edges P -> Q for the options Q of x = P + tau(P) H with
 *          g(Q) = g(P) + e(x): the cycles of this graph are exactly Min's best responses of density F
 *          against tau (every cycle of density exactly F against tau is tight).
 * mode up: nodes = pairs; tight edges P -> Q = target of a lift x of P under sigma with
 *          psi(Q) + e(x) = psi(P): its cycles are exactly the cycles of G_sigma of density F.
 * Both graphs are peeled (remove out-degree 0, then in-degree 0, repeatedly); the remaining core
 * (cycles and the tight paths between them) is printed as "E P Q b s" lines (b = lift, s = sign),
 * where x = P + b H and s in {+1, -1, 0 (even)}.
 * usage: procgen_seven_20260926_tight q k fn fd prefix lo|up
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

static long long q, N, H, fn, fd;
static uint8_t *tau, *sig;
static int16_t *g;
static int mode_lo;

static int exists(const char *prefix, const char *ext) {
    char f[1024]; snprintf(f, sizeof f, "%s.%s", prefix, ext);
    FILE *fp = fopen(f, "rb"); if (!fp) return 0; fclose(fp); return 1;
}

static void *load(const char *prefix, const char *ext, size_t bytes) {
    char f[1024]; snprintf(f, sizeof f, "%s.%s", prefix, ext);
    FILE *fp = fopen(f, "rb"); if (!fp) { printf("cannot open %s\n", f); exit(3); }
    void *p = malloc(bytes); if (!p || fread(p, 1, bytes, fp) != bytes) { printf("read error %s\n", f); exit(3); }
    fclose(fp); return p;
}

/* tight out-edges of P: up to 2 (lo: options of the tau-lift; up: sigma-targets of both lifts) */
static int edges(long long P, long long *Q, int *B, int *S) {
    int n = 0;
    long long e = (P & 1) ? fd - fn : -fn;
    if (g[P] < 0) return 0;
    if (mode_lo) {
        long long x = P + (tau[P] ? H : 0);
        if (!(x & 1)) { long long t = (x / 2) % H; if (g[t] >= 0 && g[t] == g[P] + e) { Q[n] = t; B[n] = tau[P]; S[n] = 0; n++; } }
        else for (int s = 1; s >= -1; s -= 2) {
            long long t = ((q * x + s) / 2) % H;
            if (g[t] >= 0 && g[t] == g[P] + e) { Q[n] = t; B[n] = tau[P]; S[n] = s; n++; }
        }
    } else {
        for (int b = 0; b < 2; b++) {
            long long x = P + (b ? H : 0), t; int s = 0;
            if (!(x & 1)) t = (x / 2) % H;
            else { s = sig[(x - 1) / 2] ? -1 : 1; t = ((q * x + s) / 2) % H; }
            if ((long long)g[t] + e == (long long)g[P]) {
                int dup = 0; for (int j = 0; j < n; j++) if (Q[j] == t) dup = 1;
                if (!dup) { Q[n] = t; B[n] = b; S[n] = s; n++; }
            }
        }
    }
    return n;
}

int main(int argc, char **argv) {
    if (argc != 7) { fprintf(stderr, "usage: %s q k fn fd prefix lo|up\n", argv[0]); return 2; }
    q = atoll(argv[1]); int k = atoi(argv[2]); fn = atoll(argv[3]); fd = atoll(argv[4]);
    N = 1LL << k; H = N / 2;
    mode_lo = !strcmp(argv[6], "lo");
    const char *pre = argv[5];
    if (exists(pre, "lo_taub") || exists(pre, "up_sigb")) {      /* compact format: expand */
        uint8_t *bits = load(pre, mode_lo ? "lo_taub" : "up_sigb", (H + 7) / 8);
        uint8_t *b = malloc(H);
        for (long long i = 0; i < H; i++) b[i] = (bits[i >> 3] >> (i & 7)) & 1;
        free(bits);
        if (mode_lo) tau = b; else sig = b;
        g = malloc(2 * H);
        const char *e8 = mode_lo ? "lo_g8" : "up_psi8", *e16 = mode_lo ? "lo_g16" : "up_psi16";
        if (exists(pre, e8)) { uint8_t *v = load(pre, e8, H); for (long long i = 0; i < H; i++) g[i] = v[i] == 255 ? -1 : v[i]; free(v); }
        else {
            uint16_t *v = load(pre, e16, 2 * H);
            for (long long i = 0; i < H; i++) { if (v[i] != 65535 && v[i] > 32000) { printf("potential too large\n"); exit(3); } g[i] = v[i] == 65535 ? -1 : (int16_t)v[i]; }
            free(v);
        }
    } else {
        int32_t *v = load(pre, mode_lo ? "lo_g" : "up_psi", 4 * H);
        if (mode_lo) tau = load(pre, "lo_tau", H); else sig = load(pre, "up_sig", H);
        g = malloc(2 * H);
        for (long long i = 0; i < H; i++) { if (v[i] > 32000) { printf("potential too large\n"); exit(3); } g[i] = (int16_t)(v[i] < 0 ? -1 : v[i]); }
        free(v);
    }
    uint8_t *od = calloc(H, 1), *id = calloc(H, 1), *dead = calloc(H, 1);
    int32_t *stk = malloc(4 * H);
    long long Q[2]; int B[2], S[2];
    for (long long P = 0; P < H; P++) {
        if (g[P] < 0) { dead[P] = 1; continue; }
        int n = edges(P, Q, B, S);
        od[P] = n;
        for (int j = 0; j < n; j++) if (id[Q[j]] < 255) id[Q[j]]++;
    }
    /* peel out-degree 0 */
    long long sp = 0;
    for (long long P = 0; P < H; P++) if (!dead[P] && od[P] == 0) stk[sp++] = P;
    long long qinv; { unsigned long long x = q; for (int i = 0; i < 6; i++) x *= 2ULL - q * x; qinv = (long long)(x & (N - 1)); }
    while (sp) {
        long long Qn = stk[--sp];
        if (dead[Qn]) continue;
        dead[Qn] = 1;
        long long pr[3] = { (2 * Qn) % H, (qinv * ((2 * Qn - 1 + H) % H)) % H, (qinv * ((2 * Qn + 1) % H)) % H };
        for (int j = 0; j < 3; j++) {
            long long P = pr[j];
            if (dead[P]) continue;
            int n = edges(P, Q, B, S);
            for (int i = 0; i < n; i++) if (Q[i] == Qn) { od[P]--; if (od[P] == 0) stk[sp++] = P; }
        }
    }
    /* recompute in-degrees among live nodes, peel in-degree 0 */
    memset(id, 0, H);
    for (long long P = 0; P < H; P++) if (!dead[P]) { int n = edges(P, Q, B, S); for (int j = 0; j < n; j++) if (!dead[Q[j]] && id[Q[j]] < 255) id[Q[j]]++; }
    sp = 0;
    for (long long P = 0; P < H; P++) if (!dead[P] && id[P] == 0) stk[sp++] = P;
    while (sp) {
        long long P = stk[--sp];
        if (dead[P]) continue;
        dead[P] = 1;
        int n = edges(P, Q, B, S);
        for (int j = 0; j < n; j++) if (!dead[Q[j]]) { id[Q[j]]--; if (id[Q[j]] == 0) stk[sp++] = Q[j]; }
    }
    long long alive = 0;
    for (long long P = 0; P < H; P++) if (!dead[P]) alive++;
    printf("CORE %lld\n", alive);
    if (alive <= 2000000)
        for (long long P = 0; P < H; P++) if (!dead[P]) {
            int n = edges(P, Q, B, S);
            for (int j = 0; j < n; j++) if (!dead[Q[j]]) printf("E %lld %lld %d %d\n", P, Q[j], B[j], S[j]);
        }
    return 0;
}
