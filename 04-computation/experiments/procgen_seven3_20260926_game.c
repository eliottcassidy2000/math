/*
 * procgen_seven3_20260926_game.c -- least fixed points of Min's and Max's operators of the level-k parity game of
 * q n +- 1 (lane "seven3", session collatz-procgen-20260922, 2026-09-26).  Written independently of the seven and
 * floor engines, as a cross-check of the arena semantics (THM-4486 section 2, Lemma L).
 *
 * Arena.  N = 2^k, H = 2^(k-1).  Node s in Z/N.  Options: s even -> pair s/2 mod H; s odd -> pairs (q s + 1)/2 and
 * (q s - 1)/2 mod H (signs + and -).  A pair P has the two lifts P, P + H.  Threshold F = fn/fd:
 * e(s) = fd - fn (s odd), -fn (s even).
 *   Min:  f(s) = max(0, e(s) + min_{options P} max(f(P), f(P+H)))       (finite lfp => upper certificate, Lemma G1)
 *   Max:  g(s) = max(0, -e(s) + max_{options P} min(g(P), g(P+H)))      (finite lfp => lower certificate, Lemma G2)
 * Chaotic (FIFO worklist) iteration from 0; a value above `cap` stops the run (status DIVERGED; caps only make the
 * run fail, never certify).  Output files (prefix):  .pot int32[N] (the fixed point),  .sig uint8[H] (Min: 1 iff the
 * argmin sign of odd node 2i+1 is '-'),  .lift uint8[H] (Max: argmin lift bit of pair P).  The caller re-checks every
 * certificate with exact integer arithmetic.
 *
 * usage: game min|max q k fn fd cap prefix [mhtie]   (mhtie: Min takes the MH sign wherever it is admissible)
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

typedef int64_t i64;
static i64 Q, N, H, QI;
static int32_t *f;
static uint8_t *inq;
static i64 *queue;
static i64 qh, qt, qn;

static void push(i64 s) { if (!inq[s]) { inq[s] = 1; queue[qt] = s; qt = (qt + 1) % qn; } }

int main(int argc, char **argv) {
    if (argc < 8) { fprintf(stderr, "usage: %s min|max q k fn fd cap prefix\n", argv[0]); return 2; }
    int ismin = !strcmp(argv[1], "min");
    Q = atoll(argv[2]); int k = atoi(argv[3]);
    i64 fn = atoll(argv[4]), fd = atoll(argv[5]), cap = atoll(argv[6]);
    const char *pre = argv[7];
    int mhtie = (argc > 8 && !strcmp(argv[8], "mhtie"));
    N = 1LL << k; H = N >> 1;
    /* inverse of Q mod N */
    QI = 1; for (int i = 0; i < 64; i++) QI = QI * (2 - Q * QI);
    f = calloc(N, sizeof(int32_t)); inq = calloc(N, 1); qn = N + 1; queue = malloc(qn * sizeof(i64));
    if (!f || !inq || !queue) { fprintf(stderr, "oom\n"); return 2; }
    for (i64 s = 0; s < N; s++) push(s);
    i64 ep = fd - fn, ee = -fn;
    int diverged = 0;
    while (qh != qt) {
        i64 s = queue[qh]; qh = (qh + 1) % qn; inq[s] = 0;
        i64 best;
        if (!(s & 1)) {
            i64 P = (s >> 1) & (H - 1);
            i64 a = f[P], b = f[P + H];
            best = ismin ? (a > b ? a : b) + ee : -ee + (a < b ? a : b);
        } else {
            i64 Pp = ((Q * s + 1) >> 1) & (H - 1), Pm = ((Q * s - 1) >> 1) & (H - 1);
            i64 ap = f[Pp], bp = f[Pp + H], am = f[Pm], bm = f[Pm + H];
            if (ismin) {
                i64 xp = ap > bp ? ap : bp, xm = am > bm ? am : bm;
                best = (xp < xm ? xp : xm) + ep;
            } else {
                i64 xp = ap < bp ? ap : bp, xm = am < bm ? am : bm;
                best = -ep + (xp > xm ? xp : xm);
            }
        }
        if (best < 0) best = 0;
        if (best > f[s]) {
            if (best > cap) { diverged = 1; break; }
            f[s] = (int32_t)best;
            /* predecessors of the pair P = s mod H: the even node 2P and the odd nodes with option P */
            i64 P = s & (H - 1);
            push((2 * P) & (N - 1));
            push((((2 * P - 1) * QI) & (N - 1)));   /* q s' + 1 = 2P  (mod N) */
            push((((2 * P + 1) * QI) & (N - 1)));   /* q s' - 1 = 2P */
        }
    }
    if (diverged) { printf("DIVERGED\n"); return 0; }
    i64 mx = 0; for (i64 s = 0; s < N; s++) if (f[s] > mx) mx = f[s];
    char fnm[4096];
    snprintf(fnm, sizeof fnm, "%s.pot", pre); FILE *o = fopen(fnm, "wb"); fwrite(f, sizeof(int32_t), N, o); fclose(o);
    if (ismin) {
        uint8_t *sg = calloc(H, 1);
        for (i64 i = 0; i < H; i++) {
            i64 s = 2 * i + 1;
            i64 Pp = ((Q * s + 1) >> 1) & (H - 1), Pm = ((Q * s - 1) >> 1) & (H - 1);
            i64 xp = f[Pp] > f[Pp + H] ? f[Pp] : f[Pp + H], xm = f[Pm] > f[Pm + H] ? f[Pm] : f[Pm + H];
            if (mhtie) {   /* the MH sign whenever it satisfies the potential inequality, else the other sign */
                int mhminus = ((s & 3) == 3);
                i64 xmh = mhminus ? xm : xp;
                sg[i] = (xmh + ep <= f[s]) ? mhminus : !mhminus;
            } else sg[i] = (xm < xp) ? 1 : 0;
        }
        snprintf(fnm, sizeof fnm, "%s.sig", pre); o = fopen(fnm, "wb"); fwrite(sg, 1, H, o); fclose(o);
    } else {
        uint8_t *lf = calloc(H, 1);
        for (i64 P = 0; P < H; P++) lf[P] = (f[P + H] < f[P]) ? 1 : 0;
        snprintf(fnm, sizeof fnm, "%s.lift", pre); o = fopen(fnm, "wb"); fwrite(lf, 1, H, o); fclose(o);
    }
    printf("FINITE %lld\n", (long long)mx);
    return 0;
}
