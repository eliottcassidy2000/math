/*
 * procgen_seven3_20260926_adv.c -- exact value of ONE lift strategy (Max) in the level-k parity game of q n +- 1
 * (lane "seven3", session collatz-procgen-20260922, 2026-09-26).
 *
 * Given tau : Z/H -> {0,1} (file of H bytes), Min's graph: node s in Z/N; s even -> t = P + tau(P) H with P = s/2 mod H;
 * s odd -> t_+ and t_- for the options P_+- = (q s +- 1)/2 mod H.  value(tau) = max over start nodes of the least
 * cycle density Min can reach = the largest F such that some nonempty W, closed under all Min moves, has only
 * cycles of density >= F (Lemma G2 of THM-4486 with this tau).
 * Method: Howard policy iteration (max mean of the EVEN indicator = 1 - min mean of the odd indicator) gives for
 * every node the least reachable cycle density; F = the maximum over nodes (an exact cycle a/p).  Certificate: W =
 * {nodes with value >= F} (float with tolerance, then checked EXACTLY: closed under all Min moves) and the least
 * fixed point f(s) = max(0, -e(s) + max_{moves s->t} f(t)) on W with exact int64 arithmetic, checked on every edge
 * (f(t) <= f(s) + e(s), e = p - a odd, -a even).  Output: "VALUE a p", "W size", "CERT ok maxf", "CYCLE ..." (a cycle
 * of density exactly a/p that Min can reach from the best start: the witness that tau does not give more).
 * usage: adv q k liftfile
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

typedef int64_t i64;

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: adv q k liftfile\n"); return 2; }
    i64 Q = atoll(argv[1]); int k = atoi(argv[2]);
    i64 N = 1LL << k, H = N >> 1;
    uint8_t *tau = malloc(H);
    FILE *fi = fopen(argv[3], "rb"); if (!fi || fread(tau, 1, H, fi) != (size_t)H) { fprintf(stderr, "bad liftfile\n"); return 2; } fclose(fi);
    int32_t *s0 = malloc(N * 4), *s1 = malloc(N * 4);   /* successors (s1 = -1 for even) */
    for (i64 s = 0; s < N; s++) {
        if (!(s & 1)) { i64 P = (s >> 1) & (H - 1); s0[s] = (int32_t)(P + tau[P] * H); s1[s] = -1; }
        else {
            i64 Pp = ((Q * s + 1) >> 1) & (H - 1), Pm = ((Q * s - 1) >> 1) & (H - 1);
            s0[s] = (int32_t)(Pp + tau[Pp] * H); s1[s] = (int32_t)(Pm + tau[Pm] * H);
        }
    }
    /* Howard: maximize mean of g = 1 - [odd] */
    int32_t *pol = malloc(N * 4); double *eta = malloc(N * 8), *h = malloc(N * 8);
    int8_t *st = malloc(N); int32_t *path = malloc(N * 4);
    for (i64 s = 0; s < N; s++) pol[s] = s0[s];
    for (int it = 0; it < 100000; it++) {
        memset(st, 0, N);
        for (i64 a0 = 0; a0 < N; a0++) {
            if (st[a0]) continue;
            i64 pn = 0; int32_t u = (int32_t)a0;
            while (!st[u]) { st[u] = 1; path[pn++] = u; u = pol[u]; }
            if (st[u] == 1) {
                i64 j = 0; while (path[j] != u) j++;
                i64 L = pn - j, g = 0; for (i64 t = j; t < pn; t++) g += !(path[t] & 1);
                double m = (double)g / (double)L;
                h[u] = 0; eta[u] = m;
                for (i64 t = pn - 1; t > j; t--) { int32_t v = path[t]; eta[v] = m; h[v] = (double)!(v & 1) - m + h[pol[v]]; }
                for (i64 t = j; t < pn; t++) st[path[t]] = 2;
                pn = j;
            }
            for (i64 t = pn - 1; t >= 0; t--) { int32_t v = path[t]; eta[v] = eta[pol[v]]; h[v] = (double)!(v & 1) - eta[v] + h[pol[v]]; st[v] = 2; }
        }
        int changed = 0; const double eps = 1e-10;
        for (i64 v = 0; v < N; v++) if (s1[v] >= 0) {
            int32_t o = (pol[v] == s0[v]) ? s1[v] : s0[v];
            if (eta[o] > eta[pol[v]] + eps) { pol[v] = o; changed = 1; }
        }
        if (changed) continue;
        for (i64 v = 0; v < N; v++) if (s1[v] >= 0) {
            int32_t o = (pol[v] == s0[v]) ? s1[v] : s0[v];
            if (fabs(eta[o] - eta[pol[v]]) < eps && h[o] > h[pol[v]] + eps) { pol[v] = o; changed = 1; }
        }
        if (!changed) break;
    }
    /* value = max over nodes of (1 - eta) ; the node with the least eta */
    i64 best = 0;
    for (i64 v = 1; v < N; v++) if (eta[v] < eta[best]) best = v;
    /* exact cycle from best under the policy */
    i64 pn = 0; int32_t u = (int32_t)best; memset(st, 0, N);
    while (!st[u]) { st[u] = 1; path[pn++] = u; u = pol[u]; }
    i64 j = 0; while (path[j] != u) j++;
    i64 L = pn - j, a = 0; for (i64 t = j; t < pn; t++) a += (path[t] & 1);
    i64 A = a, P = L; { i64 g = A, b = P; while (b) { i64 t = g % b; g = b; b = t; } A /= g; P /= g; }
    double Fd = (double)A / (double)P;
    /* W = {1 - eta >= F - tol} */
    uint8_t *inW = malloc(N); i64 wsz = 0;
    for (i64 v = 0; v < N; v++) { inW[v] = (1.0 - eta[v] >= Fd - 1e-9); wsz += inW[v]; }
    for (i64 v = 0; v < N; v++) if (inW[v]) { if (!inW[s0[v]] || (s1[v] >= 0 && !inW[s1[v]])) { printf("W NOT CLOSED\n"); return 1; } }
    /* lfp of f(s) = max(0, -e(s) + max f(t)) on W; worklist with predecessor scan via reverse lists */
    int32_t *roff = calloc(N + 1, 4), *rdst = malloc(2 * N * 4);
    for (i64 v = 0; v < N; v++) if (inW[v]) { roff[s0[v] + 1]++; if (s1[v] >= 0) roff[s1[v] + 1]++; }
    for (i64 v = 0; v < N; v++) roff[v + 1] += roff[v];
    { int32_t *fill = calloc(N, 4); for (i64 v = 0; v < N; v++) if (inW[v]) { rdst[roff[s0[v]] + fill[s0[v]]++] = (int32_t)v; if (s1[v] >= 0) rdst[roff[s1[v]] + fill[s1[v]]++] = (int32_t)v; } free(fill); }
    i64 *f = calloc(N, 8); int32_t *q = malloc((N + 1) * 4); uint8_t *inq = calloc(N, 1);
    i64 qh = 0, qt = 0, qn = N + 1;
    for (i64 v = 0; v < N; v++) if (inW[v]) { q[qt] = (int32_t)v; qt = (qt + 1) % qn; inq[v] = 1; }
    i64 bound = (i64)wsz * (A > 0 ? A : 1) + 1;
    while (qh != qt) {
        int32_t v = q[qh]; qh = (qh + 1) % qn; inq[v] = 0;
        i64 e = (v & 1) ? P - A : -A;
        i64 m = f[s0[v]]; if (s1[v] >= 0 && f[s1[v]] > m) m = f[s1[v]];
        i64 val = -e + m;
        if (val > f[v]) {
            f[v] = val;
            if (val > bound) { printf("CERT FAIL (diverges)\n"); return 1; }
            for (int32_t r = roff[v]; r < roff[v + 1]; r++) { int32_t x = rdst[r]; if (!inq[x]) { inq[x] = 1; q[qt] = x; qt = (qt + 1) % qn; } }
        }
    }
    i64 mx = 0;
    for (i64 v = 0; v < N; v++) if (inW[v]) {
        i64 e = (v & 1) ? P - A : -A;
        if (f[s0[v]] > f[v] + e || (s1[v] >= 0 && f[s1[v]] > f[v] + e)) { printf("CERT FAIL\n"); return 1; }
        if (f[v] > mx) mx = f[v];
    }
    /* maximality: from EVERY node Min can reach a cycle of density <= A/P.  Least fixed point of
       h(s) = max(0, e(s) + min_{moves s->t} h(t)) over all nodes (exact int64); finite iff Min can keep the e-sum
       bounded from every start, i.e. every start reaches a cycle of density <= A/P.  Bound: a finite lfp is at most
       N * max|e|, so exceeding it proves divergence (then this claim fails). */
    {
        int32_t *roff2 = calloc(N + 1, 4), *rdst2 = malloc(2 * N * 4);
        for (i64 v = 0; v < N; v++) { roff2[s0[v] + 1]++; if (s1[v] >= 0) roff2[s1[v] + 1]++; }
        for (i64 v = 0; v < N; v++) roff2[v + 1] += roff2[v];
        { int32_t *fill = calloc(N, 4); for (i64 v = 0; v < N; v++) { rdst2[roff2[s0[v]] + fill[s0[v]]++] = (int32_t)v; if (s1[v] >= 0) rdst2[roff2[s1[v]] + fill[s1[v]]++] = (int32_t)v; } free(fill); }
        i64 *hh = calloc(N, 8); i64 qh2 = 0, qt2 = 0;
        memset(inq, 0, N);
        for (i64 v = 0; v < N; v++) { q[qt2++] = (int32_t)v; inq[v] = 1; }
        i64 bound2 = N * (P > A ? P - A : A) + 1;
        int bad = 0;
        while (qh2 != qt2) {
            int32_t v = q[qh2]; qh2 = (qh2 + 1) % qn; inq[v] = 0;
            i64 e = (v & 1) ? P - A : -A;
            i64 m = hh[s0[v]]; if (s1[v] >= 0 && hh[s1[v]] < m) m = hh[s1[v]];
            i64 val = e + m;
            if (val > hh[v]) {
                hh[v] = val;
                if (val > bound2) { bad = 1; break; }
                for (int32_t r = roff2[v]; r < roff2[v + 1]; r++) { int32_t x = rdst2[r]; if (!inq[x]) { inq[x] = 1; q[qt2] = x; qt2 = (qt2 + 1) % qn; } }
            }
        }
        if (!bad) {
            for (i64 v = 0; v < N; v++) {
                i64 e = (v & 1) ? P - A : -A;
                i64 m = hh[s0[v]]; if (s1[v] >= 0 && hh[s1[v]] < m) m = hh[s1[v]];
                if (e + m > hh[v]) { bad = 1; break; }
            }
        }
        printf("MAXIMAL %s\n", bad ? "FAIL" : "ok");
    }
    printf("VALUE %lld %lld\n", (long long)A, (long long)P);
    printf("W %lld\n", (long long)wsz);
    printf("CERT ok %lld\n", (long long)mx);
    printf("CYCLE");
    for (i64 t = j; t < pn && t < j + 100000; t++) printf(" %d", path[t]);
    printf("\n");
    return 0;
}
