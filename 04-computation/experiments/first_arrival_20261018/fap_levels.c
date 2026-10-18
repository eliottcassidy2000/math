/* fap_levels.c -- backward-level energy ladder of the stopped ROOT transport into a certified forward-closed frontier V.
   mac-mini-2026-10-18.  dep(y) = number of P-steps from seed y to its first entry into V (0 on V).
   V^(k) = {y : dep(y) <= k} is forward closed (P lowers dep by one), and stopping on first entry into V^(k) freezes each
   source x at the point of its path with dep = k (or at x itself if dep(x) < k).  Hence

       E^(k)(N) = sum_{dep(y)=k} J(y)^2 + sum_{dep(x)<k} c_0(x)^2,   J(y) = sum_{x in tree(y)} sigma(x->y) c_0(x),

   E^(0) = E_V (frozen energy on V), E^(1) = E_port (port-diagonal energy), E^(maxdep) = E_0 = sum_x c_0(x)^2.
   One pass: push the level-(k+1) currents one transport step to level k and add the level-k sources.
   Also reports, per level, the merge cross term Delta_k = E^(k) - E^(k+1) and its sibling-pair part.
   Usage: ./fap_levels Vfile N1 [N2 ...] */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int nV; static int64_t Vs[8192]; static int8_t *inVflag; static int64_t maxV;
static inline int inV(int64_t s) { return s <= maxV && inVflag[s]; }
static inline int64_t seedf(int64_t u) { while ((u & 7) == 5) u = (u - 1) >> 2; return u; }
static inline int wfun(int64_t m) { return (m % 3 == 2) ? 2 : -1; }
static inline int thr(int64_t m) { int r = (int)(m % 3); return r == 0 ? 1 : (r == 1 ? -2 : 1); }
static inline int64_t Pstep(int64_t x, int *sgn) { int64_t t = 3 * x + 1; int a = 0; while (!(t & 1)) { t >>= 1; a++; } *sgn = (a == 2) ? 1 : -1; return seedf(t); }
static int cmp64(const void *a, const void *b) { int64_t x = *(const int64_t *)a, y = *(const int64_t *)b; return x < y ? -1 : x > y; }

/* open-addressing map int64 key -> int64 value */
typedef struct { int64_t *k; int64_t *v; uint64_t cap; uint64_t n; } Map;
static void minit(Map *m, uint64_t cap) { m->cap = cap; m->n = 0; m->k = malloc(cap * 8); m->v = malloc(cap * 8); memset(m->k, 0, cap * 8); }
static void mclear(Map *m) { memset(m->k, 0, m->cap * 8); m->n = 0; }
static inline uint64_t hh(int64_t x) { uint64_t z = (uint64_t)x * 0x9E3779B97F4A7C15ULL; return z ^ (z >> 29); }
static void madd(Map *m, int64_t key, int64_t val) {
    uint64_t i = hh(key) & (m->cap - 1);
    while (m->k[i] != 0 && m->k[i] != key) i = (i + 1) & (m->cap - 1);
    if (m->k[i] == 0) { m->k[i] = key; m->v[i] = 0; m->n++; if (m->n * 2 > m->cap) { fprintf(stderr, "map full\n"); exit(5); } }
    m->v[i] += val;
}

int main(int argc, char **argv) {
    FILE *fv = fopen(argv[1], "r"); char line[256];
    while (fgets(line, sizeof line, fv)) { if (line[0] == '#') continue; long long v; if (sscanf(line, "%lld", &v) == 1) Vs[nV++] = v; }
    fclose(fv); qsort(Vs, nV, 8, cmp64); maxV = Vs[nV - 1];
    inVflag = calloc(maxV + 1, 1); for (int i = 0; i < nV; i++) inVflag[Vs[i]] = 1;
    for (int ai = 2; ai < argc; ai++) {
        int64_t N = atoll(argv[ai]), M = 4 * N;
        /* initial coefficients on seeds x <= M */
        int64_t nsrc = (M + 1) / 2;
        int32_t *c0 = calloc(nsrc, 4);
        for (int64_t m = 1; m <= M; m += 2) { int k = wfun(m) + (m < 2 * N ? thr(m) : 0); if (k) c0[(seedf(m) - 1) / 2] += k; }
        /* depth of every source, bucketed */
        int maxdep = 0; uint16_t *dep = malloc(2 * nsrc);
        for (int64_t i = 0; i < nsrc; i++) {
            if (!c0[i]) { dep[i] = 0xFFFF; continue; }
            int64_t x = 2 * i + 1; int d = 0, sg; while (!inV(x)) { x = Pstep(x, &sg); d++; }
            dep[i] = d; if (d > maxdep) maxdep = d;
        }
        double E0 = 0; double *below = calloc(maxdep + 2, sizeof(double)); /* below[k] = sum_{dep<k} c0^2 */
        double *bydep = calloc(maxdep + 1, sizeof(double));
        for (int64_t i = 0; i < nsrc; i++) if (c0[i]) { double c = c0[i]; E0 += c * c; bydep[dep[i]] += c * c; }
        for (int k = 1; k <= maxdep + 1; k++) below[k] = below[k - 1] + bydep[k - 1];
        Map cur, nxt; uint64_t cap = 1; while (cap < 4 * (uint64_t)nsrc) cap <<= 1; minit(&cur, cap); minit(&nxt, cap);
        double *Ek = calloc(maxdep + 1, sizeof(double)); double *sibk = calloc(maxdep + 1, sizeof(double));
        /* bucket lists by dep */
        int64_t *cntd = calloc(maxdep + 2, 8); for (int64_t i = 0; i < nsrc; i++) if (c0[i]) cntd[dep[i] + 1]++;
        for (int d = 1; d <= maxdep + 1; d++) cntd[d] += cntd[d - 1];
        int64_t *order = malloc(8 * cntd[maxdep + 1]); int64_t *pos = calloc(maxdep + 1, 8);
        for (int64_t i = 0; i < nsrc; i++) if (c0[i]) { int d = dep[i]; order[cntd[d] + pos[d]++] = i; }
        for (int k = maxdep; k >= 0; k--) {
            mclear(&nxt);
            /* push level k+1 currents one step: siblings (s, 2s+1) land on the same y with signs +, - */
            double sibcross = 0;
            for (uint64_t j = 0; j < cur.cap; j++) if (cur.k[j]) {
                int sg; int64_t y = Pstep(cur.k[j], &sg); madd(&nxt, y, sg * cur.v[j]);
            }
            /* sibling cross term at this level: 2 * sum over pairs (s,2s+1) at level k+1 of (+J(s))(-J(2s+1)) */
            for (uint64_t j = 0; j < cur.cap; j++) if (cur.k[j]) {
                int64_t s = cur.k[j]; if ((3 * s + 1) % 4 != 0) continue;   /* s is the a=2 member */
                int64_t t = 2 * s + 1; uint64_t i = hh(t) & (cur.cap - 1);
                while (cur.k[i] != 0 && cur.k[i] != t) i = (i + 1) & (cur.cap - 1);
                if (cur.k[i] == t) sibcross += 2.0 * (double)cur.v[j] * (double)(-cur.v[i]);
            }
            for (int64_t q = cntd[k]; q < cntd[k + 1]; q++) { int64_t i = order[q]; madd(&nxt, 2 * i + 1, c0[i]); }
            double e = 0; for (uint64_t j = 0; j < nxt.cap; j++) if (nxt.k[j]) { double c = nxt.v[j]; e += c * c; }
            Ek[k] = e + below[k]; sibk[k] = sibcross;
            Map t = cur; cur = nxt; nxt = t;
        }
        printf("N=%lld |V|=%d maxdep=%d  E_0/N=%.4f  E_V/N=%.4f\n", (long long)N, nV, maxdep, E0 / N, Ek[0] / N);
        printf("  k : E^(k)/N   Delta_k/N=E^(k)-E^(k+1)   sibling part/N\n");
        for (int k = 0; k <= maxdep; k++) {
            if (k <= 25 || k % 10 == 0 || k == maxdep)
                printf("  %3d: %.4f   %+.4f   %+.4f\n", k, Ek[k] / N, k < maxdep ? (Ek[k] - Ek[k + 1]) / N : 0.0, sibk[k] / N);
        }
        fflush(stdout);
        free(c0); free(dep); free(below); free(bydep); free(cur.k); free(cur.v); free(nxt.k); free(nxt.v); free(Ek); free(sibk); free(cntd); free(order); free(pos);
    }
    return 0;
}
