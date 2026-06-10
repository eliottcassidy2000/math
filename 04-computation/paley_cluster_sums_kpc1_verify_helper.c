/* paley_cluster_sums_kpc1_verify_helper.c
 * ADVERSARIAL VERIFIER kind-pasteur-2026-06-10-S1 (independent of worker code).
 *
 * Computes distinct-tuple Legendre character sums for Paley primes by DIRECT
 * DFS ENUMERATION (no Moebius/partition inversion -- different method from the
 * worker's engine):
 *
 *   SUM over injective tuples (x_0..x_{N-1}) in F_p of prod over run-edges
 *   chi(x_{i+1}-x_i),  runs given by prev[] arrays.
 *
 * Symmetry reductions (proved, classical):
 *   - translation x -> x+t: weight invariant; fix x_0 = 0, multiply by p.
 *   - scaling x -> lambda*x (lambda != 0): weight *= chi(lambda)^{#edges};
 *     for EVEN #edges weight invariant; affine group acts FREELY on injective
 *     tuples, so additionally fix x_1 = 1 and multiply by p(p-1).
 *     (x_1 is adjacent to x_0 in run 1; chi(1-0)=chi(1)=+1, no weight lost.)
 *   - last level: sum over the final point c of chi(c - x_prev) over c not yet
 *     used equals MINUS the sum over used points (since sum over all c is 0).
 *
 * Output lines:  NAME p VALUE   (VALUE = full unreduced integer sum)
 */
#include <stdio.h>
#include <stdint.h>
#include <string.h>

static int P;
static int chid[40][40];     /* chid[a][b] = chi((b-a) mod p) */
static int NP;
static const int *PREV;
static int X[16];
static unsigned MASK;
static long long TOTAL;

static void rec(int i, int sign) {
    if (i == NP - 1) {
        /* last point: has exactly one edge, to X[PREV[i]] (true for all
           structures here: every run ends the placement order). */
        int xp = X[PREV[i]];
        long long s = 0;
        for (int j = 0; j < i; j++) s += chid[xp][X[j]];
        TOTAL += (long long)sign * (-s);
        return;
    }
    int pr = PREV[i];
    if (pr < 0) {
        for (int c = 0; c < P; c++) {
            if (MASK & (1u << c)) continue;
            X[i] = c; MASK |= (1u << c);
            rec(i + 1, sign);
            MASK &= ~(1u << c);
        }
    } else {
        int xp = X[pr];
        for (int c = 0; c < P; c++) {
            if (MASK & (1u << c)) continue;
            int w = chid[xp][c];           /* +-1 (c != xp since xp used) */
            X[i] = c; MASK |= (1u << c);
            rec(i + 1, w > 0 ? sign : -sign);
            MASK &= ~(1u << c);
        }
    }
}

typedef struct { const char *name; int npts; int prev[16]; int nedges; } S;

int main(void) {
    static const S structs[] = {
        {"A2",   3, {-1,0,1}, 2},
        {"A4",   5, {-1,0,1,2,3}, 4},
        {"A5",   6, {-1,0,1,2,3,4}, 5},
        {"A6",   7, {-1,0,1,2,3,4,5}, 6},
        {"A8",   9, {-1,0,1,2,3,4,5,6,7}, 8},
        {"J22",  6, {-1,0,1,-1,3,4}, 4},
        {"J23",  7, {-1,0,1,-1,3,4,5}, 5},
        {"J24",  8, {-1,0,1,-1,3,4,5,6}, 6},
        {"J33",  8, {-1,0,1,2,-1,4,5,6}, 6},
        {"J222", 9, {-1,0,1,-1,3,4,-1,6,7}, 6},
    };
    static const int primes[] = {3, 7, 11, 19, 23, 31};
    for (int pi = 0; pi < 6; pi++) {
        P = primes[pi];
        int qr[40]; memset(qr, 0, sizeof qr);
        for (int x = 1; x < P; x++) qr[(x * x) % P] = 1;
        for (int a = 0; a < P; a++)
            for (int b = 0; b < P; b++) {
                int d = ((b - a) % P + P) % P;
                chid[a][b] = (d == 0) ? 0 : (qr[d] ? 1 : -1);
            }
        for (unsigned si = 0; si < sizeof structs / sizeof structs[0]; si++) {
            const S *s = &structs[si];
            if (s->npts > P) { printf("%s %d 0\n", s->name, P); fflush(stdout); continue; }
            NP = s->npts; PREV = s->prev;
            TOTAL = 0;
            long long mult;
            X[0] = 0; MASK = 1u;
            if (s->nedges % 2 == 0) {        /* affine reduction valid */
                X[1] = 1; MASK |= 2u;        /* prev[1]==0, chi(1)=+1 */
                mult = (long long)P * (P - 1);
                rec(2, 1);
            } else {                         /* translation only */
                mult = P;
                rec(1, 1);
            }
            printf("%s %d %lld\n", s->name, P, mult * TOTAL);
            fflush(stdout);
        }
    }
    return 0;
}
