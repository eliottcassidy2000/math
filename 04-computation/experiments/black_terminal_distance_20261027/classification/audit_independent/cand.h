/* cand.h -- exact enumeration of descent candidates (auditor's own code).
 * A descent candidate is a Hamiltonian path Q with q_1 in X, q_n notin X for some minimum out-cut X.
 * With sep[a] = union over minimum cuts X containing a of (V \ X), Q is a candidate iff q_n in sep[q_1].
 * For each distinct end set E = sep[a], a subset DP  R[mask] = {v in mask : T[mask] has a Ham path from v ending in E}
 * prunes a DFS so that every leaf is a candidate; each candidate is produced exactly once (its start is unique).
 */
#ifndef CAND_H
#define CAND_H
#include "at.h"

typedef int (*cand_cb)(const int *Q, void *ctx);   /* return nonzero to stop */

static uint32_t *cand_ws = NULL; static int cand_wsn = -1;
static void cand_ws_alloc(int n) {
    if (cand_wsn == n) return;
    free(cand_ws); cand_ws = (uint32_t *)malloc(sizeof(uint32_t) << n); cand_wsn = n;
    if (!cand_ws) { fprintf(stderr, "oom\n"); exit(1); }
}
static void cand_dp(const tour *T, uint32_t E, uint32_t *R) {
    int n = T->n; uint32_t N = 1u << n;
    R[0] = 0;
    for (uint32_t mask = 1; mask < N; mask++) {
        uint32_t r = 0;
        if ((mask & (mask - 1)) == 0) { r = mask & E; }
        else for (uint32_t c = mask; c; c &= c - 1) { int v = lowbit(c); if (T->o[v] & R[mask ^ (1u << v)]) r |= 1u << v; }
        R[mask] = r;
    }
}
static int cand_dfs(const tour *T, const uint32_t *R, int *Q, int d, uint32_t mask, cand_cb cb, void *ctx) {
    /* Q[d-1] is current vertex v, mask = unvisited vertices including v */
    int v = Q[d - 1];
    if (mask == (1u << v)) return cb(Q, ctx);
    uint32_t rest = mask ^ (1u << v);
    for (uint32_t c = T->o[v] & R[rest]; c; c &= c - 1) {
        int w = lowbit(c); Q[d] = w;
        if (cand_dfs(T, R, Q, d + 1, rest, cb, ctx)) return 1;
    }
    return 0;
}
/* returns 1 if stopped by callback */
static int cand_enum(const tour *T, const uint32_t *sep, cand_cb cb, void *ctx) {
    int n = T->n; cand_ws_alloc(n);
    uint32_t done = 0; uint32_t F = fullmask(n); int Q[MAXV];
    for (int a = 0; a < n; a++) {
        if (done >> a & 1u) continue;
        uint32_t E = sep[a]; if (!E) { done |= 1u << a; continue; }
        cand_dp(T, E, cand_ws);
        for (int b = a; b < n; b++) if (!(done >> b & 1u) && sep[b] == E) {
            done |= 1u << b;
            if (cand_ws[F] >> b & 1u) { Q[0] = b; if (cand_dfs(T, cand_ws, Q, 1, F, cb, ctx)) return 1; }
        }
    }
    return 0;
}
static void make_sep(const tour *T, const uint32_t *cuts, int nc, uint32_t *sep) {
    uint32_t F = fullmask(T->n);
    for (int a = 0; a < T->n; a++) sep[a] = 0;
    for (int i = 0; i < nc; i++) for (int a = 0; a < T->n; a++) if (cuts[i] >> a & 1u) sep[a] |= F & ~cuts[i];
}

/* black-descent search: returns 1 and fills Qout if a black candidate exists */
typedef struct { const tour *T; int found; int Q[MAXV]; long long ncand, nblue; int single; int lambda; int needSingle; } bd_ctx;
static int bd_cb(const int *Q, void *vc) {
    bd_ctx *c = (bd_ctx *)vc; c->ncand++;
    if (t_blue(c->T, Q)) { c->nblue++; return 0; }
    int n = c->T->n;
    int single = (c->T->s[Q[0]] == c->lambda) || (n - 1 - c->T->s[Q[n - 1]] == c->lambda);
    if (c->needSingle && !single) return 0;
    c->found = 1; c->single = single; memcpy(c->Q, Q, sizeof(int) * n);
    return 1;
}
/* lambda and min cuts computed inside; returns found flag */
static int black_descent(const tour *T, int needSingle, int *Qout, int *lambdaOut, long long *ncand, long long *nblue) {
    static uint32_t cuts[1 << 16]; int nc;
    int lam = t_mincuts(T, cuts, 1 << 16, &nc);
    if (nc > (1 << 16)) { fprintf(stderr, "too many cuts\n"); exit(1); }
    uint32_t sep[MAXV]; make_sep(T, cuts, nc, sep);
    bd_ctx c; memset(&c, 0, sizeof c); c.T = T; c.lambda = lam; c.needSingle = needSingle;
    cand_enum(T, sep, bd_cb, &c);
    if (lambdaOut) *lambdaOut = lam;
    if (ncand) *ncand = c.ncand;
    if (nblue) *nblue = c.nblue;
    if (c.found && Qout) memcpy(Qout, c.Q, sizeof(int) * T->n);
    return c.found;
}
#endif
