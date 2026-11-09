/* AUDIT2 (independent of egcp_lib.py / e0search.c): exhaustive census of ALL labelled graphs on exactly n vertices.
 * For every graph: d = number of dominating sets (subset DP), universal-vertex flag, connectivity, and the FULL
 * automorphism group by backtracking, with the von Broemssen sign of every automorphism (parity of the number of
 * edges u<v whose images satisfy p(u) > p(v)).  The graph is odd iff some automorphism has odd sign.
 * Validation built in: sum over labelled graphs of |Aut|/n! = number of isomorphism classes (A000088), and the same sum
 * over even graphs must equal the number of tournaments (A000568, RPGFD 2023).
 * Output (stdout): per value k: labelled counts and iso-class counts (as sum |Aut|/n!, printed as exact rationals
 * via the integer sum of |Aut|) for odd graphs, even graphs, even graphs without universal vertex, and connected
 * odd/even graphs; plus one witness mask per (k, class).
 * usage: a2_census n            (n <= 8)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#define MAXN 8
#define MAXV 256

typedef struct {
    long long aut_odd[MAXV], aut_even[MAXV], aut_even0[MAXV], aut_codd[MAXV], aut_ceven[MAXV];
    long long lab_odd[MAXV], lab_even[MAXV], lab_even0[MAXV], lab_codd[MAXV], lab_ceven[MAXV];
    long long wit_odd[MAXV], wit_even[MAXV], wit_even0[MAXV];
    long long aut_all, aut_alleven, aut_conn_even, aut_conn;
} acc_t;

static int N, M;
static int PU[64], PV[64];

typedef struct {
    int adj[MAXN], deg[MAXN], ord[MAXN], perm[MAXN], used[MAXN];
    long long naut; int nodd;
} st_t;

static void bt(st_t *s, int i) {
    if (i == N) {
        int rev = 0;
        for (int e = 0; e < M; e++) {
            int u = PU[e], v = PV[e];
            if (((s->adj[u] >> v) & 1) && s->perm[u] > s->perm[v]) rev++;
        }
        s->naut++;
        if (rev & 1) s->nodd++;
        return;
    }
    int v = s->ord[i];
    for (int w = 0; w < N; w++) {
        if (s->used[w] || s->deg[w] != s->deg[v]) continue;
        int ok = 1;
        for (int j = 0; j < i; j++) {
            int u = s->ord[j];
            if (((s->adj[v] >> u) & 1) != ((s->adj[w] >> s->perm[u]) & 1)) { ok = 0; break; }
        }
        if (!ok) continue;
        s->perm[v] = w; s->used[w] = 1;
        bt(s, i + 1);
        s->used[w] = 0; s->perm[v] = -1;
    }
}

static int connected(const int *adj) {
    int seen = 1, frontier = 1;
    while (frontier) {
        int nf = 0;
        for (int v = 0; v < N; v++) if ((frontier >> v) & 1) nf |= adj[v];
        nf &= ~seen; seen |= nf; frontier = nf;
    }
    return seen == (1 << N) - 1;
}

int main(int argc, char **argv) {
    N = atoi(argv[1]);
    M = 0;
    for (int u = 0; u < N; u++) for (int v = u + 1; v < N; v++) { PU[M] = u; PV[M] = v; M++; }
    long long total = 1LL << M;
    int full = (1 << N) - 1;
    int nth = omp_get_max_threads();
    acc_t *acc = calloc(nth, sizeof(acc_t));
    for (int t = 0; t < nth; t++) for (int k = 0; k < MAXV; k++) acc[t].wit_odd[k] = acc[t].wit_even[k] = acc[t].wit_even0[k] = -1;
#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        acc_t *a = &acc[tid];
        st_t s;
        int dom[1 << MAXN];
#pragma omp for schedule(dynamic, 4096)
        for (long long mask = 0; mask < total; mask++) {
            memset(s.adj, 0, sizeof(s.adj));
            for (int e = 0; e < M; e++) if ((mask >> e) & 1) { s.adj[PU[e]] |= 1 << PV[e]; s.adj[PV[e]] |= 1 << PU[e]; }
            /* dominating sets by subset DP */
            int d = 0; dom[0] = 0;
            for (int S = 1; S <= full; S++) {
                int v = __builtin_ctz(S);
                dom[S] = dom[S & (S - 1)] | s.adj[v] | (1 << v);
                if (dom[S] == full) d++;
            }
            int univ = 0;
            for (int v = 0; v < N; v++) { s.deg[v] = __builtin_popcount(s.adj[v]); if (s.adj[v] == (full ^ (1 << v))) univ = 1; }
            if (N == 1) univ = 1; /* K_1 is its own universal vertex */
            /* vertex order: by degree descending, ties by index */
            for (int v = 0; v < N; v++) s.ord[v] = v;
            for (int x = 0; x < N; x++) for (int y = x + 1; y < N; y++)
                if (s.deg[s.ord[y]] > s.deg[s.ord[x]]) { int t = s.ord[x]; s.ord[x] = s.ord[y]; s.ord[y] = t; }
            for (int v = 0; v < N; v++) { s.perm[v] = -1; s.used[v] = 0; }
            s.naut = 0; s.nodd = 0;
            bt(&s, 0);
            int odd = s.nodd > 0;
            int conn = connected(s.adj);
            a->aut_all += s.naut;
            if (conn) a->aut_conn += s.naut;
            if (odd) {
                a->aut_odd[d] += s.naut; a->lab_odd[d]++;
                if (a->wit_odd[d] < 0) a->wit_odd[d] = mask;
                if (conn) { a->aut_codd[d] += s.naut; a->lab_codd[d]++; }
            } else {
                a->aut_alleven += s.naut;
                a->aut_even[d] += s.naut; a->lab_even[d]++;
                if (a->wit_even[d] < 0) a->wit_even[d] = mask;
                if (!univ) { a->aut_even0[d] += s.naut; a->lab_even0[d]++; if (a->wit_even0[d] < 0) a->wit_even0[d] = mask; }
                if (conn) { a->aut_ceven[d] += s.naut; a->lab_ceven[d]++; a->aut_conn_even += s.naut; }
            }
        }
    }
    /* merge */
    acc_t R; memset(&R, 0, sizeof(R));
    for (int k = 0; k < MAXV; k++) R.wit_odd[k] = R.wit_even[k] = R.wit_even0[k] = -1;
    for (int t = 0; t < nth; t++) {
        R.aut_all += acc[t].aut_all; R.aut_alleven += acc[t].aut_alleven; R.aut_conn_even += acc[t].aut_conn_even; R.aut_conn += acc[t].aut_conn;
        for (int k = 0; k < MAXV; k++) {
            R.aut_odd[k] += acc[t].aut_odd[k]; R.aut_even[k] += acc[t].aut_even[k]; R.aut_even0[k] += acc[t].aut_even0[k];
            R.aut_codd[k] += acc[t].aut_codd[k]; R.aut_ceven[k] += acc[t].aut_ceven[k];
            R.lab_odd[k] += acc[t].lab_odd[k]; R.lab_even[k] += acc[t].lab_even[k]; R.lab_even0[k] += acc[t].lab_even0[k];
            R.lab_codd[k] += acc[t].lab_codd[k]; R.lab_ceven[k] += acc[t].lab_ceven[k];
            if (acc[t].wit_odd[k] >= 0 && (R.wit_odd[k] < 0 || acc[t].wit_odd[k] < R.wit_odd[k])) R.wit_odd[k] = acc[t].wit_odd[k];
            if (acc[t].wit_even[k] >= 0 && (R.wit_even[k] < 0 || acc[t].wit_even[k] < R.wit_even[k])) R.wit_even[k] = acc[t].wit_even[k];
            if (acc[t].wit_even0[k] >= 0 && (R.wit_even0[k] < 0 || acc[t].wit_even0[k] < R.wit_even0[k])) R.wit_even0[k] = acc[t].wit_even0[k];
        }
    }
    long long nf = 1; for (int i = 2; i <= N; i++) nf *= i;
    printf("N %d labelled %lld nfact %lld\n", N, total, nf);
    printf("ISO all %lld (sumAut %lld) even %lld (sumAut %lld) connected %lld connected_even %lld\n",
           R.aut_all / nf, R.aut_all, R.aut_alleven / nf, R.aut_alleven, R.aut_conn / nf, R.aut_conn_even / nf);
    printf("# k lab_odd iso_odd lab_even iso_even lab_even0 iso_even0 lab_codd iso_codd lab_ceven iso_ceven wit_odd wit_even wit_even0\n");
    for (int k = 0; k <= full; k++) {
        if (!(R.lab_odd[k] || R.lab_even[k])) continue;
        printf("V %d %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld\n", k,
               R.lab_odd[k], R.aut_odd[k] / nf, R.lab_even[k], R.aut_even[k] / nf, R.lab_even0[k], R.aut_even0[k] / nf,
               R.lab_codd[k], R.aut_codd[k] / nf, R.lab_ceven[k], R.aut_ceven[k] / nf,
               R.wit_odd[k], R.wit_even[k], R.wit_even0[k]);
        if (R.aut_odd[k] % nf || R.aut_even[k] % nf || R.aut_even0[k] % nf || R.aut_codd[k] % nf || R.aut_ceven[k] % nf)
            printf("WARNING non-integral iso count at k=%d\n", k);
    }
    return 0;
}
