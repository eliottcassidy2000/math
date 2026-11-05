/* aud_blocks.c -- independent audit code (uh_deletion_law_20261104_audit_S10), shares no code with the session.
 *
 * Enumerates every labelled tournament W of order 3..NMAX (default 7) and every A subset of W, keeps the pairs
 * (W, A) that are blocks of type k <= KMAX in the sense of S75's Lemma B4 / THM-4869 B4:
 *   W strong; |A| >= 2; sum_{x in A} end_x(W) = k; every vertex of W - A has out-degree >= 2 in W;
 *   no non-empty Z inside W with (#arcs Z -> W - Z) + |Z & A| = 1,
 * reduces them to isomorphism classes of pairs (W, A) by brute force over all n! relabellings, and prints for each class:
 *   rows, A, scores, end counts, hc(W), |Aut(W, A)|, delta(v) = sum_{x in A - v} end_x(W - v) for every v,
 *   #(delta = 1), #(delta = 2), h0 = #{Hamiltonian paths of W from W - A to A},
 * and whether the class is (AT_m, {m-1, m}) for its order m (AT_m: 1..m, all arcs forward except m -> 1).
 * Usage: aud_blocks [NMAX [KMAX]]   (NMAX <= 7) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef unsigned long long u64;
#define MAXN 8

static int popc(unsigned x) { return __builtin_popcount(x); }

static int is_strong(int n, const unsigned *out) {
    unsigned full = (1u << n) - 1;
    /* forward reach from 0 */
    unsigned seen = 1, frontier = 1;
    while (frontier) { unsigned nf = 0; for (int v = 0; v < n; v++) if ((frontier >> v) & 1) nf |= out[v]; nf &= ~seen; seen |= nf; frontier = nf; }
    if (seen != full) return 0;
    /* backward reach to 0 */
    unsigned in[MAXN]; for (int v = 0; v < n; v++) { in[v] = 0; for (int w = 0; w < n; w++) if ((out[w] >> v) & 1) in[v] |= 1u << w; }
    seen = 1; frontier = 1;
    while (frontier) { unsigned nf = 0; for (int v = 0; v < n; v++) if ((frontier >> v) & 1) nf |= in[v]; nf &= ~seen; seen |= nf; frontier = nf; }
    return seen == full;
}

/* all-pairs Hamiltonian path counts P[s][t] of the tournament on vertex set given by 'mask' (other vertices ignored) */
static void ham_paths(int n, const unsigned *out, unsigned mask, u64 P[MAXN][MAXN]) {
    static u64 dp[1 << MAXN][MAXN];
    for (int s = 0; s < n; s++) for (int t = 0; t < n; t++) P[s][t] = 0;
    for (int s = 0; s < n; s++) {
        if (!((mask >> s) & 1)) continue;
        for (unsigned S = 0; S <= mask; S++) for (int v = 0; v < n; v++) dp[S][v] = 0;
        dp[1u << s][s] = 1;
        for (unsigned S = 1; S <= mask; S++) {
            if ((S & ~mask) || !((S >> s) & 1)) continue;
            for (int v = 0; v < n; v++) {
                u64 c = dp[S][v]; if (!c) continue;
                unsigned nb = out[v] & mask & ~S;
                while (nb) { int w = __builtin_ctz(nb); nb &= nb - 1; dp[S | (1u << w)][w] += c; }
            }
        }
        for (int t = 0; t < n; t++) if ((mask >> t) & 1) P[s][t] = dp[mask][t];
    }
}

static void end_counts(int n, const unsigned *out, unsigned mask, u64 *e) {
    u64 P[MAXN][MAXN]; ham_paths(n, out, mask, P);
    for (int t = 0; t < n; t++) { e[t] = 0; for (int s = 0; s < n; s++) e[t] += P[s][t]; }
    if (popc(mask) == 1) for (int t = 0; t < n; t++) e[t] = ((mask >> t) & 1) ? 1 : 0;
}

static u64 hcycles(int n, const unsigned *out) {
    u64 P[MAXN][MAXN]; ham_paths(n, out, (1u << n) - 1, P);
    u64 h = 0;   /* cycles through vertex 0: paths 0 -> t then t -> 0 */
    for (int t = 1; t < n; t++) if ((out[t] >> 0) & 1) h += P[0][t];
    return h;
}

static int is_block(int n, const unsigned *out, unsigned A) {
    unsigned full = (1u << n) - 1;
    if (popc(A) < 2) return 0;
    for (int v = 0; v < n; v++) if (!((A >> v) & 1) && popc(out[v]) < 2) return 0;
    for (unsigned Z = 1; Z <= full; Z++) {
        int leave = 0;
        for (int v = 0; v < n; v++) if ((Z >> v) & 1) leave += popc(out[v] & ~Z & full);
        if (leave + popc(Z & A) == 1) return 0;
    }
    return 1;
}

typedef struct { u64 code; unsigned A; } Key;

static Key key_of(int n, const unsigned *out, unsigned A, const int *perm) {
    /* perm[old] = new label */
    Key k; k.code = 0; k.A = 0;
    unsigned nr[MAXN]; for (int i = 0; i < n; i++) nr[i] = 0;
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if ((out[i] >> j) & 1) nr[perm[i]] |= 1u << perm[j];
    for (int i = 0; i < n; i++) k.code = (k.code << n) | nr[i];
    for (int i = 0; i < n; i++) if ((A >> i) & 1) k.A |= 1u << perm[i];
    return k;
}
static int keycmp(Key a, Key b) { if (a.code != b.code) return a.code < b.code ? -1 : 1; if (a.A != b.A) return a.A < b.A ? -1 : 1; return 0; }

static int nperm; static int perms[5040][MAXN];
static void gen_perms(int n) {
    int p[MAXN]; for (int i = 0; i < n; i++) p[i] = i; nperm = 0;
    for (;;) {
        memcpy(perms[nperm++], p, sizeof(int) * n);
        int i = n - 2; while (i >= 0 && p[i] > p[i + 1]) i--; if (i < 0) break;
        int j = n - 1; while (p[j] < p[i]) j--; int t = p[i]; p[i] = p[j]; p[j] = t;
        for (int a = i + 1, b = n - 1; a < b; a++, b--) { t = p[a]; p[a] = p[b]; p[b] = t; }
    }
}
static Key canon_key(int n, const unsigned *out, unsigned A, int *aut) {
    Key best = key_of(n, out, A, perms[0]); int cnt = 1;
    for (int q = 1; q < nperm; q++) {
        Key k = key_of(n, out, A, perms[q]); int c = keycmp(k, best);
        if (c < 0) { best = k; cnt = 1; } else if (c == 0) cnt++;
    }
    if (aut) *aut = cnt;
    return best;
}

typedef struct { int n; unsigned out[MAXN]; unsigned A; Key key; int k; } Block;
static Block found[100000]; static int nfound = 0;

int main(int argc, char **argv) {
    int NMAX = argc > 1 ? atoi(argv[1]) : 7, KMAX = argc > 2 ? atoi(argv[2]) : 5;
    if (NMAX > 7) NMAX = 7;
    printf("aud_blocks: blocks (W, A) of type k <= %d, orders 3..%d, all labelled tournaments\n", KMAX, NMAX);
    long labelled_blocks_total = 0;
    for (int n = 3; n <= NMAX; n++) {
        gen_perms(n);
        int np = n * (n - 1) / 2; int pi[64], pj[64], c = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { pi[c] = i; pj[c] = j; c++; }
        long nstrong = 0, nlab = 0; int before = nfound;
        for (u64 bits = 0; bits < (1ull << np); bits++) {
            unsigned out[MAXN]; for (int v = 0; v < n; v++) out[v] = 0;
            for (int b = 0; b < np; b++) { if ((bits >> b) & 1) out[pi[b]] |= 1u << pj[b]; else out[pj[b]] |= 1u << pi[b]; }
            if (!is_strong(n, out)) continue;
            nstrong++;
            u64 e[MAXN]; end_counts(n, out, (1u << n) - 1, e);
            /* two smallest end counts */
            u64 m1 = ~0ull, m2 = ~0ull; for (int v = 0; v < n; v++) { if (e[v] < m1) { m2 = m1; m1 = e[v]; } else if (e[v] < m2) m2 = e[v]; }
            if (m1 + m2 > (u64)KMAX) continue;
            for (unsigned A = 0; A < (1u << n); A++) {
                if (popc(A) < 2) continue;
                u64 k = 0; for (int v = 0; v < n; v++) if ((A >> v) & 1) k += e[v];
                if (k > (u64)KMAX) continue;
                if (!is_block(n, out, A)) continue;
                nlab++;
                Key key = canon_key(n, out, A, NULL);
                int f = -1; for (int q = before; q < nfound; q++) if (!keycmp(found[q].key, key)) { f = q; break; }
                if (f < 0) { Block *B = &found[nfound++]; B->n = n; memcpy(B->out, out, sizeof(out)); B->A = A; B->key = key; B->k = (int)k; }
            }
        }
        labelled_blocks_total += nlab;
        printf("order %d: %ld strong labelled tournaments, %ld labelled blocks of type <= %d, %d isomorphism classes\n",
               n, nstrong, nlab, KMAX, nfound - before);
    }
    /* AT_m keys */
    for (int k = 2; k <= KMAX; k++) {
        int cnt = 0; for (int q = 0; q < nfound; q++) if (found[q].k == k) cnt++;
        printf("\ntype %d: %d blocks\n", k, cnt);
        for (int q = 0; q < nfound; q++) {
            Block *B = &found[q]; if (B->k != k) continue;
            int n = B->n; gen_perms(n);
            int aut; canon_key(n, B->out, B->A, &aut);
            /* AT_n with A = {n-2, n-1} (0-indexed) */
            unsigned at[MAXN]; for (int i = 0; i < n; i++) at[i] = 0;
            for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (i == 0 && j == n - 1) at[n - 1] |= 1u; else at[i] |= 1u << j; }
            Key atk = canon_key(n, at, (1u << (n - 2)) | (1u << (n - 1)), NULL);
            int isAT = !keycmp(atk, B->key);
            u64 e[MAXN]; end_counts(n, B->out, (1u << n) - 1, e);
            u64 P[MAXN][MAXN]; ham_paths(n, B->out, (1u << n) - 1, P);
            u64 h0 = 0; for (int s = 0; s < n; s++) for (int t = 0; t < n; t++) if (s != t && !((B->A >> s) & 1) && ((B->A >> t) & 1)) h0 += P[s][t];
            u64 hcw = hcycles(n, B->out);
            int delta[MAXN], d1 = 0, d2 = 0;
            for (int v = 0; v < n; v++) {
                unsigned mask = ((1u << n) - 1) & ~(1u << v);
                u64 e2[MAXN]; end_counts(n, B->out, mask, e2);
                u64 s = 0; for (int x = 0; x < n; x++) if (x != v && ((B->A >> x) & 1)) s += e2[x];
                delta[v] = (int)s; if (s == 1) d1++; if (s == 2) d2++;
            }
            printf("  order %d, rows(out-nbr lists):", n);
            for (int v = 0; v < n; v++) { printf(" %d:{", v); int first = 1; for (int w = 0; w < n; w++) if ((B->out[v] >> w) & 1) { printf(first ? "%d" : ",%d", w); first = 0; } printf("}"); }
            printf("\n    A = {"); { int first = 1; for (int v = 0; v < n; v++) if ((B->A >> v) & 1) { printf(first ? "%d" : ",%d", v); first = 0; } } printf("}");
            printf(" (A = W: %s); scores", B->A == (1u << n) - 1 ? "yes" : "no");
            for (int v = 0; v < n; v++) printf(" %d", popc(B->out[v]));
            printf("; end counts"); for (int v = 0; v < n; v++) printf(" %llu", e[v]);
            printf("; hc(W) = %llu; |Aut(W,A)| = %d; is (AT_%d, {%d,%d}): %s\n", hcw, aut, n, n - 1, n, isAT ? "YES" : "no");
            printf("    delta(v), v = 0..%d:", n - 1); for (int v = 0; v < n; v++) printf(" %d", delta[v]);
            printf("   #(delta=1) = %d, #(delta=2) = %d, h0 = %llu\n", d1, d2, h0);
        }
    }
    return 0;
}
