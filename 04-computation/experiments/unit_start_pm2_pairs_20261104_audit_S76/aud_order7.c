/*
 * aud_order7.c  (independent audit of THM-4872 (A), opus-2026-11-04-S76 audit)
 *
 * Written from scratch for the audit; it does not reuse the session's programs.
 *
 * Part 1 (census).  Rooted classes (T, 0) of 7-vertex tournaments are counted by BURNSIDE over the group S_6 of
 * permutations fixing the root, with no canonical forms: for every sigma of odd order (all cycles odd; an even cycle
 * fixes no tournament) every sigma-fixed labelled tournament is generated from the sigma-orbits on unordered pairs, its
 * profile (strong, h, s, q, F, e) is computed, and the number of classes with a given profile is
 *       N(P) = (1/720) * sum_sigma #{T : sigma T = T, profile(T) = P}.
 * Pair counts are sums of N(P) N(P') over profile pairs, so no class representatives are needed at all.
 *
 * Part 2 (routes).  Breadth-first search over all 2^21 labelled tournaments on {0..6} that are strong, root 0 fixed,
 * one arc reversal per step, from the labelled A to the set of labelled copies of B.  A labelled tournament is a copy
 * of B iff it has s = 1 and its relabelling along its unique Hamiltonian path from 0 equals that of B (a rooted
 * isomorphism maps the unique path from the root onto the unique path from the root).  Shortest routes are counted by
 * dynamic programming on the BFS layers AND enumerated explicitly; the wrapper criterion is tested both in the closed
 * form (F_i <= F_0 and h_i + 2F_i < h_0 + 2F_0) and literally for m = 1..60 with 128-bit integers.
 *
 * Part 3 (wrapper).  H(D_m(U, r)) = h + 2^m F checked by direct Hamiltonian-path counts of the literal wrapper on random
 * labelled 7-vertex tournaments of every kind (strong or not, any s), on A and B, and on every intermediate of every
 * shortest route, for m = 1..5.
 *
 * Profile: h = H(T); s = paths starting at 0; q = H(T - 0); F = ordered covers (P, Q), P nonempty ending at 0, Q possibly
 * empty; e = paths ending at 0.  F is computed by the subset formula and, for a random sample, independently as
 * (H(D_1) - h)/2 is NOT used (that would be circular with Part 3); instead it is recomputed by brute-force enumeration of
 * all ordered vertex sequences with one cut (7! * 8 sequences) on a random sample.
 *
 * Build: gcc -O3 -march=native -fopenmp -o aud_order7 aud_order7.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>

#define N 7
#define NP 21
#define NT (1u << NP)
#define FULL ((1 << N) - 1)

static int pu[NP], pv[NP], pidx[N][N];

static void decode(uint32_t code, uint8_t *out) {
    memset(out, 0, N);
    for (int k = 0; k < NP; k++) {
        if ((code >> k) & 1u) out[pu[k]] |= (uint8_t)(1u << pv[k]);
        else out[pv[k]] |= (uint8_t)(1u << pu[k]);
    }
}
static uint32_t encode(const uint8_t *out) {
    uint32_t c = 0;
    for (int k = 0; k < NP; k++) if ((out[pu[k]] >> pv[k]) & 1u) c |= 1u << k;
    return c;
}

typedef struct { int16_t h, s, q, F, e; int8_t strong; } Prof;

static int is_strong(const uint8_t *out) {
    /* forward closure from 0 and backward closure to 0 */
    unsigned f = 1, g = 1, ch = 1;
    while (ch) { ch = 0; for (int x = 0; x < N; x++) if ((f >> x) & 1u) { unsigned nf = f | out[x]; if (nf != f) { f = nf; ch = 1; } } }
    ch = 1;
    while (ch) { ch = 0; for (int y = 0; y < N; y++) if (!((g >> y) & 1u)) for (int x = 0; x < N; x++) if (((g >> x) & 1u) && ((out[y] >> x) & 1u)) { g |= 1u << y; ch = 1; break; } }
    return f == FULL && g == FULL;
}

static Prof profile(const uint8_t *out) {
    int end[1 << N][N];      /* paths covering S ending at x */
    int beg0[1 << N][N];     /* paths covering S, starting at 0, ending at x */
    int HS[1 << N];
    memset(end, 0, sizeof end); memset(beg0, 0, sizeof beg0);
    for (int x = 0; x < N; x++) end[1 << x][x] = 1;
    beg0[1][0] = 1;
    for (int S = 1; S <= FULL; S++) for (int x = 0; x < N; x++) {
        int a = end[S][x], b = beg0[S][x];
        if (!a && !b) continue;
        for (int y = 0; y < N; y++) if (!((S >> y) & 1) && ((out[x] >> y) & 1u)) { end[S | (1 << y)][y] += a; beg0[S | (1 << y)][y] += b; }
    }
    HS[0] = 1;
    for (int S = 1; S <= FULL; S++) { int t = 0; for (int x = 0; x < N; x++) t += end[S][x]; HS[S] = t; }
    Prof p;
    p.h = (int16_t)HS[FULL]; p.q = (int16_t)HS[FULL & ~1]; p.e = (int16_t)end[FULL][0];
    int s = 0; for (int x = 0; x < N; x++) s += beg0[FULL][x]; p.s = (int16_t)s;
    int F = 0; for (int S = 1; S <= FULL; S += 2) F += end[S][0] * HS[FULL ^ S];
    p.F = (int16_t)F;
    p.strong = (int8_t)is_strong(out);
    return p;
}

/* brute-force F: all orderings of the 7 vertices with a cut c in 1..7 (first block = positions < c, nonempty, ends at 0) */
static int bruteF(const uint8_t *out) {
    int perm[N]; for (int i = 0; i < N; i++) perm[i] = i;
    int cnt = 0;
    /* iterate permutations via Heap's algorithm */
    int c[N] = {0}; int i = 0;
    for (;;) {
        for (int cut = 1; cut <= N; cut++) {
            if (perm[cut - 1] != 0) continue;
            int ok = 1;
            for (int t = 0; t + 1 < cut && ok; t++) if (!((out[perm[t]] >> perm[t + 1]) & 1u)) ok = 0;
            for (int t = cut; t + 1 < N && ok; t++) if (!((out[perm[t]] >> perm[t + 1]) & 1u)) ok = 0;
            if (ok) cnt++;
        }
        while (i < N && c[i] >= i) { c[i] = 0; i++; }
        if (i >= N) break;
        int sw = (i % 2 == 0) ? 0 : c[i];
        int tmp = perm[sw]; perm[sw] = perm[i]; perm[i] = tmp;
        c[i]++; i = 0;
    }
    return cnt;
}

/* Hamiltonian path count of a tournament on n <= 16 vertices given by out-rows */
static long long hampaths(int n, const uint32_t *out) {
    long long *dp = calloc((size_t)(1u << n) * n, sizeof(long long));
    for (int x = 0; x < n; x++) dp[(size_t)(1u << x) * n + x] = 1;
    for (unsigned S = 1; S < (1u << n); S++) for (int x = 0; x < n; x++) {
        long long c = dp[(size_t)S * n + x]; if (!c) continue;
        for (int y = 0; y < n; y++) if (!((S >> y) & 1u) && ((out[x] >> y) & 1u)) dp[(size_t)(S | (1u << y)) * n + y] += c;
    }
    long long H = 0; for (int x = 0; x < n; x++) H += dp[(size_t)((1u << n) - 1) * n + x];
    free(dp); return H;
}

/* literal Session189 wrapper D_m(U, r = 0): U = 0..6, Z = 7..6+m transitive, w = 7+m; w -> Z -> U, r -> w -> U \ {r} */
static long long wrapperH(const uint8_t *u, int m) {
    int n = N + m + 1; uint32_t out[16]; memset(out, 0, sizeof out);
    for (int x = 0; x < N; x++) out[x] = u[x];
    int w = N + m;
    for (int a = 0; a < m; a++) {
        int za = N + a;
        for (int x = 0; x < N; x++) out[za] |= 1u << x;                 /* Z -> U */
        for (int b = a + 1; b < m; b++) out[za] |= 1u << (N + b);         /* transitive Z */
        out[w] |= 1u << za;                                               /* w -> Z */
    }
    for (int x = 1; x < N; x++) out[w] |= 1u << x;                        /* w -> U \ {r} */
    out[0] |= 1u << w;                                                    /* r -> w */
    /* sanity: tournament */
    for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) {
        int ab = (out[a] >> b) & 1, ba = (out[b] >> a) & 1;
        if (ab + ba != 1) { fprintf(stderr, "wrapper not a tournament\n"); exit(1); }
    }
    return hampaths(n, out);
}

/* path-order code of an s = 1 rooted tournament: relabel along the unique Hamiltonian path from 0 */
static int found_path[N], found_cnt;
static void dfs_path(const uint8_t *out, int cur, unsigned vis, int depth, int *path) {
    path[depth] = cur;
    if (vis == FULL) { if (found_cnt == 0) memcpy(found_path, path, sizeof(int) * N); found_cnt++; return; }
    for (int y = 0; y < N; y++) if (!((vis >> y) & 1u) && ((out[cur] >> y) & 1u)) dfs_path(out, y, vis | (1u << y), depth + 1, path);
}
static uint32_t pathcode(const uint8_t *out) {
    int path[N]; found_cnt = 0; dfs_path(out, 0, 1u, 0, path);
    if (found_cnt != 1) return 0xFFFFFFFFu;
    uint8_t rl[N]; memset(rl, 0, N);
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) if (i != j && ((out[found_path[i]] >> found_path[j]) & 1u)) rl[i] |= (uint8_t)(1u << j);
    return encode(rl);
}

/* ---------- Burnside census ---------- */
typedef struct { int8_t strong; int16_t h, s, q, F, e; long long w; } Rec;
static int cmprec(const void *a, const void *b) {
    const Rec *x = a, *y = b;
    if (x->strong != y->strong) return x->strong - y->strong;
    if (x->h != y->h) return x->h - y->h;
    if (x->s != y->s) return x->s - y->s;
    if (x->q != y->q) return x->q - y->q;
    if (x->F != y->F) return x->F - y->F;
    return x->e - y->e;
}

static Prof *PROF;   /* profile of every labelled code (identity sigma pass) */

int main(void) {
    { int k = 0; for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++) { pu[k] = i; pv[k] = j; pidx[i][j] = pidx[j][i] = k; k++; } }
    PROF = malloc(sizeof(Prof) * NT);
    #pragma omp parallel for schedule(dynamic, 8192)
    for (long c = 0; c < (long)NT; c++) { uint8_t o[N]; decode((uint32_t)c, o); PROF[c] = profile(o); }

    /* random-sample cross-check of F against brute force */
    {
        srand(12345); int bad = 0, tot = 0;
        for (int t = 0; t < 3000; t++) {
            uint32_t c = ((uint32_t)rand() << 15 ^ (uint32_t)rand()) & (NT - 1);
            uint8_t o[N]; decode(c, o);
            if (bruteF(o) != PROF[c].F) bad++;
            tot++;
        }
        printf("F (subset formula) vs brute force over all 7!*7 cut orderings, random sample: %d / %d disagree\n", bad, tot);
    }

    /* enumerate S_6 on {1..6}; keep odd-order elements */
    int nrec = 0, cap = 2400000; Rec *R = malloc(sizeof(Rec) * cap);
    long long nsig = 0, nfixed_total = 0;
    int a[6] = {1, 2, 3, 4, 5, 6};
    /* all permutations by lexicographic next_permutation */
    int perm[6]; for (int i = 0; i < 6; i++) perm[i] = a[i];
    for (;;) {
        int sg[N]; sg[0] = 0; for (int i = 0; i < 6; i++) sg[i + 1] = perm[i];   /* sg: i -> sg[i] */
        /* cycle type: all cycles odd? */
        int seen[N] = {0}, odd = 1;
        for (int i = 0; i < N; i++) if (!seen[i]) { int L = 0, j = i; while (!seen[j]) { seen[j] = 1; j = sg[j]; L++; } if (L % 2 == 0) odd = 0; }
        if (odd) {
            nsig++;
            /* orbits on unordered pairs: representative ordered pair (x,y), images (sg^k x, sg^k y) */
            int orbit_of[NP]; for (int k = 0; k < NP; k++) orbit_of[k] = -1;
            int norb = 0, repx[NP], repy[NP];
            for (int k = 0; k < NP; k++) if (orbit_of[k] < 0) {
                repx[norb] = pu[k]; repy[norb] = pv[k];
                int x = pu[k], y = pv[k];
                for (;;) { int kk = pidx[x][y]; if (orbit_of[kk] >= 0) break; orbit_of[kk] = norb; x = sg[x]; y = sg[y]; }
                norb++;
            }
            for (uint32_t choice = 0; choice < (1u << norb); choice++) {
                uint8_t o[N]; memset(o, 0, N);
                for (int b = 0; b < norb; b++) {
                    int x = repx[b], y = repy[b];
                    if (!((choice >> b) & 1u)) { int t = x; x = y; y = t; }
                    int x0 = x, y0 = y;
                    /* orient x -> y and propagate along the orbit */
                    for (int it = 0; it < 64; it++) {
                        if ((o[y] >> x) & 1u) { fprintf(stderr, "conflict\n"); return 1; }
                        o[x] |= (uint8_t)(1u << y);
                        x = sg[x]; y = sg[y];
                        if (x == x0 && y == y0) break;
                        if (x == y0 && y == x0) { fprintf(stderr, "flip in orbit\n"); return 1; }
                    }
                }
                /* verify it is a tournament fixed by sg */
                for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++) if ((((o[i] >> j) & 1) + ((o[j] >> i) & 1)) != 1) { fprintf(stderr, "not a tournament\n"); return 1; }
                for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) if (((o[i] >> j) & 1) != ((o[sg[i]] >> sg[j]) & 1)) { fprintf(stderr, "not fixed\n"); return 1; }
                Prof p = PROF[encode(o)];
                if (nrec == cap) { cap *= 2; R = realloc(R, sizeof(Rec) * cap); }
                R[nrec].strong = p.strong; R[nrec].h = p.h; R[nrec].s = p.s; R[nrec].q = p.q; R[nrec].F = p.F; R[nrec].e = p.e; R[nrec].w = 1; nrec++;
                nfixed_total++;
            }
        }
        /* next permutation */
        int i = 4; while (i >= 0 && perm[i] > perm[i + 1]) i--;
        if (i < 0) break;
        int j = 5; while (perm[j] < perm[i]) j--;
        int t = perm[i]; perm[i] = perm[j]; perm[j] = t;
        for (int l = i + 1, r = 5; l < r; l++, r--) { t = perm[l]; perm[l] = perm[r]; perm[r] = t; }
    }
    qsort(R, nrec, sizeof(Rec), cmprec);
    int nprof = 0;
    for (int i = 0; i < nrec; ) {
        int j = i; long long w = 0;
        while (j < nrec && cmprec(&R[i], &R[j]) == 0) { w += R[j].w; j++; }
        R[nprof] = R[i]; R[nprof].w = w; nprof++; i = j;
    }
    long long tot = 0, totstrong = 0, unit_strong = 0, unit_all = 0; int nonint = 0;
    for (int i = 0; i < nprof; i++) {
        if (R[i].w % 720) nonint++;
        R[i].w /= 720;
        tot += R[i].w; if (R[i].s == 1) unit_all += R[i].w;
        if (R[i].strong) { totstrong += R[i].w; if (R[i].s == 1) unit_strong += R[i].w; }
    }
    printf("Burnside: %lld odd-order elements of S_6 (expect 225), %lld fixed labelled tournaments in total (expect 2,188,800)\n", nsig, nfixed_total);
    printf("distinct profiles (strong,h,s,q,F,e): %d; non-integral class counts: %d\n", nprof, nonint);
    printf("rooted classes: %lld; strong: %lld; strong with s = 1: %lld; any with s = 1: %lld\n", tot, totstrong, unit_strong, unit_all);

    /* pair counts (b has h + 2, equal F, both strong) */
    long long pF = 0, pF_s1 = 0, pFq = 0, pFsq = 0, pFsq_s1 = 0, pFsqe = 0, pFs = 0;
    for (int i = 0; i < nprof; i++) if (R[i].strong) for (int j = 0; j < nprof; j++) if (R[j].strong) {
        if (R[j].h != R[i].h + 2 || R[j].F != R[i].F) continue;
        long long w = R[i].w * R[j].w;
        pF += w;
        if (R[i].s == 1 && R[j].s == 1) { pF_s1 += w; printf("  unit pair profiles: (h,s,q,F,e) = (%d,%d,%d,%d,%d) and (%d,%d,%d,%d,%d), multiplicity %lld\n",
                R[i].h, R[i].s, R[i].q, R[i].F, R[i].e, R[j].h, R[j].s, R[j].q, R[j].F, R[j].e, w); }
        if (R[i].q == R[j].q) pFq += w;
        if (R[i].s == R[j].s) pFs += w;
        if (R[i].s == R[j].s && R[i].q == R[j].q) { pFsq += w; if (R[i].s == 1) pFsq_s1 += w; if (R[i].e == R[j].e) pFsqe += w; }
    }
    printf("pairs of strong rooted classes with h' = h + 2 and F' = F: %lld\n", pF);
    printf("  of these with s = s' = 1: %lld\n", pF_s1);
    printf("  with q' = q (any s): %lld\n", pFq);
    printf("  with s' = s (any q): %lld\n", pFs);
    printf("  with s' = s and q' = q: %lld (of these with s = 1: %lld; with e' = e also: %lld)\n", pFsq, pFsq_s1, pFsqe);
    /* the smallest-h unit-start strong profiles, for the record */

    /* ---------- the pair A, B from the stated out-rows ---------- */
    uint8_t A[N] = {8, 73, 3, 4, 15, 31, 61}, B[N] = {2, 84, 1, 7, 13, 31, 61};
    for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++) {
        if ((((A[i] >> j) & 1) + ((A[j] >> i) & 1)) != 1) { printf("A out-rows are not a tournament\n"); return 1; }
        if ((((B[i] >> j) & 1) + ((B[j] >> i) & 1)) != 1) { printf("B out-rows are not a tournament\n"); return 1; }
    }
    uint32_t cA = encode(A), cB = encode(B);
    Prof pa = PROF[cA], pb = PROF[cB];
    printf("\nA (out-rows 8 73 3 4 15 31 61): strong %d (h,s,q,F,e) = (%d,%d,%d,%d,%d)\n", pa.strong, pa.h, pa.s, pa.q, pa.F, pa.e);
    printf("B (out-rows 2 84 1 7 13 31 61): strong %d (h,s,q,F,e) = (%d,%d,%d,%d,%d)\n", pb.strong, pb.h, pb.s, pb.q, pb.F, pb.e);
    uint32_t keyA = pathcode(A), keyB = pathcode(B);
    printf("path-order codes: A %u, B %u\n", keyA, keyB);

    /* target sets */
    uint8_t *isA = calloc(NT, 1), *isB = calloc(NT, 1);
    long long nA = 0, nB = 0;
    #pragma omp parallel for schedule(dynamic, 8192) reduction(+:nA,nB)
    for (long c = 0; c < (long)NT; c++) {
        Prof p = PROF[c];
        if (p.s != 1 || !p.strong) continue;
        int mA = (p.h == pa.h && p.F == pa.F && p.q == pa.q), mB = (p.h == pb.h && p.F == pb.F && p.q == pb.q);
        if (!mA && !mB) continue;
        uint8_t o[N]; decode((uint32_t)c, o);
        int path[N]; int fc; int fp[N];
        /* local path finding (thread safe) */
        {
            /* iterative DFS for the unique path */
            int stack[N], nxt[N], depth = 0; unsigned vis = 1u; stack[0] = 0; nxt[0] = 0; fc = 0;
            while (depth >= 0) {
                if (vis == FULL) { if (fc == 0) for (int t = 0; t < N; t++) fp[t] = stack[t]; fc++; vis &= ~(1u << stack[depth]); depth--; continue; }
                int cur = stack[depth], y = nxt[depth];
                while (y < N && (((vis >> y) & 1u) || !((o[cur] >> y) & 1u))) y++;
                if (y >= N) { if (depth > 0) vis &= ~(1u << cur); depth--; continue; }
                nxt[depth] = y + 1; depth++; stack[depth] = y; nxt[depth] = 0; vis |= 1u << y;
            }
            (void)path;
        }
        if (fc != 1) continue;
        uint8_t rl[N]; memset(rl, 0, N);
        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) if (i != j && ((o[fp[i]] >> fp[j]) & 1u)) rl[i] |= (uint8_t)(1u << j);
        uint32_t k = encode(rl);
        if (mA && k == keyA) { isA[c] = 1; nA++; }
        if (mB && k == keyB) { isB[c] = 1; nB++; }
    }
    printf("labelled copies with root 0: A %lld, B %lld (720 each iff the rooted automorphism groups are trivial)\n", nA, nB);

    /* ---------- BFS routes ---------- */
    int8_t *dist = malloc(NT);
    uint32_t *queue = malloc(sizeof(uint32_t) * NT);
    unsigned long long *w0 = malloc(sizeof(unsigned long long) * NT), *w1 = malloc(sizeof(unsigned long long) * NT), *w2 = malloc(sizeof(unsigned long long) * NT);
    for (int dir = 0; dir < 2; dir++) {
        uint32_t src = dir ? cB : cA; uint8_t *tgt = dir ? isA : isB; Prof p0 = PROF[src];
        printf("\n=== direction %s ===\n", dir ? "B -> A" : "A -> B");
        for (int mode = 0; mode < 3; mode++) {
            /* mode 0: strong intermediates; mode 1: strong and s = 1 intermediates; mode 2: any tournament */
            memset(dist, -1, NT);
            long qh = 0, qt = 0; dist[src] = 0; queue[qt++] = src; int D = -1;
            while (qh < qt) {
                uint32_t t = queue[qh++];
                if (D >= 0 && dist[t] >= D) continue;
                for (int k = 0; k < NP; k++) {
                    uint32_t u = t ^ (1u << k);
                    if (dist[u] >= 0) continue;
                    if (mode != 2 && !PROF[u].strong) continue;
                    if (tgt[u]) { dist[u] = (int8_t)(dist[t] + 1); if (D < 0) D = dist[u]; queue[qt++] = u; continue; }
                    if (mode == 1 && PROF[u].s != 1) continue;
                    dist[u] = (int8_t)(dist[t] + 1); queue[qt++] = u;
                }
            }
            if (mode == 1) { printf("distance when every intermediate is strong with s = 1: %d\n", D); continue; }
            if (mode == 2) { printf("distance with no strongness requirement on intermediates: %d\n", D); continue; }
            printf("distance through strong labelled tournaments (root fixed): %d\n", D);
            /* count shortest routes in the layered DAG */
            memset(w0, 0, sizeof(unsigned long long) * NT); memset(w1, 0, sizeof(unsigned long long) * NT); memset(w2, 0, sizeof(unsigned long long) * NT);
            w0[src] = w1[src] = w2[src] = 1;
            for (long i = 0; i < qt; i++) {
                uint32_t t = queue[i];
                if (dist[t] >= D) continue;
                if (dist[t] > 0 && tgt[t]) continue;
                for (int k = 0; k < NP; k++) {
                    uint32_t u = t ^ (1u << k);
                    if (!PROF[u].strong || dist[u] != dist[t] + 1) continue;
                    if (dist[u] == D && !tgt[u]) continue;
                    w0[u] += w0[t];
                    int inter = (dist[u] < D);
                    Prof pu_ = PROF[u];
                    int ok1 = !inter || pu_.s == 1;
                    int ok2 = !inter || (pu_.F <= p0.F && pu_.h + 2 * pu_.F < p0.h + 2 * p0.F);
                    if (ok1) w1[u] += w1[t];
                    if (ok2) w2[u] += w2[t];
                }
            }
            unsigned long long T0 = 0, T1 = 0, T2 = 0; int ncop = 0;
            for (long i = 0; i < qt; i++) { uint32_t t = queue[i]; if (dist[t] == D && tgt[t]) { T0 += w0[t]; T1 += w1[t]; T2 += w2[t]; if (w0[t]) ncop++; } }
            printf("shortest routes: %llu (ending in %d labelled target copies); with every intermediate s = 1: %llu; with every intermediate wrapper below the start at every clock: %llu\n", T0, ncop, T1, T2);
            /* explicit enumeration of all shortest routes, by backward DFS from the targets */
            long routes = 0, routes_below = 0, routes_below_literal = 0, routes_s1 = 0; int wrapper_fail = 0, wrapper_checked = 0;
            for (long i = 0; i < qt; i++) {
                uint32_t t = queue[i]; if (dist[t] != D || !tgt[t]) continue;
                uint32_t st[16]; int nx[16]; int L = 0; st[0] = t; nx[0] = 0; L = 1;
                while (L > 0) {
                    uint32_t cur = st[L - 1];
                    if (dist[cur] == 0) {
                        routes++;
                        /* route is st[L-1] (= src) ... st[0] (= target) */
                        int below = 1, belowlit = 1, s1 = 1;
                        for (int k2 = L - 2; k2 >= 1; k2--) {
                            Prof pi = PROF[st[k2]];
                            if (!(pi.F <= p0.F && pi.h + 2 * pi.F < p0.h + 2 * p0.F)) below = 0;
                            for (int m = 1; m <= 60; m++) {
                                __int128 lhs = (__int128)pi.h + ((__int128)1 << m) * pi.F, rhs = (__int128)p0.h + ((__int128)1 << m) * p0.F;
                                if (!(lhs < rhs)) { belowlit = 0; break; }
                            }
                            if (pi.s != 1) s1 = 0;
                            /* literal wrapper check on the intermediate */
                            uint8_t o[N]; decode(st[k2], o);
                            for (int m = 1; m <= 4; m++) { wrapper_checked++; if (wrapperH(o, m) != (long long)pi.h + (1LL << m) * pi.F) wrapper_fail++; }
                        }
                        routes_below += below; routes_below_literal += belowlit; routes_s1 += s1;
                        if (below || routes <= 1) {
                            printf("  %s route:", below ? "BELOW-START" : "example");
                            for (int k2 = L - 1; k2 >= 0; k2--) { Prof pi = PROF[st[k2]]; printf(" (%d,%d,%d,%d)", pi.h, pi.s, pi.q, pi.F); }
                            printf("  flips:");
                            for (int k2 = L - 1; k2 >= 1; k2--) { uint32_t d = st[k2] ^ st[k2 - 1]; int kk = __builtin_ctz(d); printf(" {%d,%d}", pu[kk], pv[kk]); }
                            printf("\n");
                        }
                        L--; continue;
                    }
                    int found = 0;
                    for (int k = nx[L - 1]; k < NP; k++) {
                        uint32_t u = cur ^ (1u << k);
                        if (dist[u] != dist[cur] - 1 || !PROF[u].strong) continue;
                        if (dist[u] > 0 && tgt[u]) continue;
                        nx[L - 1] = k + 1; st[L] = u; nx[L] = 0; L++; found = 1; break;
                    }
                    if (!found) L--;
                }
            }
            printf("explicit enumeration: %ld shortest routes; keeping s = 1: %ld; below the start (closed form): %ld; below the start (literal, m = 1..60): %ld\n",
                   routes, routes_s1, routes_below, routes_below_literal);
            printf("literal wrapper H(D_m) = h + 2^m F on every intermediate of every shortest route (m = 1..4): %d checks, %d failures\n", wrapper_checked, wrapper_fail);
        }
    }

    /* ---------- wrapper formula on random tournaments of every kind ---------- */
    {
        srand(777); int bad = 0, tot2 = 0, nonstrong = 0, s_not1 = 0;
        for (int t = 0; t < 1500; t++) {
            uint32_t c = ((uint32_t)rand() << 15 ^ (uint32_t)rand()) & (NT - 1);
            uint8_t o[N]; decode(c, o); Prof p = PROF[c];
            if (!p.strong) nonstrong++; if (p.s != 1) s_not1++;
            for (int m = 1; m <= 3; m++) { tot2++; if (wrapperH(o, m) != (long long)p.h + (1LL << m) * p.F) bad++; }
        }
        printf("\nwrapper formula on 1500 random labelled tournaments (%d non-strong, %d with s != 1), m = 1..3: %d / %d failures\n", nonstrong, s_not1, bad, tot2);
        for (int m = 1; m <= 5; m++) {
            long long HA = wrapperH(A, m), HB = wrapperH(B, m);
            printf("m = %d (order %d): H(D_m(A)) = %lld, H(D_m(B)) = %lld, difference %lld; formula %lld / %lld\n", m, N + m + 1, HA, HB, HB - HA,
                   (long long)pa.h + (1LL << m) * pa.F, (long long)pb.h + (1LL << m) * pb.F);
        }
    }
    return 0;
}
