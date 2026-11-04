/*
 * plus_minus_two_port_routes_20261104.c  (opus-2026-11-04-S76)
 *
 * Checks the claim that two strong seven-vertex rooted components have the same port response but Hamilton path counts
 * differing by two, joined by a four-arc reversal route through strong tournaments.
 *
 * Rooted profile of (T, v), as in adaptive_metagraph_reformulation_20261103/ADAPTIVE.md section 2:
 *   h = H(T), s = Hamiltonian paths starting at v ("root start"), q = H(T - v),
 *   F = ordered spanning two-slot covers (P, Q) with P non-empty ending at v, Q possibly empty  (= E_v(T, 2) of Session189).
 * Session189's dyadic wrapper D_m(U, r) has H(D_m) = h + 2^m F (returns/PROOF.md (1)), so equal F and h' = h +- 2 shift every
 * wrapper of the clocked family by exactly +-2.
 *
 * Exhaustive over all 2^21 labelled tournaments on {0..6} with root 0:
 *   1. rooted classes (canonical form with the root individualised) and their profiles;
 *   2. pairs of strong rooted classes with |h - h'| = 2 and equal port data, under several choices of "port data";
 *   3. for each pair with equal (s, q, F) [the full root statistics other than h]: breadth-first search over strong labelled
 *      tournaments with the root fixed, one arc reversal per step, giving the exact distance, and the distance when every
 *      intermediate must also keep s = 1 or keep the wrapper below the start (h_i + 2^m F_i < h_0 + 2^m F_0 for all m >= 1).
 *
 * Usage: pm2 [MAXPAIRS]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>

#define N 7
#define NARC 21
#define NT (1u << NARC)
typedef unsigned long long u64;

static int arc_u[NARC], arc_v[NARC];

static void build(uint32_t code, uint32_t *out) {
    for (int i = 0; i < N; i++) out[i] = 0;
    for (int k = 0; k < NARC; k++) {
        if ((code >> k) & 1) out[arc_u[k]] |= 1u << arc_v[k]; else out[arc_v[k]] |= 1u << arc_u[k];
    }
}

static int strong(const uint32_t *out) {
    uint32_t full = (1u << N) - 1, seen = 1, fr = 1;
    while (fr) { int x = __builtin_ctz(fr); fr &= fr - 1; uint32_t nb = out[x] & ~seen; seen |= nb; fr |= nb; }
    if (seen != full) return 0;
    seen = 1; fr = 1;
    while (fr) {
        int x = __builtin_ctz(fr); fr &= fr - 1;
        uint32_t nb = 0; for (int y = 0; y < N; y++) if (!((seen >> y) & 1) && ((out[y] >> x) & 1)) nb |= 1u << y;
        seen |= nb; fr |= nb;
    }
    return seen == full;
}

/* profile with root 0 */
typedef struct { int h, s, q, F, e; } Prof;
static Prof profile(const uint32_t *out) {
    static __thread int dpe[1 << N][N];   /* Hamilton paths of T[S] ending at x */
    int Hs[1 << N];
    for (int S = 0; S < (1 << N); S++) for (int x = 0; x < N; x++) dpe[S][x] = 0;
    for (int x = 0; x < N; x++) dpe[1 << x][x] = 1;
    for (int S = 1; S < (1 << N); S++) for (int x = 0; x < N; x++) {
        int c = dpe[S][x]; if (!c) continue;
        uint32_t nb = out[x] & ~(uint32_t)S;
        while (nb) { int y = __builtin_ctz(nb); nb &= nb - 1; dpe[S | (1 << y)][y] += c; }
    }
    Hs[0] = 1;
    for (int S = 1; S < (1 << N); S++) { int t = 0; for (int x = 0; x < N; x++) t += dpe[S][x]; Hs[S] = t; }
    Prof p; int full = (1 << N) - 1;
    p.h = Hs[full]; p.q = Hs[full & ~1]; p.e = dpe[full][0];
    /* s: paths starting at 0 = paths of the converse ending at 0; compute directly */
    {
        int dps[1 << N][N]; memset(dps, 0, sizeof dps); dps[1][0] = 1;
        for (int S = 1; S < (1 << N); S++) if (S & 1) for (int x = 0; x < N; x++) {
            int c = dps[S][x]; if (!c) continue;
            uint32_t nb = out[x] & ~(uint32_t)S;
            while (nb) { int y = __builtin_ctz(nb); nb &= nb - 1; dps[S | (1 << y)][y] += c; }
        }
        int t = 0; for (int x = 0; x < N; x++) t += dps[full][x]; p.s = t;
    }
    int F = 0;
    for (int S = 1; S <= full; S++) if (S & 1) F += dpe[S][0] * Hs[full & ~S];
    p.F = F;
    return p;
}

/* rooted canonical code: minimum over permutations fixing 0 of the 21-bit code (brute force over 6! = 720) */
static int perms[720][N];
static uint32_t rooted_canon(const uint32_t *out) {
    uint32_t best = 0xFFFFFFFFu;
    for (int k = 0; k < 720; k++) {
        const int *p = perms[k];   /* p[new] = old */
        uint32_t c = 0;
        for (int a = 0; a < NARC; a++) {
            int i = arc_u[a], j = arc_v[a];        /* new labels i < j */
            if ((out[p[i]] >> p[j]) & 1) c |= 1u << a;
        }
        if (c < best) best = c;
    }
    return best;
}

int main(int argc, char **argv) {
    int maxpairs = argc > 1 ? atoi(argv[1]) : 50;
    { int k = 0; for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++) { arc_u[k] = i; arc_v[k] = j; k++; } }
    { /* permutations of {1..6}, 0 fixed */
        int a[6] = {1, 2, 3, 4, 5, 6}, c = 0;
        /* Heap's algorithm */
        int idx[6] = {0}; int i = 0;
        perms[c][0] = 0; for (int t = 0; t < 6; t++) perms[c][t + 1] = a[t]; c++;
        while (i < 6) {
            if (idx[i] < i) { int sw = (i % 2 == 0) ? 0 : idx[i]; int tmp = a[sw]; a[sw] = a[i]; a[i] = tmp;
                perms[c][0] = 0; for (int t = 0; t < 6; t++) perms[c][t + 1] = a[t]; c++; idx[i]++; i = 0; }
            else { idx[i] = 0; i++; }
        }
        if (c != 720) { fprintf(stderr, "perm count %d\n", c); return 1; }
    }
    uint32_t *canon_of = malloc(sizeof(uint32_t) * NT);
    unsigned char *isstrong = malloc(NT);
    Prof *prof = malloc(sizeof(Prof) * NT);
    #pragma omp parallel for schedule(dynamic, 4096)
    for (long t = 0; t < (long)NT; t++) {
        uint32_t out[N]; build((uint32_t)t, out);
        isstrong[t] = (unsigned char)strong(out);
        prof[t] = profile(out);
        canon_of[t] = rooted_canon(out);
    }
    /* rooted classes = labelled tournaments equal to their own canonical code */
    long nclass = 0, nstrong = 0;
    for (long t = 0; t < (long)NT; t++) if (canon_of[t] == (uint32_t)t) { nclass++; if (isstrong[t]) nstrong++; }
    printf("# plus_minus_two_port_routes_20261104: order 7, root 0\n");
    printf("rooted classes: %ld (strong %ld)\n", nclass, nstrong);
    uint32_t *cls = malloc(sizeof(uint32_t) * nclass); long nc = 0;
    for (long t = 0; t < (long)NT; t++) if (canon_of[t] == (uint32_t)t) cls[nc++] = (uint32_t)t;

    /* pair statistics under several port-data choices */
    long pairs_F = 0, pairs_F_s1 = 0, pairs_sqF = 0, pairs_sqF_s1 = 0, pairs_all_e = 0;
    uint32_t *pa = malloc(sizeof(uint32_t) * 200000), *pb = malloc(sizeof(uint32_t) * 200000); long npairs = 0;
    for (long i = 0; i < nc; i++) {
        uint32_t a = cls[i]; if (!isstrong[a]) continue;
        for (long j = 0; j < nc; j++) {
            uint32_t b = cls[j]; if (!isstrong[b]) continue;
            if (prof[b].h != prof[a].h + 2) continue;
            if (prof[a].F != prof[b].F) continue;
            pairs_F++;
            if (prof[a].s == 1 && prof[b].s == 1) pairs_F_s1++;
            if (prof[a].s == prof[b].s && prof[a].q == prof[b].q) {
                pairs_sqF++;
                if (prof[a].s == 1) pairs_sqF_s1++;
                if (prof[a].e == prof[b].e) pairs_all_e++;
                if (npairs < 200000) { pa[npairs] = a; pb[npairs] = b; npairs++; }
            }
        }
    }
    printf("strong rooted pairs with h' = h + 2 and equal F: %ld; of these with s = s' = 1: %ld\n", pairs_F, pairs_F_s1);
    printf("  with equal (s, q, F): %ld; of these with s = 1: %ld; with equal (s, q, F, e): %ld\n", pairs_sqF, pairs_sqF_s1, pairs_all_e);

    /* direct check of the wrapper formula H(D_m) = h + 2^m F on every strong rooted class with s = 1, m = 1..4 */
    {
        long bad = 0, checked = 0;
        for (long i = 0; i < nc; i++) {
            uint32_t a = cls[i]; if (!isstrong[a] || prof[a].s != 1) continue;
            uint32_t oa[N]; build(a, oa);
            for (int m = 1; m <= 4; m++) {
                int n2 = N + m + 1;                 /* vertices: U = 0..6, Z = 7..7+m-1, w = 7+m */
                uint32_t o2[16]; memset(o2, 0, sizeof o2);
                for (int v = 0; v < N; v++) o2[v] = oa[v];
                int w = N + m;
                for (int z = N; z < N + m; z++) { for (int v = 0; v < N; v++) o2[z] |= 1u << v; for (int z2 = z + 1; z2 < N + m; z2++) o2[z] |= 1u << z2; }
                for (int z = N; z < N + m; z++) o2[w] |= 1u << z;          /* w -> Z */
                for (int v = 1; v < N; v++) o2[w] |= 1u << v;              /* w -> U \ {r} */
                o2[0] |= 1u << w;                                          /* r -> w */
                /* Hamilton path count by subset DP */
                static long long dp[1 << 12][12];
                int full = (1 << n2) - 1;
                for (int S = 0; S <= full; S++) for (int x = 0; x < n2; x++) dp[S][x] = 0;
                for (int x = 0; x < n2; x++) dp[1 << x][x] = 1;
                for (int S = 1; S <= full; S++) for (int x = 0; x < n2; x++) { long long c = dp[S][x]; if (!c) continue;
                    uint32_t nb = o2[x] & ~(uint32_t)S; while (nb) { int y = __builtin_ctz(nb); nb &= nb - 1; dp[S | (1 << y)][y] += c; } }
                long long H = 0; for (int x = 0; x < n2; x++) H += dp[full][x];
                checked++;
                if (H != (long long)prof[a].h + (1LL << m) * prof[a].F) bad++;
            }
        }
        printf("wrapper formula H(D_m) = h + 2^m F checked directly on %ld (class, m) cases (strong rooted classes with s = 1, m = 1..4): %ld failures\n", checked, bad);
    }
    /* the unit-root-start pairs with equal F: full route analysis in both directions */
    for (long i = 0; i < nc; i++) {
        uint32_t a = cls[i]; if (!isstrong[a] || prof[a].s != 1) continue;
        for (long j = 0; j < nc; j++) {
            uint32_t b = cls[j]; if (!isstrong[b] || prof[b].s != 1) continue;
            if (prof[b].h != prof[a].h + 2 || prof[a].F != prof[b].F) continue;
            uint32_t oa[N], ob[N]; build(a, oa); build(b, ob);
            printf("\nUNIT PAIR: A (h,s,q,F,e) = (%d,%d,%d,%d,%d)  B = (%d,%d,%d,%d,%d)\n",
                   prof[a].h, prof[a].s, prof[a].q, prof[a].F, prof[a].e, prof[b].h, prof[b].s, prof[b].q, prof[b].F, prof[b].e);
            printf("  A out-rows (root 0):"); for (int v = 0; v < N; v++) printf(" %u", oa[v]); printf("   scores:");
            for (int v = 0; v < N; v++) printf(" %d", __builtin_popcount(oa[v])); printf("\n");
            printf("  B out-rows (root 0):"); for (int v = 0; v < N; v++) printf(" %u", ob[v]); printf("   scores:");
            for (int v = 0; v < N; v++) printf(" %d", __builtin_popcount(ob[v])); printf("\n");
            for (int dir = 0; dir < 2; dir++) {
                uint32_t src = dir ? b : a, dst = dir ? a : b;
                /* BFS layers from src over strong labelled tournaments (root fixed); targets = labelled copies of dst */
                int *dist = malloc(sizeof(int) * NT);
                for (long t = 0; t < (long)NT; t++) dist[t] = -1;
                uint32_t *q = malloc(sizeof(uint32_t) * NT); long qh = 0, qt = 0;
                dist[src] = 0; q[qt++] = src; int D = -1;
                while (qh < qt) {
                    uint32_t t = q[qh++];
                    if (D >= 0 && dist[t] >= D) break;
                    for (int arc = 0; arc < NARC; arc++) {
                        uint32_t u = t ^ (1u << arc);
                        if (dist[u] >= 0 || !isstrong[u]) continue;
                        dist[u] = dist[t] + 1;
                        if (canon_of[u] == dst && D < 0) D = dist[u];
                        q[qt++] = u;
                    }
                }
                /* count shortest routes into each target copy, and those that keep s = 1 / keep wrappers below the start */
                /* dynamic programme over BFS layers: ways[t] (all), ways1[t] (intermediates with s = 1), waysw[t] (intermediates below start) */
                u64 *w0 = calloc(NT, sizeof(u64)), *w1 = calloc(NT, sizeof(u64)), *w2 = calloc(NT, sizeof(u64));
                w0[src] = w1[src] = w2[src] = 1;
                for (long idx = 0; idx < qt; idx++) {
                    uint32_t t = q[idx]; if (dist[t] >= D) continue;
                    for (int arc = 0; arc < NARC; arc++) {
                        uint32_t u = t ^ (1u << arc);
                        if (!isstrong[u] || dist[u] != dist[t] + 1) continue;
                        int last = (dist[u] == D);
                        if (last && canon_of[u] != dst) continue;
                        if (!last && canon_of[u] == dst) continue;
                        w0[u] += w0[t];
                        int ok1 = last || prof[u].s == 1;
                        int ok2 = last || (prof[u].F <= prof[src].F && prof[u].h + 2 * prof[u].F < prof[src].h + 2 * prof[src].F);
                        if (ok1) w1[u] += w1[t];
                        if (ok2) w2[u] += w2[t];
                    }
                }
                u64 T0 = 0, T1 = 0, T2 = 0; long ncopies = 0;
                for (long idx = 0; idx < qt; idx++) { uint32_t t = q[idx]; if (dist[t] == D && canon_of[t] == dst) { T0 += w0[t]; T1 += w1[t]; T2 += w2[t]; ncopies++; } }
                printf("  route %s: shortest distance %d through strong tournaments; shortest routes %llu (into %ld labelled copies);"
                       " keeping s = 1 throughout: %llu; keeping every intermediate wrapper below the start at every clock m >= 1 (F_i <= F_0 and h_i + 2F_i < h_0 + 2F_0): %llu\n",
                       dir ? "B->A" : "A->B", D, T0, ncopies, T1, T2);
                /* print one shortest route with intermediate profiles */
                {
                    uint32_t tgt = 0; for (long idx = 0; idx < qt; idx++) { uint32_t t = q[idx]; if (dist[t] == D && canon_of[t] == dst && w0[t]) { tgt = t; break; } }
                    uint32_t path[16]; int L = 0; path[L++] = tgt;
                    while (dist[path[L - 1]] > 0) {
                        uint32_t t = path[L - 1];
                        for (int arc = 0; arc < NARC; arc++) { uint32_t u = t ^ (1u << arc); if (isstrong[u] && dist[u] == dist[t] - 1 && w0[u]) { path[L++] = u; break; } }
                    }
                    printf("    one shortest route (h,s,q,F):");
                    for (int k2 = L - 1; k2 >= 0; k2--) { Prof pp = prof[path[k2]]; printf(" (%d,%d,%d,%d)", pp.h, pp.s, pp.q, pp.F); }
                    printf("\n");
                }
                /* print every shortest route that keeps the wrappers below the start */
                {
                    for (long idx = 0; idx < qt; idx++) {
                        uint32_t t = q[idx]; if (dist[t] != D || canon_of[t] != dst || !w2[t]) continue;
                        /* depth-first enumeration of predecessor chains with w2 > 0 */
                        uint32_t stack[16]; int arcs[16]; int L = 0; stack[0] = t; arcs[0] = -1; L = 1;
                        while (L > 0) {
                            uint32_t cur = stack[L - 1];
                            if (dist[cur] == 0) {
                                printf("    below-start route (h,s,q,F):");
                                for (int k2 = L - 1; k2 >= 0; k2--) { Prof pp = prof[stack[k2]]; printf(" (%d,%d,%d,%d)", pp.h, pp.s, pp.q, pp.F); }
                                printf("   reversed arcs:");
                                for (int k2 = L - 1; k2 >= 1; k2--) { uint32_t diff = stack[k2] ^ stack[k2 - 1]; int a2 = __builtin_ctz(diff); printf(" {%d,%d}", arc_u[a2], arc_v[a2]); }
                                printf("\n");
                                L--; continue;
                            }
                            int found2 = 0;
                            for (int arc = arcs[L - 1] + 1; arc < NARC; arc++) {
                                uint32_t u = cur ^ (1u << arc);
                                if (!isstrong[u] || dist[u] != dist[cur] - 1 || !w2[u]) continue;
                                int okk = (dist[u] == 0) || (prof[u].F <= prof[src].F && prof[u].h + 2 * prof[u].F < prof[src].h + 2 * prof[src].F);
                                if (!okk) continue;
                                arcs[L - 1] = arc; stack[L] = u; arcs[L] = -1; L++; found2 = 1; break;
                            }
                            if (!found2) L--;
                        }
                    }
                }
                /* minimum distance if every intermediate keeps s = 1 */
                {
                    for (long t = 0; t < (long)NT; t++) dist[t] = -1;
                    qh = qt = 0; dist[src] = 0; q[qt++] = src; int D1 = -1;
                    while (qh < qt && D1 < 0) {
                        uint32_t t = q[qh++];
                        for (int arc = 0; arc < NARC; arc++) {
                            uint32_t u = t ^ (1u << arc);
                            if (dist[u] >= 0 || !isstrong[u]) continue;
                            if (canon_of[u] == dst) { D1 = dist[t] + 1; break; }
                            if (prof[u].s != 1) continue;
                            dist[u] = dist[t] + 1; q[qt++] = u;
                        }
                    }
                    printf("    distance when every intermediate keeps s = 1: %d\n", D1);
                }
                free(w0); free(w1); free(w2); free(q); free(dist);
            }
        }
    }
    return 0;
}
