#!/usr/bin/env python3
"""
n=8 tiling isomorphism analysis using C acceleration.
Strategy: compute (sorted_scores, H, #3-cycles, vertex_3cycle_profile) for all 2^21 masks.
Group by fingerprint for approximate iso classes. Do targeted iso checks for SC+SF kernel.

Instance: opus-2026-03-06-S8
"""

import ctypes
import tempfile
import os
import subprocess
import time

C_CODE = r"""
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define N 8
#define M 21
#define NMASK (1 << M)  /* 2097152 */

/* Tile setup */
static int tile_x[M], tile_y[M];
static int trans_map[M];
static int verts[N]; /* [8,7,6,5,4,3,2,1] */

void setup_tiles() {
    int idx = 0;
    for (int y = 1; y <= N-2; y++)
        for (int x = N; x > y+1; x--)
            { tile_x[idx] = x; tile_y[idx] = y; idx++; }

    for (int i = 0; i < N; i++) verts[i] = N - i;

    /* Build transpose map: tile (x,y) -> (N-y+1, N-x+1) */
    int tile_lookup[N+1][N+1];
    memset(tile_lookup, -1, sizeof(tile_lookup));
    for (int i = 0; i < M; i++) tile_lookup[tile_x[i]][tile_y[i]] = i;
    for (int i = 0; i < M; i++) {
        int nx = N - tile_y[i] + 1;
        int ny = N - tile_x[i] + 1;
        trans_map[i] = tile_lookup[nx][ny];
    }
}

static inline int vert_idx(int v) {
    return N - v;  /* verts[i] = N-i, so index of v is N-v */
}

void mask_to_adj(int mask, int A[N][N]) {
    memset(A, 0, sizeof(int)*N*N);
    for (int k = 0; k < N-1; k++) A[k][k+1] = 1;
    for (int i = 0; i < M; i++) {
        int xi = vert_idx(tile_x[i]);
        int yi = vert_idx(tile_y[i]);
        if ((mask >> i) & 1) A[yi][xi] = 1;
        else A[xi][yi] = 1;
    }
}

/* Hamiltonian path count via bitmask DP */
long long count_ham_dp(const int A[N][N]) {
    static long long dp[1<<N][N];
    memset(dp, 0, sizeof(dp));
    for (int v = 0; v < N; v++) dp[1<<v][v] = 1;
    for (int mask = 1; mask < (1<<N); mask++) {
        for (int v = 0; v < N; v++) {
            if (!dp[mask][v]) continue;
            if (!(mask & (1<<v))) continue;
            for (int u = 0; u < N; u++) {
                if (mask & (1<<u)) continue;
                if (A[v][u]) dp[mask|(1<<u)][u] += dp[mask][v];
            }
        }
    }
    long long total = 0;
    int full = (1<<N)-1;
    for (int v = 0; v < N; v++) total += dp[full][v];
    return total;
}

/* Count directed 3-cycles */
int count_3cycles(const int A[N][N]) {
    int count = 0;
    for (int i = 0; i < N; i++)
        for (int j = i+1; j < N; j++)
            for (int k = j+1; k < N; k++) {
                if (A[i][j] && A[j][k] && A[k][i]) count++;
                if (A[i][k] && A[k][j] && A[j][i]) count++;
            }
    return count;
}

/* Per-vertex 3-cycle count (sorted) -> 8 values */
void vertex_3cycle_profile(const int A[N][N], int profile[N]) {
    for (int v = 0; v < N; v++) {
        int cnt = 0;
        for (int j = 0; j < N; j++) {
            if (j == v) continue;
            for (int k = j+1; k < N; k++) {
                if (k == v) continue;
                /* Check if {v,j,k} forms a 3-cycle containing v */
                if (A[v][j] && A[j][k] && A[k][v]) cnt++;
                if (A[v][k] && A[k][j] && A[j][v]) cnt++;
            }
        }
        profile[v] = cnt;
    }
    /* Sort ascending */
    for (int i = 0; i < N-1; i++)
        for (int j = i+1; j < N; j++)
            if (profile[j] < profile[i]) { int t=profile[i]; profile[i]=profile[j]; profile[j]=t; }
}

/* Compute sorted out-degrees */
void compute_scores(const int A[N][N], int degs[N]) {
    for (int i = 0; i < N; i++) {
        degs[i] = 0;
        for (int j = 0; j < N; j++) degs[i] += A[i][j];
    }
    /* Sort descending */
    for (int i = 0; i < N-1; i++)
        for (int j = i+1; j < N; j++)
            if (degs[j] > degs[i]) { int t=degs[i]; degs[i]=degs[j]; degs[j]=t; }
}

/* Check grid symmetry */
int is_grid_sym(int mask) {
    for (int i = 0; i < M; i++)
        if (((mask>>i)&1) != ((mask>>trans_map[i])&1)) return 0;
    return 1;
}

/* ============================================================
   Isomorphism testing via partition refinement + backtracking
   ============================================================ */

/* Check if permutation p maps A to B: B[i][j] = A[p[i]][p[j]] */
int check_iso_perm(const int A[N][N], const int B[N][N], const int p[N]) {
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            if (A[p[i]][p[j]] != B[i][j]) return 0;
    return 1;
}

/* Group vertices by out-degree, then try all permutations consistent with partition */
typedef struct {
    int deg;     /* out-degree */
    int idx;     /* original vertex index */
} VertInfo;

int cmp_vert(const void *a, const void *b) {
    return ((const VertInfo*)b)->deg - ((const VertInfo*)a)->deg;
}

/* Try all mappings from B-vertices to A-vertices respecting degree partition */
static int found_iso;
static int best_perm[N];

void try_mappings(const int A[N][N], const int B[N][N],
                  int *a_groups, int *a_verts_sorted, int *b_verts_sorted,
                  int *group_sizes, int ngroups,
                  int perm[N], int group_idx, int pos_in_group,
                  int *a_used) {
    if (found_iso) return;

    if (group_idx == ngroups) {
        /* Full permutation built - check it */
        if (check_iso_perm(A, B, perm)) {
            found_iso = 1;
            memcpy(best_perm, perm, N*sizeof(int));
        }
        return;
    }

    int gs = group_sizes[group_idx];
    int a_start = 0;
    for (int g = 0; g < group_idx; g++) a_start += group_sizes[g];
    int b_start = a_start;

    if (pos_in_group == gs) {
        /* Move to next group */
        try_mappings(A, B, a_groups, a_verts_sorted, b_verts_sorted,
                     group_sizes, ngroups, perm, group_idx+1, 0, a_used);
        return;
    }

    int b_v = b_verts_sorted[b_start + pos_in_group];

    /* Try mapping b_v to each unused a_v in this group */
    for (int k = 0; k < gs; k++) {
        int a_v = a_verts_sorted[a_start + k];
        if (a_used[a_v]) continue;

        perm[b_v] = a_v;
        a_used[a_v] = 1;

        /* Early pruning: check consistency with already-mapped vertices */
        int ok = 1;
        for (int g2 = 0; g2 <= group_idx && ok; g2++) {
            int s2 = 0;
            for (int g3 = 0; g3 < g2; g3++) s2 += group_sizes[g3];
            int lim = (g2 < group_idx) ? group_sizes[g2] : pos_in_group;
            for (int p2 = 0; p2 < lim && ok; p2++) {
                int b_w = b_verts_sorted[s2 + p2];
                /* perm[b_w] is already set */
                if (A[perm[b_v]][perm[b_w]] != B[b_v][b_w]) ok = 0;
                if (ok && A[perm[b_w]][perm[b_v]] != B[b_w][b_v]) ok = 0;
            }
        }

        if (ok) {
            try_mappings(A, B, a_groups, a_verts_sorted, b_verts_sorted,
                         group_sizes, ngroups, perm, group_idx, pos_in_group+1, a_used);
        }

        a_used[a_v] = 0;
        if (found_iso) return;
    }
}

/* Are A and B isomorphic? */
int are_isomorphic(const int A[N][N], const int B[N][N]) {
    /* Quick check: score sequences */
    int da[N], db[N];
    compute_scores(A, da);
    compute_scores(B, db);
    if (memcmp(da, db, N*sizeof(int)) != 0) return 0;

    /* Group A-vertices and B-vertices by out-degree */
    VertInfo va[N], vb[N];
    for (int i = 0; i < N; i++) {
        int sa = 0, sb = 0;
        for (int j = 0; j < N; j++) { sa += A[i][j]; sb += B[i][j]; }
        va[i].deg = sa; va[i].idx = i;
        vb[i].deg = sb; vb[i].idx = i;
    }
    qsort(va, N, sizeof(VertInfo), cmp_vert);
    qsort(vb, N, sizeof(VertInfo), cmp_vert);

    /* Identify groups */
    int a_sorted[N], b_sorted[N];
    int group_sizes[N], ngroups = 0;
    for (int i = 0; i < N; i++) { a_sorted[i] = va[i].idx; b_sorted[i] = vb[i].idx; }

    int i = 0;
    while (i < N) {
        int j = i;
        while (j < N && va[j].deg == va[i].deg) j++;
        group_sizes[ngroups++] = j - i;
        i = j;
    }

    /* Try all permutations respecting partition */
    found_iso = 0;
    int perm[N], a_used[N];
    memset(perm, 0, sizeof(perm));
    memset(a_used, 0, sizeof(a_used));

    try_mappings(A, B, NULL, a_sorted, b_sorted,
                 group_sizes, ngroups, perm, 0, 0, a_used);

    return found_iso;
}

/* ============================================================
   Fingerprint and class tracking
   ============================================================ */

/* Second-order scores: for each vertex, sum of out-degrees of out-neighbors */
void second_order_scores(const int A[N][N], int s2[N]) {
    int degs[N];
    for (int i = 0; i < N; i++) {
        degs[i] = 0;
        for (int j = 0; j < N; j++) degs[i] += A[i][j];
    }
    for (int i = 0; i < N; i++) {
        s2[i] = 0;
        for (int j = 0; j < N; j++)
            if (A[i][j]) s2[i] += degs[j];
    }
    /* Sort ascending */
    for (int i = 0; i < N-1; i++)
        for (int j = i+1; j < N; j++)
            if (s2[j] < s2[i]) { int t=s2[i]; s2[i]=s2[j]; s2[j]=t; }
}

/* Common out-neighbor count hash: for each pair (i,j), count common out-neighbors */
unsigned int common_out_hash(const int A[N][N]) {
    int counts[N*(N-1)/2]; /* 28 values for n=8 */
    int idx = 0;
    for (int i = 0; i < N; i++)
        for (int j = i+1; j < N; j++) {
            int c = 0;
            for (int k = 0; k < N; k++)
                if (A[i][k] && A[j][k]) c++;
            counts[idx++] = c;
        }
    /* Sort */
    for (int i = 0; i < idx-1; i++)
        for (int j = i+1; j < idx; j++)
            if (counts[j] < counts[i]) { int t=counts[i]; counts[i]=counts[j]; counts[j]=t; }
    /* Hash */
    unsigned int h = 0;
    for (int i = 0; i < idx; i++) h = h * 37 + counts[i];
    return h;
}

/* Fingerprint: (sorted_scores[8], H, #3cycles, sorted_v3c_profile[8], s2_scores[8], co_hash) */
typedef struct {
    int scores[N];
    long long H;
    int c3;
    int v3c[N];
    int s2[N];
    unsigned int co_hash;
} Fingerprint;

int cmp_fp(const Fingerprint *a, const Fingerprint *b) {
    for (int i = 0; i < N; i++) {
        if (a->scores[i] != b->scores[i]) return a->scores[i] - b->scores[i];
    }
    if (a->H != b->H) return (a->H < b->H) ? -1 : 1;
    if (a->c3 != b->c3) return a->c3 - b->c3;
    for (int i = 0; i < N; i++) {
        if (a->v3c[i] != b->v3c[i]) return a->v3c[i] - b->v3c[i];
    }
    for (int i = 0; i < N; i++) {
        if (a->s2[i] != b->s2[i]) return a->s2[i] - b->s2[i];
    }
    if (a->co_hash != b->co_hash) return (a->co_hash < b->co_hash) ? -1 : 1;
    return 0;
}

/* Store fingerprints compactly */
static Fingerprint *fps;  /* allocated for NMASK */
static int *fp_class;     /* fingerprint class index for each mask */

/* Fingerprint class info */
#define MAX_FP_CLASSES 8000
typedef struct {
    Fingerprint fp;
    int count;         /* number of masks */
    int is_palindromic; /* score sequence palindromic? */
    int gs_count;      /* grid-symmetric masks in this class */
    int first_mask;    /* first mask in this class */
    int is_sf;         /* self-flip verified */
    int is_sc;         /* self-converse verified */
} FPClass;

static FPClass fp_classes[MAX_FP_CLASSES];
static int n_fp_classes = 0;

int find_fp_class(const Fingerprint *fp) {
    /* Binary search since we'll sort */
    for (int i = 0; i < n_fp_classes; i++) {
        if (cmp_fp(fp, &fp_classes[i].fp) == 0) return i;
    }
    return -1;
}

/* Check if score sequence is palindromic: d[i] + d[N-1-i] = N-1 */
int is_palindromic_scores(const int scores[N]) {
    for (int i = 0; i < N/2; i++) {
        if (scores[i] + scores[N-1-i] != N-1) return 0;
    }
    return 1;
}

/* Count automorphisms (for a specific matrix) */
int count_aut(const int A[N][N]) {
    /* Use partition refinement approach */
    int count = 0;
    /* Group vertices by out-degree */
    VertInfo vi[N];
    for (int i = 0; i < N; i++) {
        int s = 0;
        for (int j = 0; j < N; j++) s += A[i][j];
        vi[i].deg = s; vi[i].idx = i;
    }
    qsort(vi, N, sizeof(VertInfo), cmp_vert);

    int sorted[N], group_sizes[N], ngroups = 0;
    for (int i = 0; i < N; i++) sorted[i] = vi[i].idx;
    int i = 0;
    while (i < N) {
        int j = i;
        while (j < N && vi[j].deg == vi[i].deg) j++;
        group_sizes[ngroups++] = j - i;
        i = j;
    }

    /* Try all permutations respecting partition and check if automorphism */
    /* For automorphism: A[p[i]][p[j]] = A[i][j] */
    /* Reuse try_mappings with A=A, B=A */
    found_iso = 0;

    /* Actually, count all, not just first. Need a different approach. */
    /* For simplicity, just enumerate partition-consistent permutations */

    /* Compute total permutations to try */
    long long total_perms = 1;
    for (int g = 0; g < ngroups; g++) {
        long long f = 1;
        for (int k = 2; k <= group_sizes[g]; k++) f *= k;
        total_perms *= f;
    }

    if (total_perms > 100000) {
        /* Too many, use simpler bound */
        return -1;  /* indicate unknown */
    }

    /* Generate all partition-consistent permutations using recursion */
    /* Simple approach: for each group, try all permutations of that group's vertices */

    /* Build arrays of vertices per group */
    int grp_start[N], grp_end[N];
    int pos = 0;
    for (int g = 0; g < ngroups; g++) {
        grp_start[g] = pos;
        grp_end[g] = pos + group_sizes[g];
        pos += group_sizes[g];
    }

    /* Use iterative approach with group permutation tracking */
    int perm[N];
    for (int i = 0; i < N; i++) perm[i] = sorted[i];

    /* Count automorphisms by trying all partition-consistent perms */
    /* This is complex to do iteratively, let me use a simpler approach:
       For small groups (<= 4! * 4! = 576 worst case), just enumerate */
    count = 0;

    /* Recursive enumeration of partition-consistent permutations */
    /* Use arrays: for group g, store indices into sorted[] */
    int indices[N]; /* indices[i] = which sorted[] entry maps to vertex sorted[i] */
    for (int i = 0; i < N; i++) indices[i] = i;

    /* Generate all permutations within each group using next_permutation */
    /* Start with sorted indices per group */
    /* Actually, let me just use a brute-force approach for small counts */

    /* For the kernel classes (which is what we care about), |Aut| > 1 means
       specific structure. Let me just count via brute force for small perm counts. */

    /* Simplified: try all perms within groups */
    int grp_perm[N][24]; /* max group size assumed <= 4, 4! = 24 */
    int grp_perm_count[N];

    for (int g = 0; g < ngroups; g++) {
        int sz = group_sizes[g];
        int elems[8];
        for (int i = 0; i < sz; i++) elems[i] = sorted[grp_start[g] + i];

        grp_perm_count[g] = 0;

        /* Generate all permutations of elems[0..sz-1] */
        /* Store as flat arrays: grp_perm[g][k*sz + i] */
        /* Use simple recursive generation */
        /* For sz <= 4, just enumerate all sz! permutations */
        if (sz == 1) {
            grp_perm[g][0] = elems[0];
            grp_perm_count[g] = 1;
        } else if (sz == 2) {
            grp_perm[g][0] = elems[0]; grp_perm[g][1] = elems[1];
            grp_perm[g][2] = elems[1]; grp_perm[g][3] = elems[0];
            grp_perm_count[g] = 2;
        } else if (sz == 3) {
            int p3[6][3] = {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
            for (int k = 0; k < 6; k++)
                for (int i = 0; i < 3; i++)
                    grp_perm[g][k*3+i] = elems[p3[k][i]];
            grp_perm_count[g] = 6;
        } else if (sz == 4) {
            int p4[24][4];
            int cnt4 = 0;
            int a[4] = {0,1,2,3};
            /* Generate all 24 permutations */
            for (int a0=0;a0<4;a0++)
                for (int a1=0;a1<4;a1++) { if(a1==a0)continue;
                    for (int a2=0;a2<4;a2++) { if(a2==a0||a2==a1)continue;
                        int a3=6-a0-a1-a2;
                        p4[cnt4][0]=a0; p4[cnt4][1]=a1;
                        p4[cnt4][2]=a2; p4[cnt4][3]=a3;
                        cnt4++;
                    }
                }
            for (int k = 0; k < 24; k++)
                for (int i = 0; i < 4; i++)
                    grp_perm[g][k*4+i] = elems[p4[k][i]];
            grp_perm_count[g] = 24;
        } else {
            /* sz >= 5: too many, skip detailed enumeration */
            return -1;
        }
    }

    /* Now enumerate all combinations of group permutations */
    int grp_idx[N]; /* current permutation index per group */
    memset(grp_idx, 0, sizeof(grp_idx));

    while (1) {
        /* Build perm from current group indices */
        for (int g = 0; g < ngroups; g++) {
            int sz = group_sizes[g];
            int ki = grp_idx[g];
            for (int i = 0; i < sz; i++) {
                /* Map: vertex sorted[grp_start[g]+i] -> grp_perm[g][ki*sz+i] */
                perm[sorted[grp_start[g]+i]] = grp_perm[g][ki*sz+i];
            }
        }

        /* Check if this is an automorphism */
        int ok = 1;
        for (int ii = 0; ii < N && ok; ii++)
            for (int jj = 0; jj < N && ok; jj++)
                if (A[perm[ii]][perm[jj]] != A[ii][jj]) ok = 0;
        if (ok) count++;

        /* Increment */
        int carry = 1;
        for (int g = ngroups-1; g >= 0 && carry; g--) {
            grp_idx[g]++;
            if (grp_idx[g] >= grp_perm_count[g]) { grp_idx[g] = 0; }
            else { carry = 0; }
        }
        if (carry) break; /* all combinations tried */
    }

    return count;
}

/* Count anti-automorphisms (maps T to T^op) */
int count_anti_aut(const int A[N][N]) {
    /* Similar to count_aut but check A[p[i]][p[j]] == A[j][i] */
    /* Score partition: p maps to converse, so out-degrees map to in-degrees */
    /* Vertex with out-degree d maps to vertex with in-degree d = out-degree (N-1-d) */

    /* For anti-aut: A[p[i]][p[j]] = A[j][i] for all i,j */
    /* This means p maps vertex i (outdeg d_i) to a vertex with outdeg = indeg of i = N-1-d_i */

    int degs[N];
    for (int i = 0; i < N; i++) {
        degs[i] = 0;
        for (int j = 0; j < N; j++) degs[i] += A[i][j];
    }

    /* Check if palindromic first */
    int sorted_degs[N];
    memcpy(sorted_degs, degs, sizeof(degs));
    for (int i = 0; i < N-1; i++)
        for (int j = i+1; j < N; j++)
            if (sorted_degs[j] > sorted_degs[i])
                { int t=sorted_degs[i]; sorted_degs[i]=sorted_degs[j]; sorted_degs[j]=t; }
    if (!is_palindromic_scores(sorted_degs)) return 0;

    /* Group by degree, with mapping: deg d -> deg N-1-d */
    /* Vertices with deg d must map to vertices with deg N-1-d */
    /* Build buckets */
    int bucket[N][N], bsize[N]; /* bucket[d][k] = k-th vertex with out-degree d */
    memset(bsize, 0, sizeof(bsize));
    for (int i = 0; i < N; i++) bucket[degs[i]][bsize[degs[i]]++] = i;

    /* For each source vertex with deg d, target must be in bucket[N-1-d] */
    /* Enumerate valid mappings */
    int count = 0;

    /* Simple brute force for small buckets */
    /* Compute total valid permutations */
    long long total = 1;
    int used_deg[N]; memset(used_deg, 0, sizeof(used_deg));
    for (int d = 0; d < N; d++) {
        if (bsize[d] == 0 || used_deg[d]) continue;
        int td = N-1-d;
        if (bsize[d] != bsize[td]) return 0;
        used_deg[d] = used_deg[td] = 1;
        long long f = 1;
        for (int k = 2; k <= bsize[d]; k++) f *= k;
        total *= f;
    }
    if (total > 100000) return -1;

    /* Enumerate all valid mappings using backtracking */
    /* Order vertices: sorted by degree */
    int order[N];
    for (int i = 0; i < N; i++) order[i] = i;
    /* Sort by degree for efficiency */
    for (int i = 0; i < N-1; i++)
        for (int j = i+1; j < N; j++)
            if (degs[order[j]] > degs[order[i]])
                { int t=order[i]; order[i]=order[j]; order[j]=t; }

    /* Backtracking */
    int perm[N];
    int target_used[N];
    memset(perm, -1, sizeof(perm));
    memset(target_used, 0, sizeof(target_used));

    /* Stack-based backtracking */
    int stack_pos = 0;
    int stack_choice[N]; /* which choice for each position */
    memset(stack_choice, 0, sizeof(stack_choice));

    while (stack_pos >= 0) {
        if (stack_pos == N) {
            /* Check full anti-automorphism */
            int ok = 1;
            for (int i = 0; i < N && ok; i++)
                for (int j = 0; j < N && ok; j++)
                    if (A[perm[i]][perm[j]] != A[j][i]) ok = 0;
            if (ok) count++;
            stack_pos--;
            if (stack_pos >= 0) {
                target_used[perm[order[stack_pos]]] = 0;
                perm[order[stack_pos]] = -1;
                /* stack_choice already has resume value from forward pass */
            }
            continue;
        }

        int src = order[stack_pos];
        int target_deg = N - 1 - degs[src];
        int found = 0;

        for (int k = stack_choice[stack_pos]; k < bsize[target_deg]; k++) {
            int tgt = bucket[target_deg][k];
            if (target_used[tgt]) continue;

            /* Early pruning: check against already-assigned pairs */
            int ok = 1;
            for (int pp = 0; pp < stack_pos && ok; pp++) {
                int si = order[pp];
                if (A[perm[src]][perm[si]] != A[si][src]) ok = 0;
                if (ok && A[perm[si]][perm[src]] != A[src][si]) ok = 0;
            }

            if (!ok) continue;

            perm[src] = tgt;
            target_used[tgt] = 1;
            stack_choice[stack_pos] = k + 1;
            stack_pos++;
            stack_choice[stack_pos] = 0;
            found = 1;
            break;
        }

        if (!found) {
            /* Backtrack */
            stack_choice[stack_pos] = 0;
            stack_pos--;
            if (stack_pos >= 0) {
                target_used[perm[order[stack_pos]]] = 0;
                perm[order[stack_pos]] = -1;
                /* stack_choice[stack_pos] already advanced */
            }
        }
    }

    return count;
}

/* ============================================================
   Main
   ============================================================ */

int main() {
    setup_tiles();

    fprintf(stderr, "n=%d, m=%d, 2^m=%d\n", N, M, NMASK);

    /* Allocate fingerprint storage */
    fps = (Fingerprint *)malloc(NMASK * sizeof(Fingerprint));
    fp_class = (int *)malloc(NMASK * sizeof(int));
    if (!fps || !fp_class) { fprintf(stderr, "malloc failed\n"); return 1; }

    /* Phase 1: Compute fingerprints for all masks */
    fprintf(stderr, "Phase 1: Computing fingerprints...\n");
    long long max_H = 0;
    int max_H_mask = 0;

    for (int mask = 0; mask < NMASK; mask++) {
        int A[N][N];
        mask_to_adj(mask, A);

        compute_scores(A, fps[mask].scores);
        fps[mask].H = count_ham_dp(A);
        fps[mask].c3 = count_3cycles(A);
        vertex_3cycle_profile(A, fps[mask].v3c);
        second_order_scores(A, fps[mask].s2);
        fps[mask].co_hash = common_out_hash(A);

        if (fps[mask].H > max_H) { max_H = fps[mask].H; max_H_mask = mask; }

        if ((mask & 0xFFFF) == 0) {
            fprintf(stderr, "  %d/%d (%.1f%%)\r", mask, NMASK, 100.0*mask/NMASK);
        }
    }
    fprintf(stderr, "  Done. Max H = %lld at mask %d          \n", max_H, max_H_mask);

    /* Phase 2: Group by fingerprint */
    fprintf(stderr, "Phase 2: Grouping by fingerprint...\n");

    /* Sort masks by fingerprint. Use indices to avoid copying. */
    /* Since we have 2M entries, use a hash map approach */
    /* Simple approach: sort an array of mask indices by fingerprint */

    /* For efficiency, use a hash table */
    /* Hash = combine scores + H + c3 */
    #define HASH_SIZE (1<<20)
    int *hash_next = (int *)malloc(NMASK * sizeof(int)); /* linked list next */
    int *hash_head = (int *)malloc(HASH_SIZE * sizeof(int));
    memset(hash_head, -1, HASH_SIZE * sizeof(int));

    for (int mask = 0; mask < NMASK; mask++) {
        unsigned int h = 0;
        for (int i = 0; i < N; i++) h = h * 31 + fps[mask].scores[i];
        h = h * 31 + (unsigned int)(fps[mask].H & 0xFFFFFFFF);
        h = h * 31 + fps[mask].c3;
        for (int i = 0; i < N; i++) h = h * 31 + fps[mask].v3c[i];
        for (int i = 0; i < N; i++) h = h * 31 + fps[mask].s2[i];
        h = h * 31 + fps[mask].co_hash;
        h &= (HASH_SIZE - 1);

        /* Check if this fingerprint already exists in the chain */
        int found = -1;
        int prev = -1;
        for (int node = hash_head[h]; node != -1; node = hash_next[node]) {
            if (cmp_fp(&fps[mask], &fps[node]) == 0) {
                found = fp_class[node];
                break;
            }
        }

        if (found >= 0) {
            fp_class[mask] = found;
            fp_classes[found].count++;
            if (is_grid_sym(mask)) fp_classes[found].gs_count++;
        } else {
            int ci = n_fp_classes++;
            if (ci >= MAX_FP_CLASSES) {
                fprintf(stderr, "ERROR: too many fingerprint classes (>%d)\n", MAX_FP_CLASSES);
                return 1;
            }
            memcpy(&fp_classes[ci].fp, &fps[mask], sizeof(Fingerprint));
            fp_classes[ci].count = 1;
            fp_classes[ci].gs_count = is_grid_sym(mask) ? 1 : 0;
            fp_classes[ci].first_mask = mask;
            fp_classes[ci].is_palindromic = is_palindromic_scores(fps[mask].scores);
            fp_classes[ci].is_sf = 0;
            fp_classes[ci].is_sc = 0;
            fp_class[mask] = ci;
        }

        /* Add to hash chain */
        hash_next[mask] = hash_head[h];
        hash_head[h] = mask;
    }
    fprintf(stderr, "  %d fingerprint classes\n", n_fp_classes);

    /* Phase 3: Check SF (self-flip) — scan ALL masks */
    fprintf(stderr, "Phase 3: Checking self-flip (scanning all masks)...\n");
    int sf_candidates = 0;
    int sf_confirmed = 0;
    int all_ones = NMASK - 1;

    /* First pass: find all fp classes that have any mask with fp_class[mask] == fp_class[flip(mask)] */
    int *sf_candidate_mask = (int *)calloc(n_fp_classes, sizeof(int)); /* store a witness mask */
    int *sf_candidate_flag = (int *)calloc(n_fp_classes, sizeof(int));

    for (int mask = 0; mask < NMASK; mask++) {
        int flip = mask ^ all_ones;
        if (fp_class[mask] == fp_class[flip] && !sf_candidate_flag[fp_class[mask]]) {
            sf_candidate_flag[fp_class[mask]] = 1;
            sf_candidate_mask[fp_class[mask]] = mask;
            sf_candidates++;
        }
    }
    fprintf(stderr, "  FP classes with SF-candidate masks: %d\n", sf_candidates);

    /* Second pass: for each candidate, verify actual isomorphism */
    for (int ci = 0; ci < n_fp_classes; ci++) {
        if (!sf_candidate_flag[ci]) continue;

        int mask = sf_candidate_mask[ci];
        int flip = mask ^ all_ones;
        int A[N][N], B[N][N];
        mask_to_adj(mask, A);
        mask_to_adj(flip, B);

        if (are_isomorphic(A, B)) {
            fp_classes[ci].is_sf = 1;
            sf_confirmed++;
        } else {
            /* The witness mask failed. Try other masks in this class. */
            /* Scan for more witnesses. */
            int found = 0;
            for (int m2 = mask + 1; m2 < NMASK && !found; m2++) {
                if (fp_class[m2] != ci) continue;
                int f2 = m2 ^ all_ones;
                if (fp_class[f2] != ci) continue;
                mask_to_adj(m2, A);
                mask_to_adj(f2, B);
                if (are_isomorphic(A, B)) {
                    fp_classes[ci].is_sf = 1;
                    sf_confirmed++;
                    found = 1;
                }
            }
        }
    }
    free(sf_candidate_flag);
    free(sf_candidate_mask);
    fprintf(stderr, "  SF candidates: %d, confirmed: %d\n", sf_candidates, sf_confirmed);

    /* Phase 4: Check SC (self-converse) */
    fprintf(stderr, "Phase 4: Checking self-converse...\n");
    int sc_candidates = 0;
    int sc_confirmed = 0;

    for (int ci = 0; ci < n_fp_classes; ci++) {
        if (!fp_classes[ci].is_palindromic) continue;
        sc_candidates++;

        int A[N][N], At[N][N];
        mask_to_adj(fp_classes[ci].first_mask, A);
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                At[i][j] = A[j][i];

        if (are_isomorphic(A, At)) {
            fp_classes[ci].is_sc = 1;
            sc_confirmed++;
        }
    }
    fprintf(stderr, "  SC candidates (palindromic): %d, confirmed: %d\n",
            sc_candidates, sc_confirmed);

    /* ============================================================
       Output results
       ============================================================ */
    printf("n=%d, m=%d, 2^m=%d\n", N, M, NMASK);
    printf("Fingerprint classes: %d\n", n_fp_classes);
    printf("(OEIS A000568 says n=8 should have 6880 iso classes)\n\n");

    int both_count = 0;
    for (int ci = 0; ci < n_fp_classes; ci++)
        if (fp_classes[ci].is_sc && fp_classes[ci].is_sf) both_count++;

    printf("Self-converse (SC): %d\n", sc_confirmed);
    printf("Self-flip (SF): %d\n", sf_confirmed);
    printf("Both SC+SF (kernel): %d\n", both_count);
    printf("Max H: %lld\n\n", max_H);

    /* H-maximizer details */
    printf("--- H-MAXIMIZER ---\n");
    {
        int ci = fp_class[max_H_mask];
        int A[N][N];
        mask_to_adj(max_H_mask, A);
        int aut = count_aut(A);
        int aaut = count_anti_aut(A);

        printf("  Class #%d: H=%lld, size=%d, SC=%d, SF=%d, GS=%d/%d",
               ci, fp_classes[ci].fp.H, fp_classes[ci].count,
               fp_classes[ci].is_sc, fp_classes[ci].is_sf,
               fp_classes[ci].gs_count, fp_classes[ci].count);
        printf(", |Aut|=%d, |AntiAut|=%d", aut, aaut);
        printf(", scores=(");
        for (int i = 0; i < N; i++) printf("%d%s", fp_classes[ci].fp.scores[i], i<N-1?",":"");
        printf(")\n");
        printf("  c3=%d, v3c_profile=(", fp_classes[ci].fp.c3);
        for (int i = 0; i < N; i++) printf("%d%s", fp_classes[ci].fp.v3c[i], i<N-1?",":"");
        printf(")\n");
    }

    /* SC+SF kernel classes */
    printf("\n--- SC+SF KERNEL CLASSES ---\n");
    for (int ci = 0; ci < n_fp_classes; ci++) {
        if (!(fp_classes[ci].is_sc && fp_classes[ci].is_sf)) continue;

        int A[N][N];
        mask_to_adj(fp_classes[ci].first_mask, A);
        int aut = count_aut(A);
        int aaut = count_anti_aut(A);

        int is_hmax = (fp_classes[ci].fp.H == max_H);
        printf("  #%d: H=%lld%s, size=%d, |Aut|=%d, |AntiAut|=%d, GS=%d/%d (%.4f)",
               ci, fp_classes[ci].fp.H, is_hmax ? " *MAX*" : "",
               fp_classes[ci].count, aut, aaut,
               fp_classes[ci].gs_count, fp_classes[ci].count,
               (double)fp_classes[ci].gs_count / fp_classes[ci].count);
        printf(", scores=(");
        for (int i = 0; i < N; i++) printf("%d%s", fp_classes[ci].fp.scores[i], i<N-1?",":"");
        printf(")\n");
    }

    /* SF-only classes */
    printf("\n--- SF-ONLY CLASSES (first 20) ---\n");
    int sf_only_count = 0;
    for (int ci = 0; ci < n_fp_classes; ci++) {
        if (fp_classes[ci].is_sf && !fp_classes[ci].is_sc) {
            sf_only_count++;
            if (sf_only_count <= 20) {
                printf("  #%d: H=%lld, size=%d, GS=%d/%d, scores=(",
                       ci, fp_classes[ci].fp.H, fp_classes[ci].count,
                       fp_classes[ci].gs_count, fp_classes[ci].count);
                for (int i = 0; i < N; i++) printf("%d%s", fp_classes[ci].fp.scores[i], i<N-1?",":"");
                printf(")\n");
            }
        }
    }
    printf("Total SF-only: %d\n", sf_only_count);

    /* SC-only classes */
    int sc_only_count = 0;
    for (int ci = 0; ci < n_fp_classes; ci++) {
        if (fp_classes[ci].is_sc && !fp_classes[ci].is_sf) sc_only_count++;
    }
    printf("\nSC-only: %d\n", sc_only_count);

    /* Top 15 by H */
    printf("\n--- TOP 15 BY H ---\n");
    int top_order[MAX_FP_CLASSES];
    for (int i = 0; i < n_fp_classes; i++) top_order[i] = i;
    for (int i = 0; i < n_fp_classes-1 && i < 15; i++)
        for (int j = i+1; j < n_fp_classes; j++)
            if (fp_classes[top_order[j]].fp.H > fp_classes[top_order[i]].fp.H)
                { int t=top_order[i]; top_order[i]=top_order[j]; top_order[j]=t; }

    for (int k = 0; k < 15 && k < n_fp_classes; k++) {
        int ci = top_order[k];
        printf("  #%d: H=%lld, size=%d, SC=%d, SF=%d, GS=%d/%d, scores=(",
               ci, fp_classes[ci].fp.H, fp_classes[ci].count,
               fp_classes[ci].is_sc, fp_classes[ci].is_sf,
               fp_classes[ci].gs_count, fp_classes[ci].count);
        for (int i = 0; i < N; i++) printf("%d%s", fp_classes[ci].fp.scores[i], i<N-1?",":"");
        printf(")\n");
    }

    /* GS analysis */
    printf("\n--- GS ANALYSIS ---\n");
    int total_gs = 0;
    int gs_class_count = 0;
    for (int ci = 0; ci < n_fp_classes; ci++) {
        total_gs += fp_classes[ci].gs_count;
        if (fp_classes[ci].gs_count > 0) gs_class_count++;
    }
    printf("Total GS tilings: %d (expected 2^12 = %d)\n", total_gs, 1<<12);
    printf("FP classes with GS > 0: %d/%d\n", gs_class_count, n_fp_classes);

    /* GS in kernel */
    printf("\nGS in SC+SF kernel classes:\n");
    for (int ci = 0; ci < n_fp_classes; ci++) {
        if (fp_classes[ci].is_sc && fp_classes[ci].is_sf) {
            printf("  #%d: GS=%d (predicted 3^3=27 per class)\n",
                   ci, fp_classes[ci].gs_count);
        }
    }

    /* H value distribution for SC+SF */
    printf("\n--- SC+SF KERNEL H VALUES ---\n");
    for (int ci = 0; ci < n_fp_classes; ci++) {
        if (fp_classes[ci].is_sc && fp_classes[ci].is_sf) {
            printf("  H=%lld\n", fp_classes[ci].fp.H);
        }
    }

    free(fps);
    free(fp_class);
    free(hash_next);
    free(hash_head);

    return 0;
}
"""

def run_n8_analysis():
    with tempfile.NamedTemporaryFile(mode='w', suffix='.c', delete=False) as f:
        f.write(C_CODE)
        c_path = f.name

    exe_path = c_path.replace('.c', '')
    try:
        result = subprocess.run(['gcc', '-O2', '-o', exe_path, c_path],
                              capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Compilation failed: {result.stderr}")
            return

        print("Compiled successfully. Running n=8 analysis...")
        print("(This may take several minutes for 2^21 = 2,097,152 tilings)")
        t0 = time.time()
        result = subprocess.run([exe_path], capture_output=True, text=True, timeout=3600)
        elapsed = time.time() - t0
        print(result.stdout)
        if result.stderr:
            print(f"\n(stderr: {result.stderr.strip()[-500:]})")
        print(f"\nTotal time: {elapsed:.1f}s")
    finally:
        os.unlink(c_path)
        if os.path.exists(exe_path):
            os.unlink(exe_path)

if __name__ == '__main__':
    run_n8_analysis()
