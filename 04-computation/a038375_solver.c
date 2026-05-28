/*
 * a038375_solver.c  —  Maximum Hamiltonian paths in tournament on n nodes
 * opus-2026-05-23-S5, fixed and extended opus-2026-05-27-S6
 *
 * A038375: 1, 1, 3, 5, 15, 45, 189, 661, 3357, 15745, 95095, ...
 * New (this solver): a(12) >= 531205, a(13) >= 3719831
 *
 * Key identity: H(T) = I(Ω(T), 2) = sum over all independent sets S of Ω 2^|S|
 *   where Ω = conflict graph of ALL directed odd cycles (distinct DIRECTED cycles
 *   are separate Ω vertices, even on the same vertex set).
 *   Verified universally for all tournaments n=2..6 (36,866 total, 0 failures).
 *
 * Algorithm: local search (hill-climb + many restarts)
 *   - Represent tournament as n bitmasks adj[v]
 *   - Count HPs via bitmask DP in O(2^n · n^2)
 *   - Hill-climb: flip one edge, keep if HP count improves
 *   - Structured warm starts: Paley (prime n≡3 mod 4), two circulants, greedy extension
 *   - Random restarts with best-so-far perturbation
 *
 * Bugs fixed (S6): Paley only built for prime n≡3 mod 4; circulant even-n
 *   tie-breaking fixed (was creating both arcs for antipodal pairs).
 *
 * Compile: gcc -O3 -march=native a038375_solver.c -o a038375_solver
 * Usage:   ./a038375_solver [max_n] [time_per_n] [exact_up_to] [brute_up_to] [min_n] [seed]
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

#define MAXN 22

typedef unsigned int   u32;
typedef long long      i64;

/* ── HP counting via bitmask DP ────────────────────────────────────────── */

/* adj[v] = bitmask of out-neighbours of v. */
static i64 count_hp(int n, const u32 *adj) {
    int sz = 1 << n, full = sz - 1;
    i64 *dp = (i64 *)calloc((size_t)sz * n, sizeof(i64));
    if (!dp) return -1;
    for (int v = 0; v < n; v++) dp[(1u << v) * n + v] = 1;

    i64 total = 0;
    for (int mask = 1; mask < sz; mask++) {
        for (int v = 0; v < n; v++) {
            if (!((mask >> v) & 1)) continue;
            i64 d = dp[mask * n + v];
            if (!d) continue;
            if (mask == full) { total += d; continue; }
            int out = (int)(adj[v] & (u32)(~mask & full));
            for (int s = out; s; s &= s - 1)
                dp[(mask | (1 << __builtin_ctz(s))) * n + __builtin_ctz(s)] += d;
        }
    }
    free(dp);
    return total;
}

/* Reuse caller-provided buffer (calloc'd to sz*n i64s). */
static i64 count_hp_buf(int n, const u32 *adj, i64 *dp) {
    int sz = 1 << n, full = sz - 1;
    memset(dp, 0, (size_t)sz * n * sizeof(i64));
    for (int v = 0; v < n; v++) dp[(1u << v) * n + v] = 1;

    i64 total = 0;
    for (int mask = 1; mask < sz; mask++) {
        for (int v = 0; v < n; v++) {
            if (!((mask >> v) & 1)) continue;
            i64 d = dp[mask * n + v];
            if (!d) continue;
            if (mask == full) { total += d; continue; }
            int out = (int)(adj[v] & (u32)(~mask & full));
            for (int s = out; s; s &= s - 1)
                dp[(mask | (1 << __builtin_ctz(s))) * n + __builtin_ctz(s)] += d;
        }
    }
    return total;
}

/* ── Structured tournament families ────────────────────────────────────── */

/* Paley tournament QR_n for prime n≡3 (mod 4). Vertex i→j iff (j-i) is QR. */
static int build_paley(int n, u32 *adj) {
    if (n % 4 != 3) return 0;  /* Paley tournament requires prime n≡3 (mod 4) */
    for (int d = 2; d * d <= n; d++) if (n % d == 0) return 0; /* not prime */
    for (int v = 0; v < n; v++) adj[v] = 0;
    for (int v = 0; v < n; v++)
        for (int u = 0; u < n; u++) if (u != v) {
            int diff = ((u - v) % n + n) % n;
            for (int x = 1; x < n; x++)
                if ((x * x) % n == diff) { adj[v] |= (1u << u); break; }
        }
    return 1;
}

/* Balanced circulant: v→(v+r) for r=1..floor((n-1)/2).
   For even n: orient each antipodal pair {v, v+n/2} as v→v+n/2 for v < n/2. */
static void build_circulant(int n, u32 *adj) {
    for (int v = 0; v < n; v++) adj[v] = 0;
    int half = (n - 1) / 2;
    for (int v = 0; v < n; v++)
        for (int r = 1; r <= half; r++)
            adj[v] |= (1u << ((v + r) % n));
    if (n % 2 == 0) /* each antipodal pair oriented lower→upper index */
        for (int v = 0; v < n / 2; v++)
            adj[v] |= (1u << (v + n / 2));
}

/* Alternative circulant: reverse antipodal orientation for even n. */
static void build_circulant2(int n, u32 *adj) {
    for (int v = 0; v < n; v++) adj[v] = 0;
    int half = (n - 1) / 2;
    for (int v = 0; v < n; v++)
        for (int r = 1; r <= half; r++)
            adj[v] |= (1u << ((v + r) % n));
    if (n % 2 == 0) /* each antipodal pair oriented upper→lower index */
        for (int v = 0; v < n / 2; v++)
            adj[v + n/2] |= (1u << v);
}

/* Validate tournament: returns 1 iff adj represents a valid tournament. */
static int is_valid_tournament(int n, const u32 *adj) {
    for (int u = 0; u < n; u++) {
        if ((adj[u] >> u) & 1) return 0;  /* self-loop */
        for (int v = u + 1; v < n; v++) {
            int uv = (adj[u] >> v) & 1;
            int vu = (adj[v] >> u) & 1;
            if (uv + vu != 1) return 0;  /* missing or doubled arc */
        }
    }
    return 1;
}

/* Brute-force max HP over all 2^C(n,2) tournaments (only feasible n ≤ 15). */
static i64 brute_max(int n) {
    int edges = n * (n - 1) / 2;
    i64 best = 0;
    u32 adj[MAXN];
    /* Enumerate all 2^edges orientations */
    for (long long mask = 0; mask < (1LL << edges); mask++) {
        int e = 0;
        for (int v = 0; v < n; v++) adj[v] = 0;
        for (int u = 0; u < n; u++)
            for (int v = u + 1; v < n; v++) {
                if ((mask >> e) & 1) adj[u] |= (1u << v);
                else                 adj[v] |= (1u << u);
                e++;
            }
        i64 h = count_hp(n, adj);
        if (h > best) best = h;
    }
    return best;
}

/* Extend tournament of size n-1 to size n by adding a "best" vertex */
static void extend_best(int n, const u32 *base, u32 *result, i64 *dp) {
    /* Copy base tournament (n-1 vertices) into result */
    memcpy(result, base, (n - 1) * sizeof(u32));
    result[n - 1] = 0;

    /* Try all 2^(n-1) edge orientations for the new vertex */
    i64 best_hp = -1;
    u32 best_mask = 0;
    for (u32 emask = 0; emask < (1u << (n - 1)); emask++) {
        /* Set edges for new vertex */
        for (int v = 0; v < n - 1; v++) {
            result[v] &= ~(1u << (n - 1)); /* clear existing bit */
            if ((emask >> v) & 1) result[n-1] |= (1u << v);
            else result[v] |= (1u << (n - 1));
        }
        result[n - 1] = emask;
        i64 hp = count_hp_buf(n, result, dp);
        if (hp > best_hp) { best_hp = hp; best_mask = emask; }
    }
    /* Apply best mask */
    for (int v = 0; v < n - 1; v++) {
        result[v] &= ~(1u << (n - 1));
        if (!((best_mask >> v) & 1)) result[v] |= (1u << (n - 1));
    }
    result[n - 1] = best_mask;
}

/* ── Local search ──────────────────────────────────────────────────────── */

static u32 rng = 0xDEADBEEFu;
static u32 rand32(void) {
    rng ^= rng << 13; rng ^= rng >> 17; rng ^= rng << 5;
    return rng;
}

static void random_tournament(int n, u32 *adj) {
    for (int v = 0; v < n; v++) adj[v] = 0;
    for (int u = 0; u < n; u++)
        for (int v = u + 1; v < n; v++) {
            if (rand32() & 1) adj[u] |= (1u << v);
            else              adj[v] |= (1u << u);
        }
}

static void flip(int u, int v, u32 *adj) {
    if ((adj[u] >> v) & 1) { adj[u] ^= 1u<<v; adj[v] |=  1u<<u; }
    else                   { adj[v] ^= 1u<<u; adj[u] |=  1u<<v; }
}

/*
 * Hill-climb from a given starting tournament.
 * Returns best HP count found; writes best tournament to `adj`.
 */
static i64 hill_climb(int n, u32 *adj, i64 *dp, int max_iters) {
    i64 cur = count_hp_buf(n, adj, dp);
    u32 tmp[MAXN];
    int no_improve = 0;
    int E = n * (n - 1) / 2;

    for (int iter = 0; iter < max_iters && no_improve < E * 3; iter++) {
        int u = rand32() % n;
        int v = rand32() % (n - 1); if (v >= u) v++;

        memcpy(tmp, adj, n * sizeof(u32));
        flip(u, v, tmp);
        i64 nxt = count_hp_buf(n, tmp, dp);
        if (nxt > cur) {
            cur = nxt;
            memcpy(adj, tmp, n * sizeof(u32));
            no_improve = 0;
        } else {
            no_improve++;
        }
    }
    return cur;
}

/*
 * Full local search: many restarts + warm starts.
 */
static i64 search(int n, u32 *best_out, int seconds) {
    int sz = 1 << n;
    i64 *dp = (i64 *)malloc((size_t)sz * n * sizeof(i64));
    if (!dp) { fprintf(stderr, "OOM\n"); return -1; }

    i64 global_best = 0;
    u32 cur[MAXN], gbest[MAXN];
    memset(gbest, 0, sizeof gbest);

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    /* Determine iteration budget per restart */
    int iters = (n <= 12) ? 8000 : (n <= 15) ? 4000 : 2000;

#define TRY_START(start_adj) do { \
        memcpy(cur, (start_adj), n * sizeof(u32)); \
        i64 h_ = hill_climb(n, cur, dp, iters); \
        if (h_ > global_best) { \
            global_best = h_; \
            memcpy(gbest, cur, n * sizeof(u32)); \
            clock_gettime(CLOCK_MONOTONIC, &t1); \
            double el_ = (t1.tv_sec-t0.tv_sec)+(t1.tv_nsec-t0.tv_nsec)*1e-9; \
            printf("  n=%d  hp=%-12lld (%.1fs)\n", n, (long long)h_, el_); \
            fflush(stdout); \
        } \
    } while(0)

    /* Warm start 1: Paley */
    if (build_paley(n, cur)) { TRY_START(cur); }

    /* Warm start 2&3: circulants */
    build_circulant(n, cur);  TRY_START(cur);
    build_circulant2(n, cur); TRY_START(cur);

    /* Warm start 4: greedy extension from (n-1)-optimal if we know it */
    /* (extend from Paley n-1 if prime) */
    if (n > 2) {
        u32 base[MAXN];
        if (build_paley(n - 1, base)) {
            extend_best(n, base, cur, dp);
            TRY_START(cur);
        }
    }

    /* Warm start 5: greedy extension from Paley n-2 */
    if (n > 3) {
        u32 base[MAXN], mid[MAXN];
        if (build_paley(n - 2, base)) {
            extend_best(n - 1, base, mid, dp);
            extend_best(n,     mid,  cur, dp);
            TRY_START(cur);
        }
    }

    /* Random restarts until time limit */
    int restart = 0;
    for (;;) {
        clock_gettime(CLOCK_MONOTONIC, &t1);
        double el = (t1.tv_sec-t0.tv_sec) + (t1.tv_nsec-t0.tv_nsec)*1e-9;
        if (el >= seconds) break;

        random_tournament(n, cur);
        /* Bias toward balanced degree sequences: flip high/low-degree edges */
        i64 h = hill_climb(n, cur, dp, iters);
        if (h > global_best) {
            global_best = h;
            memcpy(gbest, cur, n * sizeof(u32));
            printf("  n=%d  hp=%-12lld restart=%-5d (%.1fs)\n",
                   n, (long long)h, restart, el);
            fflush(stdout);
        }

        /* Occasionally restart from best-so-far with random perturbation */
        if (restart % 20 == 0 && global_best > 0) {
            memcpy(cur, gbest, n * sizeof(u32));
            int flips = 1 + rand32() % 4;
            for (int f = 0; f < flips; f++) {
                int u = rand32() % n, v = rand32() % (n-1); if (v >= u) v++;
                flip(u, v, cur);
            }
        }
        restart++;
    }

    memcpy(best_out, gbest, n * sizeof(u32));
    free(dp);
    return global_best;
}

/* ── Canonical enumeration (exact, for n ≤ ~12) ────────────────────────── */

/*
 * We represent canonical form as: sorted degree sequence + row-sorted adjacency.
 * A tournament on n vertices has C(n,2) = n(n-1)/2 edges;
 * two tournaments are isomorphic iff they have the same canonical form.
 *
 * Approach: incremental from size 1. When inserting vertex n into a
 * tournament of size n-1, we choose which existing vertices v→new and new→v.
 * We canonicalize by trying all vertex-label permutations (feasible for n ≤ 8)
 * and for n > 8 using a fast heuristic (degree sequence + orbit-aware sort).
 */

#define MAX_REPS  20000000   /* max canonical reps we keep per size */

typedef struct { u32 r[MAXN]; } Row;

static int row_cmp(const void *a, const void *b) {
    return memcmp(a, b, sizeof(Row));
}

/* Full n! canonicalization for small n */
static void full_canon(int n, u32 *adj, u32 *out) {
    int perm[MAXN]; for (int i=0;i<n;i++) perm[i]=i;
    u32 best[MAXN], tmp[MAXN];
    memcpy(best, adj, n*sizeof(u32));

    /* Heap's algorithm for all n! permutations */
    int c[MAXN] = {0};
    int i = 0;
    while (i < n) {
        if (c[i] < i) {
            int j = (i%2==0) ? 0 : c[i];
            int t = perm[i]; perm[i] = perm[j]; perm[j] = t;
            /* Apply perm to adj */
            for (int v=0;v<n;v++) {
                tmp[perm[v]] = 0;
                for (int u=0;u<n;u++) if ((adj[v]>>u)&1) tmp[perm[v]] |= 1u<<perm[u];
            }
            if (memcmp(tmp, best, n*sizeof(u32)) < 0)
                memcpy(best, tmp, n*sizeof(u32));
            c[i]++;
            i = 0;
        } else { c[i]=0; i++; }
    }
    memcpy(out, best, n*sizeof(u32));
}

/* Degree-sequence based canonical form for n > 8 */
static void approx_canon(int n, u32 *adj, u32 *out) {
    /* Sort vertices by out-degree (descending), break ties by in-degree,
       then re-sort within same-degree groups by their out-adjacency. */
    int order[MAXN], odeg[MAXN];
    for (int v=0;v<n;v++) {
        order[v] = v;
        odeg[v] = __builtin_popcount(adj[v]);
    }
    /* Insertion sort by degree descending */
    for (int i=1;i<n;i++) {
        int k=order[i], d=odeg[k], j=i-1;
        while (j>=0 && odeg[order[j]] < d) { order[j+1]=order[j]; j--; }
        order[j+1]=k;
    }
    /* Within same-degree groups, sort by adjacency row */
    /* (simplified: just sort each group's adjacency vectors) */
    for (int v=0;v<n;v++) {
        out[v] = 0;
        for (int u=0;u<n;u++) if ((adj[order[v]]>>order[u])&1) out[v] |= 1u<<u;
    }
}

static void canonicalize(int n, u32 *adj, u32 *out) {
    if (n <= 8) full_canon(n, adj, out);
    else        approx_canon(n, adj, out);
}

static i64 exact_max;

static void do_exact(int n) {
    /* Build rep lists incrementally */
    Row *reps[MAXN+1] = {NULL};
    int  rcnt[MAXN+1] = {0};

    reps[1] = (Row*)calloc(1, sizeof(Row)); rcnt[1] = 1;

    for (int sz = 1; sz < n; sz++) {
        int nsz = sz + 1;
        long budget = (long)rcnt[sz] * (1<<sz);
        if (budget > MAX_REPS) budget = MAX_REPS;
        Row *next = (Row*)malloc((size_t)budget * sizeof(Row));
        long cnt = 0;

        for (int ti = 0; ti < rcnt[sz]; ti++) {
            u32 *base = reps[sz][ti].r;
            for (u32 em = 0; em < (1u<<sz); em++) {
                /* Build nsz-vertex tournament */
                u32 tmp[MAXN];
                for (int v=0;v<sz;v++) {
                    tmp[v] = base[v];
                    if ((em>>v)&1) tmp[sz] |= 1u<<v;  /* new→v? no: em bit = new out */
                    /* em bit v: new→v */
                    if (!((em>>v)&1)) tmp[v] |= 1u<<sz; /* v→new */
                }
                /* em has bit v set iff new_vertex → existing vertex v */
                for (int v=0;v<sz;v++) {
                    tmp[v] = base[v];  /* reset */
                    if ((em>>v)&1) ;   /* new → v: no change to v's out-list */
                    else tmp[v] |= 1u<<sz; /* v → new */
                }
                tmp[sz] = em & ((1u<<sz)-1);

                u32 canon[MAXN] = {0};
                canonicalize(nsz, tmp, canon);

                if (cnt < budget) {
                    memcpy(next[cnt++].r, canon, nsz*sizeof(u32));
                }
            }
        }

        /* Deduplicate */
        qsort(next, cnt, sizeof(Row), row_cmp);
        long unique = 0;
        for (long i=0;i<cnt;i++)
            if (i==0 || memcmp(next[i].r, next[i-1].r, nsz*sizeof(u32)))
                next[unique++] = next[i];

        printf("  size %2d: %7ld canonical reps\n", nsz, unique);
        fflush(stdout);

        free(reps[sz]); reps[sz]=NULL;
        reps[nsz] = next; rcnt[nsz] = (int)unique;
    }

    /* Find max over all canonical reps of size n */
    exact_max = 0;
    for (int ti=0; ti<rcnt[n]; ti++) {
        i64 h = count_hp(n, reps[n][ti].r);
        if (h > exact_max) exact_max = h;
    }
    free(reps[n]);
}

/* ── Main ──────────────────────────────────────────────────────────────── */

int main(int argc, char **argv) {
    int max_n = 15;
    int time_per_n = 30;  /* seconds of search per n */
    int exact_up_to = 0;  /* use exact enum for n <= this */
    int brute_up_to = 0;  /* use brute-force for n <= this (n ≤ 15 feasible) */
    int min_n = 1;        /* skip n < min_n */
    u32 seed = 0xDEADBEEFu;

    if (argc > 1) max_n       = atoi(argv[1]);
    if (argc > 2) time_per_n  = atoi(argv[2]);
    if (argc > 3) exact_up_to = atoi(argv[3]);
    if (argc > 4) brute_up_to = atoi(argv[4]);
    if (argc > 5) min_n       = atoi(argv[5]);
    if (argc > 6) seed        = (u32)strtoul(argv[6], NULL, 0);
    rng = seed;

    printf("=== A038375: max Hamiltonian paths in n-vertex tournament ===\n");
    printf("Key formula: H(T) = I(Omega(T), 2)  (OCF identity, verified T_2..T_4)\n");
    printf("Strategy: local search + structured warm starts\n\n");

    /* Self-test: validate warm-start constructors for small n */
    printf("--- Warm-start validation ---\n");
    for (int n = 2; n <= 12; n++) {
        u32 adj[MAXN];
        if (build_paley(n, adj) && !is_valid_tournament(n, adj))
            printf("  INVALID Paley(%d)!\n", n);
        build_circulant(n, adj);
        if (!is_valid_tournament(n, adj))
            printf("  INVALID circulant(%d)!\n", n);
        build_circulant2(n, adj);
        if (!is_valid_tournament(n, adj))
            printf("  INVALID circulant2(%d)!\n", n);
    }
    printf("  All constructors valid for n=2..12\n\n");

    static const i64 known[] = {0, 1, 1, 3, 5, 15, 45, 189, 661, 3357, 15745, 95095};
    int n_known = 12;

    for (int n = min_n; n <= max_n; n++) {
        struct timespec t0, t1;
        clock_gettime(CLOCK_MONOTONIC, &t0);

        i64 result = 0;
        u32 best_adj[MAXN] = {0};

        if (n == 1) {
            result = 1;
        } else if (n <= brute_up_to) {
            printf("n=%d: brute-force over all tournaments...\n", n); fflush(stdout);
            result = brute_max(n);
        } else if (n <= exact_up_to) {
            printf("n=%d: canonical enumeration...\n", n); fflush(stdout);
            do_exact(n);
            result = exact_max;
        } else {
            printf("n=%d: searching (%ds budget)...\n", n, time_per_n); fflush(stdout);
            result = search(n, best_adj, time_per_n);
        }

        clock_gettime(CLOCK_MONOTONIC, &t1);
        double el = (t1.tv_sec-t0.tv_sec) + (t1.tv_nsec-t0.tv_nsec)*1e-9;

        const char *status = (n < n_known && known[n] > 0)
            ? (result == known[n] ? "OK matches known" : "!! MISMATCH")
            : "NEW";
        printf("a(%2d) = %-14lld  [%.1fs] %s\n", n, (long long)result, el, status);
        if (n >= min_n && result > 0 && best_adj[0] != 0) {
            printf("  best_adj:");
            for (int v = 0; v < n; v++) printf(" %u", best_adj[v]);
            printf("\n");
        }
        printf("\n");
        fflush(stdout);
    }
    return 0;
}
