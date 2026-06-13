/*
 * circulant_search.c — Find the circulant tournament on Z_n maximizing H.
 *
 * For ODD n: a circulant tournament on Z_n is defined by a connection set
 * S ⊆ {1,...,n-1} satisfying: i ∈ S ↔ (n-i) ∉ S (antisymmetry), |S| = (n-1)/2.
 * There are 2^{(n-1)/2} such connection sets (one from each complementary pair).
 *
 * We use the circulant-reduced bitmask DP (from paley_hp_counter.c):
 *   Memory: 2^{n-1} × 8 bytes.
 *   Time: n × 2^{n-1} per tournament evaluation.
 *
 * Strategy: exhaustive search for small n; simulated annealing / local search for larger n.
 * For n=25: 2^12 = 4096 choices. With ~2s each: ~8000s. Too slow for exhaustive search.
 *   Use local search starting from Paley-like structures.
 *
 * Memory requirements:
 *   n=25: 2^24 × 8 = 128 MB  ✓ easy
 *   n=27: 2^26 × 8 = 512 MB  ✓ easy
 *   n=29: 2^28 × 8 = 2 GB    ✓ feasible
 *   n=31: 2^30 × 8 = 8 GB    ✓ feasible (large machine)
 *
 * Session: opus-2026-05-27-S7
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

typedef uint64_t u64;
typedef uint32_t u32;

static u64 *dp_arr = NULL;
static u64 dp_sz = 0;

/* Allocate the DP array (reuse if big enough) */
int alloc_dp(int n) {
    u64 need = 1ULL << (n - 1);
    if (need > dp_sz) {
        free(dp_arr);
        dp_arr = calloc(need, sizeof(u64));
        if (!dp_arr) {
            fprintf(stderr, "OOM: %.1f GB\n", (double)need * 8 / 1e9);
            return 0;
        }
        dp_sz = need;
    }
    return 1;
}

/*
 * Count HP using circulant-reduced DP.
 * S[i] for i=0..S_size-1: outgoing connection set from vertex 0.
 * (vertices j where j ∈ S means 0→j)
 */
u64 count_hp_circulant(int p, int *S, int S_size) {
    u64 sz = 1ULL << (p - 1);
    memset(dp_arr, 0, sz * sizeof(u64));
    dp_arr[0] = p;  /* initial: visited={0}, weight=p */

    u64 total = 0;
    for (int len = 1; len <= p; len++) {
        int target_pc = len - 1;
        for (u64 mask_hi = 0; mask_hi < sz; mask_hi++) {
            if (__builtin_popcountll(mask_hi) != target_pc) continue;
            u64 cnt = dp_arr[mask_hi];
            if (!cnt) continue;
            if (len == p) { total += cnt; continue; }
            u32 mask = (u32)(mask_hi << 1) | 1u;
            for (int si = 0; si < S_size; si++) {
                int d = S[si];
                if ((mask >> d) & 1u) continue;
                u32 new_mask = 1u;
                u32 m = mask;
                while (m) {
                    int v = __builtin_ctz(m);
                    int nv = (v - d + p) % p;
                    if (nv) new_mask |= (1u << nv);
                    m &= m - 1;
                }
                dp_arr[new_mask >> 1] += cnt;
            }
        }
    }
    return total;
}

/* Build Paley connection set for prime p */
int build_paley_S(int p, int *S) {
    int qr[128] = {0};
    for (int i = 1; i < p; i++) qr[(i*i)%p] = 1;
    int sz = 0;
    for (int i = 1; i < p; i++) if (qr[i]) S[sz++] = i;
    return sz;
}

/*
 * For non-prime n: build the "balanced" connection set {1, 2, ..., (n-1)/2}.
 * This is the "lower half" circulant — not optimal, but a good baseline.
 */
int build_lower_half_S(int n, int *S) {
    int sz = 0;
    for (int i = 1; i <= (n-1)/2; i++) S[sz++] = i;
    return sz;
}

/*
 * For non-prime n: try to build a "Paley-like" connection set using
 * quadratic residues mod n (treating n as if it were prime, ignoring
 * the fact that it's composite). Only valid when gcd(i, n) = 1 for all i.
 */
int build_pseudo_paley_S(int n, int *S) {
    /* For n not prime: use QR mod n for the multiplicative group */
    /* Fall back to lower_half if n is even or composite has bad structure */
    int sz = 0;
    /* Compute QR via squaring (elements with gcd=1 only) */
    int qr[128] = {0};
    for (int i = 1; i < n; i++) {
        /* Only consider i with gcd(i,n)=1 */
        int g = i, h = n;
        while (h) { int t = h; h = g % h; g = t; }
        if (g == 1) qr[(i*i)%n] = 1;
    }
    for (int i = 1; i < n; i++)
        if (qr[i] && !qr[(n-i)%n]) S[sz++] = i;
        else if (qr[i] && qr[(n-i)%n] && i < n-i) S[sz++] = i;  /* both QR: pick smaller */
    /* Pad to exactly (n-1)/2 if needed */
    /* This may not always work cleanly for composite n */
    return sz;
}

/* Local search: flip one pair to maximize H */
u64 local_search_circulant(int n, int *S, int *S_size) {
    u64 best = count_hp_circulant(n, S, *S_size);
    int improved = 1;
    while (improved) {
        improved = 0;
        for (int i = 1; i <= (n-1)/2; i++) {
            /* Check if i is in S */
            int in_S = 0, pos = -1;
            for (int j = 0; j < *S_size; j++) if (S[j] == i) { in_S = 1; pos = j; break; }
            /* Swap i with (n-i): if i ∈ S, try removing i, adding n-i; vice versa */
            int S2[64];
            memcpy(S2, S, *S_size * sizeof(int));
            if (in_S) {
                /* Remove i, add n-i */
                S2[pos] = n - i;
            } else {
                /* Find n-i in S and remove it, add i */
                for (int j = 0; j < *S_size; j++) {
                    if (S2[j] == n-i) { S2[j] = i; break; }
                }
            }
            u64 h = count_hp_circulant(n, S2, *S_size);
            if (h > best) {
                best = h;
                memcpy(S, S2, *S_size * sizeof(int));
                improved = 1;
                printf("  swap %d↔%d → H=%llu\n", i, n-i, best);
                break;  /* restart */
            }
        }
    }
    return best;
}

void print_S(int *S, int sz, int n) {
    printf("S={");
    for (int i = 0; i < sz; i++) printf("%d%s", S[i], i+1<sz?",":"");
    printf("} (complement: {");
    for (int i = 0; i < sz; i++) printf("%d%s", n-S[i], i+1<sz?",":"");
    printf("})");
}

int cmp_int(const void *a, const void *b) { return *(int*)a - *(int*)b; }

int main(int argc, char **argv) {
    int n_list[] = {13, 15, 17, 19, 21, 23, 25, 27, 29};
    int n_count = sizeof(n_list) / sizeof(n_list[0]);

    int target_n = 25;
    if (argc > 1) target_n = atoi(argv[1]);

    printf("=== Circulant Tournament HP Search, n=%d ===\n", target_n);
    printf("Memory: 2^%d × 8 = %.1f MB\n\n", target_n-1,
           (double)(1ULL << (target_n-1)) * 8 / 1e6);

    if (!alloc_dp(target_n)) return 1;

    int n = target_n;
    int S[64], S_size;

    /* 1. Try lower half baseline */
    S_size = build_lower_half_S(n, S);
    clock_t t0 = clock();
    u64 h_lower = count_hp_circulant(n, S, S_size);
    double elapsed = (double)(clock()-t0)/CLOCKS_PER_SEC;
    printf("Lower-half circulant ({1..%d}): H=%llu  [%.2fs]\n", (n-1)/2, h_lower, elapsed);

    /* 2. Try Paley if n is prime */
    int is_prime = 1;
    for (int i = 2; i * i <= n; i++) if (n % i == 0) { is_prime = 0; break; }
    if (is_prime && n % 4 == 3) {
        int S_paley[64];
        int sp = build_paley_S(n, S_paley);
        t0 = clock();
        u64 h_paley = count_hp_circulant(n, S_paley, sp);
        elapsed = (double)(clock()-t0)/CLOCKS_PER_SEC;
        printf("Paley QR_%d: H=%llu  [%.2fs]\n", n, h_paley, elapsed);
        if (h_paley > h_lower) { memcpy(S, S_paley, sp * sizeof(int)); S_size = sp; }
    }

    /* 3. Try a few specific connection sets for non-prime n */
    /* For n=25=5²: try "Paley-like" based on Z_25 */
    if (!is_prime) {
        /* Try several hand-crafted S sets */
        /* Strategy: for n=p², use the "spread" of QR_{p^2} */
        int candidates[8][32];
        int cand_sizes[8];
        int n_cands = 0;

        /* Candidate: upper half {(n+1)/2, ..., n-1} - not valid since |upper|=(n-1)/2 */
        /* Candidate: alternating {1,3,5,...} vs {2,4,6,...} */
        {
            int sc[32]; int szc = 0;
            for (int i = 1; i <= (n-1)/2; i += 2) sc[szc++] = i; /* odds in lower half */
            for (int i = 2; i <= (n-1)/2; i += 2) sc[szc++] = n-i; /* evens from upper half */
            if (szc == (n-1)/2) { memcpy(candidates[n_cands], sc, szc*4); cand_sizes[n_cands++]=szc; }
        }
        /* Candidate: {2,4,...,n-1} — even lower half + odd upper half */
        {
            int sc[32]; int szc = 0;
            for (int i = 2; i <= (n-1)/2; i += 2) sc[szc++] = i;
            for (int i = 1; i <= (n-1)/2; i += 2) sc[szc++] = n-i;
            if (szc == (n-1)/2) { memcpy(candidates[n_cands], sc, szc*4); cand_sizes[n_cands++]=szc; }
        }

        u64 best = h_lower;
        for (int ci = 0; ci < n_cands; ci++) {
            t0 = clock();
            u64 h = count_hp_circulant(n, candidates[ci], cand_sizes[ci]);
            elapsed = (double)(clock()-t0)/CLOCKS_PER_SEC;
            printf("Candidate %d: H=%llu  [%.2fs]\n", ci+1, h, elapsed);
            if (h > best) { best = h; memcpy(S, candidates[ci], cand_sizes[ci]*4); S_size=cand_sizes[ci]; }
        }
    }

    /* 4. Local search from best starting point */
    printf("\nLocal search from best found (H=%llu):\n", count_hp_circulant(n, S, S_size));
    u64 h_best = local_search_circulant(n, S, &S_size);

    qsort(S, S_size, sizeof(int), cmp_int);
    printf("\nBest circulant on Z_%d: H=%llu\n", n, h_best);
    print_S(S, S_size, n);
    printf("\n");
    printf("Lower bound for a(%d) ≥ %llu\n", n, h_best);

    free(dp_arr);
    return 0;
}
/* This section intentionally left blank -- see below for random restart version */
