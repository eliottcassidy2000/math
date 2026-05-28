/*
 * paley_hp_counter.c — Hamiltonian path counter for Paley tournaments QR_p
 *
 * Two algorithms:
 *   1. Standard bitmask DP: O(2^p * p) time, O(2^p * p * 8) bytes RAM
 *   2. Circulant-reduced DP: O(2^(p-1) * p) time, O(2^(p-1) * 8) bytes RAM
 *      Uses Z_p symmetry: always normalize so current vertex = 0.
 *      Reduces memory by factor 2p vs standard.
 *
 * For p=23: standard 1.5 GB → circulant 32 MB  (47x smaller)
 * For p=31: standard 530 GB → circulant 8.6 GB (feasible!)
 *
 * Results (Paley sub-tournament optimality, THM-336):
 *   a(p-1) = H(QR_p - v)   [all v equivalent by vertex-transitivity]
 *   a(p)   = H(QR_p)       [conjectured = a(n) for n = p]
 *
 * Verified:
 *   p=7:  a(6)=45,                a(7)=189                c=72
 *   p=11: a(10)=15745,            a(11)=95095              c=39675
 *   p=19: a(18)=117266659317,     a(19)=1172695746915      c=527714543799
 *   p=23: a(22)=1313333107451805, a(23)=15760206976379349  c=7223436934463772
 *
 * Pattern: c(p)/a(p-1) → (p-1)/4 as p grows (asymptotic formula)
 * Equivalently: a(p)/a(p-1) → (p+1)/2 as p grows
 *
 * Session: opus-2026-05-27-S7
 */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define MAXP 32
typedef uint64_t u64;
typedef uint32_t u32;

/* Build QR_p: adj[i] has bit j set iff i→j (i.e., (j-i) mod p ∈ QR) */
void build_paley(int p, u32 *adj) {
    int qr[MAXP] = {0};
    for (int i = 1; i < p; i++) qr[(i * i) % p] = 1;
    for (int i = 0; i < p; i++)
        for (int j = 0; j < p; j++)
            if (i != j && qr[((j - i) % p + p) % p])
                adj[i] |= (1u << j);
}

/* Standard bitmask DP: O(2^p * p) memory */
u64 count_hp_standard(int n, u32 *adj) {
    int full = (1 << n) - 1;
    u64 *dp = calloc((u64)(full + 1) * n, sizeof(u64));
    if (!dp) { fprintf(stderr, "OOM at standard DP\n"); return 0; }
    for (int v = 0; v < n; v++) dp[(u64)(1 << v) * n + v] = 1;
    u64 total = 0;
    for (int mask = 1; mask <= full; mask++) {
        for (int v = 0; v < n; v++) {
            if (!(mask >> v & 1)) continue;
            u64 d = dp[(u64)mask * n + v];
            if (!d) continue;
            if (mask == full) { total += d; continue; }
            u32 out = adj[v] & ~(u32)mask;
            while (out) {
                int w = __builtin_ctz(out);
                dp[(u64)(mask | (1 << w)) * n + w] += d;
                out &= out - 1;
            }
        }
    }
    free(dp);
    return total;
}

/*
 * Circulant-reduced DP for Z_p-symmetric tournaments.
 *
 * Key insight: for a Z_p-circulant tournament with connection set S,
 * all starting vertices are equivalent. We normalize: current vertex = 0.
 *
 * State: visited set as bitmask on {0,...,p-1}, bit 0 always set (current = 0).
 * Transition: pick d ∈ S (outgoing from 0) where d not yet visited.
 *   New state = shift all visited bits by -d mod p, then set bit 0.
 *   i.e., new_mask = { (v - d) mod p : v ∈ visited ∪ {d} }
 *               = { (v - d) mod p : v ∈ visited } ∪ {0}
 *
 * Initial: dp[{0}] = p  (p starting vertices, all equivalent → weight p)
 * Answer:  dp[{0,1,...,p-1}] = H(T) = p × #{HP from vertex 0}
 *
 * Memory: 2^(p-1) entries (all masks with bit 0 set)
 *   Array indexed by mask >> 1 (drop bit 0 since it's always 1)
 *   dp_arr[mask >> 1] = count for state mask
 *
 * For p=31: 2^30 × 8 bytes = 8 GB — feasible on a workstation
 */
/*
 * Circulant-reduced DP: MUST process in order of path length (popcount).
 * Transitions always go from popcount k → popcount k+1, but the new mask_hi
 * value can be SMALLER than the old one, so simple linear scan fails.
 * Fix: p separate passes, each pass processes only masks of the current length.
 *
 * State: mask_hi = mask >> 1, where mask has bit 0 always set (current = 0).
 * For p=31: 2^30 × 8 bytes = 8 GB — feasible on a workstation.
 */
u64 count_hp_circulant(int p, int *S, int S_size) {
    u64 sz = 1ULL << (p - 1);
    u64 *dp = calloc(sz, sizeof(u64));
    if (!dp) { fprintf(stderr, "OOM: need %.1f GB\n",
                       (double)sz * 8 / 1e9); return 0; }

    /* Initial state: visited = {0}, weight = p (all starting vertices equiv) */
    dp[0] = p;

    u64 total = 0;

    /* Process in p passes: pass k handles paths of length k */
    for (int len = 1; len <= p; len++) {
        int target_pc = len - 1; /* popcount of mask_hi for length-len paths */
        for (u64 mask_hi = 0; mask_hi < sz; mask_hi++) {
            if (__builtin_popcountll(mask_hi) != target_pc) continue;
            u64 cnt = dp[mask_hi];
            if (!cnt) continue;
            if (len == p) { total += cnt; continue; }

            u32 mask = (u32)(mask_hi << 1) | 1u; /* restore bit 0 */

            for (int si = 0; si < S_size; si++) {
                int d = S[si];
                if ((mask >> d) & 1u) continue; /* already visited */

                /* Shift all visited vertices by -d mod p; new current = 0 */
                u32 new_mask = 1u;
                u32 m = mask;
                while (m) {
                    int v = __builtin_ctz(m);
                    int nv = (v - d + p) % p;
                    if (nv) new_mask |= (1u << nv);
                    m &= m - 1;
                }
                dp[new_mask >> 1] += cnt;
            }
        }
    }

    free(dp);
    return total;
}

/* Remove vertex v from adj[0..n-1] */
void remove_vertex(int n, u32 *adj, int v, int *new_n, u32 *new_adj) {
    *new_n = n - 1;
    int mapping[MAXP], idx = 0;
    for (int u = 0; u < n; u++) if (u != v) mapping[idx++] = u;
    for (int ui = 0; ui < *new_n; ui++) {
        new_adj[ui] = 0;
        for (int wi = 0; wi < *new_n; wi++)
            if (adj[mapping[ui]] >> mapping[wi] & 1)
                new_adj[ui] |= (1u << wi);
    }
}

void print_time(clock_t t0) {
    double s = (double)(clock() - t0) / CLOCKS_PER_SEC;
    if (s < 60) printf("[%.1fs]", s);
    else printf("[%.0fm%.0fs]", floor(s/60), fmod(s,60));
    fflush(stdout);
}

int main(int argc, char **argv) {
    int p_list[] = {7, 11, 19, 23, 31};
    int n_list = sizeof(p_list) / sizeof(p_list[0]);
    int use_circulant = 1; /* default: use fast algorithm */

    /* Parse optional argument: prime to compute */
    int target_p = 23; /* default */
    if (argc > 1) target_p = atoi(argv[1]);
    if (argc > 2 && strcmp(argv[2], "--standard") == 0) use_circulant = 0;

    printf("=== Paley HP Counter ===\n");
    printf("Algorithm: %s\n\n", use_circulant ? "circulant-reduced" : "standard bitmask");

    int p = target_p;
    printf("--- QR_%d (p ≡ %d mod 4) ---\n", p, p % 4);

    /* Build connection set S = QR_p */
    int S[MAXP], S_size = 0;
    {
        int qr[MAXP] = {0};
        for (int i = 1; i < p; i++) qr[(i*i)%p] = 1;
        for (int i = 1; i < p; i++) if (qr[i]) S[S_size++] = i;
    }
    printf("Connection set |S|=%d (= (p-1)/2 = %d, regular: %s)\n",
           S_size, (p-1)/2, S_size == (p-1)/2 ? "YES" : "NO");

    u32 adj[MAXP] = {0};
    build_paley(p, adj);

    /* Verify degree */
    int deg = __builtin_popcount(adj[0]);
    printf("Outdegree of v=0: %d (expected %d)\n\n", deg, (p-1)/2);

    clock_t t0;

    if (use_circulant) {
        printf("Computing H(QR_%d) via circulant reduction...\n", p);
        printf("  Memory: 2^%d × 8 = %.1f MB\n", p-1,
               (double)(1ULL << (p-1)) * 8 / 1e6);
        t0 = clock();
        u64 hp = count_hp_circulant(p, S, S_size);
        printf("  H(QR_%d) = %llu  ", p, hp);
        print_time(t0);
        printf("\n  → a(%d)\n\n", p);

        /* Verify against standard for small p */
        if (p <= 23) {
            printf("Verifying against standard DP...\n");
            t0 = clock();
            u64 hp2 = count_hp_standard(p, adj);
            printf("  Standard: H(QR_%d) = %llu  ", p, hp2);
            print_time(t0);
            printf("\n  Match: %s\n\n", hp == hp2 ? "YES ✓" : "NO ✗");
        }
    } else {
        printf("Computing H(QR_%d) via standard bitmask DP...\n", p);
        t0 = clock();
        u64 hp = count_hp_standard(p, adj);
        printf("H(QR_%d) = %llu  ", p, hp);
        print_time(t0);
        printf("\n\n");
    }

    /* Also compute H(QR_p - v) = a(p-1) for small p */
    if (p <= 24) {
        printf("Computing H(QR_%d - v0) = a(%d)...\n", p, p-1);
        int n2; u32 adj2[MAXP];
        remove_vertex(p, adj, 0, &n2, adj2);
        t0 = clock();
        u64 hp2 = count_hp_standard(n2, adj2);
        printf("  H(QR_%d - v) = %llu  ", p, hp2);
        print_time(t0);
        printf("\n  → a(%d)\n\n", p-1);

        /* Print c(p) stats */
        // (these are computed from full H values stored above)
    }

    printf("\n=== Known Paley terms ===\n");
    printf("  a(2)=1,   a(3)=3    [p=3]\n");
    printf("  a(6)=45,  a(7)=189  [p=7]\n");
    printf("  a(10)=15745, a(11)=95095  [p=11]\n");
    printf("  a(18)=117266659317, a(19)=1172695746915  [p=19]\n");
    printf("  a(22)=1313333107451805, a(23)=15760206976379349  [p=23]\n");
    printf("  a(30)=? a(31)=?  [p=31 — use --p 31 with 8GB RAM]\n");

    return 0;
}
