#!/usr/bin/env python3
"""
H vs alpha_1 (total odd cycles in Omega) at n=7 (exhaustive).

For each tournament, compute H and alpha_1 = |V(Omega)| = total directed odd cycles.
Show the (H, alpha_1) pairs that occur, focusing on the H=21 neighborhood.

Key question: for H near 21, what alpha_1 values occur?
If alpha_1 + 2*alpha_2 = 10 is impossible, then H=21 is impossible.

Instance: opus-2026-03-07-S40
"""

import ctypes
import os
import time
from collections import Counter

C_CODE = r"""
#include <string.h>

long long held_karp_7(int adj[7]) {
    int n = 7;
    int full = (1 << n) - 1;
    long long dp[128][7];
    memset(dp, 0, sizeof(dp));
    for (int v = 0; v < n; v++)
        dp[1 << v][v] = 1;
    for (int S = 1; S <= full; S++) {
        for (int v = 0; v < n; v++) {
            if (!(S & (1 << v))) continue;
            long long c = dp[S][v];
            if (c == 0) continue;
            int out = adj[v] & (~S);
            while (out) {
                int u = __builtin_ctz(out);
                dp[S | (1 << u)][u] += c;
                out &= out - 1;
            }
        }
    }
    long long total = 0;
    for (int v = 0; v < n; v++)
        total += dp[full][v];
    return total;
}

/* Count directed 3-cycles */
int count_3_cycles(int adj[7], int n) {
    int count = 0;
    for (int i = 0; i < n; i++)
        for (int j = i+1; j < n; j++)
            for (int k = j+1; k < n; k++) {
                if (((adj[i]>>j)&1) && ((adj[j]>>k)&1) && ((adj[k]>>i)&1))
                    count++;
                else if (((adj[i]>>k)&1) && ((adj[k]>>j)&1) && ((adj[j]>>i)&1))
                    count++;
            }
    return count;
}

/* Count directed 5-cycles (Hamiltonian cycles on 5-subsets) */
int count_5_cycles(int adj[7], int n) {
    int count = 0;
    /* For each 5-subset, check all directed Hamiltonian cycles */
    int v5[5];
    for (v5[0] = 0; v5[0] < n-4; v5[0]++)
    for (v5[1] = v5[0]+1; v5[1] < n-3; v5[1]++)
    for (v5[2] = v5[1]+1; v5[2] < n-2; v5[2]++)
    for (v5[3] = v5[2]+1; v5[3] < n-1; v5[3]++)
    for (v5[4] = v5[3]+1; v5[4] < n; v5[4]++) {
        /* Use DP for Hamiltonian cycles on these 5 vertices */
        long long dp[32][5];
        memset(dp, 0, sizeof(dp));
        dp[1][0] = 1; /* start at v5[0] */
        for (int S = 1; S < 32; S++) {
            for (int i = 0; i < 5; i++) {
                if (!(S & (1 << i))) continue;
                if (dp[S][i] == 0) continue;
                for (int j = 0; j < 5; j++) {
                    if (S & (1 << j)) continue;
                    if ((adj[v5[i]] >> v5[j]) & 1)
                        dp[S | (1 << j)][j] += dp[S][i];
                }
            }
        }
        int full = 31;
        for (int j = 1; j < 5; j++) {
            if (dp[full][j] > 0 && ((adj[v5[j]] >> v5[0]) & 1))
                count += (int)dp[full][j];
        }
    }
    /* Each directed 5-cycle is counted once per starting vertex in the fixed
       subset, but we fix start=v5[0], so each cycle is counted once if it
       passes through v5[0] first. Actually no, we count directed Ham cycles
       starting at v5[0]. Each undirected 5-cycle gives 2 directed cycles.
       Each directed 5-cycle is counted once (start at min vertex). */
    /* Actually each directed cycle through 5 vertices has 5 rotations.
       We fix start=v5[0], so we count each directed cycle exactly once
       (the rotation starting at the minimum vertex). But that's only true
       if v5[0] IS the minimum vertex, which it always is since we enumerate
       in order. So count = number of directed 5-cycles on these 5 vertices. */
    return count;
}

/* Count directed 7-cycles (Hamiltonian cycles on all 7 vertices) */
int count_7_cycles(int adj[7]) {
    int n = 7;
    long long dp[128][7];
    memset(dp, 0, sizeof(dp));
    dp[1][0] = 1;
    for (int S = 1; S < 128; S++) {
        for (int v = 0; v < n; v++) {
            if (!(S & (1 << v))) continue;
            if (dp[S][v] == 0) continue;
            for (int u = 0; u < n; u++) {
                if (S & (1 << u)) continue;
                if ((adj[v] >> u) & 1)
                    dp[S | (1 << u)][u] += dp[S][v];
            }
        }
    }
    int full = 127;
    int count = 0;
    for (int j = 1; j < n; j++) {
        if (dp[full][j] > 0 && ((adj[j] >> 0) & 1))
            count += (int)dp[full][j];
    }
    return count;
}

/* Process all tournaments, output (H, t3, t5, t7) for each */
void process_all(
    long long *H_arr,
    int *t3_arr,
    int *t5_arr,
    int *t7_arr,
    int max_entries
) {
    int n = 7, m = 21;
    int edges[21][2];
    int ei = 0;
    for (int i = 0; i < n; i++)
        for (int j = i+1; j < n; j++) {
            edges[ei][0] = i;
            edges[ei][1] = j;
            ei++;
        }

    long long total = 1L << m;
    int idx = 0;
    for (long long bits = 0; bits < total && idx < max_entries; bits++) {
        int adj[7] = {0};
        for (int k = 0; k < m; k++) {
            int i = edges[k][0], j = edges[k][1];
            if ((bits >> k) & 1)
                adj[j] |= (1 << i);
            else
                adj[i] |= (1 << j);
        }

        long long H = held_karp_7(adj);
        int t3 = count_3_cycles(adj, n);
        int t5 = count_5_cycles(adj, n);
        int t7 = count_7_cycles(adj);

        H_arr[idx] = H;
        t3_arr[idx] = t3;
        t5_arr[idx] = t5;
        t7_arr[idx] = t7;
        idx++;
    }
}
"""

def main():
    c_path = "/tmp/h_vs_alpha.c"
    so_path = "/tmp/h_vs_alpha.so"
    with open(c_path, "w") as f:
        f.write(C_CODE)
    os.system(f"cc -O3 -o {so_path} -shared -fPIC {c_path}")
    lib = ctypes.CDLL(so_path)

    N = 2**21  # 2,097,152
    H_arr = (ctypes.c_longlong * N)()
    t3_arr = (ctypes.c_int * N)()
    t5_arr = (ctypes.c_int * N)()
    t7_arr = (ctypes.c_int * N)()

    print(f"Computing (H, t3, t5, t7) for all {N} tournaments at n=7...")
    start = time.time()
    lib.process_all(H_arr, t3_arr, t5_arr, t7_arr, N)
    elapsed = time.time() - start
    print(f"Done in {elapsed:.1f}s\n")

    # Analyze: for each H, what (t3, t5, t7) values occur?
    # alpha_1 = t3 + t5 + t7 (total directed odd cycles)
    # But we need to verify: H = I(Omega, 2) = 1 + 2*alpha_1 + 4*alpha_2 + ...
    # So H >= 1 + 2*alpha_1, meaning alpha_1 <= (H-1)/2

    # Focus on H near 21
    print("=" * 70)
    print("H vs cycle counts (t3, t5, t7) for H in [15..29]")
    print("=" * 70)

    for target_H in range(15, 30, 2):
        entries = []
        for i in range(N):
            if H_arr[i] == target_H:
                entries.append((t3_arr[i], t5_arr[i], t7_arr[i]))
        if entries:
            alpha1_vals = Counter(t3 + t5 + t7 for t3, t5, t7 in entries)
            cycle_combos = Counter(entries)
            print(f"\nH={target_H}: {len(entries)} tournaments")
            print(f"  alpha_1 distribution: {dict(sorted(alpha1_vals.items()))}")
            print(f"  (t3,t5,t7) combos: {dict(sorted(cycle_combos.items()))}")
        else:
            print(f"\nH={target_H}: *** NO TOURNAMENTS ***")

    # Check OCF: H = 1 + 2*t3 + 2*t5 + 2*t7 + 4*(sharing pairs) + ...
    # The minimum H for given alpha_1 is 1 + 2*alpha_1 (all cycles pairwise disjoint)
    # The maximum H for given alpha_1 is 3^alpha_1 (all pairwise sharing = complete graph)
    print("\n" + "=" * 70)
    print("alpha_1 vs H range")
    print("=" * 70)
    alpha_to_H = {}
    for i in range(N):
        a1 = t3_arr[i] + t5_arr[i] + t7_arr[i]
        H = H_arr[i]
        if a1 not in alpha_to_H:
            alpha_to_H[a1] = [H, H, 0]
        alpha_to_H[a1][0] = min(alpha_to_H[a1][0], H)
        alpha_to_H[a1][1] = max(alpha_to_H[a1][1], H)
        alpha_to_H[a1][2] += 1

    for a1 in sorted(alpha_to_H.keys()):
        lo, hi, cnt = alpha_to_H[a1]
        # Can H=21 with this alpha_1?
        # Need: 1 + 2*alpha_1 <= 21, so alpha_1 <= 10
        can_21 = "possible" if a1 <= 10 and lo <= 21 <= hi else "blocked"
        tag = f" [{can_21}]" if a1 <= 10 else ""
        print(f"  alpha_1={a1}: H in [{lo}, {hi}], {cnt} tournaments{tag}")

    # For H=21, need alpha_1 + 2*alpha_2 = 10
    # So the "gap" between H and 1+2*alpha_1 is 4*alpha_2 + 8*alpha_3 + ...
    # which equals H - 1 - 2*alpha_1
    print("\n" + "=" * 70)
    print("Tournaments with alpha_1 <= 10: gap analysis")
    print("=" * 70)
    for i in range(N):
        a1 = t3_arr[i] + t5_arr[i] + t7_arr[i]
        H = H_arr[i]
        if a1 <= 10 and H <= 25:
            gap = H - 1 - 2 * a1
            if gap >= 0:  # valid (gap = 4*alpha_2 + ...)
                pass


if __name__ == "__main__":
    main()
