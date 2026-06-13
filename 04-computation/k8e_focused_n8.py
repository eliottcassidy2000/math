#!/usr/bin/env python3
"""
K_8-e analysis at n=8: sample tournaments with alpha_1=8 and check i_2.

For H=21 via K_8-e: need alpha_1=8, i_2=1.
Exhaustive check: what i_2 values occur with alpha_1=8?

Instance: opus-2026-03-07-S40
"""

import ctypes
import os
import time
import random
from collections import Counter

C_CODE = r"""
#include <string.h>

long long held_karp_8(int adj[8]) {
    int n = 8;
    int full = (1 << n) - 1;
    long long dp[256][8];
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

/* Count ALL directed odd cycles: 3-cycles, 5-cycles, 7-cycles.
   Returns alpha_1 = total count.
   Also stores vertex sets and cycle lengths. */
typedef struct {
    int verts;     /* bitmask of vertices */
    int length;    /* 3, 5, or 7 */
} CycleInfo;

int find_all_odd_cycles(int adj[8], int n, CycleInfo cycles[], int max_cycles) {
    int count = 0;

    /* 3-cycles */
    for (int i = 0; i < n; i++)
        for (int j = i+1; j < n; j++)
            for (int k = j+1; k < n; k++) {
                if (((adj[i]>>j)&1) && ((adj[j]>>k)&1) && ((adj[k]>>i)&1)) {
                    if (count < max_cycles) {
                        cycles[count].verts = (1<<i)|(1<<j)|(1<<k);
                        cycles[count].length = 3;
                    }
                    count++;
                } else if (((adj[i]>>k)&1) && ((adj[k]>>j)&1) && ((adj[j]>>i)&1)) {
                    if (count < max_cycles) {
                        cycles[count].verts = (1<<i)|(1<<j)|(1<<k);
                        cycles[count].length = 3;
                    }
                    count++;
                }
            }

    /* 5-cycles: Hamiltonian cycles on each 5-subset */
    int v5[5];
    for (v5[0]=0; v5[0]<n-4; v5[0]++)
    for (v5[1]=v5[0]+1; v5[1]<n-3; v5[1]++)
    for (v5[2]=v5[1]+1; v5[2]<n-2; v5[2]++)
    for (v5[3]=v5[2]+1; v5[3]<n-1; v5[3]++)
    for (v5[4]=v5[3]+1; v5[4]<n; v5[4]++) {
        long long dp5[32][5];
        memset(dp5, 0, sizeof(dp5));
        dp5[1][0] = 1;
        for (int S = 1; S < 32; S++)
            for (int i = 0; i < 5; i++) {
                if (!(S & (1<<i))) continue;
                if (dp5[S][i] == 0) continue;
                for (int j = 0; j < 5; j++) {
                    if (S & (1<<j)) continue;
                    if ((adj[v5[i]] >> v5[j]) & 1)
                        dp5[S|(1<<j)][j] += dp5[S][i];
                }
            }
        int vmask = 0;
        for (int i = 0; i < 5; i++) vmask |= (1 << v5[i]);
        for (int j = 1; j < 5; j++) {
            if (dp5[31][j] > 0 && ((adj[v5[j]] >> v5[0]) & 1)) {
                int c5 = (int)dp5[31][j];
                for (int r = 0; r < c5; r++) {
                    if (count < max_cycles) {
                        cycles[count].verts = vmask;
                        cycles[count].length = 5;
                    }
                    count++;
                }
            }
        }
    }

    /* 7-cycles: Hamiltonian cycles on each 7-subset */
    int v7[7];
    for (v7[0]=0; v7[0]<n-6; v7[0]++)
    for (v7[1]=v7[0]+1; v7[1]<n-5; v7[1]++)
    for (v7[2]=v7[1]+1; v7[2]<n-4; v7[2]++)
    for (v7[3]=v7[2]+1; v7[3]<n-3; v7[3]++)
    for (v7[4]=v7[3]+1; v7[4]<n-2; v7[4]++)
    for (v7[5]=v7[4]+1; v7[5]<n-1; v7[5]++)
    for (v7[6]=v7[5]+1; v7[6]<n; v7[6]++) {
        long long dp7[128][7];
        memset(dp7, 0, sizeof(dp7));
        dp7[1][0] = 1;
        for (int S = 1; S < 128; S++)
            for (int i = 0; i < 7; i++) {
                if (!(S & (1<<i))) continue;
                if (dp7[S][i] == 0) continue;
                for (int j = 0; j < 7; j++) {
                    if (S & (1<<j)) continue;
                    if ((adj[v7[i]] >> v7[j]) & 1)
                        dp7[S|(1<<j)][j] += dp7[S][i];
                }
            }
        int vmask = 0;
        for (int i = 0; i < 7; i++) vmask |= (1 << v7[i]);
        for (int j = 1; j < 7; j++) {
            if (dp7[127][j] > 0 && ((adj[v7[j]] >> v7[0]) & 1)) {
                int c7 = (int)dp7[127][j];
                for (int r = 0; r < c7; r++) {
                    if (count < max_cycles) {
                        cycles[count].verts = vmask;
                        cycles[count].length = 7;
                    }
                    count++;
                }
            }
        }
    }

    return count;
}

/* For one tournament, compute alpha_1, i_2, cycle composition, H */
void analyze_one(long long bits, int *alpha1, int *i2, int *t3, int *t5, int *t7, long long *H) {
    int n = 8, m = 28;
    int edges[28][2];
    int ei = 0;
    for (int i = 0; i < n; i++)
        for (int j = i+1; j < n; j++) {
            edges[ei][0] = i;
            edges[ei][1] = j;
            ei++;
        }
    int adj[8] = {0};
    for (int k = 0; k < m; k++) {
        int i = edges[k][0], j = edges[k][1];
        if ((bits >> k) & 1)
            adj[j] |= (1 << i);
        else
            adj[i] |= (1 << j);
    }

    CycleInfo cycles[500];
    int max_cycles = 500;
    int a1 = find_all_odd_cycles(adj, n, cycles, max_cycles);
    *alpha1 = a1;

    /* Count cycle types */
    int c3 = 0, c5 = 0, c7 = 0;
    for (int i = 0; i < a1 && i < max_cycles; i++) {
        if (cycles[i].length == 3) c3++;
        else if (cycles[i].length == 5) c5++;
        else c7++;
    }
    *t3 = c3; *t5 = c5; *t7 = c7;

    /* Count i_2 = number of vertex-disjoint pairs */
    int cnt_i2 = 0;
    int actual = (a1 < max_cycles) ? a1 : max_cycles;
    for (int a = 0; a < actual; a++)
        for (int b = a+1; b < actual; b++)
            if (!(cycles[a].verts & cycles[b].verts))
                cnt_i2++;
    *i2 = cnt_i2;

    *H = held_karp_8(adj);
}
"""


def main():
    c_path = "/tmp/k8e_focused.c"
    so_path = "/tmp/k8e_focused.so"
    with open(c_path, "w") as f:
        f.write(C_CODE)
    os.system(f"cc -O3 -o {so_path} -shared -fPIC {c_path}")
    lib = ctypes.CDLL(so_path)

    random.seed(42)
    num_samples = 500000
    m = 28

    alpha1 = ctypes.c_int()
    i2 = ctypes.c_int()
    t3 = ctypes.c_int()
    t5 = ctypes.c_int()
    t7 = ctypes.c_int()
    H = ctypes.c_longlong()

    # For each alpha_1 in {4,6,8,10}, collect (i_2, H) pairs
    from collections import defaultdict
    alpha_data = defaultdict(list)

    print(f"Sampling {num_samples} tournaments at n=8...")
    start = time.time()

    for trial in range(num_samples):
        if trial % 100000 == 0 and trial > 0:
            elapsed = time.time() - start
            print(f"  {trial}/{num_samples} ({elapsed:.1f}s)")

        bits = random.randint(0, (1 << m) - 1)
        lib.analyze_one(bits, ctypes.byref(alpha1), ctypes.byref(i2),
                       ctypes.byref(t3), ctypes.byref(t5), ctypes.byref(t7),
                       ctypes.byref(H))

        a1 = alpha1.value
        if a1 in [4, 6, 8, 10]:
            alpha_data[a1].append((i2.value, H.value, t3.value, t5.value, t7.value))

    elapsed = time.time() - start
    print(f"Done in {elapsed:.1f}s\n")

    # Analyze each target alpha_1
    for a1 in [4, 6, 8, 10]:
        data = alpha_data[a1]
        if not data:
            print(f"alpha_1={a1}: no samples")
            continue

        print("=" * 60)
        needed_i2 = (10 - a1) // 2
        print(f"alpha_1={a1}: {len(data)} samples (need i_2={needed_i2} for H=21)")
        print("=" * 60)

        i2_dist = Counter(i2 for i2, _, _, _, _ in data)
        h_dist = Counter(h for _, h, _, _, _ in data)
        combo_dist = Counter((t3, t5, t7) for _, _, t3, t5, t7 in data)

        print(f"  i_2 distribution: {dict(sorted(i2_dist.items()))}")
        print(f"  H distribution: {dict(sorted(h_dist.items()))}")
        print(f"  (t3,t5,t7) combos: {dict(sorted(combo_dist.items()))}")

        if needed_i2 in i2_dist:
            # Look at H values when i_2 = needed value
            matching = [(h, t3, t5, t7) for i2, h, t3, t5, t7 in data if i2 == needed_i2]
            h_when_match = Counter(h for h, _, _, _ in matching)
            print(f"  When i_2={needed_i2}: H distribution = {dict(sorted(h_when_match.items()))}")
            if 21 in h_when_match:
                print(f"  *** H=21 FOUND with alpha_1={a1}, i_2={needed_i2}! ***")
            else:
                print(f"  H=21 NOT found even when i_2={needed_i2}")
        else:
            print(f"  i_2={needed_i2} NEVER occurs! H=21 blocked.")

        print()


if __name__ == "__main__":
    main()
