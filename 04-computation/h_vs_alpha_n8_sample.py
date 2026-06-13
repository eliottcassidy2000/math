#!/usr/bin/env python3
"""
H vs alpha_1 at n=8 (sampling) - check if H=21 decompositions are blocked.

At n=7 (exhaustive):
  alpha_1=10 -> H=29 always (i_2=2, blocks (10,0) decomposition)
  alpha_1=8  -> H=17 always (i_2=0, blocks (8,1) decomposition)
  alpha_1=6  -> H in {13,17} (i_2 in {0,1}, blocks (6,2) decomposition)

Does this persist at n=8?

Instance: opus-2026-03-07-S40
"""

import ctypes
import os
import time
import random
from collections import Counter, defaultdict

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

/* Count directed 3-cycles */
int count_3_cycles(int adj[], int n) {
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

/* Count directed 5-cycles on each 5-subset */
int count_5_cycles(int adj[], int n) {
    int count = 0;
    int v5[5];
    for (v5[0]=0; v5[0]<n-4; v5[0]++)
    for (v5[1]=v5[0]+1; v5[1]<n-3; v5[1]++)
    for (v5[2]=v5[1]+1; v5[2]<n-2; v5[2]++)
    for (v5[3]=v5[2]+1; v5[3]<n-1; v5[3]++)
    for (v5[4]=v5[3]+1; v5[4]<n; v5[4]++) {
        long long dp[32][5];
        memset(dp, 0, sizeof(dp));
        dp[1][0] = 1;
        for (int S = 1; S < 32; S++)
            for (int i = 0; i < 5; i++) {
                if (!(S & (1<<i))) continue;
                if (dp[S][i] == 0) continue;
                for (int j = 0; j < 5; j++) {
                    if (S & (1<<j)) continue;
                    if ((adj[v5[i]] >> v5[j]) & 1)
                        dp[S|(1<<j)][j] += dp[S][i];
                }
            }
        for (int j = 1; j < 5; j++)
            if (dp[31][j] > 0 && ((adj[v5[j]] >> v5[0]) & 1))
                count += (int)dp[31][j];
    }
    return count;
}

/* Count directed 7-cycles on each 7-subset */
int count_7_cycles(int adj[], int n) {
    int count = 0;
    int v7[7];
    for (v7[0]=0; v7[0]<n-6; v7[0]++)
    for (v7[1]=v7[0]+1; v7[1]<n-5; v7[1]++)
    for (v7[2]=v7[1]+1; v7[2]<n-4; v7[2]++)
    for (v7[3]=v7[2]+1; v7[3]<n-3; v7[3]++)
    for (v7[4]=v7[3]+1; v7[4]<n-2; v7[4]++)
    for (v7[5]=v7[4]+1; v7[5]<n-1; v7[5]++)
    for (v7[6]=v7[5]+1; v7[6]<n; v7[6]++) {
        long long dp[128][7];
        memset(dp, 0, sizeof(dp));
        dp[1][0] = 1;
        for (int S = 1; S < 128; S++)
            for (int i = 0; i < 7; i++) {
                if (!(S & (1<<i))) continue;
                if (dp[S][i] == 0) continue;
                for (int j = 0; j < 7; j++) {
                    if (S & (1<<j)) continue;
                    if ((adj[v7[i]] >> v7[j]) & 1)
                        dp[S|(1<<j)][j] += dp[S][i];
                }
            }
        for (int j = 1; j < 7; j++)
            if (dp[127][j] > 0 && ((adj[v7[j]] >> v7[0]) & 1))
                count += (int)dp[127][j];
    }
    return count;
}

/* Process one tournament. Returns H, t3, t5, t7. */
void process_one(long long bits, long long *H, int *t3, int *t5, int *t7) {
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
    *H = held_karp_8(adj);
    *t3 = count_3_cycles(adj, n);
    *t5 = count_5_cycles(adj, n);
    *t7 = count_7_cycles(adj, n);
}
"""

def main():
    c_path = "/tmp/h_vs_alpha_n8.c"
    so_path = "/tmp/h_vs_alpha_n8.so"
    with open(c_path, "w") as f:
        f.write(C_CODE)
    os.system(f"cc -O3 -o {so_path} -shared -fPIC {c_path}")
    lib = ctypes.CDLL(so_path)

    random.seed(42)
    # Fewer samples since 7-cycle counting is expensive (C(8,7)=8 subsets of size 7)
    num_samples = 200000
    m = 28

    alpha_to_H = defaultdict(lambda: [float('inf'), 0, 0])  # min, max, count
    h_near_21 = []

    H_out = ctypes.c_longlong()
    t3_out = ctypes.c_int()
    t5_out = ctypes.c_int()
    t7_out = ctypes.c_int()

    print(f"Sampling {num_samples} tournaments at n=8...")
    start = time.time()

    for trial in range(num_samples):
        if trial % 50000 == 0 and trial > 0:
            elapsed = time.time() - start
            print(f"  {trial}/{num_samples} ({elapsed:.1f}s)")

        bits = random.randint(0, (1 << m) - 1)
        lib.process_one(bits, ctypes.byref(H_out), ctypes.byref(t3_out),
                       ctypes.byref(t5_out), ctypes.byref(t7_out))

        H = H_out.value
        t3 = t3_out.value
        t5 = t5_out.value
        t7 = t7_out.value
        a1 = t3 + t5 + t7

        rec = alpha_to_H[a1]
        rec[0] = min(rec[0], H)
        rec[1] = max(rec[1], H)
        rec[2] += 1

        if 15 <= H <= 29:
            h_near_21.append((H, t3, t5, t7, a1))

    elapsed = time.time() - start
    print(f"Done in {elapsed:.1f}s\n")

    print("=" * 70)
    print("alpha_1 vs H range at n=8 (sampling)")
    print("=" * 70)
    for a1 in sorted(alpha_to_H.keys()):
        lo, hi, cnt = alpha_to_H[a1]
        tag = ""
        if a1 <= 10:
            # Check if H=21 is in the observed range
            if lo <= 21 <= hi:
                tag = " *** H=21 IN RANGE ***"
            elif hi < 21:
                tag = f" [max H={hi} < 21]"
            elif lo > 21:
                tag = f" [min H={lo} > 21]"
        print(f"  alpha_1={a1}: H in [{lo}, {hi}], {cnt} samples{tag}")

    # Focus on H near 21
    print("\n" + "=" * 70)
    print("Tournaments with H near 21")
    print("=" * 70)
    for target_H in range(15, 30, 2):
        matches = [(t3, t5, t7, a1) for (H, t3, t5, t7, a1) in h_near_21 if H == target_H]
        if matches:
            alpha_dist = Counter(a1 for _, _, _, a1 in matches)
            combo_dist = Counter((t3, t5, t7) for t3, t5, t7, _ in matches)
            print(f"\n  H={target_H}: {len(matches)} tournaments")
            print(f"    alpha_1 dist: {dict(sorted(alpha_dist.items()))}")
            print(f"    (t3,t5,t7) top combos: {dict(sorted(combo_dist.most_common(10)))}")
        else:
            print(f"\n  H={target_H}: *** NOT SEEN ***")

    # Key decomposition check
    print("\n" + "=" * 70)
    print("H=21 DECOMPOSITION CHECK")
    print("=" * 70)
    # For H=21, need alpha_1 + 2*i_2 = 10 with i_3 = 0
    # Possible: (alpha_1, i_2) = (10,0), (8,1), (6,2), (4,3)
    for target_a1 in [4, 6, 8, 10]:
        needed_i2 = (10 - target_a1) // 2
        needed_H = 1 + 2*target_a1 + 4*needed_i2
        rec = alpha_to_H.get(target_a1, None)
        if rec:
            lo, hi, cnt = rec
            if needed_H < lo:
                status = f"BLOCKED (need H={needed_H} but min={lo})"
            elif needed_H > hi:
                status = f"BLOCKED (need H={needed_H} but max={hi})"
            elif needed_H >= lo and needed_H <= hi:
                status = f"*** POSSIBLY ACHIEVABLE *** (H={needed_H} in [{lo},{hi}])"
            else:
                status = "UNKNOWN"
            print(f"  ({target_a1}, {needed_i2}): need H={needed_H}. Observed H in [{lo},{hi}] ({cnt} samples). {status}")
        else:
            print(f"  ({target_a1}, {needed_i2}): alpha_1={target_a1} not seen in sample")


if __name__ == "__main__":
    main()
