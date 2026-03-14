"""
alpha_mod7_constraint.py -- kind-pasteur-2026-03-14-S67

EXHAUSTIVE n=6: alpha_1 + 2*alpha_2 mod 7 NEVER equals 3.
This means H mod 7 = 0 is impossible.

Questions:
1. Is this a consequence of the 3-cycle count c3 mod 7 structure?
2. At n=6, c3 can be 0,1,2,3,4,5,6 -- what values of c3 occur?
3. Directed cycles: each 3-vertex subset has 0 or 1 directed 3-cycles.
   So alpha_1 = c3 + c5 + c7 (weighted by directed cycle count per subset).
   At n=6: only 3-cycles and 5-cycles exist (no 7-cycles).
   3-vertex subset: 0 or 1 directed 3-cycle.
   5-vertex subset: 0, 1, 2, or 3 directed Hamiltonian cycles.

   Actually wait: count_directed_hamcycles counts directed Hamiltonian cycles
   on a subset. For a 3-vertex tournament, there is exactly 1 or 0 directed 3-cycles
   (depending on whether the tournament is cyclic or transitive on those 3 vertices).
   For a 5-vertex tournament on 5 vertices, we count directed Hamiltonian cycles
   (directed 5-cycles). This can be 0, 1, 2, or 3.

   So alpha_1 = sum of (directed Hamiltonian cycle count) over all odd-sized subsets.
   At n=6: alpha_1 = (c3-directed) + (c5-directed)
   where c3-directed = number of 3-vertex subsets with a directed 3-cycle (= c3)
   and c5-directed = sum over 5-vertex subsets of (number of directed Hamiltonian 5-cycles)

Wait, c3 counts the number of 3-vertex subsets that ARE cyclic tournaments
(i.e., have a directed 3-cycle). Each such subset has EXACTLY 1 directed 3-cycle
(in one direction). But wait, a directed 3-cycle has 2 orientations.
Actually, a cyclic tournament 0->1->2->0 has EXACTLY ONE Hamiltonian cycle
(the cycle 0->1->2->0). As a directed Hamiltonian cycle, there is exactly 1
(we count the cycle starting from vertex 0 and going through all vertices back to 0).

Hmm, but count_directed_hamcycles uses a DP that counts Hamiltonian PATHS
from vertex 0 that return to vertex 0. For 3 vertices {0,1,2}:
- If 0->1->2->0 exists: path 0->1->2, return 2->0. Count = 1.
- If 0->2->1->0 exists: path 0->2->1, return 1->0. Count = 1.
- A cyclic tournament has EXACTLY ONE of these (not both).
So each cyclic 3-vertex subset contributes exactly 1 to alpha_1.

For 5-vertex subsets, the function counts all directed Hamiltonian cycles through
vertex 0. So the count can be more than 1.

Let me verify: what are the possible directed Hamiltonian cycle counts on 5 vertices?
"""

import numpy as np
from itertools import combinations
from collections import Counter

def count_directed_hamcycles(A, vertices):
    k = len(vertices)
    if k < 3 or k % 2 == 0:
        return 0
    vlist = list(vertices)
    sub = [[0]*k for _ in range(k)]
    for i in range(k):
        for j in range(k):
            if i != j:
                sub[i][j] = int(A[vlist[i]][vlist[j]])
    full = (1 << k) - 1
    dp = [[0]*k for _ in range(1 << k)]
    dp[1][0] = 1
    for mask in range(1, 1 << k):
        for v in range(k):
            if dp[mask][v] == 0:
                continue
            for u in range(1, k):
                if mask & (1 << u):
                    continue
                if sub[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    total = 0
    for v in range(1, k):
        if dp[full][v] and sub[v][0]:
            total += dp[full][v]
    return total

def main():
    print("=" * 70)
    print("ALPHA MOD 7 CONSTRAINT — STRUCTURAL ANALYSIS")
    print("=" * 70)

    # 1. Exhaustive 5-vertex Hamiltonian cycle counts
    print("\n--- DIRECTED HAMILTONIAN CYCLE COUNTS ON 5-VERTEX TOURNAMENTS ---")
    n = 5
    hc_counter = Counter()
    for bits in range(1 << 10):  # C(5,2) = 10 edges
        A = np.zeros((n, n), dtype=np.int8)
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                if bits & (1 << idx):
                    A[i][j] = 1
                else:
                    A[j][i] = 1
                idx += 1
        cnt = count_directed_hamcycles(A, list(range(n)))
        hc_counter[cnt] += 1
    print(f"  Possible directed Hamiltonian cycle counts: {sorted(hc_counter.keys())}")
    for cnt in sorted(hc_counter.keys()):
        print(f"    count={cnt}: {hc_counter[cnt]} tournaments")

    # 2. Exhaustive n=5: what are alpha_1 values?
    print(f"\n--- ALPHA ANALYSIS AT n=5 (exhaustive) ---")
    n = 5
    a1_counter = Counter()
    a2_counter = Counter()
    pair_counter = Counter()
    H_counter = Counter()

    for bits in range(1 << 10):
        A = np.zeros((n, n), dtype=np.int8)
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                if bits & (1 << idx):
                    A[i][j] = 1
                else:
                    A[j][i] = 1
                idx += 1

        # Count 3-cycles
        c3_directed = 0
        for subset in combinations(range(n), 3):
            cnt = count_directed_hamcycles(A, list(subset))
            c3_directed += cnt

        # Count 5-cycles
        c5_directed = count_directed_hamcycles(A, list(range(n)))

        alpha_1 = c3_directed + c5_directed
        # At n=5, max IS size = 1 (can't fit two disjoint odd cycles in 5 vertices)
        # Wait: 3+3=6>5. So alpha_2=0 at n=5.
        alpha_2 = 0
        H = 1 + 2*alpha_1

        a1_counter[alpha_1] += 1
        H_counter[H] += 1
        pair_counter[(c3_directed, c5_directed)] += 1

    print(f"  Achievable alpha_1: {sorted(a1_counter.keys())}")
    print(f"  Achievable H: {sorted(H_counter.keys())}")
    print(f"\n  (c3_directed, c5_directed) -> alpha_1:")
    for (c3d, c5d), cnt in sorted(pair_counter.items()):
        print(f"    c3d={c3d}, c5d={c5d}: alpha_1={c3d+c5d}, H={1+2*(c3d+c5d)}, count={cnt}")

    # 3. Exhaustive n=6: alpha structure and mod 7
    print(f"\n--- ALPHA STRUCTURE AT n=6 (exhaustive) ---")
    n = 6
    # Track c3, c5 separately
    detail_counter = Counter()

    for bits in range(1 << 15):
        A = np.zeros((n, n), dtype=np.int8)
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                if bits & (1 << idx):
                    A[i][j] = 1
                else:
                    A[j][i] = 1
                idx += 1

        # Count directed 3-cycles
        c3d = 0
        for subset in combinations(range(n), 3):
            c3d += count_directed_hamcycles(A, list(subset))

        # Count directed 5-cycles
        c5d = 0
        for subset in combinations(range(n), 5):
            c5d += count_directed_hamcycles(A, list(subset))

        # Disjoint pairs (alpha_2)
        cycles = []
        for size in range(3, n+1, 2):
            for subset in combinations(range(n), size):
                cnt = count_directed_hamcycles(A, list(subset))
                if cnt > 0:
                    cycles.append((frozenset(subset), cnt, size))

        alpha_1 = sum(cnt for _, cnt, _ in cycles)
        alpha_2 = 0
        for i in range(len(cycles)):
            for j in range(i+1, len(cycles)):
                if len(cycles[i][0] & cycles[j][0]) == 0:
                    alpha_2 += cycles[i][1] * cycles[j][1]

        H = 1 + 2*alpha_1 + 4*alpha_2

        detail_counter[(c3d, c5d, alpha_2)] += 1

        if bits % 8192 == 0:
            print(f"  {100*bits//(1<<15)}%...", flush=True)

    # Analyze the (c3d, c5d, alpha_2) triples
    print(f"\n  (c3d, c5d, alpha_2) -> H, mod 7:")
    mod7_analysis = Counter()
    for (c3d, c5d, a2), cnt in sorted(detail_counter.items()):
        a1 = c3d + c5d
        H = 1 + 2*a1 + 4*a2
        s = (a1 + 2*a2) % 7
        mod7_analysis[s] += cnt
        if cnt >= 240:  # Only show significant ones
            print(f"    c3d={c3d:2d} c5d={c5d:2d} a2={a2:1d}: a1={a1:2d} H={H:3d} mod7={H%7} s={s} ({cnt})")

    print(f"\n  (a1 + 2*a2) mod 7 summary:")
    for s in range(7):
        pct = 100 * mod7_analysis[s] / 32768
        note = " *** NEEDED FOR H=0 mod 7" if s == 3 else ""
        print(f"    s = {s}: {mod7_analysis[s]:6d} ({pct:.1f}%){note}")

    # 4. c3 mod 7 distribution at n=6
    print(f"\n  c3 (directed) mod 7:")
    c3_mod7 = Counter()
    for (c3d, c5d, a2), cnt in detail_counter.items():
        c3_mod7[c3d % 7] += cnt
    for r in range(7):
        print(f"    c3d mod 7 = {r}: {c3_mod7.get(r, 0)}")

    # 5. Check: at n=6, is alpha_1 mod 7 = 3 achievable?
    print(f"\n  alpha_1 mod 7:")
    a1_mod7 = Counter()
    for (c3d, c5d, a2), cnt in detail_counter.items():
        a1_mod7[(c3d + c5d) % 7] += cnt
    for r in range(7):
        print(f"    a1 mod 7 = {r}: {a1_mod7.get(r, 0)}")

    # 6. Key question: which (a1 mod 7, a2 mod 7) pairs occur?
    print(f"\n  (a1 mod 7, a2 mod 7) that occur:")
    pair_mod7 = Counter()
    for (c3d, c5d, a2), cnt in detail_counter.items():
        a1 = c3d + c5d
        pair_mod7[(a1 % 7, a2 % 7)] += cnt
    for (a1m, a2m) in sorted(pair_mod7.keys()):
        s = (a1m + 2*a2m) % 7
        print(f"    (a1 mod 7, a2 mod 7) = ({a1m},{a2m}): {pair_mod7[(a1m,a2m)]:6d}  s = {s}")

if __name__ == "__main__":
    main()
