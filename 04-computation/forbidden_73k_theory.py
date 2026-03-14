"""
forbidden_73k_theory.py -- kind-pasteur-2026-03-14-S67

Theoretical analysis: WHY is {7, 21, 63} = {7*3^k : k=0,1,2} forbidden?

Key structural constraint: H = I(Omega(T), 2) where I is independence polynomial.
H is ALWAYS odd (THM-002), and H = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ...

For H to equal 7*3^k, we need I(Omega(T), 2) = 7*3^k.

Observation: 7*3^k - 1 = 2*(7*3^k - 1)/2, and we need this to decompose as
2*alpha_1 + 4*alpha_2 + ... = 7*3^k - 1.

Let's check the modular constraints.

Also: I(G, x) at x=2 has a 3-adic structure related to the graph G.
Can we show v_3(I(Omega(T), 2)) != k for specific k values?

Key fact: for K_3 (triangle), I(K_3, 2) = 1 + 3*2 = 7.
So 7 = I(K_3, 2). This is the independence polynomial of the conflict graph
of a 3-vertex tournament (which always has Omega = K_3... wait, does it?

Actually for a 3-vertex tournament, there is exactly ONE directed 3-cycle
(in ONE direction). So Omega has exactly 1 vertex and no edges.
I(Omega, 2) = 1 + 2 = 3. But H(3-vertex tournament) = 3, not 7. Wait...

Let me reconsider. For n=3: single tournament T_3, H = 3.
alpha_1 = number of directed Hamiltonian cycles = 1 (the unique 3-cycle has 2 directed versions? No, a tournament on 3 vertices has exactly 1 directed Hamiltonian cycle).

Hmm. Let me recalculate properly.
"""

from itertools import combinations
from collections import Counter
import numpy as np

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

def compute_alpha(A, n):
    """Compute all alpha values."""
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

    alpha_3 = 0
    for i in range(len(cycles)):
        for j in range(i+1, len(cycles)):
            if len(cycles[i][0] & cycles[j][0]) > 0:
                continue
            for k_idx in range(j+1, len(cycles)):
                if len(cycles[i][0] & cycles[k_idx][0]) == 0 and \
                   len(cycles[j][0] & cycles[k_idx][0]) == 0:
                    alpha_3 += cycles[i][1] * cycles[j][1] * cycles[k_idx][1]

    return alpha_1, alpha_2, alpha_3, cycles

def main():
    print("=" * 70)
    print("FORBIDDEN 7*3^k — THEORETICAL ANALYSIS")
    print("=" * 70)

    # 1. What is alpha_1 at small n?
    print("\n--- ALPHA STRUCTURE AT n=3 ---")
    # Only tournament on 3 vertices (up to iso): 0->1->2->0
    n = 3
    A = np.zeros((n, n), dtype=np.int8)
    A[0][1] = 1; A[1][2] = 1; A[2][0] = 1
    a1, a2, a3, cycles = compute_alpha(A, n)
    H = 1 + 2*a1 + 4*a2 + 8*a3
    print(f"  Cyclic tournament: alpha_1={a1}, alpha_2={a2}, alpha_3={a3}, H={H}")
    print(f"  Cycles: {[(list(s), c, sz) for s,c,sz in cycles]}")

    # Transitive tournament on 3 vertices: 0->1, 0->2, 1->2
    A2 = np.zeros((n, n), dtype=np.int8)
    A2[0][1] = 1; A2[0][2] = 1; A2[1][2] = 1
    a1, a2, a3, cycles = compute_alpha(A2, n)
    H = 1 + 2*a1 + 4*a2 + 8*a3
    print(f"  Transitive tournament: alpha_1={a1}, alpha_2={a2}, alpha_3={a3}, H={H}")

    # 2. 3-adic valuations of achievable H at n=6
    print("\n--- 3-ADIC STRUCTURE OF H VALUES AT n=6 ---")
    n = 6
    H_counter = Counter()
    alpha_data = {}

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

        a1, a2, a3, cycles = compute_alpha(A, n)
        H = 1 + 2*a1 + 4*a2 + 8*a3
        H_counter[H] += 1
        if H not in alpha_data:
            alpha_data[H] = []
        if len(alpha_data[H]) < 3:  # Store a few examples
            alpha_data[H].append((a1, a2, a3))

        if bits % 8192 == 0:
            print(f"  {100*bits//(1<<15)}%...", flush=True)

    print(f"\n  H values and 3-adic valuations:")
    for H in sorted(H_counter.keys()):
        v3 = 0
        h = H
        while h % 3 == 0:
            v3 += 1
            h //= 3
        residue = H // (3**v3)
        ex = alpha_data[H][0]
        print(f"  H={H:4d} = {residue}*3^{v3}  (alpha=({ex[0]},{ex[1]},{ex[2]}))  count={H_counter[H]}")

    # 3. Check what alpha decompositions are near 7, 21, 63
    print(f"\n--- ALPHA DECOMPOSITIONS NEAR FORBIDDEN VALUES ---")
    for target in [7, 21, 63]:
        print(f"\n  Target H={target} = 7*3^{int(np.log(target/7)/np.log(3)+0.5)}:")
        # What (alpha_1, alpha_2) pairs give H=target?
        # H = 1 + 2*a1 + 4*a2
        remainder = target - 1
        print(f"    Need 2*a1 + 4*a2 [+ 8*a3 ...] = {remainder}")
        for a2 in range(remainder // 4 + 1):
            for a3 in range(remainder // 8 + 1):
                leftover = remainder - 4*a2 - 8*a3
                if leftover >= 0 and leftover % 2 == 0:
                    a1 = leftover // 2
                    print(f"    (a1={a1}, a2={a2}, a3={a3})")

    # 4. What alpha_1 values are achievable at n=6?
    print(f"\n--- ACHIEVABLE alpha_1 VALUES AT n=6 ---")
    a1_vals = set()
    a1_counter = Counter()
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

        a1, a2, a3, _ = compute_alpha(A, n)
        a1_vals.add(a1)
        a1_counter[a1] += 1

    print(f"  Achievable alpha_1: {sorted(a1_vals)}")

    # For each achievable alpha_1, what alpha_2 values are possible?
    print(f"\n--- (alpha_1, alpha_2) PAIRS AT n=6 ---")
    pair_counter = Counter()
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

        a1, a2, a3, _ = compute_alpha(A, n)
        pair_counter[(a1, a2)] += 1

    # Show all pairs and which H they give
    print(f"  {'a1':>4s} {'a2':>4s} {'a3':>4s}  {'H':>5s}  count")
    for (a1, a2), cnt in sorted(pair_counter.items()):
        H = 1 + 2*a1 + 4*a2
        print(f"  {a1:4d} {a2:4d}    0  {H:5d}  {cnt:5d}")

    # 5. Which (a1, a2) pairs would give forbidden H?
    print(f"\n--- MISSING PAIRS THAT WOULD GIVE FORBIDDEN H ---")
    for target in [7, 21, 35, 39]:
        remainder = target - 1
        for a2 in range(remainder // 4 + 1):
            leftover = remainder - 4*a2
            if leftover >= 0 and leftover % 2 == 0:
                a1 = leftover // 2
                exists = (a1, a2) in pair_counter
                print(f"  H={target}: (a1={a1}, a2={a2}) -> {'EXISTS' if exists else 'MISSING'}")

    # 6. H mod 7 structure
    print(f"\n--- H mod 7 AT n=6 ---")
    mod7 = Counter()
    for H, cnt in H_counter.items():
        mod7[H % 7] += cnt
    for r in range(7):
        pct = 100 * mod7[r] / 32768
        print(f"  H mod 7 = {r}: {mod7[r]:6d} ({pct:.1f}%)")

    # 7. Which H = 0 mod 7 are achievable?
    print(f"\n  Achievable H divisible by 7:")
    for H in sorted(H_counter.keys()):
        if H % 7 == 0:
            print(f"    H={H} = 7*{H//7}")

    print(f"\n  Missing H divisible by 7 (up to max H={max(H_counter.keys())}):")
    max_H = max(H_counter.keys())
    for h in range(7, max_H + 1, 14):  # 7 mod 14 gives odd multiples of 7
        if h not in H_counter:
            print(f"    H={h} = 7*{h//7}")

if __name__ == "__main__":
    main()
