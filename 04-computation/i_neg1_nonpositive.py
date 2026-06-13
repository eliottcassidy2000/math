#!/usr/bin/env python3
"""
i_neg1_nonpositive.py — opus-2026-03-14-S71d

POTENTIAL THEOREM: I(Ω(T), -1) ≤ 0 for all tournaments T.

I(Ω, -1) = 1 - α₁ + α₂ - α₃ + ...
         = Σ_{k≥0} (-1)^k α_k

This is the "signed independence number" — alternating count
of independent sets by size.

If I(-1) ≤ 0 always, then:
  α₁ - α₂ + α₃ - ... ≥ 1

At n≤8 (α₃=0): this becomes α₁ ≥ α₂ + 1.
Since α₁ = #(directed odd cycles) and α₂ = #(disjoint pairs):
  #(cycles) ≥ #(disjoint pairs) + 1

This is NOT trivial! The disjoint pairs are counted among C(α₁,2),
but not all pairs are disjoint. The bound says there's always
"at least one more cycle than there are disjoint pairs."
"""

import sys
import numpy as np
from itertools import combinations
from collections import Counter
sys.stdout.reconfigure(line_buffering=True)

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for ms in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != ms: continue
            for v in range(n):
                if not (mask & (1 << v)): continue
                pm = mask ^ (1 << v)
                t = 0
                for u in range(n):
                    if (pm & (1 << u)) and A[u][v]:
                        t += dp.get((pm, u), 0)
                if t:
                    dp[(mask, v)] = t
    return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))

def count_ham_cycles(A, n):
    if n < 3: return 0
    full_mask = (1 << n) - 1
    dp = {(1 << 0, 0): 1}
    for ms in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != ms: continue
            if not (mask & 1): continue
            for v in range(n):
                if not (mask & (1 << v)): continue
                if v == 0 and ms < n: continue
                pm = mask ^ (1 << v)
                if not (pm & 1): continue
                t = 0
                for u in range(n):
                    if (pm & (1 << u)) and A[u][v]:
                        t += dp.get((pm, u), 0)
                if t:
                    dp[(mask, v)] = t
    total = 0
    for v in range(1, n):
        if A[v][0] and (full_mask, v) in dp:
            total += dp[(full_mask, v)]
    return total

def count_directed_k_cycles(A, n, k):
    if k > n: return 0
    total = 0
    for combo in combinations(range(n), k):
        verts = list(combo)
        sub = np.zeros((k, k), dtype=int)
        for i in range(k):
            for j in range(k):
                sub[i][j] = A[verts[i]][verts[j]]
        total += count_ham_cycles(sub, k)
    return total

# ======================================================================
# EXHAUSTIVE CHECK AT n=3,4,5,6
# ======================================================================
print("=" * 70)
print("I(Ω, -1) ≤ 0: EXHAUSTIVE CHECK")
print("=" * 70)

for n in [3, 4, 5, 6]:
    tb = n*(n-1)//2
    violations = 0
    zeros = 0
    min_val = float('inf')
    max_val = float('-inf')

    for bits in range(1 << tb):
        A = bits_to_adj(bits, n)
        H = count_ham_paths(A, n)
        # Compute α₁
        a1 = 0
        for k in range(3, n+1, 2):
            a1 += count_directed_k_cycles(A, n, k)
        a2 = (H - 1 - 2*a1) // 4  # valid since α₃=0 for n≤8

        I_neg1 = 1 - a1 + a2
        min_val = min(min_val, I_neg1)
        max_val = max(max_val, I_neg1)
        if I_neg1 > 0:
            violations += 1
        if I_neg1 == 0:
            zeros += 1

    total = 1 << tb
    print(f"\n  n={n}: {total} tournaments")
    print(f"    I(-1) range: [{min_val}, {max_val}]")
    print(f"    I(-1) > 0: {violations} (violations)")
    print(f"    I(-1) = 0: {zeros}")
    print(f"    I(-1) < 0: {total - violations - zeros}")

# ======================================================================
# SAMPLE CHECK AT n=7,8
# ======================================================================
print("\n" + "=" * 70)
print("I(Ω, -1) ≤ 0: SAMPLE CHECK AT n=7,8")
print("=" * 70)

for n in [7, 8]:
    tb = n*(n-1)//2
    np.random.seed(42)
    N = 2000 if n == 7 else 500

    violations = 0
    zeros = 0
    min_val = float('inf')
    max_val = float('-inf')

    for trial in range(N):
        bits = np.random.randint(0, 1 << tb)
        A = bits_to_adj(bits, n)
        H = count_ham_paths(A, n)
        a1 = 0
        for k in range(3, n+1, 2):
            a1 += count_directed_k_cycles(A, n, k)
        a2 = (H - 1 - 2*a1) // 4

        I_neg1 = 1 - a1 + a2
        min_val = min(min_val, I_neg1)
        max_val = max(max_val, I_neg1)
        if I_neg1 > 0:
            violations += 1
        if I_neg1 == 0:
            zeros += 1

    print(f"\n  n={n}: {N} samples")
    print(f"    I(-1) range: [{min_val}, {max_val}]")
    print(f"    I(-1) > 0: {violations} (violations)")
    print(f"    I(-1) = 0: {zeros}")
    print(f"    I(-1) < 0: {N - violations - zeros}")

# ======================================================================
# WHEN IS I(-1) = 0?
# ======================================================================
print("\n" + "=" * 70)
print("WHEN IS I(-1) = 0?")
print("=" * 70)

n = 5
tb = n*(n-1)//2
for bits in range(1 << tb):
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)
    dc3 = count_directed_k_cycles(A, n, 3)
    dc5 = count_ham_cycles(A, n)
    a1 = dc3 + dc5
    a2 = (H - 1 - 2*a1) // 4
    I_neg1 = 1 - a1 + a2
    if I_neg1 == 0:
        scores = tuple(sorted([sum(A[i]) for i in range(n)]))
        print(f"  n=5, bits={bits}: I(-1)=0, dc3={dc3}, dc5={dc5}, α₁={a1}, α₂={a2}, H={H}, score={scores}")

n = 6
tb = n*(n-1)//2
count_zero = 0
for bits in range(1 << tb):
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)
    a1 = 0
    for k in range(3, n+1, 2):
        a1 += count_directed_k_cycles(A, n, k)
    a2 = (H - 1 - 2*a1) // 4
    I_neg1 = 1 - a1 + a2
    if I_neg1 == 0:
        count_zero += 1
        if count_zero <= 5:
            scores = tuple(sorted([sum(A[i]) for i in range(n)]))
            print(f"  n=6, bits={bits}: I(-1)=0, α₁={a1}, α₂={a2}, H={H}, score={scores}")

print(f"  n=6: {count_zero} total with I(-1)=0")

# ======================================================================
# PROOF SKETCH
# ======================================================================
print("\n" + "=" * 70)
print("PROOF SKETCH: I(Ω, -1) ≤ 0")
print("=" * 70)

print("""
  CLAIM: For any tournament T on n vertices, I(Ω(T), -1) ≤ 0.

  Equivalent to: α₁ ≥ α₂ + 1 when α₃ = 0.
  Or: #(directed odd cycles) > #(vertex-disjoint cycle pairs).

  PROOF IDEA (for n ≤ 8 where α₃ = 0):

  α₁ = total number of directed odd cycles.
  α₂ = number of vertex-disjoint pairs of directed odd cycles.

  Each disjoint pair {C, C'} uses two distinct cycles C and C'.
  The number of pairs is at most C(α₁, 2) = α₁(α₁-1)/2.
  But we need the STRONGER bound: α₂ ≤ α₁ - 1.

  This would follow if we could show that the "disjointness graph"
  (where vertices = directed cycles, edges = disjoint pairs)
  has at most |V|-1 = α₁-1 edges.

  But this is NOT always true! At n=7 with α₁=80 and α₂=14,
  we have α₂ << α₁ - 1 = 79.

  Actually, the bound α₂ ≤ α₁ - 1 IS true (from exhaustive data):
  - n=5: α₂=0, α₁≥0. Need α₁≥1 when α₂=0: TRUE except transitive
    (transitive: α₁=0, I(-1)=1 > 0!) Wait — let me recheck.

  WAIT: The transitive tournament at n=3 has dc3=0, α₁=0, α₂=0.
  I(-1) = 1 - 0 + 0 = 1 > 0!

  So I(-1) CAN be positive (=1) when the tournament is transitive!
  Let me re-verify...
""")

# Re-check n=3
print("  RECHECK n=3:")
n = 3
tb = 3
for bits in range(1 << tb):
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)
    dc3 = count_directed_k_cycles(A, n, 3)
    a1 = dc3
    a2 = 0
    I_neg1 = 1 - a1 + a2
    scores = tuple(sorted([sum(A[i]) for i in range(n)]))
    print(f"    bits={bits}: dc3={dc3}, α₁={a1}, I(-1)={I_neg1}, H={H}, score={scores}")

print("\n  Hmm, at n=3: transitive (dc3=0) has I(-1)=1 > 0!")
print("  At n=3: cyclic (dc3=1) has I(-1)=0.")
print("  So I(-1) > 0 IS possible — for transitive tournaments!")
print("  The claim I(-1) ≤ 0 is FALSE at n=3 and n=4.")

print("\n  Let's check: at what n does I(-1) become ≤ 0 for ALL tournaments?")

print("\nDone.")
