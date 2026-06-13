#!/usr/bin/env python3
"""
i3_mod3_proof.py — opus-2026-03-13-S71c

DISCOVERY: I(CG, 3) ≡ 1 (mod 3) for ALL tournaments.

Since I(CG, 3) = 1 + 3α₁ + 9α₂ + 27α₃ + ...,
this is TRIVIALLY true: I(3) = 1 + 3(α₁ + 3α₂ + 9α₃ + ...) ≡ 1 (mod 3).

But wait — what about mod 9? mod 27?
I(3) mod 9 = 1 + 3α₁ mod 9 = 1 + 3α₁ (mod 9), depends on α₁ mod 3.
α₁ = total number of odd directed cycles = c3 + c5 + c7 + ...

So: I(3) mod 9 is determined by α₁ mod 3.

More interesting questions:
1. Is α₁ mod 3 a tournament invariant with nice properties?
2. What does I(3) mod higher powers of 3 encode?
3. Compare: H ≡ 1 (mod 2) always (OCF). Is there a "mod 3 OCF"?
4. Does I(CG, x) ≡ 1 (mod x) for all x? YES trivially!

Actually the deep question: what about I(CG, -1)?
I(CG, -1) = 1 - α₁ + α₂ - α₃ + ... (alternating sum)
This is the EULER CHARACTERISTIC of the independence complex!

And I(CG, 0) = 1 (trivial).
I(CG, 1) = 1 + α₁ + α₂ + ... = total independent sets in CG.

THE REAL INSIGHT: Vandermonde at x=2,3 gives α₁, α₂.
Can we use x=-1 and x=2 instead?
I(-1) = 1 - α₁ + α₂ - α₃ + ...
I(2) = 1 + 2α₁ + 4α₂ + 8α₃ + ...

From these: 4·I(-1) = 4 - 4α₁ + 4α₂ - 4α₃ + ...
           I(2) = 1 + 2α₁ + 4α₂ + 8α₃ + ...
4·I(-1) + I(2) = 5 - 2α₁ + 8α₂ + 4α₃ + ... hmm messy

Better: I(2) + 2·I(-1) = 3 + 0·α₁ + 8α₂ + 6α₃ + ...
        = 3 + 2(4α₂ + 3α₃ + ...)
        ⟹ α₁ is NOT directly extractable from I(2) and I(-1) alone!

But I(2) - 2·I(-1) = -1 + 4α₁ + 0·α₂ + 12α₃ + ...
So (I(2) - 2·I(-1) + 1)/4 = α₁ + 3α₃ + ... ≡ α₁ (mod 3)

Interesting but messy. Let me compute I(-1) for tournaments.
"""

import sys, time
import numpy as np
from itertools import combinations, permutations
from collections import defaultdict
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
            if bin(mask).count('1') != ms:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                pm = mask ^ (1 << v)
                t = 0
                for u in range(n):
                    if (pm & (1 << u)) and A[u][v]:
                        t += dp.get((pm, u), 0)
                if t:
                    dp[(mask, v)] = t
    return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))

def count_disjoint_cycles(A, n):
    cycles = []
    seen = set()
    def find_cycles(path, start, used):
        last = path[-1]
        length = len(path)
        if length >= 3 and length % 2 == 1 and A[last][start]:
            normalized = min(path[i:] + path[:i] for i in range(length))
            key = tuple(normalized)
            if key not in seen:
                seen.add(key)
                cycles.append(frozenset(path))
        if length >= n:
            return
        for v in range(n):
            if v not in used and A[last][v]:
                find_cycles(path + [v], start, used | {v})
    for start in range(n):
        find_cycles([start], start, {start})

    alpha = defaultdict(int)
    def count_collections(idx, used, k):
        if k > 0:
            alpha[k] += 1
        for j in range(idx, len(cycles)):
            if not (cycles[j] & used):
                count_collections(j+1, used | cycles[j], k+1)
    count_collections(0, frozenset(), 0)
    return dict(alpha)

def I_at_x(alpha, x):
    return 1 + sum(count * x**k for k, count in alpha.items())

# =====================================================================
print("=" * 70)
print("I(CG, -1) — THE EULER CHARACTERISTIC")
print("=" * 70)

for n in [5, 6, 7]:
    print(f"\nn={n}:")
    tb = n*(n-1)//2
    np.random.seed(42)

    i_neg1_vals = defaultdict(int)
    i_neg1_mod = defaultdict(int)

    count = min(200, 1 << tb) if n <= 6 else 500

    for trial in range(count):
        if n <= 5:
            bits = trial if trial < (1 << tb) else np.random.randint(0, 1 << tb)
        else:
            bits = np.random.randint(0, 1 << tb)

        A = bits_to_adj(bits, n)
        alpha = count_disjoint_cycles(A, n)
        I_neg1 = I_at_x(alpha, -1)
        H = I_at_x(alpha, 2)
        I3 = I_at_x(alpha, 3)

        i_neg1_vals[I_neg1] += 1
        i_neg1_mod[I_neg1 % 2] += 1

    print(f"  I(-1) distribution: {dict(sorted(i_neg1_vals.items()))}")
    print(f"  I(-1) mod 2: {dict(sorted(i_neg1_mod.items()))}")

# =====================================================================
print(f"\n{'='*70}")
print("DEEPER: I(3) mod 9, mod 27 STRUCTURE")
print("=" * 70)

n = 7
tb = n*(n-1)//2
np.random.seed(42)

i3_mod9 = defaultdict(int)
i3_mod27 = defaultdict(int)
a1_mod3 = defaultdict(int)
a1_i3mod9 = defaultdict(set)

for trial in range(5000):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    alpha = count_disjoint_cycles(A, n)
    I3 = I_at_x(alpha, 3)
    a1 = alpha.get(1, 0)

    i3_mod9[I3 % 9] += 1
    i3_mod27[I3 % 27] += 1
    a1_mod3[a1 % 3] += 1
    a1_i3mod9[a1 % 3].add(I3 % 9)

print(f"\nn=7, 5000 samples:")
print(f"  I(3) mod 9: {dict(sorted(i3_mod9.items()))}")
print(f"  I(3) mod 27: {dict(sorted(i3_mod27.items()))}")
print(f"  α₁ mod 3: {dict(sorted(a1_mod3.items()))}")
print(f"\n  Correspondence α₁ mod 3 → I(3) mod 9:")
for k, v in sorted(a1_i3mod9.items()):
    print(f"    α₁ ≡ {k} (mod 3) → I(3) mod 9 ∈ {sorted(v)}")
    # Expected: I(3) mod 9 = (1 + 3·α₁) mod 9 if α₂ ≡ 0 mod 9
    # But 9α₂ ≡ 0 mod 9 always. So I(3) mod 9 = (1 + 3α₁) mod 9.
    for a1val in sorted(v):
        expected = (1 + 3*k) % 9
        print(f"      actual={a1val}, expected from (1+3·{k}) mod 9 = {expected}")

# =====================================================================
print(f"\n{'='*70}")
print("THE VANDERMONDE EXTRACTION: α₁ FROM (H, I(3))")
print("=" * 70)

print("\nα₂ = (2·I(3) - 3·H + 1) / 6")
print("α₁ = (H - 1 - 4·α₂) / 2")
print("\nAt n=7, this extracts α₁ = c3 + c5 + c7 (total odd cycles)")
print("α₂ = total disjoint pairs of odd cycles")

# Check what α₁ and α₂ encode
n = 7
tb = n*(n-1)//2
np.random.seed(42)

a1_to_c3c5c7 = defaultdict(set)
a2_to_config = defaultdict(set)

for trial in range(2000):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    A_np = np.array(A) if not isinstance(A, np.ndarray) else A
    c3 = int(np.trace(A_np @ A_np @ A_np)) // 3
    c5 = int(np.trace(np.linalg.matrix_power(A_np, 5))) // 5

    alpha = count_disjoint_cycles(A_np, n)
    a1 = alpha.get(1, 0)
    a2 = alpha.get(2, 0)

    # At n=7: c7 = number of Hamiltonian cycles
    # α₁ = c3 + c5 + c7, α₂ = disjoint pairs
    c7 = a1 - c3 - c5  # since only odd cycles of length 3,5,7 exist at n=7
    a1_to_c3c5c7[a1].add((c3, c5, c7))

    H = I_at_x(alpha, 2)
    I3 = I_at_x(alpha, 3)
    a2_recovered = (2*I3 - 3*H + 1) / 6
    a1_recovered = (H - 1 - 4*a2_recovered) / 2

    if abs(a1_recovered - a1) > 0.01 or abs(a2_recovered - a2) > 0.01:
        print(f"  MISMATCH! trial={trial}")

print(f"\nα₁ determines (c3,c5,c7) individually? {sum(1 for v in a1_to_c3c5c7.values() if len(v) > 1)} ambiguous / {len(a1_to_c3c5c7)} total")

# Does (α₁, α₂) = (H, I3) determine c3 individually?
h_i3_to_c3 = defaultdict(set)
h_i3_to_c5 = defaultdict(set)
h_i3_to_c7 = defaultdict(set)

for trial in range(2000):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    A_np = np.array(A) if not isinstance(A, np.ndarray) else A
    c3 = int(np.trace(A_np @ A_np @ A_np)) // 3
    c5 = int(np.trace(np.linalg.matrix_power(A_np, 5))) // 5
    alpha = count_disjoint_cycles(A_np, n)
    a1 = alpha.get(1, 0)
    a2 = alpha.get(2, 0)
    c7 = a1 - c3 - c5
    key = (a1, a2)
    h_i3_to_c3[key].add(c3)
    h_i3_to_c5[key].add(c5)
    h_i3_to_c7[key].add(c7)

print(f"\n(α₁,α₂) determines c3? Ambiguous: {sum(1 for v in h_i3_to_c3.values() if len(v) > 1)} / {len(h_i3_to_c3)}")
print(f"(α₁,α₂) determines c5? Ambiguous: {sum(1 for v in h_i3_to_c5.values() if len(v) > 1)} / {len(h_i3_to_c5)}")
print(f"(α₁,α₂) determines c7? Ambiguous: {sum(1 for v in h_i3_to_c7.values() if len(v) > 1)} / {len(h_i3_to_c7)}")

# =====================================================================
print(f"\n{'='*70}")
print("I(CG, x) AT SPECIAL POINTS")
print("=" * 70)

n = 7
tb = n*(n-1)//2
np.random.seed(42)

print(f"\nn=7, 500 samples:")
print(f"{'H':>5} {'I(-1)':>6} {'I(1)':>6} {'I(2)':>6} {'I(3)':>6} {'I(4)':>6} {'α₁':>4} {'α₂':>4} {'α₃':>4}")

# Count how many have α₃ > 0
a3_count = 0
for trial in range(500):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    alpha = count_disjoint_cycles(A, n)
    H = I_at_x(alpha, 2)
    vals = {x: I_at_x(alpha, x) for x in [-1, 1, 2, 3, 4]}
    a1 = alpha.get(1, 0)
    a2 = alpha.get(2, 0)
    a3 = alpha.get(3, 0)
    if a3 > 0:
        a3_count += 1
    if trial < 15:
        print(f"{H:5d} {vals[-1]:6d} {vals[1]:6d} {vals[2]:6d} {vals[3]:6d} {vals[4]:6d} {a1:4d} {a2:4d} {a3:4d}")

print(f"\n  α₃ > 0 count: {a3_count}/500")
print(f"  (α₃ nonzero only possible if 3 disjoint odd cycles fit in n=7,")
print(f"   requiring 3+3+3=9 > 7 vertices. So α₃ = 0 ALWAYS at n=7!)")

# =====================================================================
print(f"\n{'='*70}")
print("SYNTHESIS: THE 2-3 PRINCIPLE")
print("=" * 70)

print("""
THEOREM (Vandermonde Extraction):
  For any tournament on n vertices with max packing number m = ⌊n/3⌋:

  The values I(CG, 2) and I(CG, 3) uniquely determine α₁ and α₂.

  More generally, evaluations at x = 2, 3, ..., m+1 determine
  ALL α_k (k = 1, ..., m) via Vandermonde inversion.

  For n ≤ 8: m ≤ 2, so (2, 3) suffices.
  For n = 9: first possibility of α₃ > 0 (three disjoint 3-cycles).
    Need x = 2, 3, 4 to fully resolve.

CONGRUENCES:
  I(x) ≡ 1 (mod x) for all x (trivially, since I(x) = 1 + x·(...))
  I(2) ≡ 1 (mod 2)  ← this is the OCF parity theorem (H is odd)
  I(3) ≡ 1 (mod 3)  ← trivial analog
  I(3) mod 9 = 1 + 3·α₁ (mod 9) ← determines α₁ mod 3

DEEP CONNECTION:
  k-nacci constant → 2: the OCF evaluation point
  weighted k-nacci constant → 3: the "second evaluation point"

  Together they form the MINIMAL VANDERMONDE BASIS for
  extracting the first two levels of the odd-cycle hierarchy.

  The Jacobsthal number J(n) = (2^n - (-1)^n)/3 bridges 2 and 3:
  - Numerator: powers of 2 (OCF world)
  - Denominator: 3 (normalization from weighted world)

  "If you know 2 and 3, you know the universe" because:
  2 and 3 are the evaluation points that decompose H into its
  constituent cycle-packing levels. From H alone you get α₁ + 2α₂.
  Adding I(3) gives you α₁ and α₂ separately.
""")

print("Done.")
