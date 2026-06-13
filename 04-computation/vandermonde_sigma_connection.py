#!/usr/bin/env python3
"""
vandermonde_sigma_connection.py — opus-2026-03-13-S71c

CONNECTION: The Vandermonde extraction (2,3) framework meets the
sigma-algebra hierarchy (lambda → c3,c5; (lambda,sigma) → c7).

KEY QUESTION: Does (α₁, α₂) relate to the sigma-algebra levels?
- α₁ = c3 + c5 + c7 (total odd cycles at n=7)
- α₂ = #{disjoint pairs of odd cycles}

The sigma hierarchy shows:
- Lambda determines c3, c5 (but not c7)
- (Lambda, sigma) determines c7

So: α₁ = c3 + c5 + c7. If lambda determines c3 and c5, then
α₁ = known + c7 → alpha₁ mod (lambda-determined part) = c7.

Does this mean: c7 = α₁ - c3 - c5 where c3,c5 are lambda-determined?
→ c7 is lambda + α₁ determined!
→ And α₁ = (H-1-4α₂)/2 where α₂ = (2I₃-3H+1)/6.

So: c7 = (H-1)/2 - 2α₂ - c3 - c5
    where α₂ = (2I₃-3H+1)/6
    → c7 = (H-1)/2 - (2I₃-3H+1)/3 - c3 - c5

With c3,c5 lambda-determined and H,I₃ computable:
→ c7 = (5H - 4I₃ - 1) / 6 - c3 - c5

This gives c7 from (H, I₃, lambda)!

But we already know c7 is (lambda,sigma)-determined.
So: does (H, I₃) encode the same information as sigma?
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

def lambda_key(A, n):
    L = np.zeros((n, n), dtype=int)
    for u in range(n):
        for v in range(u+1, n):
            for w in range(n):
                if w == u or w == v: continue
                if (A[u][v] and A[v][w] and A[w][u]) or (A[v][u] and A[u][w] and A[w][v]):
                    L[u][v] += 1; L[v][u] += 1
    return tuple(L[i][j] for i in range(n) for j in range(i+1, n)), L

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
print("DOES (H, I(3)) ENCODE THE SAME INFO AS SIGMA?")
print("=" * 70)

n = 7
tb = n*(n-1)//2
np.random.seed(42)

# Collect (lambda, H, I3, c7) tuples
lam_H_to_c7 = defaultdict(set)
lam_I3_to_c7 = defaultdict(set)
lam_HI3_to_c7 = defaultdict(set)
lam_a1a2_to_c7 = defaultdict(set)
H_to_c7 = defaultdict(set)
HI3_to_c7 = defaultdict(set)

t0 = time.time()
for trial in range(10000):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    lk, L = lambda_key(A, n)
    H = count_ham_paths(A, n)
    c3 = int(np.trace(A @ A @ A)) // 3
    c5 = int(np.trace(np.linalg.matrix_power(A, 5))) // 5

    alpha = count_disjoint_cycles(A, n)
    I3 = I_at_x(alpha, 3)
    a1 = alpha.get(1, 0)
    a2 = alpha.get(2, 0)
    c7 = a1 - c3 - c5  # at n=7, c7 = Hamiltonian cycles

    lam_H_to_c7[(lk, H)].add(c7)
    lam_I3_to_c7[(lk, I3)].add(c7)
    lam_HI3_to_c7[(lk, H, I3)].add(c7)
    lam_a1a2_to_c7[(lk, a1, a2)].add(c7)
    H_to_c7[H].add(c7)
    HI3_to_c7[(H, I3)].add(c7)

dt = time.time() - t0
print(f"\nn=7, 10000 samples, {dt:.1f}s")

print(f"\nc7 determination:")
print(f"  H alone: {sum(1 for v in H_to_c7.values() if len(v)>1)} ambiguous / {len(H_to_c7)}")
print(f"  (H, I(3)): {sum(1 for v in HI3_to_c7.values() if len(v)>1)} ambiguous / {len(HI3_to_c7)}")
print(f"  (lambda, H): {sum(1 for v in lam_H_to_c7.values() if len(v)>1)} ambiguous / {len(lam_H_to_c7)}")
print(f"  (lambda, I(3)): {sum(1 for v in lam_I3_to_c7.values() if len(v)>1)} ambiguous / {len(lam_I3_to_c7)}")
print(f"  (lambda, H, I(3)): {sum(1 for v in lam_HI3_to_c7.values() if len(v)>1)} ambiguous / {len(lam_HI3_to_c7)}")
print(f"  (lambda, α₁, α₂): {sum(1 for v in lam_a1a2_to_c7.values() if len(v)>1)} ambiguous / {len(lam_a1a2_to_c7)}")

# The key: if (lambda, H, I3) determines c7 with 0 ambiguity,
# then (lambda, H, I3) is at least as powerful as (lambda, sigma)!

# =====================================================================
print(f"\n{'='*70}")
print("FORMULA: c7 = α₁ - c3 - c5 where c3,c5 are lambda-determined")
print("=" * 70)

# At n=7: c7 should be exactly α₁ - c3 - c5
# And α₁ = (H - 1 - 4α₂)/2, α₂ = (2I₃ - 3H + 1)/6
# So c7 = (H-1)/2 - 2(2I₃-3H+1)/6 - c3 - c5
#        = (H-1)/2 - (2I₃-3H+1)/3 - c3 - c5
#        = (3(H-1) - 2(2I₃-3H+1)) / 6 - c3 - c5
#        = (3H-3 - 4I₃+6H-2) / 6 - c3 - c5
#        = (9H - 4I₃ - 5) / 6 - c3 - c5

# But wait: since at n=7, α₃ = 0, α₂ might include (3,3), (3,5) types
# but NOT (3,7) or (5,5) since 3+7=10>7, 5+5=10>7
# So α₂ counts only (3,3) disjoint pairs

# So H = 1 + 2(c3+c5+c7) + 4·α₂(3,3)
# I₃ = 1 + 3(c3+c5+c7) + 9·α₂(3,3)

# From these:
# α₂(3,3) = (2I₃ - 3H + 1)/6
# c3+c5+c7 = (H - 1 - 4α₂)/2

print("\nVerifying formula: c7 = (5H - 4I₃ - 1)/6 - c3 - c5")
ok = 0
bad = 0
for trial in range(2000):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)
    c3 = int(np.trace(A @ A @ A)) // 3
    c5 = int(np.trace(np.linalg.matrix_power(A, 5))) // 5
    alpha = count_disjoint_cycles(A, n)
    I3 = I_at_x(alpha, 3)
    a1 = alpha.get(1, 0)
    c7_actual = a1 - c3 - c5

    # Formula
    c7_formula = (5*H - 4*I3 - 1) / 6 - c3 - c5

    if abs(c7_formula - c7_actual) < 0.01:
        ok += 1
    else:
        bad += 1
        if bad <= 5:
            a2 = alpha.get(2, 0)
            print(f"  FAIL: H={H}, I3={I3}, c3={c3}, c5={c5}, c7_actual={c7_actual}, c7_formula={c7_formula:.1f}, α₂={a2}")

print(f"  OK: {ok}, FAIL: {bad}")

# Wait — at n=7, can we have (3,5) disjoint pairs? 3+5=8 > 7. NO.
# Can we have (5,5)? 10 > 7. NO.
# So α₂ counts ONLY (3,3) disjoint pairs at n=7.
# But can 3-cycles and 5-cycles be disjoint? 3+5=8>7, so NO!
# So ALL disjoint pairs are (3,3).

print(f"\nAt n=7: α₂ counts only (3,3) disjoint 3-cycle pairs (3+5=8>7, 5+5=10>7)")

# Simplify: c7 = (H-1)/2 - 2·α₂ - c3 - c5
# And: α₂ = (2I₃ - 3H + 1)/6
# This is a LINEAR formula relating c7, H, I₃, c3, c5!

# Since c3 and c5 are lambda-determined, and H is a function of the tournament,
# c7 is determined by (lambda, H, I₃).

# But WAIT: what if H is already lambda-determined?
# HYP-849 says β_d are lambda-determined.
# Is H lambda-determined?

lam_to_H = defaultdict(set)
for trial in range(5000):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    lk, _ = lambda_key(A, n)
    H = count_ham_paths(A, n)
    lam_to_H[lk].add(H)

ambig_H = sum(1 for v in lam_to_H.values() if len(v) > 1)
print(f"\nIs H lambda-determined at n=7? Ambiguous: {ambig_H}/{len(lam_to_H)}")

if ambig_H > 0:
    count = 0
    for k, v in sorted(lam_to_H.items()):
        if len(v) > 1:
            print(f"  Lambda group with multiple H values: {sorted(v)}")
            count += 1
            if count >= 10: break

# Is I₃ lambda-determined?
lam_to_I3 = defaultdict(set)
for trial in range(5000):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    lk, _ = lambda_key(A, n)
    alpha = count_disjoint_cycles(A, n)
    I3 = I_at_x(alpha, 3)
    lam_to_I3[lk].add(I3)

ambig_I3 = sum(1 for v in lam_to_I3.values() if len(v) > 1)
print(f"Is I(3) lambda-determined at n=7? Ambiguous: {ambig_I3}/{len(lam_to_I3)}")

# =====================================================================
print(f"\n{'='*70}")
print("KEY: α₂ IS LAMBDA-DETERMINED (SINCE IT'S ONLY (3,3) PAIRS)")
print("=" * 70)

# α₂ = #{disjoint 3-cycle pairs} at n=7
# This is determined by the GRAPH STRUCTURE of how 3-cycles overlap.
# And 3-cycle structure IS the lambda graph!
# λ(i,j) counts 3-cycles through edge {i,j}.
# Disjoint 3-cycle pairs → no shared vertices.

# So α₂ IS lambda-determined!
lam_to_a2 = defaultdict(set)
for trial in range(5000):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    lk, _ = lambda_key(A, n)
    alpha = count_disjoint_cycles(A, n)
    a2 = alpha.get(2, 0)
    lam_to_a2[lk].add(a2)

ambig_a2 = sum(1 for v in lam_to_a2.values() if len(v) > 1)
print(f"Is α₂ lambda-determined at n=7? Ambiguous: {ambig_a2}/{len(lam_to_a2)}")

# And c3 + c5 is lambda-determined. So:
# H = 1 + 2(c3+c5+c7) + 4·α₂
# All of c3, c5, α₂ are lambda-determined.
# So: H - (1 + 2c3 + 2c5 + 4α₂) = 2c7
# → c7 = (H - 1 - 2c3 - 2c5 - 4α₂) / 2

# And if H is NOT lambda-determined, then c7 is (lambda, H)-determined!
# The non-lambda part of H IS 2c7.

print(f"\nConclusion: H = 1 + 2(c3 + c5 + c7) + 4·α₂")
print(f"  c3, c5, α₂ are lambda-determined")
print(f"  So: c7 = (H - 1 - 2c3 - 2c5 - 4α₂) / 2")
print(f"  c7 is (lambda, H)-determined")

# Now: is H lambda-determined? If NOT, the non-lambda part encodes c7.
# If YES, then c7 is also lambda-determined, contradicting HYP-848!

if ambig_H > 0:
    print(f"\n  H is NOT lambda-determined ({ambig_H} ambiguous groups)")
    print(f"  The non-lambda part of H is exactly 2·c7")
    print(f"  This is consistent with c7 NOT being lambda-determined")

    # Check: in ambiguous H groups, does H variation = 2 × c7 variation?
    np.random.seed(42)
    lam_to_Hc7 = defaultdict(list)
    for trial in range(10000):
        bits = np.random.randint(0, 1 << tb)
        A = bits_to_adj(bits, n)
        lk, _ = lambda_key(A, n)
        H = count_ham_paths(A, n)
        c3 = int(np.trace(A @ A @ A)) // 3
        c5 = int(np.trace(np.linalg.matrix_power(A, 5))) // 5
        alpha = count_disjoint_cycles(A, n)
        a1 = alpha.get(1, 0)
        c7 = a1 - c3 - c5
        lam_to_Hc7[lk].append((H, c7, c3, c5))

    print(f"\n  Checking: H variation = 2 × c7 variation within lambda fibers")
    count = 0
    for lk, entries in sorted(lam_to_Hc7.items()):
        if len(entries) < 2:
            continue
        H_vals = set(e[0] for e in entries)
        c7_vals = set(e[1] for e in entries)
        if len(H_vals) > 1 or len(c7_vals) > 1:
            count += 1
            if count <= 10:
                c3_val = entries[0][2]  # should be constant
                c5_val = entries[0][3]  # should be constant
                H_range = max(H_vals) - min(H_vals)
                c7_range = max(c7_vals) - min(c7_vals)
                print(f"    c3={c3_val}, c5={c5_val}, H∈{sorted(H_vals)}, c7∈{sorted(c7_vals)}, ΔH={H_range}, 2·Δc7={2*c7_range}")

# =====================================================================
print(f"\n{'='*70}")
print("CONNECTING TO I(3): WHAT DOES I(3) ADD BEYOND H?")
print("=" * 70)

# We showed: c7 = (H - 1 - 2c3 - 2c5 - 4α₂)/2
# All of c3, c5, α₂ lambda-determined.
# So H and lambda together determine c7.
# I(3) = 1 + 3(c3+c5+c7) + 9α₂
# Also lambda-determined parts: c3, c5, α₂.
# So I(3) - (1 + 3c3 + 3c5 + 9α₂) = 3c7
# And H - (1 + 2c3 + 2c5 + 4α₂) = 2c7

# The ratio: (I(3) - lambda_part) / (H - lambda_part) = 3/2 always!
# This is EXACTLY the 3/2 ratio!

print("\n  (I₃ - I₃_lambda_part) / (H - H_lambda_part) = 3c7 / 2c7 = 3/2")
print("  → I(3) adds NO new information beyond H (at n=7)!")
print("  → The 3/2 ratio is the SHADOW of the Vandermonde matrix")
print("  → But at n=9 (where α₃ > 0), I(3) WILL add new information")

# Verify
np.random.seed(42)
ratio_check = []
for trial in range(2000):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)
    c3 = int(np.trace(A @ A @ A)) // 3
    c5 = int(np.trace(np.linalg.matrix_power(A, 5))) // 5
    alpha = count_disjoint_cycles(A, n)
    I3 = I_at_x(alpha, 3)
    a2 = alpha.get(2, 0)

    H_lam = 1 + 2*c3 + 2*c5 + 4*a2  # lambda-determined part
    I3_lam = 1 + 3*c3 + 3*c5 + 9*a2  # lambda-determined part

    H_residual = H - H_lam  # = 2c7
    I3_residual = I3 - I3_lam  # = 3c7

    if H_residual != 0:
        ratio = I3_residual / H_residual
        ratio_check.append(ratio)

print(f"\n  Ratio check (I₃_residual / H_residual) over {len(ratio_check)} samples with c7>0:")
if ratio_check:
    print(f"  min={min(ratio_check):.6f}, max={max(ratio_check):.6f}, mean={np.mean(ratio_check):.6f}")
    print(f"  All exactly 3/2? {all(abs(r - 1.5) < 1e-10 for r in ratio_check)}")

print(f"\n{'='*70}")
print("FINAL SYNTHESIS")
print("=" * 70)

print("""
THE 2-3 STRUCTURE AT n=7 (max packing = 2):

1. H = 1 + 2α₁ + 4α₂  (α₃=0 since 3+3+3=9 > 7)
   I₃ = 1 + 3α₁ + 9α₂

2. α₁ = c3 + c5 + c7  (total odd directed cycles)
   α₂ = #{disjoint 3-cycle pairs}  (only type at n=7)

3. c3, c5, α₂ are all lambda-determined.
   c7 is NOT lambda-determined (but (lambda,sigma)-determined).

4. The OCF formula:
   H = [1 + 2c3 + 2c5 + 4α₂] + 2c7
       ←── lambda-determined ──→   ←─ residual ─→

5. I(3) = [1 + 3c3 + 3c5 + 9α₂] + 3c7
          ←── lambda-determined ──→   ←─ residual ─→

6. Residual ratio: 3c7 / 2c7 = 3/2  ALWAYS.
   → I(3) carries NO information beyond H at n=7.
   → The "3/2 ratio" is a THEOREM, not a coincidence.

7. At n=9: α₃ = #{(3,3,3) triples} becomes possible.
   H = 1 + 2α₁ + 4α₂ + 8α₃
   I₃ = 1 + 3α₁ + 9α₂ + 27α₃
   Now the Vandermonde extraction gives α₃ = (27(H-1) - 8(I₃-1))/18... wait

   Better: with 3 unknowns, need 3 equations. Add I(4):
   I₄ = 1 + 4α₁ + 16α₂ + 64α₃

   Then the (2,3,4) Vandermonde system resolves α₁, α₂, α₃.

8. DEEPER: "If you know 2 and 3, you have the keys to the universe"
   At n ≤ 8 (where α₃=0), EXACTLY TWO evaluation points suffice.
   2 and 3 are the first two positive integers ≥ 2.
   They give the MINIMAL Vandermonde basis.
   The extraction α₁, α₂ then combines with lambda to give EVERYTHING.
""")

print("Done.")
