#!/usr/bin/env python3
"""
convergence_rate_23.py — opus-2026-03-14-S71d

THE DEEP 2-3 CONVERGENCE STRUCTURE

Discoveries so far:
  φ_k (k-nacci root) → 2 at rate (1/2)^k
  ψ_k (weighted k-nacci root) → 3 at rate (2/3)^k

Questions:
1. Is rate = |subdominant/dominant| exactly (r-1)/r?
2. What is (φ_k - 2)/(ψ_k - 3) as k→∞?
3. Does I(CG, φ_k) have special meaning for n=k tournaments?
4. Fix 7-cycle enumeration and verify I(CG,2) = H properly
5. The ratio gap₃/gap₂ — does it converge to 3/2?
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

def count_directed_cycles(A, n, length):
    """Count directed cycles of given length in tournament A.
    Returns the number of VERTEX SETS that support at least one directed cycle."""
    count = 0
    for combo in combinations(range(n), length):
        verts = list(combo)
        # Check all cyclic orderings
        found = False
        # Fix first vertex, permute rest
        for perm in permutations(verts[1:]):
            order = [verts[0]] + list(perm)
            if all(A[order[i]][order[(i+1) % length]] for i in range(length)):
                found = True
                break
        if found:
            count += 1
    return count

def find_all_directed_cycles(A, n, max_length=None):
    """Find ALL directed cycles as frozensets of vertex indices."""
    if max_length is None:
        max_length = n
    cycles = []
    for length in range(3, max_length + 1, 2):  # odd lengths only
        if length > n:
            break
        for combo in combinations(range(n), length):
            verts = list(combo)
            found = False
            for perm in permutations(verts[1:]):
                order = [verts[0]] + list(perm)
                if all(A[order[i]][order[(i+1) % length]] for i in range(length)):
                    found = True
                    break
            if found:
                cycles.append(frozenset(combo))
    return cycles

def independence_polynomial_exact(A, n, x_vals):
    """Compute I(CG, x) exactly using all odd cycles."""
    cycles = find_all_directed_cycles(A, n)
    nc = len(cycles)

    # Build conflict graph
    adj = [[False]*nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if cycles[i] & cycles[j]:
                adj[i][j] = True
                adj[j][i] = True

    # Count independent sets by size
    max_size = nc
    alpha = [0] * (max_size + 1)

    for mask in range(1 << nc):
        bits_set = []
        for i in range(nc):
            if mask & (1 << i):
                bits_set.append(i)
        indep = True
        for i in range(len(bits_set)):
            for j in range(i+1, len(bits_set)):
                if adj[bits_set[i]][bits_set[j]]:
                    indep = False
                    break
            if not indep:
                break
        if indep:
            alpha[len(bits_set)] += 1

    results = {}
    for x in x_vals:
        val = sum(alpha[k] * x**k for k in range(max_size + 1))
        results[x] = val

    return results, alpha, cycles

# ======================================================================
# PART 1: Verify I(CG, 2) = H at n=5 (exhaustive) and n=7 (sample)
# ======================================================================
print("=" * 70)
print("PART 1: VERIFY I(CG, 2) = H WITH PROPER CYCLE ENUMERATION")
print("=" * 70)

# n=5 exhaustive
n = 5
tb = n*(n-1)//2
mismatches = 0
t0 = time.time()
for bits in range(1 << tb):
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)
    ix, alpha, cyc = independence_polynomial_exact(A, n, [2])
    if ix[2] != H:
        mismatches += 1
        if mismatches <= 3:
            print(f"  MISMATCH at bits={bits}: H={H}, I(2)={ix[2]}, α={alpha[:5]}, #cycles={len(cyc)}")
dt = time.time() - t0
print(f"  n=5: mismatches = {mismatches}/{1 << tb}, time = {dt:.1f}s")

# n=7 sample
n = 7
tb = n*(n-1)//2
np.random.seed(42)
mismatches = 0
t0 = time.time()
for trial in range(50):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)
    ix, alpha, cyc = independence_polynomial_exact(A, n, [2, 3, 6])
    if ix[2] != H:
        mismatches += 1
        c3 = sum(1 for c in cyc if len(c) == 3)
        c5 = sum(1 for c in cyc if len(c) == 5)
        c7 = sum(1 for c in cyc if len(c) == 7)
        print(f"  MISMATCH trial {trial}: H={H}, I(2)={ix[2]}, c3={c3}, c5={c5}, c7={c7}")
    if trial % 10 == 0:
        dt = time.time() - t0
        print(f"  trial {trial}: {dt:.1f}s, I(2)={ix[2]}, H={H}")

dt = time.time() - t0
print(f"  n=7: mismatches = {mismatches}/50, time = {dt:.1f}s")

# ======================================================================
# PART 2: PRECISE CONVERGENCE RATES
# ======================================================================
print("\n" + "=" * 70)
print("PART 2: PRECISE CONVERGENCE RATE ANALYSIS")
print("=" * 70)

def knacci_root(k):
    coeffs = [-1] * k + [1]
    roots = np.roots(coeffs[::-1])
    real_pos = [r.real for r in roots if abs(r.imag) < 1e-10 and r.real > 1]
    return max(real_pos) if real_pos else None

def knacci_subdominant(k):
    coeffs = [-1] * k + [1]
    roots = np.roots(coeffs[::-1])
    # Sort by magnitude, descending
    roots_sorted = sorted(roots, key=lambda r: abs(r), reverse=True)
    # Dominant is roots_sorted[0]
    return roots_sorted[1] if len(roots_sorted) > 1 else None

def weighted_knacci_root(k):
    coeffs = [1]
    for j in range(k):
        coeffs.append(-2**j)
    roots = np.roots(coeffs)
    real_pos = [r.real for r in roots if abs(r.imag) < 1e-10 and r.real > 1]
    return max(real_pos) if real_pos else None

def weighted_knacci_subdominant(k):
    coeffs = [1]
    for j in range(k):
        coeffs.append(-2**j)
    roots = np.roots(coeffs)
    roots_sorted = sorted(roots, key=lambda r: abs(r), reverse=True)
    return roots_sorted[1] if len(roots_sorted) > 1 else None

print(f"\n  {'k':>3s}  {'φ_k-2':>14s}  {'(φ_{k}-2)·2^k':>14s}  {'|sub/dom|':>10s}  {'ψ_k-3':>14s}  {'(ψ_k-3)·(3/2)^k':>16s}  {'|sub/dom|':>10s}")
print(f"  {'─'*3}  {'─'*14}  {'─'*14}  {'─'*10}  {'─'*14}  {'─'*16}  {'─'*10}")

for k in range(3, 20):
    phi_k = knacci_root(k)
    psi_k = weighted_knacci_root(k)
    sub_k = knacci_subdominant(k)
    wsub_k = weighted_knacci_subdominant(k)

    if phi_k is None or psi_k is None:
        continue

    g2 = phi_k - 2
    g3 = psi_k - 3
    # Normalized: multiply by (dominant/subdominant)^k to get the limit constant
    norm2 = g2 * 2**k
    norm3 = g3 * (3/2)**k

    sub_ratio = abs(sub_k) / phi_k if sub_k is not None else float('nan')
    wsub_ratio = abs(wsub_k) / psi_k if wsub_k is not None else float('nan')

    mark = " ←" if k in [7, 8] else ""
    print(f"  {k:3d}  {g2:14.10f}  {norm2:14.6f}  {sub_ratio:10.6f}  {g3:14.10f}  {norm3:16.6f}  {wsub_ratio:10.6f}{mark}")

# ======================================================================
# PART 3: THE RATIO gap₃/gap₂ — CONVERGES TO?
# ======================================================================
print("\n" + "=" * 70)
print("PART 3: RATIO gap₃/gap₂ CONVERGENCE")
print("=" * 70)

for k in range(3, 25):
    phi_k = knacci_root(k)
    psi_k = weighted_knacci_root(k)
    if phi_k is None or psi_k is None:
        continue
    g2 = phi_k - 2
    g3 = psi_k - 3
    ratio = g3 / g2 if abs(g2) > 1e-15 else float('nan')
    # What is this converging to?
    # If g2 ~ C2*(1/2)^k and g3 ~ C3*(2/3)^k, then
    # g3/g2 ~ (C3/C2) * (4/3)^k → ∞
    # So ratio should GROW, not converge!
    print(f"  k={k:3d}: g3/g2 = {ratio:12.4f}, g3/g2 * (3/4)^k = {ratio * (3/4)**k:.6f}")

# ======================================================================
# PART 4: I(CG, φ_k) at special evaluation points
# ======================================================================
print("\n" + "=" * 70)
print("PART 4: I(CG, x) AT k-NACCI ROOTS vs INTEGERS")
print("=" * 70)

n = 5
tb = n*(n-1)//2

phi_5 = knacci_root(5)
psi_5 = weighted_knacci_root(5)
print(f"\n  φ_5 = {phi_5:.10f} (k-nacci root for k=5)")
print(f"  ψ_5 = {psi_5:.10f} (weighted k-nacci root for k=5)")

x_vals_float = [1.0, phi_5, 2.0, psi_5, 3.0, 6.0]

# Sample tournaments
for bits in [0b0001000, 0b0011100, 0b1001010]:
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)
    cycles = find_all_directed_cycles(A, n)
    nc = len(cycles)

    adj = [[False]*nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if cycles[i] & cycles[j]:
                adj[i][j] = True
                adj[j][i] = True

    alpha = [0] * (nc + 1)
    for mask in range(1 << nc):
        bits_list = [i for i in range(nc) if mask & (1 << i)]
        indep = all(not adj[bits_list[a]][bits_list[b]]
                    for a in range(len(bits_list))
                    for b in range(a+1, len(bits_list)))
        if indep:
            alpha[len(bits_list)] += 1

    print(f"\n  bits={bits:010b}, H={H}, α={alpha[:4]}")
    for x in x_vals_float:
        val = sum(alpha[k] * x**k for k in range(nc + 1))
        label = {1.0: "I(1)", phi_5: f"I(φ₅)", 2.0: "I(2)=H",
                 psi_5: f"I(ψ₅)", 3.0: "I(3)", 6.0: "I(6)"}
        print(f"    {label.get(x, f'I({x:.4f})')}: {val:.4f}")

# ======================================================================
# PART 5: THE DEEP NUMEROLOGY — WHY 2 AND 3
# ======================================================================
print("\n" + "=" * 70)
print("PART 5: THE DEEP STRUCTURE — WHY 2 AND 3")
print("=" * 70)

print("""
  OBSERVATION: In the Jacobsthal tower, x = k(k-1) gives integer root k.
  The DENOMINATORS of J_{k(k-1)}(n) = (k^n - (1-k)^n)/(2k-1) are:
    k=2: denom=3    ↔ 3-cycles
    k=3: denom=5    ↔ 5-cycles
    k=4: denom=7    ↔ 7-cycles

  So the denominator 2k-1 matches the odd cycle length!

  QUESTION: Is there a structural reason?

  Consider: the number of directed k-cycles in K_n is (k-1)!/2 per vertex set.
  For 3-cycles: (3-1)!/2 = 1    → each triple has 0 or 1 3-cycle
  For 5-cycles: (5-1)!/2 = 12   → each 5-set can have 0-12 5-cycles
  For 7-cycles: (7-1)!/2 = 360  → each 7-set can have 0-360 7-cycles

  The Jacobsthal denominator 2k-1 is also the number of ADJACENT pairs
  in a (2k-1)-gon. For a k-cycle, 2k-1 edges... no, that's not right.

  Actually: for a DIRECTED cycle of length m, the sign in the cycle cover
  interpretation is (-1)^{m-1}. For ODD m: sign = +1. For EVEN m: sign = -1.

  The OCF says: H = Σ_S 2^|S| · #{independent sets of size |S|}
  The evaluation at x=2 means: each cycle in an independent set contributes
  a factor of 2. The factor of 2 is EXACTLY the golden ratio cubed divided by...
  no, 2 = 2.

  THE SIMPLEST EXPLANATION:
    2 is special because a tournament edge either points LEFT or RIGHT.
    This binary choice gives factor 2 per independent cycle.
    OCF: H = I(Ω, 2) counts "labeled" cycle packings with 2 labels.

    3 is special because:
    - It's 2+1 (the three states: left, right, neutral)
    - It's the Jacobsthal denominator for x=2
    - J_2(n) = (2^n-(-1)^n)/3 makes 3 the bridge between 2 and -1

    The 3/2 ratio in lambda fibers: I₃/H = 3/2 + O(α₂)
    For α₂=0: I₃/H = (1+3α₁)/(1+2α₁) → 3/2 as α₁ → ∞

  So 3/2 is the ASYMPTOTIC ratio of the two evaluation points!
""")

# Verify: I(3)/I(2) → 3/2 for large α₁, α₂=0
print("  I₃/H ratio as function of α₁ (with α₂=0):")
for a1 in [1, 2, 5, 10, 20, 50, 100]:
    ratio = (1 + 3*a1) / (1 + 2*a1)
    print(f"    α₁={a1:3d}: I₃/H = {ratio:.6f} (3/2 - ε = {1.5 - ratio:.6f})")

print(f"\n  I₃/H → 3/2 = {3/2} with error ~ 1/(4α₁)")

# Now: I(6)/I(2) as function of α₁ (α₂=0)
print("\n  I₆/H ratio as function of α₁ (with α₂=0):")
for a1 in [1, 2, 5, 10, 20, 50]:
    ratio = (1 + 6*a1) / (1 + 2*a1)
    print(f"    α₁={a1:3d}: I₆/H = {ratio:.6f} → 6/2=3")

print(f"\n  I₆/H → 6/2 = 3 as α₁ → ∞ (α₂=0)")
print(f"  I₆/H → 36/4 = 9 as α₂ dominates (α₁=0)")

# ======================================================================
# PART 6: THE n=7,8 THRESHOLD IN k-NACCI TERMS
# ======================================================================
print("\n" + "=" * 70)
print("PART 6: THE n=7,8 THRESHOLD — k-NACCI INTERPRETATION")
print("=" * 70)

# At n=7: α₁ = c3 + c5 + c7 (sum of all odd cycle counts)
# α₂ = #{disjoint pairs} = only (3,3) type
# I(CG, 2) = 1 + 2α₁ + 4α₂

# At n=8: α₂ gets (3,5) type
# I(CG, 2) = 1 + 2α₁ + 4(α₂^{33} + α₂^{35})

# The k-nacci at k=7 approaches 2 with gap 0.008
# This means: at k=7, the recurrence x_n = x_{n-1} + ... + x_{n-7}
# has dominant root 1.992, which is "almost 2"

# BUT: for tournaments at n=7, the "almost" matters because
# the gap 2 - φ_7 ≈ 0.008 is the correction term

# At k=8, the gap halves to 0.004
# The HALVING is exact in the limit — each new k-nacci term
# adds one more "digit" of precision to the approximation of 2

# For weighted k-nacci: each new weighted term adds one more
# "base-3/2 digit" of precision to approximation of 3

print(f"""
  AT k=7 (n=7):
    φ_7 = {knacci_root(7):.10f} (gap = {2 - knacci_root(7):.10f})
    ψ_7 = {weighted_knacci_root(7):.10f} (gap = {3 - weighted_knacci_root(7):.10f})

  AT k=8 (n=8):
    φ_8 = {knacci_root(8):.10f} (gap = {2 - knacci_root(8):.10f})
    ψ_8 = {weighted_knacci_root(8):.10f} (gap = {3 - weighted_knacci_root(8):.10f})

  TRANSITION RATIOS:
    gap₂(8)/gap₂(7) = {(2-knacci_root(8))/(2-knacci_root(7)):.10f} → 1/2
    gap₃(8)/gap₃(7) = {(3-weighted_knacci_root(8))/(3-weighted_knacci_root(7)):.10f} → 2/3

  THE CONNECTION:
    At n=7: α₂ is single-type → I(CG,2) is "almost determined" by α₁
    At n=8: α₂ splits → I(CG,2) gains a new degree of freedom

    The k-nacci gap halving at k=7→8 mirrors this:
    the 7→8 transition adds one "level" of structure
    (from single-type to cross-type packing)
    just as k=7→8 adds one "level" of recurrence terms
    (halving the approximation error to 2).

  NUMEROLOGICAL: 2³-1=7 (Mersenne prime), 2³=8
    The exponent 3 = number of vertices in smallest odd cycle.
    Is this coincidence or deeper?
""")

# ======================================================================
# PART 7: CONVERGENCE CONSTANTS
# ======================================================================
print("=" * 70)
print("PART 7: EXACT CONVERGENCE CONSTANTS")
print("=" * 70)

# If φ_k - 2 ~ C₂ · (1/2)^k, what is C₂?
print("\n  Estimating C₂ from φ_k - 2 = C₂ · (1/2)^k:")
for k in range(5, 20):
    phi_k = knacci_root(k)
    c2_est = (phi_k - 2) * 2**k
    print(f"    k={k:2d}: C₂ ≈ {c2_est:.8f}")

print("\n  C₂ → ? (appears to converge)")

# Same for weighted
print("\n  Estimating C₃ from ψ_k - 3 = C₃ · (2/3)^k:")
for k in range(5, 20):
    psi_k = weighted_knacci_root(k)
    c3_est = (psi_k - 3) * (3/2)**k
    print(f"    k={k:2d}: C₃ ≈ {c3_est:.8f}")

print("\n  C₃ → ? (appears to converge)")

print("\nDone.")
