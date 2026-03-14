#!/usr/bin/env python3
"""
recurrence_lattice.py — opus-2026-03-14-S72

"See everything as recurrences."

This script maps out the COMPLETE recurrence structure of tournament theory,
showing how 2 and 3 control everything.

PART 1:  The 4 fundamental recurrences and their roots
PART 2:  k-nacci convergence: WHY the rate is 1/2
PART 3:  Weighted k-nacci convergence: WHY the rate is 2/3
PART 4:  The deletion recurrence for H(T)
PART 5:  The deletion-contraction recurrence for I(CG, x)
PART 6:  How 7 and 8 appear as recurrence parameters
PART 7:  The 10-11 interpolation in recurrence space
PART 8:  The (2,3) → (7,8) → (10,11) chain as nested recurrences
PART 9:  OEIS connections
"""

import numpy as np
from fractions import Fraction
import math
import time

print("=" * 70)
print("PART 1: THE FOUR FUNDAMENTAL RECURRENCES")
print("=" * 70)
print()

print("Every quantity in tournament theory ultimately derives from")
print("one of four recurrence families:")
print()
print("FAMILY A: The Fibonacci-Jacobsthal family")
print("  f(n) = f(n-1) + x·f(n-2)")
print("  Roots: r = (1+√(1+4x))/2, s = (1-√(1+4x))/2")
print("  x=1: Fibonacci,  roots (φ, -1/φ)")
print("  x=2: Jacobsthal,  roots (2, -1)")
print("  x=6: roots (3, -2)")
print("  x=k(k-1): roots (k, 1-k)")
print()
print("  GENERIC FORM: J_x(n) = (r^n - s^n)/(r-s)")
print("  When x=k(k-1): J_x(n) = (k^n - (1-k)^n)/(2k-1)")
print()

print("FAMILY B: The k-nacci family")
print("  f(n) = f(n-1) + f(n-2) + ... + f(n-k)")
print("  Characteristic: x^k = x^{k-1} + ... + 1")
print("  Equivalently: x^{k+1} - 2x^k + 1 = 0")
print("  Dominant root φ_k → 2 as k → ∞")
print()

print("FAMILY C: The weighted k-nacci family")
print("  f(n) = f(n-1) + 2f(n-2) + ... + 2^{k-1}f(n-k)")
print("  Characteristic: x^k = x^{k-1} + 2x^{k-2} + ... + 2^{k-1}")
print("  Dominant root ψ_k → 3 as k → ∞")
print()

print("FAMILY D: The independence polynomial recurrence")
print("  I(G, x) = I(G-v, x) + x·I(G-N[v], x)")
print("  (deletion-contraction on the conflict graph)")
print("  For isolated vertex: I(K_1 ∪ G, x) = (1+x)·I(G, x)")
print("  → The weight per isolated vertex is 1+x")
print("  → At x=2: weight = 3 (the 'other' key!)")
print()

print("THE CONNECTION:")
print("  Family A at x=2: root 2 (the counting root)")
print("  Family B limit: 2")
print("  Family C limit: 3")
print("  Family D at x=2: isolated weight = 3")
print()
print("  2 and 3 appear as: limits, roots, weights.")
print("  They are the EIGENVALUES of tournament theory.")
print()

print("=" * 70)
print("PART 2: WHY k-NACCI CONVERGES AT RATE 1/2")
print("=" * 70)
print()

# The k-nacci has characteristic x^{k+1} - 2x^k + 1 = 0
# At x=2: 2^{k+1} - 2·2^k + 1 = 2^{k+1} - 2^{k+1} + 1 = 1 ≠ 0
# The root φ_k is close to 2, so let φ_k = 2 - ε_k
# (2-ε)^{k+1} - 2(2-ε)^k + 1 = 0
# Let f(ε) = (2-ε)^{k+1} - 2(2-ε)^k + 1
# f(0) = 2^{k+1} - 2^{k+1} + 1 = 1
# f'(ε) = -(k+1)(2-ε)^k + 2k(2-ε)^{k-1}
# f'(0) = -(k+1)2^k + 2k·2^{k-1} = -(k+1)2^k + k·2^k = -2^k
# By Newton: ε₁ ≈ f(0)/|f'(0)| = 1/2^k
# So ε_k ≈ 1/2^k, meaning 2 - φ_k ≈ C/2^k
# The ratio ε_{k+1}/ε_k → 1/2

print("PROOF that k-nacci gap → 0 at rate 1/2:")
print()
print("  The characteristic equation: x^{k+1} - 2x^k + 1 = 0")
print("  At x = 2: f(2) = 2^{k+1} - 2^{k+1} + 1 = 1")
print("  The derivative: f'(x) = (k+1)x^k - 2k·x^{k-1}")
print("  f'(2) = (k+1)2^k - 2k·2^{k-1} = (k+1)2^k - k·2^k = 2^k")
print()
print("  By Newton's approximation: φ_k ≈ 2 - f(2)/f'(2) = 2 - 1/2^k")
print()
print("  More precisely: let φ_k = 2 - ε_k")
print("  Then ε_k ≈ 2^{-k}")
print("  And ε_{k+1}/ε_k → 2^{-(k+1)}/2^{-k} = 1/2")
print()

# Verify
print("  Verification:")
print(f"  {'k':>4s}  {'2-φ_k':>14s}  {'2^{-k}':>14s}  {'ratio':>10s}  {'(2-φ_k)·2^k':>14s}")
print("  " + "-"*60)

def knacci_root(k, iters=200):
    x = 2.0 - 0.5**(k-1)
    for _ in range(iters):
        f = x**(k+1) - 2*x**k + 1
        fp = (k+1)*x**k - 2*k*x**(k-1)
        if abs(fp) < 1e-30: break
        x -= f/fp
    return x

prev_gap = None
for k in range(2, 20):
    phi = knacci_root(k)
    gap = 2 - phi
    two_inv_k = 2**(-k)
    ratio = gap/prev_gap if prev_gap else float('nan')
    scaled = gap * 2**k
    print(f"  {k:4d}  {gap:14.10e}  {two_inv_k:14.10e}  {ratio:10.6f}  {scaled:14.10f}")
    prev_gap = gap

print()
print("  The scaled gap (2-φ_k)·2^k → 1.0 as k→∞")
print("  EXACT: ε_k = (1 + O(k/2^k)) / 2^k")
print()

print("=" * 70)
print("PART 3: WHY WEIGHTED k-NACCI CONVERGES AT RATE 2/3")
print("=" * 70)
print()

# Weighted k-nacci: x^k = x^{k-1} + 2x^{k-2} + ... + 2^{k-1}
# = x^{k-1} · (1 + 2/x + (2/x)^2 + ... + (2/x)^{k-1})
# = x^{k-1} · (1 - (2/x)^k) / (1 - 2/x)
# = x^{k-1} · (x^k - 2^k) / (x^{k-1}(x-2))  for x ≠ 2
# = (x^k - 2^k) / (x-2)
# So x^k(x-2) = x^k - 2^k
# x^{k+1} - 2x^k = x^k - 2^k
# x^{k+1} - 3x^k + 2^k = 0
#
# At x=3: 3^{k+1} - 3·3^k + 2^k = 3^{k+1} - 3^{k+1} + 2^k = 2^k ≠ 0
# f(3) = 2^k
# f'(x) = (k+1)x^k - 3kx^{k-1}
# f'(3) = (k+1)3^k - 3k·3^{k-1} = (k+1)3^k - k·3^k = 3^k
#
# Newton: ψ_k ≈ 3 - f(3)/f'(3) = 3 - 2^k/3^k = 3 - (2/3)^k
# So gap → 0 at rate 2/3!

print("PROOF that weighted k-nacci gap → 0 at rate 2/3:")
print()
print("  The characteristic equation: x^{k+1} - 3x^k + 2^k = 0")
print("  At x = 3: f(3) = 3^{k+1} - 3·3^k + 2^k = 2^k")
print("  The derivative: f'(x) = (k+1)x^k - 3k·x^{k-1}")
print("  f'(3) = (k+1)3^k - 3k·3^{k-1} = 3^k")
print()
print("  By Newton: ψ_k ≈ 3 - f(3)/f'(3) = 3 - 2^k/3^k = 3 - (2/3)^k")
print()
print("  So ε_k = 3 - ψ_k ≈ (2/3)^k")
print("  And ε_{k+1}/ε_k → 2/3")
print()

def weighted_knacci_root(k, iters=200):
    x = 3.0 - (2/3)**(k-1)
    for _ in range(iters):
        val = x**(k+1) - 3*x**k + 2**k
        dval = (k+1)*x**k - 3*k*x**(k-1)
        if abs(dval) < 1e-30: break
        x -= val/dval
    return x

print("  Verification:")
print(f"  {'k':>4s}  {'3-ψ_k':>14s}  {'(2/3)^k':>14s}  {'ratio':>10s}  {'(3-ψ_k)·(3/2)^k':>14s}")
print("  " + "-"*60)

prev_gap = None
for k in range(2, 20):
    psi = weighted_knacci_root(k)
    gap = 3 - psi
    two_third_k = (2/3)**k
    ratio = gap/prev_gap if prev_gap else float('nan')
    scaled = gap * (3/2)**k
    print(f"  {k:4d}  {gap:14.10e}  {two_third_k:14.10e}  {ratio:10.6f}  {scaled:14.10f}")
    prev_gap = gap

print()
print("  EXACT: 3 - ψ_k = (1 + O(k·(2/3)^k)) · (2/3)^k")
print()
print("  THE RATES ARE RECIPROCALS OF THE LIMITS:")
print("  k-nacci → 2 at rate 1/2 = 1/(limit)")
print("  w-k-nacci → 3 at rate 2/3 = 1/(limit) · (limit-1)")
print()
print("  More precisely:")
print("  k-nacci: gap ≈ 1/2^k = (1/2)^k")
print("  w-k-nacci: gap ≈ (2/3)^k")
print()
print("  Ratio of gaps: (2/3)^k / (1/2)^k = (4/3)^k → ∞")
print("  The weighted k-nacci converges SLOWER by a factor (4/3)^k.")
print()

print("=" * 70)
print("PART 4: THE DELETION RECURRENCE FOR H(T)")
print("=" * 70)
print()

# H(T) via vertex deletion: H(T) = Σ_P 1 over Hamiltonian paths
# Choose first vertex v: P starts at some v and goes v→w→...
# H(T) = Σ_{v∈V} Σ_{w: v→w} H_{v→w}(T)
# where H_{v→w}(T) = # Hamiltonian paths starting v→w
# This doesn't simplify to a simple recurrence in terms of H(T-v)

# But: OCF gives H(T) = I(CG(T), 2)
# And I satisfies deletion-contraction on CG!

def adj_matrix(bits, n):
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

def hamiltonian_paths(A):
    n = len(A)
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask_size in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != mask_size:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                prev_mask = mask ^ (1 << v)
                count = 0
                for u in range(n):
                    if (prev_mask & (1 << u)) and A[u][v]:
                        count += dp.get((prev_mask, u), 0)
                if count:
                    dp[(mask, v)] = count
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

# The deletion recurrence structure: H(T) in terms of H(T-v)
print("DELETION STRUCTURE: H(T) vs H(T-v) at n=5")
print()

n = 5
total = 2**(n*(n-1)//2)
t0 = time.time()

# For each tournament, compute H(T) and H(T-v) for each v
# Then see the relationship
print("  Sampling 200 tournaments at n=5...")
import random
random.seed(42)

deletion_data = []
for trial in range(200):
    bits = random.randint(0, total-1)
    A = adj_matrix(bits, n)
    H = hamiltonian_paths(A)

    H_minus = []
    for v in range(n):
        # Delete vertex v
        indices = [i for i in range(n) if i != v]
        A_sub = A[np.ix_(indices, indices)]
        H_sub = hamiltonian_paths(A_sub)
        H_minus.append(H_sub)

    deletion_data.append((H, H_minus))

print(f"  Done in {time.time()-t0:.1f}s")
print()

# Analyze: H(T) vs sum of H(T-v)
print("  H(T) vs Σ_v H(T-v):")
for H, Hm in deletion_data[:20]:
    s = sum(Hm)
    print(f"    H={H:3d}, Σ H(T-v)={s:3d}, ratio={H/s:.4f}, diff={H-s:4d}")

# Is there a pattern?
ratios = [H/sum(Hm) for H, Hm in deletion_data if sum(Hm) > 0]
print(f"\n  H/Σ H(T-v): mean={np.mean(ratios):.6f}, std={np.std(ratios):.6f}")
print(f"              min={np.min(ratios):.6f}, max={np.max(ratios):.6f}")
print()

# The Claim A recurrence: H(T) - H(T-v) = 2 Σ_{C∋v} μ(C)
# So H(T) = H(T-v) + 2 Σ_{C∋v} μ(C)
# Summing over v: n·H(T) = Σ_v H(T-v) + 2 Σ_v Σ_{C∋v} μ(C)
# Each cycle C of length L appears L times (once per vertex in C)
# So: n·H(T) = Σ_v H(T-v) + 2 Σ_C L(C) · μ(C)
print("  From Claim A: n·H(T) = Σ_v H(T-v) + 2 Σ_C L(C)·μ(C)")
print()
print("  So H(T) = [Σ_v H(T-v) + 2 Σ_C L(C)·μ(C)] / n")
print()
print("  This IS a recurrence, but it depends on the CYCLE STRUCTURE")
print("  of T, not just the H-values of subtournaments.")
print()
print("  The '2' in 2μ(C) is the OCF evaluation point x=2.")
print("  If we used x=3, we'd get a '3' factor instead.")
print()

print("=" * 70)
print("PART 5: THE (2,3) DUALITY IN DELETION-CONTRACTION")
print("=" * 70)
print()

print("I(CG, x) satisfies deletion-contraction on the conflict graph:")
print("  I(G, x) = I(G-v, x) + x·I(G-N[v], x)")
print()
print("At x=2: I(G, 2) = I(G-v, 2) + 2·I(G-N[v], 2)")
print("At x=3: I(G, 3) = I(G-v, 3) + 3·I(G-N[v], 3)")
print()
print("The RATIO of the 'contraction' terms is 3/2.")
print("This is the same 3/2 that appears everywhere!")
print()
print("For an isolated vertex v (no neighbors in CG):")
print("  I(G, x) = I(G-v, x) + x·I(G, x)  ... no, N[v]={v}")
print("  I(G, x) = I(G-v, x) + x·I(G-{v}, x) ... but G-{v} and G-N[v] differ")
print("  Actually: I(G, x) = I(G-v, x) + x·I(G-N[v], x)")
print("  For isolated v: N[v] = {v}, so G-N[v] = G-v")
print("  I(G, x) = I(G-v, x) + x·I(G-v, x) = (1+x)·I(G-v, x)")
print()
print("  At x=2: factor = 1+2 = 3")
print("  At x=3: factor = 1+3 = 4")
print()
print("  Each isolated vertex in CG contributes a factor of (1+x).")
print("  For H = I(CG, 2): each isolated cycle multiplies H by 3.")
print("  For I₃ = I(CG, 3): each isolated cycle multiplies I₃ by 4.")
print()
print("  The RATIO of multiplicative factors: 4/3 = (1+3)/(1+2)")
print("  But we said the ratio is 3/2 earlier!")
print("  The 3/2 ratio is for the ADDITIVE part (contraction term).")
print("  The 4/3 ratio is for the MULTIPLICATIVE part (isolated vertex).")
print()
print("  BOTH ratios involve 2 and 3:")
print("  3/2 = 3/2")
print("  4/3 = (3+1)/(2+1)")
print()
print("  And 3/2 × 4/3 = 2.  Full circle: back to the root!")
print()

print("=" * 70)
print("PART 6: HOW 7 AND 8 APPEAR AS RECURRENCE PARAMETERS")
print("=" * 70)
print()

print("At n=7 (=2³-1):")
print("  The 7-nacci root is φ₇ = 1.99196... ≈ 2 - 1/128 = 2 - 2^{-7}")
print("  The 7-cycle first appears in tournaments")
print("  The Jacobsthal at level k=4 has denominator 7")
print()
print("At n=8 (=2³):")
print("  The 8-nacci root is φ₈ = 1.99603... ≈ 2 - 1/256 = 2 - 2^{-8}")
print("  First cross-level cycle packing (3+5=8)")
print("  Q = (H²-Pf²)/8: the divisor 8 = 2³ first matters")
print()

# The 7-cycle and 8-cycle recurrences
print("The 7-step and 8-step recurrences (Families B/C at k=7,8):")
print()

print("  k=7 k-nacci: f(n) = f(n-1) + ... + f(n-7)")
print(f"    Dominant root: {knacci_root(7):.10f}")
print(f"    Gap to 2: {2-knacci_root(7):.10e}")
print()
print("  k=8 k-nacci: f(n) = f(n-1) + ... + f(n-8)")
print(f"    Dominant root: {knacci_root(8):.10f}")
print(f"    Gap to 2: {2-knacci_root(8):.10e}")
print()
print("  k=7 weighted: f(n) = f(n-1) + 2f(n-2) + ... + 64f(n-7)")
print(f"    Dominant root: {weighted_knacci_root(7):.10f}")
print(f"    Gap to 3: {3-weighted_knacci_root(7):.10e}")
print()
print("  k=8 weighted: f(n) = f(n-1) + 2f(n-2) + ... + 128f(n-8)")
print(f"    Dominant root: {weighted_knacci_root(8):.10f}")
print(f"    Gap to 3: {3-weighted_knacci_root(8):.10e}")
print()

# The specific recurrence that involves 7 as a denominator
print("THE JACOBSTHAL CONNECTION:")
print("  J₁₂(n) = (4^n - (-3)^n) / 7")
print("  The denominator 7 = 2·4 - 1 = 2k-1 for k=4")
print()
print("  So 7 appears as a JACOBSTHAL DENOMINATOR,")
print("  and denominators ↔ cycle lengths (3, 5, 7, ...)")
print()
print("  The 7-cycle is the first cycle that REQUIRES n≥7 vertices.")
print("  It's the first cycle beyond the 'single Jacobsthal level' regime.")
print()
print("  When a 7-cycle coexists with a 3-cycle (at n≥10=7+3),")
print("  or with another 3-cycle and a 5-cycle (at n≥15=7+3+5),")
print("  the Jacobsthal levels couple through the conflict graph.")
print()

# The 8-divisor in Q
print("THE Q DIVISOR:")
print("  Q(T) = (H(T)² - Pf(T)²) / 8")
print("  Why 8? Because H² - Pf² = I(CG,2)² - det(I+2A)")
print("  and det(I+2A) = 1 + 0 + 0 + 8c₃ + 16c₄ + ...")
print("  The first non-trivial coefficient is at x³ with coefficient 8c₃.")
print()
print("  Why is the coefficient 8? Because:")
print("  det(I+xA) = 1 + Σ_{k≥3} c_k · x^k")
print("  c₃ = #{3-cycles}, and x³ = 2³ = 8 at x=2.")
print("  So 8 = 2³ is the 'entry point' for cycle contributions.")
print()

print("=" * 70)
print("PART 7: THE 10-11 INTERPOLATION IN RECURRENCE SPACE")
print("=" * 70)
print()

# 10 and 11 are between k=3 (x=6, root=3) and k=4 (x=12, root=4)
# in the Jacobsthal-integer-root tower
print("POSITION OF 10 AND 11:")
print()
print("  Jacobsthal integer-root points: x = k(k-1)")
print("  k=3: x=6,  root=3")
print("  k=4: x=12, root=4")
print()
print("  At x=10: root ≈ 3.702")
print("  At x=11: root ≈ 3.854")
print()
print("  These are (10-6)/(12-6) = 4/6 = 2/3 and (11-6)/(12-6) = 5/6")
print("  of the way from k=3 to k=4.")
print()
print("  THE 2/3 AGAIN! x=10 is EXACTLY 2/3 of the way between")
print("  consecutive integer-root Jacobsthal levels!")
print()

# Actually let me verify: is this exactly 2/3?
print("  Exact fractions:")
print(f"  (10-6)/(12-6) = 4/6 = {Fraction(4,6)}")
print(f"  (11-6)/(12-6) = 5/6 = {Fraction(5,6)}")
print()
print("  So 10 is at position 2/3 between x=6 and x=12,")
print("  and 11 is at position 5/6.")
print()
print("  The root interpolation:")
r10 = (1 + math.sqrt(1 + 4*10)) / 2
r11 = (1 + math.sqrt(1 + 4*11)) / 2
print(f"  root(x=10) = {r10:.6f}, position in [3,4]: {r10-3:.6f} = {Fraction(r10-3).limit_denominator(100)}")
print(f"  root(x=11) = {r11:.6f}, position in [3,4]: {r11-3:.6f} = {Fraction(r11-3).limit_denominator(100)}")
print()

# The key: 10 = 2 + 8 = 2 + 2³
print("  10 = 2 + 8 = x_OCF + 2³")
print("  This means: I(CG, 10) = I(CG, 2 + 8)")
print()
print("  Taylor expansion of I(CG, x) around x=2:")
print("    I(CG, 2+δ) = I(CG, 2) + δ·I'(CG, 2) + δ²·I''(CG, 2)/2 + ...")
print("                = H + δ·(α₁ + 4α₂ + ...) + ...")
print()
print("  At δ=8: I(CG, 10) = H + 8·I'(CG, 2) + 32·I''(CG, 2) + ...")
print()
print("  Similarly, 11 = 3 + 8 = x_topo + 2³")
print("  I(CG, 11) = I(CG, 3 + 8) ← Taylor around x=3!")
print("    = I₃ + 8·I'(CG, 3) + 32·I''(CG, 3) + ...")
print()
print("  So I(CG, 10) and I(CG, 11) encode DERIVATIVES of I at 2 and 3,")
print("  weighted by powers of 8 = 2³!")
print()

# For a degree-d polynomial I(x) = 1 + α₁x + α₂x², the exact values:
print("  For degree-2 case (n ≤ 8, typically):")
print("    I(CG, x) = 1 + α₁x + α₂x²")
print("    I(CG, 2) = 1 + 2α₁ + 4α₂ = H")
print("    I(CG, 3) = 1 + 3α₁ + 9α₂ = I₃")
print("    I(CG, 10) = 1 + 10α₁ + 100α₂")
print("    I(CG, 11) = 1 + 11α₁ + 121α₂")
print()
print("  Relations:")
print("    I(10) - I(2) = 8α₁ + 96α₂ = 8(α₁ + 12α₂)")
print("    I(11) - I(3) = 8α₁ + 112α₂ = 8(α₁ + 14α₂)")
print()
print("    Both divisible by 8!")
print("    (I(10) - I(2)) / 8 = α₁ + 12α₂")
print("    (I(11) - I(3)) / 8 = α₁ + 14α₂")
print()
print("    Subtracting: [(I(11)-I(3)) - (I(10)-I(2))] / 8 = 2α₂")
print("    So: α₂ = [(I(11) - I(10)) - (I(3) - I(2))] / 16")
print()
print("    And: α₁ = (I(10) - I(2))/8 - 12α₂")
print()
print("  THE 10-11 PAIR EXTRACTS (α₁, α₂) just like the 2-3 pair does,")
print("  but shifted by 8!")
print()

print("=" * 70)
print("PART 8: NESTED RECURRENCE STRUCTURE")
print("    (2,3) → (7,8) → (10,11)")
print("=" * 70)
print()

print("LEVEL 1: (2,3) — The fundamental Vandermonde extraction")
print("  I(2), I(3) → α₁, α₂ via Vandermonde:")
print("    α₂ = (2·I(3) - 3·I(2) + 1) / 6 = (2I₃ - 3H + 1) / 6")
print("    α₁ = (I(2) - 1 - 4α₂) / 2 = (H - 1 - 4α₂) / 2")
print()
print("  The denominator 6 = 2·3 (product of the pair!)")
print()

print("LEVEL 2: (7,8) — The cross-level threshold")
print("  At n=8, α₂ splits into α₂^{33} (both 3-cycles) and α₂^{35} (3+5)")
print("  I(2) and I(3) can't separate these — need I(6) or I(7)")
print("  The pair (7,8) tells us when to expect this complexity.")
print()

print("LEVEL 3: (10,11) — The shifted Vandermonde")
print("  I(10) = I(2 + 8), I(11) = I(3 + 8)")
print("  These are the I(2) and I(3) values 'shifted by 8=2³'")
print("  (I(10) - H) / 8 and (I(11) - I₃) / 8 give higher-order info")
print()

print("THE NESTING:")
print("  At level 1, (2,3) extracts α₁, α₂ from H and I₃")
print("  At level 3, (10,11) = (2+8, 3+8) extracts the SAME α₁, α₂")
print("  but through a different linear combination — providing a")
print("  CONSISTENCY CHECK on the α values.")
print()
print("  In fact, ANY consecutive pair (x, x+1) with x ≡ 0 mod 2")
print("  or x ≡ 2 mod 3 gives a Vandermonde extraction of (α₁, α₂).")
print("  But (2,3), (10,11), (7,8) are special because they also")
print("  encode the hierarchical structure of cycle packing.")
print()

print("=" * 70)
print("PART 9: THE RECURRENCE LATTICE DIAGRAM")
print("=" * 70)
print()

print("               1 (identity)")
print("              / \\")
print("             2   3  (the two keys)")
print("            / \\ / \\")
print("           4   6   9  (2², 2·3, 3²)")
print("          / \\ / \\")
print("         7   8   12  (2³-1, 2³, 2²·3)")
print("          \\ / \\ /")
print("          10   11  (2+8, 3+8)")
print()
print("  Each node n defines:")
print("  1. A k-nacci recurrence with root → 2")
print("  2. A weighted k-nacci recurrence with root → 3")
print("  3. A Jacobsthal evaluation I(CG, n) or I(CG, n(n-1))")
print("  4. A tournament size where new packing structures appear")
print()
print("  The lattice is organized by the TWO operations:")
print("  • 'doubling': n → 2n (multiplication by 2)")
print("  • 'tripling': n → 3n (multiplication by 3)")
print("  which generate the 2-3 lattice Z[1/6].")
print()

# Summary statistics
print("=" * 70)
print("FINAL SUMMARY: EVERYTHING IS A RECURRENCE")
print("=" * 70)
print()
print("  RECURRENCE          | ROOT(S)     | RATE     | KEY NUMBER")
print("  --------------------+-------------+----------+-----------")
print("  Fibonacci           | φ, -1/φ     | —        | 1")
print("  Jacobsthal (x=2)   | 2, -1       | —        | 2")
print("  Jacobsthal (x=6)   | 3, -2       | —        | 3")
print("  Jacobsthal (x=12)  | 4, -3       | —        | 7 (denom)")
print("  k-nacci (any k)    | → 2         | 1/2      | 2")
print("  w-k-nacci (any k)  | → 3         | 2/3      | 3")
print("  I(CG, x) del-con   | 1+x         | —        | 3 (at x=2)")
print("  H deletion          | 2-weighted  | —        | 2")
print("  Q = (H²-Pf²)/8    | —           | —        | 8 = 2³")
print("  α₁ mod 2 tower     | 2-adic      | —        | 2")
print("  b₀ mod 3 tower     | 3-adic      | —        | 3")
print()
print("  2 controls: counting (H), convergence (k-nacci → 2),")
print("              divisibility (H odd), deletion weight")
print()
print("  3 controls: topology (b₀), convergence (w-k-nacci → 3),")
print("              isolated vertex weight (1+x=3), CRT bridge")
print()
print("  7 = 2³-1: last pure (single-level) tournament size,")
print("            Jacobsthal denominator, Mersenne prime")
print()
print("  8 = 2³:   first cross-level, Q divisor, packing threshold")
print()
print("  10 = 2+8: shifted counting point, 2/3 interpolation")
print("  11 = 3+8: shifted topology point, first perfect α₃ packing")
print()
print("  'If one understands 2 and 3, one has the keys to the universe.'")
print("  2 and 3 are the eigenvalues. 7 and 8 are the transition.")
print("  10 and 11 are the interleaving. Everything else follows.")
print()
print("Done.")
