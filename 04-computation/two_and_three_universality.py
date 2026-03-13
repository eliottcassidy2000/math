#!/usr/bin/env python3
"""
two_and_three_universality.py — opus-2026-03-13-S71c

"If one knows 2 and 3, then we have the keys to the universe."
k-nacci constants → 2, weighted k-nacci constants → 3.

QUESTIONS:
1. What is I(CG, 3) = Σ α_k 3^k? Does it have nice properties?
2. How does I(3) relate to H = I(2)?
3. Is I(3) always odd? Always ≡ something mod powers of 3?
4. det(I + 3A) — what is it?
5. The ratio I(3)/I(2) — is it meaningful?
6. Does I(x) at x=2 and x=3 together determine more than either alone?

The weighted k-nacci recurrence: f(n) = 1·f(n-1) + 2·f(n-2) + 4·f(n-3) + ...
with weights [2^0, 2^1, ..., 2^{k-1}].

The dominant root of this approaches 3. WHY?
Characteristic equation: x^k = 1 + 2x^{-1} + 4x^{-2} + ... + 2^{k-1}x^{-(k-1)}
× x^{k-1}: x^{2k-1} = x^{k-1} + 2x^{k-2} + ... + 2^{k-1}
= Σ_{j=0}^{k-1} 2^j x^{k-1-j}
= x^{k-1} · Σ_{j=0}^{k-1} (2/x)^j
= x^{k-1} · (1 - (2/x)^k) / (1 - 2/x)  [geometric sum]
= x^{k-1} · (x^k - 2^k) / (x^k · (x-2)/x)
= x^{k-1} · x · (x^k - 2^k) / (x^k · (x-2))
= (x^k - 2^k) / (x-2)     ... wait

So x^{2k-1} = (x^k - 2^k)/(x-2) · x^... hmm let me redo this.

Actually: f(n) = Σ_{j=1}^{k} 2^{j-1} f(n-j)
Char poly: x^k = 1 + 2/x + 4/x^2 + ... + 2^{k-1}/x^{k-1}
Multiply by x^{k-1}:
x^{2k-1} = x^{k-1} + 2x^{k-2} + 4x^{k-3} + ... + 2^{k-1}

RHS = Σ_{j=0}^{k-1} 2^j x^{k-1-j}

At large k, for x=3: RHS = Σ 2^j 3^{k-1-j} = 3^{k-1} Σ (2/3)^j → 3^{k-1}/(1-2/3) = 3^k

LHS = 3^{2k-1}. So 3^{2k-1} = 3^k → only at k→∞.

Actually this shows the exact equation is:
x^{2k-1} = Σ_{j=0}^{k-1} 2^j x^{k-1-j}

For the geometric sum: RHS = x^{k-1} (1 - (2/x)^k)/(1 - 2/x)
= x^{k-1} · (x^k - 2^k)/(x^k) · x/(x-2)
= (x^k - 2^k)/(x-2)

So x^{2k-1} = (x^k - 2^k)/(x-2)
x^{2k-1}(x-2) = x^k - 2^k
x^{2k} - 2x^{2k-1} = x^k - 2^k
x^{2k} - 2x^{2k-1} - x^k + 2^k = 0

For large k, x → 3: 3^{2k} - 2·3^{2k-1} - 3^k + 2^k = 0
3^{2k}(1 - 2/3) - 3^k + 2^k = 0
3^{2k}/3 - 3^k + 2^k = 0
3^{2k-1} - 3^k + 2^k = 0

At k→∞: 3^{2k-1} >> 3^k >> 2^k, so this is NOT zero.
So x=3 is not an exact root. Let me check numerically.
"""

import sys, time
import numpy as np
from itertools import combinations, permutations
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

# =====================================================================
print("=" * 70)
print("SECTION 1: WHY WEIGHTED k-NACCI → 3")
print("=" * 70)

def weighted_knacci_constant(k):
    """Dominant root of f(n) = Σ_{j=1}^{k} 2^{j-1} f(n-j)"""
    M = np.zeros((k, k))
    for j in range(k):
        M[0, j] = 2**j
    for i in range(1, k):
        M[i, i-1] = 1
    eigs = np.linalg.eigvals(M)
    return max(abs(e) for e in eigs)

print("\nWeighted k-nacci dominant root:")
for k in range(2, 20):
    r = weighted_knacci_constant(k)
    gap = 3 - r
    print(f"  k={k:2d}: λ = {r:.10f}, 3 - λ = {gap:.10f}, gap ≈ {gap:.2e}")

print("\nThe char equation at k→∞:")
print("  x^{2k-1} = (x^k - 2^k)/(x-2)")
print("  Rearranging: x^{2k}(1-2/x) = x^k - 2^k")
print("  At x=3: 3^{2k}/3 = 3^k - 2^k → FALSE for finite k")
print("  But the RATIO 3^{2k}/3 / 3^k = 3^{k-1} → ∞")
print("  So x=3 is NOT a root. But λ_dom → 3 from below.")

# Check: what is the EXACT limit?
# The characteristic polynomial at order k is:
# p(x) = x^k - 1 - 2x^{-1} - 4x^{-2} - ... - 2^{k-1}x^{-(k-1)}
# Multiply by x^{k-1}:
# q(x) = x^{2k-1} - x^{k-1} - 2x^{k-2} - ... - 2^{k-1}
# = x^{2k-1} - Σ_{j=0}^{k-1} 2^j x^{k-1-j}

# At x=3: q(3) = 3^{2k-1} - Σ_{j=0}^{k-1} 2^j 3^{k-1-j}
# = 3^{2k-1} - 3^{k-1}(1-(2/3)^k)/(1-2/3)
# = 3^{2k-1} - 3^k(1-(2/3)^k)
# = 3^{2k-1} - 3^k + 2^k

print("\nq(3) for various k:")
for k in range(2, 15):
    q3 = 3**(2*k-1) - 3**k + 2**k
    print(f"  k={k:2d}: q(3) = {q3}")

# =====================================================================
print(f"\n{'='*70}")
print("SECTION 2: I(CG, 3) FOR TOURNAMENTS")
print("=" * 70)

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
    """Count α_k = #{collections of k disjoint directed odd cycles}."""
    # First find ALL directed odd cycles
    cycles = []
    seen = set()

    def find_cycles(path, start, used):
        last = path[-1]
        length = len(path)
        if length >= 3 and length % 2 == 1 and A[last][start]:
            # Found odd cycle
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

    # Now count disjoint collections
    alpha = defaultdict(int)
    alpha[0] = 1  # empty collection

    def count_collections(idx, used, k):
        alpha[k] += 1
        for j in range(idx, len(cycles)):
            if not (cycles[j] & used):
                count_collections(j+1, used | cycles[j], k+1)

    count_collections(0, frozenset(), 0)
    del alpha[0]
    return dict(alpha)

def I_at_x(alpha, x):
    """Compute I(CG, x) = Σ α_k x^k."""
    return 1 + sum(count * x**k for k, count in alpha.items())

# Compute for small n
for n in [5, 6, 7]:
    print(f"\nn={n}:")
    tb = n*(n-1)//2
    np.random.seed(42)

    results = []
    for trial in range(min(100, 1 << tb)):
        if n <= 6:
            bits = trial
        else:
            bits = np.random.randint(0, 1 << tb)

        A = bits_to_adj(bits, n)
        H = count_ham_paths(A, n)
        alpha = count_disjoint_cycles(A, n)

        I2 = I_at_x(alpha, 2)
        I3 = I_at_x(alpha, 3)

        # Also compute det(I + xA)
        I_mat = np.eye(n)
        det2 = int(round(np.linalg.det(I_mat + 2*A)))
        det3 = int(round(np.linalg.det(I_mat + 3*A)))

        results.append((H, I2, I3, det2, det3, alpha))

    # Show first 10
    print(f"  {'H':>5} {'I(2)':>6} {'I(3)':>6} {'det(I+2A)':>10} {'det(I+3A)':>10} {'I3/I2':>8} {'I3 mod 3':>8} {'alpha':>20}")
    for H, I2, I3, det2, det3, alpha in results[:20]:
        ratio = I3/I2 if I2 > 0 else 0
        print(f"  {H:5d} {I2:6d} {I3:6d} {det2:10d} {det3:10d} {ratio:8.3f} {I3%3:8d} {dict(alpha)}")

    # Statistics
    I3_vals = [r[2] for r in results]
    I2_vals = [r[1] for r in results]
    ratios = [r[2]/r[1] for r in results if r[1] > 0]

    print(f"\n  I(2) range: [{min(I2_vals)}, {max(I2_vals)}]")
    print(f"  I(3) range: [{min(I3_vals)}, {max(I3_vals)}]")
    print(f"  I(3)/I(2) range: [{min(ratios):.3f}, {max(ratios):.3f}], mean={np.mean(ratios):.3f}")
    print(f"  I(3) mod 3 distribution: {dict(sorted(defaultdict(int, {k: sum(1 for r in results if r[2]%3==k) for k in range(3)}).items()))}")
    print(f"  I(3) always odd? {all(r[2] % 2 == 1 for r in results)}")

    # Check: is I(3) = det(I+3A) somehow?
    det3_match = sum(1 for r in results if r[2] == r[4])
    print(f"  I(3) = det(I+3A)? matches: {det3_match}/{len(results)}")

# =====================================================================
print(f"\n{'='*70}")
print("SECTION 3: THE (2,3) DETERMINATION QUESTION")
print("=" * 70)

print("\nDo (I(2), I(3)) together determine MORE than either alone?")
n = 7
tb = n*(n-1)//2
np.random.seed(42)

# Check if (H, I3) determines more structure
h_to_c3 = defaultdict(set)
h_i3_to_c3 = defaultdict(set)
h_to_c5 = defaultdict(set)
h_i3_to_c5 = defaultdict(set)

t0 = time.time()
for trial in range(5000):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)

    A2 = A @ A
    c3 = int(np.trace(A2 @ A)) // 3
    c5 = int(np.trace(np.linalg.matrix_power(A, 5))) // 5

    alpha = count_disjoint_cycles(A, n)
    I3 = I_at_x(alpha, 3)

    h_to_c3[H].add(c3)
    h_i3_to_c3[(H, I3)].add(c3)
    h_to_c5[H].add(c5)
    h_i3_to_c5[(H, I3)].add(c5)

dt = time.time() - t0
print(f"  n=7, 5000 samples, {dt:.1f}s")

h_ambig_c3 = sum(1 for v in h_to_c3.values() if len(v) > 1)
hi3_ambig_c3 = sum(1 for v in h_i3_to_c3.values() if len(v) > 1)
h_ambig_c5 = sum(1 for v in h_to_c5.values() if len(v) > 1)
hi3_ambig_c5 = sum(1 for v in h_i3_to_c5.values() if len(v) > 1)

print(f"  H determines c3? Ambiguous: {h_ambig_c3}/{len(h_to_c3)}")
print(f"  (H,I(3)) determines c3? Ambiguous: {hi3_ambig_c3}/{len(h_i3_to_c3)}")
print(f"  H determines c5? Ambiguous: {h_ambig_c5}/{len(h_to_c5)}")
print(f"  (H,I(3)) determines c5? Ambiguous: {hi3_ambig_c5}/{len(h_i3_to_c5)}")

# =====================================================================
print(f"\n{'='*70}")
print("SECTION 4: THE JACOBSTHAL CONNECTION")
print("=" * 70)

print("\nJacobsthal numbers J(n): recurrence J(n) = J(n-1) + 2*J(n-2)")
print("This is EXACTLY the weighted 2-nacci with weights [1, 2]!")
print("Dominant root = 2 (exact!)")
print()

# Jacobsthal: J(0)=0, J(1)=1, J(n)=J(n-1)+2J(n-2)
J = [0, 1]
for i in range(2, 15):
    J.append(J[-1] + 2*J[-2])
print(f"Jacobsthal: {J}")
print(f"J(n) = (2^n - (-1)^n)/3")
for i in range(len(J)):
    formula = (2**i - (-1)**i) // 3
    print(f"  J({i}) = {J[i]}, formula = {formula}, match = {J[i] == formula}")

print(f"\nNote: J(n) = (2^n - (-1)^n)/3")
print(f"  At x=2 (OCF point): 2^n is the DOMINANT term")
print(f"  At x=(-1): (-1)^n is the ALTERNATING correction")
print(f"  The RATIO is 1/3 — there's the 3!")
print(f"  J(n) = 2^n/3 - (-1)^n/3")
print(f"  The 3 in the Jacobsthal formula IS the weighted k-nacci limit!")

# =====================================================================
print(f"\n{'='*70}")
print("SECTION 5: 2 AND 3 AS EIGENVALUES")
print("=" * 70)

print("\nThe Jacobsthal transfer matrix is [[1,2],[1,0]]")
print("Eigenvalues: 2 and -1")
print("Sum of eigenvalues: 2 + (-1) = 1")
print("Product of eigenvalues: 2 × (-1) = -2")
print()

# For the weighted k-nacci with weights [1, 2, 4, ..., 2^{k-1}]:
print("Weighted k-nacci eigenvalues:")
for k in range(2, 8):
    weights = [2**j for j in range(k)]
    M = np.zeros((k, k))
    M[0, :] = weights
    for i in range(1, k):
        M[i, i-1] = 1
    eigs = sorted(np.linalg.eigvals(M), key=lambda x: -abs(x))
    dom = eigs[0].real
    print(f"  k={k}: dominant = {dom:.6f}, all eigenvalues: {[f'{e:.3f}' for e in eigs]}")

# =====================================================================
print(f"\n{'='*70}")
print("SECTION 6: I(x) AS POLYNOMIAL — WHAT ARE THE ROOTS?")
print("=" * 70)

# For each tournament, I(CG, x) = 1 + α₁x + α₂x² + ...
# What are the ROOTS of I(CG, x)?
n = 5
tb = n*(n-1)//2

print(f"\nn={n}: Roots of I(CG, x) for each iso class")
classes = defaultdict(list)
for bits in range(1 << tb):
    A = bits_to_adj(bits, n)
    # Canonical form
    best = None
    for perm in permutations(range(n)):
        form = tuple(int(A[perm[i]][perm[j]]) for i in range(n) for j in range(i+1, n))
        if best is None or form < best:
            best = form
    classes[best].append(bits)

print(f"  {len(classes)} iso classes")
for cf in sorted(classes.keys()):
    A = bits_to_adj(classes[cf][0], n)
    H = count_ham_paths(A, n)
    alpha = count_disjoint_cycles(A, n)

    # Build polynomial coefficients: [1, α₁, α₂, ...]
    max_k = max(alpha.keys()) if alpha else 0
    coeffs = [1] + [alpha.get(k, 0) for k in range(1, max_k + 1)]

    # Roots of I(CG, x) = 0 (polynomial in x)
    if len(coeffs) > 1:
        # coeffs[0] + coeffs[1]*x + coeffs[2]*x^2 + ...
        # numpy.roots wants highest degree first
        poly_coeffs = list(reversed(coeffs))
        roots = np.roots(poly_coeffs)
        roots_str = ', '.join(f'{r:.3f}' for r in sorted(roots, key=lambda x: abs(x)))
    else:
        roots_str = "none"

    print(f"  H={H:3d}, α={dict(alpha)}, poly={coeffs}, roots=[{roots_str}]")

# =====================================================================
print(f"\n{'='*70}")
print("SECTION 7: THE DEEP IDENTITY — 2 AND 3 IN THE PARTITION FUNCTION")
print("=" * 70)

print("""
THE k-NACCI INSIGHT:
  Standard k-nacci constant → 2   (the OCF evaluation point)
  Weighted k-nacci constant → 3   (the "shadow" evaluation point)

JACOBSTHAL BRIDGE:
  J(n) = (2^n - (-1)^n) / 3

  The formula DIVIDES by 3 to connect the 2-world to the integer world.
  2 is the growth rate, 3 is the normalization.

PARTITION FUNCTION VIEW:
  H = I(CG, 2) = Σ α_k 2^k   (evaluate at x=2)
  I(CG, 3) = Σ α_k 3^k         (evaluate at x=3)

  The ratio I(3)/I(2) ≈ 3^ν(G)/2^ν(G) = (3/2)^ν(G) for dominant α_ν

  But ν(G) is the max independent set size in CG (max packing number).
  So the ratio I(3)/I(2) encodes the DIMENSION of the dominant packing!

THE 2-3 DUALITY:
  x=2: ALL k-nacci levels are active (critical point from BELOW)
  x=3: ALL weighted k-nacci levels are active (critical point from BELOW)

  The INTERVAL [2,3] is where the unweighted and weighted towers
  transition. It's the "critical strip" of tournament theory.

  This is analogous to the critical strip Re(s) ∈ [0,1] in the
  Riemann zeta function, where the partition function Σ n^{-s}
  transitions between convergent and divergent behavior.
""")

# Check: I(3)/I(2) vs max packing number
n = 5
print(f"\nn={n}: I(3)/I(2) vs max packing number:")
for cf in sorted(classes.keys()):
    A = bits_to_adj(classes[cf][0], n)
    H = count_ham_paths(A, n)
    alpha = count_disjoint_cycles(A, n)
    I3 = I_at_x(alpha, 3)
    max_k = max(alpha.keys()) if alpha else 0
    ratio = I3/H if H > 0 else 0
    expected = (3/2)**max_k if max_k > 0 else 1
    print(f"  H={H:3d}, I(3)={I3:4d}, ratio={ratio:.4f}, max_k={max_k}, (3/2)^max_k={expected:.4f}")

# =====================================================================
print(f"\n{'='*70}")
print("SECTION 8: PARITY OF I(3)")
print("=" * 70)

n = 7
tb = n*(n-1)//2
np.random.seed(42)

i3_mod2 = defaultdict(int)
i3_mod3 = defaultdict(int)
i3_mod6 = defaultdict(int)

for trial in range(2000):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    alpha = count_disjoint_cycles(A, n)
    I3 = I_at_x(alpha, 3)

    i3_mod2[I3 % 2] += 1
    i3_mod3[I3 % 3] += 1
    i3_mod6[I3 % 6] += 1

print(f"\nn=7, 2000 samples:")
print(f"  I(3) mod 2: {dict(sorted(i3_mod2.items()))}")
print(f"  I(3) mod 3: {dict(sorted(i3_mod3.items()))}")
print(f"  I(3) mod 6: {dict(sorted(i3_mod6.items()))}")

# =====================================================================
print(f"\n{'='*70}")
print("SECTION 9: I(2) AND I(3) TOGETHER — INFORMATION CONTENT")
print("=" * 70)

# From OCF: H = 1 + 2α₁ + 4α₂ + 8α₃ + ...
# From I(3): I₃ = 1 + 3α₁ + 9α₂ + 27α₃ + ...
# System: [H-1, I₃-1] = [[2,4,8,...],[3,9,27,...]] [α₁,α₂,α₃,...]
# With 2 equations, we can solve for α₁ and bound higher terms!
# α₁ = (9(H-1) - 4(I₃-1)) / (9·2 - 4·3) = (9H - 4I₃ - 5) / 6
# α₂ = (3(H-1) - 2(I₃-1)) / (3·4 - 2·9) ... wait, need to check

# If only α₁ and α₂ are nonzero:
# H = 1 + 2α₁ + 4α₂
# I₃ = 1 + 3α₁ + 9α₂
# α₂ = (3(H-1) - 2(I₃-1))/(3·4 - 2·9) = (3H-2I₃+2-3)/(12-18) = (3H-2I₃-1)/(-6)
# α₂ = (2I₃ - 3H + 1)/6
# α₁ = (H - 1 - 4α₂)/2

print("\nIf only α₁, α₂ nonzero (valid for n ≤ 8):")
print("  α₂ = (2·I(3) - 3·H + 1) / 6")
print("  α₁ = (H - 1 - 4·α₂) / 2")
print()

n = 5
for cf in sorted(classes.keys()):
    A = bits_to_adj(classes[cf][0], n)
    H = count_ham_paths(A, n)
    alpha = count_disjoint_cycles(A, n)
    I3 = I_at_x(alpha, 3)

    a2_formula = (2*I3 - 3*H + 1) / 6
    a1_formula = (H - 1 - 4*a2_formula) / 2

    a1_actual = alpha.get(1, 0)
    a2_actual = alpha.get(2, 0)

    match = "✓" if abs(a1_formula - a1_actual) < 0.01 and abs(a2_formula - a2_actual) < 0.01 else "✗"
    print(f"  H={H:3d}, I3={I3:4d}: α₁={a1_actual} (formula:{a1_formula:.1f}), α₂={a2_actual} (formula:{a2_formula:.1f}) {match}")

print(f"\nKey insight: knowing BOTH I(2) and I(3) lets us EXTRACT α₁ and α₂!")
print(f"This is WHY 2 and 3 are the 'keys to the universe' —")
print(f"they're the minimal pair of evaluation points that resolve the")
print(f"first two levels of the odd-cycle hierarchy!")

# =====================================================================
print(f"\n{'='*70}")
print("SECTION 10: GENERAL PRINCIPLE — VANDERMONDE INVERSION")
print("=" * 70)

print("""
Given evaluations at x = 2, 3, 4, ..., m+1:
  I(k) = Σ_{j=0}^{m} α_j k^j    for k = 2, 3, ..., m+1

This is a VANDERMONDE SYSTEM:
  [I(2)]     [2^0  2^1  2^2  ...  2^m ] [α_0]
  [I(3)]   = [3^0  3^1  3^2  ...  3^m ] [α_1]
  [ ⋮  ]     [ ⋮    ⋮    ⋮   ...   ⋮  ] [ ⋮ ]
  [I(m+1)]   [(m+1)^0 ... (m+1)^m     ] [α_m]

Vandermonde matrix with nodes 2,3,...,m+1 is INVERTIBLE (all distinct).
So m evaluations at x=2,3,...,m+1 determine ALL α_0,...,α_{m-1} EXACTLY!

MINIMUM INFORMATION:
  x=2 alone: determines H but not α_k individually
  x=2,3: determines α₁ and α₂ (resolves 3-cycle and 5-cycle counts)
  x=2,3,4: determines α₁, α₂, α₃ (adds 7-cycle info)

The PAIR (2,3) is the SMALLEST pair of positive integers ≥ 2 that
gives nontrivial information. And for n ≤ 8, α₃=0 so (2,3) suffices!

BEAUTIFUL: 2 and 3 are the first two primes, and they resolve
the first two levels of the tournament hierarchy.
""")

# Verify Vandermonde inversion for n=7
n = 7
tb = n*(n-1)//2
np.random.seed(42)
print(f"n=7 Vandermonde verification (100 samples):")
ok = 0
bad = 0
for trial in range(100):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    alpha = count_disjoint_cycles(A, n)

    H = I_at_x(alpha, 2)
    I3 = I_at_x(alpha, 3)

    # Recover α₁, α₂ from (H, I3) assuming α₃ = 0
    a2_rec = (2*I3 - 3*H + 1) / 6
    a1_rec = (H - 1 - 4*a2_rec) / 2

    a1_actual = alpha.get(1, 0)
    a2_actual = alpha.get(2, 0)
    a3_actual = alpha.get(3, 0)

    if a3_actual == 0 and abs(a1_rec - a1_actual) < 0.01 and abs(a2_rec - a2_actual) < 0.01:
        ok += 1
    else:
        bad += 1
        if bad <= 5:
            print(f"  MISMATCH: α=({a1_actual},{a2_actual},{a3_actual}), recovered=({a1_rec:.1f},{a2_rec:.1f})")

print(f"  OK: {ok}, mismatches: {bad}")

print("\nDone.")
