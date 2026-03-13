#!/usr/bin/env python3
"""
jacobsthal_23_deep.py — opus-2026-03-14-S71d

DEEP EXPLORATION: The 2-3 Jacobsthal duality

KEY INSIGHT from generalized_tournament_x.py:
  J_x(n) = (r^n - s^n)/(r-s) where r = (1+sqrt(1+4x))/2

  x=2: r=2, s=-1, J_2(n) = (2^n - (-1)^n)/3
  x=6: r=3, s=-2, J_6(n) = (3^n - (-2)^n)/5

  Root k requires x = k(k-1):
    k=2 → x=2: OCF lives here! H = I(CG, 2)
    k=3 → x=6: I(CG, 6) is the "3-spectrum"
    k=4 → x=12, k=5 → x=20, etc.

QUESTIONS:
1. What does I(CG, 6) count for tournaments?
2. Is there an OCF-like identity at x=6?
3. How does the (2,3) → (6,5) shift relate to (7,8)?
4. The denominators 3,5 are the odd cycle lengths — coincidence?
5. k-nacci gap ratios: does gap_3/gap_2 approach 3/2?
"""

import sys, time
import numpy as np
from itertools import combinations
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

def independence_poly(A, n, max_x=10):
    """Compute I(CG, x) for the conflict graph of tournament A."""
    # Find all odd cycles (3,5,7-cycles)
    cycles = []

    # 3-cycles
    for triple in combinations(range(n), 3):
        a, b, c = triple
        if (A[a][b] and A[b][c] and A[c][a]) or (A[a][c] and A[c][b] and A[b][a]):
            cycles.append(frozenset(triple))

    # 5-cycles (vertex-disjoint check needs actual cycles)
    if n >= 5:
        for combo in combinations(range(n), 5):
            from itertools import permutations
            for perm in permutations(combo):
                if all(A[perm[i]][perm[(i+1)%5]] for i in range(5)):
                    cycles.append(frozenset(combo))
                    break

    # 7-cycles
    if n >= 7:
        for combo in combinations(range(n), 7):
            found = False
            # Check a subset of permutations (exhaustive is too slow)
            verts = list(combo)
            for _ in range(100):
                np.random.shuffle(verts)
                if all(A[verts[i]][verts[(i+1)%7]] for i in range(7)):
                    cycles.append(frozenset(combo))
                    found = True
                    break

    cycles = list(set(cycles))
    nc = len(cycles)

    # Build conflict graph
    adj = [[False]*nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if cycles[i] & cycles[j]:
                adj[i][j] = True
                adj[j][i] = True

    # Independence polynomial by inclusion-exclusion / direct enumeration
    alpha = [0] * (nc + 1)
    alpha[0] = 1

    # Enumerate all independent sets
    for mask in range(1 << nc):
        bits_set = []
        for i in range(nc):
            if mask & (1 << i):
                bits_set.append(i)
        # Check independence
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

    # Evaluate at various x
    results = {}
    for x in range(max_x + 1):
        val = sum(alpha[k] * x**k for k in range(nc + 1))
        results[x] = val

    return results, alpha

def pfaffian_squared(A, n):
    """Compute det(I + 2A) = Pf(S)^2 where S = 2A - I."""
    M = np.eye(n) + 2 * A.astype(float)
    return int(round(np.linalg.det(M)))

print("=" * 70)
print("THE 2-3 JACOBSTHAL DUALITY: DEEP STRUCTURE")
print("=" * 70)

# Part 1: I(CG, x) at x=2,3,6 for all n=5 tournaments
print("\n" + "=" * 70)
print("PART 1: I(CG, x) spectrum at n=5 (all 1024 tournaments)")
print("=" * 70)

n = 5
tb = n*(n-1)//2
x_vals = [1, 2, 3, 4, 6, 8, 12]

spectrum = defaultdict(list)
t0 = time.time()

for bits in range(1 << tb):
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)
    ix, alpha = independence_poly(A, n, max_x=12)
    key = tuple(alpha[k] for k in range(len(alpha)) if alpha[k] > 0 or k <= 3)
    spectrum[key].append({
        'H': H, 'bits': bits,
        **{f'I({x})': ix[x] for x in x_vals}
    })

dt = time.time() - t0
print(f"  Computed in {dt:.1f}s")

# Show the spectrum
print(f"\n  {'alpha':>20s} | {'I(1)':>5s} {'I(2)':>5s} {'I(3)':>5s} {'I(6)':>5s} | {'H':>5s} | Pf²  | H²-Pf²")
print(f"  {'-'*20}-+-{'-'*23}-+-{'-'*5}-+------+--------")

for key in sorted(spectrum.keys()):
    group = spectrum[key]
    ex = group[0]
    H = ex['H']
    A = bits_to_adj(ex['bits'], n)
    pf2 = pfaffian_squared(A, n)
    gap = H*H - pf2
    alpha_str = str(list(key))
    print(f"  {alpha_str:>20s} | {ex['I(1)']:5d} {ex['I(2)']:5d} {ex['I(3)']:5d} {ex['I(6)']:5d} | {H:5d} | {pf2:5d} | {gap:7d}   (×{len(group)})")

# Part 2: The k(k-1) law — integer-root evaluations
print("\n" + "=" * 70)
print("PART 2: INTEGER ROOT LAW — I(CG, k(k-1)) spectrum")
print("=" * 70)

print("\n  k=2: x=2, I(CG,2) = H (OCF)")
print("  k=3: x=6, I(CG,6) = '6-spectrum'")
print("  k=4: x=12, I(CG,12) = '12-spectrum'")

# Check: does I(CG,6) have nice properties?
print("\n  I(CG,6) properties at n=5:")
i6_vals = sorted(set(group[0]['I(6)'] for group in spectrum.values()))
print(f"    Distinct values: {i6_vals}")
print(f"    All ≡ 1 (mod 6)? {all(v % 6 == 1 for v in i6_vals)}")
print(f"    All ≡ 1 (mod 3)? {all(v % 3 == 1 for v in i6_vals)}")
print(f"    All odd? {all(v % 2 == 1 for v in i6_vals)}")

# Part 3: Vandermonde at x=2,6 vs x=2,3
print("\n" + "=" * 70)
print("PART 3: VANDERMONDE EXTRACTION — (x=2,6) vs (x=2,3)")
print("=" * 70)

print("\n  I(CG, x) = 1 + α₁x + α₂x² + α₃x³ + ...")
print("  At n=5: max independent set size is 2 (or maybe 3)")
print("  So I = 1 + α₁x + α₂x²")
print()

# For each alpha pattern, check Vandermonde extraction
for key in sorted(spectrum.keys()):
    group = spectrum[key]
    ex = group[0]
    a0 = key[0]  # should be 1
    a1 = key[1] if len(key) > 1 else 0
    a2 = key[2] if len(key) > 2 else 0

    # Verify: I(2) = 1 + 2a1 + 4a2
    i2_check = 1 + 2*a1 + 4*a2
    # I(6) = 1 + 6a1 + 36a2
    i6_check = 1 + 6*a1 + 36*a2
    # I(3) = 1 + 3a1 + 9a2
    i3_check = 1 + 3*a1 + 9*a2

    print(f"  α=({a1},{a2}): I(2)={ex['I(2)']}={i2_check}, I(3)={ex['I(3)']}={i3_check}, I(6)={ex['I(6)']}={i6_check}")

# Part 4: The denominator pattern 3,5,7,... — odd cycles!
print("\n" + "=" * 70)
print("PART 4: JACOBSTHAL DENOMINATORS = ODD NUMBERS → CYCLE LENGTHS")
print("=" * 70)

print("""
  J_x(n) = (r^n - s^n)/(r-s) where r-s = √(1+4x)

  x = k(k-1):  r = k,  s = 1-k,  r-s = 2k-1

  k=1: x=0,  r=1, s=0,  r-s=1   → J_0(n) = 1 for all n
  k=2: x=2,  r=2, s=-1, r-s=3   → J_2(n) = (2^n-(-1)^n)/3
  k=3: x=6,  r=3, s=-2, r-s=5   → J_6(n) = (3^n-(-2)^n)/5
  k=4: x=12, r=4, s=-3, r-s=7   → J_12(n) = (4^n-(-3)^n)/7
  k=5: x=20, r=5, s=-4, r-s=9   → J_20(n) = (5^n-(-4)^n)/9

  Denominators: 1, 3, 5, 7, 9, ...  ← ODD NUMBERS!
  But cycle lengths in tournaments: 3, 5, 7, ...

  THE CORRESPONDENCE:
    k=2 (OCF level): denominator 3 ↔ 3-cycles
    k=3 (6-spectrum): denominator 5 ↔ 5-cycles
    k=4 (12-spectrum): denominator 7 ↔ 7-cycles

  Is this a coincidence, or is there a structural reason?
""")

# Part 5: I(CG,6) at n=7 — does it "see" 5-cycles more than I(CG,2)?
print("=" * 70)
print("PART 5: I(CG,x) SENSITIVITY TO CYCLE STRUCTURE AT n=7")
print("=" * 70)

n = 7
tb = n*(n-1)//2
np.random.seed(42)

print(f"\nn=7: 500 random tournaments — correlations between I(x) and cycle counts")
data = []
t0 = time.time()

for trial in range(500):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)

    H = count_ham_paths(A, n)
    c3 = int(np.trace(A @ A @ A)) // 3
    A5 = np.linalg.matrix_power(A, 5)
    c5_trace = int(np.trace(A5)) // 5  # This overcounts at n=7 (but only for 5-cycles sharing vertices)

    # Independence polynomial coefficients
    # At n=7, need cycles
    cycles_3 = []
    for triple in combinations(range(n), 3):
        a, b, c = triple
        if (A[a][b] and A[b][c] and A[c][a]) or (A[a][c] and A[c][b] and A[b][a]):
            cycles_3.append(frozenset(triple))

    cycles_5 = []
    for combo in combinations(range(n), 5):
        verts = list(combo)
        for perm in __import__('itertools').permutations(verts):
            if all(A[perm[i]][perm[(i+1)%5]] for i in range(5)):
                cycles_5.append(frozenset(combo))
                break

    all_cycles = list(set(cycles_3 + cycles_5))
    nc = len(all_cycles)

    # Build CG adjacency
    cg_adj = [[False]*nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if all_cycles[i] & all_cycles[j]:
                cg_adj[i][j] = True
                cg_adj[j][i] = True

    # Count independent sets by size (up to size 3)
    alpha = [0] * 4
    alpha[0] = 1
    for i in range(nc):
        alpha[1] += 1
    for i in range(nc):
        for j in range(i+1, nc):
            if not cg_adj[i][j]:
                alpha[2] += 1
    for i in range(nc):
        for j in range(i+1, nc):
            if cg_adj[i][j]: continue
            for k in range(j+1, nc):
                if not cg_adj[i][k] and not cg_adj[j][k]:
                    alpha[3] += 1

    I2 = sum(alpha[k] * 2**k for k in range(4))
    I3 = sum(alpha[k] * 3**k for k in range(4))
    I6 = sum(alpha[k] * 6**k for k in range(4))

    data.append({
        'H': H, 'I2': I2, 'I3': I3, 'I6': I6,
        'c3': len(cycles_3), 'c5': len(cycles_5),
        'a1': alpha[1], 'a2': alpha[2], 'a3': alpha[3]
    })

    if trial % 100 == 0:
        dt = time.time() - t0
        print(f"  trial {trial}: {dt:.1f}s")

dt = time.time() - t0
print(f"  Done: {dt:.1f}s")

# Verify I(2) = H
mismatches_h = sum(1 for d in data if d['I2'] != d['H'])
print(f"\n  I(CG,2) = H? Mismatches: {mismatches_h}/500")

# Correlations
import numpy as np
H_arr = np.array([d['H'] for d in data])
I3_arr = np.array([d['I3'] for d in data])
I6_arr = np.array([d['I6'] for d in data])
c3_arr = np.array([d['c3'] for d in data])
c5_arr = np.array([d['c5'] for d in data])
a1_arr = np.array([d['a1'] for d in data])
a2_arr = np.array([d['a2'] for d in data])

def corr(x, y):
    if np.std(x) == 0 or np.std(y) == 0:
        return 0
    return np.corrcoef(x, y)[0, 1]

print(f"\n  Correlations:")
print(f"    corr(H, c3)  = {corr(H_arr, c3_arr):.4f}")
print(f"    corr(H, c5)  = {corr(H_arr, c5_arr):.4f}")
print(f"    corr(H, α₁)  = {corr(H_arr, a1_arr):.4f}")
print(f"    corr(H, α₂)  = {corr(H_arr, a2_arr):.4f}")
print(f"    corr(I6, c3) = {corr(I6_arr, c3_arr):.4f}")
print(f"    corr(I6, c5) = {corr(I6_arr, c5_arr):.4f}")
print(f"    corr(I6, α₁) = {corr(I6_arr, a1_arr):.4f}")
print(f"    corr(I6, α₂) = {corr(I6_arr, a2_arr):.4f}")

# KEY: does I(6) - I(2) isolate α₂?
# I(6) - I(2) = (6-2)α₁ + (36-4)α₂ + (216-8)α₃ = 4α₁ + 32α₂ + 208α₃
diff_62 = I6_arr - H_arr
print(f"\n  I(6) - I(2) = 4α₁ + 32α₂ + 208α₃")
print(f"    corr(I6-I2, α₁) = {corr(diff_62, a1_arr):.4f}")
print(f"    corr(I6-I2, α₂) = {corr(diff_62, a2_arr):.4f}")

# Part 6: The 3/2 ratio in the Jacobsthal framework
print("\n" + "=" * 70)
print("PART 6: THE 3/2 RATIO — JACOBSTHAL PERSPECTIVE")
print("=" * 70)

print("""
  J_2(n) = (2^n - (-1)^n) / 3
  J_6(n) = (3^n - (-2)^n) / 5

  Ratio J_6(n)/J_2(n) → (3/2)^n · 3/5 as n → ∞

  At the OCF level:
    H = I(CG, 2) = 1 + 2α₁ + 4α₂ + 8α₃
    I₃ = I(CG, 3) = 1 + 3α₁ + 9α₂ + 27α₃
    I₆ = I(CG, 6) = 1 + 6α₁ + 36α₂ + 216α₃

  Ratio I₃/H within lambda fibers ≈ 3/2 (proved!)

  What about I₆/H ratio?
""")

# I6/H ratio
ratios = I6_arr / H_arr
print(f"  I₆/H ratio: mean={np.mean(ratios):.4f}, std={np.std(ratios):.4f}")
print(f"  I₃/H ratio: mean={np.mean(I3_arr/H_arr):.4f}, std={np.std(I3_arr/H_arr):.4f}")

# Part 7: Deeper k-nacci gap analysis
print("\n" + "=" * 70)
print("PART 7: k-NACCI GAPS — PRECISE CONVERGENCE RATES")
print("=" * 70)

def knacci_root(k):
    """Root of x^k = x^{k-1} + x^{k-2} + ... + 1"""
    from numpy.polynomial import polynomial as P
    # x^k - x^{k-1} - ... - 1 = 0
    coeffs = [-1] * k + [1]  # constant term first in numpy
    roots = np.roots(coeffs[::-1])
    real_pos = [r.real for r in roots if abs(r.imag) < 1e-10 and r.real > 1]
    return max(real_pos) if real_pos else None

def weighted_knacci_root(k):
    """Root of x^k = 1·x^{k-1} + 2·x^{k-2} + 4·x^{k-3} + ... + 2^{k-1}"""
    # x^k - sum_{j=0}^{k-1} 2^j x^{k-1-j} = 0
    coeffs = [1]
    for j in range(k):
        coeffs.append(-2**j)
    roots = np.roots(coeffs)
    real_pos = [r.real for r in roots if abs(r.imag) < 1e-10 and r.real > 1]
    return max(real_pos) if real_pos else None

print(f"  {'k':>3s}  {'φ_k':>12s}  {'gap_2':>12s}  {'ratio':>8s}  {'ψ_k':>12s}  {'gap_3':>12s}  {'ratio':>8s}  {'gap3/gap2':>10s}")
print(f"  {'─'*3}  {'─'*12}  {'─'*12}  {'─'*8}  {'─'*12}  {'─'*12}  {'─'*8}  {'─'*10}")

prev_g2, prev_g3 = None, None
for k in range(2, 25):
    phi_k = knacci_root(k)
    psi_k = weighted_knacci_root(k)
    if phi_k is None or psi_k is None:
        continue
    g2 = phi_k - 2
    g3 = psi_k - 3
    r2 = g2 / prev_g2 if prev_g2 else float('nan')
    r3 = g3 / prev_g3 if prev_g3 else float('nan')
    ratio_32 = g3/g2 if g2 > 1e-15 else float('nan')

    mark = ""
    if k == 7:
        mark = " ← k=7"
    elif k == 8:
        mark = " ← k=8"

    print(f"  {k:3d}  {phi_k:12.8f}  {g2:12.8e}  {r2:8.4f}  {psi_k:12.8f}  {g3:12.8e}  {r3:8.4f}  {ratio_32:10.4f}{mark}")
    prev_g2, prev_g3 = g2, g3

# Part 8: The MAGICAL n=7,8 interplay
print("\n" + "=" * 70)
print("PART 8: WHY 7 AND 8 — THE STRUCTURAL TRANSITION")
print("=" * 70)

print("""
  The user's insight: "if one knows 2 and 3, then we have the keys"

  Tournament reading:
    2 = I(CG, 2) root → H (Hamiltonian path count)
    3 = I(CG, 6) root → "6-spectrum" (amplified cycle information)

  But also:
    2 = dominant root of Jacobsthal recurrence (x=2)
    3 = dominant root of x=6 Jacobsthal
    3/2 = ratio I₃_res/H_res within lambda fibers (proved)

  The 7→8 transition:
    n=7: α₂ is purely (3,3) type → J_2 level only
    n=8: α₂ gets (3,5) type → J_2 × J_6 cross-coupling!

  So n=8 is where the TWO Jacobsthal levels first INTERACT.
  The "keys to the universe" (2 and 3) become coupled at n=2³=8.

  n=9: α₃ appears (3+3+3=9) → cubic independence terms
    This is where I(CG,x) becomes degree ≥ 3 in some cases
    The third Jacobsthal level (x=12, root=4) starts to matter
""")

# Check: at n=7, does knowing I(2) and I(6) together determine everything?
print("DETERMINATION POWER at n=7:")
print("  (H alone) vs (H, I₃) vs (H, I₆) vs (H, I₃, I₆)")
print()

# Group by (H), (H,I3), (H,I6), (H,I3,I6)
groups_H = defaultdict(set)
groups_HI3 = defaultdict(set)
groups_HI6 = defaultdict(set)
groups_all = defaultdict(set)

for d in data:
    c = (d['c3'], d['c5'])
    groups_H[d['H']].add(c)
    groups_HI3[(d['H'], d['I3'])].add(c)
    groups_HI6[(d['H'], d['I6'])].add(c)
    groups_all[(d['H'], d['I3'], d['I6'])].add(c)

amb_H = sum(1 for g in groups_H.values() if len(g) > 1)
amb_HI3 = sum(1 for g in groups_HI3.values() if len(g) > 1)
amb_HI6 = sum(1 for g in groups_HI6.values() if len(g) > 1)
amb_all = sum(1 for g in groups_all.values() if len(g) > 1)

print(f"  H alone: {amb_H}/{len(groups_H)} ambiguous for (c3,c5)")
print(f"  (H, I₃): {amb_HI3}/{len(groups_HI3)} ambiguous")
print(f"  (H, I₆): {amb_HI6}/{len(groups_HI6)} ambiguous")
print(f"  (H, I₃, I₆): {amb_all}/{len(groups_all)} ambiguous")

print("\nDone.")
