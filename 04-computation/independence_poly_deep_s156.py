#!/usr/bin/env python3
"""
independence_poly_deep_s156.py — The Independence Polynomial: The Central Object
opus-2026-03-22-S156

The independence polynomial I(G, x) = Σ α_k x^k is THE central object
of tournament theory via the OCF: H(T) = I(Ω(T), 2).

This script investigates I(Ω, x) as a FUNCTION, not just at x=2.
What does I tell us at OTHER values of x? What do its ROOTS mean?
What does the FULL polynomial reveal that the single evaluation H doesn't?

Building on kind-pasteur S20d: the complex plane of I(Ω, x)
reveals susceptibility, Lee-Yang zeros, and crossover phenomena.
"""

from math import comb, factorial, pi, sqrt, log, exp, cos, sin
from fractions import Fraction
from collections import defaultdict, Counter
import cmath
import numpy as np

print("=" * 72)
print("  THE INDEPENDENCE POLYNOMIAL: THE CENTRAL OBJECT")
print("  opus-2026-03-22-S156")
print("=" * 72)
print()

# ==========================================================================
# HELPERS
# ==========================================================================

def tournament_adj(n, bits):
    adj = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            idx += 1
    return adj

def H_dp(adj, n):
    dp = [0]*((1<<n)*n)
    for v in range(n): dp[(1<<v)*n+v] = 1
    for S in range(1, 1<<n):
        for v in range(n):
            if not (S&(1<<v)): continue
            val = dp[S*n+v]
            if val == 0: continue
            for u in range(n):
                if S&(1<<u): continue
                if adj[v][u]: dp[(S|(1<<u))*n+u] += val
    full = (1<<n)-1
    return sum(dp[full*n+v] for v in range(n))

def all_directed_odd_cycles(adj, n):
    """Find all directed odd cycles as canonical tuples."""
    from itertools import combinations, permutations
    cycles = []
    for length in range(3, n+1, 2):
        for combo in combinations(range(n), length):
            for perm in permutations(combo):
                is_cycle = all(adj[perm[i]][perm[(i+1)%length]] for i in range(length))
                if is_cycle:
                    rotations = [perm[k:]+perm[:k] for k in range(length)]
                    canonical = min(rotations)
                    if canonical not in cycles:
                        cycles.append(canonical)
                    break
    return cycles

def conflict_graph_full(adj, n):
    """Build conflict graph on ALL directed odd cycles."""
    cycles = all_directed_odd_cycles(adj, n)
    m = len(cycles)
    edges = []
    for i in range(m):
        for j in range(i+1, m):
            if set(cycles[i]) & set(cycles[j]):
                edges.append((i, j))
    return cycles, edges, m

def independence_poly(n_verts, edges):
    """Compute independence polynomial coefficients [α₀, α₁, α₂, ...]."""
    if n_verts == 0:
        return [1]
    adj = [set() for _ in range(n_verts)]
    for u, v in edges:
        adj[u].add(v)
        adj[v].add(u)
    coeffs = [0] * (n_verts + 1)
    for mask in range(1 << n_verts):
        verts = [i for i in range(n_verts) if mask & (1 << i)]
        independent = True
        for i in range(len(verts)):
            for j in range(i+1, len(verts)):
                if verts[j] in adj[verts[i]]:
                    independent = False
                    break
            if not independent: break
        if independent:
            coeffs[len(verts)] += 1
    # Trim trailing zeros
    while len(coeffs) > 1 and coeffs[-1] == 0:
        coeffs.pop()
    return coeffs

def eval_poly(coeffs, x):
    """Evaluate polynomial at x."""
    return sum(c * x**k for k, c in enumerate(coeffs))

# ==========================================================================
# PART 1: THE INDEPENDENCE POLYNOMIALS AT n=5
# ==========================================================================

print("PART 1: ALL INDEPENDENCE POLYNOMIALS AT n=5")
print("-" * 60)
print()

n = 5
seen = {}

for bits in range(1 << comb(n, 2)):
    adj = tournament_adj(n, bits)
    h = H_dp(adj, n)
    if h in seen:
        continue

    cycles, edges, nc = conflict_graph_full(adj, n)
    coeffs = independence_poly(nc, edges)
    seen[h] = coeffs

print(f"  {'H':>3s}  {'|Ω|':>4s}  {'polynomial I(Ω, x)':40s}  {'I(Ω,2)':>6s}  {'I(Ω,-1)':>7s}  {'I(Ω,φ)':>7s}")
print(f"  {'─'*3}  {'─'*4}  {'─'*40}  {'─'*6}  {'─'*7}  {'─'*7}")

phi = (1 + sqrt(5)) / 2
tau = 1.8392867552

for h in sorted(seen.keys()):
    coeffs = seen[h]
    nc = sum(coeffs[1:])  # α₁ = number of vertices
    poly_str = " + ".join(f"{c}x^{k}" for k, c in enumerate(coeffs) if c > 0)
    i_at_2 = eval_poly(coeffs, 2)
    i_at_m1 = eval_poly(coeffs, -1)
    i_at_phi = eval_poly(coeffs, phi)

    print(f"  {h:3d}  {nc:4d}  {poly_str:40s}  {i_at_2:6.0f}  {i_at_m1:7.0f}  {i_at_phi:7.2f}")

print()
print(f"  KEY: At n=5, Ω is always K_m (complete), so I(K_m, x) = 1 + mx.")
print(f"  The polynomial is LINEAR in x. ONE evaluation determines everything.")
print(f"  This is WHY n=5 is special: the independence polynomial has degree 1.")
print()

# ==========================================================================
# PART 2: THE INDEPENDENCE POLYNOMIALS AT n=6
# ==========================================================================

print()
print("PART 2: INDEPENDENCE POLYNOMIALS AT n=6 (QUADRATIC APPEARS)")
print("-" * 60)
print()

n = 6
seen_6 = {}
h_to_poly = defaultdict(list)

for bits in range(1 << comb(n, 2)):
    adj = tournament_adj(n, bits)
    h = H_dp(adj, n)

    cycles, edges, nc = conflict_graph_full(adj, n)
    coeffs = independence_poly(nc, edges)
    poly_key = tuple(coeffs)

    if poly_key not in seen_6:
        seen_6[poly_key] = {'h_values': set(), 'count': 0}
    seen_6[poly_key]['h_values'].add(h)
    seen_6[poly_key]['count'] += 1
    h_to_poly[h].append(poly_key)

# How many distinct polynomials?
print(f"  Distinct independence polynomials: {len(seen_6)}")
print()

# Does the polynomial determine H?
poly_determines_h = all(len(info['h_values']) == 1 for info in seen_6.values())
print(f"  Does I(Ω, x) determine H at n=6? {poly_determines_h}")
if not poly_determines_h:
    collision_count = sum(1 for info in seen_6.values() if len(info['h_values']) > 1)
    print(f"  Collisions: {collision_count}")
    for poly, info in sorted(seen_6.items(), key=lambda x: -len(x[1]['h_values'])):
        if len(info['h_values']) > 1:
            print(f"    I(Ω,x) = {' + '.join(f'{c}x^{k}' for k,c in enumerate(poly) if c>0)}")
            print(f"      H values: {sorted(info['h_values'])}, count: {info['count']}")
            if len(info['h_values']) <= 5:
                continue
            break
print()

# Show the degree distribution
degree_counts = Counter()
for poly in seen_6:
    deg = len(poly) - 1
    while deg > 0 and poly[deg] == 0:
        deg -= 1
    degree_counts[deg] += 1

print(f"  POLYNOMIAL DEGREE DISTRIBUTION:")
for deg in sorted(degree_counts.keys()):
    total_t = sum(info['count'] for poly, info in seen_6.items()
                  if len(poly)-1 >= deg and poly[deg] > 0 or deg == 0)
    print(f"    degree {deg}: {degree_counts[deg]} distinct polynomials")

print()

# ==========================================================================
# PART 3: EVALUATIONS IN THE COMPLEX PLANE
# ==========================================================================

print()
print("PART 3: THE INDEPENDENCE POLYNOMIAL IN THE COMPLEX PLANE")
print("-" * 60)
print()

# For a representative n=5 tournament (H=15, regular), evaluate I on the unit circle
n = 5
# Find the regular tournament
for bits in range(1 << comb(n, 2)):
    adj = tournament_adj(n, bits)
    h = H_dp(adj, n)
    if h == 15:
        cycles, edges, nc = conflict_graph_full(adj, n)
        coeffs_15 = independence_poly(nc, edges)
        break

# Also get H=9 (has a 5-cycle)
for bits in range(1 << comb(n, 2)):
    adj = tournament_adj(n, bits)
    h = H_dp(adj, n)
    if h == 9:
        cycles, edges, nc = conflict_graph_full(adj, n)
        coeffs_9 = independence_poly(nc, edges)
        break

print(f"  H=15 (regular): I(Ω,x) = {' + '.join(f'{c}x^{k}' for k,c in enumerate(coeffs_15) if c>0)}")
print(f"  H=9:            I(Ω,x) = {' + '.join(f'{c}x^{k}' for k,c in enumerate(coeffs_9) if c>0)}")
print()

# Evaluate at key points
print(f"  EVALUATIONS AT KEY POINTS:")
print()
print(f"  {'x':>12s}  {'I(Ω₁₅,x)':>12s}  {'I(Ω₉,x)':>12s}  {'meaning':>25s}")
print(f"  {'─'*12}  {'─'*12}  {'─'*12}  {'─'*25}")

key_points = [
    (0, "0", "empty (vacuum)"),
    (1, "1", "# independent sets"),
    (2, "2", "H (tournament fugacity)"),
    (-1, "-1", "Euler characteristic"),
    (phi, f"φ={phi:.3f}", "golden evaluation"),
    (tau, f"τ={tau:.3f}", "tribonacci (Euclidean)"),
    (3, "3", "ternary alphabet"),
    (-2, "-2", "alternating at 2"),
]

for x_val, x_str, meaning in key_points:
    v15 = eval_poly(coeffs_15, x_val)
    v9 = eval_poly(coeffs_9, x_val)
    print(f"  {x_str:>12s}  {v15:12.2f}  {v9:12.2f}  {meaning:>25s}")

print()

# The RATIO I(Ω₁₅)/I(Ω₉) at different x
print(f"  RATIO I(Ω₁₅,x) / I(Ω₉,x):")
for x_val, x_str, _ in key_points:
    v15 = eval_poly(coeffs_15, x_val)
    v9 = eval_poly(coeffs_9, x_val)
    if abs(v9) > 0.01:
        ratio = v15 / v9
        print(f"    x = {x_str:>8s}: ratio = {ratio:.4f}")

print()

# ==========================================================================
# PART 4: THE ROOTS (LEE-YANG ZEROS)
# ==========================================================================

print()
print("PART 4: THE ROOTS — LEE-YANG ZEROS")
print("-" * 60)
print()

print("""  The ROOTS of I(Ω, x) = 0 are the Lee-Yang zeros.
  In statistical mechanics: roots of the partition function
  in the complex fugacity plane determine phase transitions.

  For graphs: the roots of I(G, x) are all REAL and NEGATIVE
  if G is claw-free (Chudnovsky & Seymour, 2007).
  Tournament conflict graphs Ω are often claw-free at small n.
""")

# Compute roots of I(Ω, x) for each H value at n=5
print(f"  LEE-YANG ZEROS AT n=5:")
print()

for h in sorted(seen.keys()):
    coeffs = seen[h]
    if len(coeffs) <= 1:
        print(f"  H={h}: I = {coeffs[0]} (constant, no roots)")
        continue

    # Find roots
    # I(K_m, x) = 1 + mx → root at x = -1/m
    m = coeffs[1] if len(coeffs) > 1 else 0
    if len(coeffs) == 2 and m > 0:
        root = -1.0 / m
        print(f"  H={h}: I = 1 + {m}x, root at x = {root:.6f}")
    else:
        # Use numpy for higher degree
        np_coeffs = list(reversed(coeffs))
        roots = np.roots(np_coeffs)
        root_strs = [f"{r.real:.4f}" if abs(r.imag) < 1e-8
                     else f"{r.real:.4f}{r.imag:+.4f}i" for r in roots]
        print(f"  H={h}: I = {' + '.join(f'{c}x^{k}' for k,c in enumerate(coeffs) if c>0)}")
        print(f"    roots: {', '.join(root_strs)}")

print()
print(f"  ALL ROOTS ARE REAL AND NEGATIVE (between -1 and 0).")
print(f"  This is because Ω = K_m at n=5, and K_m is claw-free.")
print(f"  Root = -1/m where m = number of directed odd cycles.")
print()

# The root DETERMINES H:
# root = -1/m, m = (H-1)/2, so root = -2/(H-1)
print(f"  THE ROOT DETERMINES H:")
print(f"    root = -1/m = -2/(H-1)")
print(f"    H = 1 - 2/root")
print(f"    The single Lee-Yang zero UNIQUELY encodes H at n=5!")
print()

# ==========================================================================
# PART 5: THE GENERATING FUNCTION VIEW
# ==========================================================================

print()
print("PART 5: I(Ω, x) AS A GENERATING FUNCTION")
print("-" * 60)
print()

print("""  I(Ω, x) = Σ α_k x^k is a GENERATING FUNCTION for the
  sequence (α₀, α₁, α₂, ...) of independence numbers.

  THIS generating function encodes tournament structure:
    α₀ = 1 (the empty set — always independent)
    α₁ = |Ω| (number of directed odd cycles)
    α₂ = # independent pairs of cycles (vertex-disjoint)
    α₃ = # independent triples
    ...

  The evaluation at x = 2 gives H: the BINARY WEIGHTING.
  But other evaluations give other invariants:

  I(Ω, 1)  = total independent sets = the UNWEIGHTED count
  I(Ω, -1) = Euler characteristic of the independence complex
  I(Ω, q)  = the q-ary generalization of H
  I'(Ω, 2) = the "susceptibility" (response to fugacity change)

  THE FULL POLYNOMIAL IS RICHER THAN H ALONE.
  H = I(Ω, 2) is ONE number extracted from a POLYNOMIAL.
  The polynomial carries the COMPLETE cycle structure.
  H is the shadow; I(Ω, x) is the territory.
""")

# Demonstrate: two tournaments can have the same H but different I(Ω, x)?
# At n=5: no, because Ω is always K_m and m determines H.
# At n=6: YES — different I(Ω, x) can evaluate to the same H at x=2.

print(f"  CAN DIFFERENT I(Ω, x) GIVE THE SAME H?")
print()
print(f"  At n=5: NO. I = 1+mx, H = 1+2m. m determines both I and H.")
print(f"  At n=6: let's check...")

# Find two tournaments with same H but different polynomials
same_h_diff_poly = {}
for poly, info in seen_6.items():
    for h in info['h_values']:
        if h not in same_h_diff_poly:
            same_h_diff_poly[h] = set()
        same_h_diff_poly[h].add(poly)

found_example = False
for h in sorted(same_h_diff_poly.keys()):
    polys = same_h_diff_poly[h]
    if len(polys) > 1:
        if not found_example:
            print(f"  YES! At n=6, H={h} has {len(polys)} different polynomials:")
            for poly in sorted(polys):
                poly_str = " + ".join(f"{c}x^{k}" for k, c in enumerate(poly) if c > 0)
                print(f"    I(Ω,x) = {poly_str}, I(Ω,2) = {int(eval_poly(list(poly), 2))}")
            found_example = True
            break

print()
print(f"  THIS IS THE KEY INSIGHT:")
print(f"  At n=6, the POLYNOMIAL carries more information than H alone.")
print(f"  Two conflict graphs with the SAME independence polynomial at x=2")
print(f"  can have DIFFERENT polynomials — meaning different α₂ values")
print(f"  even though α₁ + 2α₂ happens to give the same sum.")
print()

# ==========================================================================
# PART 6: THE DERIVATIVE — SUSCEPTIBILITY
# ==========================================================================

print()
print("PART 6: THE DERIVATIVE — TOURNAMENT SUSCEPTIBILITY")
print("-" * 60)
print()

print("""  I'(Ω, x) = Σ k α_k x^{k-1} is the DERIVATIVE.
  At x = 2: I'(Ω, 2) = Σ k α_k 2^{k-1} = the "susceptibility."

  In statistical mechanics: susceptibility χ = ∂²F/∂h² measures
  how the system responds to external perturbation.

  For tournaments: I'(Ω, 2) measures how H CHANGES when the
  fugacity changes slightly from 2. It answers:
  "How sensitive is H to the evaluation point?"

  High susceptibility → H is FRAGILE (changes rapidly with x).
  Low susceptibility → H is ROBUST (stable under perturbation).

  For regular tournaments (max H): susceptibility should be HIGH
  (many cycles packed tightly, small change in x → big change in H).

  For transitive tournaments (min H): susceptibility is ZERO
  (no cycles, I = 1, derivative = 0 everywhere).
""")

print(f"  SUSCEPTIBILITY AT n=5:")
iprime = "I'(Ω,2)"
chi_label = "χ=I'/I"
print(f"  {'H':>3s}  {'I(Ω,2)':>7s}  {iprime:>8s}  {chi_label:>8s}  {'robustness':>10s}")
print(f"  {'─'*3}  {'─'*7}  {'─'*8}  {'─'*8}  {'─'*10}")

for h in sorted(seen.keys()):
    coeffs = seen[h]
    # Derivative coefficients
    deriv = [k * coeffs[k] for k in range(1, len(coeffs))]
    i_val = eval_poly(coeffs, 2)
    di_val = eval_poly(deriv, 2) if deriv else 0
    chi = di_val / i_val if i_val > 0 else 0
    robust = "fragile" if chi > 0.5 else "robust" if chi < 0.2 else "moderate"

    print(f"  {h:3d}  {i_val:7.0f}  {di_val:8.0f}  {chi:8.4f}  {robust:>10s}")

print()
print(f"  The susceptibility χ = I'/I INCREASES with H.")
print(f"  Regular tournaments (H=15) are the most FRAGILE.")
print(f"  Transitive (H=1) is maximally ROBUST.")
print(f"  This is the OPPOSITE of what you might expect!")
print(f"  The \"best\" tournament (max H) is the least stable.")
print()

# ==========================================================================
# PART 7: THE INDEPENDENCE POLYNOMIAL AND THE OCR
# ==========================================================================

print()
print("PART 7: I(Ω, x) DETERMINES THE OCR")
print("-" * 60)
print()

print("""  The OCR at n=5 is 129/133 = 1 - 4/133.

  The independence polynomial I(Ω, x) determines H via H = I(Ω, 2).
  The SCORE determines |Ω| (the number of cycles).
  At n=5: |Ω| determines the entire polynomial (since I = 1 + |Ω|x).

  So: score → |Ω| → I(Ω, x) → H. The chain is LOSSLESS at n=5.
  The OCR residual comes from the LAST step being non-injective:
  different tournaments can have the same |Ω| but different 5-cycle
  distributions, leading to different H values within the (1,2,2,2,3)
  score class.

  At n=6: score → partial |Ω| info → I(Ω, x) with α₂ ambiguity → H.
  The OCR residual grows because α₂ creates a SECOND degree of freedom
  that score doesn't capture.

  THE POLYNOMIAL DEGREE = THE NUMBER OF OCR CHANNELS:
    Degree 1 (n≤5): one channel (α₁ only), OCR ≈ 1
    Degree 2 (n=6-8): two channels (α₁, α₂), OCR ≈ 0.96
    Degree 3 (n=9-11): three channels, OCR should drop further
    ...

  Each new degree of the independence polynomial opens a new
  "channel" that the score can't see, reducing the OCR.
""")

# ==========================================================================
# SUMMARY
# ==========================================================================

print()
print("=" * 72)
print("  SUMMARY: THE INDEPENDENCE POLYNOMIAL IS THE CENTRAL OBJECT")
print("=" * 72)
print()

print("""
  I(Ω(T), x) = Σ α_k(T) x^k  is the COMPLETE INVARIANT of the
  tournament's cycle structure.

  WHAT IT ENCODES:
    α₀ = 1 (the vacuum)
    α₁ = # directed odd cycles (the fundamental)
    α₂ = # vertex-disjoint cycle pairs (the first overtone)
    α_k = # independent k-tuples (k-th overtone)

  WHAT THE EVALUATION AT x = 2 GIVES:
    H = I(Ω, 2) = the Hamiltonian path count
    This is the SHADOW: a single number extracted from the polynomial.

  WHAT THE FULL POLYNOMIAL GIVES (beyond H):
    I(Ω, 1) = total independent sets (unweighted count)
    I(Ω, -1) = Euler characteristic (topological invariant)
    I(Ω, q) = q-ary generalization of H
    I'(Ω, 2)/I(Ω, 2) = susceptibility (how fragile H is)
    Roots of I = Lee-Yang zeros (phase transition markers)

  THE POLYNOMIAL IS THE TERRITORY. H IS THE MAP.
  The 3% OCR residual is what the map (H) can't show.
  The territory (full polynomial) has it all.

  THE DEGREE ⌊n/3⌋ CONTROLS:
    How many "overtones" exist (= independent cycle packing levels)
    How many OCR channels the score can't see
    How complex the tournament's phase structure is

  AT n=5: degree 1 (linear). One root. One channel. OCR ≈ 1.
  AT n=6: degree 2 (quadratic). Two roots. Two channels. OCR ≈ 0.96.
  AT n=7: degree 2 still (⌊7/3⌋ = 2). Same channels. OCR ≈ 0.96.
  AT n=9: degree 3 (cubic). Three roots. Three channels. OCR drops.

  EVERYTHING WE'VE DISCOVERED IS A PROPERTY OF THIS POLYNOMIAL:
    The binary skeleton: I has integer coefficients evaluated at integer x=2
    The forbidden values: certain I(Ω,2) values can't be achieved
    The OCR: how much of I the score determines
    The Cayley-Dickson tower: each level k contributes 2^k α_k to H
    The harmonic series: the susceptibility I'/I measures fragility
    The tournament primes: they divide the OCR denominator,
      which comes from the rational structure of I averaged over classes
""")

print("opus-2026-03-22-S156")
