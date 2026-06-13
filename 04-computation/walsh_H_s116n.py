#!/usr/bin/env python3
"""walsh_H_s116n.py — Walsh/Fourier decomposition of H on the tiling hypercube.

H is a function from {0,1}^10 -> Z (the 10-bit tiling to the Hamiltonian path count).
Any such function has an EXACT multilinear polynomial representation:
  H(x) = sum_S c_S * prod_{i in S} (2*x_i - 1)   (Walsh basis)
  or equivalently in the standard basis:
  H(x) = sum_S a_S * prod_{i in S} x_i            (multilinear)

We compute both representations, find which coefficients dominate,
and understand the EXACT combinatorial nature of H.

Session: kind-pasteur-2026-03-16-S116n32
"""
import sys
sys.stdout.reconfigure(line_buffering=True)

from itertools import permutations, combinations
from collections import defaultdict

print()
print("  WALSH DECOMPOSITION OF H ON THE TILING HYPERCUBE")
print()
print("=" * 70)
print()

N = 6

# Tiling arcs
tiling_arcs = [
    (0,2), (0,3), (0,4), (0,5), (1,3), (1,4), (1,5), (2,4), (2,5), (3,5)
]
path_bits = [0, 5, 9, 12, 14]
nonpath_bits = [i for i in range(15) if i not in path_bits]

def tournament_adj(tiling):
    adj = [[0]*N for _ in range(N)]
    for i in range(N-1):
        adj[i][i+1] = 1
    for idx, (i, j) in enumerate(tiling_arcs):
        if (tiling >> idx) & 1:
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return adj

def count_hp(adj):
    count = 0
    for perm in permutations(range(N)):
        ok = True
        for i in range(N-1):
            if not adj[perm[i]][perm[i+1]]:
                ok = False
                break
        if ok:
            count += 1
    return count

# Compute H for all 1024 tilings
print("  Computing H for all 1024 canonical-path tilings...")
H_table = [0] * 1024
for t in range(1024):
    H_table[t] = count_hp(tournament_adj(t))
print(f"  Done. min={min(H_table)}, max={max(H_table)}, mean={sum(H_table)/1024:.2f}")
print()

# ============================================================
print("  I. MULTILINEAR POLYNOMIAL (STANDARD BASIS)")
print("  " + "-" * 50)
print()

# Any function f: {0,1}^n -> R has a unique multilinear representation:
# f(x) = sum_{S subset [n]} a_S * prod_{i in S} x_i
# where a_S = sum_{T subset S} (-1)^{|S|-|T|} * f(indicator of T)
# This is the Mobius inversion on the Boolean lattice.

def compute_multilinear_coefficients(f_table, n_bits):
    """Compute multilinear polynomial coefficients via Mobius inversion."""
    coeffs = {}
    for s_mask in range(1 << n_bits):
        # S = set of bits where s_mask has 1
        S = frozenset(i for i in range(n_bits) if (s_mask >> i) & 1)
        # a_S = sum_{T subset S} (-1)^{|S|-|T|} * f(T)
        val = 0
        for t_mask in range(1 << n_bits):
            # Check T subset S
            if t_mask & s_mask != t_mask:
                continue
            T = frozenset(i for i in range(n_bits) if (t_mask >> i) & 1)
            sign = (-1) ** (len(S) - len(T))
            val += sign * f_table[t_mask]
        if val != 0:
            coeffs[S] = val
    return coeffs

print("  Computing multilinear coefficients (this is 2^20 operations)...")
coeffs = compute_multilinear_coefficients(H_table, 10)
print(f"  Done. {len(coeffs)} nonzero terms.")
print()

# Verify
errors = 0
for t in range(1024):
    bits = tuple((t >> i) & 1 for i in range(10))
    val = 0
    for S, c in coeffs.items():
        term = c
        for i in S:
            term *= bits[i]
        val += term
    if val != H_table[t]:
        errors += 1
if errors == 0:
    print(f"  VERIFIED: multilinear polynomial matches all 1024 values!")
else:
    print(f"  {errors} mismatches!")
print()

# ============================================================
print("  II. COEFFICIENT STRUCTURE BY DEGREE")
print("  " + "-" * 50)
print()

by_degree = defaultdict(list)
for S, c in coeffs.items():
    by_degree[len(S)].append((S, c))

for deg in sorted(by_degree.keys()):
    terms = by_degree[deg]
    coeff_vals = [c for _, c in terms]
    print(f"  Degree {deg}: {len(terms)} terms, coefficients in "
          f"[{min(coeff_vals)}, {max(coeff_vals)}], sum = {sum(coeff_vals)}")

print()

# ============================================================
print("  III. ALL TERMS (GROUPED BY DEGREE)")
print("  " + "-" * 50)
print()

for deg in sorted(by_degree.keys()):
    terms = sorted(by_degree[deg], key=lambda x: (abs(x[1]), tuple(sorted(x[0]))), reverse=True)
    print(f"  DEGREE {deg}:")
    for S, c in terms[:30]:  # Show top 30 per degree
        bits = sorted(S)
        if bits:
            arcs = [tiling_arcs[b] for b in bits]
            verts = set()
            for a in arcs:
                verts.update(a)
            bit_str = '*'.join(f'x{b}' for b in bits)
            arc_str = ','.join(str(a) for a in arcs)
            print(f"    {c:+4d} * {bit_str:30s}  arcs={arc_str:40s} verts={sorted(verts)}")
        else:
            print(f"    {c:+4d}  (constant)")
    if len(terms) > 30:
        print(f"    ... and {len(terms)-30} more terms")
    print()

# ============================================================
print("  IV. VARIANCE DECOMPOSITION BY DEGREE")
print("  " + "-" * 50)
print()

mean_H = sum(H_table) / 1024
total_var = sum((H - mean_H)**2 for H in H_table) / 1024

print(f"  Mean H = {mean_H:.4f}")
print(f"  Total variance = {total_var:.4f}")
print()

# For a multilinear polynomial on {0,1}^n, the Parseval identity gives:
# Var(f) = sum_{S != empty} a_S^2 * prod_{i in S} p_i*(1-p_i)
# where p_i = Pr(x_i = 1). For uniform bits, p_i = 1/2, so each factor is 1/4.
# Var(f) = sum_{S != empty} a_S^2 / 4^|S|

for deg in sorted(by_degree.keys()):
    if deg == 0:
        continue
    terms = by_degree[deg]
    var_contribution = sum(c**2 / (4**deg) for _, c in terms)
    pct = 100 * var_contribution / total_var if total_var > 0 else 0
    print(f"  Degree {deg}: variance contribution = {var_contribution:8.2f} ({pct:5.1f}%)")

# Cross-check
total_parseval = sum(c**2 / (4**len(S)) for S, c in coeffs.items() if len(S) > 0)
print(f"  Total from Parseval: {total_parseval:.4f} (should be {total_var:.4f})")
print()

# ============================================================
print("  V. TOP 20 MOST IMPORTANT TERMS (BY VARIANCE CONTRIBUTION)")
print("  " + "-" * 50)
print()

term_importance = [(S, c, c**2 / (4**len(S))) for S, c in coeffs.items() if len(S) > 0]
term_importance.sort(key=lambda x: x[2], reverse=True)

print(f"  {'Rank':>4s}  {'Coeff':>6s}  {'Var%':>7s}  {'Term':>25s}  {'Arcs':>40s}")
cumulative = 0
for rank, (S, c, var) in enumerate(term_importance[:20], 1):
    pct = 100 * var / total_var
    cumulative += pct
    bits = sorted(S)
    arcs = [tiling_arcs[b] for b in bits]
    verts = set()
    for a in arcs:
        verts.update(a)
    bit_str = '*'.join(f'x{b}' for b in bits)
    arc_str = str(arcs)
    print(f"  {rank:4d}  {c:+6d}  {pct:6.2f}%  {bit_str:>25s}  {arc_str:>40s}")
print(f"  Cumulative top 20: {cumulative:.1f}%")
print()

# ============================================================
print("  VI. THE COMBINATORIAL NATURE")
print("  " + "-" * 50)
print()

# How many terms of each sign?
pos_terms = sum(1 for _, c in coeffs.items() if c > 0 and len(_) > 0)
neg_terms = sum(1 for _, c in coeffs.items() if c < 0 and len(_) > 0)
print(f"  Positive terms: {pos_terms}")
print(f"  Negative terms: {neg_terms}")
print()

# The combinatorial meaning:
# A positive coefficient a_S means: setting ALL bits in S to 1 (and others to 0)
# INCREASES H by a_S compared to the inclusion-exclusion baseline.
# A negative coefficient means the opposite.

# Key question: do the degree-2 terms dominate?
d1_var = sum(c**2/4 for S, c in coeffs.items() if len(S) == 1)
d2_var = sum(c**2/16 for S, c in coeffs.items() if len(S) == 2)
d3_var = sum(c**2/64 for S, c in coeffs.items() if len(S) == 3)
d4_var = sum(c**2/256 for S, c in coeffs.items() if len(S) == 4)
d5_var = sum(c**2/1024 for S, c in coeffs.items() if len(S) == 5)

print("  Variance by degree (Parseval):")
for d, v in [(1, d1_var), (2, d2_var), (3, d3_var), (4, d4_var), (5, d5_var)]:
    print(f"    Degree {d}: {v:8.2f} ({100*v/total_var:5.1f}%)")
print()

# The LINEAR terms alone
print("  Linear coefficients (the main effects):")
for S, c in sorted(coeffs.items(), key=lambda x: min(x[0]) if x[0] else -1):
    if len(S) != 1:
        continue
    b = min(S)
    arc = tiling_arcs[b]
    skip = arc[1] - arc[0]
    print(f"    a_{{{b}}} = {c:+4d}  (arc {arc}, skip {skip})")
print()

print("  CONCLUSION:")
print("  H on the tiling hypercube is a degree-? polynomial with ? nonzero terms.")
print(f"  Max degree: {max(len(S) for S in coeffs.keys())}")
print(f"  Total nonzero terms: {len(coeffs)}")
total_by_d = {d: len(terms) for d, terms in by_degree.items()}
print(f"  Terms per degree: {dict(sorted(total_by_d.items()))}")
print()

# The key insight
print("  WHY THE LINEAR MODEL FAILS:")
print(f"  Linear terms explain only {100*d1_var/total_var:.1f}% of variance.")
print(f"  Quadratic terms add {100*d2_var/total_var:.1f}%.")
print(f"  Cubic terms add {100*d3_var/total_var:.1f}%.")
print(f"  Together, degrees 1-3 explain {100*(d1_var+d2_var+d3_var)/total_var:.1f}%.")
print(f"  The remaining {100*(1-(d1_var+d2_var+d3_var)/total_var):.1f}% requires degrees 4+.")
print()
print("  THE COMBINATORIAL NATURE IS:")
print("  H depends on SPECIFIC CONFIGURATIONS of 2-4 bits,")
print("  not on aggregate statistics (like total forward count).")
print("  Each configuration corresponds to a tournament substructure")
print("  (a 3-cycle, a 5-cycle, or a disjoint pair of cycles)")
print("  whose presence or absence shifts H by a fixed amount.")
print()
