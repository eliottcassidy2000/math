#!/usr/bin/env python3
"""tiling_polynomial_s116n.py — The EXACT polynomial for H in tiling bits.

At n=6, H = 1 + 2*c3 + 4*alpha_2 (OCF truncation).
c3 = sum of 20 triple indicators (each depends on <= 3 bits).
alpha_2 = sum of 10 disjoint-pair products (each depends on <= 6 bits).

H is therefore a degree-6 multilinear polynomial in the 10 tiling bits.
This script computes the polynomial EXACTLY, reveals its structure,
and explains the 94% of variance the linear model misses.

Session: kind-pasteur-2026-03-16-S116n32
"""
import sys
sys.stdout.reconfigure(line_buffering=True)

from itertools import permutations, combinations
from collections import Counter, defaultdict

print()
print("  THE EXACT POLYNOMIAL FOR H(tiling)")
print()
print("=" * 70)
print()

N = 6

# ============================================================
# Define the tiling structure
# ============================================================

# Path: 0->1->2->3->4->5 (always forward)
# Non-path arcs (the 10 tiling bits):
# Index -> (i, j) pair
tiling_arcs = [
    (0, 2),  # bit 0, skip 2
    (0, 3),  # bit 1, skip 3
    (0, 4),  # bit 2, skip 4
    (0, 5),  # bit 3, skip 5
    (1, 3),  # bit 4, skip 2
    (1, 4),  # bit 5, skip 3
    (1, 5),  # bit 6, skip 4
    (2, 4),  # bit 7, skip 2
    (2, 5),  # bit 8, skip 3
    (3, 5),  # bit 9, skip 2
]

# Map (i,j) -> bit index (for i<j)
arc_to_bit = {}
for idx, (i, j) in enumerate(tiling_arcs):
    arc_to_bit[(i, j)] = idx

# Path arcs: (0,1), (1,2), (2,3), (3,4), (4,5) — always forward (i->j)
path_arcs = {(0,1), (1,2), (2,3), (3,4), (4,5)}

# ============================================================
print("  I. THE 20 TRIPLES AND THEIR CYCLICITY CONDITIONS")
print("  " + "-" * 50)
print()

# For each triple {a,b,c} (a<b<c), determine:
# - Which arcs are path arcs (fixed forward) and which are tiling bits
# - The condition for the triple to be cyclic

def triple_cycle_condition(a, b, c):
    """
    For triple {a,b,c} with a<b<c, determine the cyclicity condition.
    Returns a list of (tiling_bit_conditions, direction) pairs.
    A tiling bit condition maps bit_index -> required_value (0 or 1).
    bit=1 means i->j (forward), bit=0 means j->i (backward).

    The 3 arcs are (a,b), (b,c), (a,c).
    Arc (i,j) with i<j: path arc iff j=i+1.
    Path arc is always forward (i->j), so: a->b iff (a,b) is path or bit=1.

    For cycle a->b->c->a: need a->b, b->c, c->a.
      a->b: path iff b=a+1 (always true), else bit[(a,b)]=1
      b->c: path iff c=b+1 (always true), else bit[(b,c)]=1
      c->a: NOT path (since a<c), so bit[(a,c)]=0 (meaning c->a)

    For cycle a->c->b->a: need a->c, c->b, b->a.
      a->c: NOT path if c>a+1, so bit[(a,c)]=1
             if c=a+1... no, a<b<c so c>=a+2.
      c->b: NOT path if c>b+1. If c=b+1: path gives b->c (wrong direction!).
             So need c->b: if path (c=b+1), then it's b->c (forward), so c->b is IMPOSSIBLE.
             If not path: bit[(b,c)]=0.
      b->a: NOT path if b>a+1. If b=a+1: path gives a->b, so b->a is IMPOSSIBLE.
             If not path: bit[(a,b)]=0.
    """
    arcs = [(a,b), (b,c), (a,c)]
    is_path = [(i,j) in path_arcs for (i,j) in arcs]

    conditions = []

    # Cycle 1: a->b->c->a
    cond1 = {}
    possible1 = True
    # a->b: path gives a->b (good), or bit=1 (good)
    if not is_path[0]:  # (a,b) not path
        cond1[arc_to_bit[(a,b)]] = 1  # need bit=1 for a->b
    # b->c: path gives b->c (good), or bit=1 (good)
    if not is_path[1]:  # (b,c) not path
        cond1[arc_to_bit[(b,c)]] = 1
    # c->a: NOT a path arc (since a<c and c>=a+2). Need bit=0 for c->a.
    cond1[arc_to_bit[(a,c)]] = 0
    conditions.append(cond1)

    # Cycle 2: a->c->b->a
    cond2 = {}
    possible2 = True
    # a->c: not path (c>=a+2). Need bit=1 for a->c.
    cond2[arc_to_bit[(a,c)]] = 1
    # c->b: if (b,c) is path, then b->c is forced, so c->b IMPOSSIBLE.
    if is_path[1]:
        possible2 = False
    else:
        cond2[arc_to_bit[(b,c)]] = 0  # need c->b
    # b->a: if (a,b) is path, then a->b is forced, so b->a IMPOSSIBLE.
    if is_path[0]:
        possible2 = False
    else:
        cond2[arc_to_bit[(a,b)]] = 0  # need b->a

    result = [cond1]
    if possible2:
        result.append(cond2)

    return result

print(f"  {'Triple':>10s}  {'#path':>6s}  {'Cycle 1 condition':>30s}  {'Cycle 2':>20s}")
print(f"  {'------':>10s}  {'-----':>6s}  {'-'*30:>30s}  {'-'*20:>20s}")

all_triple_conditions = []
for a in range(N):
    for b in range(a+1, N):
        for c in range(b+1, N):
            n_path = sum(1 for (i,j) in [(a,b),(b,c)] if (i,j) in path_arcs)
            conditions = triple_cycle_condition(a, b, c)

            cond1_str = str(conditions[0])
            cond2_str = str(conditions[1]) if len(conditions) > 1 else "IMPOSSIBLE"

            print(f"  {{{a},{b},{c}}}  {n_path:6d}  {cond1_str:>30s}  {cond2_str:>20s}")

            all_triple_conditions.append((a, b, c, conditions))
print()

# ============================================================
print("  II. c3 AS A SUM OF BIT-CONDITIONS")
print("  " + "-" * 50)
print()

# For each triple, the cyclicity indicator is:
# I(triple) = indicator(all bits match condition_1) OR indicator(all bits match condition_2)
# Since the two cycles are mutually exclusive (they require opposite orientations of
# the longest arc), exactly one or zero can be satisfied.

# Actually: condition_1 requires bit[(a,c)]=0, condition_2 requires bit[(a,c)]=1.
# So they are MUTUALLY EXCLUSIVE. The triple is cyclic iff one of them holds.

# For a triple with 2 path arcs (consecutive: {i,i+1,i+2}):
#   Cycle 1: bit[(i,i+2)]=0 (backward skip-2)
#   Cycle 2: IMPOSSIBLE (path arcs force direction)
#   So: I = (1 - bit_k) where k is the skip-2 bit for (i,i+2).

# For a triple with 1 path arc:
#   More complex — depends on which arc is the path arc.

# For a triple with 0 path arcs:
#   Both cycles possible. I = indicator(cond1) + indicator(cond2).

# Let me express each indicator as a polynomial in the bits.
# bit values are 0 or 1, so:
# indicator(bit_k = v) = (1-bit_k) if v=0, or bit_k if v=1.
# Product of such indicators over all bits in a condition.

def condition_to_polynomial(cond):
    """Convert a bit-condition dict to a monomial term.
    Returns (coefficient, frozenset of bit indices) where the monomial is:
    coefficient * prod_{k in pos_bits} x_k * prod_{k in neg_bits} (1 - x_k)
    Expanded, this gives a sum of monomials with +/- coefficients.
    """
    # Each factor: x_k (if cond[k]=1) or (1-x_k) (if cond[k]=0).
    # Product = sum over subsets S of neg_bits: (-1)^|S| * prod(x_k for k in pos_bits ∪ S)

    pos_bits = frozenset(k for k, v in cond.items() if v == 1)
    neg_bits = frozenset(k for k, v in cond.items() if v == 0)

    terms = []  # list of (coefficient, frozenset of bit indices)
    for r in range(len(neg_bits) + 1):
        for subset in combinations(neg_bits, r):
            coeff = (-1) ** r
            bits = pos_bits | frozenset(subset)
            terms.append((coeff, bits))

    return terms

# Build the full polynomial for c3
c3_polynomial = defaultdict(int)  # monomial (frozenset of bits) -> coefficient

for a, b, c, conditions in all_triple_conditions:
    for cond in conditions:
        terms = condition_to_polynomial(cond)
        for coeff, bits in terms:
            c3_polynomial[bits] += coeff

# Simplify: remove zero terms
c3_polynomial = {k: v for k, v in c3_polynomial.items() if v != 0}

# Sort by degree
by_degree = defaultdict(list)
for bits, coeff in c3_polynomial.items():
    by_degree[len(bits)].append((bits, coeff))

print(f"  c3 polynomial has {len(c3_polynomial)} terms:")
for deg in sorted(by_degree.keys()):
    terms = by_degree[deg]
    print(f"    Degree {deg}: {len(terms)} terms")
    if len(terms) <= 15:
        for bits, coeff in sorted(terms, key=lambda x: tuple(sorted(x[0]))):
            bit_str = '*'.join(f'x{b}' for b in sorted(bits)) if bits else '1'
            print(f"      {coeff:+3d} * {bit_str}")
print()

# ============================================================
print("  III. VERIFY THE POLYNOMIAL")
print("  " + "-" * 50)
print()

# Check against brute-force c3 computation
def eval_polynomial(poly, bits_tuple):
    """Evaluate multilinear polynomial at a point."""
    result = 0
    for monomial_bits, coeff in poly.items():
        term = coeff
        for b in monomial_bits:
            term *= bits_tuple[b]
        result += term
    return result

def compute_c3_brute(bits_tuple):
    """Compute c3 by checking all triples."""
    # Reconstruct adjacency
    adj = [[0]*N for _ in range(N)]
    for i in range(N-1):
        adj[i][i+1] = 1  # path
    for idx, (i, j) in enumerate(tiling_arcs):
        if bits_tuple[idx]:
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    c3 = 0
    for i in range(N):
        for j in range(i+1, N):
            for k in range(j+1, N):
                if adj[i][j] and adj[j][k] and adj[k][i]:
                    c3 += 1
                elif adj[i][k] and adj[k][j] and adj[j][i]:
                    c3 += 1
    return c3

errors = 0
for tiling in range(1024):
    bits = tuple((tiling >> i) & 1 for i in range(10))
    poly_val = eval_polynomial(c3_polynomial, bits)
    brute_val = compute_c3_brute(bits)
    if poly_val != brute_val:
        errors += 1
        if errors <= 3:
            print(f"  ERROR at {bits}: poly={poly_val}, brute={brute_val}")

if errors == 0:
    print(f"  VERIFIED: c3 polynomial matches brute force for all 1024 tilings!")
else:
    print(f"  {errors} mismatches found!")
print()

# ============================================================
print("  IV. THE STRUCTURE OF THE c3 POLYNOMIAL")
print("  " + "-" * 50)
print()

# Constant term = c3 when all bits = 0 (all non-path arcs backward)
constant = c3_polynomial.get(frozenset(), 0)
print(f"  Constant term (all backward): {constant}")
print(f"  This means: with ALL non-path arcs backward, c3 = {constant}")
print()

# Linear terms: coefficient of each individual bit
print("  Linear terms (effect of flipping bit k from 0 to 1):")
for k in range(10):
    coeff = c3_polynomial.get(frozenset([k]), 0)
    arc = tiling_arcs[k]
    skip = arc[1] - arc[0]
    print(f"    x{k} (arc {arc}, skip {skip}): coefficient = {coeff:+d}")
print()

# The linear terms tell us: flipping bit k from 0 to 1 changes c3 by this amount
# (PLUS higher-order interaction corrections).

# Quadratic terms
print("  Quadratic terms (pairwise interactions):")
quad_terms = [(bits, coeff) for bits, coeff in c3_polynomial.items() if len(bits) == 2]
for bits, coeff in sorted(quad_terms, key=lambda x: tuple(sorted(x[0]))):
    b1, b2 = sorted(bits)
    a1, a2 = tiling_arcs[b1], tiling_arcs[b2]
    # Do these arcs share a vertex?
    shared = set(a1) & set(a2)
    sharing = f"share {shared}" if shared else "disjoint"
    print(f"    x{b1}*x{b2} ({a1}x{a2}, {sharing}): {coeff:+d}")
print()

# Cubic terms
print("  Cubic terms (3-way interactions):")
cubic_terms = [(bits, coeff) for bits, coeff in c3_polynomial.items() if len(bits) == 3]
for bits, coeff in sorted(cubic_terms, key=lambda x: tuple(sorted(x[0]))):
    bs = sorted(bits)
    arcs = [tiling_arcs[b] for b in bs]
    # These 3 arcs form a triangle on 3 vertices
    vertices = set()
    for a in arcs:
        vertices.update(a)
    print(f"    x{'*x'.join(str(b) for b in bs)} (arcs {arcs}): {coeff:+d}, "
          f"vertices={sorted(vertices)}")
print()

# ============================================================
print("  V. NOW THE FULL H POLYNOMIAL")
print("  " + "-" * 50)
print()

# H = 1 + 2*c3 + 4*alpha_2
# We have c3 as a polynomial. Now we need alpha_2.

# alpha_2 = number of vertex-disjoint cyclic triple pairs
# The partitions of {0,1,2,3,4,5} into two groups of 3:
partitions = []
for triple1 in combinations(range(N), 3):
    triple2 = tuple(v for v in range(N) if v not in triple1)
    if triple1 < triple2:  # avoid double counting
        partitions.append((triple1, triple2))

print(f"  {len(partitions)} vertex-disjoint triple partitions:")
for t1, t2 in partitions:
    print(f"    {t1} | {t2}")
print()

# For each partition, alpha_2 gets +1 if BOTH triples are cyclic.
# indicator(t1 cyclic) * indicator(t2 cyclic)

# Build alpha_2 polynomial
alpha2_polynomial = defaultdict(int)

for t1, t2 in partitions:
    a1, b1, c1 = t1
    a2, b2, c2 = t2

    # Get conditions for each triple
    conds1 = [c for _, _, _, conds in all_triple_conditions
              if set([_]) == set() or True  # placeholder
              for c in conds]  # This is wrong, let me fix

    # Actually, find the conditions for triple t1
    conds_t1 = None
    conds_t2 = None
    for a, b, c, conds in all_triple_conditions:
        if (a, b, c) == (a1, b1, c1):
            conds_t1 = conds
        if (a, b, c) == (a2, b2, c2):
            conds_t2 = conds

    # indicator(t1 cyclic) = sum of product indicators over conditions
    # indicator(t2 cyclic) = sum of product indicators over conditions
    # Product = sum over (c1, c2) pairs

    for c1 in conds_t1:
        terms1 = condition_to_polynomial(c1)
        for c2 in conds_t2:
            terms2 = condition_to_polynomial(c2)
            # Multiply
            for coeff1, bits1 in terms1:
                for coeff2, bits2 in terms2:
                    # Since triples are vertex-disjoint, their arcs don't overlap
                    combined_bits = bits1 | bits2
                    alpha2_polynomial[combined_bits] += coeff1 * coeff2

# Simplify
alpha2_polynomial = {k: v for k, v in alpha2_polynomial.items() if v != 0}

# H polynomial = 1 + 2*c3 + 4*alpha_2
H_polynomial = defaultdict(int)
H_polynomial[frozenset()] = 1  # constant 1
for bits, coeff in c3_polynomial.items():
    H_polynomial[bits] += 2 * coeff
for bits, coeff in alpha2_polynomial.items():
    H_polynomial[bits] += 4 * coeff

H_polynomial = {k: v for k, v in H_polynomial.items() if v != 0}

print(f"  H polynomial has {len(H_polynomial)} terms:")
H_by_degree = defaultdict(list)
for bits, coeff in H_polynomial.items():
    H_by_degree[len(bits)].append((bits, coeff))

for deg in sorted(H_by_degree.keys()):
    terms = H_by_degree[deg]
    print(f"    Degree {deg}: {len(terms)} terms")
print()

# Verify H polynomial
print("  Verifying H polynomial...")
def count_hp_brute(bits_tuple):
    adj = [[0]*N for _ in range(N)]
    for i in range(N-1):
        adj[i][i+1] = 1
    for idx, (i, j) in enumerate(tiling_arcs):
        if bits_tuple[idx]:
            adj[i][j] = 1
        else:
            adj[j][i] = 1
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

errors = 0
for tiling in range(1024):
    bits = tuple((tiling >> i) & 1 for i in range(10))
    poly_val = eval_polynomial(H_polynomial, bits)
    brute_val = count_hp_brute(bits)
    if poly_val != brute_val:
        errors += 1
        if errors <= 3:
            print(f"  ERROR at {bits}: poly={poly_val}, brute={brute_val}")

if errors == 0:
    print(f"  VERIFIED: H polynomial matches brute force for all 1024 tilings!")
else:
    print(f"  {errors} mismatches found!")
print()

# ============================================================
print("  VI. THE H POLYNOMIAL COEFFICIENTS")
print("  " + "-" * 50)
print()

# Print all terms by degree
print("  DEGREE 0 (constant):")
for bits, coeff in H_by_degree[0]:
    print(f"    {coeff}")
print()

print("  DEGREE 1 (linear effects of each bit):")
for bits, coeff in sorted(H_by_degree[1], key=lambda x: min(x[0])):
    b = min(bits)
    arc = tiling_arcs[b]
    skip = arc[1] - arc[0]
    print(f"    {coeff:+3d} * x{b} (arc {arc}, skip {skip})")

# Sum of linear coefficients
lin_sum = sum(coeff for _, coeff in H_by_degree[1])
print(f"  Sum of linear coefficients: {lin_sum}")
print()

print("  DEGREE 2 (pairwise interactions):")
for bits, coeff in sorted(H_by_degree[2], key=lambda x: tuple(sorted(x[0]))):
    bs = sorted(bits)
    arcs = [tiling_arcs[b] for b in bs]
    shared = set(arcs[0]) & set(arcs[1])
    print(f"    {coeff:+3d} * x{bs[0]}*x{bs[1]} ({arcs[0]}x{arcs[1]}, "
          f"{'share ' + str(shared) if shared else 'disjoint'})")
print()

if 3 in H_by_degree:
    print(f"  DEGREE 3: {len(H_by_degree[3])} terms")
    for bits, coeff in sorted(H_by_degree[3], key=lambda x: tuple(sorted(x[0]))):
        bs = sorted(bits)
        vertices = set()
        for b in bs:
            vertices.update(tiling_arcs[b])
        print(f"    {coeff:+3d} * x{'*x'.join(str(b) for b in bs)} (verts {sorted(vertices)})")
    print()

for deg in range(4, 7):
    if deg in H_by_degree:
        print(f"  DEGREE {deg}: {len(H_by_degree[deg])} terms")
        # Just show coefficient distribution
        coeffs = [c for _, c in H_by_degree[deg]]
        print(f"    Coefficients: {sorted(set(coeffs))}")
        print(f"    Sum: {sum(coeffs)}")
        print()

# ============================================================
print("  VII. VARIANCE DECOMPOSITION BY DEGREE")
print("  " + "-" * 50)
print()

# How much of H's variance is explained by each degree?
# Decompose H into: H = H_0 + H_1(x) + H_2(x) + H_3(x) + ...
# where H_k contains only degree-k terms.

# Compute variance contributed by each degree
mean_H = sum(count_hp_brute(tuple((t >> i) & 1 for i in range(10)))
             for t in range(1024)) / 1024

total_var = 0
for t in range(1024):
    bits = tuple((t >> i) & 1 for i in range(10))
    H = eval_polynomial(H_polynomial, bits)
    total_var += (H - mean_H) ** 2
total_var /= 1024

print(f"  Mean H = {mean_H:.4f}")
print(f"  Total variance = {total_var:.4f}")
print()

# For each degree, compute the variance of that degree's contribution
for deg in sorted(H_by_degree.keys()):
    if deg == 0:
        continue  # constant doesn't contribute to variance
    # Build sub-polynomial for this degree
    sub_poly = {bits: coeff for bits, coeff in H_polynomial.items() if len(bits) == deg}
    # Compute variance of this component
    vals = []
    for t in range(1024):
        bits = tuple((t >> i) & 1 for i in range(10))
        vals.append(eval_polynomial(sub_poly, bits))
    mean_sub = sum(vals) / 1024
    var_sub = sum((v - mean_sub)**2 for v in vals) / 1024
    pct = 100 * var_sub / total_var if total_var > 0 else 0

    # Covariance with H
    cov_with_H = 0
    for t in range(1024):
        bits = tuple((t >> i) & 1 for i in range(10))
        H = eval_polynomial(H_polynomial, bits)
        sub_val = eval_polynomial(sub_poly, bits)
        cov_with_H += (H - mean_H) * (sub_val - mean_sub)
    cov_with_H /= 1024

    print(f"  Degree {deg}: variance = {var_sub:8.2f} ({pct:5.1f}% of total), "
          f"cov(H, H_{deg}) = {cov_with_H:8.2f}")
print()

# ============================================================
print("  VIII. THE KEY INSIGHT: WHAT DRIVES H?")
print("  " + "-" * 50)
print()

# The degree-1 (linear) terms explain very little.
# The degree-2 (quadratic) terms should explain most of the variance.
# The degree-3 (cubic) terms come from the multilinear structure of c3.
# Degrees 4-6 come from alpha_2 (disjoint pairs).

# What are the LARGEST coefficients in the H polynomial?
all_terms = sorted(H_polynomial.items(), key=lambda x: abs(x[1]), reverse=True)
print("  Top 20 largest-magnitude coefficients:")
for bits, coeff in all_terms[:20]:
    bs = sorted(bits)
    arcs = [tiling_arcs[b] for b in bs]
    deg = len(bits)
    vertices = set()
    for b in bs:
        vertices.update(tiling_arcs[b])
    bit_str = '*'.join(f'x{b}' for b in bs) if bs else '1'
    print(f"    {coeff:+4d} * {bit_str:20s}  deg={deg}, verts={sorted(vertices)}")
print()

# ============================================================
print("  IX. THE COMBINATORIAL NATURE REVEALED")
print("  " + "-" * 50)
print()

# Count how many terms involve each bit
bit_involvement = Counter()
for bits, coeff in H_polynomial.items():
    for b in bits:
        bit_involvement[b] += 1

print("  How many H-polynomial terms involve each bit:")
for b in range(10):
    arc = tiling_arcs[b]
    skip = arc[1] - arc[0]
    print(f"    x{b} (arc {arc}, skip {skip}): in {bit_involvement[b]} terms")
print()

# The "effective degree" of each bit: weighted by coefficient magnitude
bit_weight = defaultdict(float)
for bits, coeff in H_polynomial.items():
    for b in bits:
        bit_weight[b] += abs(coeff)

print("  Total coefficient weight involving each bit:")
for b in range(10):
    arc = tiling_arcs[b]
    skip = arc[1] - arc[0]
    print(f"    x{b} (arc {arc}, skip {skip}): weight = {bit_weight[b]:.0f}")
print()
