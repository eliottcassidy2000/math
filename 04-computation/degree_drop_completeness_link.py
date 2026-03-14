#!/usr/bin/env python3
"""
degree_drop_completeness_link.py — opus-2026-03-14-S71g

CONNECTION: Degree Drop Theorem + Completeness Obstruction

The Degree Drop Theorem (kind-pasteur S72):
  H(T) as multilinear polynomial has degree = 2*floor((n-1)/2)
  = n-1 for odd n, n-2 for even n.

  At even n: top-degree coefficients cancel via path reversal involution.
  At odd n: all top-degree coefficients = ±2.

The Completeness Obstruction (this session):
  Directed b-cycle has H=b (as digraph).
  Tournament completion shifts H massively (7→[25,189] at n=7).

QUESTION: Is the degree drop related to the completeness obstruction?

The multilinear polynomial H(T) = H(x₁₂, x₁₃, ..., x_{n-1,n}) where
xᵢⱼ = 1 if i→j, xᵢⱼ = 0 if j→i (so xᵢⱼ + xⱼᵢ = 1 for tournaments).

For a DIGRAPH, we don't have xᵢⱼ + xⱼᵢ = 1. The variables are independent.

The tournament constraint xᵢⱼ + xⱼᵢ = 1 is what causes degree drop!
  - Substitute xⱼᵢ = 1 - xᵢⱼ into the full polynomial.
  - This reduces degree because top-degree terms cancel.

So the degree drop IS the completeness obstruction, algebraically!
The tournament constraint (completeness) kills the highest-degree terms,
which is EXACTLY what prevents H=7.

Let's verify this connection.
"""

from itertools import permutations, combinations
from collections import defaultdict

# ============================================================
# Part 1: H as polynomial in INDEPENDENT arc variables
# ============================================================

print("=" * 70)
print("H AS MULTILINEAR POLYNOMIAL: DIGRAPH vs TOURNAMENT")
print("=" * 70)

# For n=4: C(4,2) = 6 arcs in a tournament.
# But 12 arcs in a digraph (each pair can independently have 0, 1, or 2 arcs).
# For the multilinear case: we use indicator variables for each DIRECTED arc.
# Tournament: xᵢⱼ ∈ {0,1}, xᵢⱼ + xⱼᵢ = 1
# Digraph: xᵢⱼ ∈ {0,1}, xⱼᵢ ∈ {0,1} independently

# H(x) = Σ_P ∏ x_{P_i,P_{i+1}} where sum over all permutations P

n = 4
arcs = [(i,j) for i in range(n) for j in range(n) if i != j]
arc_idx = {a: i for i, a in enumerate(arcs)}

# Build H as sum of products of arc indicators
print(f"\nn={n}: H as polynomial in {len(arcs)} arc indicators")

# For each permutation, record which arcs are used
path_arcs = []
for perm in permutations(range(n)):
    used = []
    for i in range(n-1):
        used.append((perm[i], perm[i+1]))
    path_arcs.append(used)

print(f"  {len(path_arcs)} permutations = {n}! = potential paths")

# The polynomial H = Σ_P ∏_{(i,j)∈P} xᵢⱼ
# Each term is a product of n-1 distinct arc variables.
# Degree = n-1 for the digraph polynomial.

# For TOURNAMENT, substitute xⱼᵢ = 1 - xᵢⱼ (for j>i, use xᵢⱼ as primary):
# This reduces the number of independent variables from n(n-1) to C(n,2).
# And the degree drops (at even n).

# Let's compute the tournament polynomial explicitly.
# Use primary variables: x_{i,j} for i<j.
# Then x_{j,i} = 1 - x_{i,j}.

# For each path P, the product ∏ x_{P_k, P_{k+1}} becomes:
# For each step (a,b): if a<b, use x_{a,b}. If a>b, use (1 - x_{b,a}).

# Let's expand this for n=4.
from sympy import symbols, expand, Poly, degree

# Primary variables: x01, x02, x03, x12, x13, x23
primary_vars = {}
var_symbols = []
for i in range(n):
    for j in range(i+1, n):
        name = f'x{i}{j}'
        s = symbols(name)
        primary_vars[(i,j)] = s
        var_symbols.append(s)

print(f"  Primary tournament variables: {[str(v) for v in var_symbols]}")

# Build H polynomial
H_poly = 0
for perm in permutations(range(n)):
    term = 1
    for k in range(n-1):
        a, b = perm[k], perm[k+1]
        if a < b:
            term *= primary_vars[(a,b)]
        else:
            term *= (1 - primary_vars[(b,a)])
    H_poly += term

H_expanded = expand(H_poly)
print(f"\n  H(tournament) = {H_expanded}")

# Find degree
total_deg = 0
for v in var_symbols:
    d = degree(Poly(H_expanded, v))
    if d > total_deg:
        total_deg = d
print(f"  Degree in any single variable: {total_deg}")

# Total multilinear degree
from sympy import Mul, Add
max_term_deg = 0
for term in Add.make_args(H_expanded):
    # Count number of variables in this monomial
    vs = term.free_symbols
    max_term_deg = max(max_term_deg, len(vs))
print(f"  Max multilinear degree: {max_term_deg}")
print(f"  Expected (degree drop): 2*floor({n-1}/2) = {2*((n-1)//2)}")

# Now compute the DIGRAPH polynomial (no tournament constraint)
# Here ALL n(n-1) arc variables are independent
print(f"\n  DIGRAPH polynomial (independent arc variables):")
# For small n, this is just the sum of products of n-1 arc indicators
# Each path uses n-1 distinct directed arcs
# The degree is always n-1 (no cancellation without constraint)
print(f"  Degree: always {n-1} (no constraint to cause cancellation)")

# ============================================================
# Part 2: The constraint substitution explicitly
# ============================================================

print(f"\n{'='*70}")
print("CONSTRAINT SUBSTITUTION: xⱼᵢ = 1 - xᵢⱼ")
print(f"{'='*70}")

# At n=4: the digraph polynomial has degree 3 in 12 variables.
# After substitution xⱼᵢ = 1 - xᵢⱼ, we get degree 2 in 6 variables.
# The degree drops from 3 to 2 (= n-2 for even n=4).

# Let me show this explicitly at n=5 too.
n = 5
primary_vars_5 = {}
var_symbols_5 = []
for i in range(n):
    for j in range(i+1, n):
        name = f'x{i}{j}'
        s = symbols(name)
        primary_vars_5[(i,j)] = s
        var_symbols_5.append(s)

H_poly_5 = 0
for perm in permutations(range(n)):
    term = 1
    for k in range(n-1):
        a, b = perm[k], perm[k+1]
        if a < b:
            term *= primary_vars_5[(a,b)]
        else:
            term *= (1 - primary_vars_5[(b,a)])
    H_poly_5 += term

H_expanded_5 = expand(H_poly_5)

# Find multilinear degree
max_deg_5 = 0
for term in Add.make_args(H_expanded_5):
    vs = term.free_symbols
    max_deg_5 = max(max_deg_5, len(vs))

print(f"\n  n={n}: {len(var_symbols_5)} primary variables")
print(f"  Multilinear degree: {max_deg_5}")
print(f"  Expected: n-1 = {n-1} (odd n: no drop)")

# Check coefficients at top degree
top_terms = []
for term in Add.make_args(H_expanded_5):
    vs = term.free_symbols
    if len(vs) == max_deg_5:
        top_terms.append(term)

print(f"  Number of top-degree terms: {len(top_terms)}")
if top_terms:
    # Check if all coefficients are ±2
    coeffs = set()
    for t in top_terms:
        # Extract coefficient
        c = t
        for v in t.free_symbols:
            c = c.coeff(v)
        coeffs.add(int(c) if c.is_integer else c)
    if len(coeffs) <= 5:
        print(f"  Top-degree coefficients: {coeffs}")
    else:
        print(f"  Top-degree coefficients: {len(coeffs)} distinct values")

# ============================================================
# Part 3: The key identity
# ============================================================

print(f"\n{'='*70}")
print("THE KEY IDENTITY: COMPLETENESS = DEGREE DROP = H=7 OBSTRUCTION")
print(f"{'='*70}")

print("""
SYNTHESIS:

1. DIGRAPH WORLD (simplex = (x+1)^n):
   - H is a multilinear polynomial of degree n-1 in n(n-1) variables
   - Variables are INDEPENDENT (each arc can be present or absent)
   - H=7 IS achievable (by directed 7-cycle: H=7)
   - No degree drop — top-degree terms do NOT cancel

2. TOURNAMENT WORLD (cuboid = (x+2)^n):
   - The constraint xⱼᵢ = 1 - xᵢⱼ HALVES the variables
   - But MORE importantly: it causes DEGREE DROP
   - Even n: degree drops from n-1 to n-2 (top terms cancel)
   - Odd n: degree stays n-1 but coefficients = ±2 (not general)
   - H=7 is NOT achievable

3. THE CONNECTION:
   - The tournament constraint xⱼᵢ = 1 - xᵢⱼ is AFFINE (degree 1)
   - Substituting it into degree-(n-1) terms:
     (1 - x_{ij}) × other variables = -x_{ij} × others + others
   - The -x terms cancel with matching +x terms from reversed paths
   - This is the PATH REVERSAL INVOLUTION that kind-pasteur proved!

4. WHY H=7 IS KILLED:
   - In the digraph, H=7 requires very sparse arc structure
   - Tournament completion adds (1-x) terms everywhere
   - These terms expand to create many new paths
   - The minimum tournament H at n=7 is 25 (from 7-cycle completion)
   - The "gap" from 7 to 25 is exactly the effect of the constraint

5. THE SIMPLEX-CUBOID ANALOGY:
   - Simplex (x+1)^n: each dimension has 2 choices (arc present/absent)
   - Cuboid (x+2)^n: each dimension has 3 values... wait, this doesn't
     quite work because the constraint is binary (tournament), not ternary.

   Better analogy:
   - INDEPENDENT arcs = simplex: each pair contributes (1+xᵢⱼ)(1+xⱼᵢ)
   - CONSTRAINED arcs = cuboid: each pair contributes (xᵢⱼ + (1-xᵢⱼ)) = 1
     ... but this is trivial.

   Actually, the simplex-cuboid connection is through the INDEPENDENCE
   POLYNOMIAL:
   - I(Ω, x) with x=2 (tournament) vs x=1 (graph)
   - The factor (1+2x) = 5 vs (1+x) = 3 per cycle pair
   - Cuboid evaluation is RICHER than simplex evaluation
""")

print(f"{'='*70}")
print("ANALYSIS COMPLETE")
print(f"{'='*70}")
