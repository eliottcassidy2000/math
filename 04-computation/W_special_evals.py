#!/usr/bin/env python3
"""
Special evaluations of W(r) — what does the W-polynomial encode?

Known:
  W(1/2) = H(T) (forward Ham paths)
  W(-1/2) = (-1)^{n-1} H(T) (backward Ham paths with sign)
  W(0) = "tangent evaluation" (involves tangent numbers)

What about:
  W(1) = ?
  W(-1) = ?
  W(i/2) where i = sqrt(-1)?
  W at roots of F_f(r)?

Also: the DERIVATIVE W'(r) at special points.

And: what is the RESULTANT of W(r) and its flip?

kind-pasteur-2026-03-07-S26
"""
from itertools import permutations, combinations
from fractions import Fraction
from collections import defaultdict

def tournament_from_tiling(n, tiling_bits):
    A = [[0]*n for _ in range(n)]
    for i in range(n-1):
        A[i][i+1] = 1
    idx = 0
    for i in range(n):
        for j in range(i+2, n):
            if (tiling_bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def num_tiling_bits(n):
    return n*(n-1)//2 - (n-1)

def compute_W_poly(A, n):
    result = defaultdict(lambda: Fraction(0))
    for perm in permutations(range(n)):
        poly = {0: Fraction(1)}
        for i in range(n-1):
            s = Fraction(A[perm[i]][perm[i+1]]) - Fraction(1, 2)
            new_poly = {}
            for power, coeff in poly.items():
                new_poly[power+1] = new_poly.get(power+1, Fraction(0)) + coeff
                new_poly[power] = new_poly.get(power, Fraction(0)) + coeff * s
            poly = new_poly
        for power, coeff in poly.items():
            result[power] += coeff
    return dict(result)

def poly_eval(poly, r):
    return sum(c * r**p for p, c in poly.items())

def count_3cycles(A, n):
    t3 = 0
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]: t3 += 1
        if A[i][k] and A[k][j] and A[j][i]: t3 += 1
    return t3

n = 5
m = num_tiling_bits(n)

print(f"SPECIAL W(r) EVALUATIONS at n={n}")
print("="*70)

# Collect W polynomials for all tournaments
evals = {}
for bits in range(2**m):
    A = tournament_from_tiling(n, bits)
    W = compute_W_poly(A, n)
    t3 = count_3cycles(A, n)
    H = poly_eval(W, Fraction(1,2))

    # Special evaluations
    W0 = poly_eval(W, Fraction(0))
    W1 = poly_eval(W, Fraction(1))
    Wm1 = poly_eval(W, Fraction(-1))
    W2 = poly_eval(W, Fraction(2))
    W_third = poly_eval(W, Fraction(1,3))

    if bits < 10 or bits == 2**m - 1:
        print(f"  bits={bits}: t3={t3}, H={H}, W(0)={W0}, W(1)={W1}, W(-1)={Wm1}, W(2)={W2}")

    key = (int(t3),)
    if key not in evals:
        evals[key] = []
    evals[key].append({
        'bits': bits, 'H': H, 'W0': W0, 'W1': W1, 'Wm1': Wm1, 'W2': W2
    })

print(f"\n{'='*70}")
print(f"W(1) and W(-1) by t3")
print(f"{'='*70}")

for t3_val in sorted(evals.keys()):
    entries = evals[t3_val]
    W1_vals = sorted(set(int(e['W1']) for e in entries))
    Wm1_vals = sorted(set(int(e['Wm1']) for e in entries))
    print(f"  t3={t3_val[0]}: W(1) in {W1_vals}, W(-1) in {Wm1_vals}")

# W(1) interpretation:
# At r=1, each factor is (1 + s_i) where s_i in {+1/2, -1/2}
# = 3/2 if forward edge, 1/2 if backward edge
# So W(1) = sum_perm prod (3/2 or 1/2)
# = (1/2)^{n-1} sum_perm prod (3 or 1) based on edge direction
# = (1/2)^{n-1} sum_perm 3^{ascents} where ascent = forward edge

print(f"\n  W(1) = (1/2)^{{{n-1}}} * sum_perm 3^{{forward_edges}}")
print(f"  This is a WEIGHTED Hamiltonian count: 3^k * (number of perms with k forward edges)")

# Let's decompose W(1) by number of forward edges
print(f"\n{'='*70}")
print(f"W(r) as Eulerian-weighted sum: sum_k a_k * (r+1/2)^k * (r-1/2)^{{n-1-k}}")
print(f"{'='*70}")
print(f"  This is exactly the THM-059 Eulerian number interpretation!")
print(f"  a_k = #{'{'}perms with k forward edges{'}'} = Eulerian number of the TOURNAMENT")

# For a specific tournament, compute the Eulerian decomposition
for bits in [0, 1, 2**m - 1, 42]:
    A = tournament_from_tiling(n, bits)
    t3 = count_3cycles(A, n)
    H = poly_eval(compute_W_poly(A, n), Fraction(1,2))

    forward_counts = defaultdict(int)
    for perm in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if A[perm[i]][perm[i+1]])
        forward_counts[fwd] += 1

    print(f"\n  bits={bits}, t3={t3}, H={H}:")
    print(f"    Forward edge distribution: {dict(sorted(forward_counts.items()))}")
    print(f"    Sum of forward_counts = {sum(forward_counts.values())} = {n}!")
    print(f"    Sum where fwd=n-1 = {forward_counts.get(n-1, 0)} = H(T)")

    # W(r) = sum_k a_k * (r+1/2)^{n-1-k} * (r-1/2)^k  [k = descents = backward edges]
    # This is the "Eulerian polynomial" of the tournament!
    print(f"    Tournament Eulerian poly: sum a_k * p^{{n-1-k}} * q^k")
    for k in sorted(forward_counts.keys()):
        print(f"      k={n-1-k} descents: a={forward_counts[k]}")

# Grand sum: sum over all tournaments of W(r)
print(f"\n{'='*70}")
print(f"GRAND SUM: sum over ALL n={n} tournaments of W(r)")
print(f"{'='*70}")

grand_sum = defaultdict(lambda: Fraction(0))
for bits in range(2**m):
    A = tournament_from_tiling(n, bits)
    W = compute_W_poly(A, n)
    for p, c in W.items():
        grand_sum[p] += c

print(f"  Grand sum polynomial: ", end="")
terms = []
for p in sorted(grand_sum.keys(), reverse=True):
    c = grand_sum[p]
    if c != 0:
        terms.append(f"{c}*r^{p}")
print(" + ".join(terms))

# Expected: 2^m * F_{n-1}(r) (the C_0 contribution, since invariants average to specific values)
print(f"  2^m = {2**m}")
print(f"  Grand sum / 2^m:")
for p in sorted(grand_sum.keys(), reverse=True):
    print(f"    r^{p}: {grand_sum[p] / 2**m} = {float(grand_sum[p] / 2**m)}")
