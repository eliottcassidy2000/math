"""
c-TOURNAMENT SYMMETRY: PROOF ROUTE EXPLORATION

The discovery: M[a,b] = M[b,a] whenever t_ij + t_ji = c (uniform constant).

Key algebraic reformulation:
  If A is the adjacency matrix with A + A^T = c*(J-I), then M is symmetric.

  Setting t_ij = c/2 + s_ij where s_ij = -s_ji (skew-symmetric perturbation),
  every arc weight is c/2 + s_ij. The c/2 part is "base" and s_ij is "skew".

This script explores the proof route through this decomposition.
"""

from itertools import permutations
from sympy import symbols, expand, Rational
from collections import defaultdict

# ==========================================================
# Part 1: The c/2 + s decomposition at n=4
# ==========================================================
print("=" * 70)
print("c-TOURNAMENT PROOF ROUTE")
print("=" * 70)

print("\n--- Part 1: Skew decomposition at n=4 ---")

n = 4
c = symbols('c')

# Create skew-symmetric variables: s_ij = -s_ji
# Use s_ij for i < j; s_ji = -s_ij
svars = {}
for i in range(n):
    for j in range(i+1, n):
        svars[(i,j)] = symbols(f's{i}{j}')

def t_val(i, j):
    """Arc weight t_ij = c/2 + s_ij."""
    if i == j: return 0
    if i < j:
        return c/2 + svars[(i,j)]
    else:
        return c/2 - svars[(j,i)]  # s_ji = -s_ij

def hp_paths(vertex_set, start=None, end=None):
    vl = sorted(vertex_set)
    k = len(vl)
    if k == 0: return 0
    if k == 1:
        if start is not None and vl[0] != start: return 0
        if end is not None and vl[0] != end: return 0
        return 1
    total = 0
    for perm in permutations(vl):
        if start is not None and perm[0] != start: continue
        if end is not None and perm[-1] != end: continue
        prod = 1
        for i in range(len(perm)-1):
            prod *= t_val(perm[i], perm[i+1])
        total += prod
    return expand(total)

a, b = 0, 1
U = [v for v in range(n) if v != a and v != b]

# Compute M[a,b]
M_ab = 0
for mask in range(1 << len(U)):
    S = [U[i] for i in range(len(U)) if mask & (1 << i)]
    R = [U[i] for i in range(len(U)) if not (mask & (1 << i))]
    sign = (-1) ** len(S)
    ea = hp_paths(set(S)|{a}, end=a)
    bb = hp_paths(set(R)|{b}, start=b)
    M_ab = expand(M_ab + sign * ea * bb)

print(f"  M[0,1] in (c, s_ij) coordinates:")
print(f"  = {M_ab}")

# Under the skew constraint, M[b,a] should equal M[a,b]
# The swap (a <-> b) means: s_02 <-> s_12, s_03 <-> s_13, s_23 stays
# AND the skew symmetry under swap: paths through a become paths through b
# Let's compute M[1,0] directly
M_ba = 0
for mask in range(1 << len(U)):
    S = [U[i] for i in range(len(U)) if mask & (1 << i)]
    R = [U[i] for i in range(len(U)) if not (mask & (1 << i))]
    sign = (-1) ** len(S)
    eb = hp_paths(set(S)|{b}, end=b)
    ba = hp_paths(set(R)|{a}, start=a)
    M_ba = expand(M_ba + sign * eb * ba)

diff = expand(M_ab - M_ba)
print(f"\n  M[0,1] - M[1,0] = {diff}")

# ==========================================================
# Part 2: Structure of M[0,1] in skew coordinates
# ==========================================================
print(f"\n--- Part 2: Monomial structure in (c, s) ---")

# Collect by powers of c and degree in s
poly_dict = M_ab.as_coefficients_dict()

by_c_power = defaultdict(list)
for monom, coeff in poly_dict.items():
    monom_str = str(monom)
    # Count power of c
    c_power = 0
    remaining = monom
    while remaining != expand(remaining / c) * c or False:
        break
    # Use a different approach: separate c from s terms
    pass

# Let's use Poly to separate
from sympy import Poly as SPoly

# Express M_ab as polynomial in c
p = SPoly(M_ab, c)
print(f"  M[0,1] as polynomial in c (degree {p.degree()}):")
for power in range(p.degree() + 1):
    coeff = expand(p.nth(power))
    if coeff != 0:
        print(f"    c^{power}: {coeff}")

# ==========================================================
# Part 3: The key insight — what happens at c=0?
# ==========================================================
print(f"\n--- Part 3: Special values of c ---")

# At c=0: t_ij = s_ij (pure skew-symmetric)
M_c0 = expand(M_ab.subs(c, 0))
print(f"  M[0,1] at c=0 (pure skew): {M_c0}")

# At c=1: tournament
M_c1 = expand(M_ab.subs(c, 1))
print(f"  M[0,1] at c=1 (tournament): {M_c1}")

# At c=2:
M_c2 = expand(M_ab.subs(c, 2))
print(f"  M[0,1] at c=2: {M_c2}")

# ==========================================================
# Part 4: Parity analysis of M[a,b] under s -> -s (i.e., T -> T^op)
# ==========================================================
print(f"\n--- Part 4: Parity under s -> -s ---")
print("  Under T^op: all s_ij -> -s_ij (skew part flips).")
print("  The c/2 + s_ij -> c/2 - s_ij = 1 - (c/2 + s_ij) when c=1.")

# Substitute s_ij -> -s_ij
subs_neg = {svars[k]: -svars[k] for k in svars}
M_neg = expand(M_ab.subs(subs_neg))
print(f"\n  M[0,1](c, -s) = {M_neg}")

ratio = expand(M_ab - (-1)**(n-2) * M_neg)
print(f"  M(c,s) - (-1)^(n-2)*M(c,-s) = {ratio}")
print(f"  T^op equivalence holds? {ratio == 0}")

# Separate M into even and odd parts in s
even_part = expand((M_ab + M_neg) / 2)
odd_part = expand((M_ab - M_neg) / 2)
print(f"\n  Even part (in s): {even_part}")
print(f"  Odd part (in s): {odd_part}")
print(f"  M = even + odd: {expand(even_part + odd_part - M_ab) == 0}")

# At n=4: (-1)^{n-2} = 1, so M(c,-s) = M(c,s).
# This means M is an EVEN function of s!
# Which means: only even-degree s-monomials survive.
print(f"\n  At n={n}: (-1)^(n-2) = {(-1)**(n-2)}")
if (-1)**(n-2) == 1:
    print(f"  => M is EVEN in s: only s-degree 0 and 2 contribute")
    print(f"  => odd part should be 0: {odd_part == 0}")
else:
    print(f"  => M is ODD in s: only s-degree 1 and 3 contribute")
    print(f"  => even part should be 0: {even_part == 0}")

# ==========================================================
# Part 5: n=5 skew decomposition
# ==========================================================
print(f"\n--- Part 5: n=5 parity analysis ---")
n5 = 5
svars5 = {}
for i in range(n5):
    for j in range(i+1, n5):
        svars5[(i,j)] = symbols(f'u{i}{j}')

def t_val5(i, j):
    if i == j: return 0
    if i < j: return c/2 + svars5[(i,j)]
    return c/2 - svars5[(j,i)]

def hp5(vset, start=None, end=None):
    vl = sorted(vset)
    if len(vl) <= 1:
        if start is not None and (len(vl)==0 or vl[0] != start): return 0
        if end is not None and (len(vl)==0 or vl[0] != end): return 0
        return 1 if len(vl)==1 else 0
    total = 0
    for perm in permutations(vl):
        if start is not None and perm[0] != start: continue
        if end is not None and perm[-1] != end: continue
        prod = 1
        for ii in range(len(perm)-1):
            prod *= t_val5(perm[ii], perm[ii+1])
        total += prod
    return expand(total)

a5, b5 = 0, 1
U5 = [v for v in range(n5) if v != a5 and v != b5]
M5_ab = 0
for mask in range(1 << len(U5)):
    S = [U5[i] for i in range(len(U5)) if mask & (1 << i)]
    R = [U5[i] for i in range(len(U5)) if not (mask & (1 << i))]
    sign = (-1) ** len(S)
    ea = hp5(set(S)|{a5}, end=a5)
    bb = hp5(set(R)|{b5}, start=b5)
    M5_ab = expand(M5_ab + sign * ea * bb)

subs_neg5 = {svars5[k]: -svars5[k] for k in svars5}
M5_neg = expand(M5_ab.subs(subs_neg5))

print(f"  At n=5: (-1)^(n-2) = {(-1)**3}")
ratio5 = expand(M5_ab - (-1)**3 * M5_neg)
print(f"  M(c,s) + M(c,-s) = {ratio5}")

even5 = expand((M5_ab + M5_neg) / 2)
odd5 = expand((M5_ab - M5_neg) / 2)
print(f"  Even part = 0? {even5 == 0}")
print(f"  => M is ODD in s: {even5 == 0}")

if even5 == 0:
    print(f"\n  M[0,1] (n=5) = odd part = {odd5}")
    # As polynomial in c
    p5 = SPoly(odd5, c)
    print(f"\n  As polynomial in c (degree {p5.degree()}):")
    for power in range(p5.degree() + 1):
        coeff = expand(p5.nth(power))
        if coeff != 0:
            print(f"    c^{power}: {coeff}")

# ==========================================================
# Part 6: The determinantal perspective
# ==========================================================
print(f"\n--- Part 6: Determinantal / Pfaffian connection ---")
print("""
Key observation: the arc matrix A = (c/2)J_off + S where S is skew-symmetric.
  A + A^T = c*(J - I)  (by construction)
  A - A^T = 2*S          (skew-symmetric part)

The transfer matrix M[a,b] is symmetric iff M is even in S (n even) or odd in S (n odd).

For a SKEW-SYMMETRIC matrix S, the Pfaffian Pf(S) is a natural invariant.
But M[a,b] involves PATH PRODUCTS, not minors.

Connection: the permanent per(A) counts cycle covers (= Hamiltonian decompositions).
The determinant det(A) involves signs. The Hamiltonian path count H(T) is related
to per(A_hat) where A_hat has modified diagonal.

Question: can M[a,b] be expressed as a determinant or Pfaffian of some matrix
derived from A = (c/2)J + S?

If so, the symmetry M[a,b] = M[b,a] would follow from det/Pfaffian properties.
""")

# ==========================================================
# Part 7: Testing if M[a,b] is a polynomial in the INNER PRODUCTS of rows
# ==========================================================
print(f"\n--- Part 7: Symmetry via row inner products ---")
print("  If M[a,b] depends on rows of A only through inner products,")
print("  and A+A^T = cJ gives specific inner product structure,")
print("  then symmetry might follow.")

# At n=4 in (c,s) coordinates:
p4 = SPoly(M_ab, c)
print(f"\n  n=4: M[0,1] by c-degree:")
for power in range(p4.degree() + 1):
    coeff = expand(p4.nth(power))
    if coeff != 0:
        # Further analyze: degree in s-variables
        s_vars_list = list(svars.values())
        if s_vars_list:
            ps = SPoly(coeff, *s_vars_list)
            print(f"    c^{power} (s-degree {ps.total_degree()}): {coeff}")
        else:
            print(f"    c^{power}: {coeff}")

# ==========================================================
# Part 8: The CRUCIAL question — M at c=0
# ==========================================================
print(f"\n--- Part 8: M at c=0 (pure skew-symmetric) ---")
print("  When c=0, arcs are t_ij = s_ij, t_ji = -s_ij.")
print("  This is a 'signed tournament' where arc weights are ±w_ij.")

M_c0_explicit = expand(M_ab.subs(c, 0))
print(f"\n  n=4: M[0,1](c=0) = {M_c0_explicit}")

# At c=0, the even/odd decomposition tells us:
# n=4 (even in s): M(0,s) should be even in s
# The s-degree 0 term:
print(f"  s-degree 0 (constant): {M_c0_explicit.subs({v: 0 for v in svars.values()})}")

# Check if M(c=0) is always 0 for n even?
if M_c0_explicit == 0:
    print(f"  M(c=0) = 0 identically! Pure skew gives M=0.")
else:
    print(f"  M(c=0) nonzero.")

if even5 == 0:
    M5_c0 = expand(odd5.subs(c, 0))
    print(f"\n  n=5: M[0,1](c=0) = {M5_c0}")

print(f"\n" + "=" * 70)
print("KEY FINDINGS")
print("=" * 70)
print(f"""
1. M[a,b] = M[b,a] for ALL c-tournaments (t_ij + t_ji = c, uniform c).
2. The constraint must be UNIFORM (same c for all pairs) — non-uniform fails.
3. ALL pairs need constraining — no proper subset suffices.
4. In skew coordinates (t_ij = c/2 + s_ij, s skew-symmetric):
   - n even: M is EVEN in s (invariant under T -> T^op)
   - n odd: M is ODD in s (changes sign under T -> T^op)
5. This is exactly the T^op equivalence: M(T^op) = (-1)^(n-2) M(T).
6. The proof must leverage that A + A^T = c*(J-I) uniformly.
""")
