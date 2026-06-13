"""
Exploring the Pfaffian structure of M[a,b] in skew coordinates.

At c=0, M[a,b] is a homogeneous polynomial of degree n-2 in the
skew-symmetric variables s_ij = -s_ji.

Key question: is M[a,b](c=0) expressible as a Pfaffian, determinant,
or other algebraic object with known symmetry properties?

At n=4: M[0,1](c=0) = s02*s13 + s02*s23 + s03*s12 - s03*s23 + s12*s23 - s13*s23
At n=5: M[0,1](c=0) = 24-term cubic in u-variables
"""

from itertools import permutations
from sympy import symbols, expand, Matrix, zeros
from collections import defaultdict

print("=" * 70)
print("SKEW-PFAFFIAN TRANSFER MATRIX ANALYSIS")
print("=" * 70)

# ==========================================================
# Part 1: Can M(c=0) be expressed via the skew matrix S?
# ==========================================================
print("\n--- Part 1: M(c=0) at n=4 ---")

n = 4
svars = {}
for i in range(n):
    for j in range(i+1, n):
        svars[(i,j)] = symbols(f's{i}{j}')

def s_val(i, j):
    """Skew-symmetric: s_ij = -s_ji."""
    if i == j: return 0
    if i < j: return svars[(i,j)]
    return -svars[(j,i)]

# Build the full skew matrix S
S_mat = zeros(n, n)
for i in range(n):
    for j in range(n):
        S_mat[i,j] = s_val(i, j)

print(f"  S = {S_mat}")
print(f"  det(S) = {expand(S_mat.det())}")

# Pfaffian of S (4x4)
pf = expand(S_mat[0,1]*S_mat[2,3] - S_mat[0,2]*S_mat[1,3] + S_mat[0,3]*S_mat[1,2])
print(f"  Pf(S) = {pf}")
print(f"  Pf(S)^2 = {expand(pf**2)}")
print(f"  det(S) = Pf^2? {expand(S_mat.det() - pf**2) == 0}")

# The c=0 transfer matrix entry
# M[0,1](c=0) = s02*s13 + s02*s23 + s03*s12 - s03*s23 + s12*s23 - s13*s23
M01_c0 = expand(svars[(0,2)]*svars[(1,3)] + svars[(0,2)]*svars[(2,3)]
              + svars[(0,3)]*svars[(1,2)] - svars[(0,3)]*svars[(2,3)]
              + svars[(1,2)]*svars[(2,3)] - svars[(1,3)]*svars[(2,3)])

print(f"\n  M[0,1](c=0) = {M01_c0}")

# Check if it's related to cofactors of S
# The (i,j)-cofactor of S is (-1)^{i+j} * det(S with row i, col j removed)
# For a skew matrix, cofactors relate to Pfaffian minors

# Sub-Pfaffians: Pf_{ij} = Pfaffian of S with rows/cols i,j removed
# For 4x4: removing 2 rows/cols gives 2x2, whose Pf is just the entry

print(f"\n  Sub-Pfaffians of S:")
for i in range(n):
    for j in range(i+1, n):
        # Remove rows/cols i and j
        remaining = [k for k in range(n) if k != i and k != j]
        sub = Matrix(2, 2, lambda r, c: S_mat[remaining[r], remaining[c]])
        sub_pf = sub[0, 1]  # Pf of 2x2 skew matrix is the (0,1) entry
        print(f"    Pf_{{{i}{j}}} = {sub_pf}")

# Check: is M[0,1](c=0) a linear combination of sub-Pfaffians?
# Sub-Pfs: Pf_01 = s23, Pf_02 = -s13, Pf_03 = s12, Pf_12 = -s03, Pf_13 = s02, Pf_23 = s01
# M = s02*s13 + s02*s23 + s03*s12 - s03*s23 + s12*s23 - s13*s23

# M = s02*(s13+s23) + s03*(s12-s23) + s23*(s12-s13)
# = s02*(s13+s23) + (s12-s13)*s23 + s03*(s12-s23)
# Hmm, doesn't simplify to sub-Pfaffians easily.

# Try: M = Pf(S') where S' is S with s01 set to something?
# Pf(S) = s01*s23 - s02*s13 + s03*s12
# So Pf - s01*s23 = -s02*s13 + s03*s12
# Our M has s02*s13 + s03*s12 + more terms.
# M - Pf + s01*s23 = s02*s13 + s03*s12 + s02*s23 - s03*s23 + s12*s23 - s13*s23
#                   - (-s02*s13 + s03*s12)
#                   = 2*s02*s13 + s02*s23 - s03*s23 + s12*s23 - s13*s23

# No clean relationship. Let me try something different.

# ==========================================================
# Part 2: Direct computation — what IS M[a,b](c=0)?
# ==========================================================
print(f"\n--- Part 2: Direct computation of M at c=0 ---")

def hp_skew(vertex_set, start=None, end=None):
    """Path count with skew arc weights."""
    vl = sorted(vertex_set)
    if len(vl) <= 1:
        if start is not None and (len(vl)==0 or vl[0] != start): return 0
        if end is not None and (len(vl)==0 or vl[0] != end): return 0
        return 1 if len(vl) == 1 else 0
    total = 0
    for perm in permutations(vl):
        if start is not None and perm[0] != start: continue
        if end is not None and perm[-1] != end: continue
        prod = 1
        for i in range(len(perm)-1):
            prod *= s_val(perm[i], perm[i+1])
        total += prod
    return expand(total)

a, b = 0, 1
U = [v for v in range(n) if v != a and v != b]

# Compute each E and B at c=0
print(f"  Path counts with skew weights:")
for mask in range(1 << len(U)):
    S = [U[i] for i in range(len(U)) if mask & (1 << i)]
    R = [U[i] for i in range(len(U)) if not (mask & (1 << i))]

    ea = hp_skew(set(S)|{a}, end=a)
    bb = hp_skew(set(R)|{b}, start=b)

    print(f"    S={str(S):10s}: E_0(S+0)={str(ea):20s}, B_1(R+1)={bb}")

# ==========================================================
# Part 3: The row/column structure of S
# ==========================================================
print(f"\n--- Part 3: Row/column sums and products ---")

# For a c-tournament at c=0: each row of A sums to 0 (since sum_j t_ij = sum_j s_ij,
# and sum_j s_ij is NOT zero in general for a skew matrix).
# But for the COMPLETE skew matrix, row sums are zero by antisymmetry
# (sum_{j!=i} s_ij = -sum_{j!=i} s_ji, and if we sum BOTH we get 0).
# But individual row sums are NOT zero.

for i in range(n):
    row_sum = sum(s_val(i, j) for j in range(n) if j != i)
    print(f"  Row {i} sum: {expand(row_sum)}")

# ==========================================================
# Part 4: M(c=0) as a sum over "2-path covers" with skew weights
# ==========================================================
print(f"\n--- Part 4: 2-path cover interpretation ---")
print("""
M[a,b](c=0) = sum_S (-1)^|S| * (product of s_ij along path ending at a in S+a)
                               * (product of s_ij along path starting from b in R+b)

Each term is a product of n-2 skew variables (one per arc used).
The total degree is n-2 (confirmed: degree 2 at n=4, degree 3 at n=5).

Since each s_ij appears at most once (paths are vertex-disjoint from each other
and Hamiltonian within their part), these are SQUARE-FREE monomials.
""")

# Let's enumerate all monomials at n=4
print(f"  Monomials in M[0,1](c=0):")
poly = M01_c0.as_coefficients_dict()
for monom, coeff in sorted(poly.items(), key=lambda x: str(x[0])):
    if coeff != 0:
        print(f"    {'+' if int(coeff) > 0 else ''}{int(coeff)} * {monom}")

# ==========================================================
# Part 5: n=5 at c=0 — look for patterns
# ==========================================================
print(f"\n--- Part 5: n=5 at c=0 ---")

n5 = 5
sv5 = {}
for i in range(n5):
    for j in range(i+1, n5):
        sv5[(i,j)] = symbols(f'u{i}{j}')

def s5(i, j):
    if i == j: return 0
    if i < j: return sv5[(i,j)]
    return -sv5[(j,i)]

def hp5_skew(vset, start=None, end=None):
    vl = sorted(vset)
    if len(vl) <= 1:
        if start is not None and (len(vl)==0 or vl[0] != start): return 0
        if end is not None and (len(vl)==0 or vl[0] != end): return 0
        return 1 if len(vl) == 1 else 0
    total = 0
    for perm in permutations(vl):
        if start is not None and perm[0] != start: continue
        if end is not None and perm[-1] != end: continue
        prod = 1
        for ii in range(len(perm)-1):
            prod *= s5(perm[ii], perm[ii+1])
        total += prod
    return expand(total)

a5, b5 = 0, 1
U5 = [v for v in range(n5) if v != a5 and v != b5]
M5_c0 = 0
for mask in range(1 << len(U5)):
    S = [U5[i] for i in range(len(U5)) if mask & (1 << i)]
    R = [U5[i] for i in range(len(U5)) if not (mask & (1 << i))]
    sign = (-1) ** len(S)
    ea = hp5_skew(set(S)|{a5}, end=a5)
    bb = hp5_skew(set(R)|{b5}, start=b5)
    M5_c0 = expand(M5_c0 + sign * ea * bb)

print(f"  M[0,1](c=0, n=5) has {len(M5_c0.as_ordered_terms())} terms")

# Group by which variables appear
poly5 = M5_c0.as_coefficients_dict()
# Check: each monomial is a product of 3 s-variables
degrees = set()
for monom, coeff in poly5.items():
    if coeff != 0 and monom != 1:
        d = sum(1 for v in sv5.values() if monom.has(v))
        degrees.add(d)
print(f"  Degrees of monomials: {degrees}")

# Count +1 and -1 coefficients
plus_count = sum(1 for c in poly5.values() if c == 1)
minus_count = sum(1 for c in poly5.values() if c == -1)
print(f"  Coefficients: {plus_count} terms with +1, {minus_count} terms with -1")

# ==========================================================
# Part 6: Symmetry under a<->b swap at c=0
# ==========================================================
print(f"\n--- Part 6: Symmetry check of M(c=0) ---")

# At n=4: swap a=0 <-> b=1 means s02<->s12, s03<->s13, s23 stays
# AND signs: s_02 -> s_12, s_12 -> s_02, etc. (just relabeling)
# s_20 = -s_02 -> s_21 = -s_12, etc.

# Check: is M[0,1](c=0) symmetric under 0<->1?
subs_swap4 = {
    svars[(0,2)]: svars[(1,2)],
    svars[(0,3)]: svars[(1,3)],
    svars[(1,2)]: svars[(0,2)],
    svars[(1,3)]: svars[(0,3)],
    # s23 stays
}
M01_swapped = expand(M01_c0.subs(subs_swap4))
print(f"  n=4: M[0,1](c=0) under 0<->1 swap: {M01_swapped}")
print(f"  Equal to M[0,1]? {M01_swapped == M01_c0}")
print(f"  => M[0,1](c=0) IS invariant under swapping a and b (= symmetry)")

# At n=5:
subs_swap5 = {
    sv5[(0,2)]: sv5[(1,2)],
    sv5[(0,3)]: sv5[(1,3)],
    sv5[(0,4)]: sv5[(1,4)],
    sv5[(1,2)]: sv5[(0,2)],
    sv5[(1,3)]: sv5[(0,3)],
    sv5[(1,4)]: sv5[(0,4)],
}
M5_swapped = expand(M5_c0.subs(subs_swap5))
# For odd n, M is ODD in s, so under swap we expect M -> M (symmetry)
# But wait: the swap doesn't flip signs. It just relabels.
# M[0,1] = M[1,0] means M(s) with 0<->1 = M(s) with 0<->1.
# But M[1,0] is computed with roles reversed.
# Actually M[0,1] being symmetric under 0<->1 IS M[0,1] = M[1,0].

# But we need to be careful: swapping 0<->1 in the formula for M[0,1]
# gives M[1,0] (because the same subsets S of U are summed over,
# but now E_1(S+1) and B_0(R+0) are used).

# The substitution is the right check.
print(f"\n  n=5: M[0,1](c=0) symmetric under 0<->1? {expand(M5_swapped - M5_c0) == 0}")

# ==========================================================
# Part 7: Can M(c=0) be expressed as a DETERMINANT?
# ==========================================================
print(f"\n--- Part 7: Determinantal test at n=4 ---")

# At n=4, M[0,1](c=0) is a degree-2 polynomial in 5 variables.
# A 2x2 determinant of linear forms would give exactly degree 2.
# det([[a1, a2], [a3, a4]]) = a1*a4 - a2*a3

# We need to find linear forms L1, L2, L3, L4 in {s02, s03, s12, s13, s23}
# such that L1*L4 - L2*L3 = M[0,1](c=0).

# M = s02*s13 + s02*s23 + s03*s12 - s03*s23 + s12*s23 - s13*s23

# Try: factor M as a product or determinant.
# Let's check if it factors.
from sympy import factor
M_factored = factor(M01_c0)
print(f"  factor(M[0,1](c=0)) = {M_factored}")

# Try collecting by s23:
# M = s23*(s02 - s03 + s12 - s13) + s02*s13 + s03*s12
# = s23*(s02 + s12 - s03 - s13) + (s02*s13 + s03*s12)
# The first part: s23 * [(s02+s12) - (s03+s13)]
# The second part: s02*s13 + s03*s12 (symmetric in 0<->1!)

# Note: s02*s13 + s03*s12 = partial Pfaffian (sum of two matching terms)
# And s23*(...) involves the U-internal edge

# Try as det: det([[s02+s12, s03+s13], [-1, s23]]) = (s02+s12)*s23 + (s03+s13)
# No, that gives wrong terms.

# Try det([[s02, -s03], [s12-s23, s13+s23]])
# = s02*(s13+s23) - (-s03)*(s12-s23) = s02*s13 + s02*s23 + s03*s12 - s03*s23
# = M - s12*s23 + s13*s23
# Close but not quite!

# det([[s02, -(s03-s23)], [s12-s23, s13]])
# = s02*s13 + (s03-s23)*(s12-s23)
# = s02*s13 + s03*s12 - s03*s23 - s12*s23 + s23^2
# No, has s23^2.

# Actually, let me try: det([[s02+s23, -s03+s23], [-s13+s23, s12+s23]])
# Hmm this is getting complicated. Let me just check if it's even possible
# by checking the "discriminant" of the quadratic form.

# M as bilinear form in two groups:
# Group A = {s02, s03}: edges from vertex 0 to U
# Group B = {s12, s13}: edges from vertex 1 to U
# "Internal" s23: edge within U

# M = s02*s13 + s03*s12 + s23*(s02-s03+s12-s13)
# = s02*(s13+s23) + s03*(s12-s23) + s23*(s12-s13)
# = s02*(s13+s23) + (s12-s13)*s23 + s03*(s12-s23)

# Alternative grouping:
# M = (s02+s23)*(s13+s23) - s23*(s13+s23) + s03*(s12-s23) + s23*(s12-s13)
# This is messy.

# Let me try a 3x3 permanent or something.
# Actually, for n=4, |U|=2, and the paths have length 0, 1, or 2.
# The path of length 0 is just {a} (trivial).
# Maybe M has a MATRIX structure indexed by U.

# Define P[u,v] = path-product of u->v (if exists) in the skew graph on U
# For |U|=2: P[2,3] = s23, P[3,2] = -s23

# E_0({0,u}) = s_u0 = -s_{0u} (path u->0)
# B_1({1,v}) = s_1v (path 1->v)
# E_0({0,2,3}) = s_20*s_32 + s_30*s_23 = -s02*(-s23) + (-s03)*s23 = s02*s23 - s03*s23

# Hmm, let me think about this as a matrix product.
# For each u in U, define:
#   e_a(u) = s_{ua} = -s_{au} (1-step path: u -> a)
#   b_b(u) = s_{bu} (1-step path: b -> u)

print(f"\n  Edge vectors:")
for u in U:
    ea_u = s_val(u, a)
    bb_u = s_val(b, u)
    print(f"    u={u}: e_a(u) = s_{u}{a} = {ea_u}, b_b(u) = s_{b}{u} = {bb_u}")

# The full paths E_0({0,2,3}):
# Orderings of {0,2,3} ending at 0: (2,3,0) and (3,2,0)
# (2,3,0): s23 * s30 = s23 * (-s03)
# (3,2,0): s32 * s20 = (-s23) * (-s02) = s23 * s02
# Total: s23*(-s03) + s23*s02 = s23*(s02 - s03)
e023 = hp_skew({0,2,3}, end=0)
print(f"\n  E_0({{0,2,3}}) = {e023}")
print(f"  = s23 * (s02 - s03) ? {expand(e023 - svars[(2,3)]*(svars[(0,2)] - svars[(0,3)])) == 0}")

# Similarly:
b123 = hp_skew({1,2,3}, start=1)
print(f"  B_1({{1,2,3}}) = {b123}")

# The key formula: at n=4, M[0,1](c=0) =
# S=empty: 1 * B_1({1,2,3})
# S={2}: -1 * E_0({0,2}) * B_1({1,3})
# S={3}: -1 * E_0({0,3}) * B_1({1,2})
# S={2,3}: +1 * E_0({0,2,3}) * 1
# = B_1({1,2,3}) - E_0({0,2})*B_1({1,3}) - E_0({0,3})*B_1({1,2}) + E_0({0,2,3})

# ==========================================================
# Part 8: The SIGNED PERMANENT interpretation
# ==========================================================
print(f"\n--- Part 8: Signed permanent interpretation ---")
print("""
At c=0, every path of length k through a set of vertices has weight
that is a product of k skew variables, with each variable's sign
determined by the direction of traversal.

For the 2-path cover (path1 through S+a ending at a, path2 through R+b
starting from b), the weight is the product of all arc weights along both paths.

The key: the SAME SET of n-2 edges is used, but the SIGNS depend on
which direction each edge is traversed. The alternating sum over S
selects specific combinations of edge directions.
""")

# ==========================================================
# Part 9: Is M(c=0) the (a,b)-cofactor of something?
# ==========================================================
print(f"\n--- Part 9: Cofactor interpretation ---")

# For a skew matrix S of order n (even), det(S) = Pf(S)^2.
# The cofactor matrix cof(S) has cof_ij = (-1)^{i+j} det(S_ij).
# For skew S, the cofactor matrix has special structure.

# Does M[0,1](c=0) equal some cofactor of S?
S4 = S_mat  # our 4x4 skew matrix

# (0,1)-cofactor: (-1)^{0+1} * det(S with row 0, col 1 removed)
S_01_minor = S4.minor_submatrix(0, 1)
cof_01 = -1 * S_01_minor.det()
print(f"  Cofactor(S, 0, 1) = {expand(cof_01)}")

# (1,0)-cofactor:
S_10_minor = S4.minor_submatrix(1, 0)
cof_10 = S_10_minor.det()
print(f"  Cofactor(S, 1, 0) = {expand(cof_10)}")

# Compare with M[0,1](c=0)
print(f"  M[0,1](c=0) = {M01_c0}")
print(f"  Equal to cof(0,1)? {expand(cof_01 - M01_c0) == 0}")
print(f"  Equal to -cof(0,1)? {expand(-cof_01 - M01_c0) == 0}")

# Maybe a different matrix? Try A = J_off + S (with c=2, so t_ij = 1 + s_ij)
A_mat = zeros(n, n)
for i in range(n):
    for j in range(n):
        if i != j:
            A_mat[i,j] = 1 + s_val(i, j)

print(f"\n  A = J_off + S:")
for i in range(n):
    print(f"    row {i}: {[A_mat[i,j] for j in range(n)]}")
print(f"  det(A) = {expand(A_mat.det())}")

# M at c=2
M01_c2 = expand(M01_c0 + 2)  # from earlier: M = c^2/2 + Q(s), so at c=2: Q(s)+2
print(f"  M[0,1](c=2) = {M01_c2}")

# Cofactor of A at (0,1):
A_01_minor = A_mat.minor_submatrix(0, 1)
cof_A_01 = -1 * A_01_minor.det()
print(f"  Cofactor(A, 0, 1) = {expand(cof_A_01)}")
print(f"  Equal to M(c=2)? {expand(cof_A_01 - M01_c2) == 0}")

# What about cofactor / (n-2)! ?
import math
print(f"  Cofactor/(n-2)! = {expand(cof_A_01 / math.factorial(n-2))}")

# ==========================================================
# Part 10: The INVERSE matrix connection
# ==========================================================
print(f"\n--- Part 10: Inverse/adjugate connection ---")

# adj(A) = cofactor matrix^T. For a matrix with A+A^T = c*(J-I):
# If M[a,b] = adj(A)[b,a] / (n-2)!, that would make M symmetric!
# Because adj(A) is always symmetric when A+A^T is symmetric.
# Actually adj(A)^T_ij = cof(A)_ij, so adj(A)_ij = cof(A)_ji.

# If A + A^T = c*(J-I), then adj(A) + adj(A^T) = ?
# For any matrix, adj(A^T) = adj(A)^T.
# So adj(A) + adj(A)^T comes from A + A^T = c*(J-I).

# Let's just test: adj(A)[1,0] vs M[0,1]
adj_A = A_mat.adjugate()
print(f"  adj(A)[1,0] = {expand(adj_A[1,0])}")
print(f"  adj(A)[0,1] = {expand(adj_A[0,1])}")
print(f"  M[0,1](c=2) = {M01_c2}")
print(f"  adj(A)[1,0] = M(c=2)? {expand(adj_A[1,0] - M01_c2) == 0}")
print(f"  adj(A)[0,1] = M(c=2)? {expand(adj_A[0,1] - M01_c2) == 0}")
print(f"  adj(A)[1,0]/(n-2)! = {expand(adj_A[1,0] / math.factorial(n-2))}")

# Also try with the matrix (c/2)*J_off + S (general c)
c_sym = symbols('c')
A_gen = zeros(n, n)
for i in range(n):
    for j in range(n):
        if i != j:
            A_gen[i,j] = c_sym/2 + s_val(i, j)

det_gen = expand(A_gen.det())
print(f"\n  det(A_gen) with c/2+s = {det_gen}")

adj_gen = A_gen.adjugate()
adj_10 = expand(adj_gen[1,0])
adj_01 = expand(adj_gen[0,1])
print(f"  adj(A_gen)[1,0] = {adj_10}")
print(f"  adj(A_gen)[0,1] = {adj_01}")
print(f"  adj[1,0] = adj[0,1]? {expand(adj_10 - adj_01) == 0}")

# M at general c = c^2/2 + Q(s)
M01_gen = expand(c_sym**2 / 2 + M01_c0)
print(f"\n  M[0,1](c,s) = {M01_gen}")
print(f"  adj(A_gen)[1,0] / (n-2)! = {expand(adj_10 / math.factorial(n-2))}")
print(f"  Equal to M? {expand(adj_10 / math.factorial(n-2) - M01_gen) == 0}")

print(f"\n" + "=" * 70)
print("FINDINGS")
print("=" * 70)
