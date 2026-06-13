#!/usr/bin/env python3
"""
Search for a matrix identity expressing M[a,b] in a manifestly symmetric form.

Key insight from audit: M[a,b](-r,s) = M[b,a](r,s) is PROVED.
So symmetry <=> M is even in r.

Can we find an expression for M that is MANIFESTLY even in r?

IDEA 1: M might be a determinant of a matrix whose entries are even in r.
IDEA 2: M might be a Pfaffian of a skew matrix with entries even in r.
IDEA 3: M might be expressible as f(r^2, s) for some function f.
IDEA 4: M might relate to eigenvalues of A = r(J-I) + S.

opus-2026-03-06-S21
"""

from itertools import permutations, combinations
from sympy import (symbols, expand, Poly, Matrix, eye, ones, det,
                   sqrt, Rational, factor, collect, simplify)

def setup(n):
    r = symbols('r')
    sv = {}
    for i in range(n):
        for j in range(i+1, n):
            sv[(i,j)] = symbols(f's{i}{j}')
    def s(i, j):
        if i == j: return 0
        if i < j: return sv[(i,j)]
        return -sv[(j,i)]
    def t(i, j):
        if i == j: return 0
        return r + s(i, j)
    return r, sv, s, t

def hp(t_fn, vset, start=None, end=None):
    vl = sorted(vset)
    k = len(vl)
    if k == 0: return 0
    if k == 1:
        if start is not None and vl[0] != start: return 0
        if end is not None and vl[0] != end: return 0
        return 1
    total = 0
    for p in permutations(vl):
        if start is not None and p[0] != start: continue
        if end is not None and p[-1] != end: continue
        prod = 1
        for i in range(len(p)-1):
            prod *= t_fn(p[i], p[i+1])
        total += prod
    return expand(total)

def transfer_M(t_fn, n, a, b):
    U = [v for v in range(n) if v != a and v != b]
    result = 0
    for mask in range(1 << len(U)):
        S = [U[i] for i in range(len(U)) if mask & (1 << i)]
        R = [U[i] for i in range(len(U)) if not (mask & (1 << i))]
        sign = (-1)**len(S)
        ea = hp(t_fn, set(S)|{a}, end=a)
        bb = hp(t_fn, set(R)|{b}, start=b)
        result += sign * ea * bb
    return expand(result)

print("=" * 70)
print("SEARCHING FOR MATRIX IDENTITY")
print("=" * 70)

# ============================================================
# Part 1: Try M = det(X) or M = Pf(X) for various X
# ============================================================
print("\n--- Part 1: Determinantal / Pfaffian search at n=4 ---")

n = 4
r, sv, s, t = setup(n)
a, b = 0, 1
M01 = transfer_M(t, n, a, b)
print(f"  M[0,1] = {M01}")
print(f"  M[0,1] = 2*r^2 + Q where Q = {expand(M01 - 2*r**2)}")

# Is M[0,1] a 2x2 determinant?
# det([[a,b],[c,d]]) = ad - bc
# Need: ad - bc = 2r^2 + s02*s13 + s02*s23 + s03*s12 - s03*s23 + s12*s23 - s13*s23

# Try: a = r + alpha, d = r + delta, b = beta, c = gamma
# Then ad - bc = r^2 + (alpha+delta)*r + alpha*delta - beta*gamma
# For r^1 to vanish: alpha + delta = 0, so delta = -alpha.
# Then ad - bc = r^2 - alpha^2 - beta*gamma
# And the r^0 part: -alpha^2 - beta*gamma = Q/2 (since overall factor is 2)

# Actually M = 2r^2 + Q. So M/2 = r^2 + Q/2.
# As a 2x2 det: r^2 - alpha^2 - beta*gamma = r^2 + Q/2
# So alpha^2 + beta*gamma = -Q/2

Q = expand(M01 - 2*r**2)
print(f"\n  Q = {Q}")
print(f"  -Q/2 = {expand(-Q/2)}")

# -Q/2 = (-s02*s13 - s02*s23 - s03*s12 + s03*s23 - s12*s23 + s13*s23) / 2
# Can this be written as alpha^2 + beta*gamma?

# If alpha = 0: need beta*gamma = -Q/2. So M/2 = r^2 - beta*gamma.
# That would mean M/2 = det([[r, beta], [gamma, r]]).
# And M = 2*det([[r, beta], [gamma, r]]) = 2*(r^2 - beta*gamma).
# So -beta*gamma = Q/2.

# Let's try: beta = s_{something}, gamma = s_{something else}

# Actually, since Q has 6 terms and is degree 2, beta*gamma must be
# a product of two degree-1 polynomials in s. Let's try:
# beta = a1*s02 + a2*s03 + a3*s12 + a4*s13 + a5*s23
# gamma = b1*s02 + b2*s03 + b3*s12 + b4*s13 + b5*s23

# beta*gamma should equal -Q/2. This is an underdetermined system.
# Let's just check if -Q/2 factors as a product of two linear forms.

from sympy import factor
Q_half = expand(-Q/2)
Q_factored = factor(Q_half)
print(f"  factor(-Q/2) = {Q_factored}")

# If it doesn't factor, try with a nonzero alpha
# M/2 = r^2 + Q/2 = det([[r+alpha, beta], [-gamma, r-alpha]])
# = (r+alpha)(r-alpha) + beta*gamma = r^2 - alpha^2 + beta*gamma
# So -alpha^2 + beta*gamma = Q/2

# What if alpha is a linear polynomial in s?
# Then alpha^2 is degree 2, beta*gamma is degree 2.
# We need: beta*gamma - alpha^2 = Q/2.

# ============================================================
# Part 2: Try Schur complement / block determinant
# ============================================================
print("\n--- Part 2: Block determinant approach ---")
print("  Can M be expressed as a block determinant?")

# det([[A11, A12], [A21, A22]]) = det(A11) * det(A22 - A21 * A11^{-1} * A12)
# If we can write M as a Schur complement...

# Actually let me try: does M relate to the (a,b) entry of adj(rI - S)?
# where S is the skew matrix

S_mat = Matrix(n, n, lambda i, j: s(i,j) if i != j else 0)
print(f"\n  Skew matrix S:\n{S_mat}")

rI_minus_S = r * eye(n) - S_mat
adj_rS = rI_minus_S.adjugate()

print(f"\n  adj(rI - S)[0,1] = {expand(adj_rS[0,1])}")
print(f"  M[0,1] = {M01}")
print(f"  Equal? {expand(adj_rS[0,1] - M01) == 0}")

# Try (rI + S)
rI_plus_S = r * eye(n) + S_mat
adj_rS_plus = rI_plus_S.adjugate()
print(f"\n  adj(rI + S)[0,1] = {expand(adj_rS_plus[0,1])}")
print(f"  Equal to M? {expand(adj_rS_plus[0,1] - M01) == 0}")

# Try r(J-I) - S
J = ones(n,n)
I_n = eye(n)
A_mat = r * (J - I_n) + S_mat
IminA = I_n - A_mat
adj_IminA = IminA.adjugate()
print(f"\n  adj(I - A)[0,1] = {expand(adj_IminA[0,1])}")
# This involves s01 which doesn't appear in M. Expected.

# The issue: M uses only arcs WITHIN S+a and R+b, not arc (a,b).
# So s01 = s(0,1) should NOT appear in M. Let's verify:
from sympy import Symbol
s01_sym = sv.get((0,1), None)
if s01_sym is not None:
    M01_has_s01 = M01.coeff(s01_sym)
    print(f"\n  Does M[0,1] contain s01? {M01_has_s01 != 0}")
    print(f"  Coefficient of s01 in M[0,1]: {M01_has_s01}")

# ============================================================
# Part 3: The Laplacian approach
# ============================================================
print("\n--- Part 3: Laplacian / Matrix-Tree connection ---")

# For A = r(J-I) + S, the Laplacian is L = D - A where D = diag(row sums).
# Row sum of A_i = sum_{j≠i} (r + s(i,j)) = (n-1)r + sigma_i
# where sigma_i = sum_{j≠i} s(i,j)

sigma = [sum(s(i,j) for j in range(n) if j != i) for i in range(n)]
print(f"  Row sums sigma_i: {sigma}")

D_mat = Matrix(n, n, lambda i, j: (n-1)*r + sigma[i] if i == j else 0)
L = D_mat - A_mat
print(f"\n  Laplacian L = D - A:")
for i in range(n):
    print(f"    row {i}: {[expand(L[i,j]) for j in range(n)]}")

# Matrix-Tree theorem: all cofactors of L are equal (for strongly connected)
# cofactor(L, i, j) = (-1)^{i+j} * det(L with row i, col j deleted)
# For a tournament, this counts spanning arborescences rooted at j.

# But M counts something related to Hamiltonian paths, not spanning trees.
# The relationship is: M[a,b] = sum_{arborescences rooted at a with...} ??
# Not obvious.

# Let's compute the (0,1) cofactor of L
minor_01 = L.minor_submatrix(0, 1)
cof_01 = (-1)**(0+1) * det(minor_01)
cof_01 = expand(cof_01)
print(f"\n  cofactor(L, 0, 1) = {cof_01}")
print(f"  M[0,1] = {M01}")
print(f"  Equal? {expand(cof_01 - M01) == 0}")
ratio = None
if M01 != 0:
    # Check if it's a constant multiple
    p_cof = Poly(cof_01, r)
    p_M = Poly(M01, r)
    if p_M.degree() > 0:
        ratio = expand(p_cof.nth(p_cof.degree()) / p_M.nth(p_M.degree()))
        print(f"  Leading coeff ratio: {ratio}")

# ============================================================
# Part 4: Characteristic polynomial approach
# ============================================================
print("\n--- Part 4: Characteristic polynomial of S ---")

# S is skew-symmetric. Its characteristic polynomial has the form:
# det(lambda*I - S) = lambda^n + c_{n-2} lambda^{n-2} + c_{n-4} lambda^{n-4} + ...
# (only even powers if n even, only odd+even if n odd)

char_S = expand(det(symbols('lam') * I_n - S_mat))
lam = symbols('lam')
char_S = expand(det(lam * I_n - S_mat))
p_char = Poly(char_S, lam)
print(f"  char poly of S: degree {p_char.degree()}")
for k in range(p_char.degree() + 1):
    c_k = expand(p_char.nth(k))
    if c_k != 0:
        print(f"    lambda^{k}: {c_k}")

# For n=4 skew: det(lam*I - S) = lam^4 + c2*lam^2 + c0
# c2 = sum of 2x2 principal Pfaffians squared
# c0 = Pf(S)^2

# The coefficients c2, c0 are EVEN functions of S (invariant under S -> -S)
# because S and -S are conjugate by a permutation matrix... no, that's not right.
# Actually det(lam*I - S) = det(lam*I + S^T) = det(lam*I + (-S)) = char poly of -S.
# Hmm, that gives the SAME polynomial because eigenvalues of -S are negatives of eigs of S,
# and char poly of -S at lam is char poly of S at -lam.
# So det(lam*I - S) at -lam = det(-lam*I - S) = (-1)^n det(lam*I + S) = (-1)^n det(lam*I - (-S)).
# For S skew, -S = S^T, and det(lam*I - S^T) = det(lam*I - S)^T = det(lam*I - S).
# So: char_S(-lam) = (-1)^n char_S(lam). For n=4: char(-lam) = char(lam), all even.

# ============================================================
# Part 5: The key question — can M be a polynomial in r^2 and
# the coefficients of char(S)?
# ============================================================
print("\n--- Part 5: M in terms of r^2 and invariants of S ---")

# At n=4: M = 2r^2 + Q. And char(S) = lam^4 + c2*lam^2 + c0.
# c2 and c0 are symmetric functions of the s-variables.

c2_val = expand(p_char.nth(2))
c0_val = expand(p_char.nth(0))
print(f"  c2 (of char poly) = {c2_val}")
print(f"  c0 (of char poly) = {c0_val}")
print(f"  Q (of M) = {Q}")

# Is Q a linear combination of c2 and c0?
# Q = s02*s13 + s02*s23 + s03*s12 - s03*s23 + s12*s23 - s13*s23
# c2 = s02^2 + s03^2 + s12^2 + s13^2 + s23^2 + ... (sum of squares of principal 2x2 pfaffians)

# c2 is a sum of SQUARES, Q has cross-terms. These are different.
# So M is NOT a simple function of the characteristic polynomial coefficients.

# But maybe M is a coefficient of a MODIFIED characteristic polynomial?

# ============================================================
# Part 6: The (a,b)-restricted characteristic polynomial
# ============================================================
print("\n--- Part 6: Restricted characteristic polynomial ---")

# What if we consider the characteristic polynomial of the
# submatrix of S obtained by deleting rows/cols a,b?
# Or the "bordered" matrix?

# Delete row a, col b from (rI - S):
restricted = (r * I_n - S_mat).minor_submatrix(a, b)
det_restricted = expand(det(restricted))
print(f"  det(rI-S with row {a}, col {b} deleted) = {det_restricted}")
# This is the (a,b) cofactor of (rI-S), which is adj(rI-S)[b,a].

# Hmm, we already checked adj(rI-S). Let's look at it differently.
# adj(rI-S)[b,a] = (-1)^{a+b} det(minor) ... same thing.

# What about the (a,b) entry of (rI-S)^{-1} * det(rI-S)?
# That IS the adjugate, which we've checked.

# ============================================================
# Part 7: DIRECT r^2 substitution
# ============================================================
print("\n--- Part 7: Is M expressible as f(r^2, s)? ---")

# Since M has only even powers of r, we should be able to write
# M(r,s) = g(r^2, s) for some polynomial g.

# At n=4: M = 2r^2 + Q. So g(u,s) = 2u + Q where u = r^2.

# At n=5: let's compute
n5 = 5
r5, sv5, s5, t5 = setup(n5)
M01_5 = transfer_M(t5, n5, 0, 1)
p5 = Poly(M01_5, r5)
print(f"\n  n=5: M[0,1] as polynomial in r:")
for k in range(p5.degree() + 1):
    ck = expand(p5.nth(k))
    if ck != 0:
        print(f"    r^{k}: {ck}")

# At n=5: M = r^2 * Q2 + Q0 (only even powers)
Q2_5 = expand(p5.nth(2))
Q0_5 = expand(p5.nth(0))
print(f"\n  n=5: M = r^2 * ({Q2_5}) + ({Q0_5})")

# So g(u,s) = u * Q2 + Q0 for n=5.

# ============================================================
# Part 8: NEW IDEA — the "signed permanent" connection
# ============================================================
print("\n--- Part 8: Signed permanent / immanant ---")
print("""
The transfer matrix M[a,b] = sum_S (-1)^|S| E_a B_b involves a sum
over all 2-path-covers (partitions of V into two directed paths).

This is reminiscent of the PERMANENT of a matrix:
per(A) = sum over permutations sigma, product A[i, sigma(i)]

A permutation decomposes into cycles; a 2-PATH-COVER is not a
permutation but rather two paths covering all vertices.

However, the HAMILTONIAN PATH count itself is:
HP(a,b) = sum over Hamiltonian paths a -> ... -> b of product of weights
         = sum over certain permutations of products of A entries

The transfer matrix combines Hamiltonian paths with inclusion-exclusion.
""")

# ============================================================
# Part 9: Is M related to a HAFNIAN or permanent?
# ============================================================
print("\n--- Part 9: Hafnian / permanent connection ---")

# The hafnian of a symmetric matrix B is:
# haf(B) = sum over perfect matchings M, product B[i,j] for {i,j} in M

# A 2x2 "hafnian" is just B[0,1].
# A 4x4 hafnian has 3 terms.

# Our Q at n=4 has 6 terms. Not a hafnian.

# What about: M = permanent of some matrix?
# per([[a,b],[c,d]]) = ad + bc.
# This has 2 terms. M at n=3 has 2 terms: s02 + s12. Hmm!

# n=3: M[0,1] = s02 + s12 = per([[s02, 0], [0, s12]])? No, per of diagonal is product.
# per([[0, s02], [s12, 0]]) = 0*0 + s02*s12. Not right.
# per([[1, 1], [s02, s12]]) = 1*s12 + 1*s02 = s02 + s12. YES!

# So at n=3, M[0,1] = per([[1, 1], [s02, s12]]).
# This is a 2x1 permanent... or rather a 2x2 permanent with a row of 1s.

# At n=4: M = 2r^2 + Q. Not obviously a permanent of a small matrix.

# ============================================================
# Part 10: Summary of n=3 structure
# ============================================================
print("\n--- Part 10: n=3 deep look ---")
n3 = 3
r3, sv3, s3, t3 = setup(n3)
M01_3 = transfer_M(t3, n3, 0, 1)
M10_3 = transfer_M(t3, n3, 1, 0)
M02_3 = transfer_M(t3, n3, 0, 2)
M20_3 = transfer_M(t3, n3, 2, 0)
print(f"  M[0,1] = {M01_3}")
print(f"  M[1,0] = {M10_3}")
print(f"  M[0,2] = {M02_3}")
print(f"  M[2,0] = {M20_3}")

# M[0,1] = s02 + s12. This is the sum of s_v2 for v in {0,1}.
# M[0,2] = s02 - s12. Wait, is it? Let me check:
# Actually with just s02 and s12 (since n=3, there's also s01).
# But s01 shouldn't appear... hmm.

# At n=3, U = {2}.
# M[0,1]: S=empty, R={2}.
#   E_0({0}, end=0) = 1, B_1({2,1}, start=1) = t(1,2) = r + s12.
#   S={2}: sign=-1. E_0({2,0}, end=0) = t(2,0) = r + s(2,0) = r - s02.
#   B_1({1}, start=1) = 1.
#   M[0,1] = 1*(r + s12) + (-1)*(r - s02)*1 = r + s12 - r + s02 = s02 + s12.
# The r cancels! And s01 doesn't appear. Good.

# M[0,2]: U = {1}.
#   S=empty, R={1}. E_0({0}, end=0) = 1. B_2({1,2}, start=2) = t(2,1) = r + s(2,1) = r - s12.
#   S={1}: sign=-1. E_0({1,0}, end=0) = t(1,0) = r + s(1,0) = r - s01.
#   B_2({2}, start=2) = 1.
#   M[0,2] = 1*(r - s12) + (-1)*(r - s01)*1 = r - s12 - r + s01 = s01 - s12.

print(f"\n  M[0,1] = s02 + s12 (doesn't depend on r or s01)")
print(f"  M[0,2] = s01 - s12 (doesn't depend on r or s02)")
print(f"  M is degree 0 in r at n=3 (trivially even)")

# Entire transfer matrix at n=3:
print(f"\n  Full transfer matrix at n=3:")
for aa in range(n3):
    for bb in range(n3):
        if aa == bb: continue
        M_val = transfer_M(t3, n3, aa, bb)
        print(f"    M[{aa},{bb}] = {M_val}")

# ============================================================
# Part 11: n=4 full matrix
# ============================================================
print(f"\n--- Part 11: Full transfer matrix at n=4 ---")
n4 = 4
r4, sv4, s4, t4 = setup(n4)
print(f"  Computing all M[a,b] for n=4...")
for aa in range(n4):
    for bb in range(n4):
        if aa == bb: continue
        M_val = transfer_M(t4, n4, aa, bb)
        p = Poly(M_val, r4)
        r2_coeff = expand(p.nth(2))
        r0_coeff = expand(p.nth(0))
        r1_coeff = expand(p.nth(1))
        print(f"    M[{aa},{bb}] = {r2_coeff}*r^2 + {r1_coeff}*r + {r0_coeff}")

print("\n" + "=" * 70)
print("FINDINGS")
print("=" * 70)
print("""
1. M is NOT a simple determinant, cofactor, or adjugate of obvious matrices.

2. At n=3, M is independent of r (degree 0 in r), so trivially even.

3. At n=4, M = 2*r^2 + Q_0(s) where Q_0 is degree 2 in s.
   The r^2 coefficient is the universal constant 2 = (n-2)!/2^{n-3}.
   Wait: for n=4, (n-2)!/2^{n-3} = 2!/2 = 1. But we get 2. Hmm.

4. The r^2 coefficient of M[a,b] at n=4 is always 2 regardless of (a,b).
   This means: the leading coefficient of M is (a,b)-independent!
   This is a STRONG constraint that should have a clean proof.

5. M[a,b](-r) = M[b,a](r) is proved and holds for ALL entries.
""")
