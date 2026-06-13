#!/usr/bin/env python3
"""
Does M[a,b] relate to adj(r^2*I - S^2) or similar?

Since S is skew-symmetric, S^2 is symmetric negative-semidefinite.
And r^2*I - S^2 is a symmetric matrix with only r^2 (not r).
Its adjugate/determinant/cofactors are polynomials in r^2 — manifestly even!

If M equals or relates to a cofactor of r^2*I - S^2, that would
explain the even-r-powers property.

opus-2026-03-06-S21
"""

from itertools import permutations
from sympy import (symbols, expand, Matrix, eye, ones, det, simplify,
                   Poly, factor, sqrt)

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

print("=" * 70)
print("S-SQUARED CONNECTION")
print("=" * 70)

# ============================================================
# Part 1: Compute M and adj(r^2*I - S^2) at n=4
# ============================================================
print("\n--- Part 1: n=4 comparison ---")

n = 4
r, sv, s, t = setup(n)
a, b = 0, 1

# Build S matrix
S_mat = Matrix(n, n, lambda i, j: s(i,j) if i != j else 0)
I_n = eye(n)

# S^2
S2 = S_mat * S_mat
S2_simplified = Matrix(n, n, lambda i, j: expand(S2[i,j]))
print(f"  S^2 diagonal [0,0] = {S2_simplified[0,0]}")
print(f"  S^2 off-diag [0,1] = {S2_simplified[0,1]}")
print(f"  S^2 is symmetric? {all(expand(S2[i,j] - S2[j,i]) == 0 for i in range(n) for j in range(n))}")

# r^2*I - S^2
X = r**2 * I_n - S2
adj_X = X.adjugate()

M01 = transfer_M(t, n, a, b)

print(f"\n  M[0,1] = {M01}")
print(f"  adj(r^2*I - S^2)[0,1] = {expand(adj_X[0,1])}")
print(f"  adj(r^2*I - S^2)[1,0] = {expand(adj_X[1,0])}")
# Since r^2*I - S^2 is symmetric, adj should be symmetric
print(f"  adj symmetric? {expand(adj_X[0,1] - adj_X[1,0]) == 0}")

# Check if M relates to adj by a scalar
print(f"\n  Checking if M = c * adj(r^2*I - S^2)[0,1] for some constant c...")
ratio = None
p_M = Poly(M01, r)
p_adj = Poly(expand(adj_X[0,1]), r)
if p_M.degree() > 0 and p_adj.degree() > 0:
    # Compare leading coefficients
    r_M = p_M.nth(p_M.degree())
    r_adj = p_adj.nth(p_adj.degree())
    if r_adj != 0:
        ratio_lead = expand(r_M / r_adj)
        print(f"  Leading coefficient ratio: {ratio_lead}")

        # Check if this ratio works for all coefficients
        check = expand(M01 - ratio_lead * adj_X[0,1])
        print(f"  M - ratio*adj = {check}")

# ============================================================
# Part 2: Try adj(r^2*I + S^2) — note S^2 has NEGATIVE eigenvalues
# ============================================================
print("\n--- Part 2: Try r^2*I + S^2 (adding, not subtracting) ---")

Y = r**2 * I_n + S2
adj_Y = Y.adjugate()
print(f"  adj(r^2*I + S^2)[0,1] = {expand(adj_Y[0,1])}")
check2 = expand(M01 - adj_Y[0,1])
if check2 != 0:
    # Try various scalar multiples
    for c_val in [1, -1, 2, -2, 1/2, -1/2]:
        check = expand(M01 - c_val * adj_Y[0,1])
        if check == 0:
            print(f"  M = {c_val} * adj(r^2*I + S^2)[0,1] ✓")
            break
    else:
        print(f"  M ≠ c * adj(r^2*I + S^2)[0,1] for simple c")

# ============================================================
# Part 3: Try the Cayley-Hamilton type identities
# ============================================================
print("\n--- Part 3: Cayley-Hamilton approach ---")

# For S skew, the Cayley-Hamilton theorem gives:
# S^n + c_{n-2} S^{n-2} + ... = 0
# For n=4: S^4 + c_2 S^2 + c_0 I = 0

char_coeffs = Poly(det(symbols('x') * I_n - S_mat), symbols('x'))
c2 = expand(char_coeffs.nth(2))
c0 = expand(char_coeffs.nth(0))
print(f"  Char poly: x^4 + ({c2})*x^2 + ({c0})")

# So S^4 + c2*S^2 + c0*I = 0, meaning S^4 = -c2*S^2 - c0*I
S4 = S_mat**4
S4_check = expand((S4 + c2 * S2 + c0 * I_n)[0,0])
print(f"  Cayley-Hamilton check [0,0]: {S4_check}")

# ============================================================
# Part 4: Build the "augmented" skew matrix approach
# ============================================================
print("\n--- Part 4: Augmented matrix ---")
print("  Build a (n+2) x (n+2) skew matrix that encodes M...")

# Idea: extend S with two new rows/columns for endpoints a,b
# and set up r-dependent entries that give M as a Pfaffian.

# For a (2k)x(2k) skew matrix, Pf^2 = det.
# If we extend our 4x4 S to a 6x6 skew matrix with r in two rows...

# At n=4, M = 2r^2 + Q. Can this be a Pfaffian?
# Pf of 6x6 skew has 15 terms. M has 7 terms. So possibly as a
# Pfaffian with structured zeros.

# Let me try a specific construction. Add two "auxiliary" vertices
# x, y (representing endpoints a, b) to the skew matrix:

# X = [[0, s01, s02, s03, ?, ?],
#      [-s01, 0, s12, s13, ?, ?],
#      [-s02, -s12, 0, s23, ?, ?],
#      [-s03, -s13, -s23, 0, ?, ?],
#      [?, ?, ?, ?, 0, ?],
#      [?, ?, ?, ?, ?, 0]]

# For M[0,1]: a=0, b=1. Let x=4, y=5 be auxiliary.
# We want Pf(X) or some function of X to give M[0,1].

# Since M is manifestly even in r, and Pfaffians of skew matrices
# with r-dependent entries CAN be even, this might work.

# ============================================================
# Part 5: Direct numerical test — is M a cofactor of (rI-S)(rI+S)?
# ============================================================
print("\n--- Part 5: Products (rI-S)(rI+S) = r^2*I - S^2 ---")

# (rI-S)(rI+S) = r^2 I - S^2 (since SI = IS)
# So adj((rI-S)(rI+S)) = adj(r^2 I - S^2).
# But det((rI-S)(rI+S)) = det(rI-S)*det(rI+S) = Pf(S)^2 * ... hmm.
# For skew S: det(rI-S) = r^4 + c2*r^2 + c0 for n=4.
# det(rI+S) = (-1)^4 det(-rI-S) = det(rI+S).
# Actually det(rI+S) = det(rI - (-S)) = det(rI-S^T) = det((rI-S)^T) = det(rI-S).
# So det((rI-S)(rI+S)) = det(rI-S)^2.

# And adj((rI-S)(rI+S)) = adj(rI+S) * adj(rI-S).
# = adj((rI-S)^T) * adj(rI-S) = adj(rI-S)^T * adj(rI-S).

# The (a,b) entry of adj(rI-S)^T * adj(rI-S) is:
# sum_k adj(rI-S)[k,a] * adj(rI-S)[k,b]
# = sum_k adj[a,k]^T * adj[k,b] = (adj^T * adj)[a,b]

# This is an inner product of the a-th and b-th columns of adj(rI-S).
# For a symmetric positive definite matrix, this would be related to
# the Gram matrix. Since adj(rI-S) involves r, this is a polynomial.

# Let's compute adj(rI-S) and check
rI_minus_S = r * I_n - S_mat
adj_rS = rI_minus_S.adjugate()

inner_01 = 0
for k in range(n):
    inner_01 += adj_rS[k, a] * adj_rS[k, b]
inner_01 = expand(inner_01)

print(f"  sum_k adj(rI-S)[k,0]*adj(rI-S)[k,1] = {inner_01}")
print(f"  M[0,1] = {M01}")
print(f"  Equal? {expand(inner_01 - M01) == 0}")

# Try the cross-column inner product with just rows not a,b
inner_01_restricted = 0
for k in range(n):
    if k == a or k == b: continue
    inner_01_restricted += adj_rS[k, a] * adj_rS[k, b]
inner_01_restricted = expand(inner_01_restricted)
print(f"\n  sum_{k not in {a,b}} adj[k,0]*adj[k,1] = {inner_01_restricted}")

# What about just the (a,b) cofactor of (rI-S)(rI+S) = r^2*I - S^2?
# The cofactor is (-1)^{a+b} * det(minor)
minor_ab = (r**2 * I_n - S2).minor_submatrix(a, b)
cof_ab = (-1)**(a+b) * det(minor_ab)
cof_ab = expand(cof_ab)
print(f"\n  cofactor(r^2*I - S^2, {a}, {b}) = {cof_ab}")
print(f"  M[0,1] = {M01}")
print(f"  Equal? {expand(cof_ab - M01) == 0}")

# Since r^2*I - S^2 is symmetric, cof(a,b) = cof(b,a)
cof_ba = (-1)**(b+a) * det((r**2 * I_n - S2).minor_submatrix(b, a))
cof_ba = expand(cof_ba)
print(f"  cofactor(r^2*I - S^2, {b}, {a}) = {cof_ba}")
print(f"  cof(a,b) = cof(b,a)? {expand(cof_ab - cof_ba) == 0}")

# ============================================================
# Part 6: Try restricted matrix — delete rows/cols a,b
# ============================================================
print("\n--- Part 6: Restricted to U = {2,...,n-1} ---")

# The "internal" skew matrix S_U = S restricted to U
S_U = Matrix(len(U_l := list(range(2,n))), len(U_l),
             lambda i, j: s(U_l[i], U_l[j]) if U_l[i] != U_l[j] else 0)
print(f"  S_U = {S_U}")

# det(rI - S_U) for the 2x2 case
det_rI_SU = expand(det(r * eye(len(U_l)) - S_U))
print(f"  det(rI - S_U) = {det_rI_SU}")

# Does M relate to det(rI - S_U) somehow?
print(f"  M[0,1] = {M01}")
print(f"  det(rI-S_U) = {det_rI_SU}")

# ============================================================
# Part 7: Try M = f(adj(rI-S)) for some function f
# ============================================================
print("\n--- Part 7: Trying various combinations of adj(rI-S) ---")

# adj(rI-S)[i,j] is a degree-(n-1) polynomial in r.
# M is degree n-2. So adj/r might work?

# adj(rI-S)[0,1] / r = ?
adj_01 = expand(adj_rS[0,1])
p_adj01 = Poly(adj_01, r)
print(f"  adj(rI-S)[0,1] = {adj_01}")
print(f"  Degree: {p_adj01.degree()}")
print(f"  Coefficients:")
for k in range(p_adj01.degree() + 1):
    print(f"    r^{k}: {expand(p_adj01.nth(k))}")

# What about adj(rI-S)[0,1] + adj(rI-S)[1,0]?
adj_10 = expand(adj_rS[1,0])
sym_adj = expand(adj_01 + adj_10)
print(f"\n  adj[0,1] + adj[1,0] = {sym_adj}")
p_sym = Poly(sym_adj, r)
print(f"  Degree: {p_sym.degree()}")
for k in range(p_sym.degree() + 1):
    ck = expand(p_sym.nth(k))
    if ck != 0:
        print(f"    r^{k}: {ck}")

# Since S is skew, (rI-S)^T = rI+S, so adj(rI-S)^T = adj(rI+S).
# adj(rI+S)[j,i] = adj(rI-S)[i,j] ... is this right?
# adj(M^T) = adj(M)^T. So adj(rI+S) = adj((rI-S)^T) = adj(rI-S)^T.
# So adj(rI+S)[i,j] = adj(rI-S)[j,i].

# adj(rI-S)[0,1] + adj(rI+S)[0,1] = adj(rI-S)[0,1] + adj(rI-S)[1,0]
# = symmetric part of adj.
print(f"\n  This is the symmetric part of adj(rI-S).")
print(f"  For S skew: adj(rI-S)[0,1] + adj(rI-S)[1,0] has only even r-powers")
print(f"    if adj(rI-S) satisfies adj(-r-S) = adj(rI-S)^T = adj(rI+S)")

# Check: adj(rI-S)[0,1] at -r vs adj(rI-S)[1,0] at r
adj_01_neg = expand(adj_01.subs(r, -r))
print(f"\n  adj(rI-S)[0,1](-r) = {adj_01_neg}")
print(f"  adj(rI-S)[1,0](r) = {adj_10}")
print(f"  Equal? {expand(adj_01_neg - adj_10) == 0}")
# adj(-rI-S) = adj(-(rI+S)) = (-1)^{n-1} adj(rI+S) = (-1)^{n-1} adj(rI-S)^T
# So adj(-rI-S)[0,1] = (-1)^{n-1} adj(rI-S)[1,0]
# = (-1)^3 adj[1,0] = -adj[1,0] for n=4
print(f"  Expected: adj(-rI-S)[0,1] = (-1)^(n-1) adj(rI-S)[1,0] = {(-1)**(n-1)} * adj[1,0]")
check = expand(adj_01_neg - (-1)**(n-1) * adj_10)
print(f"  Check: {check == 0}")

# So: adj(rI-S)[0,1] + (-1)^{n-1} * adj(rI-S)[0,1](-r) ... hmm.
# This means the SYMMETRIC part adj[0,1] + adj[1,0] has:
# = adj[0,1](r) + adj[0,1](-r) * (-1)^{n-1} / (-1)^{n-1}
# Actually: adj[1,0] = (-1)^{n-1} adj[0,1](-r)... not quite.

# Let me just check the parity of the symmetric part
p_sym_check = Poly(sym_adj, r)
print(f"\n  Parity check of adj[0,1] + adj[1,0]:")
for k in range(p_sym_check.degree() + 1):
    ck = expand(p_sym_check.nth(k))
    if ck != 0:
        parity = "even" if k % 2 == 0 else "ODD"
        print(f"    r^{k} [{parity}]: {ck}")

# ============================================================
# Part 8: n=5 check
# ============================================================
print("\n--- Part 8: n=5 --- ")
n = 5
r, sv, s, t = setup(n)
a, b = 0, 1

S_mat5 = Matrix(n, n, lambda i, j: s(i,j) if i != j else 0)
I_5 = eye(n)
S2_5 = S_mat5 * S_mat5

M01_5 = transfer_M(t, n, a, b)

# Cofactor of r^2*I - S^2
X5 = r**2 * I_5 - S2_5
minor_ab5 = X5.minor_submatrix(a, b)
cof_ab5 = (-1)**(a+b) * det(minor_ab5)
cof_ab5 = expand(cof_ab5)

print(f"  M[0,1] = {M01_5}")
print(f"  cofactor(r^2*I - S^2, 0, 1) = {cof_ab5}")
print(f"  Equal? {expand(cof_ab5 - M01_5) == 0}")

# Check ratio
if M01_5 != 0 and cof_ab5 != 0:
    p_M5 = Poly(M01_5, r)
    p_cof5 = Poly(cof_ab5, r)
    for k in range(max(p_M5.degree(), p_cof5.degree()) + 1):
        mval = expand(p_M5.nth(k))
        cval = expand(p_cof5.nth(k))
        if mval != 0 and cval != 0:
            print(f"    r^{k}: M={mval}, cof={cval}")

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
