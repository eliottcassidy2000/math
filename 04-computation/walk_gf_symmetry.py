"""
Irving-Omar walk generating function and transfer matrix symmetry.

W(z) = det(I + z*A^T) / det(I - z*A)

For a c-tournament, A + A^T = c*(J-I), so A^T = c*(J-I) - A.
Thus I + z*A^T = I + z*(c*(J-I) - A) = (1-cz)*I + cz*J - z*A

Question: does W(z) encode M[a,b] in a way that makes symmetry manifest?

Also: the "pointed" version of W — extracting the (a,b) coefficient —
might directly relate to M[a,b].
"""

from sympy import symbols, expand, Matrix, eye, ones, det, Rational, factor, collect
import random
from itertools import permutations

print("=" * 70)
print("WALK GENERATING FUNCTION AND SYMMETRY")
print("=" * 70)

# ==========================================================
# Part 1: Symbolic W(z) at n=3
# ==========================================================
print("\n--- Part 1: n=3 symbolic ---")
n = 3
z = symbols('z')
c = symbols('c')

# Skew variables
svars = {}
for i in range(n):
    for j in range(i+1, n):
        svars[(i,j)] = symbols(f's{i}{j}')

# Arc matrix: A_ij = c/2 + s_ij
A = Matrix(n, n, lambda i, j: 0 if i == j else c/2 + (svars[(i,j)] if i < j else -svars[(j,i)]))
print(f"  A = {A}")
print(f"  A + A^T = {expand(A + A.T)}")

I_n = eye(n)
J = ones(n, n)

# W(z) = det(I + z*A^T) / det(I - z*A)
numer = expand(det(I_n + z * A.T))
denom = expand(det(I_n - z * A))
print(f"\n  det(I + z*A^T) = {numer}")
print(f"  det(I - z*A) = {denom}")

# Factor in z
from sympy import Poly
numer_poly = Poly(numer, z)
denom_poly = Poly(denom, z)
print(f"\n  Numerator coefficients by z-power:")
for i in range(numer_poly.degree() + 1):
    print(f"    z^{i}: {expand(numer_poly.nth(i))}")
print(f"  Denominator coefficients by z-power:")
for i in range(denom_poly.degree() + 1):
    print(f"    z^{i}: {expand(denom_poly.nth(i))}")

# ==========================================================
# Part 2: The adjugate / inverse of (I - zA) at n=3
# ==========================================================
print(f"\n--- Part 2: (I - zA)^(-1) structure ---")

# adj(I - zA) / det(I - zA) = (I - zA)^(-1)
# The (a,b) entry of adj(I-zA) is the (b,a)-cofactor of (I-zA).
# This is a polynomial in z.

IzA = I_n - z * A
adj_IzA = IzA.adjugate()

print(f"  adj(I - zA):")
for i in range(n):
    for j in range(n):
        entry = expand(adj_IzA[i,j])
        print(f"    [{i},{j}]: {entry}")

# The coefficient of z^{n-1} in adj(I-zA)[a,b] counts...
# what exactly?

# For the INVERSE: (I-zA)^{-1} = sum_{k>=0} z^k A^k
# So (I-zA)^{-1}[a,b] = sum_k z^k A^k[a,b]
# = sum_k z^k (number of walks of length k from a to b)
# But A^k[a,b] counts weighted walks (with arc weights c/2 + s_ij).

# adj(I-zA)[a,b] = det(I-zA) * (I-zA)^{-1}[a,b]
# So adj(I-zA)[a,b] / det(I-zA) = (I-zA)^{-1}[a,b]

# The HAMILTONIAN path count from a to b is related to the
# coefficient of z^{n-1} in det(I-zA)^{-1} * ... hmm, this is the
# permanent formulation.

# ==========================================================
# Part 3: Connection to transfer matrix M
# ==========================================================
print(f"\n--- Part 3: M[a,b] vs matrix entries ---")

# Compute M[a,b] directly
def hp(T_func, vset, start=None, end=None):
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
            prod *= T_func(perm[ii], perm[ii+1])
        total += prod
    return expand(total)

def arc(i, j):
    if i == j: return 0
    if i < j: return c/2 + svars[(i,j)]
    return c/2 - svars[(j,i)]

for a in range(n):
    for b in range(n):
        if a == b: continue
        U = [v for v in range(n) if v != a and v != b]
        M_val = 0
        for mask in range(1 << len(U)):
            S = [U[ii] for ii in range(len(U)) if mask & (1 << ii)]
            R = [U[ii] for ii in range(len(U)) if not (mask & (1 << ii))]
            sign = (-1)**len(S)
            ea = hp(arc, set(S)|{a}, end=a)
            bb = hp(arc, set(R)|{b}, start=b)
            M_val += sign * ea * bb
        M_val = expand(M_val)

        # Compare with adj(I-zA) entries
        adj_entry = expand(adj_IzA[a,b])
        # Get coefficient of z^{n-1}
        adj_coeff = expand(Poly(adj_entry, z).nth(n-1))

        print(f"  M[{a},{b}] = {M_val}")
        print(f"  adj(I-zA)[{a},{b}] z^{n-1} coeff = {adj_coeff}")
        print(f"  Match? {expand(M_val - adj_coeff) == 0}")
        print()

# ==========================================================
# Part 4: n=4 numerical test of adj connection
# ==========================================================
print(f"\n--- Part 4: n=4 numerical test ---")
random.seed(42)

n4 = 4
for trial in range(5):
    c_val = random.uniform(0.5, 2.0)
    s_vals = {}
    for i in range(n4):
        for j in range(i+1, n4):
            s_vals[(i,j)] = random.uniform(-1, 1)

    # Build A
    A_num = [[0.0]*n4 for _ in range(n4)]
    for i in range(n4):
        for j in range(n4):
            if i == j: continue
            if i < j:
                A_num[i][j] = c_val/2 + s_vals[(i,j)]
            else:
                A_num[i][j] = c_val/2 - s_vals[(j,i)]

    # Build I - zA as matrix at z=1... no, we need symbolic z
    # Let's compute adj(I-A) = adj at z=1
    import numpy as np
    A_np = np.array(A_num)
    IminA = np.eye(n4) - A_np

    # adj(I-A) = det(I-A) * inv(I-A)
    det_IminA = np.linalg.det(IminA)
    if abs(det_IminA) > 1e-10:
        inv_IminA = np.linalg.inv(IminA)
        adj_IminA = det_IminA * inv_IminA
    else:
        adj_IminA = None

    # Compute M[a,b] directly
    def hp_num(vset, start=None, end=None):
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
                prod *= A_num[perm[ii]][perm[ii+1]]
            total += prod
        return total

    a_t, b_t = 0, 1
    U_t = [v for v in range(n4) if v != a_t and v != b_t]
    m_val = 0
    for mask in range(1 << len(U_t)):
        S = [U_t[ii] for ii in range(len(U_t)) if mask & (1 << ii)]
        R = [U_t[ii] for ii in range(len(U_t)) if not (mask & (1 << ii))]
        sign = (-1)**len(S)
        ea = hp_num(set(S)|{a_t}, end=a_t)
        bb = hp_num(set(R)|{b_t}, start=b_t)
        m_val += sign * ea * bb

    # Compare
    if adj_IminA is not None:
        adj_val = adj_IminA[a_t, b_t]
        print(f"  Trial {trial+1}: M[0,1]={m_val:.6f}, adj(I-A)[0,1]={adj_val:.6f}, "
              f"ratio={m_val/adj_val:.6f}" if abs(adj_val) > 1e-10 else f"adj=~0")
    else:
        print(f"  Trial {trial+1}: det(I-A)~0, skipping")

# ==========================================================
# Part 5: The (I - zA)^{-1} power series at n=4
# ==========================================================
print(f"\n--- Part 5: Power series analysis ---")
print("""
(I - zA)^{-1} = I + zA + z²A² + z³A³ + ...

The (a,b) entry at z^k gives A^k[a,b] = sum of weighted walks of length k from a to b.

For Hamiltonian paths, we want walks that visit each vertex exactly once.
The permanent of A (restricted to avoid revisits) gives the count.

The transfer matrix M[a,b] = sum_S (-1)^|S| E_a(S+a) B_b(R+b) is an
INCLUSION-EXCLUSION formula that extracts the Hamiltonian contribution.

Key identity (Irving-Omar): ham(D) = sum_S det(Abar[S]) * per(A[S^c])
where Abar = complement adjacency.

For c-tournaments: Abar[i,j] = 1 - A[i,j] = 1 - c/2 - s_ij = (2-c)/2 - s_ij = (2-c)/2 + s_ji.
So Abar is ALSO a c'-tournament with c' = 2-c!

This is beautiful: the complement of a c-tournament is a (2-c)-tournament.
At c=1: complement is also c=1 (tournaments are self-complementary in this sense).
At c=0: complement is c=2.
""")

# Verify: Abar + Abar^T = ?
# Abar_ij = 1 - A_ij = 1 - c/2 - s_ij (for i≠j), 0 for i=j
# Abar_ji = 1 - A_ji = 1 - c/2 + s_ij
# Abar_ij + Abar_ji = 2 - c = (2-c) ✓

print(f"  Verification: A + A^T = c*(J-I), (J-I-A) + (J-I-A)^T = (2-c)*(J-I) ✓")

# ==========================================================
# Part 6: The key algebraic insight
# ==========================================================
print(f"\n--- Part 6: Algebraic structure ---")
print("""
For a c-tournament A with A + A^T = c*(J-I):

  A = (c/2)*(J-I) + S    where S is skew-symmetric (S = -S^T)

  det(I - zA) = det((1-cz/2+cz/(2n))I - z*S - (cz/2)*(J-I))

Actually, since J has eigenvalue n on the all-ones vector and 0 on
its orthogonal complement:

  A = (c/2)*(J-I) + S

The eigenvalues of (c/2)*(J-I) are:
  - (c/2)*(n-1) on the all-ones vector
  - -c/2 on the (n-1)-dimensional complement

So A has the spectral structure of S perturbed by a rank-1 shift.

For S skew-symmetric: eigenvalues are purely imaginary ±iλ_k (plus 0 if n odd).
For A = (c/2)(J-I) + S: eigenvalues are -c/2 ± iλ_k (on the complement)
plus c(n-1)/2 + 0 = c(n-1)/2 on the all-ones vector.

The SYMMETRY of M might follow from the spectral structure:
  - S is skew-symmetric => its "transfer matrix" has definite parity
  - The rank-1 perturbation (c/2)(J-I) is symmetric, adding only to the even part

This is consistent with: M = c^{n-2} * f(n) + polynomial in S (even/odd).
""")

# Verify spectral structure numerically
print(f"\n  Numerical eigenvalue check (n=5, c=1):")
n5 = 5
random.seed(42)
T5 = [[0]*n5 for _ in range(n5)]
for i in range(n5):
    for j in range(i+1, n5):
        if random.random() < 0.5:
            T5[i][j] = 1
        else:
            T5[j][i] = 1

A5 = np.array(T5, dtype=float)
S5 = A5 - A5.T  # skew part (times 2)
eig_A = np.linalg.eigvals(A5)
eig_S = np.linalg.eigvals(S5 / 2)

print(f"  Eigenvalues of A: {sorted(eig_A, key=lambda x: -x.real)}")
print(f"  Eigenvalues of S/2: {sorted(eig_S, key=lambda x: -abs(x.imag))}")
print(f"  A = (1/2)(J-I) + S/2")
print(f"  Expected: eig(A) ≈ -1/2 + eig(S/2), plus (n-1)/2 = {(n5-1)/2}")

# Check
J5 = np.ones((n5,n5)) - np.eye(n5)
A5_check = 0.5 * J5 + S5/2
print(f"  A = (1/2)(J-I) + S/2? {np.allclose(A5, A5_check)}")

print(f"\n" + "=" * 70)
print("KEY INSIGHT")
print("=" * 70)
print("""
The c-tournament structure A = (c/2)(J-I) + S decomposes A into:
  1. A SYMMETRIC rank-1 perturbation: (c/2)(J-I)
  2. A SKEW-SYMMETRIC part: S

The transfer matrix M[a,b] = sum_S (-1)^|S| E_a B_b is a polynomial
in the entries of A. Under the decomposition:
  M = (polynomial in c and S-entries)

The symmetry M[a,b] = M[b,a] is equivalent to M being EVEN in S (n even)
or ODD in S (n odd), which follows from:
  M(c, -S) = (-1)^{n-2} M(c, S)   [T^op equivalence]

The CHALLENGE: prove this identity for the specific polynomial M.

The MOST PROMISING ROUTE: use the spectral structure.
Since S is skew and J-I is symmetric with known spectrum,
the Hamiltonian path generating function has specific parity properties.
""")
