#!/usr/bin/env python3
"""
fibonacci_determinant.py — opus-2026-03-12-S67d

CAN H BE EXPRESSED AS A FIBONACCI-TYPE DETERMINANT?

Key idea: The Cassini identity B_{m-1}·B_{m+1} - B_m² = x
is equivalent to det([[B_m, B_{m-1}], [B_{m+1}, B_m]]) = -x.
This 2×2 determinant pattern is the MV transfer matrix.

QUESTION: Can we build an m×m matrix whose DETERMINANT or PERMANENT
gives H(Int, p) where p = 2m+1?

APPROACH 1: Tridiagonal determinant
  The eigenvalues Q_k of the Interval tournament come from a
  tridiagonal (Jacobi) matrix J with entries related to sin/cos.
  det(xI - J) = B_m(x) (characteristic polynomial)
  Could det(f(J)) = H for some function f?

APPROACH 2: Permanent of circulant
  H counts Hamiltonian paths = restricted permanent.
  For circulant matrices, permanent formulas exist (van der Waerden-like).

APPROACH 3: Fibonacci matrix construction
  Build an m×m matrix M(p) such that det(M) = H/p or perm(M) = H/p.
  Use Morgan-Voyce structure as guide.
"""

import numpy as np
from math import comb, factorial, gcd
from collections import Counter

def morgan_voyce_B(m, x):
    return sum(comb(m+j, 2*j) * x**j for j in range(m+1))

def fib(n):
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a+b
    return a

def interval_Q(p):
    m = (p-1)//2
    return np.array([
        (np.sin(m*np.pi*k/p)/np.sin(np.pi*k/p))**2
        for k in range(1, m+1)
    ])

def legendre(a, p):
    a = a % p
    if a == 0: return 0
    ls = pow(a, (p-1)//2, p)
    return -1 if ls == p-1 else ls

print("=" * 70)
print("FIBONACCI DETERMINANT — opus-2026-03-12-S67d")
print("=" * 70)
print()

H_int = {3: 3, 5: 25, 7: 175, 11: 93027, 13: 3711175}

# ============================================================
# PART 1: JACOBI MATRIX AND CHARACTERISTIC POLYNOMIAL
# ============================================================

print("PART 1: JACOBI MATRIX FOR Q-SPECTRUM")
print("=" * 70)
print()

# The Q_k are eigenvalues of an m×m matrix.
# For the Interval tournament, Q_k = (sin(mπk/p)/sin(πk/p))²
# These are eigenvalues of D_m^T D_m where D_m is the Dirichlet kernel matrix.
# But they are also eigenvalues of a TRIDIAGONAL matrix!

# The tridiagonal matrix with eigenvalues Q_k:
# Since Q_k = |D_m(2πk/p)|², and D_m(θ) = sin(mθ/2)/sin(θ/2) (Dirichlet kernel),
# the matrix that has these as eigenvalues relates to a Chebyshev problem.

# Actually, the matrix is: M = 2I + N + N^T where N is a nilpotent shift
# with specific entries. For a circulant S = {1,...,m}, the autocorrelation
# matrix r_S gives rise to a Toeplitz structure.

for p in [5, 7, 11, 13]:
    m = (p-1)//2
    Q = interval_Q(p)
    H = H_int[p]
    h = H // p

    # Construct the m×m matrix whose eigenvalues are Q_k
    # Using the DFT relationship: Q_k = |Σ_{s∈S} ω^{sk}|²
    # For S = {1,...,m}, this is |Σ_{j=1}^m ω^{jk}|² = |D_m(2πk/p)|²

    # The matrix whose eigenvalues are Q_k can be constructed as
    # A A^* where A is the m×p DFT restricted to rows 1..m

    # But let's try a different approach: find a NICE m×m matrix
    # whose determinant relates to H.

    print(f"p={p}, m={m}: H/p = {h}, F_p = {fib(p)}")

    # Tridiagonal Toeplitz with constant diagonal a and off-diagonal b:
    # eigenvalues = a + 2b·cos(kπ/(m+1)) for k=1,...,m
    # det = prod eigenvalues

    # Try: what values of (a, b) make det = h?
    # For m=1: trivially a = h
    # For m=2: a² - b² = h
    # For m=3: a(a² - 2b²) = h

    # More interesting: can we build a matrix from the Q_k?
    # det(Q_diag) = prod(Q_k) = 1 for p≡3 mod 4 (when all Q_k are real)

    prod_Q = np.prod(Q)
    print(f"  prod(Q_k) = {prod_Q:.6f}")
    print(f"  det(diag(Q)) = prod(Q_k) = {prod_Q:.6f}")

    # det(I + Q_diag) = prod(1 + Q_k) = F_p
    det_IQ = np.prod(1 + Q)
    print(f"  det(I + diag(Q)) = {det_IQ:.1f} = F_p")

    # det(2I + Q_diag) = prod(2 + Q_k) = 2^m B_m(1/2)
    det_2Q = np.prod(2 + Q)
    Bm_half = 2**m * morgan_voyce_B(m, 0.5)
    print(f"  det(2I + diag(Q)) = {det_2Q:.1f}, 2^m B_m(1/2) = {Bm_half:.1f}")

    # What value of 'a' gives det(aI + diag(Q)) = h = H/p?
    # This is a^m B_m(1/a) = h
    # We need to solve for a.
    # At a = h^{1/m} approximately:
    a_approx = h ** (1/m)
    det_a = np.prod(a_approx + Q)
    print(f"  a s.t. det(aI+Q)=H/p: a ≈ {a_approx:.4f}, det = {det_a:.1f}")

    # Actually, is there a MATRIX whose eigenvalues give H?
    # Try: matrix M with eigenvalues Q_k, then some function of M...
    # det(2M + I) = prod(2Q_k + 1) = ?

    det_2M1 = np.prod(2*Q + 1)
    print(f"  det(2·diag(Q)+I) = prod(2Q+1) = {det_2M1:.1f}")
    print(f"  H/p / prod(2Q+1) = {h / det_2M1:.6f}")
    print()

# ============================================================
# PART 2: TOEPLITZ DETERMINANTS
# ============================================================

print("=" * 70)
print("PART 2: TOEPLITZ DETERMINANT APPROACH")
print("=" * 70)
print()

# For a circulant tournament, the autocorrelation r_S(t) = |{(a,b)∈S²: a+b≡t}|
# forms a Toeplitz structure.
# The Toeplitz determinant D_m(f) = det of m×m Toeplitz matrix with
# entries c_{i-j} = Fourier coefficients of f(θ).
# By Szegő's theorem: D_m(f) → exp(m·∫log f dθ + ...)

# For the Interval tournament, the "symbol" is
# f(θ) = |D_m(θ)|² = Q(θ)
# The Toeplitz determinant is det T_m where T_m[i,j] = ∫f(θ)e^{i(i-j)θ}dθ

# Let's compute: the Toeplitz matrix T with symbol |D_m(θ)|²

for p in [5, 7, 11]:
    m = (p-1)//2
    H = H_int[p]
    h = H // p

    # Compute Fourier coefficients of |D_m(θ)|²
    # |D_m(θ)|² = Σ r_S(n) e^{inθ} where r_S is the autocorrelation
    S = set(range(1, m+1))

    r_S = [0] * p
    for a in S:
        for b in S:
            r_S[(a+b) % p] += 1

    print(f"p={p}, m={m}: r_S = {r_S}")

    # Build Toeplitz matrix T using r_S as the symbol
    # T[i,j] = r_S[(j-i) mod p]
    T = np.zeros((m, m))
    for i in range(m):
        for j in range(m):
            T[i,j] = r_S[(j-i) % p]

    det_T = np.linalg.det(T)
    print(f"  det(T_{m}) = {det_T:.4f}")
    print(f"  H/p = {h}")
    print(f"  H/p / det = {h / det_T:.6f}" if abs(det_T) > 0.01 else "  det ≈ 0")

    # Try p×p Toeplitz (circulant)
    C = np.zeros((p, p))
    for i in range(p):
        for j in range(p):
            C[i,j] = r_S[(j-i) % p]

    det_C = np.linalg.det(C)
    print(f"  det(C_{p}) = {det_C:.4f}")

    # The circulant C has eigenvalues = p · Q_k (plus Q_0 = m²)
    # So det(C) = p^p · prod(Q_k) · Q_0 ... hmm
    eigs_C = np.sort(np.linalg.eigvalsh(C))[::-1]
    print(f"  eigenvalues of C: {np.round(eigs_C, 4)}")
    print(f"  Q_k from formula: {np.round(np.sort(interval_Q(p))[::-1], 4)}")
    print()

# ============================================================
# PART 3: FIBONACCI MATRIX AND TOURNAMENT PERMANENT
# ============================================================

print("=" * 70)
print("PART 3: FIBONACCI MATRICES AND PERMANENTS")
print("=" * 70)
print()

# Fibonacci permanent theorem:
# For the n×n matrix F_n = [[1,1,0,...], [1,1,1,0,...], [0,1,1,1,...], ...]
# (tridiagonal with all 1s on diagonal and super/sub-diagonal),
# perm(F_n) = F_{n+1} (Fibonacci number!)

# Does a related permanent give H?

for p in [5, 7]:
    m = (p-1)//2
    H = H_int[p]
    h = H // p

    # Build the adjacency matrix of Interval tournament
    A = np.zeros((p, p), dtype=int)
    S = set(range(1, m+1))
    for i in range(p):
        for j in range(p):
            if i != j and (j-i) % p in S:
                A[i][j] = 1

    # Compute permanent of A (expensive!)
    # For small p, we can use the formula
    # perm(A) = sum over all permutations sigma of prod A[i,sigma(i)]
    from itertools import permutations

    perm_A = 0
    for sigma in permutations(range(p)):
        prod = 1
        for i in range(p):
            prod *= A[i][sigma[i]]
        perm_A += prod

    print(f"p={p}: perm(A) = {perm_A}")
    print(f"  H = {H}")
    print(f"  perm(A) / H = {perm_A / H:.6f}")
    print(f"  perm(A) / p! = {perm_A / factorial(p):.6f}")
    print()

    # The permanent counts perfect matchings in the bipartite graph
    # defined by A. Each matching assigns each row i to a column sigma(i)
    # such that A[i,sigma(i)] = 1.

    # For a tournament, perm(A) counts the number of permutations sigma
    # such that i beats sigma(i) for all i. This is different from
    # Hamiltonian paths!

    # Hamilton paths: perm of a DIFFERENT matrix
    # H = sum over permutations sigma s.t. A[sigma(0),sigma(1)] * A[sigma(1),sigma(2)] * ...
    # This is NOT perm(A) but the permanent of a path-adjacency matrix.

    # Let me compute: path matrix P where P[i,j] = A[i,j] for "step i can go to vertex j"
    # Actually H = trace of A^(p-1) in the "path" sense...
    # No, H is the number of Hamiltonian paths = sum of all directed HP

# ============================================================
# PART 4: THE HANKEL DETERMINANT
# ============================================================

print("=" * 70)
print("PART 4: HANKEL DETERMINANTS AND MOMENTS")
print("=" * 70)
print()

# A Hankel matrix has entries H[i,j] = a_{i+j}
# The Hankel determinant of the Fibonacci sequence is:
# det([[F_i+j]]) = ?

# For our Q_k, the "moments" are μ_n = Σ Q_k^n
# The Hankel matrix H_m with H[i,j] = μ_{i+j} has
# det(H_m) = prod_{i<j} (Q_i - Q_j)² (Vandermonde squared!)

for p in [5, 7, 11, 13]:
    m = (p-1)//2
    Q = interval_Q(p)
    H_val = H_int[p]
    h = H_val // p

    # Moments μ_n = Σ Q_k^n
    moments = [np.sum(Q**n) for n in range(2*m)]

    # Hankel determinant
    Hank = np.zeros((m, m))
    for i in range(m):
        for j in range(m):
            Hank[i,j] = moments[i+j]

    det_Hank = np.linalg.det(Hank)

    # Vandermonde product
    vdm = 1.0
    for i in range(m):
        for j in range(i+1, m):
            vdm *= (Q[i] - Q[j])**2

    print(f"p={p}, m={m}:")
    print(f"  Moments: {np.round(moments[:6], 4)}")
    print(f"  det(Hankel) = {det_Hank:.6f}")
    print(f"  Vandermonde² = {vdm:.6f}")
    print(f"  Ratio det/Vdm = {det_Hank/vdm:.6f}" if abs(vdm) > 1e-10 else "  Vdm ≈ 0")
    print(f"  H/p = {h}")
    print(f"  H/p / det(Hankel) = {h / det_Hank:.6f}" if abs(det_Hank) > 1e-10 else "  det ≈ 0")

    # Try Hankel of B_m coefficients
    bcoeffs = [comb(m+j, 2*j) for j in range(m+1)]
    print(f"  B_m coefficients: {bcoeffs}")

    # Hankel of B_m coefficients
    if m >= 2:
        Hb = np.zeros((m, m))
        for i in range(m):
            for j in range(m):
                idx = i + j
                Hb[i,j] = bcoeffs[idx] if idx < len(bcoeffs) else 0
        det_Hb = np.linalg.det(Hb)
        print(f"  det(Hankel(B_m coeffs)) = {det_Hb:.4f}")
    print()

# ============================================================
# PART 5: FIBONACCI SEQUENCE OPERATIONS ON H/p
# ============================================================

print("=" * 70)
print("PART 5: FIBONACCI OPERATIONS ON H/p")
print("=" * 70)
print()

# H/p values: 1, 5, 25, 8457, 285475, 805251147
# Check various Fibonacci-related operations

h_vals = [(3, 1), (5, 5), (7, 25), (11, 8457), (13, 285475), (17, 805251147)]

print("H/p mod F_p:")
for p, h in h_vals:
    Fp = fib(p)
    print(f"  p={p}: H/p = {h}, F_p = {Fp}, (H/p) mod F_p = {h % Fp}")

print()
print("H/p in terms of F_p:")
for p, h in h_vals:
    Fp = fib(p)
    q, r = divmod(h, Fp)
    print(f"  p={p}: H/p = {q}·F_p + {r}")

print()

# Check: does H/p = sum of products of Fibonacci numbers?
print("Factorizations of H/p:")
for p, h in h_vals:
    # Factor h
    n = h
    factors = []
    for d in range(2, min(n+1, 100000)):
        while n % d == 0:
            factors.append(d)
            n //= d
    if n > 1:
        factors.append(n)
    print(f"  p={p}: H/p = {h} = {' × '.join(str(f) for f in factors)}")

print()

# ============================================================
# PART 6: THE VANDERMONDE APPROACH
# ============================================================

print("=" * 70)
print("PART 6: VANDERMONDE AND Q-SPECTRUM")
print("=" * 70)
print()

# The Vandermonde determinant V(Q) = prod_{i<j}(Q_i - Q_j) measures
# how SPREAD OUT the Q_k are. For Interval, this is large (concentrated Q
# means large differences).

for p in [7, 11, 13]:
    m = (p-1)//2
    Q = np.sort(interval_Q(p))[::-1]
    H_val = H_int[p]
    h = H_val // p

    # Vandermonde
    V = 1.0
    for i in range(m):
        for j in range(i+1, m):
            V *= (Q[i] - Q[j])

    print(f"p={p}: V(Q) = {V:.4f}")
    print(f"  |V(Q)| = {abs(V):.4f}")
    print(f"  H/p = {h}")
    print(f"  H/p / |V| = {h / abs(V):.4f}")

    # Also: discriminant = V²
    disc = V**2
    print(f"  disc(Q) = V² = {disc:.4f}")
    print(f"  H/p / disc = {h / disc:.6f}" if abs(disc) > 1e-10 else "")
    print()

# ============================================================
# PART 7: MATRIX WHOSE PERMANENT IS H/p — CONSTRUCTION ATTEMPT
# ============================================================

print("=" * 70)
print("PART 7: SEEKING A MATRIX WITH perm(M) = H/p")
print("=" * 70)
print()

# For p=5, m=2: H/p = 5
# Find 2×2 matrix with perm = 5
# perm([[a,b],[c,d]]) = ad + bc = 5
# E.g., [[1,2],[1,1]] has perm = 1 + 2 = 3. No.
# [[1,1],[2,3]] has perm = 3+2 = 5. Yes!
# [[2,1],[1,3]] has perm = 6+1 = 7. No.
# The simplest: [[F_2, F_3], [F_1, F_2]] = [[1,2],[1,1]], perm = 3 ≠ 5
# Or [[F_3, F_4], [F_1, F_2]] = [[2,3],[1,1]], perm = 2+3 = 5 ✓!

print("p=5, m=2: H/p = 5")
print("  perm([[F_3, F_4], [F_1, F_2]]) = perm([[2,3],[1,1]]) = 2+3 = 5 ✓")
print()

# For p=7, m=3: H/p = 25
# Find 3×3 matrix with perm = 25
# Try Fibonacci-based matrices

from itertools import permutations

def perm(M):
    """Permanent of matrix M."""
    n = len(M)
    total = 0
    for sigma in permutations(range(n)):
        prod = 1
        for i in range(n):
            prod *= M[i][sigma[i]]
        total += prod
    return total

# Try various 3×3 matrices built from Fibonacci numbers
print("p=7, m=3: H/p = 25, seeking 3×3 matrix with perm = 25")

# Try [[F_a, F_b, F_c], [F_d, F_e, F_f], [F_g, F_h, F_i]]
# where the indices relate to the tournament structure
best_match = None
for a in range(1, 8):
    for b in range(1, 8):
        for c in range(1, 8):
            for d in range(1, 8):
                for e in range(1, 8):
                    for f in range(1, 8):
                        M = [[fib(a), fib(b), fib(c)],
                             [fib(d), fib(e), fib(f)],
                             [fib(a+d-1), fib(b+e-1), fib(c+f-1)]]  # constrained
                        p_val = perm(M)
                        if p_val == 25:
                            if best_match is None:
                                best_match = (M, (a,b,c,d,e,f))

if best_match:
    M, idx = best_match
    print(f"  Found: perm({M}) = 25, Fib indices = {idx}")
else:
    print("  No simple Fibonacci matrix found with perm = 25")

# Let's try a simpler approach: binary matrices
print()
print("  Binary 3×3 matrices with perm = 25: none (max perm of 3×3 binary = 6)")
print("  So need non-binary entries.")
print()

# Actually, let's think about this differently.
# A p×p circulant matrix C_S has perm(C_S) counting doubly stochastic matchings.
# H counts Hamiltonian paths, which is a DIFFERENT combinatorial object.

# The key insight: H is the number of Hamiltonian paths in the tournament.
# For a circulant tournament, H = p · h where h counts HP starting at vertex 0.

# A HP starting at 0 is a sequence 0, v_1, v_2, ..., v_{p-1} where
# v_i - v_{i-1} ∈ S (mod p) and all v_i are distinct.

# This is the permanent of a DIFFERENT matrix:
# Let M[step, vertex] = 1 if vertex can be visited at this step
# This is a (p-1) × (p-1) matrix (or p-1 × p with a missing column)

# For small p, let's compute this step-reachability matrix
print("Step-reachability permanent for Interval tournament:")
for p in [5, 7]:
    m = (p-1)//2
    S = set(range(1, m+1))
    H = H_int[p]
    h = H // p

    # DP: compute which vertices can be at each step
    # Step 0: vertex 0
    # Step 1: vertices reachable from 0 via S
    # etc.

    # Actually, the "permanent" approach uses:
    # Construct bipartite graph G:
    # Left = {steps 1,...,p-1}, Right = {vertices 1,...,p-1}
    # Edge (step, vertex) if vertex could be visited at that step
    # (depends on ALL previous choices, so this isn't quite right)

    # The correct formulation: h = #HP from 0
    # = sum_{v1 ∈ S} sum_{v2: v2-v1 ∈ S, v2≠0,v1} ...
    # This is the permanent of the DIRECTED adjacency matrix
    # restricted to non-starting vertices.

    # More precisely: consider the (p-1)×(p-1) matrix A' where
    # A'[i,j] = 1 if edge from vertex i to vertex j in the tournament
    # (for vertices 1,...,p-1, with vertex 0 removed)
    # Then h ≠ perm(A'), because perm counts all perfect matchings
    # while we need Hamiltonian paths.

    # But: h = sum over permutations sigma of {1,...,p-1}
    # such that sigma is a Hamiltonian path from 0:
    # 0 → sigma(1) → sigma(2) → ... → sigma(p-1)
    # with sigma(k+1) - sigma(k) ∈ S mod p for each k,
    # and sigma(1) ∈ S.

    print(f"  p={p}: H = {H}, H/p = {h}")

print()

# ============================================================
# PART 8: THE RECIPROCAL POLYNOMIAL AND H
# ============================================================

print("=" * 70)
print("PART 8: RECIPROCAL POLYNOMIAL B_m*(x) AND H")
print("=" * 70)
print()

# B_m(x) = sum C(m+j,2j) x^j has reciprocal
# B_m*(x) = x^m B_m(1/x) = sum C(m+j,2j) x^{m-j}
# The roots of B_m(x) are at x = -1/Q_k (negative reciprocals of eigenvalues)
# The roots of B_m*(x) are at x = -Q_k

# prod(a + Q_k) = (-1)^m B_m*(-a) = a^m B_m(1/a)

# KEY: What is B_m*(2)?
# B_m*(2) = 2^m B_m(1/2)

for m in range(1, 9):
    p = 2*m + 1
    Bm2 = morgan_voyce_B(m, 2)
    Bm_half = morgan_voyce_B(m, 0.5)
    Bm_star_2 = 2**m * Bm_half

    Fp = fib(p)
    h = H_int.get(p, None)

    print(f"  m={m}, p={p}: B_m(2) = {Bm2}, B_m*(2) = 2^m·B_m(1/2) = {Bm_star_2:.0f}", end="")
    if h is not None:
        print(f", H/p = {h // p if isinstance(h, int) else h}", end="")
        print(f", H/(p·B_m(2)) = {(h)/(p*Bm2):.6f}" if h else "", end="")
    print()

print()

# ============================================================
# PART 9: SYNTHESIS — TOWARDS THE DETERMINANTAL FORMULA
# ============================================================

print("=" * 70)
print("SYNTHESIS: THE DETERMINANTAL FORMULA FOR H")
print("=" * 70)
print()
print("WHAT WE FOUND:")
print()
print("1. The Q_k spectrum determines H (for circulant tournaments)")
print("2. det(aI + diag(Q)) = a^m B_m(1/a) — parametric identity")
print("3. det(I + diag(Q)) = F_p — Fibonacci")
print("4. The Cassini identity is det([[B_m, B_{m-1}],[B_{m+1}, B_m]]) = -x")
print("5. H/p does NOT equal any simple determinant of diag(Q)")
print("6. H/p is NOT perm of a simple Fibonacci matrix")
print("7. The Hankel determinant of Q-moments = Vandermonde² (classical)")
print()
print("KEY INSIGHT: H depends on MORE than just the Q_k.")
print("H depends on the PHASES of the eigenvalues λ_k = -1/2 + i·D_k")
print("where D_k = Σ σ_i sin(2πk·a_i/p).")
print("The Q_k = 1/4 + D_k² lose the phase information.")
print()
print("For the Interval tournament, the phases are FIXED (all σ_i = +1).")
print("So the eigenvalues λ_k = -1/2 + i·D_k(+1,...,+1) are determined")
print("by the FULL eigenvalue spectrum, not just Q_k = |λ_k|².")
print()
print("CONCLUSION: To express H as a determinant, we need a matrix")
print("built from the FULL eigenvalues λ_k (including phases), not just Q_k.")
print()
print("det(zI - Λ) where Λ = diag(λ_1,...,λ_{p-1}) at z = ??? gives H?")
print()

# Check: det(I - A/z) for the tournament adjacency matrix at special z
for p in [5, 7]:
    m = (p-1)//2
    S = set(range(1, m+1))
    H_val = H_int[p]

    # Eigenvalues
    lambdas = []
    for k in range(p):
        lk = sum(np.exp(2j*np.pi*s*k/p) for s in S)
        lambdas.append(lk)

    lambdas = np.array(lambdas)
    print(f"p={p}: eigenvalues λ_k:")
    for k in range(p):
        print(f"  k={k}: λ = {lambdas[k]:.4f}, |λ|² = {abs(lambdas[k])**2:.4f}")

    # det(I - A) = prod(1 - λ_k) = characteristic polynomial at z=1
    det_ImA = np.prod(1 - lambdas)
    print(f"  det(I - A) = prod(1-λ_k) = {det_ImA:.4f}")
    print(f"  H = {H_val}")

    # Characteristic polynomial det(zI - A) evaluated at various z
    for z in [0, 1, 2, 3, -1, m, m+1]:
        char_val = np.prod(z - lambdas)
        print(f"  det({z}I - A) = {char_val:.4f}, H/det = {H_val/char_val.real:.4f}" if abs(char_val) > 0.01 else f"  det({z}I - A) = {char_val:.4f}")

    print()
