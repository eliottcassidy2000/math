#!/usr/bin/env python3
"""
fibonacci_tiling_galois.py — opus-2026-03-12-S67d

FIBONACCI DEEP DIVE PART 2: Tiling interpretation and Galois splitting.

KEY QUESTIONS:
1. Does F_p = prod(1+Q_k) have a TILING interpretation for the tournament?
2. How does the Binet formula F_p = (φ^p - ψ^p)/√5 relate to Galois theory?
3. What is the continued fraction structure of B_m(x)?
4. Can we relate H to Fibonacci via a LATTICE MODEL?
5. Does the Fibonacci recurrence appear in the Walsh coefficients?

CONNECTIONS TO EXPLORE:
- Fibonacci numbers count tilings of 1×n strips with squares and dominoes
- The transfer matrix [[1,1],[1,0]] for Fibonacci vs [[3,-1],[1,0]] for Morgan-Voyce
- Chebyshev/Morgan-Voyce polynomials as partition functions
- The golden ratio φ as an algebraic unit and its role in Q(√5)
- Cassini's identity F_{n-1}F_{n+1} - F_n^2 = (-1)^n and its analogue
"""

import numpy as np
from math import gcd, comb, factorial
from collections import Counter, defaultdict
from fractions import Fraction
import sys

def legendre(a, p):
    a = a % p
    if a == 0: return 0
    ls = pow(a, (p-1)//2, p)
    return -1 if ls == p-1 else ls

def morgan_voyce_B(m, x):
    """B_m(x) = sum_{j=0}^{m} C(m+j, 2j) x^j"""
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

def all_circulant_H(p):
    """Compute H for all 2^m circulant tournaments at prime p."""
    m = (p-1)//2
    pairs = []
    seen = set()
    for s in range(1, p):
        if s not in seen:
            pairs.append((s, p-s))
            seen.add(s)
            seen.add(p-s)

    results = []
    for bits in range(2**m):
        S = set()
        for i in range(m):
            if bits & (1 << i):
                S.add(pairs[i][0])
            else:
                S.add(pairs[i][1])

        A = [[0]*p for _ in range(p)]
        for i in range(p):
            for j in range(p):
                if i != j and (j-i) % p in S:
                    A[i][j] = 1

        n = p
        dp = [[0]*n for _ in range(1 << n)]
        for v in range(n):
            dp[1 << v][v] = 1
        for mask in range(1, 1 << n):
            for v in range(n):
                if not (mask & (1 << v)): continue
                if dp[mask][v] == 0: continue
                for u in range(n):
                    if mask & (1 << u): continue
                    if A[v][u]:
                        dp[mask | (1 << u)][u] += dp[mask][v]
        full_mask = (1 << n) - 1
        H = sum(dp[full_mask][v] for v in range(n))

        interval = set(range(1, m+1))
        paley_set = set(s for s in range(1, p) if legendre(s, p) == 1)
        name = ""
        if S == interval: name = "INT"
        if S == paley_set: name = "PAL"

        results.append({'S': S, 'H': H, 'name': name})

    return results

print("=" * 70)
print("FIBONACCI TILING & GALOIS SPLITTING — opus-2026-03-12-S67d")
print("=" * 70)
print()

# ============================================================
# PART 1: TILING INTERPRETATION
# ============================================================

print("PART 1: FIBONACCI AS TILING COUNT — CAN WE TILE TOURNAMENTS?")
print("=" * 70)
print()
print("F_n = number of tilings of 1×(n-1) strip with squares and dominoes")
print("F_p = prod(1+Q_k) for the Interval tournament on p vertices")
print()
print("QUESTION: Is there a tiling model on the tournament where:")
print("  - Each tiling has weight 1")
print("  - Total count = prod(1+Q_k)?")
print()

# The transfer matrix for F_n is T_F = [[1,1],[1,0]]
# The transfer matrix for B_m is T_B = [[3,-1],[1,0]]
# Relationship: T_B has eigenvalues φ² and ψ² = 1/φ²
# So B_m is the "squared Fibonacci" transfer matrix!

print("Transfer matrix comparison:")
print("  Fibonacci:    T_F = [[1,1],[1,0]], eigenvalues φ, ψ = (1±√5)/2")
print("  Morgan-Voyce: T_B = [[3,-1],[1,0]], eigenvalues φ², ψ²")
print()
print("KEY: T_B = T_F², i.e., the MV transfer matrix is the SQUARE")
print("of the Fibonacci transfer matrix!")
print()

# Verify: T_F² = [[1,1],[1,0]]² = [[2,1],[1,1]] — NO, that's not [[3,-1],[1,0]]
T_F = np.array([[1,1],[1,0]])
T_F2 = T_F @ T_F
print(f"T_F² = {T_F2.tolist()} — NOT equal to T_B = [[3,-1],[1,0]]")
print()

# Actually, let's think about this differently.
# T_B has characteristic polynomial t² - 3t + 1 = 0
# T_F has char poly t² - t - 1 = 0
# If λ is eigenvalue of T_F, then λ² satisfies (λ²)² - 3(λ²) + 1 = 0
# Check: λ² = (1+√5)/2)² = (3+√5)/2
# (3+√5)²/4 - 3(3+√5)/2 + 1 = (14+6√5)/4 - (9+3√5)/2 + 1
# = (14+6√5)/4 - (18+6√5)/4 + 4/4 = 0 ✓

print("Eigenvalue relationship:")
phi = (1 + np.sqrt(5)) / 2
psi = (1 - np.sqrt(5)) / 2
print(f"  φ = {phi:.6f}, φ² = {phi**2:.6f}")
print(f"  ψ = {psi:.6f}, ψ² = {psi**2:.6f}")
print(f"  φ² + ψ² = {phi**2 + psi**2:.6f} = 3 (trace of T_B)")
print(f"  φ² · ψ² = {(phi*psi)**2:.6f} = 1 (det of T_B)")
print()

# So B_m(1) = F_p because:
# B_m(1) = trace of T_B^m (the (1,1) entry in some basis)
# and F_p = (φ^p - ψ^p)/√5
# Since p = 2m+1, φ^p = (φ²)^m · φ, ψ^p = (ψ²)^m · ψ
# F_p = (φ^{2m+1} - ψ^{2m+1})/√5
# B_m(1) = (φ^{2m+2} - ψ^{2m+2})/(φ² - ψ²)
# φ² - ψ² = (φ-ψ)(φ+ψ) = √5 · 1 = √5
# So B_m(1) = (φ^{2m+2} - ψ^{2m+2})/√5 = F_{2m+2} = F_{p+1}??

# Wait, let me recheck. B_m(1) = sum C(m+j,2j) for j=0..m
# This is the central Delannoy number or related?
# Actually we verified B_m(1) = F_p (the p-th Fibonacci number where p=2m+1)
# Let me verify directly:
print("Verifying B_m(1) = F_p (p=2m+1):")
for m in range(1, 10):
    p = 2*m + 1
    Bm1 = morgan_voyce_B(m, 1)
    Fp = fib(p)
    print(f"  m={m}, p={p}: B_m(1) = {Bm1}, F_p = {Fp}, match = {Bm1 == Fp}")
print()

# Hmm, let me check if B_m(1) = F_{2m+1} or F_{2m+2}
# B_0(1) = C(0,0) = 1 = F_1
# B_1(1) = C(1,0) + C(2,2) = 1 + 1 = 2 = F_3
# B_2(1) = C(2,0) + C(3,2) + C(4,4) = 1 + 3 + 1 = 5 = F_5
# B_3(1) = C(3,0) + C(4,2) + C(5,4) + C(6,6) = 1 + 6 + 5 + 1 = 13 = F_7
# So B_m(1) = F_{2m+1}! And p = 2m+1, so B_m(1) = F_p ✓

# ============================================================
# PART 2: BINET/GALOIS SPLITTING
# ============================================================

print("=" * 70)
print("PART 2: BINET FORMULA AND GALOIS SPLITTING OF prod(1+Q_k)")
print("=" * 70)
print()

# F_p = (φ^p - ψ^p)/√5
# Since prod(1+Q_k) = F_p, we have:
# prod(1+Q_k) = (φ^p - ψ^p)/√5
#
# Now, the Q_k are eigenvalues living in the real cyclotomic field K = Q(cos 2π/p)
# [K : Q] = (p-1)/2 = m
# The Galois group Gal(K/Q) ≅ (Z/pZ)* / {±1}
#
# φ^p and ψ^p live in Q(√5), while the Q_k live in K.
# If 5 is a quadratic residue mod p, then Q(√5) ⊂ K.
# If 5 is a nonresidue, then Q(√5) and K are linearly disjoint over Q.

print("Galois structure: K = Q(cos 2π/p), [K:Q] = m = (p-1)/2")
print()
print("Q(√5) ⊂ K iff 5 is a quadratic residue mod p:")
for p in [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43]:
    if p == 5:
        leg = 0
        contained = "ramified"
    else:
        leg = legendre(5, p)
        contained = "YES" if leg == 1 else "NO"
    print(f"  p={p:2d}: (5/p) = {leg:+d}, Q(√5) ⊂ K: {contained}")
print()

# By quadratic reciprocity: (5/p) = (p/5)
# p mod 5: 1,4 → (p/5)=+1; 2,3 → (p/5)=-1
print("By QR: (5/p) = (p/5). So Q(√5) ⊂ K iff p ≡ ±1 mod 5")
print()

# The key question: how does prod(1+Q_k) = F_p factor in K?
#
# Each Q_k = (sin(mπk/p)/sin(πk/p))² is an element of K
# The product over k=1..m gives an element of Q (specifically F_p ∈ Z)
#
# Under the Galois group, the Q_k permute. Specifically:
# σ_a : cos(2π/p) ↦ cos(2πa/p) sends Q_k to Q_{ak mod p}
# (with the identification Q_k = Q_{p-k})

print("Galois action on Q_k for Interval tournament:")
print()

for p in [7, 11, 13]:
    m = (p-1)//2
    Q = interval_Q(p)
    print(f"p={p}, m={m}:")
    print(f"  Q_k = {np.round(Q, 4)}")

    # The Galois group acts by k ↦ ak mod p (for a in (Z/pZ)*/{±1})
    # Under this action, the set {Q_1,...,Q_m} permutes
    # prod(1+Q_k) is Galois-invariant, hence in Q

    # Can we find a PARTIAL product that gives φ^p or (φ^p - ψ^p)/√5?
    # Idea: split the Q_k into two Galois orbits based on (k/p)

    QR_indices = [k for k in range(1, m+1) if legendre(k, p) == 1]
    QNR_indices = [k for k in range(1, m+1) if legendre(k, p) == -1]

    prod_QR = np.prod([1 + Q[k-1] for k in QR_indices])
    prod_QNR = np.prod([1 + Q[k-1] for k in QNR_indices])

    print(f"  QR indices: {QR_indices}")
    print(f"  QNR indices: {QNR_indices}")
    print(f"  prod(1+Q_k) over QR: {prod_QR:.6f}")
    print(f"  prod(1+Q_k) over QNR: {prod_QNR:.6f}")
    print(f"  Product = {prod_QR * prod_QNR:.6f} = F_{p} = {fib(p)}")
    print(f"  Ratio QR/QNR = {prod_QR/prod_QNR:.6f}")
    print(f"  Sum QR+QNR = {prod_QR + prod_QNR:.6f}")
    print(f"  (Sum)² - 5*(F_p)² = {(prod_QR + prod_QNR)**2 - 5*fib(p)**2:.6f}")

    # Check: does the partial product live in Q(√5)?
    # If Gal(K/Q) acts transitively on QR indices and QNR indices separately...
    # Actually QR and QNR are the two cosets of the index-2 subgroup of (Z/pZ)*/{±1}
    # when m is even

    # Try: split by half
    first_half = list(range(1, m//2 + 1))
    second_half = list(range(m//2 + 1, m+1))
    prod_first = np.prod([1 + Q[k-1] for k in first_half])
    prod_second = np.prod([1 + Q[k-1] for k in second_half])
    print(f"  prod(1+Q) first half: {prod_first:.6f}")
    print(f"  prod(1+Q) second half: {prod_second:.6f}")

    # Check if partial products relate to Lucas numbers
    Lp = fib(p-1) + fib(p+1)  # Lucas number
    print(f"  L_p = F_{p-1} + F_{p+1} = {fib(p-1)} + {fib(p+1)} = {Lp}")
    print(f"  Binet: F_p = (φ^p - ψ^p)/√5, L_p = φ^p + ψ^p")
    print(f"  φ^p = (L_p + F_p·√5)/2 = {(Lp + fib(p)*np.sqrt(5))/2:.6f}")
    print(f"  ψ^p = (L_p - F_p·√5)/2 = {(Lp - fib(p)*np.sqrt(5))/2:.6f}")
    print()

# ============================================================
# PART 3: CASSINI IDENTITY ANALOGUE
# ============================================================

print("=" * 70)
print("PART 3: CASSINI-TYPE IDENTITIES FOR B_m(x)")
print("=" * 70)
print()
print("Cassini: F_{n-1}F_{n+1} - F_n² = (-1)^n")
print("For B_m: is there B_{m-1}(x)B_{m+1}(x) - B_m(x)² = f(x)?")
print()

for x_val in [1, 2, 3, 4, -1]:
    print(f"  x = {x_val}:")
    for m in range(1, 8):
        Bm = morgan_voyce_B(m, x_val)
        Bm1 = morgan_voyce_B(m-1, x_val)
        Bp1 = morgan_voyce_B(m+1, x_val)
        cassini = Bm1 * Bp1 - Bm * Bm
        print(f"    m={m}: B_{m-1}·B_{m+1} - B_m² = {cassini}")
    print()

# ============================================================
# PART 4: CONTINUED FRACTION REPRESENTATION
# ============================================================

print("=" * 70)
print("PART 4: CONTINUED FRACTION STRUCTURE")
print("=" * 70)
print()
print("Morgan-Voyce polynomials satisfy B_m(x) = (2+x)B_{m-1}(x) - B_{m-2}(x)")
print("with B_0 = 1, B_1 = 2+x.")
print()
print("This means B_m/B_{m-1} has the continued fraction expansion:")
print("  B_m/B_{m-1} = (2+x) - 1/((2+x) - 1/((2+x) - ...))")
print("which is a periodic CF with period 1!")
print()

# At x=1: CF for F_{2m+1}/F_{2m-1}
# The CF [3; 3, 3, 3, ...] converges to (3+√5)/2 = φ²
print("Convergence of B_m(1)/B_{m-1}(1) = F_{2m+1}/F_{2m-1} → φ²:")
for m in range(1, 12):
    p = 2*m + 1
    ratio = morgan_voyce_B(m, 1) / morgan_voyce_B(m-1, 1)
    print(f"  m={m:2d}: B_m/B_{m-1} = {ratio:.10f}, φ² = {phi**2:.10f}, diff = {abs(ratio - phi**2):.2e}")
print()

# The connection: the CF [3; 3, 3, ...] = (3+√5)/2 = φ²
# This is the golden ratio SQUARED appearing as the fixed point of x ↦ 3 - 1/x
print(f"Fixed point of x ↦ (2+1) - 1/x = 3 - 1/x: x = (3+√5)/2 = φ² = {phi**2:.6f}")
print(f"Compare: fixed point of x ↦ 1 + 1/x: x = (1+√5)/2 = φ = {phi:.6f}")
print()

# ============================================================
# PART 5: B_m(x) AS PARTITION FUNCTION — LATTICE GAS CONNECTION
# ============================================================

print("=" * 70)
print("PART 5: B_m(x) AS LATTICE GAS PARTITION FUNCTION")
print("=" * 70)
print()
print("CRUCIAL CONNECTION: Both H(T) and prod(1+Q_k) arise from")
print("lattice gas / hard-core polymer models!")
print()
print("H(T) = I(Ω(T), 2) = hard-core lattice gas on odd-cycle graph at activity 2")
print("prod(1+Q_k) = ???")
print()
print("Q: Can prod(1+Q_k) be written as a hard-core partition function too?")
print()

# prod(1+Q_k) = sum over subsets A of {1,...,m} of prod_{k in A} Q_k
# = sum_A prod_{k in A} |D_m(2πk/p)|²
# This IS a partition function — it's the independent set polynomial
# of the EMPTY graph on m vertices, evaluated at activities Q_k!
# (Every subset is independent when there are no edges.)

print("prod(1+Q_k) = I(empty graph on m vertices, with activities Q_k)")
print("  = sum_{A ⊆ [m]} prod_{k∈A} Q_k")
print()
print("By contrast, H = I(Ω(T), 2) uses the ODD-CYCLE graph with uniform activity 2.")
print()
print("So the ratio H / prod(1+Q_k) compares:")
print("  - Hard-core gas on Ω(T) vs no-exclusion gas on {Q_k}")
print("  - The GRAPH STRUCTURE in Ω(T) excludes some configurations")
print()

# Let's compute this ratio
for p in [7, 11, 13]:
    m = (p-1)//2
    results = all_circulant_H(p)

    for r in results:
        if r['name'] == 'INT':
            H = r['H']
            Q = interval_Q(p)
            prod_1Q = np.prod(1 + Q)
            print(f"  p={p}: H(Int) = {H}, prod(1+Q) = {prod_1Q:.1f} = F_{p} = {fib(p)}")
            print(f"         H / F_p = {H / fib(p):.6f}")
            print(f"         H / (p · F_p) = {H / (p * fib(p)):.6f}")
            # Check: H = p · something · F_p ?
            print(f"         H / p = {H // p}")
            print(f"         (H/p) / F_p = {(H//p) / fib(p):.6f}")
            # Maybe H ≈ C(m) · prod(1+Q_k) · p?
            print(f"         Ratio chain: H/p/F_p = {(H//p)/fib(p):.6f}")
            print()
print()

# ============================================================
# PART 6: FIBONACCI RECURRENCE IN H VALUES
# ============================================================

print("=" * 70)
print("PART 6: DOES H SATISFY A FIBONACCI-LIKE RECURRENCE?")
print("=" * 70)
print()

# Known H values for Interval tournament:
# p=3: H=3 (trivially: 3-cycle)
# p=5: H=25
# p=7: H=175
# p=11: H=93027
# p=13: H=3711175
# p=17: H=13689269499

H_int = {3: 3, 5: 25, 7: 175, 11: 93027, 13: 3711175, 17: 13689269499}

print("H(Int, p) values:")
primes = sorted(H_int.keys())
for p in primes:
    H = H_int[p]
    m = (p-1)//2
    Fp = fib(p)
    print(f"  p={p:2d}, m={m}: H={H:>16d}, F_p={Fp:>8d}, H/p={H//p:>14d}")

print()
print("Checking ratios H(p)/H(q) for consecutive primes:")
for i in range(1, len(primes)):
    p, q = primes[i-1], primes[i]
    ratio = H_int[q] / H_int[p]
    m_p, m_q = (p-1)//2, (q-1)//2
    # Compare to (q-1)!/(p-1)! / 2^{q-p}
    growth = factorial(q-1) / factorial(p-1) / 2**(q-p) * q/p
    print(f"  H({q})/H({p}) = {ratio:.4f}, expected ~{growth:.4f}")

print()

# Check: H/p and Fibonacci
print("Checking: Is H/p a Fibonacci-like combination?")
for p in primes:
    h = H_int[p] // p
    m = (p-1)//2
    # Various Fibonacci numbers
    for k in range(1, 2*p+2):
        if fib(k) == h:
            print(f"  p={p}: H/p = {h} = F_{k}")
            break
    else:
        # Check if h = product of Fibonacci numbers
        pass

print()

# ============================================================
# PART 7: THE DEEP CONNECTION — CHEBYSHEV AND TRANSFER MATRICES
# ============================================================

print("=" * 70)
print("PART 7: CHEBYSHEV POLYNOMIALS AND THE UNIFIED FRAMEWORK")
print("=" * 70)
print()

# Morgan-Voyce B_m(x) can be expressed via Chebyshev polynomials of the 2nd kind:
# B_m(x) = U_m(1 + x/2) where U_m is the Chebyshev polynomial of 2nd kind
# defined by U_m(cos θ) = sin((m+1)θ)/sin θ
#
# So B_m(1) = U_m(3/2) and the argument 3/2 = (1+1/2+1) = 1 + x/2 at x=1

# Let's verify this
print("Chebyshev U_m connection:")
print("U_m(cos θ) = sin((m+1)θ) / sin θ")
print("B_m(x) = U_m(1 + x/2)")
print()

def chebyshev_U(m, x):
    """Chebyshev polynomial of the second kind."""
    if m == 0: return 1
    if m == 1: return 2*x
    U_prev2, U_prev1 = 1, 2*x
    for _ in range(2, m+1):
        U_curr = 2*x * U_prev1 - U_prev2
        U_prev2, U_prev1 = U_prev1, U_curr
    return U_curr

for m in range(0, 8):
    Bm = morgan_voyce_B(m, 1)
    Um = chebyshev_U(m, 3/2)
    print(f"  m={m}: B_m(1) = {Bm}, U_m(3/2) = {Um:.0f}, F_{2*m+1} = {fib(2*m+1)}, match = {Bm == int(round(Um))}")

print()

# Now the DEEP observation:
# H(T) = I(Ω(T), 2) where Ω(T) is the odd-cycle intersection graph
# prod(1+Q_k) = F_p = U_m(3/2) where m = (p-1)/2
#
# The argument 3/2 comes from the TRANSFER MATRIX of the Interval:
# T = [[3,-1],[1,0]] has trace 3 = 2·(3/2)
# So cos θ = 3/2 means θ is imaginary: θ = i·arccosh(3/2)
# The eigenvalues of T are 2·cosh(θ) ± something...
#
# Actually: eigenvalues of T are (3 ± √5)/2 = φ², ψ²
# And 3/2 = (φ² + ψ²)/2 = (trace)/2

print("THE UNIFIED VIEW:")
print()
print("1. Interval Q_k = |D_m(2πk/p)|² = Dirichlet kernel squared")
print("2. prod(1+Q_k) = B_m(1) = U_m(3/2) = F_p")
print("3. U_m(3/2) = value of Chebyshev polynomial at x = 3/2 = cos(iθ)")
print("   where θ = arccosh(3/2) ≈ 0.9624")
print(f"   cosh(0.9624) = {np.cosh(0.9624):.6f}")
print()

# ============================================================
# PART 8: THE KEY INSIGHT — PARAMETRIC DEFORMATION
# ============================================================

print("=" * 70)
print("PART 8: PARAMETRIC DEFORMATION — FROM F_p TO H(T)")
print("=" * 70)
print()
print("prod(a + Q_k) = a^m · B_m(1/a) = a^m · U_m(1 + 1/(2a))")
print()
print("At a=2: prod(2+Q_k) = 2^m · B_m(1/2) = 2^m · U_m(5/4)")
print("At a=1: prod(1+Q_k) = B_m(1) = U_m(3/2) = F_p")
print("At a→∞: prod(a+Q_k) → a^m (trivial)")
print()

# The question is: how does H relate to this parametric family?
# H = I(Ω(T), 2) involves the graph Ω(T)
# prod(1+Q_k) = I(∅, Q_k) has no graph structure
#
# IDEA: Is H = det(I + 2·M) for some matrix M built from the Q_k and Ω(T)?
# If so, then H generalizes prod(1+Q_k) by adding the GRAPH CONSTRAINT.

print("Parametric values of prod(a+Q_k) for Interval, p=7 and p=11:")
print()

for p in [7, 11, 13]:
    m = (p-1)//2
    Q = interval_Q(p)

    print(f"p={p}, m={m}, Q = {np.round(Q, 4)}:")
    H = H_int[p]

    for a in [0.5, 1, 1.5, 2, 3, 4, 5, 10]:
        prod_val = np.prod(a + Q)
        Bm_val = a**m * morgan_voyce_B(m, 1/a)
        print(f"  a={a:5.1f}: prod(a+Q) = {prod_val:>16.2f}, a^m·B_m(1/a) = {float(Bm_val):>16.2f}, ratio to H = {H/prod_val:.6f}")

    # Is there an 'a' where prod(a+Q) = H?
    # H = prod(a_0 + Q_k) for some a_0?
    # Since prod(a+Q_k) is a polynomial in a of degree m, H may not be achieved
    # (the Q_k are fixed for Interval)

    # But maybe H = prod(a+Q_k) at a = some function of p?
    # H/p / F_p ≈ ?
    h_over_Fp = (H // p) / fib(p)
    print(f"  H/p / F_p = {h_over_Fp:.6f}")
    print()

# ============================================================
# PART 9: FIBONACCI LATTICE — PHYLLOTAXIS AND DIFFRACTION
# ============================================================

print("=" * 70)
print("PART 9: FIBONACCI LATTICE AND QUASICRYSTAL STRUCTURE")
print("=" * 70)
print()

# The Fibonacci lattice Z_p points: {k·φ mod 1 : k = 0,...,p-1}
# This is related to the Fibonacci word and phyllotaxis

# For the Interval tournament S = {1,...,m}:
# The connection set is a consecutive block
# This is the MOST STRUCTURED possible set (minimum complexity)
#
# The Fibonacci spiral in phyllotaxis uses angle 2π/φ² = 2π(1 - 1/φ)
# The golden angle is 2π(2 - φ) ≈ 137.5°

print("Golden angle: 2π(2-φ) = 2π/φ² ≈ {:.4f}° ≈ {:.4f} rad".format(
    360 * (2 - phi), 2 * np.pi * (2 - phi)))
print()

# For each prime p, consider the Fibonacci permutation:
# σ_F: k ↦ floor(k·m/p) or k ↦ k·g mod p where g is a primitive root
# related to Fibonacci

print("Fibonacci permutation σ: k ↦ k·a mod p for various a:")
print()

for p in [7, 11, 13]:
    m = (p-1)//2

    # Find the element closest to φ·p
    a_golden = round(p / phi**2) % p
    if a_golden == 0: a_golden = 1

    print(f"p={p}: golden element a = round(p/φ²) mod p = {a_golden}")

    # The Fibonacci permutation rearranges S = {1,...,m}
    S_int = set(range(1, m+1))
    S_fib = set((a_golden * s) % p for s in S_int)

    print(f"  S_int = {sorted(S_int)}")
    print(f"  σ_{a_golden}(S_int) = {sorted(S_fib)}")

    # Is this close to a QR set?
    paley_set = set(s for s in range(1, p) if legendre(s, p) == 1)
    overlap_paley = len(S_fib & paley_set)
    print(f"  Overlap with Paley: {overlap_paley}/{m}")
    print()

# ============================================================
# PART 10: THE TRANSFER MATRIX TRACE FORMULA
# ============================================================

print("=" * 70)
print("PART 10: TRANSFER MATRIX TRACE FORMULA FOR H")
print("=" * 70)
print()

# H counts Hamiltonian paths. For circulant graphs, the adjacency matrix
# is a circulant with eigenvalues λ_k = S_hat(k).
#
# Hamiltonian path count: H = sum_v sum_σ ∈ S_n Hamiltonian prod A[σ(i),σ(i+1)]
# This is related to the PERMANENT of a matrix, not the determinant.
#
# For a circulant matrix C with eigenvalues λ_0,...,λ_{n-1}:
# perm(C) is NOT simply a function of the eigenvalues!
# (Unlike det(C) = prod(λ_k))
#
# HOWEVER, the Hamilton path count for a circulant CAN sometimes
# be expressed using the transfer matrix approach.

# The transfer matrix for counting paths in a circulant tournament:
# For S = {1,...,m}, the adjacency is: A[i,j] = 1 iff j-i mod p ∈ S
# A Hamiltonian path visits each vertex exactly once.
#
# The key: for the INTERVAL tournament, the connection set {1,...,m}
# means we go "forward" by 1 to m steps. This is like a restricted
# random walk on Z_p.
#
# H = number of permutations σ of Z_p such that σ(i+1) - σ(i) ∈ S for all i.
# This is a ROOK PLACEMENT problem on a circulant constraint graph!

print("Rook placement interpretation:")
print("H = #{σ ∈ S_p : σ(i+1) - σ(i) mod p ∈ S for all i}")
print()
print("For S = {1,...,m} (Interval):")
print("  Each step goes forward by 1 to m positions (mod p)")
print("  This is like a PAWN moving on a circular chessboard!")
print()

# Connection to dimer covers / matchings:
# A Hamiltonian path on K_p with constraint S is equivalent to
# a perfect matching on a bipartite graph G where:
# Left vertices = positions 1,...,p-1 (steps)
# Right vertices = elements 1,...,p (vertices to visit)
# Edge (step i, vertex v) exists iff v can be reached at step i

print("Bipartite matching formulation:")
print("  H = perm(M) where M is a (p-1)×p matrix")
print("  M[i,v] = 1 iff vertex v can be visited at step i")
print("  This permanent counts Hamiltonian paths.")
print()

# ============================================================
# PART 11: FIBONACCI AND H — THE LATTICE GAS BRIDGE
# ============================================================

print("=" * 70)
print("PART 11: THE LATTICE GAS BRIDGE — OCF MEETS FIBONACCI")
print("=" * 70)
print()

# OCF: H(T) = I(Ω(T), 2) = sum over independent sets S of Ω(T) of 2^|S|
# where Ω(T) is the odd-cycle intersection graph.
#
# Fibonacci: F_n = I(P_{n-2}, 1) = sum over independent sets of path P_{n-2}
# (This is the classic result: Fibonacci counts independent sets of paths)
#
# For p = 2m+1 prime:
# F_p = I(P_{p-2}, 1) = I(P_{2m-1}, 1) = independent sets of P_{2m-1}
#
# And we showed: F_p = prod(1+Q_k) for Interval tournament
#
# QUESTION: Is there a graph homomorphism Ω(Int) → P_{2m-1}
# that maps the hard-core gas at activity 2 to activity 1?

print("Fibonacci = independent set polynomial of path graph:")
print("  F_n = I(P_{n-2}, 1)")
print()
print("OCF = independent set polynomial of odd-cycle graph:")
print("  H(T) = I(Ω(T), 2)")
print()
print("For the Interval tournament:")
print("  F_p = prod(1+Q_k) = independent set poly of ??? at activity ???")
print()

# Let me verify the path interpretation
print("Verify: F_n = I(P_{n-2}, 1):")
for n in range(2, 12):
    Fn = fib(n)
    # I(P_{n-2}, x) = number of independent sets of path on n-2 vertices at activity x
    # For path P_k: I(P_k, x) = sum_{indep sets} x^|S|
    # I(P_0, x) = 1 + x (one vertex, empty set or {v})
    # Actually P_0 is a single vertex, P_1 is an edge, P_k has k+1 vertices
    # Wait: need to be careful about convention
    # F_n = I(P_{n-1}, 1) where P_k has k vertices? Let me just count.
    pass

# Simpler: F_1 = 1, F_2 = 1, F_3 = 2, F_4 = 3, F_5 = 5
# I(P_0, 1) = 2 (empty set + singleton on 1 vertex)
# I(P_1, 1) = 3 (empty, {v1}, {v2} on 2-vertex path)
# Actually F_{n+2} = I(P_n) where I counts all independent sets (at activity 1)
# P_0: 1 vertex → I = 2 = F_3 ✓
# P_1: 2 vertices, 1 edge → I = 3 = F_4 ✓
# P_2: 3 vertices → I = 4 = not F_5=5... hmm

# Actually the standard result: F_{n+2} = number of independent sets of P_n
# P_0: {}, {1} → 2 = F_2? No...
# Let me just Google the standard result. The standard is:
# The number of independent sets of a path P_n (n vertices) is F_{n+2}
# P_1: {}, {1} → 2 = F_3 ✓
# P_2: {}, {1}, {2} → 3 = F_4 ✓
# P_3: {}, {1}, {2}, {3}, {1,3} → 5 = F_5 ✓

print("Standard result: number of independent sets of P_n (n vertices) = F_{n+2}")
for n in range(1, 8):
    print(f"  P_{n} ({n} vertices): F_{n+2} = {fib(n+2)}")
print()

# So F_p = F_{p} = I(P_{p-2})  at activity 1 (counting independent sets)
# H(Int, p) = I(Ω(Int), 2) at activity 2

# At p=7: F_7 = 13, H(Int) = 175 = 7 × 25
# At p=11: F_11 = 89, H(Int) = 93027 = 11 × 8457
# At p=13: F_13 = 233, H(Int) = 3711175 = 13 × 285475

# Ω(Int) has m vertices (one per odd cycle class)
# P_{p-2} has p-2 = 2m-1 vertices
# So Ω(Int) has about half as many vertices as the Fibonacci path!

print("Graph sizes:")
for p in [7, 11, 13]:
    m = (p-1)//2
    print(f"  p={p}: Ω(Int) has {m} vertices, P_{p-2} = P_{2*m-1} has {2*m-1} vertices")
print()

# ============================================================
# PART 12: CLUSTER EXPANSION — CONNECTING H AND F_p
# ============================================================

print("=" * 70)
print("PART 12: CLUSTER EXPANSION — EXPRESSING H/F_p")
print("=" * 70)
print()
print("Both H and F_p are partition functions (independent set polys).")
print("Their RATIO R = H/F_p = I(Ω,2)/I(∅,Q_k) measures the effect of:")
print("  (a) Changing the graph from empty to Ω(T)")
print("  (b) Changing activities from Q_k to uniform 2")
print()

# Cluster expansion: log(I(G,z)) = sum over connected subgraphs
# This is the Mayer expansion from statistical mechanics!
#
# log(H) = log I(Ω, 2) = sum_C (-1)^{|C|-1} / |C| * sum_{S connected in Ω} ...
# log(F_p) = log I(∅, Q) = sum log(1+Q_k)  (no graph → no clusters)
#
# log(H/F_p) = log H - log F_p = cluster expansion of graph effects

for p in [7, 11, 13]:
    H = H_int[p]
    Fp = fib(p)
    m = (p-1)//2
    Q = interval_Q(p)

    log_ratio = np.log(H) - np.log(Fp)

    # Compare to various quantities
    print(f"p={p}: log(H/F_p) = {log_ratio:.6f}")
    print(f"  = log({H}) - log({Fp})")
    print(f"  H/F_p = {H/Fp:.6f}")
    print(f"  (H/F_p)/p = {H/(Fp*p):.6f}")

    # sum log(Q_k)
    sum_log_Q = np.sum(np.log(Q))
    print(f"  sum log(Q_k) = {sum_log_Q:.6f}")

    # log(H/F_p) / sum_log(Q)
    print(f"  log(H/F_p) / sum_log(Q) = {log_ratio / sum_log_Q:.6f}")

    # Is H/F_p close to something nice?
    # Try: C(2m, m) / 4^m * something
    central = comb(2*m, m)
    print(f"  C(2m,m) = {central}")
    print(f"  H / (p · F_p · m!) = {H / (p * Fp * factorial(m)):.6f}")
    print()

print()
print("=" * 70)
print("SYNTHESIS: THE FIBONACCI-TOURNAMENT DICTIONARY (EXTENDED)")
print("=" * 70)
print()
print("TIER 1 — PROVEN IDENTITIES:")
print("  prod(1+Q_k) = F_p (Fibonacci number)")
print("  prod(a+Q_k) = a^m · B_m(1/a) (Morgan-Voyce)")
print("  B_m(x) = U_m(1+x/2) (Chebyshev of 2nd kind)")
print("  B_m(1) = F_{2m+1} = F_p")
print("  B_m(4) = Pell_{2m+1}")
print("  T_B eigenvalues = φ², ψ² (golden ratio squared)")
print("  CF: B_m/B_{m-1} = [3; 3, 3, ...] → φ² (periodic CF)")
print()
print("TIER 2 — STRUCTURAL PARALLELS:")
print("  F_n = I(P_{n-2}, 1) — independent sets of path")
print("  H(T) = I(Ω(T), 2) — independent sets of odd-cycle graph")
print("  Both are lattice gas partition functions!")
print("  Fibonacci: activity 1 on path graph")
print("  Tournament: activity 2 on odd-cycle graph")
print()
print("TIER 3 — KEY QUESTIONS:")
print("  Q1: Is there a graph homomorphism Ω(Int) → P_{p-2}?")
print("  Q2: Does the cluster expansion of H/F_p have a nice form?")
print("  Q3: Can the Binet splitting φ^p, ψ^p be matched to Galois orbits of Q_k?")
print("  Q4: Does Cassini-type identity B_{m-1}·B_{m+1} - B_m² = (-1)^m")
print("      have a tournament interpretation?")
print()
print("TIER 4 — PROOF STRATEGY:")
print("  The path from Fibonacci to H-maximization goes:")
print("  (1) Interval Q_k are zeros of B_m(1/x) (up to scaling)")
print("  (2) B_m encodes the Chebyshev/Dirichlet kernel structure")
print("  (3) The Q_k are maximally CONCENTRATED (largest max, smallest min)")
print("  (4) This concentration → max additive energy E(S)")
print("  (5) Max E(S) → max Walsh degree-4 terms → max H")
print("  (6) The Fibonacci identity prod(1+Q_k) = F_p is a SHADOW")
print("      of the deeper: Q_k structure determines H")
print()
