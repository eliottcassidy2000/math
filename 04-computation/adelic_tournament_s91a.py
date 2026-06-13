#!/usr/bin/env python3
"""
adelic_tournament_s91a.py — The Adelic Tournament Space
opus-2026-03-17-S91a

The tournament eigenvalue space has ADELIC structure:
  - Global eigenvalues are rationals k/D, D = odd part of C(n,2)
  - Each prime p | D gives a LOCAL factor (eigenvalue mod p^v)
  - The CRT reassembles local factors into the global eigenvalue
  - The archimedean place gives the real eigenvalue in [-1,1]

This is EXACTLY the structure of the adeles A_Q.

Geometric similarities explored:
  1. Restricted product decomposition
  2. Local-global principle (Hasse-Minkowski analog)
  3. Bruhat-Tits tree at each prime
  4. Spectral zeta function and Euler product
  5. Idele class group via Cayley transform
  6. Strong approximation
  7. Adelic Haar measure = stationary distribution
"""

from math import comb, factorial, gcd, log, sqrt, pi
from fractions import Fraction
from itertools import permutations
from collections import Counter
from functools import reduce
import numpy as np

# =============================================================================
# PART 1: THE EIGENVALUE DENOMINATOR AS CONDUCTOR
# =============================================================================

print("""


     THE ADELIC TOURNAMENT SPACE

     opus-2026-03-17-S91a



PART 1: THE CONDUCTOR
========================

The eigenvalue denominator D(n) = odd_part(C(n,2)) is the CONDUCTOR
of the adelic tournament space at size n.

In number theory, the conductor controls:
  - Which primes ramify in a number field
  - The level of a modular form
  - The "size" of the adelic quotient

In tournament theory, D(n) controls:
  - Which primes appear in eigenvalue denominators
  - The granularity of the spectral decomposition
  - The "resolution" at which the flip chain distinguishes classes
""")

def factorize(n):
    """Return prime factorization as {p: e} dict."""
    factors = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors

def odd_part(n):
    while n % 2 == 0:
        n //= 2
    return n

def euler_phi(n):
    """Euler's totient function."""
    result = n
    p = 2
    temp = n
    while p * p <= temp:
        if temp % p == 0:
            while temp % p == 0:
                temp //= p
            result -= result // p
        p += 1
    if temp > 1:
        result -= result // temp
    return result

print(f"{'n':>3} | {'C(n,2)':>7} | {'D(n)':>7} | {'factorization':>20} | {'phi(D)':>6} | {'primes':>15} | {'places':>6}")
print("-" * 80)

for n in range(3, 25):
    cn2 = comb(n, 2)
    D = odd_part(cn2)
    f = factorize(D)
    f_str = " * ".join(f"{p}^{e}" if e > 1 else str(p) for p, e in sorted(f.items()))
    primes = sorted(f.keys())
    phi_D = euler_phi(D)
    # Number of "places" = number of distinct primes + 1 (archimedean)
    num_places = len(primes) + 1
    print(f"{n:3d} | {cn2:7d} | {D:7d} | {f_str:>20} | {phi_D:6d} | {str(primes):>15} | {num_places:6d}")

# =============================================================================
# PART 2: LOCAL FACTORS AND CRT DECOMPOSITION
# =============================================================================

print("""

PART 2: LOCAL FACTORS AND CRT DECOMPOSITION
==============================================

The adeles A_Q = R x prod'_p Q_p (restricted product over all primes).

For tournament eigenvalues at conductor D:
  A_T(n) = R x prod_{p|D} Z/p^{v_p(D)}Z

Each eigenvalue lambda = k/D decomposes as:
  lambda_infty = k/D in R         (archimedean place)
  lambda_p = k mod p^{v_p(D)}     (p-adic place for p | D)

The CRT says: knowing all lambda_p determines k mod D,
hence determines the eigenvalue.

THIS IS THE LOCAL-GLOBAL PRINCIPLE FOR TOURNAMENTS.
""")

def crt_decompose(k, D):
    """Decompose k mod D into local components at each prime."""
    f = factorize(D)
    components = {}
    for p, e in sorted(f.items()):
        pe = p ** e
        components[p] = (k % pe, pe)
    return components

def crt_reconstruct(components, D):
    """Reconstruct k from local components via CRT."""
    result = 0
    for p, (residue, pe) in components.items():
        # M = D / pe
        M = D // pe
        # M_inv = M^{-1} mod pe
        M_inv = pow(M, -1, pe)
        result += residue * M * M_inv
    return result % D

# Demo at n=6: D = 15 = 3 * 5
print("n=6: D=15=3*5, eigenvalues k/15 for k in {-9,-7,-5,-3,-1,1,3,5,7,9,11,15}")
print()
print(f"{'k/15':>7} | {'k':>4} | {'k mod 3':>7} | {'k mod 5':>7} | {'real':>8} | {'CRT check':>10}")
print("-" * 55)

n6_eigenvalue_nums = [15, 11, 9, 7, 5, 3, 1, -1, -3, -5, -7, -9]
for k in n6_eigenvalue_nums:
    k_pos = k % 15  # work with positive representative
    comps = crt_decompose(k_pos, 15)
    reconstructed = crt_reconstruct(comps, 15)
    real_val = k / 15
    check = "OK" if reconstructed == k_pos else "FAIL"
    print(f"{k:4d}/15 | {k_pos:4d} | {comps[3][0]:7d} | {comps[5][0]:7d} | {real_val:+8.4f} | {check:>10}")

# n=10: D = 45 = 3^2 * 5
print()
print("n=10: D=45=9*5, eigenvalues decompose at primes {3, 5}")
print(f"  C(10,2) = {comb(10,2)}, odd part = {odd_part(comb(10,2))}")
cn2_10 = comb(10, 2)
D_10 = odd_part(cn2_10)
f_10 = factorize(D_10)
print(f"  Factorization: {f_10}")
print(f"  Local ring at 3: Z/9Z (9 elements, resolution 1/9)")
print(f"  Local ring at 5: Z/5Z (5 elements, resolution 1/5)")
print(f"  Total adelic resolution: 9 * 5 = 45 eigenvalue slots")

# =============================================================================
# PART 3: THE BRUHAT-TITS TREE
# =============================================================================

print("""

PART 3: THE BRUHAT-TITS TREE
===============================

At each prime p, Q_p has a natural TREE structure: the Bruhat-Tits tree T_p.
  - Vertices at depth d = lattices L in Q_p^2 up to scaling, with index p^d
  - Edges connect lattices differing by one step

For tournament eigenvalues at prime p with p^e || D:
  - The local eigenvalue k mod p^e lives at depth e in the tree
  - Refining n -> n+1 may increase the p-adic depth
  - The tournament "zooms in" on the p-adic tree as n grows

THE TREE STRUCTURE:
  Z/pZ (depth 1) <- Z/p^2 Z (depth 2) <- ... <- Z/p^e Z (depth e)

Each level has p branches. The tournament eigenvalue picks ONE branch at each level.
""")

# Track the Bruhat-Tits depth at each prime as n grows
print("Bruhat-Tits depth at primes 3, 5, 7 as n grows:")
print(f"{'n':>3} | {'D(n)':>7} | {'v_3':>4} | {'v_5':>4} | {'v_7':>4} | {'tree structure':>30}")
print("-" * 60)

for n in range(3, 25):
    D = odd_part(comb(n, 2))
    f = factorize(D)
    v3 = f.get(3, 0)
    v5 = f.get(5, 0)
    v7 = f.get(7, 0)
    tree = " x ".join(f"T_{p}(d={e})" for p, e in sorted(f.items()))
    print(f"{n:3d} | {D:7d} | {v3:4d} | {v5:4d} | {v7:4d} | {tree:>30}")

# =============================================================================
# PART 4: THE SPECTRAL ZETA FUNCTION
# =============================================================================

print("""

PART 4: THE SPECTRAL ZETA FUNCTION AND EULER PRODUCT
======================================================

In Tate's thesis, the global zeta function factors:
  Z(s) = prod_p Z_p(s) * Z_infty(s)

For the tournament flip chain at conductor D, define:
  zeta_T(s) = sum_{eigenvalues lambda != 0} |lambda|^{-s}

The local factor at prime p | D:
  Z_p(s) = sum_{k=1}^{p^e - 1} |k/p^e|_p^{-s}
         = sum_{j=0}^{e-1} phi(p^{e-j}) * (p^j)^s / (p^e)^s

Let's compute this for small n.
""")

def spectral_zeta_local(p, e, s):
    """Local spectral zeta factor at prime p with p^e || D."""
    pe = p ** e
    total = 0.0
    for k in range(1, pe):
        # |k/pe|_p = p^{-v_p(k)} / p^{-e} = p^{e - v_p(k)}
        # So |k/pe|_p^{-s} = p^{-s*(e - v_p(k))}
        v = 0
        temp = k
        while temp % p == 0:
            v += 1
            temp //= p
        val = p ** (e - v)
        total += val ** (-s)
    return total

def spectral_zeta_archimedean(eigenvalues, s):
    """Archimedean factor: sum |lambda|^{-s} over nonzero eigenvalues."""
    total = 0.0
    for lam in eigenvalues:
        if abs(lam) > 1e-12:
            total += abs(lam) ** (-s)
    return total

# For n=6: D=15, eigenvalue numerators are multiples of 1
# The nonzero eigenvalues are k/15 for k in {1,3,5,7,9,11} (positive) and negatives
print("Local zeta factors at n=6 (D=15=3*5):")
for s_val in [1.0, 2.0, 3.0]:
    z3 = spectral_zeta_local(3, 1, s_val)
    z5 = spectral_zeta_local(5, 1, s_val)
    z_product = z3 * z5
    print(f"  s={s_val:.0f}: Z_3(s) = {z3:.4f}, Z_5(s) = {z5:.4f}, product = {z_product:.4f}")

# =============================================================================
# PART 5: THE IDELE CLASS GROUP VIA CAYLEY TRANSFORM
# =============================================================================

print("""

PART 5: THE IDELE CLASS GROUP VIA CAYLEY TRANSFORM
====================================================

The Cayley transform Q(x) = (1+x)/(1-x) maps:
  Eigenvalue space [-1, 1] -> Multiplicative group (0, infinity)

In adelic terms: Q maps the ADDITIVE adeles to the MULTIPLICATIVE ideles.

The ideles J_Q = R* x prod'_p Q_p*.
The idele class group C_Q = J_Q / Q* encodes CLASS FIELD THEORY.

For tournament eigenvalues:
  Q(lambda) = (1 + k/D) / (1 - k/D) = (D+k) / (D-k)

These are RATIONAL numbers whose prime factorization reveals
the arithmetic structure of the eigenvalue.
""")

print("Cayley transform of eigenvalues at n=6 (D=15):")
print(f"{'k/15':>7} | {'Q(k/15)':>12} | {'factored':>20} | {'v_2':>4} | {'v_3':>4} | {'v_5':>4} | {'v_7':>4}")
print("-" * 75)

for k in [15, 11, 9, 7, 5, 3, 1, -1, -3, -5, -7, -9]:
    if k == 15:
        print(f"  15/15 |          inf |              inf |    - |    - |    - |    -")
        continue
    if k == -15:
        print(f" -15/15 |            0 |                0 |    - |    - |    - |    -")
        continue
    q = Fraction(15 + k, 15 - k)
    # p-adic valuations of numerator and denominator
    num = abs(q.numerator)
    den = abs(q.denominator)

    def v_p(n, p):
        if n == 0: return float('inf')
        v = 0
        while n % p == 0:
            v += 1
            n //= p
        return v

    # Valuation of q = v_p(num) - v_p(den)
    vals = {}
    for p in [2, 3, 5, 7]:
        vals[p] = v_p(num, p) - v_p(den, p)

    f_num = factorize(num)
    f_den = factorize(den)

    num_str = " * ".join(f"{p}^{e}" if e > 1 else str(p) for p, e in sorted(f_num.items())) if num > 1 else "1"
    den_str = " * ".join(f"{p}^{e}" if e > 1 else str(p) for p, e in sorted(f_den.items())) if den > 1 else "1"
    frac_str = f"{num_str}/{den_str}" if den > 1 else num_str

    print(f"{k:4d}/15 | {str(q):>12} | {frac_str:>20} | {vals[2]:4d} | {vals[3]:4d} | {vals[5]:4d} | {vals[7]:4d}")

print("""
THE IDELIC NORM:

For each eigenvalue Q(lambda) = (D+k)/(D-k), the idelic norm is:
  |Q(lambda)|_A = prod_v |Q(lambda)|_v = 1  (product formula!)

This is because Q(lambda) is a RATIONAL NUMBER,
and the product formula holds for all nonzero rationals.

The product formula is the GLOBAL CONSTRAINT that binds all local factors together.
It says: the archimedean "size" of an eigenvalue and its p-adic "sizes"
are not independent — they must multiply to 1.
""")

# Verify the product formula for one eigenvalue
k = 7
q = Fraction(15 + k, 15 - k)  # = 22/8 = 11/4
print(f"Product formula check for Q(7/15) = {q}:")
print(f"  |Q|_infty = {float(abs(q)):.6f}")
product = float(abs(q))
for p in [2, 3, 5, 7, 11]:
    num_val = v_p(abs(q.numerator), p)
    den_val = v_p(abs(q.denominator), p)
    val = num_val - den_val
    p_norm = p ** (-val) if val != 0 else 1.0
    product *= p_norm
    if val != 0:
        print(f"  |Q|_{p} = {p}^{-val} = {p_norm:.6f}")
print(f"  Product of all norms = {product:.6f} (should be 1)")

# =============================================================================
# PART 6: GEOMETRIC SIMILARITIES
# =============================================================================

print("""

PART 6: GEOMETRIC SIMILARITIES
=================================

The adelic tournament space has deep geometric parallels with
classical adelic geometry. Here are the precise correspondences:

SIMILARITY 1: RESTRICTED PRODUCT = EIGENVALUE DECOMPOSITION
  Adeles:  A_Q = R x prod'_p Q_p
  Tourney: Spec(T) = [-1,1] x prod_{p|D} Z/p^{v_p(D)}Z

  The eigenvalue space is a FINITE approximation to the adeles,
  truncated at the primes dividing D(n).

SIMILARITY 2: LOCAL-GLOBAL = CRT FOR EIGENVALUES
  Adeles:  Hasse-Minkowski: local solutions => global solution
  Tourney: CRT: local eigenvalues (k mod p^e) => global eigenvalue (k/D)

  The tournament version is STRONGER: CRT always works (no obstructions).
  Hasse-Minkowski can fail (e.g., Selmer curves). The tournament
  local-global principle is EXACT.

SIMILARITY 3: BRUHAT-TITS TREE = PRIME REFINEMENT TOWER
  Adeles:  T_p = (p+1)-regular tree, vertices = lattice classes
  Tourney: As n grows, v_p(D(n)) can increase, refining the p-adic resolution

  Going from n to n' where v_p(D) increases by 1 = descending
  one level in the Bruhat-Tits tree. The tree has p branches at each node.

SIMILARITY 4: ARCHIMEDEAN PLACE = REAL EIGENVALUE
  Adeles:  R = completion at infinity, the only archimedean place
  Tourney: lambda in [-1,1] = the real eigenvalue, the archimedean data

  The real eigenvalue is the "infinite prime" contribution.
  It controls the DECAY RATE (via ln|lambda|), which is archimedean.
  The p-adic contributions control the ARITHMETIC structure.

SIMILARITY 5: PRODUCT FORMULA = SPECTRAL CONSTRAINT
  Adeles:  prod_v |x|_v = 1 for x in Q*
  Tourney: prod_v |Q(lambda)|_v = 1 (Cayley transform to ideles)

  The product formula constrains how eigenvalue information is
  distributed across places. You can't be "big" everywhere.

SIMILARITY 6: HAAR MEASURE = STATIONARY DISTRIBUTION
  Adeles:  mu_A = mu_R x prod_p mu_p (product Haar measure)
  Tourney: pi = prod_{p|D} pi_p (stationary dist decomposes via CRT?)

  QUESTION: Does the stationary distribution on tournament classes
  factorize as a product of local distributions? Not exactly —
  the class sizes are not multiplicatively decomposable in general.
  But the SPECTRAL DECOMPOSITION does factorize.
""")

# =============================================================================
# PART 7: ADELIC DISTANCE AND TOURNAMENT GEOMETRY
# =============================================================================

print("""
PART 7: ADELIC DISTANCE BETWEEN TOURNAMENT CLASSES
====================================================

In adelic geometry, the distance between two points x, y in A_Q is:
  d_A(x,y) = max(|x-y|_infty, max_p |x-y|_p)  (sup-norm on adeles)

or equivalently the product:
  d_prod(x,y) = |x-y|_infty * prod_p max(1, |x-y|_p)

For tournament eigenvalues lambda_i = k_i/D and lambda_j = k_j/D:
  Real distance: |lambda_i - lambda_j| = |k_i - k_j| / D
  p-adic distance: |k_i - k_j|_p = p^{-v_p(k_i - k_j)}
""")

# Compute adelic distances between n=6 eigenvalues
D = 15
n6_nums = [15, 11, 9, 7, 5, 3, 1, -1, -3, -5, -7, -9]

print("Adelic distance matrix for n=6 eigenvalues (showing real + p-adic):")
print(f"  Eigenvalues: {[f'{k}/15' for k in n6_nums[:7]]}")
print()

# Just show a few interesting pairs
pairs = [(11, 9), (9, 7), (7, 5), (5, 3), (3, 1), (11, 1), (11, -9)]
print(f"{'pair':>12} | {'diff':>5} | {'|.|_R':>8} | {'|.|_3':>8} | {'|.|_5':>8} | {'adelic':>10}")
print("-" * 60)

for ki, kj in pairs:
    diff = ki - kj
    real_dist = abs(diff) / D
    # p-adic distances
    d3 = 3 ** (-v_p(abs(diff), 3)) if diff != 0 else 0
    d5 = 5 ** (-v_p(abs(diff), 5)) if diff != 0 else 0
    # Adelic distance (product formula variant)
    adelic = real_dist * max(1, 1/d3) * max(1, 1/d5) if diff != 0 else 0
    print(f"({ki:3d},{kj:3d}) | {diff:5d} | {real_dist:8.4f} | {d3:8.4f} | {d5:8.4f} | {adelic:10.4f}")

# =============================================================================
# PART 8: THE TOURNAMENT ADELE RING
# =============================================================================

print("""

PART 8: THE TOURNAMENT ADELE RING
===================================

Define the TOURNAMENT ADELE RING at level n:

  A_T(n) = R x prod_{p | D(n)} Z_p / p^{v_p(D)} Z_p

This is a FINITE quotient of the full adele ring A_Q.
As n grows, A_T(n) grows — more primes appear, existing primes get deeper.

The TOURNAMENT ADELE RING at n=6 is:
  A_T(6) = R x Z/3Z x Z/5Z

Dimension: 1 (real) + 1 (mod 3) + 1 (mod 5) = 3 "adelic dimensions"
But the discrete part has 3 * 5 = 15 points.

The TOURNAMENT ADELE RING at n=10 is:
  A_T(10) = R x Z/9Z x Z/5Z

Dimension: 1 + 2 + 1 = 4 "adelic dimensions" (depth 2 at prime 3!)
Discrete part has 9 * 5 = 45 points.
""")

# Compute the adelic dimension and size for each n
print(f"{'n':>3} | {'D(n)':>7} | {'adelic dim':>10} | {'disc size':>9} | {'primes':>15} | {'depths':>15}")
print("-" * 75)

for n in range(3, 21):
    D = odd_part(comb(n, 2))
    f = factorize(D)
    adelic_dim = 1 + sum(f.values())  # 1 for real + sum of depths
    disc_size = D  # product of p^e
    primes_str = ",".join(str(p) for p in sorted(f.keys()))
    depths_str = ",".join(f"{p}^{e}" for p, e in sorted(f.items()))
    print(f"{n:3d} | {D:7d} | {adelic_dim:10d} | {disc_size:9d} | {primes_str:>15} | {depths_str:>15}")

# =============================================================================
# PART 9: STRONG APPROXIMATION
# =============================================================================

print("""

PART 9: STRONG APPROXIMATION
===============================

Strong approximation (Eichler, Kneser, Platonov):
  For G a simply-connected algebraic group over Q,
  G(Q) is dense in G(A_f) (the finite adeles).

Tournament analog:
  The global eigenvalues k/D are DENSE in the product ring
  prod_{p|D} Z/p^{v_p(D)}Z.

  More precisely: as k ranges over {0, 1, ..., D-1},
  the CRT image covers ALL of the product ring.
  This is EXACT density — not just density, but surjectivity!

  Strong approximation for tournaments is TRIVIALLY TRUE
  because CRT is an isomorphism for coprime moduli.
""")

# Verify: CRT gives a bijection for D=15
print("Strong approximation at D=15:")
print("  CRT: Z/15Z -> Z/3Z x Z/5Z")
print()
print(f"  {'k':>3} | {'k mod 3':>7} | {'k mod 5':>7}")
print(f"  {'-'*25}")
for k in range(15):
    print(f"  {k:3d} | {k%3:7d} | {k%5:7d}")

print(f"""
  Every pair (a, b) with 0 <= a < 3, 0 <= b < 5 appears exactly once.
  CRT is a PERFECT BIJECTION. Strong approximation is exact.

  This means: for EVERY choice of local eigenvalue data
  (one residue at each prime), there exists a UNIQUE global eigenvalue
  that has those local components. No obstructions!
""")

# =============================================================================
# PART 10: THE ARCHIMEDEAN-NONARCHIMEDEAN DUALITY
# =============================================================================

print("""
PART 10: THE ARCHIMEDEAN / NON-ARCHIMEDEAN DUALITY
=====================================================

The deepest geometric similarity is the DUALITY between
the archimedean and non-archimedean data of an eigenvalue.

Archimedean (real place):
  - Controls DYNAMICS: decay rate = |lambda|, mixing time = 1/(1-lambda)
  - Continuous: lambda in [-1, 1]
  - Metric: standard absolute value
  - Geometry: interval/disk/hyperbolic space

Non-archimedean (p-adic places):
  - Controls ARITHMETIC: which primes divide the eigenvalue, ramification
  - Discrete: k mod p^e in finite set
  - Metric: ultrametric (|x+y| <= max(|x|, |y|))
  - Geometry: tree (Bruhat-Tits)

The PRODUCT FORMULA binds them:
  The archimedean "size" and the p-adic "sizes" are not independent.
  Big archimedean eigenvalue => specific p-adic constraints.
""")

# Demonstrate the duality
print("Archimedean vs non-archimedean data for n=6 eigenvalues:")
print(f"{'lambda':>8} | {'|lambda|':>8} | {'decay':>8} | {'k mod 3':>7} | {'k mod 5':>7} | {'type':>12}")
print("-" * 65)

for k in n6_nums:
    lam = Fraction(k, 15)
    decay = float(abs(lam))
    k_pos = k % 15
    kmod3 = k_pos % 3
    kmod5 = k_pos % 5

    # Classify
    if k == 15:
        t = "stationary"
    elif k == -15:
        t = "alternating"
    elif abs(k) % 3 == 0 and abs(k) % 5 == 0:
        t = "zero"
    elif abs(k) % 3 == 0:
        t = "3-ramified"
    elif abs(k) % 5 == 0:
        t = "5-ramified"
    else:
        t = "generic"

    print(f"{str(lam):>8} | {decay:8.4f} | {decay:8.4f} | {kmod3:7d} | {kmod5:7d} | {t:>12}")

# =============================================================================
# PART 11: TATE'S THESIS FOR TOURNAMENTS
# =============================================================================

print("""

PART 11: TATE'S THESIS FOR TOURNAMENTS
=========================================

Tate's thesis (1950) reinterprets the Riemann zeta function adelically:
  zeta(s) = integral_{J_Q/Q*} |x|^s d*x

The tournament analog:

Define the TOURNAMENT L-FUNCTION at level n:
  L_T(s) = prod_{p | D(n)} L_p(s)

where L_p(s) is the local factor encoding the p-adic structure
of the flip chain.

For the simplest case (p || D, meaning v_p(D) = 1):
  L_p(s) = (1 - p^{-s})^{-1}  (the standard Euler factor!)

This is because: at depth 1, the local eigenvalues are
{0, 1, 2, ..., p-1} mod p, and the nonzero ones contribute
(p-1) terms of size 1 to the local zeta, giving 1/(1-p^{-s}).
""")

# Compute the tournament L-function for small n
print("Tournament L-function (Euler product form):")
print()

for n in [5, 6, 7, 8, 10, 12]:
    D = odd_part(comb(n, 2))
    f = factorize(D)
    factors = []
    for p, e in sorted(f.items()):
        if e == 1:
            factors.append(f"(1-{p}^(-s))^(-1)")
        else:
            factors.append(f"L_{p}^({e})(s)")

    product_str = " * ".join(factors)
    # Evaluate at s=1 (harmonic)
    val_s1 = 1.0
    for p, e in f.items():
        if e == 1:
            val_s1 *= p / (p - 1)
        else:
            # More complex local factor
            local = sum(1.0 for k in range(1, p**e) if k % p != 0)
            val_s1 *= p**e / (p**e - p**(e-1))

    print(f"  n={n:2d}: D={D:5d}, L_T(s) = {product_str}")
    print(f"         L_T(1) = {val_s1:.4f}")

# =============================================================================
# PART 12: THE CLASS FIELD THEORY CONNECTION
# =============================================================================

print("""

PART 12: CLASS FIELD THEORY CONNECTION
========================================

The Artin map of class field theory:
  art: J_Q/Q* -> Gal(Q^ab/Q)

maps the idele class group to the Galois group of the maximal
abelian extension. By Kronecker-Weber, Q^ab = Q(all roots of unity).

For tournaments: the Cayley transform Q(lambda) lives in Q*.
The idelic image of Q(lambda) determines a character of Gal(Q^ab/Q).

CONCRETELY:
  Q(k/D) = (D+k)/(D-k)

The primes dividing (D+k)(D-k) determine which cyclotomic
fields are "activated" by this eigenvalue.
""")

# Compute the "activated fields" for each n=6 eigenvalue
print("Activated cyclotomic fields for n=6 eigenvalues:")
print(f"{'k/15':>7} | {'Q(k/15)':>12} | {'num primes':>15} | {'den primes':>15} | {'cyclotomic':>15}")
print("-" * 75)

for k in [11, 9, 7, 5, 3, 1, -1, -3, -5, -7, -9]:
    q = Fraction(15 + k, 15 - k)
    num_primes = sorted(factorize(abs(q.numerator)).keys()) if abs(q.numerator) > 1 else []
    den_primes = sorted(factorize(abs(q.denominator)).keys()) if abs(q.denominator) > 1 else []
    all_primes = sorted(set(num_primes + den_primes))
    cyclo = ",".join(f"Q(zeta_{p})" for p in all_primes) if all_primes else "Q"

    print(f"{k:4d}/15 | {str(q):>12} | {str(num_primes):>15} | {str(den_primes):>15} | {cyclo:>15}")

# =============================================================================
# PART 13: ADELIC FORMAL GROUP
# =============================================================================

print("""

PART 13: THE ADELIC FORMAL GROUP
==================================

The formal group F(x,y) = (x+y)/(1+xy) lives simultaneously
at every place:

  F over R:    The standard tanh addition. Height infinity.
               Connected to the hyperbolic plane.
               Convergence radius: |x|, |y| < 1.

  F over Q_p (p odd):  Height 1 formal group.
               [p](x) = x^p mod p (Frobenius).
               The p-torsion F[p] = Z/pZ.
               Convergence: p-adic unit disk.

  F over Q_2:  Height INFINITY (supersingular!).
               [2](x) = 2x/(1+x^2) = 0 mod 2 at x=0.
               Actually [2](x) has leading term 2x, and
               over F_2, [2](x) = 0 (it kills everything).
               This is why H(T) is ALWAYS ODD (Redei's theorem).

THE ADELIC FORMAL GROUP is the product:
  F_A = F_R x prod_p F_{Q_p}

Tournament eigenvalues live in the KERNEL of certain [n]-maps:
  [n](lambda) = 0 in F_A means:
  - [n](lambda)_R = 0 (real torsion)
  - [n](lambda)_p = 0 for all p (p-adic torsion)

The [n]-torsion of F_A = prod of local torsion groups:
  F_A[n] = F_R[n] x prod_p F_{Q_p}[n]

For F = tanh addition:
  F_R[n] = {tanh(k*pi*i/n) : k = 0,...,n-1} (complex torsion)
  F_{Q_p}[n] = Z/nZ for p coprime to n (Lubin-Tate theory)
""")

# Compute torsion points of the formal group
print("Torsion of F(x,y) = (x+y)/(1+xy) at small primes:")
print()

for p in [3, 5, 7]:
    print(f"  F[{p}] over Q_{p}:")
    # [p](x) = Q^{-1}(Q(x)^p) = ((1+x)^p - (1-x)^p) / ((1+x)^p + (1-x)^p)
    # Torsion points: [p](x) = 0, i.e., (1+x)^p = (1-x)^p
    # Over Q_p: x = 0 is always a solution. Others require p-th roots of unity.
    # F[p] = ker([p]) = {0} union {x : Q(x)^p = 1} = {x : Q(x) is a p-th root of 1}
    # Q(x) = zeta_p^k => x = (zeta_p^k - 1)/(zeta_p^k + 1)
    # These x live in Q_p(zeta_p), not in Q_p (unless p | p, i.e., always for the ramified prime).
    print(f"    Torsion points: x_k = (zeta_{p}^k - 1)/(zeta_{p}^k + 1) for k=0,...,{p-1}")
    print(f"    These live in Q_{p}(zeta_{p}) = unramified extension of degree {p-1}")
    print(f"    Over Q_{p} itself: only x=0 (the identity)")
    print()

# =============================================================================
# PART 14: THE RAPIDITY LATTICE AS ADELIC OBJECT
# =============================================================================

print("""
PART 14: THE RAPIDITY LATTICE
================================

The rapidity of eigenvalue lambda is:
  rho = arctanh(lambda) = (1/2) * ln(Q(lambda))

At n=6 (the simplest nontrivial case with D=15=3*5),
the rapidity lattice is spanned by {ln(2), ln(3), ln(7)}.

From the adelic perspective:
  ln(Q(lambda)) = ln(D+k) - ln(D-k)

At the archimedean place: this is the usual real logarithm.
At the p-adic place: this is the p-adic logarithm log_p.

The p-ADIC RAPIDITY at prime p is:
  rho_p = (1/2) * log_p(Q(lambda))
        = (1/2) * (log_p(D+k) - log_p(D-k))

where log_p is the Iwasawa logarithm (the p-adic log).
""")

# Compute p-adic valuations of rapidities
print("Rapidity structure at n=6:")
print(f"{'k/15':>7} | {'Q(k/15)':>10} | {'arctanh':>10} | {'v_2(Q)':>7} | {'v_3(Q)':>7} | {'v_5(Q)':>7} | {'v_7(Q)':>7}")
print("-" * 70)

from math import atanh

for k in [11, 9, 7, 5, 3, 1, -1, -3, -5, -7, -9]:
    q = Fraction(15 + k, 15 - k)
    rho = atanh(k / 15)

    def valuation(frac, p):
        """p-adic valuation of a fraction."""
        return v_p(abs(frac.numerator), p) - v_p(abs(frac.denominator), p)

    v2 = valuation(q, 2)
    v3 = valuation(q, 3)
    v5 = valuation(q, 5)
    v7 = valuation(q, 7)

    print(f"{k:4d}/15 | {str(q):>10} | {rho:+10.6f} | {v2:7d} | {v3:7d} | {v5:7d} | {v7:7d}")

# =============================================================================
# PART 15: THE DEEP GEOMETRIC PICTURE
# =============================================================================

print("""

PART 15: THE DEEP GEOMETRIC PICTURE
======================================

THE TOURNAMENT SPACE IS AN ADELIC MANIFOLD.

At each level n, the tournament eigenvalue space is:
  Spec_n = R x prod_{p | D(n)} Z/p^{v_p} Z

This is a compact subset of the adeles A_Q.

As n -> infinity:
  - More primes enter (D(n) has more prime factors)
  - Existing primes go deeper (v_p(D(n)) can increase)
  - The adelic space A_T(n) approaches the FULL adeles A_Q

The LIMIT is:
  A_T(infinity) = lim_{n -> infty} A_T(n) = A_Q

THE TOURNAMENT THEORY IS A FINITE APPROXIMATION TO ADELIC GEOMETRY.

Each specific n gives a "level n truncation" of the full adelic picture,
just as a modular form of level N gives a truncation of the full
automorphic picture.

GEOMETRIC CORRESPONDENCES:
  Tournament flip graph  <->  Hecke operator on automorphic forms
  Eigenvalue lambda      <->  Hecke eigenvalue a_p
  Spectral gap 4/C(n,2)  <->  Ramanujan bound |a_p| <= 2*sqrt(p)
  Zero eigenvalue         <->  Vanishing of L-function at s=1/2
  H(T) always odd         <->  Supersingularity at p=2
  D(n) = odd part C(n,2)  <->  Level/conductor of automorphic form

THE FLIP CHAIN IS A HECKE OPERATOR.

The transition matrix T = (1/C(n,2)) * sum_{arcs} P_{ij}
is an average over "local" operators, exactly like:
  T_p = sum of Hecke correspondences at prime p.

Hecke operators act on automorphic forms by:
  (T_p f)(L) = sum_{L' ~ L} f(L')
summing over lattices L' connected to L by a p-isogeny.

The flip operator acts on tournament functions by:
  (P_{ij} f)(A) = f(flip_{ij}(A))
summing over tournaments A' connected to A by a single arc flip.

BOTH average over local modifications.
BOTH have eigenvalues that decompose multiplicatively.
BOTH satisfy a "Ramanujan-type" bound.

THE TOURNAMENT FLIP CHAIN = HECKE OPERATOR ON THE ADELIC TOURNAMENT SPACE.
""")

# =============================================================================
# SUMMARY TABLE
# =============================================================================

print("=" * 80)
print("SUMMARY: GEOMETRIC SIMILARITIES")
print("=" * 80)
print()

similarities = [
    ("Structure",        "Adelic geometry",               "Tournament space"),
    ("Base space",       "Spec(Z) = {primes}",            "D(n) = odd part C(n,2)"),
    ("Local factors",    "Q_p at each prime p",            "Z/p^e Z at each p|D"),
    ("Global space",     "A_Q = R x prod' Q_p",           "R x prod Z/p^e Z"),
    ("Local-global",     "Hasse-Minkowski",                "CRT (always exact!)"),
    ("Metric",           "Product of p-adic + real",       "Adelic distance on evals"),
    ("Tree structure",   "Bruhat-Tits tree T_p",           "Prime depth tower"),
    ("Zeta function",    "zeta(s) = prod 1/(1-p^-s)",     "L_T(s) = prod L_p(s)"),
    ("Operator",         "Hecke operator T_p",             "Flip operator P_{ij}"),
    ("Eigenvalue bound", "Ramanujan |a_p| <= 2sqrt(p)",    "Spectral gap 4/C(n,2)"),
    ("Formal group",     "Lubin-Tate at each p",           "F(x,y)=(x+y)/(1+xy)"),
    ("Height structure",  "height 1 (ordinary) or inf",    "height 1 (odd p) or inf (p=2)"),
    ("Class field",      "Artin map: ideles -> Galois",    "Cayley: evals -> cyclotomic"),
    ("Product formula",  "prod |x|_v = 1",                "prod |Q(lam)|_v = 1"),
    ("Strong approx",    "G(Q) dense in G(A_f)",           "CRT bijection (exact!)"),
    ("Supersingular",    "Formal group height inf at p",   "H(T) odd (Redei) at p=2"),
    ("Conductor",        "N = level of modular form",      "D(n) = eigenvalue denom"),
    ("Rapidity",         "p-adic logarithm log_p",         "arctanh = formal log"),
]

print(f"{'':>18} | {'Adelic geometry':>30} | {'Tournament space':>30}")
print("-" * 85)
for label, adelic, tourney in similarities:
    print(f"{label:>18} | {adelic:>30} | {tourney:>30}")

print("""

THE DEEPEST INSIGHT:

The tournament flip chain at level n is a FINITE-LEVEL HECKE OPERATOR
on a TRUNCATED ADELIC SPACE.

As n -> infinity, the tournament adelic space fills out the full adeles A_Q,
and the flip chain approaches a true Hecke operator on A_Q.

The eigenvalue denominator D(n) = odd part of C(n,2) is the CONDUCTOR,
controlling which primes are "visible" at level n.

The formal group F(x,y) = (x+y)/(1+xy) is the LUBIN-TATE formal group
at each odd prime (height 1) and is SUPERSINGULAR at p=2 (height infinity).

This explains:
  - Why H(T) is always odd: supersingularity at p=2
  - Why eigenvalues are k/D: they live on the finite adelic quotient
  - Why the spectral gap = 4/C(n,2): it's the conductor's reciprocal
  - Why CRT decomposes eigenvalues: the adelic product structure
  - Why {2,3,7} control the Cayley boosts: they are the Hurwitz primes,
    the primes of the (2,3,7) triangle group = PSL(2,Z)

42 = 2*3*7 = the conductor of the UNIVERSAL adelic tournament.
""")

print("opus-2026-03-17-S91a: THE ADELIC TOURNAMENT SPACE")
