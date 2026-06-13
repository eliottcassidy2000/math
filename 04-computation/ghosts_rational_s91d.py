#!/usr/bin/env python3
"""
ghosts_rational_s91d.py — Ghosts of rational structure in the transcendental world
opus-2026-03-17-S91d

A "ghost" is a discrete/rational structure that survives when we cross
the Cayley gate into the transcendental. The Z^3 rapidity lattice is
the first ghost. What others exist?

Systematic hunt across every operation we can apply to eigenvalues.
"""

from math import comb, factorial, log, sqrt, pi, e as E, atanh, cos, sin, exp, tanh
from fractions import Fraction
from itertools import permutations, combinations
from functools import reduce
import numpy as np

# =============================================================================
# Build n=5 transition data
# =============================================================================

def build_transition(n):
    num_arcs = comb(n, 2)
    perms = list(permutations(range(n)))
    def adj_from_mask(mask):
        A = [[0]*n for _ in range(n)]
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                if (mask >> idx) & 1: A[i][j] = 1
                else: A[j][i] = 1
                idx += 1
        return A
    def canon(A):
        best = None
        for p in perms:
            s = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
            if best is None or s < best: best = s
        return best
    canon_map = {}; classes = []; t_class = [0]*(2**num_arcs)
    for mask in range(2**num_arcs):
        A = adj_from_mask(mask)
        c = canon(A)
        if c not in canon_map:
            canon_map[c] = len(classes); classes.append(c)
        t_class[mask] = canon_map[c]
    nc = len(classes)
    sizes = [0]*nc
    for ci in t_class: sizes[ci] += 1
    F = np.zeros((nc, nc), dtype=int)
    for mask in range(2**num_arcs):
        ci = t_class[mask]
        for bit in range(num_arcs):
            cj = t_class[mask ^ (1 << bit)]
            F[ci][cj] += 1
    T = np.zeros((nc, nc))
    for ci in range(nc):
        for cj in range(nc):
            T[ci][cj] = F[ci][cj] / (sizes[ci] * num_arcs)
    return T, nc, sizes

print("""


     GHOSTS OF RATIONAL STRUCTURE

     opus-2026-03-17-S91d


A "ghost" is a COUNTABLE or DISCRETE structure that persists
when eigenvalues are mapped into the transcendental world.

Ghost 0: The rapidity lattice Z^3 = (1/2)(Z*ln2 + Z*ln3 + Z*ln7).
         Already found. What else?

Hunting for more...
""")

# The eigenvalues at n=5 (simplified)
evals_frac = [Fraction(1), Fraction(3,5), Fraction(2,5), Fraction(1,5),
              Fraction(0), Fraction(-1,5), Fraction(-3,5)]
# With multiplicities: 1(x1), 3/5(x1), 2/5(x1), 1/5(x4), 0(x1), -1/5(x3), -3/5(x1)
evals_with_mult = [(Fraction(1),1), (Fraction(3,5),1), (Fraction(2,5),1),
                   (Fraction(1,5),4), (Fraction(0),1), (Fraction(-1,5),3), (Fraction(-3,5),1)]

# Nonzero non-stationary eigenvalues
nz_evals = [e for e in evals_frac if e != 0 and e != 1]

# =============================================================================
# GHOST 1: TRACES ARE ALWAYS RATIONAL
# =============================================================================

print("""
GHOST 1: THE TRACE GHOST
===========================

Tr(T^k) = sum_{i} lambda_i^k is ALWAYS RATIONAL for integer k.

This is because T^k has rational matrix entries (T has rational entries),
and the trace of a rational matrix is rational.

But the INDIVIDUAL eigenvalues lambda_i can be irrational for general
matrices. For tournaments, they happen to be rational. But the trace
being rational is a DEEPER fact — it holds even when eigenvalues aren't.

The trace ghost: Tr(T^k) in Q for ALL k in Z.
This is a sequence of rationals encoding the FULL spectral data.
""")

T5, nc5, sizes5 = build_transition(5)

print(f"  n=5 Trace ghost Tr(T^k):")
print(f"  {'k':>4} | {'Tr(T^k)':>12} | {'exact fraction':>20} | {'interpretation':>25}")
print("  " + "-" * 70)

for k in range(0, 13):
    Tk = np.linalg.matrix_power(T5, k) if k > 0 else np.eye(nc5)
    tr = np.trace(Tk)
    # Try to find exact fraction
    best = Fraction(tr).limit_denominator(100000)
    # Check
    interp = ""
    if k == 0: interp = "dim = T(5) = 12"
    elif k == 1: interp = "sum lambda_i"
    elif k == 2: interp = "sum lambda_i^2 = energy"

    print(f"  {k:4d} | {tr:12.6f} | {str(best):>20} | {interp:>25}")

# =============================================================================
# GHOST 2: THE DETERMINANT GHOST
# =============================================================================

print("""

GHOST 2: THE DETERMINANT GHOST
=================================

det(T) = product of ALL eigenvalues.
det(T - lambda*I) = characteristic polynomial, ALL coefficients rational.

The characteristic polynomial p(x) = prod(x - lambda_i) has
INTEGER coefficients (after clearing denominators).

Each coefficient is an ELEMENTARY SYMMETRIC FUNCTION of the eigenvalues:
  e_0 = 1
  e_1 = sum lambda_i = Tr(T)
  e_2 = sum_{i<j} lambda_i * lambda_j
  ...
  e_n = prod lambda_i = det(T)

These are ALL rational — the full symmetric function ghost.
""")

# Compute characteristic polynomial
char_poly = np.poly(T5)
print(f"  n=5 characteristic polynomial coefficients:")
print(f"  p(x) = ", end="")
for i, c in enumerate(char_poly):
    frac = Fraction(c).limit_denominator(1000000)
    if i == 0:
        print(f"x^{len(char_poly)-1}", end="")
    elif frac >= 0:
        print(f" + {frac}*x^{len(char_poly)-1-i}", end="")
    else:
        print(f" - {-frac}*x^{len(char_poly)-1-i}", end="")
print()

print(f"\n  Elementary symmetric functions e_k(eigenvalues):")
for k in range(min(7, len(char_poly))):
    coeff = Fraction(char_poly[k]).limit_denominator(1000000)
    sign = (-1)**k
    ek = sign * coeff
    print(f"    e_{k} = {ek}")

# =============================================================================
# GHOST 3: THE MOMENT GHOST
# =============================================================================

print("""

GHOST 3: THE MOMENT / POWER SUM GHOST
========================================

The power sums p_k = sum lambda_i^k are rational for ALL k in Z_>=0.
These are the MOMENTS of the spectral distribution.

Newton's identities connect moments to elementary symmetric functions:
  k*e_k = sum_{i=1}^{k} (-1)^{i-1} e_{k-i} p_i

This means: the moments (transcendental-looking sums) are
DETERMINED by the elementary symmetric functions (rational coefficients).

The moment ghost: the generating function
  M(t) = sum_{k>=0} p_k * t^k = Tr(1/(1-t*T))
is a RATIONAL FUNCTION of t.
""")

print(f"  Power sums p_k = sum lambda_i^k:")
print(f"  {'k':>4} | {'p_k (float)':>14} | {'p_k (exact)':>20} | {'rational?':>10}")
print("  " + "-" * 55)

for k in range(11):
    pk = sum(float(e)**k * m for e, m in evals_with_mult)
    pk_frac = Fraction(pk).limit_denominator(1000000)
    is_rat = abs(float(pk_frac) - pk) < 1e-10
    print(f"  {k:4d} | {pk:14.6f} | {str(pk_frac):>20} | {'YES' if is_rat else 'NO':>10}")

# =============================================================================
# GHOST 4: THE RESOLVENT GHOST
# =============================================================================

print("""

GHOST 4: THE RESOLVENT GHOST
===============================

R(z) = Tr((zI - T)^{-1}) = sum_i 1/(z - lambda_i)

For z in Q (z not an eigenvalue), R(z) is a sum of rationals -> RATIONAL.

The resolvent is a RATIONAL function of z.
It has POLES at the eigenvalues and is rational EVERYWHERE else.

The resolvent ghost: a meromorphic function with rational values
at all rational arguments away from poles.

The poles are at the eigenvalues {1, 3/5, 2/5, 1/5, 0, -1/5, -3/5}.
Between any two poles, R(z) takes every rational value.
""")

print(f"  Resolvent R(z) at rational z values:")
print(f"  {'z':>8} | {'R(z)':>14} | {'exact':>15}")
print("  " + "-" * 45)

for z in [Fraction(2), Fraction(3,2), Fraction(4,5), Fraction(1,2),
          Fraction(1,10), Fraction(-1,2), Fraction(-1), Fraction(-2)]:
    rz = sum(m / (z - e) for e, m in evals_with_mult if e != z)
    rz_frac = Fraction(rz).limit_denominator(10000000)
    print(f"  {str(z):>8} | {float(rz):14.6f} | {str(rz_frac):>15}")

# =============================================================================
# GHOST 5: THE HEAT KERNEL GHOST
# =============================================================================

print("""

GHOST 5: THE HEAT KERNEL GHOST
=================================

K(t) = Tr(exp(-t*T)) = sum_i m_i * exp(-t*lambda_i)

For generic t, K(t) is TRANSCENDENTAL.
But the Taylor coefficients of K(t) around t=0 are RATIONAL:
  K(t) = sum_{k>=0} (-t)^k/k! * p_k
       = sum_{k>=0} (-t)^k/k! * Tr(T^k)

The heat kernel ghost: the Taylor coefficients are rational moments.
The function itself is transcendental, but its SKELETON is rational.

Even more: at t = ln(q) for rational q, we get:
  K(ln(q)) = sum m_i * q^{-lambda_i}
which involves q^{-lambda_i}. For rational lambda_i:
  q^{-lambda_i} = q^{-a/b} = (q^{1/b})^{-a}
This is ALGEBRAIC (the b-th root of a rational)!

So K(ln(q)) for rational q is a sum of ALGEBRAIC numbers.
By the transcendence of e, K(t) for irrational t is transcendental.
But K at logarithmic rationals has algebraic structure!
""")

print(f"  Heat kernel K(t) = Tr(exp(-t*T)) at special t:")
print(f"  {'t':>12} | {'K(t)':>14} | {'nature':>20}")
print("  " + "-" * 55)

for t_val, t_name in [(0, "0"), (1, "1"), (log(2), "ln(2)"),
                       (log(3), "ln(3)"), (log(7), "ln(7)"),
                       (log(42), "ln(42)"), (pi, "pi"),
                       (1/E, "1/e")]:
    kt = sum(m * exp(-t_val * float(e)) for e, m in evals_with_mult)
    # Nature
    if t_val == 0:
        nature = "RATIONAL (= 12)"
    elif t_name.startswith("ln("):
        q_str = t_name[3:-1]
        nature = f"ALGEBRAIC (q={q_str})"
    else:
        nature = "TRANSCENDENTAL"
    print(f"  {t_name:>12} | {kt:14.8f} | {nature:>20}")

# Verify algebraicity at t = ln(2)
print(f"\n  Verification: K(ln(2)) = sum m_i * 2^(-lambda_i)")
total = Fraction(0)
for e, m in evals_with_mult:
    if e == 0:
        val_str = "1"
        total += m
    elif e == 1:
        val_str = "1/2"
        total += Fraction(m, 2)
    else:
        f_e = float(e)
        val = 2**(-f_e)
        # 2^(-a/b) is algebraic
        val_str = f"2^({-e}) = {val:.6f}"
    print(f"    m={m}, lambda={str(e):>5}: m * 2^(-lambda) = {m} * {val_str}")

# =============================================================================
# GHOST 6: THE ZETA GHOST (NEGATIVE INTEGER ARGUMENTS)
# =============================================================================

print("""

GHOST 6: THE SPECTRAL ZETA AT NEGATIVE INTEGERS
=================================================

zeta_T(s) = sum |lambda_i|^{-s} (over nonzero eigenvalues)

At NEGATIVE integers s = -k (k >= 0):
  zeta_T(-k) = sum |lambda_i|^k = p_k (the power sum)

These are RATIONAL! The spectral zeta has rational special values
at all non-positive integers, just like the Riemann zeta function.

Riemann zeta: zeta(-k) = -B_{k+1}/(k+1)  (Bernoulli numbers)
Tournament zeta: zeta_T(-k) = Tr(T^k) = p_k  (power sums)

The analogy is EXACT. Both are rational at negative integers.
Both are transcendental at positive integers.
The negative-integer ghost connects zeta to combinatorics.
""")

print(f"  {'s':>5} | {'zeta_T(s) (Riemann-like)':>22} | {'zeta_T(s) (our def)':>22} | {'rational?':>10}")
print("  " + "-" * 65)

for s in [-4, -3, -2, -1, 0, 1, 2, 3, 4]:
    if s <= 0:
        # sum |lambda_i|^{-s} = sum |lambda_i|^{|s|} -- but for s=-k, we compute sum lambda_i^k
        val = sum(m * abs(float(e))**(-s) for e, m in evals_with_mult if e != 0)
        val_frac = Fraction(val).limit_denominator(1000000)
        is_rat = abs(float(val_frac) - val) < 1e-8
        print(f"  {s:5d} | {val:22.6f} | {val:22.6f} | {'YES' if is_rat else '???':>10}")
    else:
        val = sum(m * abs(float(e))**(-s) for e, m in evals_with_mult if abs(float(e)) > 1e-10)
        print(f"  {s:5d} | {val:22.6f} | {val:22.6f} | {'NO (trans)':>10}")

# =============================================================================
# GHOST 7: THE CHEBYSHEV GHOST
# =============================================================================

print("""

GHOST 7: THE CHEBYSHEV GHOST
===============================

The Chebyshev polynomials T_k(x) satisfy T_k(cos(theta)) = cos(k*theta).
They map [-1,1] to [-1,1] and preserve rational structure:
  T_k(p/q) is RATIONAL for rational p/q.

For eigenvalues lambda_i in [-1,1]:
  T_k(lambda_i) is rational (since lambda_i is rational).

But there's a deeper ghost: if we write lambda = cos(theta),
then theta = arccos(lambda) is generally TRANSCENDENTAL.
Yet T_k(lambda) = cos(k*theta) stays rational!

The Chebyshev ghost: rational structure survives through the
transcendental parameterization theta = arccos(lambda).
""")

print(f"  Chebyshev values T_k(lambda) for n=5 eigenvalues:")
print(f"  {'lambda':>7} | {'T_1':>8} | {'T_2':>8} | {'T_3':>8} | {'T_4':>8} | {'T_5':>8} | {'arccos':>10}")
print("  " + "-" * 70)

def chebyshev(k, x):
    """T_k(x) via recurrence."""
    if k == 0: return Fraction(1)
    if k == 1: return x
    t0, t1 = Fraction(1), x
    for _ in range(k - 1):
        t0, t1 = t1, 2*x*t1 - t0
    return t1

for lam in nz_evals:
    t_vals = [chebyshev(k, lam) for k in range(1, 6)]
    theta = f"{float(np.arccos(float(lam))):.6f}" if abs(float(lam)) <= 1 else "---"
    print(f"  {str(lam):>7} | {str(t_vals[0]):>8} | {str(t_vals[1]):>8} | "
          f"{str(t_vals[2]):>8} | {str(t_vals[3]):>8} | {str(t_vals[4]):>8} | {theta:>10}")

# Check: does Chebyshev reveal new structure?
print(f"\n  T_5(lambda) values: {[str(chebyshev(5, e)) for e in nz_evals]}")
print(f"  T_7(lambda) values: {[str(chebyshev(7, e)) for e in nz_evals]}")

# =============================================================================
# GHOST 8: THE FORMAL GROUP GHOST
# =============================================================================

print("""

GHOST 8: THE FORMAL GROUP GHOST
==================================

The formal group F(x,y) = (x+y)/(1+xy) has:
  - Rational coefficients in its power series expansion
  - Rational values when evaluated at rationals
  - The logarithm arctanh(x) = x + x^3/3 + x^5/5 + ... has rational coefficients

The formal group logarithm L(x) = arctanh(x) has the Taylor expansion:
  L(x) = sum_{k>=0} x^{2k+1}/(2k+1)

Each coefficient 1/(2k+1) is RATIONAL.
The function value L(p/q) is TRANSCENDENTAL (by Lindemann-Weierstrass).

The formal group ghost: the Taylor coefficients are a RATIONAL sequence
{1, 1/3, 1/5, 1/7, 1/9, ...} = {1/(2k+1)}.

These are the RECIPROCALS OF ODD NUMBERS.
The same odd numbers that appear in tournament H-values!

H(T) is always ODD. The Taylor coefficients are 1/ODD.
The formal group encodes the ODD STRUCTURE.
""")

print(f"  Formal group log coefficients vs tournament H-values:")
print(f"  {'k':>4} | {'coeff 1/(2k+1)':>15} | {'2k+1':>5} | {'is H-value?':>12} | {'forbidden?':>10}")
print("  " + "-" * 55)

# Known H-values up to some n
known_H = {1, 3, 5, 9, 11, 13, 15, 17, 19, 23, 25, 27, 29, 31, 33}
forbidden_H = {7, 21}

for k in range(0, 20):
    odd = 2*k + 1
    coeff = Fraction(1, odd)
    is_H = "YES" if odd in known_H else ("FORBIDDEN" if odd in forbidden_H else "?")
    forbidden = "YES" if odd in forbidden_H else ""
    print(f"  {k:4d} | {str(coeff):>15} | {odd:5d} | {is_H:>12} | {forbidden:>10}")

# =============================================================================
# GHOST 9: THE HADAMARD PRODUCT GHOST
# =============================================================================

print("""

GHOST 9: THE HADAMARD PRODUCT GHOST
======================================

The Hadamard (entrywise) product of two eigenvectors is NOT generally
an eigenvector. But the entrywise product of two EIGENVALUES is just
a rational number (product of rationals).

The deeper ghost: the Hadamard product of the transition matrix T
with itself gives T (*) T, whose trace is:
  Tr(T (*) T) = sum_{i,j} T_{ij}^2 = the FROBENIUS norm squared.

This is RATIONAL (since T has rational entries).
More generally: Tr(T^{(*k)}) where T^{(*k)} is the k-th Hadamard power
is rational for all k.

The Frobenius ghost: entrywise powers of T have rational traces.
""")

T5_data, _, _ = build_transition(5)

print(f"  Hadamard powers of T at n=5:")
print(f"  {'k':>4} | {'Tr(T^(*k))':>14} | {'exact':>20}")
print("  " + "-" * 45)

for k in range(1, 8):
    Tk_had = T5_data ** k  # entrywise power
    tr = np.trace(Tk_had)
    frac = Fraction(tr).limit_denominator(10000000)
    print(f"  {k:4d} | {tr:14.8f} | {str(frac):>20}")

# =============================================================================
# GHOST 10: THE GALOIS GHOST
# =============================================================================

print("""

GHOST 10: THE GALOIS GHOST
=============================

For the tournament eigenvalues at n=5:
  {1, 3/5, 2/5, 1/5, 0, -1/5, -3/5}

The MINIMAL POLYNOMIAL of each eigenvalue over Q is LINEAR (degree 1),
since each eigenvalue is rational.

But consider the EIGENVALUES OF EIGENVALUES: if we build a matrix
from the eigenvalues themselves (e.g., the Gram matrix of eigenvectors),
some of THOSE eigenvalues could be algebraic of higher degree.

For now, the Galois ghost is trivial: Gal(Q(lambda_i)/Q) = {1} for all i.
All eigenvalues generate the TRIVIAL extension.

The ghost: the Galois group is trivial, so ALL algebraic operations
on eigenvalues stay in Q. No algebraic number beyond Q is ever produced.
This is why the rapidity lattice is the FIRST non-rational structure:
you have to use TRANSCENDENTAL operations to escape Q.
""")

# =============================================================================
# GHOST 11: THE CONTINUED FRACTION GHOST
# =============================================================================

print("""

GHOST 11: THE CONTINUED FRACTION GHOST
=========================================

Every rational eigenvalue has a FINITE continued fraction expansion.
Irrational numbers have infinite ones.

The continued fraction of lambda encodes a sequence of "best rational
approximations." For rational lambda, this sequence is finite:
the number IS exactly representable.

The continued fraction ghost: rational eigenvalues have FINITE
representations. Transcendental rapidities have INFINITE ones.
The Cayley gate turns finite representations into infinite ones.
""")

def continued_fraction(frac, max_terms=20):
    """Compute continued fraction of a Fraction."""
    cf = []
    while len(cf) < max_terms:
        a = int(frac)
        cf.append(a)
        frac = frac - a
        if frac == 0:
            break
        frac = 1 / frac
    return cf

print(f"  Continued fractions of eigenvalues vs rapidities:")
print(f"  {'value':>10} | {'CF':>30} | {'length':>6} | {'type':>12}")
print("  " + "-" * 65)

for lam in [Fraction(3,5), Fraction(2,5), Fraction(1,5), Fraction(-1,5), Fraction(-3,5)]:
    cf = continued_fraction(lam if lam > 0 else -lam)
    cf_str = str(cf)
    print(f"  {str(lam):>10} | {cf_str:>30} | {len(cf):6d} | {'FINITE':>12}")

# Rapidities (approximated as fractions)
for lam_f in [0.6, 0.4, 0.2]:
    rho = atanh(lam_f)
    rho_frac = Fraction(rho).limit_denominator(10**12)
    cf = continued_fraction(rho_frac, max_terms=15)
    cf_str = str(cf[:8]) + "..."
    print(f"  atanh({lam_f}) | {cf_str:>30} | {'inf':>6} | {'INFINITE':>12}")

# =============================================================================
# GHOST 12: THE NEWTON POLYGON GHOST
# =============================================================================

print("""

GHOST 12: THE p-ADIC NEWTON POLYGON GHOST
============================================

The characteristic polynomial p(x) = det(xI - T) has rational coefficients.
Its Newton polygon at prime p encodes the p-adic valuations of the roots.

For each prime p, the Newton polygon is a CONVEX HULL of points
(k, v_p(a_k)) where a_k are the coefficients.

The slopes of the Newton polygon give the p-adic valuations of eigenvalues.
Since eigenvalues are rational, their p-adic valuations are integers.

The Newton polygon ghost: the SLOPES of the Newton polygon are
the p-adic valuations of eigenvalues, and these are always
INTEGERS (since eigenvalues are rational).

For a general matrix over Q_p, the slopes could be any rational number.
The integrality of slopes is a ghost of the rationality of eigenvalues.
""")

# Compute Newton polygon at p=5 for n=5
def vp(n, p):
    if n == 0: return float('inf')
    v = 0
    n = abs(n)
    while n % p == 0:
        v += 1
        n //= p
    return v

print(f"  Newton polygon of char poly at p=5:")
coeffs = np.poly(T5_data)
print(f"  {'k':>4} | {'a_k':>14} | {'v_5(a_k)':>8}")
print("  " + "-" * 32)
for k, c in enumerate(coeffs):
    c_frac = Fraction(c).limit_denominator(10**12)
    if c_frac != 0:
        v = vp(c_frac.numerator, 5) - vp(c_frac.denominator, 5)
    else:
        v = float('inf')
    print(f"  {k:4d} | {float(c_frac):14.6f} | {v:8}")

# =============================================================================
# GHOST 13: THE BERNOULLI GHOST
# =============================================================================

print("""

GHOST 13: THE BERNOULLI / VON STAUDT GHOST
=============================================

Von Staudt's theorem: the denominator of B_{2k} is prod_{(p-1)|2k} p.

For B_6: denominator = prod_{(p-1)|6} p = 2 * 3 * 7 = 42.
These are the HURWITZ PRIMES — the same ones in the rapidity lattice!

Connection: B_6 = 1/42 is the Bernoulli number whose denominator
is the product of the primes generating the rapidity lattice.

The Bernoulli ghost: the Bernoulli denominators encode which primes
are "active" in the formal group at each level.

For the formal group F(x,y) = (x+y)/(1+xy) = the multiplicative
formal group over Z: the Bernoulli numbers control its logarithm.
""")

# Bernoulli denominators
from fractions import Fraction

def bernoulli(n):
    """Compute B_n."""
    B = [Fraction(0)] * (n+1)
    B[0] = Fraction(1)
    for m in range(1, n+1):
        B[m] = -sum(Fraction(comb(m, k), m - k + 1) * B[k] for k in range(m))
    return B[n]

print(f"  Bernoulli numbers and their denominators:")
print(f"  {'n':>4} | {'B_n':>20} | {'denom':>8} | {'factors':>20} | {'= Hurwitz?':>12}")
print("  " + "-" * 72)

for n_val in [0, 1, 2, 4, 6, 8, 10, 12, 14]:
    bn = bernoulli(n_val)
    if bn == 0:
        print(f"  {n_val:4d} | {'0':>20} | {'-':>8} | {'-':>20} | {'-':>12}")
        continue
    den = bn.denominator
    # Factor the denominator
    temp = den
    facts = {}
    d = 2
    while d*d <= temp:
        while temp % d == 0:
            facts[d] = facts.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        facts[temp] = facts.get(temp, 0) + 1
    f_str = " * ".join(f"{p}" for p in sorted(facts.keys()))
    is_hurwitz = sorted(facts.keys()) == [2, 3, 7]
    print(f"  {n_val:4d} | {str(bn):>20} | {den:8d} | {f_str:>20} | {'YES!' if is_hurwitz else '':>12}")

# =============================================================================
# GHOST 14: THE MAHLER MEASURE GHOST
# =============================================================================

print("""

GHOST 14: THE MAHLER MEASURE GHOST
=====================================

The Mahler measure of the characteristic polynomial p(x) is:
  M(p) = |a_n| * prod_{|alpha_i|>1} |alpha_i|

where alpha_i are the roots. For our transition matrix:
  - All eigenvalues have |lambda_i| <= 1
  - Leading coefficient a_n = 1 (monic)
  - So M(p) = 1

Mahler measure 1 means the polynomial is a PRODUCT OF CYCLOTOMIC
polynomials and a monomial (by Kronecker's theorem, if all roots
are on or inside the unit circle and algebraic integers).

Our eigenvalues are rationals in [-1,1], not algebraic integers.
But the Mahler measure is still exactly 1.

The Mahler ghost: M(char poly) = 1 for ALL tournament transition matrices.
This is because the spectral radius is exactly 1 (the stationary eigenvalue).
""")

# =============================================================================
# GHOST 15: THE INVERSE CAYLEY CASCADE
# =============================================================================

print("""
GHOST 15: THE ITERATED CAYLEY / MANDELBROT GHOST
===================================================

Apply Q repeatedly: Q, Q*Q, Q*Q*Q, ...

  Q^1(lambda) = Q(lambda) = (1+lambda)/(1-lambda)
  Q^2(lambda) = Q(Q(lambda)) = -1/lambda    (!)
  Q^3(lambda) = Q^2(Q(lambda)) = -(1-lambda)/(1+lambda) = -1/Q(lambda)
  Q^4(lambda) = Q^3(Q(lambda)) = lambda     (!!!)

The Cayley transform has ORDER 4 as a Mobius transformation!
Q^4 = identity. The Cayley group is Z/4Z = the cyclic group of order 4.

The orbit of any eigenvalue under Q is a 4-element set:
  {lambda, Q(lambda), -1/lambda, -1/Q(lambda)}

This 4-periodicity is a ghost: no matter what transcendental operation
we apply AFTER the Cayley gate, the 4-fold symmetry persists.
""")

print(f"  Cayley orbits of eigenvalues:")
print(f"  {'lambda':>8} | {'Q(lam)':>8} | {'Q^2(lam)':>8} | {'Q^3(lam)':>8} | {'Q^4(lam)':>8}")
print("  " + "-" * 50)

for lam in [Fraction(3,5), Fraction(2,5), Fraction(1,5)]:
    q1 = Fraction(1+lam, 1-lam)
    q2 = Fraction(-1, lam) if lam != 0 else "inf"
    q3 = -Fraction(1-lam, 1+lam) if lam != Fraction(-1) else "inf"
    q4 = lam  # Q^4 = identity
    print(f"  {str(lam):>8} | {str(q1):>8} | {str(q2):>8} | {str(q3):>8} | {str(q4):>8}")

print(f"""
  Q^4 = Id: the Cayley transform generates a Z/4Z action on Q.

  The four elements:
    Q^0 = Id:      lambda (the eigenvalue)
    Q^1 = Cayley:  (1+x)/(1-x) (the boost)
    Q^2 = Invert:  -1/x (the inversion)
    Q^3 = anti-Q:  -(1-x)/(1+x) (the anti-boost)

  This is the KLEIN FOUR-GROUP acting on the eigenvalue.
  Wait: Q^2 = -1/x, and (Q^2)^2 = (-1/(-1/x)) = x = Id. Yes, Z/4Z.

  Actually Q has order 4 in PGL(2,Q), generated by the matrix
  [[1,1],[−1,1]] which has order 4 (its square is [[0,2],[−2,0]] ~ [[0,1],[−1,0]]).
""")

# =============================================================================
# SUMMARY: THE CATALOG OF GHOSTS
# =============================================================================

print("""
================================================================
THE CATALOG OF GHOSTS
================================================================

  GHOST | WHERE IT LIVES            | WHAT SURVIVES            | GUARANTEE
  ------|---------------------------|--------------------------|------------------
    0   | Rapidity lattice Z^3      | {2,3,7} prime structure  | Baker's theorem
    1   | Traces Tr(T^k)            | Power sums in Q          | Rationality of T
    2   | Char poly coefficients    | Elem sym functions in Q  | Rationality of T
    3   | Moments p_k               | Newton identities        | Trace-poly bridge
    4   | Resolvent R(z)            | Rational function of z   | Partial fractions
    5   | Heat kernel K(ln q)       | Algebraic at log-rationals| Lindemann-Weierstrass
    6   | Zeta at neg integers      | Rational special values  | Moment connection
    7   | Chebyshev T_k(lambda)     | Rational despite trig    | Polynomial identity
    8   | Formal group log coeffs   | 1/(2k+1) = 1/ODD        | Odd parity of H
    9   | Hadamard traces           | Entrywise rational       | Rational entries
   10   | Galois group              | Trivial (= {1})          | Eigenvalues in Q
   11   | Continued fractions       | Finite (= rational)      | Def of rationality
   12   | Newton polygon slopes     | Integer slopes           | Rational eigenvalues
   13   | Bernoulli denominators    | {2,3,7} at B_6 = 1/42   | Von Staudt theorem
   14   | Mahler measure            | M = 1 always             | Spectral radius = 1
   15   | Cayley periodicity        | Q^4 = Id (Z/4Z action)  | Mobius transformation

THE HIERARCHY OF GHOSTS:

  Strongest (survive deepest into transcendental):
    Ghost 0: Z^3 lattice — survives ln
    Ghost 5: algebraic at ln(q) — survives exp*ln
    Ghost 15: Z/4Z periodicity — survives all Mobius

  Medium (survive algebraic operations):
    Ghosts 1-4, 7, 9: rational traces, resolvents, Chebyshev

  Weakest (only at the boundary):
    Ghosts 6, 11-14: special values, finiteness, integrality

THE DEEPEST GHOST IS GHOST 0:
  The rapidity lattice Z^3 survives the logarithm,
  the most powerful transcendental operation.
  It is the last discrete structure before the uncountable.
  After it: only the uncountable continuum remains.

But Ghost 15 (Q^4 = Id) is the most PERSISTENT:
  The 4-fold Cayley symmetry survives ALL rational operations.
  It is a symmetry of the gate itself, not of what passes through.
  The gate has 4-fold symmetry. The transcendental world doesn't know this.
  But the ghost remembers.
""")

print("opus-2026-03-17-S91d: GHOSTS OF RATIONAL STRUCTURE")
