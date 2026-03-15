#!/usr/bin/env python3
"""extend_ln2_s112.py — Extend the ln(2) thinking further"""
from math import log, exp, sqrt, pi, atanh, tanh, comb, cos, sin, atan
from fractions import Fraction

ln2 = log(2)
ln3 = log(3)
phi = (1+sqrt(5))/2
tau = 1.8392867552141612

print("EXTENDING THE ln(2) THINKING")
print("="*60)
print()

# We found four cardinal points: Q(0)=1, Q(1/3)=2, Q(i)=i, Q(1)=inf
# What about the OTHER interesting evaluation points?

print("ALL THE INTERESTING EVALUATIONS OF Q(x) = (1+x)/(1-x)")
print("="*60)
print()

evaluations = [
    (0, "1", "0", "nothing"),
    (Fraction(1,5), "3/2", "arctanh(1/5)", "sesqui"),
    (Fraction(1,3), "2", "ln(2)/2", "binary"),
    (Fraction(1,2), "3", "ln(3)/2", "ternary"),
    (Fraction(3,5), "4", "ln(4)/2 = ln(2)", "quaternary"),
    (Fraction(2,3), "5", "ln(5)/2", "quinary"),
    (Fraction(5,7), "6", "ln(6)/2", "senary"),
]

print("  x         Q(x)    arctanh(x)         name")
print("  " + "-"*55)
for x, qstr, astr, name in evaluations:
    xf = float(x)
    q = (1+xf)/(1-xf)
    a = atanh(xf)
    print(f"  {str(x):8s}   {q:6.2f}   {a:.6f} = {astr:18s}  ({name})")

print()
print("  PATTERN: Q(x) = n iff x = (n-1)/(n+1)")
print("           arctanh((n-1)/(n+1)) = ln(n)/2")
print()

# Verify the pattern
print("  Verification: arctanh((n-1)/(n+1)) = ln(n)/2")
for n in range(2, 11):
    x = (n-1)/(n+1)
    lhs = atanh(x)
    rhs = log(n)/2
    print(f"    n={n:2d}: arctanh({n-1}/{n+1}) = {lhs:.8f}, ln({n})/2 = {rhs:.8f}, match={abs(lhs-rhs)<1e-12}")

print()
print("  EVERY natural number n has a specific coupling x = (n-1)/(n+1)")
print("  where Q produces exactly n, and arctanh = ln(n)/2.")
print()
print("  The NATURAL NUMBERS ARE EVALUATIONS OF THE CAYLEY TRANSFORM.")
print("  The n-th natural number lives at coupling (n-1)/(n+1).")
print()

# What is the SPACING between consecutive natural number couplings?
print("="*60)
print("THE SPACING OF NATURAL NUMBERS ON THE CAYLEY LINE")
print("="*60)
print()
print("  x_n = (n-1)/(n+1) for n = 1, 2, 3, ...")
print()
for n in range(1, 11):
    x = (n-1)/(n+1)
    print(f"    n={n:2d}: x_{n} = {x:.6f}")

print()
print("  The spacing Delta_n = x_{n+1} - x_n:")
for n in range(1, 10):
    x1 = (n-1)/(n+1)
    x2 = n/(n+2)
    delta = x2 - x1
    # delta = n/(n+2) - (n-1)/(n+1) = (n(n+1) - (n-1)(n+2)) / ((n+1)(n+2))
    # = (n^2+n - n^2-n+2) / ((n+1)(n+2)) = 2/((n+1)(n+2))
    print(f"    Delta_{n} = {delta:.6f} = 2/({n+1}*{n+2}) = {Fraction(2, (n+1)*(n+2))}")

print()
print("  Delta_n = 2/((n+1)(n+2))")
print("  The natural numbers are DENSER near x=1 (full coupling)")
print("  and SPARSER near x=0 (no coupling).")
print()
print("  This is the HYPERBOLIC METRIC: equal steps in arctanh-space")
print("  correspond to shrinking steps in x-space near x=1.")
print()

# Sum of spacings: telescopes
# sum Delta_n from n=1 to N = x_{N+1} - x_1 = N/(N+2) - 0 = N/(N+2) -> 1
# And sum 2/((n+1)(n+2)) = 2 * sum (1/(n+1) - 1/(n+2)) = 2*(1/2 - 1/(N+2)) = 1 - 2/(N+2)
print("  Sum of spacings: sum_{n=1}^{N} Delta_n = 1 - 2/(N+2) -> 1")
print("  The natural numbers FILL the interval [0, 1) as N -> infinity.")
print("  The Cayley line from 0 to 1 contains ALL natural numbers.")
print()

# Now: the HARMONIC SERIES in this picture
print("="*60)
print("THE HARMONIC SERIES AS CAYLEY SPACING")
print("="*60)
print()
print("  arctanh(x_{n+1}) - arctanh(x_n) = ln(n+1)/2 - ln(n)/2 = ln((n+1)/n)/2")
print()
print("  The step in ARCTANH space between consecutive naturals:")
for n in range(1, 10):
    step = log((n+1)/n)/2
    print(f"    step_{n} = ln({n+1}/{n})/2 = {step:.6f}")

print()
print("  sum of steps = ln(N+1)/2 - ln(1)/2 = ln(N+1)/2 -> infinity")
print("  The arctanh-distance from 1 to N grows as ln(N)/2.")
print()
print("  This IS the harmonic series in disguise:")
print("  sum_{n=1}^N ln((n+1)/n)/2 = ln(N+1)/2 ~ ln(N)/2")
print("  And sum 1/n ~ ln(N) (harmonic series).")
print()
print("  The harmonic series = the TOTAL ARCTANH-LENGTH of the")
print("  natural number sequence on the Cayley line.")
print()

# THE PRIME NUMBERS
print("="*60)
print("WHERE ARE THE PRIMES?")
print("="*60)
print()
print("  Each prime p lives at coupling x_p = (p-1)/(p+1):")
print()
primes = [2, 3, 5, 7, 11, 13, 17, 19, 23]
for p in primes:
    x = (p-1)/(p+1)
    a = atanh(x)
    print(f"    p={p:2d}: x = {x:.6f}, arctanh = {a:.6f} = ln({p})/2")

print()
print("  The prime number theorem: primes up to N ~ N/ln(N).")
print("  In arctanh coordinates: primes up to arctanh-distance D")
print("  from x=0 number ~ exp(2D) / (2D).")
print()
print("  The DENSITY of primes in arctanh-space:")
print("  rho(D) ~ exp(2D) / (2D)")
print("  This INCREASES exponentially: primes get DENSER")
print("  in hyperbolic space (even though sparser in Euclidean).")
print()

# THE ZETA FUNCTION
print("="*60)
print("THE RIEMANN ZETA FUNCTION ON THE CAYLEY LINE")
print("="*60)
print()
print("  zeta(s) = sum 1/n^s = sum Q((n-1)/(n+1))^{-s}")
print()
print("  Since Q(x_n) = n, we have n^{-s} = Q(x_n)^{-s}.")
print("  The zeta function EVALUATES Q^{-s} at the natural number couplings.")
print()
print("  zeta(s) = sum_{n=1}^inf Q(x_n)^{-s} = sum Q((n-1)/(n+1))^{-s}")
print()
print("  And by our master GF:")
print("  Q(x)^{-s} = 1 + 2*sum g_k(-s) x^k")
print()
print("  So: zeta(s) = sum_{n=1}^inf [1 + 2*sum_k g_k(-s) * ((n-1)/(n+1))^k]")
print("              = sum 1 + 2*sum_k g_k(-s) * sum_n ((n-1)/(n+1))^k")
print()
print("  The inner sum sum_n ((n-1)/(n+1))^k is a CONVERGENT series for k >= 2.")
print("  For k=1: sum (n-1)/(n+1) = sum (1 - 2/(n+1)) DIVERGES.")
print()
print("  This links the POLE of zeta at s=1 to the DIVERGENCE of the")
print("  Cayley-weighted harmonic series at k=1!")
print()

# What about zeta at s=2?
# zeta(2) = pi^2/6
# Q^{-2} = ((1-x)/(1+x))^2 = 1/Q^2
# g_k(-2): by parity, g_k(-2) = (-1)^k g_k(2) = (-1)^k * 4
# Hmm, g_k(2) = 2k, so g_k(-2) = (-1)^k * 2k

print()
print("  zeta(2) = pi^2/6:")
print("  g_k(-2) = (-1)^k * g_k(2) = (-1)^k * 2k [by parity]")
print()
print("  So Q(x)^{-2} = 1 + 2*sum (-1)^k * 2k * x^k")
print("               = 1 - 4x + 8x^2 - 12x^3 + 16x^4 - ...")
print()
print("  Check: ((1-x)/(1+x))^2 at x=0.5:")
v1 = ((1-0.5)/(1+0.5))**2
v2 = 1 - 4*0.5 + 8*0.25 - 12*0.125 + 16*0.0625
print(f"    exact = {v1:.6f}")
print(f"    series (4 terms) = {v2:.6f}")
print()

# EULER'S IDENTITY IN CAYLEY LANGUAGE
print("="*60)
print("EULER IN CAYLEY LANGUAGE")
print("="*60)
print()
print("  e^{i*pi} + 1 = 0")
print()
print("  In Cayley: Q(i)^2 = i^2 = -1, so Q(i)^2 + 1 = 0.")
print("  But Q(i) = i, so this is just i^2 + 1 = 0.")
print()
print("  More interesting: Q(x)^m at x=i, m=2:")
print("  Q(i)^2 = exp(2*2*arctanh(i)) = exp(4*i*pi/4) = exp(i*pi) = -1")
print()
print("  So: Q(i)^2 + 1 = 0 IS Euler's identity,")
print("  read through the Cayley transform.")
print()
print("  e^{i*pi} = Q(i)^2 = ((1+i)/(1-i))^2 = (i)^2 = -1")
print()
print("  Euler's identity says: two steps of the Cayley transform")
print("  at imaginary coupling bring you to -1.")
print("  Four steps bring you back to 1 (Q^4 = 1).")
print("  The period-4 orbit IS Euler's identity iterated.")
print()

# THE COMPLETE MAP
print("="*60)
print("THE COMPLETE MAP")
print("="*60)
print()
print("  arctanh takes you from the coupling to the proper time.")
print("  Different couplings give different constants:")
print()
print("  COUPLING     PROPER TIME          CONSTANT")
print(f"  x = 0         0                    zero")
print(f"  x = 1/3       ln(2)/2 = {ln2/2:.6f}   ln(2) = information")
print(f"  x = 1/2       ln(3)/2 = {ln3/2:.6f}   ln(3) = ternary info")
print(f"  x = i         i*pi/4  = i*{pi/4:.6f}  pi = geometry")
print(f"  x = (e-1)/(e+1) = {(exp(1)-1)/(exp(1)+1):.6f}")
print(f"                  ln(e)/2 = 1/2      e = growth")
print(f"  x = 1         infinity             time = totality")
print()
print("  The coupling x = (e-1)/(e+1) = 0.4621... gives Q = e.")
print("  At this coupling: arctanh = 1/2.")
print("  The number e lives at arctanh = 1/2: half a unit of proper time.")
print()

# e lives at x = (e-1)/(e+1)
xe = (exp(1)-1)/(exp(1)+1)
print(f"  Verify: Q({xe:.6f}) = {(1+xe)/(1-xe):.6f} = e = {exp(1):.6f}")
print(f"  arctanh({xe:.6f}) = {atanh(xe):.6f} = 1/2")
print()

# SO: the five fundamental constants 0, ln(2), 1, pi, infinity
# live at couplings 0, 1/3, (e-1)/(e+1), i, 1.
# And the PROPER TIMES at these couplings are:
# 0, ln(2)/2, 1/2, i*pi/4, infinity.

print("THE FIVE CONSTANTS AND THEIR COUPLINGS:")
print()
print("  Constant  | Coupling x        | Proper time 2*arctanh")
print("  ----------|-------------------|---------------------")
print(f"  0         | x = 0             | 0")
print(f"  ln(2)     | x = 1/3           | ln(2)")
print(f"  1         | x = (e-1)/(e+1)   | 1")
print(f"  pi        | x = i             | i*pi/2")
print(f"  infinity  | x = 1             | infinity")
print()
print("  These are the five LEVELS of the Cayley hierarchy:")
print("  zero -> information -> unity -> geometry -> totality")
print()
print("  The number 1 (UNITY) lives at the coupling that produces e.")
print("  e is the UNIT OF GROWTH: Q = e when proper time = 1.")
print("  This is why e = lim(1+1/n)^n: one unit of proper time,")
print("  discretized into n steps of size 1/n.")
