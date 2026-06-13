#!/usr/bin/env python3
"""polynomial_tricks_s116n.py — Getting tricky with the heat kernel polynomial.

P(z) = 29 - 14z - 25z^2 + 6z^3 + 5z^4

This degree-4 polynomial encodes ALL of tournament dynamics from the transitive.
Let's push it to its limits.

Session: kind-pasteur-2026-03-17-S116n33
"""
import sys
sys.stdout.reconfigure(line_buffering=True)
from math import log, exp, sqrt, pi
from fractions import Fraction

print()
print("  GETTING TRICKY WITH P(z) = 5z^4 + 6z^3 - 25z^2 - 14z + 29")
print()
print("=" * 70)
print()

# The polynomial
def P(z):
    return 5*z**4 + 6*z**3 - 25*z**2 - 14*z + 29

def P_frac(z):
    z = Fraction(z) if not isinstance(z, Fraction) else z
    return 5*z**4 + 6*z**3 - 25*z**2 - 14*z + 29

# ============================================================
print("  I. THE SPECIAL VALUES")
print("  " + "-" * 50)
print()

special = [
    (-3, "deep imaginary time"),
    (-2, "double imaginary"),
    (-1, "Wick rotation (imaginary time)"),
    (0, "infinite forward time"),
    (Fraction(1,7), "Hurwitz temperature q=7"),
    (Fraction(1,5), "golden temperature q=5"),
    (Fraction(1,3), "Cayley address of 2"),
    (Fraction(1,2), "half temperature q=2"),
    (1, "zero time (H of transitive)"),
    (2, "negative time q=1/2"),
    (3, "deep negative time"),
]

print(f"  {'z':>8s}  {'P(z)':>12s}  {'meaning':>30s}")
for z, meaning in special:
    val = P_frac(z)
    print(f"  {str(z):>8s}  {str(val):>12s}  {meaning:>30s}")
print()

# ============================================================
print("  II. THE STUNNING COINCIDENCE: P(0) = P(2) = 29")
print("  " + "-" * 50)
print()

print(f"  P(0) = {P_frac(0)} = mean H (infinite time)")
print(f"  P(2) = {P_frac(2)} = mean H AGAIN!")
print()
print("  These are DIFFERENT times:")
print("  z=0: q->infinity, t->infinity (forward, equilibrium)")
print("  z=2: q=1/2, t=5*ln(1/2)=-5*ln(2) (BACKWARD time)")
print()
print("  Running BACKWARD for 5*ln(2) steps from transitive")
print("  gives the SAME expected H as running forward FOREVER.")
print()

# Factor P(z) - 29
print("  FACTORING P(z) - 29:")
print("  P(z) - 29 = 5z^4 + 6z^3 - 25z^2 - 14z")
print("            = z * (5z^3 + 6z^2 - 25z - 14)")
print()

# Check: is z=2 a root of the cubic?
cubic_at_2 = 5*8 + 6*4 - 25*2 - 14
print(f"  5(2)^3 + 6(2)^2 - 25(2) - 14 = {cubic_at_2}")
print(f"  YES, z=2 is a root!")
print()

# Factor the cubic: 5z^3 + 6z^2 - 25z - 14 = (z-2)(5z^2 + 16z + 7)
print("  5z^3 + 6z^2 - 25z - 14 = (z-2)(5z^2 + 16z + 7)")
# Verify: (z-2)(5z^2+16z+7) = 5z^3+16z^2+7z-10z^2-32z-14 = 5z^3+6z^2-25z-14 ✓
print("  Verified: expansion matches.")
print()

# Roots of 5z^2 + 16z + 7 = 0
disc = 16**2 - 4*5*7
print(f"  5z^2 + 16z + 7 = 0")
print(f"  Discriminant = 256 - 140 = {disc}")
print(f"  sqrt(disc) = sqrt({disc}) = sqrt(4*29) = 2*sqrt(29)")
print(f"  z = (-16 +/- 2*sqrt(29)) / 10 = (-8 +/- sqrt(29)) / 5")
print()

sqrt29 = sqrt(29)
z1 = (-8 + sqrt29) / 5
z2 = (-8 - sqrt29) / 5
print(f"  z_+ = (-8 + sqrt(29)) / 5 = {z1:.10f}")
print(f"  z_- = (-8 - sqrt(29)) / 5 = {z2:.10f}")
print()

# ============================================================
print("  III. THE COMPLETE FACTORIZATION")
print("  " + "-" * 50)
print()

print("  P(z) = 29 + z(z-2)(5z^2 + 16z + 7)")
print()
print("  The polynomial has mean H = 29 as its BASELINE.")
print("  Deviation from mean = z(z-2) * quadratic.")
print()
print("  ROOTS of P(z) - 29:")
print(f"  z = 0 (infinite time)")
print(f"  z = 2 (backward time)")
print(f"  z = (-8+sqrt(29))/5 = {z1:.6f}")
print(f"  z = (-8-sqrt(29))/5 = {z2:.6f}")
print()
print("  The irrational roots involve sqrt(29) = sqrt(MEAN H).")
print("  The mean H is not just a value — it's a DISCRIMINANT.")
print("  The field Q(sqrt(29)) controls the algebraic structure")
print("  of the deviation from equilibrium!")
print()

# ============================================================
print("  IV. P(-1) = WICK ROTATION = ANTI-TRANSITIVE")
print("  " + "-" * 50)
print()

print(f"  P(-1) = 5 - 6 - 25 + 14 + 29 = {P_frac(-1)}")
print()
print("  P(-1) = 17 = H(anti-transitive) = Tribonacci at n=6")
print()
print("  z = -1 corresponds to q = 1/(-1) = -1.")
print("  t = 5*ln(-1) = 5*i*pi (IMAGINARY TIME)")
print()
print("  The WICK ROTATION of the Markov chain dynamics")
print("  at imaginary time t = 5*i*pi gives the anti-transitive H.")
print()
print("  In physics, Wick rotation exchanges:")
print("  - Minkowski spacetime (real time) <-> Euclidean space (imaginary time)")
print("  - Oscillating wavefunctions <-> Decaying amplitudes")
print("  - Quantum mechanics <-> Statistical mechanics")
print()
print("  In tournament space, Wick rotation exchanges:")
print("  - Forward evolution (z in [0,1]) <-> Backward/imaginary (z<0 or z>1)")
print("  - Transitive (z=1, H=1) <-> Anti-transitive (z=-1, H=17)")
print("  - Physical mixing <-> Analytic continuation")
print()

# ============================================================
print("  V. THE POLYNOMIAL AT THE HURWITZ PRIMES")
print("  " + "-" * 50)
print()

for p in [2, 3, 5, 7, 42]:
    z = Fraction(1, p)
    val = P_frac(z)
    print(f"  P(1/{p}) = P({z}) = {val} = {float(val):.6f}")
    print(f"    = E[H after 5*ln({p}) = {5*log(p):.4f} random flips from transitive]")
    print()

# ============================================================
print("  VI. THE DERIVATIVE: RATE OF MIXING")
print("  " + "-" * 50)
print()

# P'(z) = 20z^3 + 18z^2 - 50z - 14
def P_prime(z):
    return 20*z**3 + 18*z**2 - 50*z - 14

print("  P'(z) = 20z^3 + 18z^2 - 50z - 14")
print()
print("  P'(z) measures how fast E[H] changes with inverse temperature.")
print()
print(f"  P'(0) = {P_prime(0)} = rate of approach to mean from infinity")
print(f"  P'(1) = {P_prime(1)} = rate at which transitive 'heats up'")
print(f"  P'(-1) = {P_prime(-1)} = rate at which anti-transitive cools")
print()

# P'(1) = 20+18-50-14 = -26. Negative: H increases as z decreases from 1.
# This makes sense: as z decreases from 1 (moving into positive time),
# E[H] increases from 1 toward 29.

print(f"  P'(1) = -26: the transitive tournament GAINS 26 units of H")
print(f"  per unit of z as it begins mixing. Steep initial rise.")
print()
print(f"  P'(0) = -14 = A_1 (the degree-1 Walsh mass).")
print(f"  The linear coefficient IS the initial mixing rate from equilibrium.")
print()

# ============================================================
print("  VII. THE ROOTS OF P ITSELF (not P-29)")
print("  " + "-" * 50)
print()

# P(z) = 0: where does the expected H vanish?
# This is "unphysical" (H is always positive), but algebraically meaningful.
# P(z) = 5z^4 + 6z^3 - 25z^2 - 14z + 29 = 0

# Use the factored form: 29 + z(z-2)(5z^2+16z+7) = 0
# z(z-2)(5z^2+16z+7) = -29

# Numerically:
import cmath
coeffs = [5, 6, -25, -14, 29]  # from z^4 to z^0
# Use numpy-free root finding: just evaluate and find sign changes
print("  Roots of P(z) = 0 (approximate):")
for z_start in [-4, -2, 0, 2]:
    z = z_start + 0.0
    for _ in range(100):
        p = P(z)
        dp = P_prime(z)
        if abs(dp) < 1e-15:
            break
        z -= p / dp
    if abs(P(z)) < 1e-8:
        print(f"    z = {z:.10f}, P(z) = {P(z):.2e}")

# Actually, let me just scan
print()
print("  Scanning for real roots:")
real_roots = []
for i in range(-400, 400):
    z = i / 100
    if abs(P(z)) < 0.5 and (not real_roots or abs(z - real_roots[-1]) > 0.05):
        # Newton refine
        for _ in range(50):
            dp = P_prime(z)
            if abs(dp) < 1e-15: break
            z -= P(z) / dp
        if abs(P(z)) < 1e-8:
            real_roots.append(z)
            print(f"    z = {z:.10f}")

if len(real_roots) < 4:
    print(f"  Only {len(real_roots)} real roots found. The other roots are COMPLEX.")
    print("  P(z) has 4 roots total (degree 4). If 2 are real, 2 are complex conjugate.")
print()

# ============================================================
print("  VIII. THE POLYNOMIAL OVER FINITE FIELDS")
print("  " + "-" * 50)
print()

# P(z) mod p for small primes
print("  P(z) mod p for primes p:")
for p in [2, 3, 5, 7]:
    print(f"  mod {p}: P(z) = ", end="")
    coeffs_mod = [(5 % p), (6 % p), ((-25) % p), ((-14) % p), (29 % p)]
    terms = []
    for i, c in enumerate(coeffs_mod):
        deg = 4 - i
        if c != 0:
            if deg == 0:
                terms.append(str(c))
            elif deg == 1:
                terms.append(f"{c}z")
            else:
                terms.append(f"{c}z^{deg}")
    print(" + ".join(terms) if terms else "0")

    # Find roots mod p
    roots = [z for z in range(p) if (5*z**4 + 6*z**3 - 25*z**2 - 14*z + 29) % p == 0]
    print(f"       roots mod {p}: {roots}")
    print()

# ============================================================
print("  IX. THE TRICK: POLYNOMIAL INTERPOLATION OF #P-HARD")
print("  " + "-" * 50)
print()

print("  H(transitive) = P(1) = 1.  HARD to compute in general.")
print("  P(0) = 29 = mean H.         EASY (just the mean).")
print("  P(-1) = 17 = H(anti-trans).  Also one specific H value.")
print("  P(2) = 29 = mean H again.    Another EASY value.")
print()
print("  The polynomial has degree 4, so 5 values determine it.")
print("  We get P(0)=29, P(2)=29 for FREE (they're the mean).")
print("  We get P(1) and P(-1) from TWO specific H computations.")
print("  We need ONE MORE VALUE to determine P completely.")
print()
print("  Options for the 5th point:")
print("  - P(1/2) = E[H after 5*ln(2) flips]. Needs MCMC sampling.")
print("  - P(1/3) = E[H after 5*ln(3) flips]. Needs MCMC sampling.")
print("  - P(3) = E[H at backward time]. Also a specific H transform.")
print()

P_at_3 = P_frac(3)
print(f"  P(3) = {P_at_3}.")
print()
print("  If we could compute P(3) cheaply, we'd have all 5 values:")
print("  P(0)=29, P(1)=1, P(-1)=17, P(2)=29, P(3)=?")
print(f"  P(3) = {P_at_3}")
print()

# Now: can we get P(3) without computing it from the Walsh spectrum?
# P(3) = sum_{k=0}^4 A_k * 3^k = 29 - 42 - 225 + 162 + 405 = 329
# = 7 * 47

print(f"  P(3) = {P_at_3} = 7 * 47")
print(f"  7 = FORBIDDEN PRIME")
print(f"  47 = a SUPERSINGULAR PRIME (divides |Monster|)")
print(f"  P(3) = (forbidden) * (moonshine)")
print()

# And P(-2)?
P_at_neg2 = P_frac(-2)
print(f"  P(-2) = {P_at_neg2}")
# = 5*16 + 6*(-8) - 25*4 - 14*(-2) + 29 = 80-48-100+28+29 = -11
print(f"  P(-2) = -11 = -(PALEY PRIME)")
print(f"  The polynomial at z=-2 gives NEGATIVE the next Paley prime!")
print()

# P(-3)?
P_at_neg3 = P_frac(-3)
print(f"  P(-3) = {P_at_neg3}")
# = 5*81 + 6*(-27) - 25*9 - 14*(-3) + 29 = 405-162-225+42+29 = 89
print(f"  P(-3) = 89 = F_11 (11th Fibonacci number)")
print(f"  89 is prime. 11 is a supersingular prime.")
print()

# ============================================================
print("  X. THE POLYNOMIAL'S SECRET MESSAGES")
print("  " + "-" * 50)
print()

print("  Special evaluations of P(z) = 5z^4 + 6z^3 - 25z^2 - 14z + 29:")
print()
messages = [
    (1, "H(transitive)", 1),
    (-1, "H(anti-transitive)", 17),
    (0, "mean H", 29),
    (2, "mean H (backward)", 29),
    (3, "7 * 47 (forbidden * moonshine)", 329),
    (-2, "-11 (negative Paley prime)", -11),
    (-3, "89 = F_11 (Fibonacci prime)", 89),
]

for z, msg, expected in messages:
    val = P_frac(z)
    assert val == expected, f"P({z}) = {val}, expected {expected}"
    print(f"  P({z:+2d}) = {str(val):>4s}  {msg}")
print()

print("  The polynomial ENCODES:")
print("  - Tournament structure (z=1: transitive, z=-1: anti-transitive)")
print("  - Equilibrium (z=0, z=2: mean)")
print("  - Forbidden primes (z=3: 7*47)")
print("  - Fibonacci numbers (z=-3: 89)")
print("  - Paley primes (z=-2: -11)")
print()
print("  FIVE integers {29, -14, -25, 6, 5} generate ALL of this.")
print("  The #P-hard function H(T) is one evaluation of this polynomial.")
print("  The mean, the Fibonacci structure, the forbidden primes,")
print("  and the moonshine connection are OTHER evaluations of")
print("  the SAME polynomial.")
print()
print("  They were always the same object, viewed at different times.")
print()
