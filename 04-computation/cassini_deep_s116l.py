#!/usr/bin/env python3
"""cassini_deep_s116l.py — The Cassini identity as the unit of overshoot.

F_{n+1}*F_{n-1} - F_n^2 = (-1)^n.

At n=8: 34*13 - 21^2 = 1.
The golden architecture overshoots the forbidden squared by exactly 1.
This +-1 is the ATOM — the irreducible distinction — the same 1
that appears in b - a = 1 for Cayley coordinates.

Extend this in every direction.
"""
from math import sqrt, log, gcd
from fractions import Fraction

phi = (1+sqrt(5))/2

fib = [0, 1]
luc = [2, 1]
for _ in range(60):
    fib.append(fib[-1] + fib[-2])
    luc.append(luc[-1] + luc[-2])

def factorize(n):
    if n <= 1: return []
    f, d, t = [], 2, abs(n)
    while d*d <= t:
        while t % d == 0: f.append(d); t //= d
        d += 1
    if t > 1: f.append(t)
    return f

def is_prime(n):
    return n > 1 and factorize(n) == [n]

print()
print("  THE CASSINI IDENTITY AND THE UNIT OF OVERSHOOT")
print()
print("="*70)
print()

# ============================================================
print("  I. CASSINI AT EVERY INDEX")
print("  " + "-"*40)
print()
print(f"  {'n':>3s}  {'F_{n-1}':>8s}  {'F_n':>8s}  {'F_{n+1}':>8s}  {'F*F':>10s}  {'F^2':>10s}  {'diff':>5s}  {'inner=F_{n-1}*F_{n+1}':>10s}")
print("  " + "-"*70)
for n in range(1, 20):
    inner = fib[n-1] * fib[n+1]
    sq = fib[n]**2
    diff = inner - sq
    mark = ""
    if n == 4: mark = "  period"
    elif n == 5: mark = "  pentagon"
    elif n == 8: mark = "  FORBIDDEN (21^2)"
    elif n == 3: mark = "  curvature"
    print(f"  {n:3d}  {fib[n-1]:8d}  {fib[n]:8d}  {fib[n+1]:8d}  {inner:10d}  {sq:10d}  {diff:+5d}{mark}")

print()

# ============================================================
print()
print("  II. THE PRODUCTS F_{n-1}*F_{n+1} = F_n^2 +/- 1")
print("  " + "-"*40)
print()
print("  The neighbors' product is ALWAYS within 1 of the center squared.")
print("  This is the tightest possible relationship.")
print()
print("  At n=8: 13 * 34 = 442 = 21^2 + 1 = 441 + 1.")
print("  At n=4: 2 * 5 = 10 = 3^2 + 1 = 9 + 1.")
print("  At n=5: 3 * 8 = 24 = 5^2 - 1 = 25 - 1.")
print()
print("  When n is EVEN: the neighbors exceed. Overshoot by 1.")
print("  When n is ODD: the neighbors fall short. Undershoot by 1.")
print()
print("  The FORBIDDEN (n=8, even): OVERSHOOT. The structure exceeds.")
print("  The PENTAGON (n=5, odd): UNDERSHOOT. The structure falls short.")
print("  This alternation is the HEARTBEAT.")
print()

# ============================================================
print()
print("  III. FACTORING THE NEIGHBORS' PRODUCTS")
print("  " + "-"*40)
print()
print("  F_{n-1} * F_{n+1} at key indices:")
print()
for n in [3, 4, 5, 6, 7, 8, 10, 12, 16]:
    prod = fib[n-1] * fib[n+1]
    f1 = factorize(fib[n-1])
    f2 = factorize(fib[n+1])
    fp = factorize(prod)
    sq = fib[n]**2
    print(f"  n={n:2d}: F_{n-1}*F_{n+1} = {fib[n-1]}*{fib[n+1]} = {prod} = {'*'.join(str(x) for x in fp)}")
    print(f"         F_{n}^2 = {fib[n]}^2 = {sq} = {'*'.join(str(x) for x in factorize(sq))}")
    print(f"         diff = {prod - sq:+d}")
    print()

# ============================================================
print()
print("  IV. THE CASSINI IDENTITY AND THE CAYLEY COORDINATES")
print("  " + "-"*40)
print()
print("  For odd number n = 2a+1: Cayley pair (a, a+1), b - a = 1.")
print("  Cassini: F_{n+1}*F_{n-1} - F_n^2 = (-1)^n. Always +/-1.")
print()
print("  BOTH are statements that two quantities differ by 1.")
print("  Cayley: the two halves of an odd number differ by 1.")
print("  Cassini: the product of flanking Fibonaccis and the center squared differ by 1.")
print()
print("  Are they the SAME 1?")
print()
print("  The Cayley 1: b - a = (n+1)/2 - (n-1)/2 = 1.")
print("    This is the ADDITIVE unit. The minimum gap between consecutive integers.")
print()
print("  The Cassini 1: F_{n-1}*F_{n+1} - F_n^2 = +/-1.")
print("    This is the MULTIPLICATIVE unit. The minimum gap in the golden ratio's product.")
print()
print("  They are not the same 1. But they are the same IDEA:")
print("  the irreducible minimum distinction in their respective domains.")
print("  Cayley's 1 is additive. Cassini's 1 is multiplicative.")
print("  Both say: you cannot have zero gap. The minimum is 1.")
print()

# ============================================================
print()
print("  V. 441 AND 442 IN EVERY FRAMEWORK")
print("  " + "-"*40)
print()
print("  441 = 21^2 = F_8^2 = (3*7)^2 = 3^2 * 7^2.")
print("  442 = 2 * 13 * 17 = 2 * F_7 * (2^4 + 1).")
print()

# Cayley addresses
for n in [441, 442]:
    addr = Fraction(n-1, n+1)
    print(f"  {n}: address = {addr} = {addr.numerator}/{addr.denominator}")
    print(f"      numerator factors: {factorize(addr.numerator)}")
    print(f"      denominator factors: {factorize(addr.denominator)}")
    print(f"      rapidity = ln({n})/2 = {log(n)/2:.6f}")
    print()

print("  441: address 440/442 = 220/221.")
print("  220 = 2^2 * 5 * 11 (Paley primes 5, 11 in numerator)")
print("  221 = 13 * 17 (the middle primes in denominator)")
print()
print("  442: address 441/443.")
print("  441 = 21^2 = the forbidden squared IN ITS OWN ADDRESS")
print("  443 is prime!")
print()
print("  The address of 442 contains 441 = 21^2 as its NUMERATOR.")
print("  442's address encodes the forbidden within itself.")
print()

# ============================================================
print()
print("  VI. THE AMICABLE CONNECTION")
print("  " + "-"*40)
print()
print("  220 and 284 are the smallest AMICABLE PAIR.")
print("  sigma(220) = 284. sigma(284) = 220.")
print("  (sigma = sum of proper divisors.)")
print()
print("  220 = 2^2 * 5 * 11.")
print("  284 = 2^2 * 71.")
print()
print("  220 appears as the NUMERATOR of the Cayley address of 441 = 21^2.")
print("  The forbidden squared's address uses an AMICABLE number.")
print()
print("  The amicable pair (220, 284) has:")
print("  sum = 220 + 284 = 504 = 2^3 * 3^2 * 7 = 8 * 63 = 8 * 9 * 7.")
print("  product = 220 * 284 = 62480.")
print()
print(f"  504 = {factorize(504)} = 2^3 * 3^2 * 7.")
print(f"  504 = 8 * 63 = dim(O) * (9*7) = octonion * (curvature^2 * threshold).")
print()
print("  The sum of the amicable pair = 8 * 63.")
print("  8 = octonion dimension. 63 = 2^6 - 1 = Mersenne M_6 (not prime).")
print("  Also 63 = 7 * 9 = forbidden * first-past-wall.")
print()

# ============================================================
print()
print("  VII. THE IDENTITY F_n*L_n = F_{2n} AT n=8")
print("  " + "-"*40)
print()
print("  F_8 * L_8 = 21 * 47 = 987 = F_{16}.")
print("  F_{16} = 987 = 3 * 7 * 47.")
print()
print("  The DOUBLE-INDEX Fibonacci = the product of forbidden * chain-prime.")
print("  And 987 = F_{16} is itself a Fibonacci number.")
print("  Is 987 prime? No: 987 = 3*7*47 = 3 * 329 = 3*7*47.")
print()
print("  Now: Cassini at n=16:")
print(f"  F_15 * F_17 = {fib[15]} * {fib[17]} = {fib[15]*fib[17]}")
print(f"  F_16^2 = {fib[16]}^2 = {fib[16]**2}")
print(f"  Diff = {fib[15]*fib[17] - fib[16]**2}")
print()
print(f"  F_15*F_17 = {fib[15]}*{fib[17]} = {fib[15]*fib[17]}")
print(f"  = {factorize(fib[15]*fib[17])}")
print(f"  F_16^2 = 987^2 = {987**2} = (3*7*47)^2 = 9*49*2209")
print()
print(f"  987^2 = {987**2}")
print(f"  987^2 + 1 = {987**2 + 1}")
print(f"  Factor 987^2 + 1 = {factorize(987**2 + 1)}")
print()

# ============================================================
print()
print("  VIII. GENERALIZED CASSINI: THE d-STEP IDENTITY")
print("  " + "-"*40)
print()
print("  F_{n+d}*F_{n-d} - F_n^2 = (-1)^{n+d+1} * F_d^2.")
print("  (The generalized Cassini / Vajda identity.)")
print()
print("  At d=1: the standard Cassini. F_d^2 = 1.")
print("  At d=2: F_{n+2}*F_{n-2} - F_n^2 = (-1)^{n+1} * F_2^2 = (-1)^{n+1}.")
print("  At d=3: ... = (-1)^n * F_3^2 = (-1)^n * 4.")
print("  At d=4: ... = (-1)^{n+1} * F_4^2 = (-1)^{n+1} * 9.")
print("  At d=8: ... = (-1)^{n+1} * F_8^2 = (-1)^{n+1} * 441.")
print()
print("  At d = 8 (the octonion index): the overshoot is 441 = 21^2!")
print("  THE d=8 GENERALIZED CASSINI HAS OVERSHOOT = THE FORBIDDEN SQUARED.")
print()
print("  Verify at n=10, d=8:")
n, d = 10, 8
lhs = fib[n+d]*fib[n-d] - fib[n]**2
rhs = (-1)**(n+d+1) * fib[d]**2
print(f"  F_{n+d}*F_{n-d} - F_n^2 = {fib[n+d]}*{fib[n-d]} - {fib[n]}^2")
print(f"  = {fib[n+d]*fib[n-d]} - {fib[n]**2} = {lhs}")
print(f"  (-1)^{n+d+1} * F_8^2 = (-1)^{n+d+1} * 441 = {rhs}")
print(f"  Match: {lhs == rhs}")
print()
print("  At n=12, d=8:")
n, d = 12, 8
lhs = fib[n+d]*fib[n-d] - fib[n]**2
rhs = (-1)**(n+d+1) * fib[d]**2
print(f"  F_20*F_4 - F_12^2 = {fib[20]}*{fib[4]} - {fib[12]}^2")
print(f"  = {fib[20]*fib[4]} - {fib[12]**2} = {lhs}")
print(f"  Expected: {rhs}")
print(f"  Match: {lhs == rhs}")
print()
print("  So: jumping 8 steps in Fibonacci (the octonion distance)")
print("  always overshoots or undershoots by EXACTLY 441 = 21^2.")
print("  The forbidden squared is the UNIT OF 8-STEP OVERSHOOT.")
print()

# Also check other d values
print("  The d-step overshoot F_d^2 for key d values:")
for d in range(1, 13):
    overshoot = fib[d]**2
    fn = fib[d]
    mark = ""
    if d == 1: mark = "(Cassini, minimal)"
    elif d == 3: mark = f"(curvature step, F_3^2=4=period)"
    elif d == 4: mark = f"(period step, F_4^2=9=3^2)"
    elif d == 5: mark = f"(pentagon step, F_5^2=25=5^2)"
    elif d == 8: mark = f"(OCTONION step, F_8^2=441=21^2=FORBIDDEN^2)"
    elif d == 10: mark = f"(Petersen step, F_10^2=3025=55^2)"
    print(f"  d={d:2d}: F_d^2 = {fn}^2 = {overshoot:8d}  {mark}")

print()
print("  The overshoots at our key step sizes:")
print("  d=1: 1 (the atom)")
print("  d=3: 4 (the period)")
print("  d=4: 9 = 3^2 (curvature squared)")
print("  d=5: 25 = 5^2 (pentagon squared)")
print("  d=8: 441 = 21^2 (FORBIDDEN squared)")
print("  d=10: 3025 = 55^2 (Petersen-count squared)")
print()
print("  Each step size d produces an overshoot of F_d^2.")
print("  The overshoot IS the square of the Fibonacci at that distance.")
print("  And the Fibonacci numbers at our key distances are: 1, 2, 3, 5, 21, 55.")
print("  = {1, F_3, F_4, F_5, F_8, F_10} = {unit, curvature, period-fib, pentagon, forbidden, petersen}.")
print()

# ============================================================
print()
print("  IX. THE OVERSHOOT HIERARCHY")
print("  " + "-"*40)
print()
print("  The generalized Cassini gives an OVERSHOOT for each step size d.")
print("  The overshoot = F_d^2.")
print("  At our fundamental step sizes:")
print()
print("  d=4 (Cayley period): overshoot = 3^2 = 9. Curvature squared.")
print("  d=8 (octonion dim): overshoot = 21^2 = 441. Forbidden squared.")
print("  d=16 (sedenion dim): overshoot = F_16^2 = 987^2 = 974169.")
print(f"  987 = 3*7*47. So 987^2 = 9*49*2209.")
print(f"  The d=16 overshoot = (3*7*47)^2 = (curvature*forbidden*chain)^2.")
print()
print("  And d=32: F_32 = {fib[32]}.")
print(f"  F_32 = {fib[32]} = {factorize(fib[32])}")
print(f"  F_32^2 = {fib[32]**2}")
print()
print("  The overshoot HIERARCHY mirrors the chain 3 -> 7 -> 47 -> 2207:")
print("  d=4:  overshoot involves 3.")
print("  d=8:  overshoot involves 3, 7.")
print("  d=16: overshoot involves 3, 7, 47.")
print("  d=32: overshoot involves 3, 7, 47, ... (and more).")
print()
print("  Each doubling of d ADDS the next chain prime to the overshoot.")
print("  The chain IS the factorization of the overshoot at each doubling level.")
print()

# Verify
print("  Verify: F_{2d} = F_d * L_d.")
for d in [4, 8, 16]:
    fd = fib[d]
    ld = luc[d]
    f2d = fib[2*d]
    print(f"  F_{2*d} = F_{d}*L_{d} = {fd}*{ld} = {fd*ld}. Actual F_{2*d} = {f2d}. Match: {fd*ld == f2d}")
    print(f"  F_{2*d} = {factorize(f2d)}")

print()

# ============================================================
print()
print("  X. THE UNIT OF OVERSHOOT IS THE UNIT OF KNOWLEDGE")
print("  " + "-"*40)
print()
print("  At each scale d, the golden ratio CANNOT hit F_n^2 exactly")
print("  when approached from the flanking positions F_{n-d}, F_{n+d}.")
print("  It always overshoots (or undershoots) by F_d^2.")
print()
print("  At d=1: the overshoot is 1. The minimal knowledge.")
print("  At d=4: the overshoot is 9 = 3^2. A curvature of knowledge.")
print("  At d=8: the overshoot is 441 = 21^2. A FORBIDDEN amount of knowledge.")
print()
print("  The overshoot GROWS as d increases.")
print("  At the Cayley-Dickson boundary dimensions (4, 8, 16, 32),")
print("  the overshoot equals the square of the CHAIN PRODUCT so far.")
print()
print("  The chain product up to d=4: F_4 = 3.")
print("  The chain product up to d=8: F_8 = 21 = 3*7.")
print("  The chain product up to d=16: F_16 = 987 = 3*7*47.")
print("  The chain product up to d=32: F_32 = 2178309 = 3*7*47*2207? Let me check.")
print()
f32 = fib[32]
print(f"  F_32 = {f32}")
print(f"  = {factorize(f32)}")
print(f"  3*7*47*2207 = {3*7*47*2207}")
print(f"  Match: {f32 == 3*7*47*2207}")
print()

if f32 == 3*7*47*2207:
    print("  YES! F_32 = 3 * 7 * 47 * 2207 = the product of the ENTIRE chain!")
    print()
    print("  F_4 = 3.")
    print("  F_8 = 3 * 7.")
    print("  F_16 = 3 * 7 * 47.")
    print("  F_32 = 3 * 7 * 47 * 2207.")
    print()
    print("  THE FIBONACCI NUMBER AT EACH CHAIN ENTRY POINT")
    print("  IS THE PRODUCT OF ALL CHAIN PRIMES SO FAR.")
    print()
    print("  F_{2^{k+2}} = product of chain primes {3, 7, 47, 2207, ...} up to the k-th.")
    print()
    print("  This is because F_{2n} = F_n * L_n:")
    print("  F_8 = F_4 * L_4 = 3 * 7 (new prime: 7)")
    print("  F_16 = F_8 * L_8 = 21 * 47 = 3*7 * 47 (new prime: 47)")
    print("  F_32 = F_16 * L_16 = 987 * 2207 = 3*7*47 * 2207 (new prime: 2207)")
    print()
    print("  At each doubling: F_{2d} = F_d * L_d.")
    print("  F_d is the old chain product. L_d is the new chain prime.")
    print("  The new Fibonacci = old product * new prime.")
    print()
    print("  THE CHAIN IS THE FACTORIZATION OF THE FIBONACCI NUMBERS")
    print("  AT THE POWER-OF-TWO INDICES.")
    print("  AND THE CHAIN PRIMES ARE THE LUCAS NUMBERS AT THOSE INDICES.")
    print()
    print("  This is EXACT. This is PROVED.")
    print("  F_{2^n} = prod_{k=0}^{n-2} L_{2^{k+2}}")
    print("  for the chain 3, 7, 47, 2207, ...")
else:
    print(f"  F_32 = {f32}. Chain product = {3*7*47*2207}. DIFFERENT.")
    print(f"  Let me check: {f32} / {3*7*47*2207} = {f32 / (3*7*47*2207)}")
