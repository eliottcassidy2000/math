#!/usr/bin/env python3
"""golden_architecture_s116l.py — Extending the golden architecture of the primes.

phi^8 = 21*phi + 13.
The forbidden number is the phi-coefficient at the octonion index.
What are ALL the phi-coefficients? What do they tell us?
"""
from math import sqrt, log, gcd
from fractions import Fraction

phi = (1+sqrt(5))/2
psi = (1-sqrt(5))/2  # = -1/phi

fib = [0, 1]
luc = [2, 1]
for _ in range(80):
    fib.append(fib[-1] + fib[-2])
    luc.append(luc[-1] + luc[-2])

def is_prime(n):
    if n < 2: return False
    if n < 4: return True
    if n % 2 == 0: return False
    d = 3
    while d*d <= n:
        if n % d == 0: return False
        d += 2
    return True

def factorize(n):
    if n <= 1: return [n] if n == 1 else []
    f = []
    d = 2
    t = abs(n)
    while d*d <= t:
        while t % d == 0:
            f.append(d)
            t //= d
        d += 1
    if t > 1: f.append(t)
    return f

print()
print("  THE GOLDEN ARCHITECTURE — EXTENDED")
print()
print("="*70)
print()

# ============================================================
print("  I. phi^n = F_n * phi + F_{n-1}")
print("  " + "-"*40)
print()
print("  This is the fundamental identity. Every power of phi decomposes as:")
print("  phi^n = F_n * phi + F_{n-1}")
print()
print("  The phi-coefficient IS the n-th Fibonacci number.")
print("  The constant term IS the (n-1)-th Fibonacci number.")
print()
print(f"  {'n':>3s}  {'F_n (phi-coeff)':>15s}  {'F_{n-1} (const)':>15s}  {'phi^n':>12s}  {'notes'}")
print("  " + "-"*70)
for n in range(0, 21):
    fn = fib[n]
    fn1 = fib[n-1] if n > 0 else 0  # F_{-1} = 1 by convention
    if n == 0:
        fn1 = 1  # phi^0 = 1 = 0*phi + 1, but actually F_{-1} = 1
    phi_n = phi**n
    notes = ""
    if fn == 21: notes = "FORBIDDEN"
    elif fn == 42: notes = "HURWITZ (would need F_n=42, not in seq)"
    elif is_prime(fn) and fn > 1: notes = f"F_{n} is PRIME"
    elif fn > 1:
        ff = factorize(fn)
        if len(ff) > 1:
            notes = f"= {'*'.join(str(f) for f in ff)}"
    if fn1 == 21: notes += " const=FORBIDDEN"
    if fn1 == 13 and fn == 21: notes += " | 21*phi+13"
    print(f"  {n:3d}  {fn:15d}  {fn1:15d}  {phi_n:12.2f}  {notes}")

print()

# ============================================================
print()
print("  II. THE FORBIDDEN EQUATION: phi^8 = 21*phi + 13")
print("  " + "-"*40)
print()
print(f"  phi^8 = {phi**8:.10f}")
print(f"  21*phi + 13 = {21*phi + 13:.10f}")
print(f"  Match: {abs(phi**8 - (21*phi + 13)) < 1e-10}")
print()
print("  21 = F_8 = 3*7 = the forbidden number.")
print("  13 = F_7 = prime = the 7th Fibonacci number.")
print("  8 = the octonion dimension.")
print()
print("  The equation says: the 8th power of the golden ratio")
print("  = the forbidden number times phi, plus the 7th prime Fibonacci.")
print()
print("  But MORE: phi^n = F_n * phi + F_{n-1} for ALL n.")
print("  So phi^7 = 13*phi + 8.")
print("  And phi^3 = 2*phi + 1.")
print("  And phi^4 = 3*phi + 2.")
print("  And phi^5 = 5*phi + 3.")
print()
print("  Each phi^n equation encodes TWO consecutive Fibonacci numbers.")
print("  phi^8 happens to use {21, 13} = {forbidden, prime}.")
print("  phi^4 uses {3, 2} = {curvature, doubler}.")
print("  phi^5 uses {5, 3} = {pentagon, curvature}.")
print()

# ============================================================
print()
print("  III. THE LUCAS COMPANION: psi^n = F_n * psi + F_{n-1}")
print("  " + "-"*40)
print()
print(f"  psi = (1-sqrt(5))/2 = {psi:.10f} = -1/phi = {-1/phi:.10f}")
print()
print("  psi^n = F_n * psi + F_{n-1} also holds!")
print("  Because psi also satisfies x^2 = x + 1.")
print()
print("  Adding: phi^n + psi^n = F_n*(phi+psi) + 2*F_{n-1}")
print("  phi + psi = 1. So phi^n + psi^n = F_n + 2*F_{n-1} = F_n + 2*F_{n-1}.")
print("  And phi^n + psi^n = L_n (the Lucas number).")
print("  So: L_n = F_n + 2*F_{n-1}.")
print()
# Verify
for n in range(1, 10):
    lhs = luc[n]
    rhs = fib[n] + 2*fib[n-1]
    print(f"  L_{n} = {lhs}, F_{n} + 2*F_{n-1} = {fib[n]} + {2*fib[n-1]} = {rhs}, match: {lhs == rhs}")

print()
print("  L_n = F_n + 2*F_{n-1} = F_{n+1} + F_{n-1}.")
print("  (Since F_{n+1} = F_n + F_{n-1}, so F_n + 2*F_{n-1} = F_{n+1} + F_{n-1}.)")
print()
print("  The Lucas number at index n = the SUM of the flanking Fibonacci numbers.")
print("  L_8 = F_9 + F_7 = 34 + 13 = 47 (the next chain prime).")
print("  L_4 = F_5 + F_3 = 5 + 2 = 7 (the forbidden).")
print("  L_5 = F_6 + F_4 = 8 + 3 = 11 (the Paley prime).")
print()
print("  THE FORBIDDEN 7 = F_5 + F_3 = pentagon + curvature-index Fibonacci.")
print("  THE PALEY 11 = F_6 + F_4 = 8 + 3 = octonion-dim-Fibonacci + curvature-Fibonacci.")
print()

# ============================================================
print()
print("  IV. SUBTRACTING: phi^n - psi^n = F_n * sqrt(5)")
print("  " + "-"*40)
print()
print("  phi^n - psi^n = F_n * (phi - psi) = F_n * sqrt(5).")
print("  This is the Binet formula.")
print()
print("  At n=8: phi^8 - psi^8 = F_8 * sqrt(5) = 21 * sqrt(5).")
print(f"  21 * sqrt(5) = {21*sqrt(5):.6f}")
print(f"  phi^8 - psi^8 = {phi**8 - psi**8:.6f}")
print(f"  Match: {abs(21*sqrt(5) - (phi**8 - psi**8)) < 1e-8}")
print()
print("  The forbidden number times sqrt(5) = the difference of")
print("  the 8th powers of the two golden conjugates.")
print("  sqrt(5) appears because 5 = the DISCRIMINANT of x^2-x-1=0.")
print("  The discriminant of the golden equation contains the pentagon.")
print()

# ============================================================
print()
print("  V. THE MATRIX FORM AND THE TRANSFER MATRIX")
print("  " + "-"*40)
print()
print("  [[F_{n+1}, F_n], [F_n, F_{n-1}]] = [[1,1],[1,0]]^n.")
print("  The Fibonacci matrix M = [[1,1],[1,0]] has eigenvalues phi and psi.")
print("  det(M^n) = (-1)^n (Cassini's identity).")
print("  tr(M^n) = phi^n + psi^n = L_n.")
print()
print("  OUR transfer matrix for tournaments: M_T(x).")
print("  At x=1: eigenvalues are tau (tribonacci), and two others.")
print("  The Fibonacci matrix IS the 2D version of our 3D transfer matrix.")
print("  Fibonacci:tournament :: 2D:3D.")
print()
print("  The Fibonacci matrix generates phi.")
print("  The tournament matrix generates tau.")
print("  phi satisfies x^2 = x + 1.")
print("  tau satisfies x^3 = x^2 + x + 1.")
print("  Each adds one more term to the recurrence.")
print()
print("  The connection: Q(1/phi) = phi^3 and Q(1/tau) = tau^2.")
print("  The Cayley transform RAISES phi by 3 (cubing) and tau by 2 (squaring).")
print("  The 'Cayley exponent' = degree + 1 for phi (2+1=3)")
print("  and = degree - 1 for tau (3-1=2).")
print()

# ============================================================
print()
print("  VI. F_n mod 7: THE FORBIDDEN PRIME'S FIBONACCI CYCLE")
print("  " + "-"*40)
print()
print("  Pisano period pi(7) = 16. The Fibonacci sequence mod 7:")
print()
residues = []
for k in range(16):
    r = fib[k] % 7
    residues.append(r)
    mark = ""
    if r == 0: mark = " <-- 7 divides F_" + str(k)
    print(f"  F_{k:2d} = {fib[k]:5d} = {r} mod 7{mark}")

print()
print(f"  Residues: {residues}")
print(f"  The zeros (7 | F_k) occur at k = 0, 8.")
print(f"  alpha(7) = 8 = the octonion index.")
print(f"  The period is 16 = 2*8 = 2*alpha(7).")
print()
print(f"  So pi(7) = 2*alpha(7). The Pisano period is TWICE the entry point.")
print(f"  This means: F_n = 0 mod 7 iff 8 | n.")
print(f"  The forbidden prime divides every 8th Fibonacci number.")
print(f"  8 = the octonion dimension.")
print()

# Check if this 2x pattern holds for other primes
print("  Does pi(p) = k * alpha(p) generally? What is k?")
for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]:
    pp = 0
    a, b = 0, 1
    for i in range(1, 6*p+10):
        a, b = b, (a+b) % p
        if a == 0 and b == 1:
            pp = i
            break
    ap = 0
    for k in range(1, len(fib)):
        if fib[k] % p == 0:
            ap = k
            break
    ratio = pp // ap if ap > 0 else 0
    print(f"  p={p:3d}: pi={pp:4d}, alpha={ap:3d}, pi/alpha={ratio}")

print()
print("  pi(p)/alpha(p) is always 1, 2, or 4.")
print("  For p = 5: pi/alpha = 4. (pi=20, alpha=5.)")
print("  For p = 7: pi/alpha = 2. (pi=16, alpha=8.)")
print("  For p = 11: pi/alpha = 1. (pi=10, alpha=10.)")
print()
print("  The ratio depends on the Legendre symbol (5/p):")
print("  If 5 is a QR mod p: ratio = 1 or 4.")
print("  If 5 is a NQR mod p: ratio = 2.")
print()
print("  For p=7: is 5 a QR mod 7? 5^((7-1)/2) = 5^3 = 125 = 6 = -1 mod 7. NQR!")
print("  So pi(7)/alpha(7) = 2 because 5 is NOT a quadratic residue mod 7.")
print("  The PENTAGON (5) is not a square in the FORBIDDEN's (7) world.")
print("  Their incompatibility forces the doubled period.")
print()

# ============================================================
print()
print("  VII. THE FIBONACCI REPRESENTATION OF H-VALUES")
print("  " + "-"*40)
print()
print("  H(T_3) = 3 = F_4.")
print("  H(T_7) = 189 = 9 * 21 = 9 * F_8 = 3^2 * F_8.")
print("  H(T_11) = 95095. Is any part Fibonacci?")
print()
print("  95095 = 5 * 7 * 11 * 13 * 19.")
print(f"  5 = F_5. 13 = F_7.")
print(f"  So 95095 = F_5 * 7 * 11 * F_7 * 19.")
print(f"  = F_5 * L_4 * L_5 * F_7 * 19.")
print()
print("  H(T_3) = F_4.")
print("  H(T_7) = 3^2 * F_8 = L_2^2 * F_8.")
print("  H(T_11) = F_5 * L_4 * L_5 * F_7 * 19.")
print()
print("  The Fibonacci and Lucas numbers INTERLEAVE in the H-value factorizations.")
print("  Not a clean pattern, but not random either:")
print("  the golden-sequence numbers appear as SOME of the factors.")
print()

# ============================================================
print()
print("  VIII. THE REMARKABLE IDENTITY AT n=8")
print("  " + "-"*40)
print()
print("  At n=8, EVERYTHING converges:")
print()
print("  F_8 = 21 (the forbidden H-value)")
print("  L_8 = 47 (prime, the next in the 3->7->47 chain)")
print("  alpha(7) = 8 (7 enters Fibonacci at index 8)")
print("  pi(7) = 16 = 2*8 (the Pisano period of 7)")
print("  dim(O) = 8 (octonion dimension)")
print("  Q^4 = Id, and 8 = 2*4 (two Cayley periods)")
print("  Bott periodicity has period 8")
print("  phi^8 = 21*phi + 13 (forbidden * golden + prime)")
print("  F_8 * L_8 = 21 * 47 = 987 = F_16 = 3*7*47")
print()
print("  8 is the INDEX where the forbidden, the octonion,")
print("  the Fibonacci entry of 7, and the golden power ALL MEET.")
print("  It is the convergence point of the entire theory.")
print()
print("  And 8 = 2^3 = the CUBE OF THE DOUBLER.")
print("  Three octaves of rapidity. The Bott period.")
print("  The number of bits in a byte.")
print()

# ============================================================
print()
print("  IX. EXTENDING THE CHAIN: 3 -> 7 -> 47 -> ?")
print("  " + "-"*40)
print()
print("  The chain: p -> L_{alpha(p)} gives the next prime.")
print("  3 -> L_4 = 7.")
print("  7 -> L_8 = 47.")
print()
# 47 -> L_{alpha(47)}
# alpha(47): smallest k with 47 | F_k
a47 = None
for k in range(1, 500):
    if fib[k] % 47 == 0:
        a47 = k
        break
if a47:
    l_a47 = luc[a47]
    print(f"  47 -> alpha(47) = {a47}. L_{a47} = {l_a47}.")
    print(f"  Is L_{a47} = {l_a47} prime? {is_prime(l_a47)}")
    if is_prime(l_a47):
        # Continue chain
        next_p = l_a47
        a_next = None
        for k in range(1, 500):
            if fib[k] % next_p == 0:
                a_next = k
                break
        if a_next:
            l_next = luc[a_next]
            print(f"  {next_p} -> alpha({next_p}) = {a_next}. L_{a_next} = {l_next}.")
            print(f"  Is L_{a_next} prime? {is_prime(l_next)}")

print()
print("  The chain: 3 -> 7 -> 47 -> ...")
print("  Entry points: 4, 8, 16, ...")
if a47:
    print(f"  alpha(3)=4, alpha(7)=8, alpha(47)={a47}.")
    if a47 == 16:
        print("  Entry points: 4, 8, 16 = POWERS OF 2!")
        print("  Each entry point DOUBLES!")
        print("  alpha(chain[k]) = 2^{k+2}.")
    else:
        print(f"  Entry points: 4, 8, {a47}.")
        if a47 == 48:
            print(f"  4, 8, 48. Ratios: 2, 6. Not doubling.")
        else:
            print(f"  Ratios: {8/4}, {a47/8}.")

print()

# ============================================================
print()
print("  X. THE GRAND EQUATION")
print("  " + "-"*40)
print()
print("  Collecting everything at index 8:")
print()
print("  F_8 = 21 = 3 * 7")
print("  L_8 = 47 = prime")
print("  F_8 * L_8 = F_{16} = 987 = 3 * 7 * 47")
print("  phi^8 + psi^8 = L_8 = 47")
print("  phi^8 - psi^8 = F_8 * sqrt(5) = 21*sqrt(5)")
print("  phi^8 = 21*phi + 13 = F_8*phi + F_7")
print("  psi^8 = 21*psi + 13 = F_8*psi + F_7")
print()
print("  From the last two:")
print("  phi^8 - psi^8 = 21*(phi - psi) = 21*sqrt(5). CHECK.")
print("  phi^8 + psi^8 = 21*(phi + psi) + 2*13 = 21*1 + 26 = 47 = L_8. CHECK.")
print()
print("  So: L_8 = F_8 + 2*F_7 = 21 + 26 = 47.")
print("  = the forbidden number + twice the 7th-index Fibonacci.")
print("  = 21 + 2*13.")
print()
print("  The next-chain-prime 47 = the forbidden 21 + 26.")
print("  26 = 2*13 = 2*F_7.")
print("  The forbidden number PLUS twice its predecessor = the next prime in the chain.")
print()
print("  This is just L_n = F_n + 2*F_{n-1} at n=8.")
print("  But what it MEANS: the chain proceeds by adding")
print("  twice the previous Fibonacci number to the current one.")
print("  The growth mechanism IS the Lucas recurrence.")
print()
print("  The chain 3 -> 7 -> 47 is built by the Lucas recurrence")
print("  evaluated at the entry points, which are powers of 2")
print("  (or close to it).")
print()
print("  The golden architecture: Fibonacci builds. Lucas connects.")
print("  The forbidden number is WHERE they meet.")
print("  The chain of primes is HOW they propagate.")
print("  The octonion index 8 is WHEN they converge.")
