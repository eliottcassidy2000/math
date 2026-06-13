#!/usr/bin/env python3
"""cayley_coordinates_s116g.py — Every odd number n has Cayley coordinates (a, b) = ((n-1)/2, (n+1)/2).
a + b = n. b - a = 1. a/b = Cayley address. Q(a/b) = n.

Turn over every stone.
"""
from math import sqrt, log, pi, exp, atanh, gcd, factorial
from fractions import Fraction
from collections import defaultdict

def is_prime(n):
    if n < 2: return False
    for d in range(2, int(sqrt(n))+1):
        if n % d == 0: return False
    return True

def factorize(n):
    if n <= 1: return []
    factors = []
    d = 2
    while d*d <= n:
        while n % d == 0:
            factors.append(d)
            n //= d
        d += 1
    if n > 1:
        factors.append(n)
    return factors

print()
print("  CAYLEY COORDINATES OF THE ODD NUMBERS")
print("  turning over every stone")
print()
print("="*70)
print()

# ============================================================
print("  THE TABLE")
print("  " + "-"*40)
print()
print("  n     a=(n-1)/2   b=(n+1)/2   a factors       b factors       a*b")
print("  " + "-"*75)
for n in range(3, 62, 2):
    a = (n-1)//2
    b = (n+1)//2
    af = "*".join(str(f) for f in factorize(a)) if a > 1 else "1"
    bf = "*".join(str(f) for f in factorize(b)) if b > 1 else "1"
    ab = a*b
    mark = ""
    if n == 7: mark = " <-- FORBIDDEN"
    elif n == 21: mark = " <-- FORBIDDEN"
    elif is_prime(n): mark = " prime"
    print(f"  {n:3d}    {a:4d}         {b:4d}        {af:<14s}  {bf:<14s}  {ab:6d}{mark}")

print()

# ============================================================
print()
print("  STONE 1: a*b = (n-1)(n+1)/4 = (n^2-1)/4")
print("  " + "-"*40)
print()
print("  For odd n, a*b = (n-1)/2 * (n+1)/2 = (n^2-1)/4.")
print("  This is always an integer (since n is odd, n^2-1 is divisible by 8,")
print("  so (n^2-1)/4 is even).")
print()
print("  n     a*b = (n^2-1)/4")
for n in range(3, 32, 2):
    ab = (n*n - 1) // 4
    print(f"  {n:3d}    {ab:5d}", end="")
    if ab > 1:
        print(f" = {'*'.join(str(f) for f in factorize(ab))}", end="")
    print()

print()
print("  The sequence a*b: 2, 6, 12, 20, 30, 42, 56, 72, 90, 110, 132, 156, ...")
print("  These are the PRONIC NUMBERS! a*b = k*(k+1) where k = (n-1)/2.")
print("  a*b = a*(a+1) = the a-th pronic number = 2*T_a (twice triangular).")
print()
print("  WE ALREADY SAW PRONIC NUMBERS in the rapidity composition algebra!")
print("  The number of rapidity decompositions of interval c = tau(c(c+1))/2.")
print("  And here: the Cayley coordinate product = a-th pronic number.")
print()
print("  So the Cayley product of n = the rapidity decomposition count of a = (n-1)/2.")
print()

# ============================================================
print()
print("  STONE 2: WHEN IS a PRIME? WHEN IS b PRIME?")
print("  " + "-"*40)
print()
print("  a = (n-1)/2 is prime iff n = 2p+1 for some prime p.")
print("  These n are called SAFE PRIMES (if n itself is prime too,")
print("  then p is a Sophie Germain prime).")
print()
print("  b = (n+1)/2 is prime iff n = 2q-1 for some prime q.")
print()

a_prime = []
b_prime = []
both_prime = []
neither = []
for n in range(3, 100, 2):
    a, b = (n-1)//2, (n+1)//2
    ap = is_prime(a)
    bp = is_prime(b)
    if ap and bp: both_prime.append(n)
    elif ap: a_prime.append(n)
    elif bp: b_prime.append(n)
    else: neither.append(n)

print(f"  a prime only:     {a_prime[:12]}")
print(f"  b prime only:     {b_prime[:12]}")
print(f"  BOTH a,b prime:   {both_prime[:12]}")
print(f"  neither prime:    {neither[:12]}")
print()
print("  When BOTH a and b are prime:")
for n in both_prime[:10]:
    a, b = (n-1)//2, (n+1)//2
    print(f"    n={n}: a={a}, b={b} (twin primes shifted! b-a=1, so (a,b) are CONSECUTIVE integers, both prime? No, they can't both be prime unless one is 2.)")

print()
print("  Wait: a and b differ by 1. Two consecutive integers can't")
print("  both be prime unless one of them is 2.")
print("  So a=2, b=3 -> n=5. Or a=prime, b=a+1=even -> b not prime (unless b=2).")
print("  b=2 -> a=1 -> n=3.")
print()
print("  Correcting: for n >= 5, at most ONE of a,b is prime.")
print("  (Because consecutive integers > 2 can't both be prime.)")
print()
print(f"  Verify: both_prime = {both_prime}")
print("  Only n=5 (a=2,b=3) and n=3 (a=1,b=2). For n=3, a=1 is not prime.")
print("  So only n=5 has both coordinates prime (and only because 2,3 are")
print("  the unique pair of consecutive primes).")
print()

# ============================================================
print()
print("  STONE 3: THE PARITY OF a")
print("  " + "-"*40)
print()
print("  n is odd. a = (n-1)/2.")
print("  n = 1 mod 4 => n-1 = 0 mod 4 => a even.")
print("  n = 3 mod 4 => n-1 = 2 mod 4 => a odd.")
print()
print("  So the parity of a encodes n mod 4:")
print("  a even <=> n = 1 mod 4")
print("  a odd  <=> n = 3 mod 4")
print()
print("  For tournaments: H(T) is always odd (Redei).")
print("  The Cayley coordinate a of an odd H-value tells us H mod 4:")
print("  a even => H = 1 mod 4")
print("  a odd  => H = 3 mod 4")
print()
print("  The forbidden H=7: a=3 (odd) => 7 = 3 mod 4. Correct.")
print("  The forbidden H=21: a=10 (even) => 21 = 1 mod 4. Correct.")
print()
print("  The two forbidden numbers have DIFFERENT parities of a!")
print("  7 has odd a (3). 21 has even a (10).")
print("  One lives in H = 3 mod 4, the other in H = 1 mod 4.")
print("  The forbidden numbers span BOTH residue classes mod 4.")
print()

# ============================================================
print()
print("  STONE 4: TRIANGULAR NUMBERS AS CAYLEY PRODUCTS")
print("  " + "-"*40)
print()
print("  a*b = a*(a+1) = 2*T_a where T_a = a(a+1)/2 = triangular number.")
print()
print("  So a*b = 2*T_a = number of edges in K_{a+1} (complete graph).")
print("  The Cayley coordinate product of n = edges of K_{(n+1)/2}.")
print()
print("  For n=7: a*b = 3*4 = 12 = edges of K_4 = edges of the TETRAHEDRON.")
print("  For n=21: a*b = 10*11 = 110 = edges of K_11.")
print("  For n=11: a*b = 5*6 = 30 = edges of K_6 = edges of the OCTAHEDRON... ")
print("  (K_6 has 15 edges, not 30. Actually K_{a+1} = K_6 has C(6,2) = 15 edges.)")
print("  Wait, 2*T_a = a*(a+1) = 30, and T_5 = 15. So 2*T_5 = 30 = 2 * edges of K_6.")
print()
print("  Let me be precise: a*b = a*(a+1) = TWICE the a-th triangular number.")
print("  = the number of ORDERED pairs from {1,...,a+1} = a+1 choose 2 times 2.")
print()

# ============================================================
print()
print("  STONE 5: CAYLEY COORDINATES OF PRIMES")
print("  " + "-"*40)
print()
print("  If n is an odd prime p, what are its Cayley coordinates?")
print()
print("  p     a=(p-1)/2    b=(p+1)/2    a factors         a*b")
print("  " + "-"*60)
for p in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]:
    a, b = (p-1)//2, (p+1)//2
    af = "*".join(str(f) for f in factorize(a)) if a > 1 else "1"
    print(f"  {p:3d}    {a:4d}          {b:4d}         {af:<16s}  {a*b:6d}")

print()
print("  For prime p: a = (p-1)/2 is a VERY important number.")
print("  It is the number of QUADRATIC RESIDUES mod p.")
print("  (By Euler's criterion: the QRs mod p are exactly")
print("   {1^2, 2^2, ..., ((p-1)/2)^2} mod p.)")
print()
print("  So the Cayley a-coordinate of a prime = its QR count!")
print("  And the Paley tournament T_p is built from QRs.")
print("  The Cayley address 'a/b = QR_count / NQR_count_plus_1'.")
print()

# ============================================================
print()
print("  STONE 6: THE SUM a+b=n AS PARTITION")
print("  " + "-"*40)
print()
print("  n = a + b with b = a+1 is the UNIQUE partition of n")
print("  into two parts differing by exactly 1.")
print()
print("  This is the MOST BALANCED partition of n into two parts.")
print("  (For even n, the most balanced is n/2 + n/2.)")
print("  (For odd n, the most balanced is (n-1)/2 + (n+1)/2 = a + b.)")
print()
print("  The Cayley address IS the ratio of the most balanced partition.")
print("  It measures how close to 50/50 you can get.")
print("  For large n: a/b -> 1 (very close to balanced).")
print("  For n=3: a/b = 1/2 (the least balanced odd partition).")
print()
print("  The RAPIDITY of the address = how far from balanced:")
print("  arctanh(a/b) = arctanh((n-1)/(n+1)) = ln(n)/2 = rapidity(n).")
print("  Rapidity measures the DEPARTURE FROM BALANCE.")
print()

# ============================================================
print()
print("  STONE 7: WHAT a AND b SHARE AND DON'T SHARE")
print("  " + "-"*40)
print()
print("  gcd(a, b) = gcd(a, a+1) = 1. Always.")
print("  CONSECUTIVE INTEGERS ARE ALWAYS COPRIME.")
print()
print("  So a and b share NO prime factors.")
print("  The factorization of n = a+b is COMPLETELY SPLIT between a and b.")
print("  Every prime factor of n goes to EITHER a or b, never both.")
print()
print("  For n=15: a=7, b=8. 15=3*5. Neither 3 nor 5 divides both 7 and 8.")
print("  The factors of 15 are DISTRIBUTED: 15 = 7+8, 7 is prime, 8 = 2^3.")
print("  None of {3, 5} appear in {7} or {2^3}. They appear in the SUM, not the parts.")
print()
print("  This means: the prime factorization of n CANNOT be read from")
print("  the factorizations of a and b individually.")
print("  It emerges ONLY from their SUM.")
print("  The whole is truly more than the sum of its parts.")
print()

# ============================================================
print()
print("  STONE 8: THE CAYLEY COORDINATES OF PERFECT NUMBERS")
print("  " + "-"*40)
print()
print("  Even perfect numbers: 6, 28, 496, 8128, ...")
print("  These are even, so their Cayley addresses don't reduce (a+b=2n).")
print("  But what about their NEIGHBORS?")
print()
for pn in [6, 28, 496, 8128]:
    for n in [pn-1, pn+1]:
        if n % 2 == 1:
            a, b = (n-1)//2, (n+1)//2
            af = "*".join(str(f) for f in factorize(a)) if a > 1 else "1"
            bf = "*".join(str(f) for f in factorize(b)) if b > 1 else "1"
            print(f"  {pn} +/- 1 = {n}: a={a} ({af}), b={b} ({bf})")

print()

# ============================================================
print()
print("  STONE 9: CAYLEY COORDINATE ARITHMETIC")
print("  " + "-"*40)
print()
print("  If n has coordinates (a, a+1) and m has coordinates (c, c+1),")
print("  what are the coordinates of n*m?")
print()
print("  n*m is odd (product of odds). Its coordinates are ((nm-1)/2, (nm+1)/2).")
print("  = ((2a+1)(2c+1)-1)/2, ((2a+1)(2c+1)+1)/2)")
print("  = (4ac+2a+2c)/2, (4ac+2a+2c+2)/2)")
print("  = (2ac+a+c, 2ac+a+c+1)")
print()
print("  So the a-coordinate of n*m = 2ac + a + c = 2ac + a + c.")
print("  = a(2c+1) + c = a*m + c")
print("  = a*m + (m-1)/2")
print("  Hmm, let me verify:")
print()

for n, m in [(3,5), (3,7), (5,7), (3,3), (7,7), (3,21), (5,9)]:
    a_n, a_m = (n-1)//2, (m-1)//2
    prod = n*m
    a_prod = (prod-1)//2
    formula = 2*a_n*a_m + a_n + a_m
    print(f"  {n}*{m} = {prod}: a_prod = {a_prod}, 2*{a_n}*{a_m}+{a_n}+{a_m} = {formula}, match: {a_prod == formula}")

print()
print("  THEOREM: The Cayley a-coordinate of n*m = 2*a_n*a_m + a_n + a_m.")
print("  Or: a_{nm} = (2a_n + 1)*a_m + a_n = n*a_m + a_n.")
print()
print("  Even cleaner: a_{nm} = n*a_m + a_n.")
print("  The a-coordinate of a product = WEIGHTED SUM of coordinates!")
print()

# Verify
for n, m in [(3,5), (3,7), (5,7), (3,3), (7,7), (3,21)]:
    a_n, a_m = (n-1)//2, (m-1)//2
    a_prod = (n*m - 1)//2
    formula = n*a_m + a_n
    print(f"  a_({n}*{m}) = a_{n*m} = {a_prod}, n*a_m + a_n = {n}*{a_m}+{a_n} = {formula}, match: {a_prod == formula}")

print()
print("  PROOF: a_{nm} = (nm-1)/2 = (n(2a_m+1) - 1)/2 = (2n*a_m + n - 1)/2")
print("       = n*a_m + (n-1)/2 = n*a_m + a_n. QED.")
print()
print("  This is a SEMI-DIRECT PRODUCT structure!")
print("  Multiplication of odd numbers acts on Cayley coordinates")
print("  as an AFFINE MAP: a_{nm} = n * a_m + a_n.")
print("  The group of odd numbers under multiplication acts on")
print("  the a-coordinates by affine transformations x -> n*x + a_n.")
print()

# ============================================================
print()
print("  STONE 10: THE FORBIDDEN NUMBERS IN COORDINATE ARITHMETIC")
print("  " + "-"*40)
print()
print("  7 = 3 + 4, coordinates (3, 4).")
print("  21 = 7 * 3. In coordinates: a_21 = 7*a_3 + a_7 = 7*1 + 3 = 10.")
print("  Check: a_21 = (21-1)/2 = 10. Correct!")
print()
print("  So 21 = 7*3 and a_21 = 7*1 + 3 = 10.")
print("  The a-coordinate of 21 is 7 times the a-coordinate of 3")
print("  plus the a-coordinate of 7.")
print()
print("  In other words: 21's address (10/11) is built from")
print("  7's address (3/4) and 3's address (1/2) by the affine map.")
print()
print("  The two forbidden numbers are RELATED by the affine coordinate map.")
print("  21 = 7 * 3 in the numbers.")
print("  10 = 7 * 1 + 3 in the coordinates.")
print("  The multiplication STRUCTURE of the forbidden pair")
print("  is reflected in the AFFINE STRUCTURE of their coordinates.")
print()

# ============================================================
print()
print("  STONE 11: CAYLEY COORDINATES AND MODULAR ARITHMETIC")
print("  " + "-"*40)
print()
print("  a_n = (n-1)/2 mod small primes:")
print()
print("  n    a_n   a mod 2  a mod 3  a mod 5  a mod 7")
print("  " + "-"*50)
for n in range(3, 44, 2):
    a = (n-1)//2
    print(f"  {n:3d}   {a:3d}      {a%2}        {a%3}        {a%5}        {a%7}")

print()
print("  a_n mod 2: alternates 1,1,0,0,1,1,0,0,... (period 4, matching n mod 8)")
print("  a_n mod 3: cycles 1,2,0,1,2,0,... (period 6 in n)")
print("  a_n mod p: has period 2p in n.")
print()

# ============================================================
print()
print("  STONE 12: THE GOLDEN COORDINATES")
print("  " + "-"*40)
print()
print("  The 'Cayley coordinates' of phi and tau (not integers, but still):")
print()
print("  phi = {:.6f}. Not odd, not integer. But:".format((1+sqrt(5))/2))
print("  a_phi = (phi-1)/2 = 1/(2*phi) = {:.6f}".format(1/(2*(1+sqrt(5))/2)))
print("  b_phi = (phi+1)/2 = phi^2/2 = {:.6f}".format(((1+sqrt(5))/2)**2/2))
print("  a_phi + b_phi = phi. The golden ratio IS its own Cayley sum.")
print()
phi = (1+sqrt(5))/2
a_phi = (phi-1)/2
b_phi = (phi+1)/2
print(f"  a_phi = {a_phi:.10f} = 1/(2*phi) = phi-1)/2")
print(f"  b_phi = {b_phi:.10f} = phi^2/2")
print(f"  a_phi * b_phi = {a_phi*b_phi:.10f}")
print(f"  (phi^2-1)/4 = {(phi**2-1)/4:.10f}")
print(f"  = phi/4 = {phi/4:.10f}")
print(f"  Hmm: a_phi * b_phi = (phi^2-1)/4 = ((phi+1)(phi-1))/4 = phi/4")
print(f"  (since phi^2-1 = phi, the defining property!)")
print()
print("  THE CAYLEY PRODUCT OF PHI = PHI/4.")
print("  The golden ratio's self-similarity (phi^2 = phi+1)")
print("  means its Cayley product is just phi scaled by 1/4.")
print()

# ============================================================
print()
print("  STONE 13: THE COORDINATE MAP AS A FUNCTOR")
print("  " + "-"*40)
print()
print("  The map n -> a_n = (n-1)/2 sends:")
print("  - Odd integers to non-negative integers")
print("  - Multiplication to affine maps: a_{nm} = n*a_m + a_n")
print("  - The identity 1 to 0: a_1 = 0")
print("  - The forbidden 7 to 3: a_7 = 3 (the curvature quantum)")
print("  - The forbidden 21 to 10: a_21 = 10 (the pentagonal number 2*5)")
print()
print("  This map is a HOMOMORPHISM from (odd Z, *) to (Z, affine).")
print("  It 'linearizes' multiplicative structure into affine structure.")
print("  Compare: logarithm linearizes multiplication into addition.")
print("  This is a DISCRETE LOGARITHM-LIKE map, but affine, not linear.")
print()
print("  The affine structure means: products don't just ADD coordinates")
print("  (that would be a logarithm). They SCALE and SHIFT.")
print("  The scaling factor IS the original number n.")
print("  So multiplication by n is 'multiplication by n' in coordinates too,")
print("  but with an additive correction a_n.")
print()
print("  This correction a_n = (n-1)/2 is the 'logarithmic residue'")
print("  of the number n. It is what remains after you extract the")
print("  multiplicative part.")
print()

# ============================================================
print()
print("="*70)
print()
print("  WHAT THE STONES REVEAL")
print()
print("="*70)
print()
print("  1. The Cayley coordinate product a*b is always pronic = 2*triangular.")
print("  2. Consecutive integers a, a+1 are always coprime: the whole > parts.")
print("  3. For prime p: a = (p-1)/2 = number of quadratic residues mod p.")
print("  4. The a-coordinate encodes n mod 4: even a <=> n = 1 mod 4.")
print("  5. Multiplication of odds is AFFINE in coordinates: a_{nm} = n*a_m + a_n.")
print("  6. The forbidden numbers are affinely related: a_21 = 7*a_3 + a_7.")
print("  7. The golden ratio's Cayley product = phi/4 (from phi^2 = phi+1).")
print("  8. Only n=5 has both coordinates prime (the unique consecutive prime pair 2,3).")
print("  9. Every odd n is the most balanced partition of itself.")
print(" 10. The coordinate map is a discrete affine logarithm.")
print()
print("  And underneath it all:")
print("  n = a + b, b - a = 1, gcd(a,b) = 1, a*b = pronic.")
print("  Four constraints. They determine everything.")
