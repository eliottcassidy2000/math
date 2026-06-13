#!/usr/bin/env python3
"""
numthy_arithmetic_23.py — Number Theory & Arithmetic Geometry through (2,3)
opus-2026-03-14-S84

Exploring how KEY1=2 and KEY2=3 govern:
- Prime distribution and the prime counting function
- Quadratic reciprocity (KEY1 is THE prime)
- Zeta functions and special values
- p-adic numbers (at p = KEY1 and p = KEY2)
- Quadratic forms and class numbers
- Bernoulli numbers and their denominators
- Continued fractions and Diophantine approximation
- Partitions and congruences (Ramanujan)
- Arithmetic functions (phi, sigma, mu, tau)
- ABC conjecture implications
- Goldbach, twin primes, and additive number theory
- Dirichlet L-functions and characters
"""

import math
from fractions import Fraction
from functools import lru_cache
from collections import Counter

# Tournament vocabulary
KEY1 = 2
KEY2 = 3
KEY_SUM = KEY1 + KEY2  # 5
H_forb_1 = 7
V_PET = 10
h_E6 = 12
h_G2 = 6
BT = 24
BO = 48
BI = 120

print("=" * 70)
print("  NUMBER THEORY & ARITHMETIC GEOMETRY THROUGH (2,3)")
print("  The primes that rule all primes")
print("=" * 70)

# =====================================================================
# Part 1: PRIME COUNTING AND THE FIRST PRIMES
# =====================================================================
print("\n" + "=" * 70)
print("  Part 1: THE FIRST PRIMES ARE THE TOURNAMENT CONSTANTS")
print("=" * 70)

def sieve(n):
    """Sieve of Eratosthenes."""
    is_prime = [True] * (n + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(n**0.5) + 1):
        if is_prime[i]:
            for j in range(i*i, n+1, i):
                is_prime[j] = False
    return [p for p in range(n+1) if is_prime[p]]

primes = sieve(200)

print("\n  The first primes: ", primes[:15])
print(f"  p_1 = {primes[0]} = KEY1")
print(f"  p_2 = {primes[1]} = KEY2")
print(f"  p_3 = {primes[2]} = KEY_SUM")
print(f"  p_4 = {primes[3]} = H_forb_1")
print()

# Prime counting function values at tournament constants
pi_vals = {}
for n in [KEY1, KEY2, KEY_SUM, H_forb_1, V_PET, h_E6, BT, BO, BI]:
    count = sum(1 for p in primes if p <= n)
    pi_vals[n] = count

print("  Prime counting function pi(n) at tournament constants:")
for n, pi_n in sorted(pi_vals.items()):
    label = ""
    if pi_n == KEY1: label = " = KEY1"
    elif pi_n == KEY2: label = " = KEY2"
    elif pi_n == KEY_SUM: label = " = KEY_SUM"
    elif pi_n == H_forb_1: label = " = H_forb_1"
    elif pi_n == V_PET: label = " = V(Pet)"
    print(f"  pi({n}) = {pi_n}{label}")

print(f"""
  REMARKABLE:
  pi(KEY1) = 1
  pi(KEY2) = KEY1 = 2
  pi(KEY_SUM) = KEY2 = 3
  pi(H_forb_1) = KEY1^KEY1 = 4

  The prime counting function SHIFTS tournament constants!
  This is not a coincidence — KEY1, KEY2, KEY_SUM, H_forb_1
  are the first four primes, so pi permutes them.

  PRIME GAPS starting from KEY1:
  gap(KEY1, KEY2) = 1 (twin-adjacent)
  gap(KEY2, KEY_SUM) = KEY1 (the first gap of size KEY1!)
  gap(KEY_SUM, H_forb_1) = KEY1 (another gap of KEY1)
  gap(H_forb_1, 11) = KEY1^KEY1 = 4 (first gap > KEY1)
""")

# =====================================================================
# Part 2: QUADRATIC RECIPROCITY
# =====================================================================
print("\n" + "=" * 70)
print("  Part 2: QUADRATIC RECIPROCITY — KEY1 IS THE PRIME")
print("=" * 70)

def legendre(a, p):
    """Legendre symbol (a/p)."""
    if a % p == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1

print("  Legendre symbols (a/p) for small primes:")
print(f"  {'a':>4}", end="")
for p in [KEY1, KEY2, KEY_SUM, H_forb_1, 11, 13]:
    print(f"  (a/{p})", end="")
print()

for a in range(1, 13):
    print(f"  {a:>4}", end="")
    for p in [KEY1, KEY2, KEY_SUM, H_forb_1, 11, 13]:
        if p == KEY1:
            val = 0 if a % 2 == 0 else 1  # everything odd is QR mod 2
            print(f"  {val:>5}", end="")
        else:
            val = legendre(a, p)
            print(f"  {val:>5}", end="")
    print()

print(f"""
  QUADRATIC RECIPROCITY (Gauss):
  For odd primes p, q:
  (p/q)(q/p) = (-1)^((p-1)/KEY1 * (q-1)/KEY1)

  The exponent involves (p-1)/KEY1 — division by KEY1!

  SUPPLEMENTARY LAWS:
  (-1/p) = (-1)^((p-1)/KEY1)
  (KEY1/p) = (-1)^((p^KEY1-1)/KEY1^KEY2)

  KEY1 is special in quadratic reciprocity because:
  1. It's the only even prime (all others are odd = not divisible by KEY1)
  2. The formula involves (p-1)/KEY1 (half the group order)
  3. Z/KEY1 = F_KEY1 is the coefficient field for mod-2 invariants

  QUADRATIC RESIDUES mod small primes:
""")

for p in [KEY2, KEY_SUM, H_forb_1, 11, 13]:
    qrs = [a for a in range(1, p) if legendre(a, p) == 1]
    nqrs = [a for a in range(1, p) if legendre(a, p) == -1]
    print(f"  QR mod {p}: {qrs}")
    print(f"  NQR mod {p}: {nqrs}")
    print(f"    #QR = #NQR = (p-1)/KEY1 = {(p-1)//KEY1}")
    print()

# =====================================================================
# Part 3: BERNOULLI NUMBERS — THE (2,3) DENOMINATORS
# =====================================================================
print("\n" + "=" * 70)
print("  Part 3: BERNOULLI NUMBERS — DENOMINATORS ARE (2,3)")
print("=" * 70)

@lru_cache(maxsize=None)
def bernoulli(n):
    B = [Fraction(0)] * (n + 1)
    B[0] = Fraction(1)
    for m in range(1, n + 1):
        B[m] = Fraction(0)
        for k in range(m):
            B[m] -= Fraction(math.comb(m + 1, k), m + 1) * B[k]
    return B[n]

print("  Bernoulli numbers B_n (nonzero even ones):")
for n in range(0, 32, 2):
    b = bernoulli(n)
    if b != 0:
        num = b.numerator
        den = b.denominator
        # Factor denominator
        d = den
        factors = []
        for p in primes[:20]:
            while d % p == 0:
                factors.append(p)
                d //= p
            if d == 1:
                break
        factor_str = " * ".join(str(f) for f in factors) if factors else "1"
        label = ""
        if den == h_G2: label = " = h(G2)"
        elif den == 30: label = " = KEY1*KEY2*KEY_SUM"
        elif den == 42: label = " = KEY1*KEY2*H_forb_1"
        elif den == BT * 11 * 5: label = ""
        print(f"  B_{n} = {num}/{den} (den = {factor_str}){label}")

print(f"""
  VON STAUDT-CLAUSEN THEOREM:
  The denominator of B_(KEY1*k) is:
  denom(B_(KEY1*k)) = prod(p : (p-1) | KEY1*k)

  For B_KEY1 = 1/h(G2) = 1/6:
  Primes with (p-1) | KEY1: p = KEY1 (1|2), p = KEY2 (2|2)
  denom = KEY1 * KEY2 = h(G2) = 6. Correct!

  For B_(KEY1^KEY1) = B_4:
  (p-1) | 4: p = KEY1(1|4), KEY2(2|4), KEY_SUM(4|4)
  denom = KEY1 * KEY2 * KEY_SUM = 30. Correct!

  For B_h(G2) = B_6:
  (p-1) | 6: p = KEY1(1|6), KEY2(2|6), KEY_SUM... no (4 nmid 6), H_forb_1(6|6)
  denom = KEY1 * KEY2 * H_forb_1 = 42. Correct!

  PATTERN: B_(KEY1*k) always has KEY1 and KEY2 in the denominator!
  (Since (KEY1-1)=1 divides everything, and (KEY2-1)=KEY1 divides KEY1*k)

  SIGNIFICANCE:
  denom(B_(KEY1*k)/(KEY1*KEY1*k)) = order of im(J) in pi_(KEY1*KEY1*k-1)^s
  This connects Bernoulli numbers to stable homotopy theory!
""")

# =====================================================================
# Part 4: ZETA FUNCTION SPECIAL VALUES
# =====================================================================
print("\n" + "=" * 70)
print("  Part 4: ZETA FUNCTION — SPECIAL VALUES AT (2,3)")
print("=" * 70)

print("""
  RIEMANN ZETA FUNCTION zeta(s) = sum n^(-s):
""")

# Compute zeta at even positive integers via Bernoulli
print("  zeta(KEY1*k) for k = 1, 2, 3, ...:")
for k in range(1, 8):
    n = 2 * k
    # zeta(2k) = (-1)^(k+1) B_{2k} (2pi)^{2k} / (2 * (2k)!)
    b = bernoulli(n)
    coeff = abs(b) * Fraction(2**(n-1), math.factorial(n))
    # zeta(2k) = coeff * pi^{2k}
    zeta_val = float(coeff) * math.pi**(n)
    marker = ""
    if n == KEY1: marker = f" = pi^KEY1 / h(G2)"
    elif n == KEY1**2: marker = f" = pi^4 / 90"
    elif n == h_G2: marker = f" = pi^6 / 945"
    print(f"  zeta({n}) = {zeta_val:.10f} = |B_{n}| * (KEY1*pi)^{n} / (KEY1 * {n}!){marker}")

print(f"""
  KEY SPECIAL VALUES:
  zeta(KEY1) = pi^KEY1 / h(G2) = pi^2/6
  zeta(KEY1^KEY1) = pi^(KEY1^KEY1) / 90
  zeta(h(G2)) = pi^(h(G2)) / 945

  The denominators: 6, 90, 945, 9450, 93555, ...
  6 = h(G2) = KEY1 * KEY2
  90 = KEY1 * KEY2^KEY1 * V(Pet) = 2 * 9 * 5
  945 = KEY2^KEY2 * KEY_SUM * H_forb_1

  FUNCTIONAL EQUATION:
  zeta(s) = KEY1^s * pi^(s-1) * sin(pi*s/KEY1) * Gamma(1-s) * zeta(1-s)

  The KEY1 appears in:
  - KEY1^s factor
  - sin(pi*s/KEY1) — half-angle!
  - The symmetry point: s = 1/KEY1 (the critical line Re(s) = 1/KEY1)

  THE CRITICAL LINE is Re(s) = 1/KEY1!
  Riemann Hypothesis: all nontrivial zeros have Re(s) = 1/KEY1.
  KEY1 literally defines the most important unsolved problem in math.
""")

# =====================================================================
# Part 5: p-ADIC NUMBERS
# =====================================================================
print("\n" + "=" * 70)
print("  Part 5: p-ADIC NUMBERS AT p = KEY1 AND p = KEY2")
print("=" * 70)

print(f"""
  The p-adic integers Z_p complete Z at the prime p.

  KEY1-ADIC (2-adic):
  Z_KEY1 = lim Z/KEY1^n = {{...b_n b_(n-1) ... b_1 b_0}} (binary!)
  Units: Z_KEY1* = {{odd integers}} (index KEY1 subgroup)
  Z_KEY1* = {{+1, -1}} x (1 + KEY1^KEY1 * Z_KEY1)
  |Z_KEY1* / (Z_KEY1*)^KEY1| = KEY1^KEY2 = 8

  The KEY1-adic integers are the "natural" binary representation.
  Every KEY1-adic integer IS a binary string!

  KEY2-ADIC (3-adic):
  Z_KEY2 = lim Z/KEY2^n = {{...t_n t_(n-1) ... t_1 t_0}} (ternary!)
  Units: Z_KEY2* = {{not divisible by KEY2}}
  Z_KEY2* = {{+1, -1}} x (1 + KEY2 * Z_KEY2)
  |Z_KEY2* / (Z_KEY2*)^KEY2| = KEY2^KEY1 = 9

  HASSE-MINKOWSKI: A quadratic form over Q is isotropic iff
  it is isotropic over R AND over Q_p for ALL primes p.

  The local-global principle involves checking at:
  - The archimedean place (R)
  - p = KEY1 (the wildest prime)
  - p = KEY2 (often the next obstruction)
  - All other primes (usually easier)

  KEY1 is special among p-adic primes because:
  - KEY1 is the only RAMIFIED prime in Z[i] (Gaussian integers)
  - The KEY1-adic valuation detects PARITY
  - Hensel's lemma for p = KEY1 needs extra care (lifting mod KEY1^KEY2)
""")

# =====================================================================
# Part 6: QUADRATIC FORMS AND CLASS NUMBERS
# =====================================================================
print("\n" + "=" * 70)
print("  Part 6: QUADRATIC FORMS AND CLASS NUMBERS")
print("=" * 70)

def class_number_neg(d):
    """Compute class number h(-|d|) for negative fundamental discriminant."""
    D = -abs(d)
    if D >= 0:
        return 0
    # Kronecker symbol sum
    h = 0
    for a in range(1, abs(D)):
        # Kronecker symbol (D/a)
        kr = 0
        if math.gcd(a, abs(D)) > 1:
            kr = 0
        else:
            kr = pow(D % a if D % a >= 0 else a + (D % a), (a - 1) // 2, a) if a > 2 else 1
            if kr > a // 2:
                kr -= a
        h += kr
    # Simpler: use formula h(-d) = (1/d) sum_{a=1}^{|d|} (D/a) * a... no
    # Just count reduced forms
    h_count = 0
    b_max = int(math.sqrt(abs(D) / 3)) + 1
    for b in range(0, b_max + 1):
        if (b * b - D) % 4 != 0 and (b * b + abs(D)) % 4 != 0:
            # Try both signs
            pass
        disc_rem = abs(D) + b * b if D < 0 else b * b - D
        if disc_rem % 4 != 0:
            continue
        four_ac = disc_rem
        for a in range(max(1, b), int(math.sqrt(four_ac / 4)) + 2):
            if four_ac % (4 * a) == 0:
                c = four_ac // (4 * a)
                if c >= a and a >= b and b >= 0:
                    if b == 0 or b == a or a == c:
                        h_count += 1
                    else:
                        h_count += 2  # count both +b and -b
    return h_count

# Use known values instead
class_nums = {
    3: 1, 4: 1, 7: 1, 8: 1, 11: 1, 19: 1, 43: 1, 67: 1, 163: 1,  # h=1
    15: 2, 20: 2, 24: 2,  # h=2
    23: 3, 31: 3,  # h=3
    56: 4, 84: 4,  # h=4
    47: 5,  # h=5
}

print("  CLASS NUMBERS h(-d) for fundamental discriminants -d:")
print("  (h(-d) = 1 means unique factorization in Q(sqrt(-d)))")
print()

# Class number 1 discriminants (Heegner numbers!)
heegner = [1, 2, 3, 7, 11, 19, 43, 67, 163]
print("  h(-d) = 1 for d =", heegner)
print(f"  (these are the {len(heegner)} = KEY2^KEY1 Heegner numbers!)")
print()

print("  Tournament constants in Heegner numbers:")
for d in heegner:
    labels = []
    if d == 1: labels.append("1")
    elif d == KEY1: labels.append("KEY1")
    elif d == KEY2: labels.append("KEY2")
    elif d == H_forb_1: labels.append("H_forb_1")
    elif d == 11: labels.append("11")
    elif d == 19: labels.append("19")
    elif d == 43: labels.append("43")
    elif d == 67: labels.append("67")
    elif d == 163: labels.append("163")
    print(f"  d = {d}: h(-{d}) = 1" + (f" [{', '.join(labels)}]" if labels else ""))

print(f"""
  The Heegner numbers: 1, KEY1, KEY2, H_forb_1, 11, 19, 43, 67, 163

  First FOUR Heegner numbers = 1, KEY1, KEY2, H_forb_1
  = 1, 2, 3, 7 = the tournament constants!

  RAMANUJAN'S CONSTANT:
  e^(pi * sqrt(163)) = 262537412640768743.999999999999250...
  Almost an integer! (Off by < 10^(-12))

  j((-1+sqrt(-163))/2) = -640320^KEY2 = -(640320)^3
  640320 = KEY1^6 * KEY2 * KEY_SUM * 23 * 29
  640320^KEY2 = 262537412640768000

  The KEY2 = 3 power in j = -640320^KEY2 is fundamental:
  j-invariant involves CUBIC (KEY2) expressions.

  CLASS NUMBER FORMULA:
  h(-d) = w * sqrt(|d|) / (KEY1 * pi) * L(1, chi_d)
  where w = number of roots of unity.
  The KEY1 in the denominator is universal!
""")

# =====================================================================
# Part 7: PARTITION FUNCTION — RAMANUJAN CONGRUENCES
# =====================================================================
print("\n" + "=" * 70)
print("  Part 7: PARTITIONS — RAMANUJAN'S (2,3,5,7) CONGRUENCES")
print("=" * 70)

@lru_cache(maxsize=None)
def partition(n):
    """Number of partitions of n."""
    if n < 0: return 0
    if n == 0: return 1
    result = 0
    for k in range(1, n + 1):
        # Euler's pentagonal recurrence
        for sign in [1, -1]:
            m = k * (3 * k + sign) // 2  # generalized pentagonal
            if m > n:
                continue
            if k % 2 == 1:
                result += partition(n - m)
            else:
                result -= partition(n - m)
    return result

# Actually use proper pentagonal number theorem
@lru_cache(maxsize=None)
def p(n):
    if n < 0: return 0
    if n == 0: return 1
    s = 0
    k = 1
    while True:
        pent1 = k * (3*k - 1) // 2
        pent2 = k * (3*k + 1) // 2
        if pent1 > n and pent2 > n:
            break
        sign = (-1)**(k+1)
        if pent1 <= n:
            s += sign * p(n - pent1)
        if pent2 <= n:
            s += sign * p(n - pent2)
        k += 1
    return s

print("  First 30 partition numbers p(n):")
for i in range(30):
    print(f"  p({i:2d}) = {p(i)}", end="")
    notes = []
    if p(i) == KEY1: notes.append("KEY1")
    elif p(i) == KEY2: notes.append("KEY2")
    elif p(i) == KEY_SUM: notes.append("KEY_SUM")
    elif p(i) == H_forb_1: notes.append("H_forb_1")
    elif p(i) == BT - 2: notes.append("KEY1*11")
    elif p(i) == 42: notes.append("KEY1*KEY2*H_forb_1")
    if notes:
        print(f"  [{', '.join(notes)}]", end="")
    print()

print(f"""
  RAMANUJAN'S CONGRUENCES:
  p(KEY_SUM*n + KEY1^KEY1) = 0  (mod KEY_SUM)     [p(5n+4) = 0 mod 5]
  p(H_forb_1*n + KEY_SUM) = 0   (mod H_forb_1)   [p(7n+5) = 0 mod 7]
  p(11*n + h(G2)-h(G2)+6) = 0   (mod 11)          [p(11n+6) = 0 mod 11]

  The moduli: KEY_SUM = 5, H_forb_1 = 7, 11
  These are 3 consecutive primes starting from KEY_SUM!

  GENERATING FUNCTION:
  sum p(n) q^n = prod 1/(1-q^n)
  = 1/((1-q)(1-q^KEY1)(1-q^KEY2)(1-q^KEY1^KEY1)...)
""")

# Verify Ramanujan congruences
print("  Verification of p(5n+4) = 0 (mod 5):")
for n in range(8):
    val = p(5*n + 4)
    print(f"  p({5*n+4}) = {val}, mod 5 = {val % 5}")

print("\n  Verification of p(7n+5) = 0 (mod 7):")
for n in range(6):
    val = p(7*n + 5)
    print(f"  p({7*n+5}) = {val}, mod 7 = {val % 7}")

# =====================================================================
# Part 8: ARITHMETIC FUNCTIONS
# =====================================================================
print("\n" + "=" * 70)
print("  Part 8: ARITHMETIC FUNCTIONS — EULER, MOBIUS, SIGMA")
print("=" * 70)

def euler_phi(n):
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

def sigma(n, k=1):
    """Sum of k-th powers of divisors."""
    s = 0
    for d in range(1, n + 1):
        if n % d == 0:
            s += d**k
    return s

def mobius(n):
    if n == 1: return 1
    factors = []
    temp = n
    for p in range(2, n + 1):
        if temp % p == 0:
            count = 0
            while temp % p == 0:
                temp //= p
                count += 1
            if count > 1:
                return 0
            factors.append(p)
        if temp == 1:
            break
    return (-1)**len(factors)

print("  phi(n), sigma(n), mu(n) for tournament constants:")
print(f"  {'n':>4}  {'phi':>6}  {'sigma':>6}  {'mu':>4}  notes")
for n in [KEY1, KEY2, KEY_SUM, h_G2, H_forb_1, V_PET, h_E6, BT, BO, BI]:
    ph = euler_phi(n)
    sig = sigma(n)
    mu = mobius(n)
    notes = []
    if ph == KEY1: notes.append("phi=KEY1")
    elif ph == KEY1**2: notes.append("phi=KEY1^2")
    elif ph == KEY1**3: notes.append("phi=KEY1^3")
    if sig == KEY2: notes.append("sig=KEY2")
    elif sig == KEY1**2: notes.append("sig=KEY1^2")
    elif sig == BT: notes.append("sig=BT!")
    elif sig == h_E6: notes.append("sig=h(E6)")
    elif sig == 60: notes.append("sig=|A5|=BI/KEY1")
    print(f"  {n:>4}  {ph:>6}  {sig:>6}  {mu:>4}  {', '.join(notes)}")

print(f"""
  REMARKABLE VALUES:
  phi(KEY1) = 1, phi(KEY2) = KEY1, phi(KEY_SUM) = KEY1^KEY1,
  phi(H_forb_1) = h(G2)
  phi(BT) = KEY1^KEY2 = 8

  sigma(KEY1) = KEY2 (!), sigma(KEY2) = KEY1^KEY1,
  sigma(KEY_SUM) = h(G2)
  sigma(V(Pet)) = 18 = KEY1 * KEY2^KEY1

  PERFECT NUMBERS:
  sigma(n) = KEY1 * n iff n is perfect.
  P_1 = h(G2) = 6 = KEY1 * KEY2
  P_2 = 28 = KEY1^KEY1 * H_forb_1
  P_3 = 496 = KEY1^KEY1^KEY1 * 31
  P_4 = 8128 = KEY1^h(G2) * 127

  Perfect numbers = KEY1^(p-1) * (KEY1^p - 1) when KEY1^p - 1 is prime.
  They are built ENTIRELY from KEY1!

  sigma(KEY1) = KEY2 is the simplest instance:
  KEY1 + 1 = KEY2, so sigma(KEY1) = 1 + KEY1 = KEY2.
""")

# =====================================================================
# Part 9: CONTINUED FRACTIONS — (2,3) APPROXIMATION
# =====================================================================
print("\n" + "=" * 70)
print("  Part 9: CONTINUED FRACTIONS")
print("=" * 70)

def cf_expansion(x, n_terms=15):
    """Compute continued fraction expansion of x."""
    result = []
    for _ in range(n_terms):
        a = int(math.floor(x))
        result.append(a)
        frac = x - a
        if abs(frac) < 1e-12:
            break
        x = 1.0 / frac
    return result

# Key continued fractions
targets = [
    ("sqrt(KEY1) = sqrt(2)", math.sqrt(2)),
    ("sqrt(KEY2) = sqrt(3)", math.sqrt(3)),
    ("sqrt(KEY_SUM) = sqrt(5)", math.sqrt(5)),
    ("phi = (1+sqrt(5))/2", (1 + math.sqrt(5)) / 2),
    ("e", math.e),
    ("pi", math.pi),
    ("log_2(3)", math.log(3) / math.log(2)),
    ("log_3(2)", math.log(2) / math.log(3)),
]

for name, val in targets:
    cf = cf_expansion(val)
    print(f"  {name} = {val:.10f}")
    print(f"    CF = [{cf[0]}; {', '.join(str(c) for c in cf[1:])}]")

print(f"""
  GOLDEN RATIO phi has CF = [1; 1, 1, 1, ...] — all 1s!
  This makes phi the MOST IRRATIONAL number.

  sqrt(KEY1): CF = [1; KEY1, KEY1, KEY1, ...] — all KEY1!
  sqrt(KEY2): CF = [1; 1, KEY1, 1, KEY1, ...] — period (1, KEY1)!

  log_KEY1(KEY2) has CF = [1; 1, 1, KEY1, KEY1, KEY2, ...]
  This encodes the relationship between KEY1 and KEY2
  in terms of continued fraction convergents.

  BEST RATIONAL APPROXIMATIONS to log_2(3):
  The convergents give the musical temperaments!
""")

# Compute convergents
cf_log23 = cf_expansion(math.log(3) / math.log(2), 12)
h_prev, h_curr = 0, 1
k_prev, k_curr = 1, 0
print("  Convergents of log_2(3):")
for i, a in enumerate(cf_log23):
    h_new = a * h_curr + h_prev
    k_new = a * k_curr + k_prev
    error = abs(math.log(3) / math.log(2) - h_new / k_new)
    label = ""
    if k_new == KEY_SUM: label = " <- pentatonic (KEY_SUM notes)"
    elif k_new == H_forb_1: label = " <- diatonic (H_forb_1 notes)"
    elif k_new == h_E6: label = " <- chromatic (h(E6) notes)"
    elif k_new == 53: label = " <- 53-TET (excellent)"
    print(f"  {h_new}/{k_new} = {h_new/k_new:.8f} (error {error:.2e}){label}")
    h_prev, h_curr = h_curr, h_new
    k_prev, k_curr = k_curr, k_new

# =====================================================================
# Part 10: DIRICHLET L-FUNCTIONS AND CHARACTERS
# =====================================================================
print("\n" + "=" * 70)
print("  Part 10: DIRICHLET L-FUNCTIONS — CHARACTERS MOD KEY1, KEY2")
print("=" * 70)

print(f"""
  Dirichlet characters mod q:
  chi: (Z/qZ)* -> C*

  CHARACTERS MOD KEY1 = 2:
  Only the trivial character chi_0.
  (Z/KEY1Z)* = {{1}} has order 1.

  CHARACTERS MOD KEY2 = 3:
  (Z/KEY2Z)* = {{1, 2}} ~ Z/KEY1
  chi_0: trivial (1, 1)
  chi_1: nontrivial (1, -1)

  chi_1 is the LEGENDRE SYMBOL (-KEY2/.)
  L(1, chi_1) = pi / (KEY2 * sqrt(KEY2)) = pi / (3*sqrt(3))

  CHARACTERS MOD KEY1^KEY1 = 4:
  (Z/4Z)* = {{1, 3}} ~ Z/KEY1
  chi: the Dirichlet character (-4/.) = (-1)^((n-1)/2)
  L(1, chi) = pi/KEY1^KEY1 = pi/4 (Leibniz formula!)

  CHARACTERS MOD KEY1^KEY2 = 8:
  (Z/8Z)* = {{1, 3, 5, 7}} ~ Z/KEY1 x Z/KEY1
  KEY1 independent characters (plus trivial).

  PRIME NUMBER THEOREM IN ARITHMETIC PROGRESSIONS:
  pi(x; q, a) ~ x / (phi(q) * ln(x))

  For q = KEY1: pi(x; KEY1, 1) ~ x / (1 * ln(x))
  — ALL primes > KEY1 are odd!

  For q = KEY2: pi(x; KEY2, 1) ~ pi(x; KEY2, 2) ~ x / (KEY1 * ln(x))
  — primes split equally between 1 mod KEY2 and KEY1 mod KEY2.

  CHEBYSHEV BIAS:
  Primes tend to be KEY1 mod KEY2 (quadratic nonresidues)
  more often than 1 mod KEY2 (quadratic residues).
  This is the "Chebyshev bias" — controlled by zeros of L-functions!
""")

# Count primes in residue classes
primes_1mod3 = [p for p in primes if p > 3 and p % 3 == 1]
primes_2mod3 = [p for p in primes if p > 3 and p % 3 == 2]
print(f"  Primes up to 200 that are 1 mod KEY2: {len(primes_1mod3)}")
print(f"  Primes up to 200 that are KEY1 mod KEY2: {len(primes_2mod3)}")
print(f"  Bias: {len(primes_2mod3) - len(primes_1mod3)} more in KEY1 mod KEY2 class")

# =====================================================================
# Part 11: ABC CONJECTURE AND ADDITIVE NUMBER THEORY
# =====================================================================
print("\n" + "=" * 70)
print("  Part 11: ABC CONJECTURE AND ADDITIVE PROBLEMS")
print("=" * 70)

def rad(n):
    """Radical of n: product of distinct prime factors."""
    r = 1
    for p in range(2, n + 1):
        if n % p == 0:
            r *= p
            while n % p == 0:
                n //= p
        if n == 1:
            break
    return r

print(f"""
  ABC CONJECTURE (Masser-Oesterle, 1985):
  For a + b = c with gcd(a,b) = 1:
  c < K(eps) * rad(abc)^(1+eps) for any eps > 0.

  TOURNAMENT INSTANCES:
  KEY1 + KEY2 = KEY_SUM: rad(KEY1*KEY2*KEY_SUM) = KEY1*KEY2*KEY_SUM = 30
  KEY1 + KEY_SUM = H_forb_1: rad(KEY1*KEY_SUM*H_forb_1) = KEY1*KEY_SUM*H_forb_1 = 70
  KEY2 + KEY1^KEY1 = H_forb_1: rad(KEY2*KEY1^KEY1*H_forb_1) = KEY1*KEY2*H_forb_1 = 42
""")

# Best ABC triples
print("  High-quality ABC triples (q = log(c)/log(rad(abc))):")
abc_triples = [
    (1, 2, 3),
    (1, 8, 9),    # 1 + 2^3 = 3^2
    (1, 63, 64),  # 1 + 3^2*7 = 2^6
    (5, 27, 32),  # 5 + 3^3 = 2^5
    (2, 6561, 6563), # not, let me use real ones
    (1, 2, 3),
    (1, 8, 9),
    (2, 3**10, 3**10 + 2), # probably not great
]

# Just show the famous ones
for a, b in [(1, 2), (1, 8), (5, 27), (1, 63), (1, 80), (32, 49)]:
    c = a + b
    r = rad(a * b * c)
    if r > 0:
        q = math.log(c) / math.log(r) if r > 1 else 0
        labels = []
        if a in [KEY1, KEY2, KEY_SUM, H_forb_1]: labels.append(f"a=tc")
        if b in [KEY1, KEY2, KEY_SUM, H_forb_1, 8, 27, 64, 81]: labels.append(f"b=tc")
        if c in [KEY1, KEY2, KEY_SUM, H_forb_1]: labels.append(f"c=tc")
        print(f"  {a} + {b} = {c}: rad = {r}, q = {q:.4f} {labels}")

print(f"""
  GOLDBACH CONJECTURE:
  Every even number > KEY1 is a sum of KEY1 primes.
  KEY1 is BOTH the base case AND the parity constraint!

  KEY1^KEY1 = 4 = KEY1 + KEY1
  h(G2) = 6 = KEY2 + KEY2
  KEY1^KEY2 = 8 = KEY2 + KEY_SUM = KEY2 + KEY_SUM
  V(Pet) = 10 = KEY2 + H_forb_1 = KEY_SUM + KEY_SUM
  h(E6) = 12 = KEY_SUM + H_forb_1
  BT = 24 = KEY_SUM + 19 = H_forb_1 + 17 = 11 + 13

  WARING'S PROBLEM:
  g(k) = smallest s such that every n is a sum of s k-th powers.
  g(KEY1) = KEY1^KEY1 = 4 (Lagrange: 4 squares suffice!)
  g(KEY2) = 9 (KEY2^KEY1 = 9 cubes suffice)
  g(KEY1^KEY1) = 19 (19 fourth powers)

  g(KEY1) = KEY1^KEY1 = 4: EVERY integer is a sum of KEY1^KEY1 squares!
  g(KEY2) = KEY2^KEY1 = 9: EVERY integer is a sum of KEY2^KEY1 cubes!
""")

# =====================================================================
# Part 12: GRAND SYNTHESIS — NUMBER THEORY IS (2,3)
# =====================================================================
print("\n" + "=" * 70)
print("  Part 12: GRAND SYNTHESIS — NUMBER THEORY IS (2,3)")
print("=" * 70)

print("""
======================================================================
  NUMBER THEORY IS BUILT FROM KEY1 = 2 AND KEY2 = 3
======================================================================

1. PRIMES:
   p_1 = KEY1, p_2 = KEY2, p_3 = KEY_SUM, p_4 = H_forb_1
   The first four primes ARE the tournament constants!
   pi(KEY1)=1, pi(KEY2)=KEY1, pi(KEY_SUM)=KEY2, pi(H_forb_1)=KEY1^KEY1

2. QUADRATIC RECIPROCITY:
   (p/q)(q/p) = (-1)^((p-1)/KEY1 * (q-1)/KEY1)
   KEY1 divides every formula in number theory via (p-1)/KEY1
   KEY1 is the unique even prime

3. BERNOULLI NUMBERS:
   denom(B_(KEY1*k)) always contains KEY1 * KEY2
   These denominators = orders of image of J
   Von Staudt-Clausen: denom = prod(p : (p-1) | KEY1*k)

4. ZETA FUNCTION:
   zeta(KEY1) = pi^KEY1 / (KEY1 * KEY2) = pi^2/6
   Critical line: Re(s) = 1/KEY1 (Riemann Hypothesis!)
   Functional equation involves KEY1^s

5. p-ADIC NUMBERS:
   Z_KEY1 = binary (KEY1-ary) representation of integers
   Z_KEY2 = ternary (KEY2-ary) representation
   KEY1-adic is wildest (only ramified prime in Z[i])

6. CLASS NUMBERS:
   Heegner numbers: 1, KEY1, KEY2, H_forb_1, 11, 19, 43, 67, 163
   First four = tournament constants!
   j-invariant involves KEY2 powers

7. PARTITIONS:
   Ramanujan congruences mod KEY_SUM, H_forb_1, 11
   = three consecutive primes starting from KEY_SUM = KEY1 + KEY2

8. ARITHMETIC FUNCTIONS:
   sigma(KEY1) = KEY2 (simplest nontrivial!)
   phi(KEY2) = KEY1 (simplest nontrivial!)
   Perfect numbers = KEY1^(p-1) * (KEY1^p - 1)

9. CONTINUED FRACTIONS:
   sqrt(KEY1): CF = [1; KEY1, KEY1, KEY1, ...]
   log_KEY1(KEY2) encodes all of music theory

10. WARING:
    g(KEY1) = KEY1^KEY1 = 4 squares
    g(KEY2) = KEY2^KEY1 = 9 cubes

THE CROWN JEWEL:
   The Riemann Hypothesis lives on Re(s) = 1/KEY1.
   Bernoulli denominators are products of primes p with (p-1)|KEY1k.
   Perfect numbers are built from KEY1 alone.
   The first four primes are the tournament constants.

   NUMBER THEORY IS (2,3).
""")
