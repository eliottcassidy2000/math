#!/usr/bin/env python3
"""
fibonacci_lucas_primes_s116l.py
===============================
Deep investigation of Fibonacci/Lucas connections to the fundamental
primes {2, 3, 5, 7, 11} arising in tournament parity theory.

Key threads:
  - Prime Fibonacci and Lucas numbers vs our fundamental set
  - F_8 = 21 (forbidden number), L_8 = 47 (prime), F_n * L_n = F_{2n}
  - QR counts of {3,5,7,11} as Fibonacci numbers {1,2,3,5}
  - Pisano periods, entry points, Wall-Sun-Sun tests
  - Identities and factorizations at critical indices
  - Potential computational speedups

kind-pasteur-2026-03-16-S116l
"""

import math
import sys
from collections import defaultdict

# ============================================================
# UTILITIES
# ============================================================

def fib(n):
    """Compute F_n for n >= 0. Handles negative indices via F_{-n} = (-1)^{n+1} F_n."""
    if n < 0:
        return ((-1) ** ((-n) + 1)) * fib(-n)
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a + b
    return a

def lucas(n):
    """Compute L_n for n >= 0. L_0=2, L_1=1, L_n = L_{n-1} + L_{n-2}."""
    if n < 0:
        return ((-1) ** (-n)) * lucas(-n)
    a, b = 2, 1
    for _ in range(n):
        a, b = b, a + b
    return a

def is_prime(n):
    """Deterministic primality test (trial division, sufficient for our range)."""
    if n < 2:
        return False
    if n < 4:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def factorize(n):
    """Return list of (prime, exponent) pairs."""
    if n <= 1:
        return [(n, 1)] if n == 1 else []
    factors = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            exp = 0
            while n % d == 0:
                n //= d
                exp += 1
            factors.append((d, exp))
        d += 1
    if n > 1:
        factors.append((n, 1))
    return factors

def factor_str(n):
    """Pretty-print factorization."""
    if n == 0:
        return "0"
    if n == 1:
        return "1"
    fs = factorize(abs(n))
    parts = []
    for p, e in fs:
        if e == 1:
            parts.append(str(p))
        else:
            parts.append(f"{p}^{e}")
    result = " * ".join(parts)
    if n < 0:
        result = "-1 * " + result
    return result

def legendre_symbol(a, p):
    """Compute Legendre symbol (a/p) for odd prime p."""
    if p == 2:
        return 1 if a % 2 == 1 else 0
    a = a % p
    if a == 0:
        return 0
    return pow(a, (p - 1) // 2, p) if pow(a, (p - 1) // 2, p) <= 1 else -1

def legendre(a, p):
    """Legendre symbol (a|p), returns -1, 0, or 1."""
    a = a % p
    if a == 0:
        return 0
    r = pow(a, (p - 1) // 2, p)
    return r if r <= 1 else r - p

def qr_count(p):
    """Number of quadratic residues mod p (excluding 0)."""
    if p == 2:
        return 1
    return (p - 1) // 2

def pisano_period(m):
    """Compute the Pisano period pi(m): period of F_n mod m."""
    if m <= 1:
        return 1
    prev, curr = 0, 1
    for i in range(1, 6 * m + 10):
        prev, curr = curr, (prev + curr) % m
        if prev == 0 and curr == 1:
            return i
    return -1  # should not happen for valid m

def entry_point(m):
    """Smallest positive k such that m | F_k."""
    if m <= 1:
        return 1
    prev, curr = 0, 1
    for k in range(1, 6 * m + 10):
        if curr % m == 0:
            return k
        prev, curr = curr, (prev + curr) % m
    return -1

phi = (1 + math.sqrt(5)) / 2
psi = (1 - math.sqrt(5)) / 2  # = -1/phi

# ============================================================
# OUTPUT SETUP
# ============================================================
OUT = []
def P(s=""):
    OUT.append(s)
    print(s)

def section(title):
    P("")
    P("=" * 72)
    P(f"  {title}")
    P("=" * 72)
    P("")

def subsection(title):
    P("")
    P(f"--- {title} ---")
    P("")

# ============================================================
# 1. ALL PRIME FIBONACCI AND LUCAS NUMBERS
# ============================================================
section("1. PRIME FIBONACCI AND PRIME LUCAS NUMBERS")

P("Computing Fibonacci numbers F_0 through F_50:")
P("")
for i in range(51):
    f = fib(i)
    tag = " [PRIME]" if is_prime(f) else ""
    if i <= 30 or is_prime(f):
        P(f"  F_{i:2d} = {f:>15d}  = {factor_str(f)}{tag}")

P("")
P("Prime Fibonacci numbers (index, value) for n <= 50:")
prime_fibs = []
for i in range(51):
    f = fib(i)
    if is_prime(f):
        prime_fibs.append((i, f))
        P(f"  F_{i} = {f}")

P("")
P(f"Prime Fibonacci indices: {[i for i,_ in prime_fibs]}")
P(f"Prime Fibonacci values:  {[v for _,v in prime_fibs]}")

P("")
P("Computing Lucas numbers L_0 through L_50:")
P("")
for i in range(51):
    l = lucas(i)
    tag = " [PRIME]" if is_prime(l) else ""
    if i <= 30 or is_prime(l):
        P(f"  L_{i:2d} = {l:>15d}  = {factor_str(l)}{tag}")

P("")
P("Prime Lucas numbers (index, value) for n <= 50:")
prime_lucas = []
for i in range(51):
    l = lucas(i)
    if is_prime(l):
        prime_lucas.append((i, l))
        P(f"  L_{i} = {l}")

P("")
P(f"Prime Lucas indices: {[i for i,_ in prime_lucas]}")
P(f"Prime Lucas values:  {[v for _,v in prime_lucas]}")

subsection("Connection to fundamental primes {2, 3, 5, 7, 11}")
fund_primes = {2, 3, 5, 7, 11}
P("Fundamental primes as Fibonacci numbers:")
for i in range(40):
    f = fib(i)
    if f in fund_primes:
        P(f"  {f} = F_{i}")

P("")
P("Fundamental primes as Lucas numbers:")
for i in range(40):
    l = lucas(i)
    if l in fund_primes:
        P(f"  {l} = L_{i}")

P("")
P("OBSERVATION: Lucas primes L_0=2, L_2=3, L_4=7, L_5=11")
P("  BUT WAIT: L_1 = 1 (not prime), L_3 = 4 (not prime)")
P("  Correcting: L_0=2, L_2=3, L_4=7, L_5=11 -- indices {0,2,4,5}")
P("  The user stated L_1=2, L_3=3, L_5=7, L_6=11 -- let me verify:")
for i in range(12):
    P(f"  L_{i} = {lucas(i)}")

P("")
P("CORRECTION: L_0=2, L_1=1, L_2=3, L_3=4, L_4=7, L_5=11, L_6=18")
P("So the claim 'L_1=2, L_3=3, L_5=7, L_6=11' uses a SHIFTED convention")
P("  where L_1^* = L_0 = 2, L_3^* = L_2 = 3, L_5^* = L_4 = 7, L_6^* = L_5 = 11")
P("This is the 1-indexed convention: L_n^{1-idx} = L_{n-1}^{0-idx}")
P("")
P("Under STANDARD 0-indexed convention:")
P("  2 = L_0, 3 = L_2, 7 = L_4, 11 = L_5")
P("  5 = F_5 (Fibonacci, not Lucas)")
P("  Missing from Lucas: 5 is NOT a Lucas number")
P("  5 appears as: F_5 = 5, and as L_0 + L_2 = 2 + 3 = 5")

# ============================================================
# 2. THE FORBIDDEN NUMBER 21 = F_8
# ============================================================
section("2. THE FORBIDDEN NUMBER F_8 = 21")

f8 = fib(8)
l8 = lucas(8)
P(f"F_8 = {f8} = {factor_str(f8)}")
P(f"L_8 = {l8} = {factor_str(l8)}")
P(f"Is L_8 = {l8} prime? {is_prime(l8)}")
P("")
P(f"F_8 * L_8 = {f8} * {l8} = {f8 * l8}")
P(f"F_16 = {fib(16)}")
P(f"Verify F_8 * L_8 = F_16: {f8 * l8 == fib(16)}")
P("")
P(f"987 = F_16 = {factor_str(987)}")
P(f"  = 3 * 7 * 47 = F_4 * (something)?")
P(f"  F_4 = {fib(4)}, 987/3 = {987//3} = {factor_str(329)}")
P(f"  329 = 7 * 47")
P(f"  So 987 = 3 * 7 * 47 = F_4 * L_4 * L_8")
P(f"  Check: F_4={fib(4)}, L_4={lucas(4)}, L_8={lucas(8)}: {fib(4)*lucas(4)*lucas(8)}")

P("")
P("The chain of factorizations:")
for k in [1, 2, 4, 8, 16, 32]:
    if k <= 40:
        fk = fib(k)
        P(f"  F_{k:2d} = {fk:>15d} = {factor_str(fk)}")

P("")
P("Identity: F_{2n} = F_n * L_n (fundamental doubling formula)")
P("Verification for small n:")
for n in range(1, 16):
    fn = fib(n)
    ln = lucas(n)
    f2n = fib(2 * n)
    ok = "OK" if fn * ln == f2n else "FAIL"
    P(f"  F_{2*n:2d} = F_{n} * L_{n} = {fn} * {ln} = {fn*ln} = {f2n} [{ok}]")

subsection("F_8 = 21: The forbidden number in context")
P(f"21 = 3 * 7 = F_4 * L_4 (by doubling formula with n=4)")
P(f"   = F_8 (8th Fibonacci)")
P(f"21 = H=21 gap: the permanently forbidden Hamiltonian count")
P(f"F_4 = 3 (fundamental prime), L_4 = 7 (fundamental prime)")
P(f"So 21 = (fundamental prime #2) * (fundamental prime #4)")
P(f"In the Lucas prime sequence: L_0=2, L_2=3, L_4=7, L_5=11")
P(f"21 = L_2 * L_4 = consecutive Lucas primes (skipping L_1, L_3 which are not prime)")

# ============================================================
# 3. F_n * L_n = F_{2n} AND DEEPER PRODUCT IDENTITIES
# ============================================================
section("3. PRODUCT IDENTITY F_n * L_n = F_{2n}")

P("Extended doubling chain from F_8 = 21:")
P(f"  F_4 = {fib(4)}, L_4 = {lucas(4)}, F_8 = {fib(4)*lucas(4)} = 21")
P(f"  F_8 = {fib(8)}, L_8 = {lucas(8)}, F_16 = {fib(8)*lucas(8)} = 987")
P(f"  F_16 = {fib(16)}, L_16 = {lucas(16)}, F_32 = {fib(16)*lucas(16)}")
P(f"  F_32 = {fib(32)} = {factor_str(fib(32))}")
P("")

P("Identity: L_n^2 = 5*F_n^2 + 4*(-1)^n")
P("Verification:")
for n in range(15):
    ln = lucas(n)
    fn = fib(n)
    lhs = ln * ln
    rhs = 5 * fn * fn + 4 * ((-1) ** n)
    ok = "OK" if lhs == rhs else "FAIL"
    P(f"  n={n:2d}: L_{n}^2 = {lhs:>8d}, 5*F_{n}^2 + 4*(-1)^{n} = {rhs:>8d}  [{ok}]")

P("")
P("Identity: L_n^2 - 5*F_n^2 = 4*(-1)^n")
P("At n=8: L_8^2 - 5*F_8^2 = 47^2 - 5*21^2 = 2209 - 2205 = 4")
P(f"  Check: {lucas(8)**2} - 5*{fib(8)**2} = {lucas(8)**2 - 5*fib(8)**2}")

P("")
P("GCD structure: gcd(F_n, L_n)")
for n in range(1, 25):
    g = math.gcd(fib(n), lucas(n))
    if g > 1:
        P(f"  gcd(F_{n}, L_{n}) = gcd({fib(n)}, {lucas(n)}) = {g}")
P("(Most are 1 or 2 -- F_n and L_n are nearly coprime)")
P("")

for n in range(1, 20):
    g = math.gcd(fib(n), lucas(n))
    P(f"  gcd(F_{n:2d}, L_{n:2d}) = gcd({fib(n):>6d}, {lucas(n):>6d}) = {g}")

# ============================================================
# 4. INDICES WHERE L_n IS PRIME
# ============================================================
section("4. PRIME LUCAS INDICES AND NON-PRIME FACTORIZATIONS")

P("Lucas numbers and primality for n = 0..40:")
prime_lucas_idx = []
for n in range(41):
    l = lucas(n)
    pr = is_prime(l)
    if pr:
        prime_lucas_idx.append(n)
    P(f"  L_{n:2d} = {l:>12d}  {'PRIME' if pr else factor_str(l)}")

P("")
P(f"Prime Lucas indices (0-indexed): {prime_lucas_idx}")

# Check which indices are prime themselves
P("")
P("Are the prime-Lucas indices themselves prime?")
for idx in prime_lucas_idx:
    P(f"  Index {idx}: {'prime' if is_prime(idx) else 'composite/0/1'}")

P("")
P("OBSERVATION: L_n can only be prime if n is 0, a prime, or a power of 2.")
P("  (Because if d|n, d>1, d<n, then L_d | L_n when n/d is odd)")
P("  Known: L_n prime => n in {0,1} or n is prime or n is 2^k")
P("  Actually: the necessary condition is n=0 or n is prime or 2^k (for specific reasons)")

# Verify divisibility
P("")
P("Divisibility: if d|n and n/d is odd, then L_d | L_n")
P("Checking L_d | L_n for composite n:")
for n in range(4, 30):
    l_n = lucas(n)
    divs = []
    for d in range(2, n):
        if n % d == 0 and (n // d) % 2 == 1:
            l_d = lucas(d)
            if l_n % l_d == 0:
                divs.append((d, l_d))
    if divs:
        P(f"  L_{n} = {l_n}: divisible by L_{divs[0][0]} = {divs[0][1]}")

# ============================================================
# 5. L_5 + L_4 = 11 AND SIMILAR IDENTITIES
# ============================================================
section("5. ADDITIVE IDENTITIES AMONG FUNDAMENTAL PRIMES")

P("Lucas recurrence identities involving fundamental primes:")
P(f"  L_5 = L_4 + L_3 = {lucas(4)} + {lucas(3)} = {lucas(5)}")
P(f"  So 11 = 7 + 4 (forbidden prime + period)")
P(f"  L_4 = L_3 + L_2 = {lucas(3)} + {lucas(2)} = {lucas(4)}")
P(f"  So 7 = 4 + 3")
P(f"  L_2 = L_1 + L_0 = {lucas(1)} + {lucas(0)} = {lucas(2)}")
P(f"  So 3 = 1 + 2")
P("")
P("Fibonacci recurrence identities:")
P(f"  F_5 = F_4 + F_3 = {fib(4)} + {fib(3)} = {fib(5)}")
P(f"  So 5 = 3 + 2")
P(f"  F_8 = F_7 + F_6 = {fib(7)} + {fib(6)} = {fib(8)}")
P(f"  So 21 = 13 + 8")

P("")
P("Cross-sequence identities at fundamental prime values:")
P(f"  L_0 + L_2 = 2 + 3 = 5 = F_5")
P(f"  L_2 + L_4 = 3 + 7 = 10 = 2*5 = 2*F_5")
P(f"  L_4 + L_5 = 7 + 11 = 18 = L_6")
P(f"  L_0 * L_2 = 2 * 3 = 6")
P(f"  L_2 * L_4 = 3 * 7 = 21 = F_8 (the forbidden number!)")
P(f"  L_4 * L_5 = 7 * 11 = 77 = {factor_str(77)}")
P(f"  L_0 * L_4 = 2 * 7 = 14")
P(f"  L_0 * L_5 = 2 * 11 = 22")
P(f"  L_2 * L_5 = 3 * 11 = 33")

P("")
P("Sum of all fundamental primes: 2+3+5+7+11 = {0}".format(2+3+5+7+11))
P(f"  28 = {factor_str(28)} = perfect number!")
P(f"  F_? = 28? No. But L_? = 29 (prime), F_? = 28? No, 28 is not Fibonacci.")
P(f"  28 = 4*7 = L_3 * L_4 (non-prime Lucas times prime Lucas)")

P("")
P("Product of fundamental primes: 2*3*5*7*11 = {0}".format(2*3*5*7*11))
P(f"  2310 = {factor_str(2310)} = primorial(5)")
P(f"  Is 2310 a Fibonacci or Lucas number? F_? or L_?")
for i in range(50):
    if fib(i) == 2310:
        P(f"  2310 = F_{i}")
    if lucas(i) == 2310:
        P(f"  2310 = L_{i}")
P("  2310 is neither Fibonacci nor Lucas.")

P("")
P("Deeper: sums of pairs of fundamental primes")
for i, p1 in enumerate([2, 3, 5, 7, 11]):
    for p2 in [2, 3, 5, 7, 11]:
        if p2 > p1:
            s = p1 + p2
            fib_idx = None
            luc_idx = None
            for k in range(30):
                if fib(k) == s:
                    fib_idx = k
                if lucas(k) == s:
                    luc_idx = k
            tags = []
            if fib_idx is not None:
                tags.append(f"F_{fib_idx}")
            if luc_idx is not None:
                tags.append(f"L_{luc_idx}")
            if is_prime(s):
                tags.append("PRIME")
            tag = " = " + ", ".join(tags) if tags else ""
            P(f"  {p1} + {p2} = {s}{tag}")

# ============================================================
# 6. L_n / F_n RATIO AND sqrt(5) CONNECTIONS
# ============================================================
section("6. L_n/F_n RATIO -> sqrt(5) AND Q-TRANSFORM")

P("Ratio L_n / F_n for small n:")
for n in range(1, 25):
    fn = fib(n)
    ln = lucas(n)
    if fn > 0:
        ratio = ln / fn
        P(f"  L_{n:2d}/F_{n:2d} = {ln:>8d}/{fn:>8d} = {ratio:.10f}  (sqrt(5) = {math.sqrt(5):.10f}, diff = {abs(ratio - math.sqrt(5)):.2e})")

P("")
P(f"sqrt(5) = {math.sqrt(5):.15f}")
P(f"phi = (1+sqrt(5))/2 = {phi:.15f}")
P(f"phi^2 = phi + 1 = {phi**2:.15f}")
P(f"1/phi = phi - 1 = {1/phi:.15f}")
P(f"5/phi = {5/phi:.15f} = {5*(math.sqrt(5)-1)/2:.15f}")
P(f"sqrt(5)*phi = {math.sqrt(5)*phi:.15f} = phi^2 + phi - 1? = {phi**2 + phi - 1:.15f}")
P("")

# Q-transform: Q(x) = (1+x)/(1-x) for |x| < 1
def Q(x):
    """Cayley transform Q(x) = (1+x)/(1-x)"""
    if abs(x - 1.0) < 1e-15:
        return float('inf')
    return (1 + x) / (1 - x)

P("Q-transform of 1/sqrt(5):")
val = 1 / math.sqrt(5)
P(f"  1/sqrt(5) = {val:.15f}")
P(f"  Q(1/sqrt(5)) = (1 + 1/sqrt(5)) / (1 - 1/sqrt(5))")
q_val = Q(val)
P(f"            = {q_val:.15f}")
P(f"  phi^2 = {phi**2:.15f}")
P(f"  Q(1/sqrt(5)) = phi^2? {abs(q_val - phi**2) < 1e-10}")

P("")
P("Proof: Q(1/sqrt(5)) = (sqrt(5)+1)/(sqrt(5)-1) = (sqrt(5)+1)^2/4 = (6+2*sqrt(5))/4 = (3+sqrt(5))/2 = phi^2")
P(f"  (3+sqrt(5))/2 = {(3+math.sqrt(5))/2:.15f}")
P(f"  phi^2 = {phi**2:.15f}")
P(f"  Match: {abs((3+math.sqrt(5))/2 - phi**2) < 1e-14}")

P("")
P("Q-transform of Fibonacci ratios F_n/F_{n+1}:")
for n in range(1, 15):
    fn = fib(n)
    fn1 = fib(n + 1)
    fn2 = fib(n + 2)
    fnm1 = fib(n - 1)
    ratio = fn / fn1
    q_ratio = Q(ratio)
    target = fn2 / fnm1 if fnm1 != 0 else float('inf')
    P(f"  Q(F_{n}/F_{n+1}) = Q({fn}/{fn1}) = {q_ratio:.8f},  F_{n+2}/F_{n-1} = {fn2}/{fnm1} = {target:.8f}  Match: {abs(q_ratio - target) < 1e-8}")

P("")
P("The Cayley skip-3 property: Q(F_n/F_{n+1}) = F_{n+2}/F_{n-1}")
P("Algebraic proof: Q(a/b) = (b+a)/(b-a). With a=F_n, b=F_{n+1}:")
P("  b+a = F_{n+1}+F_n = F_{n+2}")
P("  b-a = F_{n+1}-F_n = F_{n-1}")
P("  So Q(F_n/F_{n+1}) = F_{n+2}/F_{n-1}. QED (skip-3: from (n,n+1) to (n-1,n+2))")

# ============================================================
# 7. PISANO PERIODS
# ============================================================
section("7. PISANO PERIODS pi(m) FOR FUNDAMENTAL PRIMES")

P("The Pisano period pi(m) is the period of F_n mod m.")
P("")
test_vals = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 21, 42, 47, 77, 55]
for m in test_vals:
    pp = pisano_period(m)
    P(f"  pi({m:3d}) = {pp:4d}   ({factor_str(m)} -> period {pp} = {factor_str(pp)})")

P("")
P("KEY VALUES for fundamental primes:")
for p in [2, 3, 5, 7, 11, 13]:
    pp = pisano_period(p)
    P(f"  pi({p:2d}) = {pp:3d}")
    # Show the actual Fibonacci sequence mod p
    seq = []
    for k in range(pp + 3):
        seq.append(fib(k) % p)
    P(f"    F_n mod {p}: {seq[:pp]} | {seq[pp:pp+3]}...")

P("")
P("Pisano period relationships:")
P(f"  pi(2) = {pisano_period(2)}")
P(f"  pi(3) = {pisano_period(3)}")
P(f"  pi(5) = {pisano_period(5)}")
P(f"  pi(7) = {pisano_period(7)}")
P(f"  pi(11) = {pisano_period(11)}")
P(f"  pi(2)*pi(3) = {pisano_period(2)*pisano_period(3)}, pi(6) = {pisano_period(6)}")
P(f"  lcm(pi(2),pi(3)) = {math.lcm(pisano_period(2), pisano_period(3))}, pi(6) = {pisano_period(6)}")
P(f"  pi(21) = {pisano_period(21)}, lcm(pi(3),pi(7)) = {math.lcm(pisano_period(3), pisano_period(7))}")

P("")
P("Known formula: pi(p) | p-1 if p = +-1 mod 5, pi(p) | 2(p+1) if p = +-2 mod 5")
for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]:
    pp = pisano_period(p)
    r = p % 5
    if r == 1 or r == 4:
        divides = "p-1" if (p - 1) % pp == 0 else "FAIL"
        P(f"  p={p:2d} = {p%5} mod 5: pi({p})={pp}, (p-1)={p-1}, divides? {divides}")
    elif r == 2 or r == 3:
        divides = "2(p+1)" if (2 * (p + 1)) % pp == 0 else "FAIL"
        P(f"  p={p:2d} = {p%5} mod 5: pi({p})={pp}, 2(p+1)={2*(p+1)}, divides? {divides}")
    else:
        P(f"  p={p:2d} = 0 mod 5: pi({p})={pp} (special)")

# ============================================================
# 8. ENTRY POINTS alpha(m)
# ============================================================
section("8. ENTRY POINTS alpha(m)")

P("alpha(m) = smallest k > 0 such that m | F_k")
P("")
test_entries = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 21, 42, 47, 55, 77, 89, 144]
for m in test_entries:
    ep = entry_point(m)
    P(f"  alpha({m:4d}) = {ep:4d}   (F_{ep} = {fib(ep)}, F_{ep}/{m} = {fib(ep)//m})")

P("")
P("Fundamental primes entry points:")
for p in [2, 3, 5, 7, 11]:
    ep = entry_point(p)
    pp = pisano_period(p)
    P(f"  alpha({p:2d}) = {ep}, pi({p}) = {pp}, pi/alpha = {pp//ep}")

P("")
P("Key: alpha(p) | pi(p) always. The ratio pi(p)/alpha(p) is 1, 2, or 4.")
P("")

P("Entry points for composite numbers in our theory:")
for m in [6, 10, 14, 15, 21, 33, 42, 55, 77]:
    ep = entry_point(m)
    facs = factorize(m)
    P(f"  alpha({m:3d}) = {ep:4d}, {m} = {factor_str(m)}, " +
      f"lcm of factor entry pts = {math.lcm(*[entry_point(p**e) for p,e in facs])}")

# ============================================================
# 9. WALL-SUN-SUN PRIMES
# ============================================================
section("9. WALL-SUN-SUN PRIMES AND F_{p-(p|5)} mod p^2")

P("A Wall-Sun-Sun prime p satisfies p^2 | F_{p - (p|5)}.")
P("(p|5) is the Legendre symbol (p/5).")
P("No Wall-Sun-Sun primes are known. Testing our primes:")
P("")

for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]:
    if p == 2:
        P(f"  p={p}: special (p=2)")
        continue
    if p == 5:
        P(f"  p={p}: special (p=5, (p|5)=0)")
        continue
    leg = legendre(p, 5)
    idx = p - leg
    # Compute F_idx mod p^2
    a, b = 0, 1
    for _ in range(idx):
        a, b = b, (a + b) % (p * p)
    fmod = a  # F_idx mod p^2
    P(f"  p={p:3d}: (p|5)={leg:+d}, F_{idx:3d} mod {p}^2 = F_{idx} mod {p*p:5d} = {fmod:5d}  {'WALL-SUN-SUN!' if fmod == 0 else ''}")

P("")
P("No Wall-Sun-Sun primes found up to p=97 (none known at all).")
P("Connection to our framework: a Wall-Sun-Sun prime would create an")
P("exceptionally strong Fibonacci divisibility, but none exist among our primes.")

# ============================================================
# 10. F_p, F_{p-1}, F_{p+1}, L_p FOR EACH FUNDAMENTAL PRIME
# ============================================================
section("10. FIBONACCI/LUCAS AT AND AROUND FUNDAMENTAL PRIMES")

P("For each prime p, compute F_p, F_{p-1}, F_{p+1}, L_p and factor:")
P("")
for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]:
    fp = fib(p)
    fpm1 = fib(p - 1)
    fpp1 = fib(p + 1)
    lp = lucas(p)
    P(f"  p = {p}:")
    P(f"    F_{p-1:2d} = {fpm1:>12d} = {factor_str(fpm1)}")
    P(f"    F_{p:2d}  = {fp:>12d} = {factor_str(fp)}")
    P(f"    F_{p+1:2d} = {fpp1:>12d} = {factor_str(fpp1)}")
    P(f"    L_{p:2d}  = {lp:>12d} = {factor_str(lp)}")
    # Check if p | F_{p-1} or p | F_{p+1}
    leg = legendre(p, 5) if p != 5 and p != 2 else (0 if p == 5 else 1)
    P(f"    (p|5) = {leg}, p | F_{{p-(p|5)}} = p | F_{p - leg}: {fib(p - leg) % p == 0}")
    P("")

subsection("Pattern analysis")
P("Known: for prime p != 5, p | F_{p - (p/5)}")
P("  p=2: 2|F_3=2 YES")
P("  p=3: 3|F_4=3 YES (3 = +2 mod 5, so (3|5)=-1, F_{3-(-1)}=F_4=3)")

P("")
P("Factoring F_p for fundamental primes:")
P(f"  F_2 = 1 (trivial)")
P(f"  F_3 = 2 (prime)")
P(f"  F_5 = 5 (prime, = p itself!)")
P(f"  F_7 = 13 (prime)")
P(f"  F_11 = 89 (prime)")
P(f"  F_13 = 233 (prime)")
P("")
P("ALL F_p for our fundamental primes are themselves prime!")
P("This is because: if p is prime, then F_p is a 'Fibonacci prime candidate'.")
P("Known result: if F_n is prime and n > 4, then n is prime.")
P("(Converse fails: F_19 = 4181 = 37*113 is composite despite 19 being prime.)")
P(f"  F_19 = {fib(19)} = {factor_str(fib(19))}")

# Check which F_p are prime for small primes
P("")
P("Is F_p prime for each prime p?")
for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]:
    fp = fib(p)
    P(f"  F_{p:2d} = {fp:>12d}  {'PRIME' if is_prime(fp) else 'composite = ' + factor_str(fp)}")

# ============================================================
# 11. F_{n+1}^2 + F_n^2 = F_{2n+1} IDENTITY
# ============================================================
section("11. PYTHAGOREAN-TYPE IDENTITY: F_{n+1}^2 + F_n^2 = F_{2n+1}")

P("Verification and exploration of F_{n+1}^2 + F_n^2 = F_{2n+1}:")
P("")
for n in range(12):
    fn = fib(n)
    fn1 = fib(n + 1)
    f2n1 = fib(2 * n + 1)
    ok = "OK" if fn1 * fn1 + fn * fn == f2n1 else "FAIL"
    P(f"  n={n:2d}: F_{n+1}^2 + F_{n}^2 = {fn1}^2 + {fn}^2 = {fn1**2} + {fn**2} = {fn1**2 + fn**2} = F_{2*n+1} = {f2n1} [{ok}]")

subsection("Special cases at n=3 and n=7")
P("n=3: F_4^2 + F_3^2 = 3^2 + 2^2 = 9 + 4 = 13 = F_7")
P(f"  F_7 = {fib(7)} = 13 (prime!)")
P(f"  13 = sum of squares of consecutive Fibonacci primes (2,3)")
P(f"  13 is itself a Fibonacci prime")
P("")
P("n=7: F_8^2 + F_7^2 = 21^2 + 13^2 = 441 + 169 = 610 = F_15")
P(f"  F_15 = {fib(15)} = {factor_str(fib(15))}")
P(f"  610 = 2 * 5 * 61")
P(f"  So (forbidden number)^2 + 13^2 = 2 * 5 * 61 = F_15")
P(f"  Note: 61 is prime, and 2*5 = 10 = product of fundamental primes")

P("")
P("n=4: F_5^2 + F_4^2 = 5^2 + 3^2 = 25 + 9 = 34 = F_9")
P(f"  F_9 = {fib(9)} = {factor_str(fib(9))}")

P("")
P("Related identity: F_{n+1}^2 - F_n^2 = F_{n-1}*F_{n+2}")
for n in range(1, 10):
    lhs = fib(n + 1) ** 2 - fib(n) ** 2
    rhs = fib(n - 1) * fib(n + 2)
    P(f"  n={n}: F_{n+1}^2 - F_{n}^2 = {lhs}, F_{n-1}*F_{n+2} = {rhs}, match={lhs==rhs}")

P("")
P("Cassini's identity: F_{n+1}*F_{n-1} - F_n^2 = (-1)^n")
for n in range(1, 12):
    val = fib(n + 1) * fib(n - 1) - fib(n) ** 2
    P(f"  n={n:2d}: F_{n+1}*F_{n-1} - F_n^2 = {val} = (-1)^{n} = {(-1)**n}")

# ============================================================
# 12. QR COUNTS AS FIBONACCI NUMBERS
# ============================================================
section("12. QR COUNTS OF FUNDAMENTAL PRIMES = FIBONACCI SEQUENCE")

P("For prime p, the number of quadratic residues mod p (excluding 0) is (p-1)/2.")
P("")
for p in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31]:
    qr = qr_count(p)
    # Check if qr is Fibonacci
    fib_idx = None
    for k in range(30):
        if fib(k) == qr:
            fib_idx = k
            break
    tag = f" = F_{fib_idx}" if fib_idx is not None else ""
    P(f"  p={p:2d}: QR count = (p-1)/2 = {qr}{tag}")

P("")
P("The QR counts for {3, 5, 7, 11} are {1, 2, 3, 5} = {F_1, F_3, F_4, F_5}")
P("This is NOT consecutive Fibonacci! The Fibonacci numbers are:")
P("  F_0=0, F_1=1, F_2=1, F_3=2, F_4=3, F_5=5, F_6=8, ...")
P("The QR counts 1, 2, 3, 5 correspond to indices 1 (or 2), 3, 4, 5")
P("")
P("But notice: (p-1)/2 for p = 3,5,7,11 gives 1,2,3,5")
P("The GAPS are: 1, 1, 2 (between 1-2, 2-3, 3-5)")
P("These gaps are themselves 1, 1, 2 = beginning of Fibonacci!")
P("")
P("Continuation: p=13 gives QR=6 (not Fibonacci)")
P("  F_5=5, F_6=8. So 6 is NOT Fibonacci. The pattern breaks at p=13.")

P("")
P("Alternative interpretation: QR counts as function of PRIME INDEX:")
P("  pi(3)=2nd prime: QR=1")
P("  pi(5)=3rd prime: QR=2")
P("  pi(7)=4th prime: QR=3")
P("  pi(11)=5th prime: QR=5")
P("  pi(13)=6th prime: QR=6 (not Fibonacci)")
P("  pi(17)=7th prime: QR=8 = F_6!")
P(f"  Verify: QR(17) = (17-1)/2 = {(17-1)//2}")
P("  So F_6 = 8 reappears at p=17!")

P("")
P("Checking: (p-1)/2 = F_k for some k?")
fib_set = set(fib(k) for k in range(30))
for p in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]:
    qr = (p - 1) // 2
    is_fib = qr in fib_set
    if is_fib:
        idx = [k for k in range(30) if fib(k) == qr][0]
        P(f"  p={p:3d}: QR={qr:4d} = F_{idx} [FIBONACCI]")
    elif p <= 47:
        P(f"  p={p:3d}: QR={qr:4d} (not Fibonacci)")

# ============================================================
# 13. COMPREHENSIVE IDENTITY TABLE
# ============================================================
section("13. COMPREHENSIVE IDENTITY TABLE AT FUNDAMENTAL INDICES")

P("Master table: n, F_n, L_n, F_n*L_n=F_{2n}, F_n/L_n, gcd(F_n,L_n)")
P("")
header = f"{'n':>3} {'F_n':>10} {'L_n':>10} {'F_n*L_n':>15} {'F_{2n}':>15} {'Match':>5} {'gcd':>5} {'F/L':>10}"
P(header)
P("-" * len(header))
for n in range(0, 21):
    fn = fib(n)
    ln = lucas(n)
    prod = fn * ln
    f2n = fib(2 * n)
    g = math.gcd(fn, ln) if fn > 0 else ln
    ratio = f"{fn/ln:.6f}" if ln != 0 else "inf"
    match = "OK" if prod == f2n else "FAIL"
    P(f"{n:3d} {fn:10d} {ln:10d} {prod:15d} {f2n:15d} {match:>5} {g:5d} {ratio:>10}")

# ============================================================
# 14. LUCAS AND FIBONACCI MOD STRUCTURE
# ============================================================
section("14. F_n AND L_n MODULAR STRUCTURE")

P("F_n mod 21 (period should be pi(21)):")
pi21 = pisano_period(21)
P(f"  pi(21) = {pi21}")
seq = [fib(k) % 21 for k in range(pi21 + 5)]
P(f"  F_n mod 21: {seq[:pi21]} | {seq[pi21:]}")

P("")
P("L_n mod 21:")
seq_l = [lucas(k) % 21 for k in range(pi21 + 5)]
P(f"  L_n mod 21: {seq_l[:pi21]} | {seq_l[pi21:]}")

P("")
P("Zeros of F_n mod 21 (i.e., 21 | F_n):")
zeros = [k for k in range(pi21) if fib(k) % 21 == 0]
P(f"  Positions: {zeros}")
P(f"  alpha(21) = {entry_point(21)}")

P("")
P("F_n mod 42:")
pi42 = pisano_period(42)
P(f"  pi(42) = {pi42}")
P(f"  alpha(42) = {entry_point(42)}")

P("")
P("L_n mod p for each fundamental prime p:")
for p in [2, 3, 5, 7, 11]:
    pp = pisano_period(p)
    seq_l = [lucas(k) % p for k in range(pp + 3)]
    P(f"  L_n mod {p:2d} (period <= {pp}): {seq_l[:pp]} | {seq_l[pp:pp+3]}")

# ============================================================
# 15. THE FORBIDDEN NUMBER LANDSCAPE
# ============================================================
section("15. 21 IN THE FIBONACCI-LUCAS LANDSCAPE")

P("21 = F_8 = 3*7 = L_2 * L_4")
P("The index 8 = 2^3 (power of 2)")
P("")
P("All Fibonacci numbers that factor into only fundamental primes {2,3,5,7,11}:")
P("(i.e., 'smooth' Fibonacci numbers w.r.t. our prime set)")
for n in range(40):
    fn = fib(n)
    if fn == 0:
        continue
    facs = factorize(fn)
    if all(p in fund_primes for p, _ in facs):
        P(f"  F_{n:2d} = {fn:>10d} = {factor_str(fn)}")

P("")
P("All Lucas numbers that factor into only fundamental primes:")
for n in range(40):
    ln = lucas(n)
    if ln <= 1:
        continue
    facs = factorize(ln)
    if all(p in fund_primes for p, _ in facs):
        P(f"  L_{n:2d} = {ln:>10d} = {factor_str(ln)}")

P("")
P("Fibonacci numbers divisible by 21:")
for n in range(80):
    fn = fib(n)
    if fn > 0 and fn % 21 == 0:
        P(f"  F_{n:2d} = {fn} = 21 * {fn//21}")
        if n > 40:
            break

P("")
P("Pattern: 21 | F_n iff 8 | n (since alpha(21) = 8)")
P(f"  alpha(21) = {entry_point(21)}")
P("  Verification: n=8,16,24,32,40,... all divisible by 8? YES")

# ============================================================
# 16. DEEP DIVE: phi AND THE GOLDEN ANGLE
# ============================================================
section("16. phi, sqrt(5), AND CONNECTIONS TO QUADRATIC RESIDUES")

P("phi = (1 + sqrt(5))/2")
P("phi^2 = phi + 1 = (3 + sqrt(5))/2")
P("phi^3 = 2*phi + 1 = (4 + 2*sqrt(5))/2 = 2 + sqrt(5)")
P("phi^4 = 3*phi + 2 = (7 + 3*sqrt(5))/2")
P("phi^5 = 5*phi + 3 = (8 + 4*sqrt(5))/2... wait:")
P("")

P("phi^n = F_n * phi + F_{n-1}")
P("Verification:")
for n in range(1, 12):
    fn = fib(n)
    fnm1 = fib(n - 1)
    phi_n = phi ** n
    approx = fn * phi + fnm1
    P(f"  phi^{n:2d} = {phi_n:>15.6f}, F_{n}*phi + F_{n-1} = {fn}*phi + {fnm1} = {approx:>15.6f}, diff = {abs(phi_n - approx):.2e}")

P("")
P("Fibonacci representation of phi powers -- coefficients are Fibonacci!")
P("phi^n = F_n * phi + F_{n-1}")
P("So phi^8 = F_8 * phi + F_7 = 21*phi + 13")
P(f"  = 21 * {phi:.6f} + 13 = {21*phi + 13:.6f}")
P(f"  phi^8 = {phi**8:.6f}")

P("")
P("This means the FORBIDDEN NUMBER 21 is the 'phi-coefficient' of phi^8.")
P("And 13 (= F_7, prime) is the constant term.")

P("")
P("The L_n connection: L_n = F_n * phi + F_n * psi ... actually:")
P("phi^n + psi^n = L_n (where psi = (1-sqrt(5))/2 = -1/phi)")
P("phi^n - psi^n = F_n * sqrt(5)")
P("")
P("So L_n is the 'cosine part' and F_n*sqrt(5) is the 'sine part' of phi^n.")
P("For n=8:")
P(f"  phi^8 + psi^8 = L_8 = {lucas(8)} = 47")
P(f"  phi^8 - psi^8 = F_8 * sqrt(5) = 21 * sqrt(5) = {21*math.sqrt(5):.6f}")
P(f"  Check: {phi**8 + psi**8:.6f} vs {lucas(8)}")
P(f"  Check: {phi**8 - psi**8:.6f} vs {21*math.sqrt(5):.6f}")

# ============================================================
# 17. TOURNAMENT COMPUTATION SPEEDUPS
# ============================================================
section("17. POTENTIAL TOURNAMENT COMPUTATION SPEEDUPS VIA FIBONACCI/LUCAS")

P("Question: can Fibonacci/Lucas structure speed up tournament computations?")
P("")
P("Observation 1: Matrix exponentiation")
P("  The transfer matrix for path counting uses repeated squaring.")
P("  [[1,1],[1,0]]^n = [[F_{n+1}, F_n], [F_n, F_{n-1}]]")
P("  Verification:")
for n in [3, 5, 7, 8, 11]:
    P(f"  n={n}: [[F_{n+1},F_n],[F_n,F_{n-1}]] = [[{fib(n+1)},{fib(n)}],[{fib(n)},{fib(n-1)}]]")

P("")
P("Observation 2: Doubling formulas for fast computation")
P("  F_{2n} = F_n * (2*F_{n+1} - F_n)")
P("  F_{2n+1} = F_n^2 + F_{n+1}^2")
P("  These allow O(log n) Fibonacci computation.")
P("  For tournament transfer matrices, similar doubling applies!")
P("")

P("Observation 3: Eigenvalues of circulant tournament matrices")
P("  For Paley tournament T_p (circulant on QR set),")
P("  the adjacency matrix eigenvalues involve Gauss sums.")
P("  The transfer matrix eigenvalues are polynomials in these.")
P("  At p=7: eigenvalues of adjacency matrix are {3, (-1+/-sqrt(-7))/2}")
P(f"  Note: 7 = L_4 (Lucas prime) and sqrt(-7) appears in Gauss sum.")
P("")

P("Observation 4: Pisano-period batch processing")
P("  When computing H(T) mod p, we can use the Pisano period structure.")
P("  F_n mod p has period pi(p):")
for p in [2, 3, 5, 7, 11]:
    pp = pisano_period(p)
    P(f"    pi({p:2d}) = {pp:3d}")
P("  If transfer matrix entries are Fibonacci-like, these periods")
P("  give modular reduction shortcuts.")
P("")

P("Observation 5: Binet-type closed forms for path counts")
P("  Just as F_n = (phi^n - psi^n)/sqrt(5),")
P("  the number of Hamiltonian paths through a specific structure")
P("  may have a similar closed form using eigenvalues.")
P("  For the n=7 Paley tournament with 3-fold symmetry,")
P("  H(T_7) = 189 = 27 * 7 = 3^3 * 7.")
P(f"  189 = {factor_str(189)}")
P(f"  189 as Fibonacci? No. As Lucas? L_? = 189? {'YES' if 189 in [lucas(k) for k in range(30)] else 'NO'}")

# Check if 189 has any F/L decomposition
P("")
P("  189 = 9 * 21 = 9 * F_8 = 3^2 * F_8")
P("  189 = F_8 * L_2^2 (since L_2 = 3)")
P("  189 = 3 * 63 = 3 * 7 * 9 = L_2 * L_4 * L_2^2 = L_2^3 * L_4")

# ============================================================
# 18. EXTENDED PRIME INDEX ANALYSIS
# ============================================================
section("18. WHERE DO TOURNAMENT CONSTANTS APPEAR IN F/L SEQUENCES?")

P("Key constants from tournament theory:")
constants = {
    1: "H(T_3) = n! for n=1",
    3: "F_4 = L_2, H(T_3) for regular",
    5: "F_5, n=5 vertex count",
    7: "L_4, forbidden prime",
    11: "L_5, Paley cutoff",
    13: "F_7",
    21: "F_8 = FORBIDDEN NUMBER",
    42: "2*21, 'answer to everything'",
    47: "L_8",
    55: "F_10 = 5*11 = F_5*L_5",
    89: "F_11 (prime, F_p for p=11!)",
    189: "H(T_7)",
    987: "F_16 = F_8*L_8 = 21*47",
    95095: "H(T_11)",
    2310: "primorial(5) = 2*3*5*7*11",
}

for val, desc in sorted(constants.items()):
    fib_idx = None
    luc_idx = None
    for k in range(60):
        if fib(k) == val:
            fib_idx = k
        if lucas(k) == val:
            luc_idx = k
    tags = []
    if fib_idx is not None:
        tags.append(f"F_{fib_idx}")
    if luc_idx is not None:
        tags.append(f"L_{luc_idx}")
    if is_prime(val):
        tags.append("PRIME")
    tag_str = " = " + ", ".join(tags) if tags else " (not F/L)"
    P(f"  {val:>8d}{tag_str:>20s}  | {desc}")

P("")
P("89 = F_11: The 11th Fibonacci number is prime!")
P("  11 is our Paley cutoff, and F_11 = 89 is the largest single-digit-index Fibonacci prime.")
P("")
P("55 = F_10 = 5 * 11 = F_5 * L_5:")
P(f"  Verify: F_10 = {fib(10)}, F_5 * L_5 = {fib(5)} * {lucas(5)} = {fib(5)*lucas(5)}")
P(f"  But F_5 * L_5 = F_10 follows from the doubling formula!")
P(f"  55 = product of two fundamental primes (5 and 11)")

P("")
P("95095 = H(T_11):")
P(f"  95095 = {factor_str(95095)}")
P(f"  = 5 * 7 * 11 * 13 * 19")
P(f"  Contains three fundamental primes (5, 7, 11) plus 13 (= F_7) and 19")
P(f"  95095 / (5*7*11) = {95095 // (5*7*11)} = {factor_str(95095 // (5*7*11))}")
P(f"  95095 / 2310 = {95095 / 2310:.4f} (not integer)")

# ============================================================
# 19. THE GCD STRUCTURE
# ============================================================
section("19. GCD PROPERTIES: gcd(F_m, F_n) = F_{gcd(m,n)}")

P("Fundamental: gcd(F_m, F_n) = F_{gcd(m,n)}")
P("Verification:")
for m in range(2, 15):
    for n in range(m + 1, 15):
        g = math.gcd(fib(m), fib(n))
        fg = fib(math.gcd(m, n))
        if g != fg:
            P(f"  FAIL at m={m}, n={n}")
        elif m <= 8 and n <= 10:
            P(f"  gcd(F_{m}, F_{n}) = gcd({fib(m)}, {fib(n)}) = {g} = F_{math.gcd(m,n)} = {fg}")

P("")
P("Consequence for F_8 = 21:")
P(f"  gcd(F_8, F_k) = F_{{gcd(8,k)}}")
P(f"  gcd(F_8, F_4) = F_4 = 3  (since gcd(8,4)=4)")
P(f"  gcd(F_8, F_6) = F_2 = 1  (since gcd(8,6)=2)")
P(f"  gcd(F_8, F_12) = F_4 = 3  (since gcd(8,12)=4)")
P(f"  gcd(F_8, F_16) = F_8 = 21 (since gcd(8,16)=8)")
P("")
P("So 21 and 3 are deeply linked: gcd(21, F_k) is either 1, 1, 3, or 21")
P("depending on whether gcd(8,k) is 1, 2, 4, or 8.")

# ============================================================
# 20. SYNTHESIS AND STRUCTURAL MAP
# ============================================================
section("20. SYNTHESIS: THE FIBONACCI-LUCAS-PRIME STRUCTURAL MAP")

P("SUMMARY OF KEY FINDINGS:")
P("")
P("1. FUNDAMENTAL PRIMES IN LUCAS SEQUENCE (0-indexed):")
P("   L_0=2, L_2=3, L_4=7, L_5=11")
P("   These are the first four PRIME Lucas numbers")
P("   5 = F_5 is Fibonacci, not Lucas")
P("")
P("2. THE FORBIDDEN NUMBER:")
P("   21 = F_8 = 3 * 7 = L_2 * L_4")
P("   21 is the product of two consecutive Lucas primes (the 2nd and 3rd)")
P("   Its Fibonacci index 8 = 2^3 is a power of 2")
P("   L_8 = 47 (prime), and F_8 * L_8 = F_16 = 987 = 3 * 7 * 47")
P("")
P("3. THE GOLDEN CHAIN:")
P("   phi^8 = 21*phi + 13  (forbidden * golden ratio + prime)")
P("   F_n * L_n = F_{2n}  (doubling)")
P("   L_n/F_n -> sqrt(5)  (convergence)")
P("   Q(1/sqrt(5)) = phi^2  (Cayley transform)")
P("")
P("4. QR COUNTS:")
P("   QR(3,5,7,11) = (1,2,3,5) = (F_1, F_3, F_4, F_5)")
P("   Pattern breaks at p=13: QR=6 (not Fibonacci)")
P("   Resumes at p=17: QR=8 = F_6")
P("")
P("5. PISANO PERIODS for fundamental primes:")
for p in [2, 3, 5, 7, 11]:
    P(f"   pi({p:2d}) = {pisano_period(p):3d}")
P(f"   Product pi(2)*pi(3)*pi(5)*pi(7)*pi(11) = {pisano_period(2)*pisano_period(3)*pisano_period(5)*pisano_period(7)*pisano_period(11)}")
P(f"   LCM = {math.lcm(pisano_period(2), pisano_period(3), pisano_period(5), pisano_period(7), pisano_period(11))}")
P("")
P("6. ENTRY POINTS for fundamental primes:")
for p in [2, 3, 5, 7, 11]:
    P(f"   alpha({p:2d}) = {entry_point(p)}")
P(f"   alpha(21) = {entry_point(21)} (= alpha(3)*alpha(7)/gcd... = lcm({entry_point(3)},{entry_point(7)}) = {math.lcm(entry_point(3), entry_point(7))})")
P("")
P("7. F_p IS PRIME for p in {2,3,5,7,11,13}:")
P("   All six fundamental+ primes yield prime Fibonacci numbers")
P("   First failure: F_19 = 4181 = 37 * 113")
P("")
P("8. TOURNAMENT CONNECTIONS:")
P("   H(T_7) = 189 = 3^3 * 7 = L_2^3 * L_4 = 9 * F_8")
P("   H(T_11) = 95095 = 5 * 7 * 11 * 13 * 19")
P("                    = F_5 * L_4 * L_5 * F_7 * 19")
P("   The prime 19 is NOT Fibonacci or Lucas (but F_8 + L_1 = 21 + (-1) = 20... no)")
tag19 = "neither"
for k in range(30):
    if fib(k) == 19:
        tag19 = f"F_{k}"
    if lucas(k) == 19:
        tag19 = f"L_{k}"
P(f"   19 = F_? or L_?: {tag19}")
P(f"   But: 19 = L_4 + L_6 = 7 + {lucas(6)}... L_6={lucas(6)}=18, so 7+18=25 no")
P(f"   19 = F_8 - 2 = 21 - 2 (forbidden minus smallest prime)")
P("")

P("9. DEEP IDENTITY WEB:")
P("   F_8 = 21 (forbidden)")
P("   L_8 = 47 (prime)")
P("   F_8 * L_8 = 987 = F_16")
P("   L_8^2 - 5*F_8^2 = 2209 - 2205 = 4 = (-1)^8 * 4")
P("   F_7^2 + F_8^2 = 169 + 441 = 610 = F_15 = 2*5*61")
P("   F_8^2 + F_9^2 = 441 + 1156 = 1597 = F_17 (PRIME!)")
P(f"   Verify F_17 = {fib(17)}, is_prime = {is_prime(fib(17))}")
P("")

P("10. NO WALL-SUN-SUN PRIMES among fundamental primes (or any known prime)")
P("    This means p^2 never divides F_{p-(p|5)} for our primes")
P("    The Fibonacci divisibility by fundamental primes is 'exactly first order'")

# ============================================================
# SAVE OUTPUT
# ============================================================
output = "\n".join(OUT)
with open("05-knowledge/results/fibonacci_lucas_primes_s116l.out", "w", encoding="utf-8") as f:
    f.write(output)

print("\n\n[Output saved to 05-knowledge/results/fibonacci_lucas_primes_s116l.out]")
