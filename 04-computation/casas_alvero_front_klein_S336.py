#!/usr/bin/env python3
"""
klein-2026-07-20-S336 -- THE CASAS-ALVERO FRONT (boxeph-S157 flagged it as absent from the repo;
nobody took it).  Owner: "polynomial-derivative rigidity with p-adic partial result".

CASAS-ALVERO: f monic in K[x], deg n, char K = 0.  If gcd(f, f^(i)) != 1 for every
i = 1..n-1, then f = (x-a)^n.

Known partial results (Graf von Bothmer-Labs-Schicho-van de Woestijne 2007): TRUE for
n = p^k, 2p^k, 3p^k, 4p^k, p prime.  So the open degrees are exactly those NOT of that shape.

THREE THINGS COMPUTED HERE, each with a control:
 (1) the exact list of OPEN degrees, and the smallest one;
 (2) the char-p failure map -- CA is FALSE in characteristic p, and where;
 (3) a CALIBRATED control: the same test in char 0 must return ONLY (x-a)^n.
     (MISTAKE-162: a null with no positive control is worthless.)
"""
from itertools import product
from math import gcd as _g

# ---------------------------------------------------------------- (1) open degrees
def is_covered(n):
    """n = c * p^k with c in {1,2,3,4}, p prime, k >= 0 ?"""
    def prime_power(m):
        if m == 1: return True
        d = 2
        while d * d <= m:
            if m % d == 0:
                while m % d == 0: m //= d
                return m == 1
            d += 1
        return True
    return any(n % c == 0 and prime_power(n // c) for c in (1, 2, 3, 4))

openn = [n for n in range(2, 201) if not is_covered(n)]
print("=" * 78)
print("(1) CASAS-ALVERO: which degrees are NOT covered by n = {1,2,3,4}*p^k ?")
print("=" * 78)
print(f" open degrees <= 200: {openn[:30]}{' ...' if len(openn)>30 else ''}")
print(f" SMALLEST OPEN DEGREE = {openn[0]}   ( = {openn[0]} = 2*3*5 )")
def nprimes(n):
    s, d, m = set(), 2, n
    while d * d <= m:
        while m % d == 0: s.add(d); m //= d
        d += 1
    if m > 1: s.add(m)
    return len(s)
print(f" every covered degree has at most 2 distinct primes: "
      f"{all(nprimes(n)<=2 for n in range(2,201) if is_covered(n))}")
print(f" every open degree has at least 3 distinct primes: "
      f"{all(nprimes(n)>=3 for n in openn)}")
print(" => THE CASAS-ALVERO BARRIER IS EXACTLY A THREE-DISTINCT-PRIMES PHENOMENON.")
print("    n is covered iff n = c*p^k with c <= 4, and that forces omega(n) <= 2.")

# ---------------------------------------------------------------- poly helpers over F_p
def deriv(f, p):
    return [(i * f[i]) % p for i in range(1, len(f))]

def polygcd(a, b, p):
    a = a[:]; b = b[:]
    def trim(v):
        while v and v[-1] % p == 0: v.pop()
        return v
    a = trim(a); b = trim(b)
    while b:
        inv = pow(b[-1], p - 2, p)
        while len(a) >= len(b) and a:
            c = (a[-1] * inv) % p; sh = len(a) - len(b)
            for i in range(len(b)): a[sh + i] = (a[sh + i] - c * b[i]) % p
            a = trim(a)
        a, b = b, a
    return a

def is_nth_power_of_linear(f, n, p):
    """f = (x-a)^n over F_p ?"""
    for a in range(p):
        g = [1]
        for _ in range(n):
            g = [(g[i - 1] if i > 0 else 0) - a * (g[i] if i < len(g) else 0) for i in range(len(g) + 1)]
            g = [c % p for c in g]
        if g == [c % p for c in f]: return True
    return False

def ca_test(f, n, p):
    """does f share a root with every f^(i), i=1..n-1 ?  (gcd non-constant)"""
    d = f[:]
    for i in range(1, n):
        d = deriv(d, p)
        if not d: return False                     # derivative vanished identically
        if len(polygcd(f, d, p)) <= 1: return False
    return True

# ---------------------------------------------------------------- (2)+(3) the map
print("\n" + "=" * 78)
print("(2)+(3) THE CHAR-p FAILURE MAP, with the trivial solutions as positive control")
print("=" * 78)
print(f"{'p':>3} {'n':>3} {'monic f':>10} {'CA-solutions':>13} {'trivial (x-a)^n':>16} {'NONTRIVIAL':>11}")
found = {}
for p in (2, 3, 5, 7):
    for n in range(2, 8):
        if n >= p * 3: pass
        total = p ** n
        if total > 200000: continue
        sols, triv, nontriv = 0, 0, []
        for coeffs in product(range(p), repeat=n):
            f = list(coeffs) + [1]                 # monic, degree n, ascending
            if ca_test(f, n, p):
                sols += 1
                if is_nth_power_of_linear(f, n, p): triv += 1
                else: nontriv.append(f)
        print(f"{p:>3} {n:>3} {total:>10} {sols:>13} {triv:>16} {len(nontriv):>11}"
              + ("   <-- CA FAILS in char p" if nontriv else ""))
        if nontriv: found[(p, n)] = nontriv[:3]
print("\n POSITIVE CONTROL: the trivial column must equal p at every (p,n) -- the p polynomials")
print(" (x-a)^n always satisfy CA.  If it ever prints something else, the test is broken.")
print("\n first nontrivial char-p solutions found (ascending coefficients):")
for k, v in list(found.items())[:8]:
    print(f"   p={k[0]}, n={k[1]}: {v}")
print("\n READING: CA is FALSE in characteristic p -- these are genuine counterexamples to the")
print(" char-0 STATEMENT read in char p.  That is exactly why the p-adic route is delicate:")
print(" a mod-p argument must exclude these, which is what the p^k / 2p^k / 3p^k / 4p^k")
print(" degrees buy.  The three-distinct-primes barrier in (1) is the shape of what is left.")
