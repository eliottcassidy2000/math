#!/usr/bin/env python3
"""Independent exact controls for simple-negative trinomial channel roots.

No repository mathematics and no floating point are imported. The first
producer uses canonical factorial rows; the second scans original charge
equations at every mass. Root counts use an elementary rational Sturm chain,
not the Laguerre factorization that proves the infinite theorem. All gates
remain active with Python -O. The finite universe is declared in main().
"""

from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from math import comb, factorial, gcd, lcm, prod


CHECKS = 0


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def rising(a, j):
    return prod((a+i for i in range(j)), start=F(1))


def falling(a, j):
    return prod((a-i for i in range(j)), start=F(1))


def generalized_binomial(a, j):
    return falling(a, j)/factorial(j)


def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def primitive(p):
    p = trim(p)
    denominator = lcm(*(c.denominator for c in p))
    integers = [int(c*denominator) for c in p]
    content = gcd(*integers)
    return [F(c//content) for c in integers] if content else [F(0)]


def remainder(p, q):
    p, q = list(p), trim(q)
    while p != [0] and len(p) >= len(q):
        shift, coefficient = len(p)-len(q), p[-1]/q[-1]
        for j, a in enumerate(q):
            p[j+shift] -= coefficient*a
        p = trim(p)
    return primitive(p)


def sign(x):
    return (x > 0)-(x < 0)


def variations(signs):
    signs = [s for s in signs if s]
    return sum(a != b for a, b in zip(signs, signs[1:]))


def root_audit(coefficients):
    """Count distinct roots in (-infinity,0) and independently test gcd(P,P')."""
    p = primitive([F(a) for a in coefficients])
    degree = len(p)-1
    require(all(a > 0 for a in p), "positive coefficients and both nonzero ends")
    if degree == 0:
        return 0
    chain = [p, primitive([F(j)*p[j] for j in range(1, len(p))])]
    while True:
        rem = remainder(chain[-2], chain[-1])
        if rem == [0]:
            break
        chain.append([-a for a in rem])
    at_minus_infinity = [sign(q[-1])*(-1)**(len(q)-1) for q in chain]
    at_zero = [sign(q[0]) for q in chain]
    negative_count = variations(at_minus_infinity)-variations(at_zero)
    require(negative_count == degree, ("all roots strictly negative", degree, negative_count, p))
    require(len(chain[-1]) == 1, "squarefree by Euclidean gcd")
    return degree


def row(A, B, h, r, z, x):
    delta = B-A
    mass = x+B*h+r+z
    vectors = [(x+delta*j, B*h+r-B*j, z+A*j) for j in range(h+1)]
    return vectors, [factorial(mass)//prod(factorial(u) for u in v) for v in vectors]


def split_audit(A, B, h, r, z, x):
    require(0 < A < B and h >= 0 and x >= 0 and 0 <= r < B and 0 <= z < A,
            "weak integer parameter hypotheses")
    vectors, coefficients = row(A, B, h, r, z, x)
    delta = B-A
    eta = [F(B*h+r-s, B) for s in range(B) if s != r]
    mu = [F(x+i, delta) for i in range(1, delta+1)]
    nu = [F(z+i, A) for i in range(1, A+1) if i != A-z]
    scale = F(B**B, delta**delta*A**A)
    require(len(eta) == len(mu)+len(nu) == B-1, "balanced factor count")
    require(all(e > h-1 for e in eta) and all(u > 0 for u in mu+nu),
            "strict finite Laguerre parameter ranges")
    for j, coefficient in enumerate(coefficients):
        split = comb(h, j)*scale**j
        split *= prod((falling(e, j) for e in eta), start=F(1))
        split /= prod((rising(u, j) for u in mu+nu), start=F(1))
        require(F(coefficient, coefficients[0]) == split, "exact split coefficient")
    root_audit(coefficients)
    return vectors, coefficients


def original_channels(a, b, c, mass):
    result = {}
    for x in range(mass+1):
        numerator = a*x-b*(mass-x)
        if numerator % (c-b):
            continue
        z = numerator//(c-b)
        y = mass-x-z
        if min(y, z) >= 0:
            result[x, y, z] = factorial(mass)//prod(factorial(t) for t in (x, y, z))
    return result


def charge_audit(a, b, c, mass):
    original = original_channels(a, b, c, mass)
    g = gcd(a+b, a+c)
    if not original:
        return None
    require(a*mass % g == 0, "charge congruence necessity")
    A, B, target = (a+b)//g, (a+c)//g, a*mass//g
    z = next(t for t in range(A) if (target-B*t) % A == 0)
    y = (target-B*z)//A
    h, r = divmod(y, B)
    x = mass-y-z
    require(h >= 0 and x > 0, "positive-charge complete canonical fibre")
    vectors, coefficients = split_audit(A, B, h, r, z, x)
    require(dict(zip(vectors, coefficients)) == original, "independent complete weighted charge row")
    return h, coefficients


def symbol_audit():
    rows = 0
    for n in range(1, 9):
        for numerator in range(1, 15):
            eta, mu = F(n-1)+F(numerator, 7), F(numerator, 7)
            falling_symbol = [comb(n, j)*falling(eta, j) for j in range(n+1)]
            inverse_symbol = [F(comb(n, j))/rising(mu, j) for j in range(n+1)]
            # Independently expand n! t^n L_n^(eta-n)(-1/t).
            reverse_laguerre = [factorial(n)*generalized_binomial(eta, j)/factorial(n-j)
                                for j in range(n+1)]
            # Independently expand n!/(mu)_n L_n^(mu-1)(-t).
            direct_laguerre = [F(factorial(n))*generalized_binomial(n+mu-1, n-j)
                              /(rising(mu, n)*factorial(j)) for j in range(n+1)]
            require(falling_symbol == reverse_laguerre, "falling/reversed Laguerre identity")
            require(inverse_symbol == direct_laguerre, "inverse-rising/Laguerre identity")
            root_audit(falling_symbol)
            root_audit(inverse_symbol)
            rows += 1
    return rows


def main():
    symbols = symbol_audit()
    canonical, noncoprime = Counter(), 0
    digest = sha256()
    # No gcd or charge positivity filter: the weak factorial theorem is broader.
    for B in range(2, 6):
        for A in range(1, B):
            for h in range(6):
                for r in range(B):
                    for z in range(A):
                        for x in range(3):
                            _, coefficients = split_audit(A, B, h, r, z, x)
                            canonical[h] += 1
                            noncoprime += gcd(A, B) > 1
                            digest.update(repr((A, B, h, r, z, x, coefficients)).encode())
    charge_rows, empty_rows, contents, degree_counts = 0, 0, Counter(), Counter()
    # All supports in the box, with no primitive-content restriction, every mass.
    for a in range(1, 6):
        for b in range(1, 5):
            for c in range(b+1, 8):
                for mass in range(1, 19):
                    result = charge_audit(a, b, c, mass)
                    if result is None:
                        empty_rows += 1
                    else:
                        h, coefficients = result
                        charge_rows += 1
                        contents[gcd(gcd(a, b), c)] += 1
                        degree_counts[h] += 1
                        digest.update(repr((a, b, c, mass, coefficients)).encode())
    controls = [(4, 1, 6, 5), (13, 1, 8, 7), (13, 1, 8, 14),
                (15, 33, 49, 16), (15, 33, 49, 32)]
    control_rows = []
    for a, b, c, mass in controls:
        h, coefficients = charge_audit(a, b, c, mass)
        control_rows.append(((-a, b, c), mass, h, coefficients))
    # Wide first fibres, including h=12, x=0, and noncoprime A,B controls.
    wide_rows = 0
    for A, B in ((1, 2), (2, 3), (3, 7), (7, 12), (2, 4)):
        for h in (6, 8, 12):
            for r, z, x in ((0, 0, 0), (B-1, A-1, 3)):
                split_audit(A, B, h, r, z, x)
                wide_rows += 1
    resonance_rows = 0
    for p in range(1, 5):
        for q in range(1, 5):
            for mass in range(19):
                h, r = divmod(mass, p+q)
                _, coefficients = split_audit(q, p+q, h, r, 0, 0)
                direct = [factorial(mass)//(factorial(p*j)*factorial(q*j)
                          *factorial(mass-(p+q)*j)) for j in range(h+1)]
                require(coefficients == direct, "all-charge proportional-resonance corollary")
                resonance_rows += 1
    # Exact boundaries: positivity alone and a noncanonical deletion fail.
    require(1-4 < 0, "1+t+t^2 is positive but not real-rooted")
    _, quadratic = row(1, 2, 2, 0, 0, 1)
    require(quadratic == [5, 30, 10], "higher-multiplicity hostile control")
    require(-4*quadratic[0]*quadratic[2] < 0,
            "deleting the middle channel gives nonreal roots")
    root_audit([1, 1])
    root_audit([1, 3, 2])
    require(sum(F(a)*F(-1)**j for j, a in enumerate([1, 1])) == 0
            and sum(F(a)*F(-1)**j for j, a in enumerate([1, 3, 2])) == 0,
            "root location and simplicity do not imply two-polynomial coprimality")
    print("FINITE-EXACT controls; infinite theorem uses the cited finite-preserver input.")
    print("Laguerre symbol pairs:", symbols)
    print("Canonical integer rows:", sum(canonical.values()), "by degree:", sorted(canonical.items()))
    print("Noncoprime canonical rows included:", noncoprime)
    print("Original charge rows:", charge_rows, "empty:", empty_rows,
          "by content:", sorted(contents.items()))
    print("Original charge rows by degree:", sorted(degree_counts.items()))
    print("Wide independent Sturm rows:", wide_rows)
    print("Direct proportional-resonance rows:", resonance_rows)
    for control in control_rows:
        print("Control:", control)
    print("Semantic SHA-256:", digest.hexdigest())
    print("Explicit gates:", CHECKS)


if __name__ == "__main__":
    main()
