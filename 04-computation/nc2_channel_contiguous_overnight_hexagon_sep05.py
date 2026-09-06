#!/usr/bin/env python3
"""Exact cross-level probes and contiguous interlacing controls.

Local differential identities are proved in the companion note. First/second
signs and adjacent-level gcds below are FINITE-EXACT, not universal claims.
Original charge equations independently reconstruct the complete moment rows.
SymPy supplies exact rational polynomial arithmetic and isolating intervals;
no floating-point root finder or repository mathematical import is used.
"""

from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from math import factorial, gcd, prod

import sympy as sp


t = sp.symbols("t")
CHECKS = 0


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def polynomial(coefficients):
    return sp.Poly.from_list(list(reversed(coefficients)), t, domain=sp.QQ)


def row(A, B, h, r, z, x):
    mass = x+B*h+r+z
    vectors = [(x+(B-A)*j, B*h+r-B*j, z+A*j) for j in range(h+1)]
    coefficients = [factorial(mass)//prod(factorial(k) for k in v) for v in vectors]
    return mass, vectors, polynomial(coefficients)


def original(a, b, c, mass):
    result = {}
    for x in range(mass+1):
        numerator = a*x-b*(mass-x)
        if numerator % (c-b):
            continue
        z = numerator//(c-b)
        y = mass-x-z
        if min(y, z) >= 0:
            result[x, y, z] = factorial(mass)//prod(factorial(k) for k in (x, y, z))
    return result


def canonical(a, b, c, mass):
    g = gcd(a+b, a+c)
    if a*mass % g:
        return None
    A, B, target = (a+b)//g, (a+c)//g, a*mass//g
    z = next(j for j in range(A) if (target-B*j) % A == 0)
    y = (target-B*z)//A
    if y < 0:
        return None
    h, r = divmod(y, B)
    x = mass-y-z
    require(x > 0, "positive-charge canonical left endpoint")
    _, vectors, p = row(A, B, h, r, z, x)
    return A, B, h, r, z, x, vectors, p


def interior_count(p, a, b):
    """SymPy count_roots includes endpoints; product isolators may share one."""
    answer = int(p.count_roots(a, b))
    if a != b:
        answer -= int(p.eval(a) == 0)+int(p.eval(b) == 0)
    return answer


def root_word(p, q):
    require(p.gcd(q).degree() == 0, "no common root before ordered-root audit")
    word = []
    for (a, b), multiplicity in (p*q).intervals():
        require(multiplicity == 1, "simple product roots")
        word.append("P" if interior_count(p, a, b) else "Q")
    require(word.count("P") == p.degree() and word.count("Q") == q.degree(),
            "complete endpoint-safe root labels")
    return "".join(word)


def interval_multiply(left, right):
    values = [a*b for a in left for b in right]
    return min(values), max(values)


def interval_evaluate(p, a, b):
    result = (F(0), F(0))
    for coefficient in p.all_coeffs():
        result = interval_multiply(result, (F(a), F(b)))
        result = result[0]+F(coefficient), result[1]+F(coefficient)
    return result


def signs_at_roots(p, q):
    require(p.gcd(q).degree() == 0, "root-value refinement has no common zero")
    remainder = q.rem(p)
    signs, certificates = [], []
    for (a, b), multiplicity in p.intervals():
        require(multiplicity == 1 and b <= 0 and p.eval(0) != 0,
                "simple negative first roots, allowing a non-root zero interval endpoint")
        for _ in range(200):
            lo, hi = interval_evaluate(remainder, a, b)
            if lo > 0 or hi < 0:
                break
            a, b = p.refine_root(a, b, eps=(b-a)/4)
        else:
            raise RuntimeError("root-value interval refinement did not terminate")
        signs.append(1 if lo > 0 else -1)
        certificates.append((str(a), str(b), str(lo), str(hi)))
    require(len(signs) == p.degree(), "all first roots included in sign query")
    return signs, certificates


def coordinate_controls():
    counts = Counter()
    for B in range(2, 5):
        for A in range(1, B):
            for h in range(1, 4):
                for r in range(B):
                    for z in range(A):
                        for x in (0, 2):
                            m, _, p = row(A, B, h, r, z, x)
                            updates = [("x", (A, B, h, r, z, x+1), x+1, B-A, "QP"*h)]
                            if z+1 < A:
                                updates.append(("z", (A, B, h, r, z+1, x), z+1, A, "QP"*h))
                            if r+1 < B:
                                updates.append(("r", (A, B, h, r+1, z, x), B*h+r+1, -B, "PQ"*h))
                            else:
                                updates.append(("r_carry", (A, B, h+1, 0, z, x), B*(h+1), -B,
                                                "QP"*h+"Q"))
                            for kind, parameters, constant, slope, expected_word in updates:
                                mm, _, q = row(*parameters)
                                require(mm == m+1, "one-coordinate mass increment")
                                require((constant*q+slope*sp.Poly(t, t)*q.diff()) == (m+1)*p,
                                        "exact Euler contiguous relation including normalization")
                                require(root_word(p, q) == expected_word, "strict directed interlacing")
                                counts[kind] += 1
    return counts


def first_second_controls():
    counts, carries = Counter(), Counter()
    digest = sha256()
    for B in range(2, 7):
        for A in range(1, B):
            if gcd(A, B) != 1:
                continue
            for h in (2, 3, 4):
                for r in range(B):
                    for z in range(A):
                        for x in (1, 2, 3):
                            a = A*(B*h+r)+B*z
                            g = x+B*h+r+z
                            b, c = g*A-a, g*B-a
                            if b <= 0 or gcd(a, g) != 1:
                                continue
                            require(c > b and gcd(gcd(a, b), c) == 1, "eligible primitive first-row support")
                            first = canonical(a, b, c, g)
                            second = canonical(a, b, c, 2*g)
                            require(first[:6] == (A, B, h, r, z, x), "first parameter recovery")
                            ey, ez = 2*r//B, 2*z//A
                            require(second[:6] == (A, B, 2*h+ey+ez, 2*r-B*ey, 2*z-A*ez,
                                                   2*x-(B-A)*ez), "both actual doubling carries")
                            for mass, result in ((g, first), (2*g, second)):
                                vectors, poly = result[-2:]
                                expected = {v: int(poly.nth(j)) for j, v in enumerate(vectors)}
                                require(original(a, b, c, mass) == expected, "independent full weighted charge row")
                            p, q = first[-1], second[-1]
                            signs, certificates = signs_at_roots(p, q)
                            require(all(s*(-1)**ez == -1 for s in signs),
                                    "FINITE-EXACT normalized second sign at every first root")
                            counts[h] += 1
                            carries[ey, ez] += 1
                            digest.update(repr((A, B, h, r, z, x, signs, certificates)).encode())
    return counts, carries, digest.hexdigest()


def adjacent_controls():
    supports, pairs, jumps = 0, 0, Counter()
    first_drop = None
    for a in range(1, 9):
        for b in range(1, 6):
            for c in range(b+1, 11):
                if gcd(gcd(a, b), c) != 1:
                    continue
                supports += 1
                previous = None
                for mass in range(1, 41):
                    result = canonical(a, b, c, mass)
                    if result is None:
                        continue
                    q = result[-1]
                    if previous is not None:
                        old_mass, p = previous
                        if min(p.degree(), q.degree()) > 0:
                            require(p.gcd(q).degree() == 0, "FINITE-EXACT adjacent-level coprimality")
                            pairs += 1
                            jump = q.degree()-p.degree()
                            jumps[jump] += 1
                            if jump < 0 and first_drop is None:
                                first_drop = ((-a, b, c), old_mass, mass, p.degree(), q.degree())
                    previous = mass, q
    return supports, pairs, jumps, first_drop


def large_offset_controls():
    counts = Counter()
    for h in (3, 4, 5, 6, 8, 10, 12):
        for x in (5, 7, 11, 101, 503, 997):
            if gcd(2*h, x) != 1:
                continue
            def normalized(H, X):
                return polynomial([F(factorial(2*H), factorial(j)*factorial(2*H-2*j)
                                     *prod(range(X+1, X+j+1))) for j in range(H+1)])
            p, q = normalized(h, x), normalized(2*h, 2*x)
            signs, _ = signs_at_roots(p, q)
            require(all(s == -1 for s in signs), "FINITE-EXACT large-offset hostile-direction test")
            counts[h] += 1
    return counts


def paired_product_hostile():
    # f=u^-4+u+t*u^6, q=f^5; compute every frequency from the original trinomial.
    by_frequency = {}
    for x in range(6):
        for z in range(6-x):
            y = 5-x-z
            exponent = -4*x+y+6*z
            coefficient = factorial(5)//(factorial(x)*factorial(y)*factorial(z))
            by_frequency.setdefault(exponent, {})[z] = coefficient
    q = {e: polynomial([co.get(j, 0) for j in range(max(co)+1)]) for e, co in by_frequency.items()}
    p = q[0]
    remainders, sign_rows = {}, {}
    for e in (5, 10, 15, 20):
        pair = q[e]*q[-e]
        remainders[e] = pair.rem(p)
        if remainders[e].is_zero:
            sign_rows[e] = [0, 0]
        else:
            sign_rows[e], _ = signs_at_roots(p, pair)
    require(sign_rows == {5: [-1, -1], 10: [0, 0], 15: [-1, 1], 20: [1, -1]},
            "individual opposite-frequency products have both signs")
    require(remainders[15] == -20*remainders[20], "exact grouped-product sidecar")
    require(remainders[5] == polynomial([560, 3220]), "first paired remainder")
    require(remainders[15] == polynomial([125, 700]), "mixed-sign paired remainder")
    second = canonical(4, 1, 6, 10)[-1]
    require((p*p+2*sum(q[e]*q[-e] for e in q if e > 0 and -e in q)) == second,
            "complete constant-term square convolution")
    return sorted(sign_rows.items())


def hostile_controls():
    p = canonical(3, 1, 5, 8)[-1]
    q = canonical(3, 1, 5, 12)[-1]
    require(p == polynomial([28, 280, 420, 56]), "literal mass-eight hostile")
    require(q == polynomial([220, 3960, 16632, 18480, 3960]), "literal mass-twelve hostile")
    require(root_word(p, q) == "PQQPQPQ", "consecutive admissible masses do not interlace")
    intervals = [("P", F(-7), F(-6)), ("Q", F(-4), F(-3)),
                 ("Q", F(-4, 5), F(-3, 4)), ("P", F(-2, 3), F(-1, 2)),
                 ("Q", F(-1, 3), F(-1, 4)), ("P", F(-1, 8), F(-1, 9)),
                 ("Q", F(-1, 12), F(-1, 13))]
    for letter, a, b in intervals:
        poly = p if letter == "P" else q
        require(poly.eval(a)*poly.eval(b) < 0, "independent rational sign bracket")
    require(all(a[2] < b[1] for a, b in zip(intervals, intervals[1:])), "ordered disjoint hostile brackets")
    # Opposite signs show CT of a square is not a positive norm, even at zero mean.
    square_values = []
    for terms in ({1: 1, -1: 1}, {1: 1, -1: -1}):
        square_values.append(sum(value*terms.get(-exponent, 0) for exponent, value in terms.items()))
    require(square_values == [2, -2], "constant-term square is an indefinite bilinear form")
    return str(p.as_expr()), str(q.as_expr())


def main():
    coordinate = coordinate_controls()
    first, carries, digest = first_second_controls()
    supports, adjacent, jumps, first_drop = adjacent_controls()
    offsets = large_offset_controls()
    paired = paired_product_hostile()
    hostile = hostile_controls()
    print("PROVED local identities; cross-mass sign/gcd universes are FINITE-EXACT only.")
    print("Exact polynomial backend: SymPy", sp.__version__)
    print("Exact local Euler/interlacing pairs:", sorted(coordinate.items()))
    print("Higher-multiplicity first/second sign rows:", sum(first.values()), sorted(first.items()))
    print("First/second carry types:", sorted(carries.items()))
    print("Root-value certificate SHA-256:", digest)
    print("Adjacent-return support universe:", supports, "coprime pairs:", adjacent)
    print("Adjacent-return degree jumps:", sorted(jumps.items()))
    print("First degree-drop control:", first_drop)
    print("Large-offset first/second sign rows:", sum(offsets.values()), sorted(offsets.items()))
    print("Opposite-frequency paired-product hostile:", paired)
    print("Interlacing hostile P:", hostile[0])
    print("Interlacing hostile Q:", hostile[1])
    print("Interlacing hostile exact root order: PQQPQPQ")
    print("Explicit gates:", CHECKS)


if __name__ == "__main__":
    main()
