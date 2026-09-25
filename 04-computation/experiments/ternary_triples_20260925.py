"""Primitive triple denominator filtration and guarded signed Collatz lift.

Source note: 05-knowledge/results/ternary_triples_20260925.md.
All claims of exhaustive coverage below refer to explicit finite universes.
No global Collatz convergence assumption or third-party package is used.
"""

from fractions import Fraction
from math import gcd, isqrt


def check(test, label):
    if not test:
        raise RuntimeError(label)


def valuation(n, prime):
    check(n > 0, "positive valuation input")
    exponent = 0
    while n % prime == 0:
        n //= prime
        exponent += 1
    return exponent


def triple(s, t):
    check(s > 0 and t > 0 and s % 2 == t % 2 == 1 and gcd(s, t) == 1,
          "primitive odd spinor domain")
    return s*t, (s*s-t*t)//2, (s*s+t*t)//2


def decode(point):
    a, b, c = point
    check(a > 0 and c > 0 and a*a+b*b == c*c and gcd(a, abs(b)) == 1,
          "oriented primitive triple")
    check(a % 2 == 1 and b % 2 == 0, "marked odd and even legs")
    s, t = isqrt(c+b), isqrt(c-b)
    check(s*s == c+b and t*t == c-b and triple(s, t) == point,
          "square-gap inverse")
    return s, t


def lift(s, t, sign=1):
    triple(s, t)
    check(sign in (-1, 1), "sign domain")
    raw = 3*s+sign*t
    check(raw > 0, "positive-sheet guard")
    k = valuation(raw, 2)
    g = gcd(raw, t)
    check(k >= 1 and g == gcd(3, t), "exact dyadic and primitive split")
    out = raw//(2**k*g), t//g
    triple(*out)
    return out, k, g


def native_rational(x, sign=1):
    z = 3*x+sign
    check(z > 0, "native positive-sheet guard")
    k = 0
    while z.numerator % 2 == 0:
        z /= 2
        k += 1
    return z, k


def inverse_parent(u, v, k, sign=1):
    check(k >= 1, "positive exponent")
    parent = (2**k*Fraction(u, v)-sign)/3
    check(parent > 0, "positive inverse parent")
    return parent.numerator, parent.denominator


def main():
    pairs = signed_edges = contractions = 0
    for s in range(1, 256, 2):
        for t in range(1, 128, 2):
            if gcd(s, t) != 1:
                continue
            point = triple(s, t)
            a, b, c = point
            check(decode(point) == (s, t), "spinor/triple inverse")
            check(a*a+b*b == c*c and gcd(a, abs(b)) == 1, "primitive identity")
            check(c % 3 != 0 and ((a % 3 == 0) != (b % 3 == 0)), "one leg divisible by3")
            check((a % 3 == 0) == (s % 3 == 0 or t % 3 == 0), "which leg carries3")
            check(b % 4 == 0, "even leg dyadic seam")
            if s != t:
                check(valuation(abs(b), 2) == valuation(abs(s-t), 2)+valuation(s+t, 2)-1,
                      "even leg exact valuation")
            pairs += 1
            for sign in (-1, 1):
                if 3*s+sign*t <= 0:
                    continue
                (u, v), k, g = lift(s, t, sign)
                native, native_k = native_rational(Fraction(s, t), sign)
                check((Fraction(u, v), k) == (native, native_k), "native rational conjugacy")
                check(inverse_parent(u, v, k, sign) == (s, t), "inverse exact edge")
                check(valuation(v, 3) == max(valuation(t, 3)-1, 0), "denominator level")
                check(v//3**valuation(v, 3) == t//3**valuation(t, 3), "3-free denominator invariant")
                # Cancelling the common3 before or after all halvings commutes.
                check((3*s+sign*t)//g//2**k == u, "coprime cancellation commutes with halving")
                if t % 3 == 0:
                    c2 = triple(u, v)[2]
                    check((3 if sign == 1 else 4)*c2 < c, "strict intrinsic hypotenuse contraction")
                    contractions += 1
                signed_edges += 1
    print("ODD COPRIME PAIRS: s<=255,t<=127; pairs", pairs,
          "; native signed edges", signed_edges, "; off-axis contractions", contractions)

    # Euclid-free primitive-triple enumeration independently checks chart coverage.
    triangles = 0
    for c in range(1, 301):
        for a in range(1, c+1, 2):
            b = isqrt(c*c-a*a)
            if b*b != c*c-a*a or b % 2 or gcd(a, b) != 1:
                continue
            for signed_b in ({b, -b} if b else {0}):
                check(triple(*decode((a, signed_b, c))) == (a, signed_b, c), "independent chart coverage")
                triangles += 1
    print("EUCLID-FREE COVERAGE: oriented primitive triples with c<=300:", triangles)

    entries = inverse_edges = 0
    for r in range(7):
        t = 3**r
        for s in range(1, 256, 2):
            if gcd(s, t) != 1:
                continue
            u, v = s, t
            for expected in range(r-1, -1, -1):
                old_c = triple(u, v)[2]
                (u, v), _, _ = lift(u, v)
                check(v == 3**expected and 3*triple(u, v)[2] < old_c, "exact plus entry clock")
            check(v == 1, "entered integer axis")
            entries += 1
            for sign in (-1, 1):
                for k in range(1, 13):
                    if 2**k*s-sign*t <= 0:
                        continue
                    p, q = inverse_parent(s, t, k, sign)
                    check(lift(p, q, sign)[:2] == ((s, t), k), "inverse exponent passport")
                    if r:
                        check(q == 3*t, "off-axis reverse raises level exactly once")
                    else:
                        check((q == 1) == ((2**k*s-sign) % 3 == 0), "integer-axis ternary guard")
                    inverse_edges += 1
    print("POWER9 GAP FILTRATION: r=0..6,s<=255; exact entries", entries,
          "; signed inverse edges k=1..12", inverse_edges)

    point = triple(13, 9)
    path = [point]
    for _ in range(2):
        spinor, _, _ = lift(*decode(path[-1]))
        path.append(triple(*spinor))
    check(path == [(117, 44, 125), (3, -4, 5), (1, 0, 1)], "orientation crossing witness")
    print("PLUS ORIENTATION CROSSING:", path)
    check(lift(3, 1)[0] == (5, 1) and lift(1, 3)[0] == (1, 1), "reciprocal collapse hostile")
    print("UNMARKED COLLISION: (3,+4,5)->(5,12,13), but (3,-4,5)->(1,0,1)")

    for k in range(3, 21):
        t = 2**k-3
        check(t % 3 != 0 and lift(1, t)[:2] == ((1, t), k), "noninteger fixed stratum")
    print("FULL-TRIPLE HOSTILE: x=1/(2^k-3) fixed for k>=3; k=3..20 checked; (5,-12,13) is first")

    for r in range(1, 21):
        numerator = 2*8**(r-1)+3**r
        check(numerator % 5 == 0, "counter-family integrality")
        s, t = numerator//5, 3**r
        for expected in range(r, 0, -1):
            (s, t), k, _ = lift(s, t)
            check(k == (1 if expected == 1 else 3), "counter-family exact exponents")
        check((s, t) == (1, 1), "counter-family terminal")
    print("COUNTER FAMILY: s_r=(2*8^(r-1)+3^r)/5,t_r=3^r; r=1..20; exact word 3^(r-1),1")

    local_hostiles = 0
    for h in range(1, 33):
        points = []
        for k in (h+3, h+4):
            s, t = 2**k-3, 9
            point = triple(s, t)
            check(lift(s, t)[1] == k, "arbitrary dyadic carry")
            points.append(tuple(x % 2**h for x in point))
        check(points[0] == points[1], "identical local triple data with different carries")
        local_hostiles += 1
    print("DYADIC LOCAL OBSTRUCTION: equal full triples mod2^h with distinct halving counts; h=1..32:", local_hostiles)

    # The positive minus domain genuinely differs on the rational extension.
    s, t = 29, 27
    minus_path = [triple(s, t)]
    for _ in range(2):
        (s, t), _, _ = lift(s, t, -1)
        minus_path.append(triple(s, t))
    check(minus_path == [(783, 56, 785), (45, -28, 53), (3, -4, 5)], "minus boundary approach")
    try:
        lift(s, t, -1)
    except RuntimeError:
        pass
    else:
        raise RuntimeError("zero was accepted as odd rational state")
    print("MINUS POSITIVITY HOSTILE:", minus_path, "then 3x-1=0, rejected")
    for seed, period in ((1, 1), (5, 2), (17, 7)):
        s, t = seed, 1
        for _ in range(period):
            (s, t), _, _ = lift(s, t, -1)
        check((s, t) == (seed, 1), "minus integer-axis cycle")
    print("MINUS INTEGER AXIS: cycles at 1,5,17 retained; off-axis contraction does not settle the boundary")
    print("ALL EXACT CHECKS PASS; universal integer-axis root generation remains OPEN")


if __name__ == "__main__":
    main()
