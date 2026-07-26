#!/usr/bin/env python3
"""Exact companion for THM-2412.

Checks the scaled Gregory--Newton lowering identity, the umbral
intertwiner, terminating discrete exponentials, the central Bernoulli-triangle
split, Catalan leakage recurrences, and the tournament Hamming-shell count.
Only standard-library exact rational arithmetic is used.
"""

from fractions import Fraction
from math import comb, factorial


def trim(p):
    q = list(p)
    while len(q) > 1 and q[-1] == 0:
        q.pop()
    return tuple(q)


def add(p, q):
    n = max(len(p), len(q))
    return trim(
        [
            (p[i] if i < len(p) else 0)
            + (q[i] if i < len(q) else 0)
            for i in range(n)
        ]
    )


def scale(p, c):
    return trim([c * x for x in p])


def multiply(p, q):
    out = [Fraction(0) for _ in range(len(p) + len(q) - 1)]
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            out[i + j] += x * y
    return trim(out)


def derivative(p):
    if len(p) == 1:
        return (Fraction(0),)
    return trim([i * p[i] for i in range(1, len(p))])


def shift(p, h):
    """Coefficients of p(x+h)."""

    out = [Fraction(0) for _ in p]
    for j, c in enumerate(p):
        for i in range(j + 1):
            out[i] += c * comb(j, i) * h ** (j - i)
    return trim(out)


def difference_quotient(p, h):
    return scale(add(shift(p, h), scale(p, -1)), 1 / h)


def falling_poly(k, h):
    out = (Fraction(1),)
    for j in range(k):
        out = multiply(out, (-j * h, Fraction(1)))
    return out


def umbral(p, h):
    """U_h(x^k)=x^(falling k,h), extended linearly."""

    out = (Fraction(0),)
    for k, c in enumerate(p):
        out = add(out, scale(falling_poly(k, h), c))
    return out


def falling_value(x, k, h):
    out = Fraction(1)
    for j in range(k):
        out *= x - j * h
    return out


def catalan(n):
    return comb(2 * n, n) // (n + 1)


def inclusive_half(n):
    return sum(comb(2 * n, k) for k in range(n + 1))


def strict_half(n):
    return sum(comb(2 * n, k) for k in range(n))


def main():
    hs = (Fraction(1), Fraction(2, 3), Fraction(-3, 5))
    checked_basis = 0
    checked_intertwiners = 0
    for h in hs:
        for k in range(13):
            lhs = difference_quotient(falling_poly(k, h), h)
            rhs = (
                (Fraction(0),)
                if k == 0
                else scale(falling_poly(k - 1, h), k)
            )
            assert lhs == rhs
            checked_basis += 1

        for degree in range(11):
            # A deterministic hostile coefficient vector with mixed signs.
            p = tuple(
                Fraction(((-1) ** k) * (degree + k + 1), k + 1)
                for k in range(degree + 1)
            )
            assert difference_quotient(umbral(p, h), h) == umbral(
                derivative(p), h
            )
            checked_intertwiners += 1

    checked_exponentials = 0
    for h in (Fraction(1), Fraction(2, 3), Fraction(-1, 4)):
        for lam in (
            Fraction(1),
            Fraction(3, 2),
            Fraction(-2, 3),
        ):
            for n in range(13):
                x = n * h
                newton = sum(
                    lam**k * falling_value(x, k, h) / factorial(k)
                    for k in range(n + 1)
                )
                assert newton == (1 + h * lam) ** n
                checked_exponentials += 1

    plus = []
    minus = []
    checked_central = 0
    for n in range(13):
        p = inclusive_half(n)
        m = strict_half(n)
        center = comb(2 * n, n)
        assert p == (4**n + center) // 2
        assert m == (4**n - center) // 2
        assert p + m == 4**n
        assert p - m == center
        if n:
            assert p == 4 * plus[-1] - catalan(n - 1)
            assert m == 4 * minus[-1] + catalan(n - 1)
        plus.append(p)
        minus.append(m)
        checked_central += 1

    # The first coefficients of the two claimed ordinary generating functions.
    # (1-4z)^(-1/2) has coefficient C(2n,n).
    for n, (p, m) in enumerate(zip(plus, minus)):
        assert 2 * p == 4**n + comb(2 * n, n)
        assert 2 * m == 4**n - comb(2 * n, n)

    checked_tournaments = 0
    for vertices in range(1, 10):
        edges = comb(vertices, 2)
        shells = [comb(edges, k) for k in range(edges + 1)]
        assert sum(shells) == 2**edges
        checked_tournaments += 1

    # Hostiles to the two most tempting overstatements.
    x_squared = (Fraction(0), Fraction(0), Fraction(1))
    assert difference_quotient(x_squared, Fraction(1)) == (
        Fraction(1),
        Fraction(2),
    )
    assert difference_quotient(x_squared, Fraction(1)) != derivative(
        x_squared
    )
    assert 2 ** comb(4, 2) == 64
    assert 2**4 == 16

    print("theorem=THM-2412")
    print("arithmetic=Fraction-only")
    print(f"scaled-falling-basis-checks={checked_basis}")
    print(f"umbral-intertwiner-checks={checked_intertwiners}")
    print(f"terminating-exponential-checks={checked_exponentials}")
    print(f"central-layer-checks={checked_central}")
    print("A032443-first=1,3,11,42,163,638")
    print("A000346-shifted-first=1,5,22,93,386,1586")
    print("central-sum=4^n")
    print("central-difference=binomial(2n,n)")
    print("Catalan-leakage=plus:-C_(n-1),minus:+C_(n-1)")
    print(f"tournament-Hamming-shell-checks={checked_tournaments}")
    print("hostile-power-basis=D_1(x^2)=2x+1")
    print("hostile-tournament-exponent=v=4 gives 2^6,not 2^4")
    print("status=PASS")


if __name__ == "__main__":
    main()
