#!/usr/bin/env python3
"""Exact certificate that the sharp 5x5 attachment point is non-torsion.

This is deliberately dependency-free.  It reconstructs the rational
constants, verifies two explicit residue places of the degree-six coordinate
field, computes the reduced point orders by the chord-and-tangent law, and
checks the incompatible prime-to-residue-characteristic order invoice.
"""

from fractions import Fraction
from hashlib import sha256
from math import gcd


TRANSCRIPT = []


def gate(label: str, condition: bool, detail: str) -> None:
    if not condition:
        raise RuntimeError(f"FAIL {label}: {detail}")
    line = f"PASS {label}: {detail}"
    TRANSCRIPT.append(line)
    print(line)


def frac_mod(q: Fraction, p: int) -> int:
    return q.numerator * pow(q.denominator, -1, p) % p


def add(P, Q, p: int):
    """Group law on E/F_p: y^2=x^3+1; None denotes O."""
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2 and (y1 + y2) % p == 0:
        return None
    if P == Q:
        slope = 3 * x1 * x1 * pow(2 * y1, -1, p) % p
    else:
        slope = (y2 - y1) * pow(x2 - x1, -1, p) % p
    x3 = (slope * slope - x1 - x2) % p
    y3 = (slope * (x1 - x3) - y1) % p
    return x3, y3


def mul(n: int, P, p: int):
    R = None
    Q = P
    while n:
        if n & 1:
            R = add(R, Q, p)
        Q = add(Q, Q, p)
        n >>= 1
    return R


def exact_order(P, p: int, claimed: int) -> bool:
    if mul(claimed, P, p) is not None:
        return False
    for ell in (2, 3, 5, 7, 11, 13, 17, 19):
        if claimed % ell == 0 and mul(claimed // ell, P, p) is None:
            return False
    return True


def valuation(n: int, ell: int) -> int:
    e = 0
    while n % ell == 0:
        n //= ell
        e += 1
    return e


def main() -> None:
    # Dimensionless coefficients: the omitted powers are
    # c30*A5^-3, c40*A5^-4, c02*A5^-3.
    c30 = Fraction(2752, 135)
    c40 = Fraction(-512, 9)
    # 105*A5^4*c40_actual + 90*A5^3*c02_actual + 11392 = 0.
    c02 = -(105 * c40 + 11392) / 90
    gate("constant-relation", 105 * c40 + 90 * c02 + 11392 == 0,
         f"c40={c40}, c02={c02}")
    gate("sharp-hostile-c02", c02 == Fraction(-8128, 135),
         "[y^2]Rtilde=-8128/(135*A5^3)")

    # gamma=-a^3/2 and A5=a^5 turn both weight-six coefficients
    # into a^-12 times the following rational scalars.
    epsilon = -c30 / 2
    kappa = -c02 / 2
    total = epsilon + kappa
    gate("raw-epsilon", epsilon == Fraction(-1376, 135),
         "[p^3]R=(-1376/135)*a^-12")
    gate("raw-kappa", kappa == Fraction(4064, 135),
         "[y^2]R=(4064/135)*a^-12")
    gate("nonresonance", total == Fraction(2688, 135) and total != 0,
         "epsilon+kappa=(2688/135)*a^-12")

    x_cube = -epsilon / total
    y_square = kappa / total
    gate("normalized-X", x_cube == Fraction(43, 84), "X^3=43/84")
    gate("normalized-Y", y_square == Fraction(127, 84), "Y^2=127/84")
    gate("curve-equation", y_square == x_cube + 1, "Y^2=X^3+1")

    # The valuations v_43(43/84)=1 and v_127(127/84)=1 prove that the
    # cubic and quadratic are irreducible.  Their coprime degrees force
    # trivial field intersection (the quadratic extension is Galois).
    vx43 = valuation(abs(x_cube.numerator), 43) - valuation(x_cube.denominator, 43)
    vy127 = valuation(abs(y_square.numerator), 127) - valuation(y_square.denominator, 127)
    gate("coordinate-field", vx43 == 1 and vy127 == 1 and gcd(3, 2) == 1,
         "[Q(alpha,beta):Q]=3*2=6 by the 43- and 127-valuations")

    places = [
        (11, (9, 2), 12),
        (17, (2, 3), 6),
    ]
    for p, P, claimed_order in places:
        x, y = P
        gate(f"place-{p}-alpha", x**3 % p == frac_mod(x_cube, p),
             f"alpha maps to {x}")
        gate(f"place-{p}-beta", y*y % p == frac_mod(y_square, p),
             f"beta maps to {y}")
        gate(f"place-{p}-curve", (y*y - x**3 - 1) % p == 0,
             f"Pbar={P} lies on E/F_{p}")
        gate(f"place-{p}-good", (-432) % p != 0,
             f"Delta(E)=-432 is a {p}-adic unit")
        gate(f"place-{p}-order", exact_order(P, p, claimed_order),
             f"ord(Pbar)={claimed_order}")

    P11 = places[0][1]
    gate("p11-multiples",
         (mul(2, P11, 11), mul(3, P11, 11), mul(4, P11, 11),
          mul(6, P11, 11)) == ((2, 8), (5, 4), (0, 10), (10, 0)),
         "2P=(2,8), 3P=(5,4), 4P=(0,10), 6P=(10,0)")
    P17 = places[1][1]
    gate("p17-multiples",
         (mul(2, P17, 17), mul(3, P17, 17)) == ((0, 1), (16, 0)),
         "2P=(0,1), 3P=(16,0)")

    # Good reduction has no prime-to-p torsion in its kernel.  Therefore a
    # hypothetical global order N would be 12*11^a and also 6*17^b.  Since
    # neither residue characteristic is 2, their 2-adic valuations must agree.
    gate("all-order-contradiction",
         2 not in (11, 17) and valuation(12, 2) != valuation(6, 2),
         "v2(N)=v2(12)=2 at 11 but v2(N)=v2(6)=1 at 17")

    semantic = sha256("\n".join(TRANSCRIPT).encode("ascii")).hexdigest()
    print(f"SEMANTIC_SHA256 {semantic}")
    print("ALL ATTACHMENT NON-TORSION CHECKS PASSED")


if __name__ == "__main__":
    main()
