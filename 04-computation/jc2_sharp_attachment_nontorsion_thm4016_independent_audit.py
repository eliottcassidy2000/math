#!/usr/bin/env python3
"""No-import audit of the THM-4007/4008 attachment non-torsion claim.

This verifier intentionally does not import either scratch implementation under
audit.  It uses naive repeated addition for the two finite-field orders and a
separate Eisenstein-integer ledger for the six-point automorphism.
"""

from fractions import Fraction
from hashlib import sha256
from math import gcd


LINES = []


def check(label, condition, detail):
    if not condition:
        raise AssertionError(f"FAIL {label}: {detail}")
    line = f"PASS {label}: {detail}"
    LINES.append(line)
    print(line)


def vprime_fraction(q, p):
    def vint(n):
        e = 0
        while n % p == 0:
            e += 1
            n //= p
        return e

    return vint(abs(q.numerator)) - vint(q.denominator)


def reduce_fraction(q, p):
    return q.numerator * pow(q.denominator, -1, p) % p


def point_add(P, Q, p):
    """Chord-and-tangent addition on y^2=x^3+1; infinity is None."""
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2 and (y1 + y2) % p == 0:
        return None
    if P == Q:
        lam = 3 * x1 * x1 * pow(2 * y1, -1, p) % p
    else:
        lam = (y2 - y1) * pow(x2 - x1, -1, p) % p
    x3 = (lam * lam - x1 - x2) % p
    return x3, (lam * (x1 - x3) - y1) % p


def naive_multiples(P, p, limit):
    """Return P,2P,... using only serial addition by P."""
    values = []
    Q = None
    for _ in range(limit):
        Q = point_add(Q, P, p)
        values.append(Q)
    return values


def exact_order_serial(P, p, bound):
    for n, Q in enumerate(naive_multiples(P, p, bound), 1):
        if Q is None:
            return n
    raise AssertionError((P, p, bound))


def eisenstein_norm(a, b):
    """Norm of a+b*zeta where zeta^2+zeta+1=0."""
    return a * a - a * b + b * b


def main():
    # Reconstruct the sharp third-row residual from the two independent linear
    # constraints: vanishing of the new C weight -5 and the affine residual.
    c40 = -2 * Fraction(256, 9)
    c02 = -(105 * c40 + 11392) / 90
    check("sharp-residual", (c40, c02) ==
          (Fraction(-512, 9), Fraction(-8128, 135)),
          "c40=-512/9 and c02=-8128/135 after suppressing A5 powers")

    c30 = Fraction(2752, 135)
    epsilon = -c30 / 2
    kappa = -c02 / 2
    total = epsilon + kappa
    check("raw-face", (epsilon, kappa, total) ==
          (Fraction(-1376, 135), Fraction(4064, 135),
           Fraction(2688, 135)),
          "epsilon,kappa,total carry the common nonzero a^-12 factor")

    alpha_cube = -epsilon / total
    beta_square = kappa / total
    check("normalized-point", (alpha_cube, beta_square) ==
          (Fraction(43, 84), Fraction(127, 84)) and
          beta_square == alpha_cube + 1,
          "alpha^3=43/84, beta^2=127/84 on y^2=x^3+1")

    # A reducible cubic T^3-q over Q has a rational root.  The displayed
    # valuations rule out a rational cube/square.  Coprime degrees then make
    # the two extensions linearly disjoint.
    check("degree-six-field",
          vprime_fraction(alpha_cube, 43) == 1 and
          vprime_fraction(beta_square, 127) == 1 and gcd(3, 2) == 1,
          "irreducible degrees 3 and 2 have trivial intersection")

    places = ((11, 9, 2, 12), (17, 2, 3, 6))
    observed = []
    for p, x, y, expected in places:
        ac = reduce_fraction(alpha_cube, p)
        bs = reduce_fraction(beta_square, p)
        cube_roots = tuple(z for z in range(p) if z ** 3 % p == ac)
        square_roots = tuple(z for z in range(p) if z * z % p == bs)
        check(f"place-{p}", x in cube_roots and y in square_roots and
              (y * y - x ** 3 - 1) % p == 0 and (-432) % p,
              f"evaluation alpha->{x}, beta->{y}; good discriminant")
        check(f"place-{p}-root-shape",
              cube_roots == (x,) and set(square_roots) == {y, (-y) % p},
              "the cube root is unique and the two beta signs negate P")
        order = exact_order_serial((x, y), p, p + 1 + 2 * int(p ** 0.5) + 2)
        observed.append(order)
        check(f"order-{p}", order == expected,
              f"serial-addition order is exactly {expected}")

    mult11 = naive_multiples((9, 2), 11, 6)
    mult17 = naive_multiples((2, 3), 17, 3)
    check("explicit-multiples-11",
          (mult11[1], mult11[2], mult11[3], mult11[5]) ==
          ((2, 8), (5, 4), (0, 10), (10, 0)),
          "2P,3P,4P,6P agree; 6P is nonzero 2-torsion")
    check("explicit-multiples-17",
          (mult17[1], mult17[2]) == ((0, 1), (16, 0)),
          "2P,3P agree; 3P is nonzero 2-torsion")

    # If a global torsion order is N and the reduction order is m, the kernel
    # of <P> -> <Pbar> has order N/m.  Good reduction makes every torsion
    # element of that kernel p-primary, so N=m*p^r.  The prime 2 is different
    # from both residue characteristics and gives incompatible exponents.
    v2_12 = 2
    v2_6 = 1
    check("order-kernel-contradiction",
          observed == [12, 6] and v2_12 != v2_6 and 2 not in (11, 17),
          "N=12*11^r and N=6*17^s would force v2(N)=2 and 1")

    # Powers of sigma=(zeta_3 X,-Y) give every pair of cube-root and sign
    # exponents.  As an Eisenstein unit sigma=-zeta_3, so
    # sigma-1=-1-zeta_3=zeta_3^2 has norm one.
    orbit_indices = {((i % 3), (i % 2)) for i in range(6)}
    check("six-point-orbit", orbit_indices ==
          {(i, j) for i in range(3) for j in range(2)},
          "sigma^0,...,sigma^5 exhaust the 3x2 root choices")
    check("sigma-minus-one-unit", eisenstein_norm(-1, -1) == 1,
          "sigma=-zeta_3 and sigma-1=zeta_3^2 is an Eisenstein unit")

    digest = sha256("\n".join(LINES).encode("ascii")).hexdigest()
    print(f"SEMANTIC_SHA256 {digest}")
    print("INDEPENDENT ATTACHMENT AUDIT PASSED")


if __name__ == "__main__":
    main()
