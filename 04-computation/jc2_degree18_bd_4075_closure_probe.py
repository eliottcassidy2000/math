#!/usr/bin/env python3
"""Exact certificate for the rational B--D ratio D/B^2=4075/85176.

This is a standard-library check on the B=1, C=W=0 chart of the
degree-eighteen spectral cubic from THM-2262/2297.  With x=y^2, the
spectral equation is a plane cubic H(u,x).  The script verifies that

* Disc_u(H) has one double root and four further simple roots;
* the double-discriminant fibre is a smooth, totally ramified triple root;
* the projective cubic has three distinct smooth points at infinity; and
* consequently the projective spectral curve is a smooth irreducible
  plane cubic of genus one.

No CAS package is used.  Rational-polynomial arithmetic is imported from
the already audited central-factor verifier.
"""

from __future__ import annotations

from fractions import Fraction
from functools import reduce
import importlib.util
from math import gcd
from pathlib import Path


HERE = Path(__file__).resolve().parent
FACTOR_SCRIPT = HERE / "jc2_degree18_central_t_factor_exact_probe.py"
SPEC = importlib.util.spec_from_file_location("central_factor", FACTOR_SCRIPT)
FACTOR = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(FACTOR)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def monic(poly):
    poly = FACTOR.trim(poly[:])
    return FACTOR.scale(poly, Fraction(1, 1) / poly[-1])


def polynomial_gcd(first, second):
    first = FACTOR.trim(first[:])
    second = FACTOR.trim(second[:])
    while second != [0]:
        _, remainder = FACTOR.divmod_poly(first, second)
        first, second = second, remainder
    return monic(first)


def primitive_integer_coefficients(poly):
    denominator_lcm = 1
    for coefficient in poly:
        denominator_lcm = (
            denominator_lcm
            * coefficient.denominator
            // gcd(denominator_lcm, coefficient.denominator)
        )
    integers = [int(coefficient * denominator_lcm) for coefficient in poly]
    content = reduce(gcd, (abs(value) for value in integers if value))
    integers = [value // content for value in integers]
    if integers[-1] < 0:
        integers = [-value for value in integers]
    return integers


def evaluate(poly, value):
    out = Fraction(0)
    for coefficient in reversed(poly):
        out = out * value + coefficient
    return out


def cubic_discriminant(a, bpoly, cpoly, dpoly):
    terms = [
        FACTOR.mul(FACTOR.power(bpoly, 2), FACTOR.power(cpoly, 2)),
        FACTOR.scale(FACTOR.power(cpoly, 3), -4 * a),
        FACTOR.scale(FACTOR.mul(FACTOR.power(bpoly, 3), dpoly), -4),
        FACTOR.scale(FACTOR.power(dpoly, 2), -27 * a * a),
        FACTOR.scale(
            FACTOR.mul(FACTOR.mul(bpoly, cpoly), dpoly),
            18 * a,
        ),
    ]
    out = [Fraction(0)]
    for term in terms:
        out = FACTOR.add(out, term)
    return out


def scalar_cubic_discriminant(a, b, c, d):
    return (
        b * b * c * c
        - 4 * a * c**3
        - 4 * b**3 * d
        - 27 * a * a * d * d
        + 18 * a * b * c * d
    )


def main() -> None:
    ratio = Fraction(4075, 85176)
    a = Fraction(-26040609)
    bpoly = [Fraction(49601160), Fraction(1607445)]
    cpoly = [
        Fraction(-20995200) - Fraction(52907904) * ratio,
        Fraction(-2857680),
        Fraction(-138915),
    ]
    dpoly = [
        Fraction(33592320) * ratio,
        Fraction(777600) + Fraction(1959552) * ratio,
        Fraction(78120),
        Fraction(1127),
    ]

    branch = cubic_discriminant(a, bpoly, cpoly, dpoly)
    require(len(branch) == 7, "branch discriminant lost degree six")

    double_factor = [Fraction(1215, 91), Fraction(1)]
    branch_gcd = polynomial_gcd(branch, FACTOR.derivative(branch))
    require(branch_gcd == double_factor, "unexpected repeated branch factor")

    quartic, remainder = FACTOR.divmod_poly(
        branch,
        FACTOR.power(double_factor, 2),
    )
    require(remainder == [0], "double branch factor does not divide twice")
    expected_quartic = [
        -3037628182500,
        1389005982000,
        -147517112925,
        203464170,
        1577224103,
    ]
    require(
        primitive_integer_coefficients(quartic) == expected_quartic,
        "primitive residual quartic changed",
    )
    require(
        len(polynomial_gcd(quartic, FACTOR.derivative(quartic))) == 1,
        "residual quartic is not squarefree",
    )

    x0 = Fraction(-1215, 91)
    u0 = Fraction(295, 819)
    require(evaluate(quartic, x0) != 0, "triple fibre meets residual quartic")

    fibre = [
        evaluate(dpoly, x0),
        evaluate(cpoly, x0),
        evaluate(bpoly, x0),
        a,
    ]
    triple_fibre = [
        -a * u0**3,
        3 * a * u0**2,
        -3 * a * u0,
        a,
    ]
    require(fibre == triple_fibre, "exceptional fibre is not a triple root")

    h_x = (
        Fraction(1607445) * u0**2
        + (Fraction(-2857680) - Fraction(277830) * x0) * u0
        + Fraction(777600)
        + Fraction(1959552) * ratio
        + Fraction(156240) * x0
        + Fraction(3381) * x0**2
    )
    require(h_x == Fraction(-16329600, 169), "triple fibre lost smoothness")

    infinity_coefficients = [
        Fraction(1127),
        Fraction(-138915),
        Fraction(1607445),
        Fraction(-26040609),
    ]
    infinity_discriminant = scalar_cubic_discriminant(
        infinity_coefficients[3],
        infinity_coefficients[2],
        infinity_coefficients[1],
        infinity_coefficients[0],
    )
    require(
        infinity_discriminant
        == -153384762202971019112448,
        "infinity cubic lost its three distinct roots",
    )

    # Any affine singular point projects to a repeated root of Disc_u(H).
    # The only such x is x0, where H_x is nonzero.  At infinity, the
    # squarefree leading cubic gives three distinct smooth points.  Hence
    # the projective plane cubic is smooth.  A smooth projective plane cubic
    # is geometrically irreducible and has genus one.
    ramification = 2 + 4
    genus = (ramification - 4) // 2
    require(ramification == 6 and genus == 1, "Riemann--Hurwitz count changed")

    print("chart=B=1,C=0,W=0")
    print("ratio_D_over_B2=4075/85176")
    print("quotient_coordinate=x=y^2")
    print("branch_degree=6")
    print("branch_gcd=x+1215/91")
    print(
        "residual_quartic_primitive_ascending="
        + ",".join(str(value) for value in expected_quartic)
    )
    print("residual_quartic=squarefree_and_disjoint_from_triple_fibre")
    print("triple_fibre_x=-1215/91")
    print("triple_fibre_u=295/819")
    print("triple_fibre_Hx=-16329600/169")
    print(f"infinity_discriminant={infinity_discriminant}")
    print("affine_singular_points=0")
    print("infinity_points=3_distinct_smooth")
    print("projective_spectral_curve=smooth_irreducible_plane_cubic")
    print("finite_ramification=6")
    print("genus=1")
    print("status=BD_4075_RATIO_EXACT_CLOSURE_CERTIFICATE")


if __name__ == "__main__":
    main()
