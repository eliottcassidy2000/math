#!/usr/bin/env python3
"""Exact rational referee for THM-2164.

The analytic input is THM-2085's ordered Selberg interval pair.  This
companion verifies the elementary sinc majorants and the two finite
whole-subtorus packet estimates used for relative rank harvesting.
"""

from fractions import Fraction


ALPHA = Fraction(6, 7)
DIMENSION = 13


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def cosine_quartic_upper(x: Fraction) -> Fraction:
    """The alternating-Taylor upper polynomial 1-x^2/2+x^4/24."""
    return 1 - x * x / 2 + x**4 / 24


def residue_sinc_bound(index: int) -> Fraction:
    """Rational upper bound for |hat(1_[1/14,13/14])(index)|."""
    require(index >= 1, "positive Fourier index required")
    residue = index % 7
    if residue == 0:
        return Fraction(0)
    if residue in (1, 6):
        return Fraction(1, 7 * index)
    if residue in (2, 5):
        return Fraction(1, 4 * index)
    return Fraction(5, 16 * index)


def packet_margin(
    height: int, support: int, residue_sensitive: bool
) -> Fraction:
    """Lower bound for the Selberg packet on a primitive relation line."""
    require(height >= 1, "positive degree required")
    require(3 <= support <= DIMENSION, "support range")
    epsilon = Fraction(1, height + 1)
    upper_constant = ALPHA + epsilon
    defect_constant = 2 * epsilon
    constant_term = upper_constant**12 * (ALPHA - 25 * epsilon)
    tail = Fraction(0)
    for index in range(1, height + 1):
        if residue_sensitive:
            interval_coefficient = residue_sinc_bound(index)
        else:
            interval_coefficient = Fraction(5, 16 * index)
        upper_coefficient = interval_coefficient + epsilon
        coefficient_bound = (
            upper_coefficient**support
            * upper_constant ** (DIMENSION - support)
            + support
            * defect_constant
            * upper_coefficient ** (support - 1)
            * upper_constant ** (DIMENSION - support)
            + (DIMENSION - support)
            * defect_constant
            * upper_coefficient**support
            * upper_constant ** (DIMENSION - support - 1)
        )
        tail += coefficient_bound
    return constant_term - 2 * tail


def main() -> None:
    # Elementary wrappers for the three residue-dependent sinc bounds.
    # pi>157/50.  For the +/-2 residues,
    # sin(2pi/7)=cos(3pi/14) and 3pi/14>471/700.
    require(
        cosine_quartic_upper(Fraction(471, 700)) < Fraction(157, 200),
        "+/-2 residue cosine bound",
    )
    # For the +/-3 residues, pi/14>3/14 and
    # cos(pi/14)<49/50<5pi/16.
    require(
        cosine_quartic_upper(Fraction(3, 14)) < Fraction(49, 50),
        "+/-3 residue cosine bound",
    )
    require(Fraction(49, 50) < Fraction(157, 160), "5pi/16 lower bound")

    signed_margins_34 = {
        support: packet_margin(34, support, True)
        for support in range(3, DIMENSION + 1)
    }
    signed_min_support = min(signed_margins_34, key=signed_margins_34.get)
    signed_minimum = signed_margins_34[signed_min_support]
    require(signed_min_support == 3, "signed packet minimizing support")
    require(signed_minimum > Fraction(1, 212), "height-34 signed packet")
    require(packet_margin(33, 3, True) < 0, "height-33 signed boundary")

    uniform_margins_43 = {
        support: packet_margin(43, support, False)
        for support in range(3, DIMENSION + 1)
    }
    uniform_min_support = min(uniform_margins_43, key=uniform_margins_43.get)
    uniform_minimum = uniform_margins_43[uniform_min_support]
    require(uniform_min_support == 3, "uniform packet minimizing support")
    require(uniform_minimum > Fraction(1, 634), "height-43 uniform packet")
    require(packet_margin(42, 3, False) < 0, "height-42 uniform boundary")

    print("THM-2164 exact relative-packet rank-harvesting referee")
    print("interval sinc table: residues 0,+/-1,+/-2,+/-3 use 0,1/(7l),1/(4l),5/(16l)")
    print("elementary pi/cosine wrappers for the residue table passed")
    print(
        "signed-unit packet H=34: "
        f"minimum support={signed_min_support}, margin={signed_minimum} > 1/212"
    )
    print("signed-unit H=33 support-3 certificate is negative")
    print(
        "arbitrary-relation packet H=43: "
        f"minimum support={uniform_min_support}, margin={uniform_minimum} > 1/634"
    )
    print("arbitrary-relation H=42 support-3 certificate is negative")
    print("ALL EXACT ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
