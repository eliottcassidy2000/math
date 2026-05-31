#!/usr/bin/env python3
"""
log_triangle_q_decomposition_s396.py

codex-2026-05-31 S396

Exact notes for the two-piece area

    f(x) = x      on 0 <= x <= 1
         = 1/x    on 1 <= x <= 2.

The prompt asks for two decompositions leaving an ln(2) remainder and for
an exploration of the correction term q(x).  The exact algebra says:

* the first decomposition is the visible triangle area 1/2 plus ln(2);
* the second decomposition, if it uses a right triangle with base [0,2]
  and still leaves ln(2), must subtract the line x/4, not a 30-60-90 line;
* after subtracting x/4, the remainder is H(x-1)/x + q(x), where q has
  total integral zero.
"""

from __future__ import annotations

from fractions import Fraction
from math import atan, degrees, log, sqrt


def fmt(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def q_integral(a: Fraction) -> Fraction:
    """Integral of q_a over [0,2].

    For triangle line a*x:
        q_a(x) = (1-a)x   on [0,1]
              = -a x      on [1,2].
    """

    return Fraction(1, 2) - 2 * a


def q_moment(a: Fraction, n: int) -> Fraction:
    """Integral x^n q_a(x) dx over [0,2]."""

    return (1 - a * (2 ** (n + 2))) / Fraction(n + 2, 1)


def qstar_moment(n: int) -> Fraction:
    return q_moment(Fraction(1, 4), n)


def qstar_potential_value(x: Fraction) -> Fraction:
    """Antiderivative Q with Q(0)=Q(2)=0 for the log-preserving q."""

    if x <= 1:
        return Fraction(3, 8) * x * x
    return Fraction(4, 8) - Fraction(1, 8) * x * x


def print_area_decompositions() -> None:
    print("1. Exact area ledger")
    triangle_01 = Fraction(1, 2)
    print("   f(x)=x on [0,1], f(x)=1/x on [1,2]")
    print(f"   area = {fmt(triangle_01)} + ln(2) = {0.5 + log(2):.12f}")
    print()
    print("   Subtract a right-triangle ramp y=a*x over [0,2].")
    print("   Its area is integral_0^2 a*x dx = 2a.")
    print("   To leave exactly ln(2), one needs 2a=1/2, hence a=1/4.")
    print(
        f"   The angle is atan(1/4) = {degrees(atan(Fraction(1, 4))):.6f} degrees,"
    )
    print("   not 30 degrees.")
    print()


def print_306090_obstruction() -> None:
    print("2. 30-60-90 obstruction")
    print(
        "   A literal Euclidean 30-60-90 triangle leaving ln(2) would need area 1/2."
    )
    hyp = 2 / (3 ** 0.25)
    print(
        "   Such a triangle has hypotenuse 2/3^(1/4) "
        f"= {hyp:.6f}, so its horizontal span is < 2."
    )
    print("   Conversely, common span-2 interpretations have the wrong area:")
    print(f"     long leg horizontal 2: area = 2/sqrt(3) = {2/sqrt(3):.12f}")
    print(f"     short leg horizontal 2: area = 2*sqrt(3) = {2*sqrt(3):.12f}")
    print(f"     hypotenuse length 2:  area = sqrt(3)/2 = {sqrt(3)/2:.12f}")
    print(
        "   So the log-preserving second carve is a 4:1 ramp triangle.  If a\n"
        "   30-60-90 triangle is intended, it must live in a transformed metric\n"
        "   or a different normalization."
    )
    print()


def print_q_family() -> None:
    print("3. The q_a family")
    print("   After subtracting y=a*x, the residual can be written")
    print("       H(x-1)/x + q_a(x)")
    print("   where")
    print("       q_a(x)=(1-a)x on [0,1],    q_a(x)=-a*x on [1,2].")
    print("   The q_a area is 1/2 - 2a, so q_a is area-neutral iff a=1/4.")
    print()
    print("   a          triangle area   integral q_a   residual area")
    for a in [
        Fraction(0),
        Fraction(1, 6),
        Fraction(1, 4),
        Fraction(1, 3),
        Fraction(1, 2),
    ]:
        residual_offset = q_integral(a)
        sign = "+" if residual_offset >= 0 else "-"
        offset_abs = abs(residual_offset)
        print(
            f"   {fmt(a):<10} {fmt(2*a):<15} {fmt(residual_offset):<13} "
            f"ln(2) {sign} {fmt(offset_abs)}"
        )
    print()


def print_qstar_structure() -> None:
    print("4. The log-preserving q = q_{1/4}")
    print("   q(x)=3x/4 on [0,1], q(x)=-x/4 on [1,2].")
    print("   It has zero area but nonzero transport moments.")
    print()
    print("   moments M_n = integral_0^2 x^n q(x) dx = (1-2^n)/(n+2)")
    for n in range(0, 8):
        print(f"     M_{n} = {fmt(qstar_moment(n))}")
    print()
    print("   A primitive Q with Q(0)=Q(2)=0 is")
    print("       Q(x)=3x^2/8      on [0,1]")
    print("       Q(x)=(4-x^2)/8   on [1,2]")
    print("   Sample values:")
    for x in [Fraction(0), Fraction(1, 2), Fraction(1), Fraction(3, 2), Fraction(2)]:
        print(f"     Q({fmt(x)}) = {fmt(qstar_potential_value(x))}")
    print(
        "   This makes q a flux term: it moves area inside the decomposition\n"
        "   but contributes no net area, leaving the logarithmic piece untouched."
    )
    print()


def print_mellin_tangent() -> None:
    print("5. Mellin/log tangent")
    print("   The Mellin transform of q_a is")
    print("       M_a(s) = integral_0^2 x^(s-1) q_a(x) dx")
    print("              = (1 - a*2^(s+1))/(s+1).")
    print("   For a=1/4:")
    print("       M(s) = (1 - 2^(s-1))/(s+1), so M(1)=0.")
    print("   The first log moment is the derivative at the zero:")
    print("       M'(1)= integral q(x) log(x) dx = -ln(2)/2.")
    print(f"   Numerically: {-log(2)/2:.12f}")
    print(
        "   So q is not noise; it is the zero-mass distribution whose first\n"
        "   Mellin variation remembers exactly half a bit."
    )
    print()


def print_synthesis() -> None:
    print("6. Synthesis")
    print(
        "   The robust identity is not '30-60-90 + ln(2)' in Euclidean area.\n"
        "   It is 'base-[0,2] ramp of slope 1/4 + a zero-area q transport + ln(2)'."
    )
    print(
        "   This suggests a useful general pattern: logarithms survive when the\n"
        "   polynomial/triangular part is adjusted until its signed correction\n"
        "   has zero mass.  The correction q then stores scale information in\n"
        "   its Mellin derivatives rather than in its area."
    )


def main() -> None:
    print("Log triangle q-decomposition probe (codex-2026-05-31 S396)")
    print()
    print_area_decompositions()
    print_306090_obstruction()
    print_q_family()
    print_qstar_structure()
    print_mellin_tangent()
    print_synthesis()


if __name__ == "__main__":
    main()
