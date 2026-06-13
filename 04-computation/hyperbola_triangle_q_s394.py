#!/usr/bin/env python3
"""
hyperbola_triangle_q_s394.py

codex-2026-05-31 S394

Explore the area under

    p(x) = x       for 0 <= x <= 1
         = 1/x     for 1 <= x <= 2,

as a logarithmic/hyperbolic renormalization problem.  The main object is the
one-sided reciprocal counterterm

    q(x) = x - 1/x    for 0 < x <= 1,
         = 0          for 1 <= x <= 2,

so that p(x) = 1/x + q(x) after the singular hyperbola is interpreted with a
cutoff at 0.
"""

from __future__ import annotations

from fractions import Fraction
from math import atan, degrees, log, radians, sqrt, tan


def f6(x: float) -> str:
    return f"{x:.6f}"


def q_bare(x: float) -> float:
    return x - 1.0 / x


def q_gate(x: float) -> float:
    return q_bare(x) if x <= 1.0 else 0.0


def p_curve(x: float) -> float:
    return x if x <= 1.0 else 1.0 / x


def print_header(title: str) -> None:
    print()
    print(title)
    print("-" * len(title))


def triangle_sanity_checks() -> None:
    print_header("1. Triangle sanity checks")
    total = 0.5 + log(2)
    print(f"area under p(x)=min(x,1/x) on [0,2] = 1/2 + ln 2 = {f6(total)}")
    print("To leave exactly ln 2, the carved triangle must have area 1/2.")
    print()

    print("Line from the origin reaching height 1 at angle theta:")
    print("  base = cot(theta), triangle area = 1/(2 tan(theta))")
    for theta in (30, 45, 60):
        area = 1.0 / (2.0 * tan(radians(theta)))
        print(
            f"  theta={theta:2d} deg: area={f6(area)}, "
            f"remainder if total fixed={f6(total - area)}"
        )
    print("  Only theta=45 deg gives the exact 1/2 triangle while gluing to (1,1).")
    print()

    print("If a triangle has horizontal span 2 and area 1/2, its height is 1/2.")
    print(f"  The corresponding base angle is atan((1/2)/2) = {f6(degrees(atan(0.25)))} deg.")
    print("If a 30-60-90 triangle has hypotenuse/horizontal span 2, its area is sqrt(3)/2.")
    print(f"  sqrt(3)/2 = {f6(sqrt(3)/2)}, not 1/2.")


def q_cutoff_table() -> None:
    print_header("2. q(x)+1/x as a cutoff identity")
    print("For epsilon>0:")
    print("  int_epsilon^2 1/x dx = ln(2/epsilon)")
    print("  int_epsilon^1 q(x) dx = 1/2 + ln(epsilon) - epsilon^2/2")
    print("  sum = 1/2 + ln 2 - epsilon^2/2 -> 1/2 + ln 2")
    print()
    print("epsilon     hyperbola       q-counterterm   finite total")
    for n in (10, 100, 1000, 10000, 100000):
        eps = 1.0 / n
        hyper = log(2.0 / eps)
        qpart = 0.5 + log(eps) - eps * eps / 2.0
        print(f"1/{n:<7d} {hyper:14.9f} {qpart:15.9f} {hyper + qpart:14.9f}")


def log_coordinate_view() -> None:
    print_header("3. Log-coordinate view")
    print("Set x=e^u.  Then dx/x = du, so 1/x is flat log-measure.")
    print("The patched curve has area density p(e^u)e^u:")
    print("  e^(2u) for u<=0, and 1 for 0<=u<=ln 2.")
    print("Hence:")
    print("  int_{-infty}^0 e^(2u) du = 1/2")
    print("  int_0^ln2 1 du = ln 2")
    print("The triangle is an exponential soft cutoff of the negative log-ray.")
    print()
    print("In the same coordinate, q(e^u)e^u = e^(2u)-1 for u<=0, and 0 for u>=0.")
    print("So q deletes the infinite negative log-ray and replaces it by a finite tail.")


def q_algebra() -> None:
    print_header("4. Algebra of the bare reciprocal defect")
    print("Bare q0(x)=x-1/x is anti-invariant under reciprocal inversion:")
    for x in (0.25, 0.5, 2.0, 4.0):
        print(
            f"  x={x:<4g}: q0(x)={q_bare(x):>9.6f}, "
            f"q0(1/x)={q_bare(1/x):>9.6f}, sum={q_bare(x)+q_bare(1/x):>9.2e}"
        )
    print()
    print("With r(x)=x+1/x:")
    print("  r^2 - q0^2 = 4")
    print("  D q0 = r and D r = q0 for D=x*d/dx")
    print("  so D^2 q0 = q0. In log-coordinate, q0(e^u)=2 sinh(u).")
    print()
    print("Mellin moments of q0 on (0,1):")
    print("  int_0^1 x^n (x-1/x) dx = 1/(n+2) - 1/n = -2/(n(n+2))")
    for n in range(1, 8):
        moment = Fraction(1, n + 2) - Fraction(1, n)
        print(f"  n={n}: {moment} = {float(moment): .9f}")
    print("The meromorphic Mellin transform is M_q(s)=1/(s+1)-1/(s-1)=-2/(s^2-1).")


def cutoff_family() -> None:
    print_header("5. A family of possible q-cutoffs")
    print("Replace the left branch by x^a while keeping p(1)=1.")
    print("Then the left mass is int_0^1 x^a dx = 1/(a+1).")
    print("The exact 1/2 mass is forced by a=1.")
    print()
    print("a       left mass       scale needed for mass 1/2       jump at x=1")
    for a in (0.25, 0.5, 1.0, 2.0, 3.0, 5.0):
        mass = 1.0 / (a + 1.0)
        scale = (a + 1.0) / 2.0
        print(f"{a:<5.2f} {mass:14.9f} {scale:30.9f} {scale - 1.0:15.9f}")
    print()
    print("So the isosceles triangle is unique if we demand both:")
    print("  (i) continuity with the hyperbola at x=1, and")
    print("  (ii) triangle/cutoff mass exactly 1/2.")
    print("Other triangles are possible only by paying a q-defect: a jump, a slope change, or an area debt.")


def q_tangents() -> None:
    print_header("6. Tangents suggested by q")
    tangents = [
        "q is a renormalization counterterm, not a positive area.",
        "q0=x-1/x is the anti-reciprocal coordinate paired with r=x+1/x.",
        "The cutoff p(x)=1/x+q(x) is a one-sided Heaviside gate at the fixed point x=1.",
        "The Mellin kernel -2/(s^2-1) makes q a two-pole residue object.",
        "In log-space, the triangle is the exponential tail e^(2u), while 1/x is flat measure.",
        "A literal 30-60-90 triangle cannot keep the ln 2 remainder without an extra q-defect.",
    ]
    for i, tangent in enumerate(tangents, start=1):
        print(f"  {i}. {tangent}")


def main() -> None:
    print("Hyperbola/triangle q exploration (codex-2026-05-31 S394)")
    print("Object: p(x)=x on [0,1], p(x)=1/x on [1,2].")
    triangle_sanity_checks()
    q_cutoff_table()
    log_coordinate_view()
    q_algebra()
    cutoff_family()
    q_tangents()


if __name__ == "__main__":
    main()
