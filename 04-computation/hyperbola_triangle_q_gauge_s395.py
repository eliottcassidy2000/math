#!/usr/bin/env python3
"""
hyperbola_triangle_q_gauge_s395.py

codex-2026-05-31 S395

Exploration of the user's hyperbola/triangle decomposition:

    area under y=x on [0,1] plus area under y=1/x on [1,2]
    = 1/2 + log(2).

There are two q-objects worth separating.

1. A triangle-density gauge q_L.  Any integrable q with total mass 1/2
   gives the same finite part

       FP int_0^2 (1/x + q(x)) dx = log(2) + 1/2.

   The equal-area line family is q_L(x)=x/L^2 on [0,L].

2. A pointwise counterterm q_ct(x)=x-1/x on (0,1).  Then the original splice
   is literally 1/x + q_ct below 1 and 1/x above 1, with logarithmic
   divergences cancelling.

The second object is reciprocal-odd: q_ct(1/x)=-q_ct(x).  In log coordinate
x=e^t it is 2*sinh(t), while 1/x dx becomes dt.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import atan, degrees, log, pi, sqrt, tan


TARGET_B = 2.0
TARGET_AREA = 0.5 + log(TARGET_B)


@dataclass(frozen=True)
class TriangleGauge:
    label: str
    L: float
    slope: float
    height: float
    angle_deg: float
    log_cut: float
    area: float
    finite_part_total: float


def triangle_from_base(label: str, L: float) -> TriangleGauge:
    slope = 1.0 / (L * L)
    height = slope * L
    return TriangleGauge(
        label=label,
        L=L,
        slope=slope,
        height=height,
        angle_deg=degrees(atan(slope)),
        log_cut=log(L),
        area=0.5,
        finite_part_total=0.5 + log(TARGET_B),
    )


def triangle_from_angle(label: str, angle_deg: float) -> TriangleGauge:
    slope = tan(angle_deg * pi / 180.0)
    L = 1.0 / sqrt(slope)
    height = slope * L
    return TriangleGauge(
        label=label,
        L=L,
        slope=slope,
        height=height,
        angle_deg=angle_deg,
        log_cut=log(L),
        area=0.5,
        finite_part_total=0.5 + log(TARGET_B),
    )


def q_ct(x: float) -> float:
    return x - 1.0 / x


def cutoff_counterterm_sum(epsilon: float) -> tuple[float, float, float]:
    hyperbola = log(TARGET_B) - log(epsilon)
    counterterm = (1.0 - epsilon * epsilon) / 2.0 + log(epsilon)
    return hyperbola, counterterm, hyperbola + counterterm


def forced_base_angle(base: float, angle_deg: float) -> tuple[float, float, float, float]:
    slope = tan(angle_deg * pi / 180.0)
    area = slope * base * base / 2.0
    finite_total = log(TARGET_B) + area
    normalized_slope = 1.0 / (base * base)
    y_scale_to_look_like_angle = slope / normalized_slope
    return area, finite_total, normalized_slope, y_scale_to_look_like_angle


def mellin_q_L(L: float, s: float) -> float:
    # int_0^L x^s * (x/L^2) dx = L^s/(s+2), valid for Re(s)>-2.
    return (L**s) / (s + 2.0)


def print_header() -> None:
    print("HYPERBOLA + TRIANGLE Q EXPLORATION")
    print("=" * 78)
    print(f"Original area int_0^1 x dx + int_1^2 dx/x = {TARGET_AREA:.12f}")
    print(f"  1/2 = {0.5:.12f}")
    print(f"  log(2) = {log(2.0):.12f}")
    print()


def print_triangle_gauges() -> None:
    print("1. Equal-area triangle gauges q_L")
    print("-" * 78)
    print("Define q_L(x)=x/L^2 on 0<=x<=L, and 0 otherwise.")
    print("Then int q_L dx = 1/2 and FP int_0^2 (1/x+q_L) dx = log(2)+1/2.")
    print()
    gauges = [
        triangle_from_base("45 degree isosceles, L=1", 1.0),
        triangle_from_base("base-two equal-area, L=2", 2.0),
        triangle_from_angle("30 degree equal-area", 30.0),
        triangle_from_angle("60 degree equal-area", 60.0),
    ]
    print("  label                         L        slope     height    angle   log L      FP total")
    for gauge in gauges:
        print(
            f"  {gauge.label:<28} "
            f"{gauge.L:>8.6f} {gauge.slope:>10.6f} {gauge.height:>9.6f} "
            f"{gauge.angle_deg:>7.3f} {gauge.log_cut:>9.6f} {gauge.finite_part_total:>11.9f}"
        )
    print()
    print("Interpretation:")
    print("  L=1 is the usual isosceles right triangle.")
    print("  L=2 is the equal-area triangle spread across x=0..2; its Euclidean")
    print("      angle is atan(1/4), not 30 degrees.")
    print("  A Euclidean 30 degree equal-area triangle has L=3^(1/4), so its")
    print("      log-coordinate cutoff is log(3)/4.")
    print()


def print_forced_triangle_warning() -> None:
    print("2. If one forces a base-two Euclidean 30-60-90 triangle")
    print("-" * 78)
    for angle in (30.0, 60.0):
        area, finite_total, normalized_slope, y_scale = forced_base_angle(2.0, angle)
        print(
            f"  base=2, angle={angle:>4.0f}: triangle area={area:.9f}, "
            f"log(2)+area={finite_total:.9f}"
        )
        print(
            f"    equal-area base-two slope would be {normalized_slope:.9f}; "
            f"vertical display scale to look like {angle:.0f} deg is {y_scale:.9f}"
        )
    print()
    print("So a literal 30-60-90 triangle with x-projection 2 changes the value.")
    print("To keep log(2)+1/2, either rescale the triangle or view the angle")
    print("after an anisotropic x/y display normalization.")
    print()


def print_counterterm() -> None:
    print("3. Pointwise q counterterm")
    print("-" * 78)
    print("On 0<x<1, write x = 1/x + (x - 1/x).")
    print("Thus q_ct(x)=x-1/x cancels the hyperbola's left divergence.")
    print()
    print("  epsilon        int_eps^2 1/x      int_eps^1 q_ct      sum")
    for epsilon in (1e-1, 1e-2, 1e-4, 1e-8):
        hyperbola, counterterm, total = cutoff_counterterm_sum(epsilon)
        print(
            f"  {epsilon:>9.1e} {hyperbola:>18.9f} "
            f"{counterterm:>18.9f} {total:>12.9f}"
        )
    print(f"  limit = {TARGET_AREA:.9f}")
    print()
    print("Reciprocal parity:")
    print("  x       q_ct(x)       q_ct(1/x)    sum")
    for x in (0.5, 1.0, 2.0, 3.0):
        a = q_ct(x)
        b = q_ct(1.0 / x)
        print(f"  {x:>3.1f} {a:>13.9f} {b:>13.9f} {a+b:>9.2e}")
    print()


def print_log_coordinate() -> None:
    print("4. Log-coordinate form")
    print("-" * 78)
    print("Set x=e^t.  Then dx/x=dt, so the hyperbola is uniform in t.")
    print("The original splice becomes")
    print("  int_{-infty}^0 e^(2t) dt + int_0^log(2) dt = 1/2 + log(2).")
    print("The counterterm is q_ct=e^t-e^(-t)=2*sinh(t).")
    print()
    print("For q_L, the normalized mass 2*q_L(x) dx is an exponential tail")
    print("in t with mean log(L)-1/2 and variance 1/4.")
    print()
    for L in (1.0, 2.0, 3.0 ** 0.25):
        mean_t = log(L) - 0.5
        print(
            f"  L={L:.6f}: log L={log(L):.9f}, "
            f"mean_t={mean_t:.9f}, variance_t=0.250000000"
        )
    print()


def print_mellin_tangent() -> None:
    print("5. Mellin-transform tangent")
    print("-" * 78)
    print("For q_L(x)=x/L^2 on [0,L],")
    print("  M_q(s)=int_0^L x^s q_L(x) dx = L^s/(s+2).")
    print("The triangle has moved the hyperbola pole at s=0 to a pole at s=-2.")
    print("Changing L only multiplies by L^s, i.e. translates in log coordinate.")
    print()
    print("  L        M_q(0)    M_q(1)    M_q(-1)")
    for L in (1.0, 2.0, 3.0 ** 0.25):
        print(
            f"  {L:>8.6f} {mellin_q_L(L,0.0):>9.6f} "
            f"{mellin_q_L(L,1.0):>9.6f} {mellin_q_L(L,-1.0):>10.6f}"
        )
    print()


def print_synthesis() -> None:
    print("SYNTHESIS")
    print("=" * 78)
    print("There are two complementary readings of q(x):")
    print("  q_L is a triangle-density gauge with fixed total mass 1/2.")
    print("  q_ct=x-1/x is the reciprocal-odd counterterm making the splice")
    print("      pointwise equal to 1/x+q_ct below x=1.")
    print()
    print("The 30-60-90 tangent is a gauge choice, not a new value by itself.")
    print("A true 30 degree equal-area q has log cutoff log(3)/4; a base-two")
    print("equal-area q has log cutoff log(2) and slope 1/4.")
    print()
    print("The deeper invariant is dx/x.  Logs are areas in rapidity/log")
    print("coordinate, while triangles are exponential cutoff packets added")
    print("to that invariant measure.")
    print()


def main() -> None:
    print_header()
    print_triangle_gauges()
    print_forced_triangle_warning()
    print_counterterm()
    print_log_coordinate()
    print_mellin_tangent()
    print_synthesis()


if __name__ == "__main__":
    main()
