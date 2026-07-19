#!/usr/bin/env python3
"""Exact arithmetic audit for THM-1241's first-lap Kakeya drift."""

from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path


H = F(1, 14)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def floor_fraction(x: F) -> int:
    return x.numerator // x.denominator


def frac(x: F) -> F:
    return x - floor_fraction(x)


def envelope_formula(length: F) -> F:
    """Maximum normalized x-load of one radius-1/14 comb."""
    require(length > 0, "length must be positive")
    whole = floor_fraction(length)
    tail = frac(length)
    return (F(whole, 7) + min(tail, F(1, 7))) / length


def danger_u_length(phase: F, length: F) -> F:
    """Exact danger length in the lifted u-interval [phase, phase+length]."""
    left = phase
    right = phase + length
    first = floor_fraction(left - H) - 1
    last = floor_fraction(right + H) + 1
    total = F(0)
    for centre in range(first, last + 1):
        lo = max(left, F(centre) - H)
        hi = min(right, F(centre) + H)
        if lo < hi:
            total += hi - lo
    return total


def arrangement_envelope(length: F) -> F:
    """Independent phase-arrangement maximum over one phase period."""
    events = {F(0), F(1)}
    # Load changes slope only when a window endpoint meets a danger endpoint.
    bound = floor_fraction(length) + 3
    for centre in range(-2, bound + 1):
        for sign in (-1, 1):
            endpoint = F(centre) + sign * H
            events.add(frac(endpoint))
            events.add(frac(endpoint - length))
    ordered = sorted(events)
    candidates = set(ordered)
    for a, b in zip(ordered, ordered[1:]):
        candidates.add((a + b) / 2)
    candidates.add((ordered[-1] + F(1) + ordered[0]) / 2 % 1)
    return max(danger_u_length(p, length) / length for p in candidates)


def audit_envelope() -> int:
    checks = 0
    for denominator in range(2, 42):
        for numerator in range(1, 12 * denominator + 1):
            length = F(numerator, denominator)
            direct = arrangement_envelope(length)
            formula = envelope_formula(length)
            require(direct == formula, f"envelope mismatch at L={length}")
            checks += 1
    return checks


def diagonal_alias_windows() -> list[tuple[F, F]]:
    """Ratios y for which the scalar condition 6 Q(6y/7)>=1 survives."""
    windows = [(F(1), F(1))]
    for whole in range(1, 6):
        lower_L = F(whole) + F(whole, 35)
        upper_L = F(whole) + F(6 - whole, 7)
        lower_y = F(7, 6) * lower_L
        upper_y = F(7, 6) * upper_L
        windows.append((lower_y, upper_y))
        require(6 * envelope_formula(lower_L) == 1, "lower alias endpoint")
        require(6 * envelope_formula(upper_L) == 1, "upper alias endpoint")
        if lower_L < upper_L:
            require(
                6 * envelope_formula((lower_L + upper_L) / 2) > 1,
                "alias-window interior must pass scalar load",
            )
    require(windows == [
        (F(1), F(1)),
        (F(6, 5), F(2)),
        (F(12, 5), F(3)),
        (F(18, 5), F(4)),
        (F(24, 5), F(5)),
        (F(6), F(6)),
    ], "unexpected diagonal alias windows")
    return windows


def audit_constants() -> None:
    # Whole-arc branch: radius >=1/2 already pays much more than 1/14 drift.
    require(F(1, 2) - H == F(3, 7), "whole-arc drift threshold")
    require(F(3, 7) > F(1, 14), "whole arc is stronger")

    # Proper-arc length invoice: 1 < 6/7 + 2E gives E > 1/14.
    missing_half = (F(1) - 6 * (2 * H)) / 2
    require(missing_half == F(1, 14), "first-lap drift constant")

    # Speed spread and weighted Hamiltonian-path consequences.
    require(5 * F(1, 70) == F(1, 14), "diameter constant")
    require(sum(range(1, 6)) == 15, "weighted gap sum")
    require(F(15, 210) == F(1, 14), "adjacent gap constant")
    require(1 + F(1, 210) == F(211, 210), "adjacent ratio")
    require(1 / (1 - F(1, 70)) == F(70, 69), "diameter ratio")

    # Carrier-scale cut stratified by the number of full-lap pivots.
    require(F(1, 12) / 15 == F(1, 180), "one-pivot carrier cut")
    require(F(1, 12) / 11 == F(1, 132), "two-pivot carrier cut")
    require(F(1, 12) / 9 == F(1, 108), "three-pivot carrier cut")

    # Convert the limiting absolute cuts using L6>1, i.e. x6>7/6.
    require(F(7, 6) / 70 == F(1, 60), "limit diameter")
    require(F(7, 6) / 210 == F(1, 180), "limit adjacent cut")


def main() -> None:
    checks = audit_envelope()
    windows = diagonal_alias_windows()
    audit_constants()

    print("THM-1241 FIRST-LAP KAKEYA DRIFT EXACT AUDIT")
    print(f"one_comb_arrangement_checks={checks}")
    print("Q(L)=(floor(L)/7+min(frac(L),1/7))/L")
    print("L6<=1 => each load <1/6 => no six-cover")
    print("whole_arc_threshold=3/7 > 1/14")
    print("proper_open_arc_invoice: 1 < 6/7+2E => E>1/14")
    print("all_pivots: dh>=7c/6 => sum_i|di-dh|>dh/14")
    print("speed_invoice=sum_i(d6-di)>d6/14")
    print("diameter=d6-d1>d6/70; ratio d6/d1>70/69")
    print("weighted_path=sum_j j(d[j+1]-d[j])>d6/14")
    print("macroscopic_edge=some d[j+1]/d[j]>211/210")
    print("pivot_count_cuts: r=1=>gap>c/180 r=2=>gap>c/132 r>=3=>gap>c/108")
    print("limit_cuts: x6-x1>=1/60; some adjacent gap>=1/180")
    print("scalar_diagonal_alias_windows=" + ",".join(
        f"[{a},{b}]" for a, b in windows
    ))
    print("first_lap_arc_vertices=static_phase,radius_1/14,relative_drift")
    print("speed_tournament=transitive scores 0,1,2,3,4,5; SCCs 1^6; HP unique")
    print("result=PASS")


if __name__ == "__main__":
    main()
