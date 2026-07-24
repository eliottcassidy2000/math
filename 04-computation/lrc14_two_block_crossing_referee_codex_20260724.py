#!/usr/bin/env python3
"""Exact referee for THM-2145.

This script has three independent jobs:

1. enumerate all 7-subsets of {1,...,13} on the exact danger-comb boundary
   arrangement and find the sharp safe-mass floor;
2. recheck the minimizing core by direct rational interval-union merging; and
3. certify the Jackson-kernel error and the two crossing margins using only
   Fraction arithmetic and the classical rational inequality pi < 355/113.
"""

from fractions import Fraction
from itertools import combinations


RADIUS = Fraction(1, 14)
PI2_UPPER = Fraction(355 * 355, 113 * 113)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def circle_fraction(x: Fraction) -> Fraction:
    return x % 1


def is_safe(speed: int, t: Fraction) -> bool:
    phase = circle_fraction(speed * t)
    return min(phase, 1 - phase) >= RADIUS


def boundary_cells(max_speed: int = 13) -> tuple[Fraction, ...]:
    points = {Fraction(0), Fraction(1)}
    for speed in range(1, max_speed + 1):
        for j in range(speed):
            points.add(circle_fraction(Fraction(14 * j - 1, 14 * speed)))
            points.add(circle_fraction(Fraction(14 * j + 1, 14 * speed)))
    return tuple(sorted(points))


def seven_core_masses() -> tuple[dict[tuple[int, ...], Fraction], tuple[Fraction, ...]]:
    points = boundary_cells(13)
    masses = {core: Fraction(0) for core in combinations(range(1, 14), 7)}
    for left, right in zip(points, points[1:]):
        midpoint = (left + right) / 2
        safe = {speed for speed in range(1, 14) if is_safe(speed, midpoint)}
        width = right - left
        for core in combinations(sorted(safe), 7):
            masses[core] += width
    return masses, points


def danger_intervals(speed: int) -> list[tuple[Fraction, Fraction]]:
    radius = Fraction(1, 14 * speed)
    pieces: list[tuple[Fraction, Fraction]] = []
    for j in range(speed):
        center = Fraction(j, speed)
        left, right = center - radius, center + radius
        if left < 0:
            pieces.append((Fraction(0), right))
            pieces.append((1 + left, Fraction(1)))
        elif right > 1:
            pieces.append((left, Fraction(1)))
            pieces.append((Fraction(0), right - 1))
        else:
            pieces.append((left, right))
    return pieces


def union_length(intervals: list[tuple[Fraction, Fraction]]) -> Fraction:
    ordered = sorted(intervals)
    total = Fraction(0)
    left, right = ordered[0]
    for next_left, next_right in ordered[1:]:
        if next_left <= right:
            right = max(right, next_right)
        else:
            total += right - left
            left, right = next_left, next_right
    return total + right - left


def jackson_coefficient(n: int, k: int) -> int:
    """Unnormalized Fourier coefficient C_k of F_n^2."""
    require(n >= 2 and 0 <= k <= 2 * n - 2, "Jackson coefficient index")
    if k <= n:
        numerator = 4 * n**3 - 6 * n * k**2 + 2 * n + 3 * k**3 - 3 * k
    else:
        d = 2 * n - k
        numerator = d**3 - d
    require(numerator % 6 == 0, "Jackson coefficient integrality")
    value = numerator // 6
    require(value > 0, "Jackson coefficient positivity")
    return value


def jackson_data(n: int) -> tuple[int, Fraction, Fraction]:
    """Return C_0, odd-mode sum S, and rigorous upper bound for eta."""
    c0 = n * (2 * n * n + 1) // 3
    require(c0 == jackson_coefficient(n, 0), "Jackson C0 formula")
    odd_sum = sum(
        (Fraction(jackson_coefficient(n, k), k * k) for k in range(1, 2 * n - 2, 2)),
        Fraction(0),
    )
    # eta = 2 int ||t|| J_n(t) dt
    #     = 1/2 - 4 S/(pi^2 C_0).
    # Since pi^2 < PI2_UPPER, this is a strict upper bound.
    eta_upper = Fraction(1, 2) - 4 * odd_sum / (PI2_UPPER * c0)
    require(eta_upper > 0, "Jackson eta upper bound")
    return c0, odd_sum, eta_upper


def crossing_margin(alpha: Fraction, beta: Fraction, eta: Fraction) -> Fraction:
    return (alpha - 6 * eta) * (beta - 7 * eta) - 13 * eta


def first_certified_n(alpha: Fraction, beta: Fraction, stop: int) -> int | None:
    for n in range(2, stop + 1):
        _, _, eta_upper = jackson_data(n)
        if alpha - 6 * eta_upper > 0 and beta - 7 * eta_upper > 0:
            if crossing_margin(alpha, beta, eta_upper) > 0:
                return n
    return None


def main() -> None:
    masses, points = seven_core_masses()
    minimum = min(masses.values())
    minimizers = [core for core, mass in masses.items() if mass == minimum]
    expected_core = (1, 5, 7, 8, 9, 11, 13)
    require(len(points) == 178, "boundary point count")
    require(len(masses) == 1716, "seven-core count")
    require(minimum == Fraction(45107, 229320), "safe-mass minimum")
    require(minimizers == [expected_core], "unique minimizing core")

    bad_pieces = [
        interval for speed in expected_core for interval in danger_intervals(speed)
    ]
    bad_union = union_length(bad_pieces)
    independent_safe = 1 - bad_union
    require(bad_union == Fraction(184213, 229320), "independent bad union")
    require(independent_safe == minimum, "independent safe mass")
    common_denominator = 5_045_040
    require(bad_union * common_denominator == 4_052_686, "raw bad numerator")
    require(independent_safe * common_denominator == 992_354, "raw safe numerator")

    alpha = Fraction(61, 273)  # arbitrary distinct six-block, THM-1234
    beta_general = Fraction(15, 154)  # arbitrary distinct seven-block, THM-1221
    beta_core = minimum

    # General 6+7 split.
    n_general = 293
    c0_general, s_general, eta_general = jackson_data(n_general)
    eta0_general = Fraction(1439, 1_000_000)
    floor_s_general = s_general.numerator // s_general.denominator
    require(eta_general < eta0_general, "general eta cutoff")
    gap_general = crossing_margin(alpha, beta_general, eta0_general)
    require(
        gap_general == Fraction(548679648961, 10510500000000000) > 0,
        "general crossing margin",
    )
    require(
        first_certified_n(alpha, beta_general, n_general) == n_general,
        "general Jackson certificate minimality",
    )

    # Defect-6 specialization with retained core E subset {1,...,13}.
    n_core = 150
    c0_core, s_core, eta_core = jackson_data(n_core)
    require(c0_core == 2_250_050, "core C0")
    require(s_core > 2_760_290, "core odd-mode floor")
    floor_bound_core = (
        Fraction(1, 2)
        - Fraction(4 * 12_769 * 2_760_290, 126_025 * c0_core)
    )
    require(
        floor_bound_core == Fraction(159340717, 56712510250),
        "core compact eta bound",
    )
    eta0_core = Fraction(281, 100_000)
    require(eta_core < floor_bound_core < eta0_core, "core eta cutoff")
    gap_core = crossing_margin(alpha, beta_core, eta0_core)
    require(
        gap_core == Fraction(322476924229, 7825545000000000) > 0,
        "core crossing margin",
    )
    require(
        first_certified_n(alpha, beta_core, n_core) == n_core,
        "core Jackson certificate minimality",
    )

    print("THM-2145 exact two-block spectral-crossing referee")
    print(f"boundary points={len(points)}, open cells={len(points)-1}, seven-cores={len(masses)}")
    print(f"unique safe-mass minimizer={expected_core}: {minimum}")
    print(f"independent union check: bad={bad_union}, safe={independent_safe}")
    print("raw common denominator 5045040: bad=4052686, safe=992354")
    print(
        f"general N={n_general}, H={2*n_general-2}, C0={c0_general}, "
        f"floor(S)={floor_s_general}, eta_upper<{eta0_general}"
    )
    print(f"general crossing gap at eta0={gap_general} > 0; N=292 certificate fails")
    print(
        f"core N={n_core}, H={2*n_core-2}, C0={c0_core}, "
        f"S>2760290, eta_upper<{eta0_core}"
    )
    print(f"core crossing gap at eta0={gap_core} > 0; N=149 certificate fails")
    print("normalized retained-core carry bounds: 584*70=40880, 298*70=20860")
    print("ALL EXACT ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
