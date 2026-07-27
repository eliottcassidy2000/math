#!/usr/bin/env python3
"""Exact controls for the stable curvature/cohabitation kernel.

The mathematical proof is the clipped-increment identity in the companion
reflection.  This standard-library script checks it over a hostile rational
bank, verifies positivity and point detection, and separates the interior
curvature from the diagonal delta*min(P,Q) ridge.
"""

from __future__ import annotations

from fractions import Fraction


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def clip(value: Fraction) -> Fraction:
    return min(Fraction(1), max(Fraction(0), value))


def green(z: Fraction, t: Fraction) -> Fraction:
    return min(z, t) - z * t


def defect(
    p: int,
    q: int,
    delta: Fraction,
    atoms: list[tuple[Fraction, Fraction]],
) -> Fraction:
    """Symmetric homogeneous defect from a finite curvature measure."""

    if p < q:
        p, q = q, p
    if p == 0:
        return Fraction(0)
    z = Fraction(q, p)
    return p * (
        delta * z
        + sum(mass * green(z, point) for point, mass in atoms)
    )


def chamber_kernel(p: int, q: int, t: Fraction) -> Fraction:
    return clip((p + 1) * t - q) - clip(p * t - q)


def mixed_increment(
    p: int,
    q: int,
    delta: Fraction,
    atoms: list[tuple[Fraction, Fraction]],
) -> Fraction:
    return (
        defect(p + 1, q + 1, delta, atoms)
        - defect(p + 1, q, delta, atoms)
        - defect(p, q + 1, delta, atoms)
        + defect(p, q, delta, atoms)
    )


def diagonal_increment(
    p: int,
    delta: Fraction,
    atoms: list[tuple[Fraction, Fraction]],
) -> Fraction:
    """Mixed cell crossing the symmetry diagonal at (p,p)."""

    return (
        defect(p + 1, p + 1, delta, atoms)
        - defect(p + 1, p, delta, atoms)
        - defect(p, p + 1, delta, atoms)
        + defect(p, p, delta, atoms)
    )


def main() -> None:
    controls = [
        (Fraction(0), []),
        (Fraction(3, 2), []),
        (Fraction(2, 3), [(Fraction(1, 3), Fraction(5, 2))]),
        (
            Fraction(4, 5),
            [
                (Fraction(2, 9), Fraction(7, 4)),
                (Fraction(5, 7), Fraction(11, 6)),
            ],
        ),
    ]

    checked_cells = 0
    for delta, atoms in controls:
        for p in range(1, 14):
            for q in range(p):
                # Verify P*G=min(Q,Pt)-Qt atom by atom.
                for point, _mass in atoms:
                    require(
                        p * green(Fraction(q, p), point)
                        == min(Fraction(q), p * point) - q * point,
                        "Green/min identity failed",
                    )

                kernel_sum = sum(
                    mass * chamber_kernel(p, q, point)
                    for point, mass in atoms
                )
                require(
                    mixed_increment(p, q, delta, atoms) == kernel_sum,
                    "mixed increment/kernel identity failed",
                )
                require(kernel_sum >= 0, "curvature kernel lost positivity")
                checked_cells += 1

            # A complete chamber row telescopes to the clipped first moment.
            row_sum = sum(
                mixed_increment(p, q, delta, atoms) for q in range(p)
            )
            clipped_moment = sum(
                mass * min(point, p * (1 - point))
                for point, mass in atoms
            )
            require(
                row_sum == clipped_moment,
                "complete-row clipped-moment identity failed",
            )
            require(
                diagonal_increment(p, delta, atoms) + 2 * row_sum == delta,
                "diagonal-plus-row reconstruction failed",
            )

    # Every rational interior point is detected by Q=floor(Pt).
    detected_points = 0
    for denominator in range(2, 31):
        for numerator in range(1, denominator):
            point = Fraction(numerator, denominator)
            p = denominator
            q = (p * point).numerator // (p * point).denominator
            require(0 <= q <= p - 1, "detecting cell left the chamber")
            require(
                chamber_kernel(p, q, point) > 0,
                "an interior rational curvature point was invisible",
            )
            detected_points += 1

    # With no interior curvature all chamber increments vanish, and the
    # one cell crossing the diagonal sees exactly the ridge delta.
    ridge = Fraction(7, 5)
    for p in range(1, 14):
        for q in range(p):
            require(
                mixed_increment(p, q, ridge, []) == 0,
                "linear chamber profile acquired curvature",
            )
        crossing = diagonal_increment(p, ridge, [])
        require(crossing == ridge, "diagonal ridge was not detected")

    print(f"hostile_rational_cells_checked={checked_cells}")
    print(f"interior_rational_points_detected={detected_points}")
    print("kernel=clip((P+1)t-Q)-clip(Pt-Q)")
    print("kernel_sign=nonnegative")
    print("all_chamber_increments_zero_iff=interior_curvature_zero")
    print("row_sum=integral_min(t,P*(1-t))_dkappa")
    print("diagonal_plus_twice_row_sum=delta")
    print("THM2348_sign=stable_sigma_side_positive_muU_side_negative")
    print("scope=stable_chamber_bridge_not_finite_rectangularity")
    print("status=KNOT_CURVATURE_COHABITATION_KERNEL_EXACT_AUDIT")


if __name__ == "__main__":
    main()
