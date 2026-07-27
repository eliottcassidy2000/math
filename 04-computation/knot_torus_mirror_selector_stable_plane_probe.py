#!/usr/bin/env python3
"""Exact arithmetic audit for torus-knot mirror selector planes.

For K_g=T(2,2g+1), g>=2, the proof uses two switched
Levine--Tristram inertia gauges.  This companion checks:

* their coefficient matrix and mixture support function;
* domination of every possible two-frequency inertia support;
* the sharp lower and upper stable-plane envelopes;
* both rigid defect endpoints; and
* the one-atom curvature extremal and its two moment budgets.

The finite loops are hostile arithmetic controls for identities proved
symbolically in the accompanying note; they are not a finite substitute
for the all-g argument.
"""

from __future__ import annotations

from fractions import Fraction


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def switched_lower(g: int, p: int, q: int) -> int:
    return max(
        abs(g * p + (2 - g) * q),
        abs((2 - g) * p + g * q),
    )


def blanchfield_floor(g: int, p: int, q: int) -> int:
    return p + q + (g - 1) * abs(p - q)


def lower_envelope(g: int, delta: Fraction, x: Fraction, y: Fraction):
    length = 2 * (g - 1)
    return max(
        (2 * g - delta) * x,
        2 * x + length * y,
        2 * g * y,
    )


def upper_envelope(g: int, delta: Fraction, x: Fraction, y: Fraction):
    return (2 * g - delta) * max(x, y) + delta * y


def main() -> None:
    for g in range(2, 31):
        determinant = g * g - (2 - g) ** 2
        require(determinant == 4 * (g - 1), "selector determinant changed")

        for p in range(9):
            for q in range(9):
                require(
                    switched_lower(g, p, q) == blanchfield_floor(g, p, q),
                    "switched selector lost the Blanchfield mixture floor",
                )

        length = 2 * (g - 1)
        for x_num in range(9):
            for y_num in range(9):
                x = Fraction(x_num, 7)
                y = Fraction(y_num, 7)

                # A general two-frequency inertia row has support
                # epsilon*X+d*Y, epsilon in {0,1,2},
                # 0<=d<=2g-epsilon.  The extremal bank below dominates it.
                bank = max(2 * x + length * y, 2 * g * y)
                for epsilon in range(3):
                    for difference in range(2 * g - epsilon + 1):
                        require(
                            epsilon * x + difference * y <= bank,
                            "full inertia bank escaped the two extremal gauges",
                        )

                for delta_num in range(0, 2 * length + 1):
                    delta = Fraction(delta_num, 2)
                    if delta > length:
                        continue
                    lower = lower_envelope(g, delta, x, y)
                    upper = upper_envelope(g, delta, x, y)
                    require(lower <= upper, "stable envelope order reversed")

                    if delta == 0 or delta == length:
                        require(
                            lower == upper,
                            "rigid endpoint retained a plane ambiguity",
                        )
                    elif x > y > 0:
                        require(
                            lower < upper,
                            "interior defect unexpectedly rigid",
                        )

                    if 0 < delta < length:
                        t0 = delta / (2 * length - delta)
                        mass = (2 * length - delta) / 2
                        require(0 < t0 < 1 and mass > 0, "bad curvature atom")
                        require(
                            mass * (1 - t0) == length - delta,
                            "left curvature budget changed",
                        )
                        require(
                            mass * t0 == delta / 2,
                            "right curvature budget changed",
                        )

    # Trefoil boundary: the switched matrix has rank one.
    require(1 * 1 - (2 - 1) ** 2 == 0, "trefoil determinant changed")

    print("family=K_g=T(2,2g+1),g>=2")
    print("selector_matrix=[[g,2-g],[2-g,g]]")
    print("selector_determinant=4(g-1)")
    print("context_floor=P+Q+(g-1)|P-Q|")
    print("context_floor=real_Blanchfield_dimension")
    print("full_two_frequency_inertia_bank=dominated_by_extreme_switched_pair_and_odd_axis")
    print("balanced_selector_ceiling=2_per_copy")
    print("defect_interval=0<=delta<=2(g-1)")
    print("profile_bounds=delta*z<=h(z)<=min(2(g-1)z,delta(1+z)/2)")
    print("curvature_budgets=A<=2(g-1)-delta,B<=delta/2")
    print("rigid_endpoints=delta=0_and_delta=2(g-1)")
    print("interior_defects=continuum_abstract_plane_ambiguity")
    print("trefoil_boundary=selector_determinant_zero")
    print("audit_range=g_2_to_30_exact_fraction_controls")
    print("status=TORUS_MIRROR_SELECTOR_PLANE_EXACT_AUDIT")


if __name__ == "__main__":
    main()
