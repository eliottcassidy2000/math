#!/usr/bin/env python3
"""Exact certificate for THM-4003's scale-two component erosion.

The symbolic proof is in THM-4003.  This companion checks the two component
widths, the odd-integer boundary, the four directed owner-event residues, and
the finite owner-relaxed certificate through t <= 1001.  It does not search
for counterexamples and does not prove LRC(14).
"""

from __future__ import annotations

from fractions import Fraction as F
from math import gcd
import sys


sys.stdout.reconfigure(newline="\n")


CHECKS = 0


def require(condition: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(label)


def positive_residue(value: int, modulus: int) -> int:
    """Return the least strictly positive directed residue."""
    residue = value % modulus
    require(residue != 0, f"nonzero odd residue value={value},modulus={modulus}")
    return residue


def owner_gaps(t: int, u: int) -> tuple[F, F, F, F]:
    """Four exact inward distances for one hypothetical owner u."""
    d = gcd(t, u)
    return (
        F(positive_residue(3 * t - 4 * u, 42 * d), 42 * u),
        F(positive_residue(9 * t + 16 * u, 126 * d), 126 * u),
        F(positive_residue(9 * t - 110 * u, 126 * d), 126 * u),
        F(positive_residue(3 * t + 38 * u, 42 * d), 42 * u),
    )


def main() -> None:
    q_height = 91**6
    require(q_height == 567869252041, "THM-3818 pair-height constant")
    boundaries = (F(2, 21), F(8, 63), F(55, 63), F(19, 21))
    widths = (boundaries[1] - boundaries[0], boundaries[3] - boundaries[2])
    require(widths == (F(2, 63), F(2, 63)), "scale-two component widths")

    # THM-3995's per-component collars retain u <= U.
    for U in range(11, 1002):
        retained = F(2, 63) - F(1, 42 * U) - F(1, 126 * U)
        require(retained == F(2 * (U - 1), 63 * U), f"retained width U={U}")

    for U in (11, q_height - 1, q_height, q_height + 1, 2 * q_height + 17):
        height = min(U, q_height)
        pair_height_component = F(2 * height - 1, 84 * U * (height - 1))
        distinctness_component = F(2 * U - 1, 84 * U * (U - 1))
        require(
            pair_height_component >= distinctness_component,
            f"pair-height upgrade U={U}",
        )

    def top_balance_bound(ratio: F) -> F:
        return max(F(0), F(3) * ratio / 154 - F(11, 154))

    require(top_balance_bound(F(11)) == F(1, 7), "all-type top-balance threshold")
    require(
        top_balance_bound(F(1001, 189)) == F(2, 63),
        "scale-two top-balance threshold",
    )
    require(top_balance_bound(F(22)) == F(5, 14), "large-ratio truncation boundary")

    # Distinct boundary owners strengthen the body-deep component to
    # 1/(84U)+1/(84(U-1)); a common boundary owner gives a still longer whole
    # tooth.  Its t-image can fit in a retained component only under the
    # quadratic inequality below.  For odd t this removes one old U-layer,
    # and a second layer when t == 1 (mod 4).
    boundary_checks = 0
    removed_rows = []
    for t in range(3, 2002, 2):
        for U in range(1, t + 1):
            body_component = F(2 * U - 1, 84 * U * (U - 1)) if U > 1 else F(1)
            geometric = t * body_component <= F(2 * (U - 1), 63 * U)
            arithmetic = 3 * t * (2 * U - 1) <= 8 * (U - 1) ** 2
            integer_threshold = U >= (3 * t) // 4 + 2 + (t % 4 == 1)
            require(
                geometric == arithmetic == integer_threshold,
                f"odd boundary t={t},U={U}",
            )
            boundary_checks += 1
        old_minimum = max(11, (3 * t) // 4 + 1)
        new_minimum = (3 * t) // 4 + 2 + (t % 4 == 1)
        for U in range(old_minimum, min(t + 1, new_minimum)):
            require(3 * t * (2 * U - 1) > 8 * (U - 1) ** 2, f"removed layer t={t},U={U}")
            removed_rows.append((t, U))
    require(removed_rows[0] == (13, 11), "first nonvacuous removed layer")

    # THM-3995's oddness is load-bearing: its first and fourth directed
    # residues both vanish in the legal even arithmetic control (t,u)=(12,9).
    even_modulus = 42 * gcd(12, 9)
    require((3 * 12 - 4 * 9) % even_modulus == 0, "even left-wall collision")
    require((3 * 12 + 38 * 9) % even_modulus == 0, "even right-wall collision")

    # Minimize each correctly oriented event distance over every hypothetical
    # owner 1 <= u <= U.  This is conservative: a real eleven-owner body can
    # only have larger collars.  Enumerate every old conditional strip cell
    # through t=1001 and record the closures beyond the symbolic one-layer cut.
    exact_cells = 0
    uniform_closed = []
    exact_extra = []
    survivors = []
    for t in range(11, 1002, 2):
        old_minimum = (3 * t) // 4 + 1
        best: list[F | None] = [None, None, None, None]
        for U in range(1, t + 1):
            candidates = owner_gaps(t, U)
            require(
                candidates[0] == candidates[3] and candidates[1] == candidates[2],
                f"reflection identities t={t},u={U}",
            )
            for index, candidate in enumerate(candidates):
                if best[index] is None or candidate < best[index]:
                    best[index] = candidate
            if U < max(11, old_minimum):
                continue
            gaps = tuple(value for value in best if value is not None)
            require(len(gaps) == 4, f"four gaps t={t},U={U}")
            residual = max(
                F(0),
                F(2, 63) - gaps[0] - gaps[1],
                F(2, 63) - gaps[2] - gaps[3],
            )
            body_image = t * F(2 * U - 1, 84 * U * (U - 1))
            simple_gate = 3 * t * (2 * U - 1) > 8 * (U - 1) ** 2
            exact_gate = body_image > residual
            require(not simple_gate or exact_gate, f"exact contains simple t={t},U={U}")
            if simple_gate:
                uniform_closed.append((t, U))
            elif exact_gate:
                exact_extra.append((t, U))
            else:
                survivors.append((t, U))
            exact_cells += 1

    require(len(exact_extra) == 77, "frozen exact-extra count")
    require(exact_extra[0] == (11, 11), "first exact-extra pair")
    require(exact_extra[-1] == (185, 141), "last exact-extra pair through 1001")

    print("THM4003_SCALE_TWO_COMPONENT_EROSION_EXACT")
    print("scope=odd_t;scale2_(2,1,9);conditional_t>=U;LRC14=OPEN")
    print("obstruction_component_width=2/63")
    print("uniform_retained_component_width=2(U-1)/(63U)")
    print("pair_height_body_component_width>=(2H-1)/(84U(H-1));H=min(U,91^6)")
    print("general_failure_gate=3t(2H-1)<=8(H-1)(U-1)")
    print("top_balance=Ulambda>=max(0,(3U/V-11)/154)_for_U/V<22;>=5/14_after")
    print("top_balance_closures=all17_if_U/V>=11;scale2_if_U/V>=1001/189")
    print("body_component_width>=(2U-1)/(84U(U-1))")
    print("failure_requires=3t(2U-1)<=8(U-1)^2")
    print("odd_t_equivalent=U>=floor(3t/4)+2+indicator(t_mod_4=1)")
    print("reflection_identity=e1=e4_and_e2=e3")
    print(f"boundary_checks={boundary_checks}")
    print(f"removed_old_boundary_rows_through_2001={len(removed_rows)}")
    print(f"first_removed_rows={tuple(removed_rows[:8])}")
    print(f"exact_old_strip_cells_through_1001={exact_cells}")
    print(f"uniform_closed_through_1001={len(uniform_closed)}")
    print(f"exact_extra_closed={len(exact_extra)}")
    print(f"exact_extra_pairs={tuple(exact_extra)}")
    print(f"owner_relaxed_survivors={len(survivors)};first={tuple(survivors[:12])}")
    print("even_wall_collision_controls=2")
    print(f"checks={CHECKS}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
