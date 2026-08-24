#!/usr/bin/env python3
"""Exact controls for THM-4025's owner-residue dilation semigroup.

The all-depth statement is proved algebraically in THM-4025.  This companion
checks scaled-owner invariance, the quantitative margin law, the transported
odd-index operation, the minimal ordinary-order revival and hostile stronger
monotonicity claims.  It does not prove LRC(14).
"""

from __future__ import annotations

from fractions import Fraction as F
from math import gcd
import sys


sys.stdout.reconfigure(newline="\n")


WIDTH = F(2, 63)
Q_HEIGHT = 91**6
CHECKS = 0


def require(condition: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(label)


def positive_residue(value: int, modulus: int) -> int:
    residue = value % modulus
    require(residue != 0, f"nonzero odd residue value={value},modulus={modulus}")
    return residue


def owner_pair(t: int, u: int) -> tuple[F, F]:
    divisor = gcd(t, u)
    return (
        F(positive_residue(3 * t - 4 * u, 42 * divisor), 42 * u),
        F(positive_residue(9 * t + 16 * u, 126 * divisor), 126 * u),
    )


def gate_details(t: int, U: int) -> tuple[bool, F, F, F, int, F, int]:
    require(t % 2 == 1 and 2 <= U <= t, f"gate domain t={t},U={U}")
    best1: tuple[F, int] | None = None
    best2: tuple[F, int] | None = None
    for owner in range(1, U + 1):
        first, second = owner_pair(t, owner)
        if best1 is None or first < best1[0]:
            best1 = (first, owner)
        if best2 is None or second < best2[0]:
            best2 = (second, owner)
    require(best1 is not None and best2 is not None, f"minima t={t},U={U}")
    demand = F(t * (2 * U - 1), 84 * U * (U - 1))
    residual = max(F(0), WIDTH - best1[0] - best2[0])
    return demand > residual, demand, residual, best1[0], best1[1], best2[0], best2[1]


def physical_old_strip(t: int, U: int) -> bool:
    return t % 2 == 1 and max(11, 3 * t // 4 + 1) <= U <= t


def phi(index: int) -> int:
    return 2 * index - 1


def star(first: int, second: int) -> int:
    return 2 * first * second - first - second + 1


def physical_states(bound: int):
    """One-pass cumulative owner minima for every old-strip cell."""
    rows = []
    for t in range(11, bound + 1, 2):
        best1: tuple[F, int] | None = None
        best2: tuple[F, int] | None = None
        for U in range(1, t + 1):
            first, second = owner_pair(t, U)
            if best1 is None or first < best1[0]:
                best1 = (first, U)
            if best2 is None or second < best2[0]:
                best2 = (second, U)
            if not physical_old_strip(t, U):
                continue
            require(best1 is not None and best2 is not None, f"profile t={t},U={U}")
            demand = F(t * (2 * U - 1), 84 * U * (U - 1))
            residual = max(F(0), WIDTH - best1[0] - best2[0])
            scale = gcd(t, U)
            rows.append(
                (
                    t,
                    U,
                    (t // scale, U // scale),
                    scale,
                    demand > residual,
                    demand,
                    residual,
                    best1,
                    best2,
                )
            )
    return rows


def main() -> None:
    require(Q_HEIGHT == 567869252041, "THM-4003 pair-height boundary")

    # Algebraic controls.  A scaled owner is available in the larger owner
    # set and has exactly the same two directed gaps.
    scaling_checks = 0
    for t in range(3, 102, 2):
        for U in range(2, t + 1):
            if 4 * U <= 3 * t:
                continue
            closed, demand, residual, epsilon1, _, epsilon2, _ = gate_details(t, U)
            margin = residual - demand
            for multiplier in range(1, 12, 2):
                scaled = gate_details(multiplier * t, multiplier * U)
                closed_k, demand_k, residual_k, eta1, _, eta2, _ = scaled
                delta = F(
                    t * (multiplier - 1),
                    84 * (U - 1) * (multiplier * U - 1),
                )
                require(demand - demand_k == delta, f"demand delta {t},{U},{multiplier}")
                require(eta1 <= epsilon1 and eta2 <= epsilon2, f"gap minima {t},{U},{multiplier}")
                require(residual_k >= residual, f"residual {t},{U},{multiplier}")
                require(residual_k - demand_k >= margin + delta, f"margin {t},{U},{multiplier}")
                if not closed:
                    require(not closed_k, f"survivor dilation {t},{U},{multiplier}")
                for owner in {1, max(1, U // 2), U}:
                    require(
                        owner_pair(multiplier * t, multiplier * owner) == owner_pair(t, owner),
                        f"owner invariance {t},{owner},{multiplier}",
                    )
                scaling_checks += 1

    # Treating the n-th odd number (or its square) as n transports
    # multiplication to star.  Check its semigroup laws independently of the
    # residue computation.
    star_checks = 0
    for first in range(1, 41):
        for second in range(1, 41):
            require(phi(star(first, second)) == phi(first) * phi(second), "star transport")
            require(star(first, second) == star(second, first), "star commutative")
            require(star(first, 1) == first, "star identity")
            for third in range(1, 8):
                require(
                    star(star(first, second), third) == star(first, star(second, third)),
                    "star associative",
                )
            star_checks += 1

    # The smallest physical ordinary-order revival.  The two adjacent odd
    # multipliers live on the primitive arithmetic ray (5,4), but are
    # incomparable in divisibility order.
    first = gate_details(45, 36)
    second = gate_details(55, 44)
    require(
        first == (False, F(71, 2352), F(215, 7084), F(1, 966), 23, F(1, 2772), 22),
        "(45,36) survivor details",
    )
    require(
        second == (True, F(145, 4816), F(8, 287), F(1, 1722), 41, F(17, 5166), 41),
        "(55,44) closure details",
    )
    require(first[1] - first[2] == F(-97, 595056), "survivor margin")
    require(second[1] - second[2] == F(63, 28208), "closure margin")
    require(11 % 9 != 0, "ordinary revival scales incomparable")

    rows = physical_states(501)
    by_ray: dict[tuple[int, int], list[tuple]] = {}
    for row in rows:
        by_ray.setdefault(row[2], []).append(row)

    revivals = []
    divisibility_checks = 0
    for ray, ray_rows in by_ray.items():
        ray_rows.sort(key=lambda row: row[3])
        seen_survivor = None
        status_by_scale = {row[3]: row[4] for row in ray_rows}
        for row in ray_rows:
            if not row[4] and seen_survivor is None:
                seen_survivor = row
            if row[4] and seen_survivor is not None:
                revivals.append((row[0], ray, seen_survivor[:5], row[:5]))
                break
        for scale, closed in status_by_scale.items():
            if closed:
                continue
            for multiplier in range(1, max(status_by_scale) // scale + 1, 2):
                target = scale * multiplier
                if target in status_by_scale:
                    require(not status_by_scale[target], f"divisor implication {ray},{scale},{target}")
                    divisibility_checks += 1
    revivals.sort()
    require(revivals and revivals[0][0] == 55 and revivals[0][1] == (5, 4), "minimal revival")

    # Hostiles to tempting stronger claims.
    close_base = gate_details(11, 11)
    survive_triple = gate_details(33, 33)
    require(close_base[0] and not survive_triple[0], "closure is not upward")
    require(
        gate_details(21, 16)[5] == gate_details(63, 48)[5] == F(1, 504),
        "coordinate minimum need not strictly decrease",
    )
    require(
        gate_details(33, 33)[3] > gate_details(11, 11)[3] / 3,
        "no inverse-linear epsilon bound",
    )

    ray_word = "".join(
        "C" if gate_details(5 * multiplier, 4 * multiplier)[0] else "S"
        for multiplier in range(3, 32, 2)
    )
    require(ray_word == "CCCSCCSSSSSSSSS", "ray (5,4) finite word")

    print("THM4025_OWNER_RESIDUE_ODD_DILATION_SEMIGROUP_EXACT")
    print("scope=THM4003_owner_relaxed_arithmetic_all_U;LRC_gate_only_when_U<91^6;LRC14=OPEN")
    print("owner_invariance=e_i(kt,ku)=e_i(t,u)_for_odd_k")
    print("margin_law=(B-D)(kt,kU)>=(B-D)(t,U)+t(k-1)/(84(U-1)(kU-1))")
    print("sequence_law=phi(n)=2n-1;n_star_m=2nm-n-m+1;closure(n_star_m)<=min(closure(n),closure(m))")
    print("minimal_ordinary_revival=ray_(5,4);k9_(45,36)=S;k11_(55,44)=C")
    print("minimal_revival_margins=-97/595056,63/28208")
    print(f"ray_(5,4)_physical_word_k3_to31={ray_word}")
    print(f"scaling_checks={scaling_checks};star_checks={star_checks}")
    print(f"physical_scan_cells_through_501={len(rows)};rays={len(by_ray)}")
    print(f"ordinary_revivals_through_501={len(revivals)};first={revivals[0]}")
    print(f"divisibility_survival_checks={divisibility_checks}")
    print("hostiles=closure_not_upward;epsilon_not_strict;no_1/k_coordinate_bound")
    print(f"checks={CHECKS}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
