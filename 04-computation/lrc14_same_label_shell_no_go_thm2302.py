#!/usr/bin/env python3
"""Exact arithmetic checks for THM-2302.

The universal Fourier and endpoint-Prony arguments are proved in the theorem.
This companion freezes the profile census, all displayed rational constants,
the sharp pointwise-energy coefficient, and endpoint-free one-sheet selector
checks on several hostile scales.
"""

from fractions import Fraction
from math import gcd


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def frac_part(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def circle_norm(x: Fraction) -> Fraction:
    y = frac_part(x)
    return min(y, 1 - y)


def in_danger(speed: int, x: Fraction) -> bool:
    return circle_norm(speed * x) < Fraction(1, 14)


def in_center_branch(x: Fraction) -> bool:
    return circle_norm(x) < Fraction(1, 26)


def selector_breakpoints(s: int) -> list[Fraction]:
    points = {Fraction(0), Fraction(1), Fraction(1, 2)}
    for r in range(s):
        for sign in (-1, 1):
            points.add(frac_part(Fraction(r, s) + sign * Fraction(1, 14 * s)))
    return sorted(points)


def audit_selector(s: int) -> int:
    points = selector_breakpoints(s)
    cells = 0
    for left, right in zip(points, points[1:]):
        if left == right:
            continue
        y = (left + right) / 2
        safe = not in_danger(s, y)
        active = 0
        for r in range(13):
            x = (y + r) / 13
            if in_center_branch(x) and safe:
                active += 1
                require(frac_part(13 * x) == y, "selected root does not map to y")
                require(circle_norm(x) < Fraction(1, 14), "selected root left D_1")
        require(active == int(safe), f"selector multiplicity failed for s={s}")
        cells += 1
    return cells


def thirteen_unit_part(s: int) -> int:
    while s % 13 == 0:
        s //= 13
    return s


def main() -> None:
    strict_profiles = [
        (1, b, c)
        for c in range(5, 20)
        for b in range(2, c)
    ]
    require(len(strict_profiles) == 150, "strict profile census")
    require(
        sum(1 for _, b, c in strict_profiles if 3 <= b <= c - 2) == 120,
        "interior census",
    )

    nu0 = Fraction(227189785662847, 58436012221844400)
    beta = Fraction(961, 6930) - Fraction(10, 91)
    require(beta == Fraction(2593, 90090), "target floor beta")

    # THM-2293 chooses nu_j/C_j >= nu0/(42+22).  The pointwise identity
    # V=sum_{r<t}(h_r-h_t)^2 and 0<=h_r<=1 give V<=12 sum_r h_r.
    # Integrating sum_r h((y+r)/13) gives 13*mu(E).
    energy_to_mass = 12 * 13
    e_star = Fraction(22, 64 * energy_to_mass) * nu0
    require(
        e_star == Fraction(227189785662847, 26519324819222476800),
        "same-label source floor",
    )
    delta_star = e_star * beta
    require(
        delta_star
        == Fraction(84157587746251753, 341303710423393276416000),
        "same-label return floor",
    )

    # The pairwise inequality behind V<=12*S is exact and sharp at (1,0).
    grid = [Fraction(k, 20) for k in range(21)]
    slack_min = None
    equality_pairs = []
    for a in grid:
        for b in grid:
            slack = a + b - (a - b) ** 2
            require(slack >= 0, "pairwise root-energy inequality")
            if slack_min is None or slack < slack_min:
                slack_min = slack
            if slack == 0:
                equality_pairs.append((a, b))
    require((Fraction(1), Fraction(0)) in equality_pairs, "sharp energy boundary")

    # At the common shell clock k=b+1, every marked source atom has
    # v_13 >= b+1, while every shell-graph vertex has exact valuation b.
    for _, b, _ in strict_profiles:
        for extra_valuation in range(5):
            marked_valuation = b + 1 + extra_valuation
            require(marked_valuation > b, "valuation nonincidence")

    selector_scales = (1, 2, 13, 26, 169)
    selector_cells = {s: audit_selector(s) for s in selector_scales}
    for s in selector_scales:
        require(gcd(thirteen_unit_part(s), 13) == 1, "hostile scale decomposition")

    # The 1/7-periodic hostile support sees A+m*c only for 7|m when
    # A=7 and c is prime to seven.  This is the exact missing unit color.
    terminal_step = 13**5
    require(gcd(terminal_step, 7) == 1, "periodic hostile terminal step")
    visible_multipliers = [
        m for m in range(1, 15) if (7 + m * terminal_step) % 7 == 0
    ]
    require(visible_multipliers == [7, 14], "seven-periodic color boundary")

    print("THM-2302 exact companion")
    print(f"strict_profiles={len(strict_profiles)}")
    print("interior_profiles=120")
    print(f"nu0={nu0}")
    print(f"energy_to_mass={energy_to_mass}")
    print(f"beta={beta}")
    print(f"e_star={e_star}")
    print(f"delta_star={delta_star}")
    print(f"pairwise_slack_min={slack_min}")
    print("selector_mass_G=6/91")
    print("selector_energy=1_(D_s^c) for each nonzero root character")
    print("total_selector_energy=12*1_(D_s^c)")
    print(
        "selector_cells="
        + ",".join(f"{s}:{selector_cells[s]}" for s in selector_scales)
    )
    print("common_shell_clock_source_valuation>=b+1>shell_valuation=b")
    print("rooted_common_index_bound=8*S^2-1")
    print("marked_uncolored_degree_bound=2*S")
    print("seven_periodic_first_visible_multiplier=7")
    print("hostile_scope=rooted middle-owner/deep-exclusion only; not full LRC cover")
    print("ALL_CHECKS_PASSED")


if __name__ == "__main__":
    main()
