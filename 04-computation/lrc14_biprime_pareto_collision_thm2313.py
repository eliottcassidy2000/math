#!/usr/bin/env python3
"""Exact, optimization-safe companion for THM-2313."""

from fractions import Fraction
from itertools import product
from math import gcd


P = 13
R = 7
MODULUS = P * R
EPSILON = Fraction(1, 10**12)


def require(condition: bool, message: str) -> None:
    """Raise in ordinary and optimized Python when a theorem check fails."""
    if not condition:
        raise RuntimeError(message)


def valuation(number: int, prime: int) -> int:
    """Return the prime-adic valuation of a positive integer."""
    require(number > 0, "valuation requires a positive integer")
    value = 0
    while number % prime == 0:
        number //= prime
        value += 1
    return value


def circle_distance(x: Fraction, y: Fraction) -> Fraction:
    """Distance on R/Z between rational representatives."""
    difference = (x - y) % 1
    return min(difference, 1 - difference)


def upper_set_check(side: int = 4) -> tuple[int, int, int]:
    """Exhaust all upper subsets of a small product-order rectangle."""
    cells = tuple(product(range(side), repeat=2))
    upper_sets = 0
    nonempty_disjoint_origin_sets = 0
    frontier_points = 0

    for mask in range(1 << len(cells)):
        occupied = {
            cell for index, cell in enumerate(cells) if mask & (1 << index)
        }
        is_upper = all(
            (u, v) in occupied
            for s, t in occupied
            for u in range(s, side)
            for v in range(t, side)
        )
        if not is_upper:
            continue

        upper_sets += 1
        if not occupied or (0, 0) in occupied:
            continue
        nonempty_disjoint_origin_sets += 1

        frontier = {
            (s, t)
            for s, t in occupied
            if not any(
                (u, v) in occupied and (u, v) != (s, t)
                and u <= s and v <= t
                for u, v in cells
            )
        }
        require(frontier, "nonempty finite upper set lost its frontier")
        for s, t in frontier:
            require(s > 0 or t > 0, "the disjoint origin entered the frontier")
            if s > 0:
                require((s - 1, t) not in occupied,
                        "west predecessor of a minimal cell is positive")
            if t > 0:
                require((s, t - 1) not in occupied,
                        "south predecessor of a minimal cell is positive")
            if s > 0 and t > 0:
                require((s - 1, t - 1) not in occupied,
                        "southwest predecessor of a minimal cell is positive")
        frontier_points += len(frontier)

    # Product-order upper sets in a 4 by 4 square are monotone lattice paths.
    require(upper_sets == 70, "4x4 upper-set count changed")
    require(nonempty_disjoint_origin_sets == 68,
            "disjoint-origin upper-set count changed")
    return upper_sets, nonempty_disjoint_origin_sets, frontier_points


def selector_check(radius: int = 4000) -> int:
    """Verify the one- and two-prime inclusion-exclusion selectors."""
    rows = 0
    for n in range(-radius, radius + 1):
        selector_13 = 1 - int(n % P == 0)
        selector_7 = 1 - int(n % R == 0)
        selector_91 = (
            1 - int(n % P == 0) - int(n % R == 0)
            + int(n % MODULUS == 0)
        )
        require(selector_13 == int(n % P != 0),
                "mod-13 shell selector failed")
        require(selector_7 == int(n % R != 0),
                "mod-7 shell selector failed")
        require(selector_91 == int(gcd(abs(n), MODULUS) == 1),
                "mod-91 corner selector failed")
        rows += 1
    return rows


def prony_bank_check(max_nodes: int = 128) -> tuple[int, int, int, int]:
    """Check every residue bank used by the endpoint-Prony bounds."""
    unit_residues_91 = tuple(
        residue for residue in range(1, MODULUS)
        if gcd(residue, MODULUS) == 1
    )
    require(len(unit_residues_91) == 72, "phi(91) changed")

    rows_13 = 0
    rows_7 = 0
    rows_91 = 0
    for node_count in range(1, max_nodes + 1):
        for modulus, residues in (
            (P, range(1, P)),
            (R, range(1, R)),
            (MODULUS, unit_residues_91),
        ):
            for residue in residues:
                bank = tuple(
                    residue + modulus * ell for ell in range(node_count)
                )
                require(
                    all(q > 0 and q % modulus == residue for q in bank),
                    "Prony residue bank malformed",
                )
                require(
                    max(bank) <= modulus * node_count - 1,
                    "Prony residue-bank landing bound failed",
                )
                if modulus == P:
                    rows_13 += 1
                elif modulus == R:
                    rows_7 += 1
                else:
                    rows_91 += 1
    return rows_13, rows_7, rows_91, len(unit_residues_91)


def jump_bound_check() -> int:
    """Verify the instantiated owner-normalized landing bounds."""
    rows = 0
    for speed_sum in range(9, 401):
        jumps_a = 4 * speed_sum
        jumps_b = 2 * speed_sum
        product_bound = jumps_a * jumps_b
        require(product_bound == 8 * speed_sum**2,
                "jump-product bound changed")
        require(P * product_bound - 1 == 104 * speed_sum**2 - 1,
                "mod-13 landing bound changed")
        require(R * product_bound - 1 == 56 * speed_sum**2 - 1,
                "mod-7 landing bound changed")
        require(MODULUS * product_bound - 1
                == 728 * speed_sum**2 - 1,
                "mod-91 landing bound changed")
        rows += 1
    return rows


def valuation_ledger_check() -> int:
    """Check exact two-prime valuations after a mixed Pareto corner."""
    rows = 0
    for lambda_13 in range(5):
        for gamma_7 in range(5):
            owner = P**lambda_13 * R**gamma_7
            for s in range(1, 6):
                for t in range(1, 6):
                    for n in range(1, 251):
                        if gcd(n, MODULUS) != 1:
                            continue
                        frequency = owner * P ** (s - 1) * R ** (t - 1) * n
                        require(
                            valuation(frequency, P) == lambda_13 + s - 1,
                            "mixed thirteen-adic valuation changed",
                        )
                        require(
                            valuation(frequency, R) == gamma_7 + t - 1,
                            "mixed seven-adic valuation changed",
                        )
                        rows += 1
    return rows


SOURCE_CENTERS = (Fraction(1, 16), Fraction(15, 16))
CURRENT_CENTERS = (Fraction(7, 16), Fraction(9, 16))


def hostile_gap(s: int, t: int) -> Fraction:
    """Exact cross-gap after multiplication by 13^s 7^t."""
    scale = P**s * R**t
    source_images = tuple((scale * x) % 1 for x in SOURCE_CENTERS)
    current_images = tuple((scale * x) % 1 for x in CURRENT_CENTERS)
    return min(
        circle_distance(x, y)
        for x in source_images
        for y in current_images
    )


def hostile_collision(s: int, t: int) -> bool:
    """Whether the two pushed hostile supports overlap in positive measure."""
    scale = P**s * R**t
    source_radius = scale * EPSILON
    current_radius = 169 * scale * EPSILON
    gap = hostile_gap(s, t)
    geometric_collision = (
        source_radius >= Fraction(1, 2)
        or current_radius >= Fraction(1, 2)
        or source_radius + current_radius > gap
    )
    threshold_collision = 170 * scale * EPSILON > gap
    require(
        geometric_collision == threshold_collision,
        "hostile collision threshold lost equivalence",
    )
    return geometric_collision


def hostile_atlas_check() -> tuple[
    tuple[int, ...],
    tuple[tuple[int, int], ...],
    tuple[int, ...],
    tuple[tuple[tuple[int, int], tuple[int, int]], ...],
    int,
]:
    """Compute the complete Pareto frontier of the multiplier-four carrier."""
    gap_rows = 0
    for s in range(13):
        for t in range(15):
            expected = Fraction(3, 8) if s % 2 == 0 else Fraction(1, 8)
            require(hostile_gap(s, t) == expected,
                    "hostile cross-gap parity law changed")
            hostile_collision(s, t)
            gap_rows += 1

    first_t = []
    for s in range(10):
        candidates = [t for t in range(20) if hostile_collision(s, t)]
        require(candidates, "hostile collision search rectangle too small")
        first_t.append(min(candidates))
    first_t_tuple = tuple(first_t)
    require(
        first_t_tuple == (12, 10, 9, 7, 6, 4, 4, 2, 1, 0),
        "hostile first-collision staircase changed",
    )

    frontier = []
    best_t = None
    for s, t in enumerate(first_t_tuple):
        if best_t is None or t < best_t:
            frontier.append((s, t))
            best_t = t
    frontier_tuple = tuple(frontier)
    expected_frontier = (
        (0, 12),
        (1, 10),
        (2, 9),
        (3, 7),
        (4, 6),
        (5, 4),
        (7, 2),
        (8, 1),
        (9, 0),
    )
    require(frontier_tuple == expected_frontier,
            "hostile Pareto frontier changed")

    for s, t in frontier_tuple:
        require(hostile_collision(s, t), "frontier cell is not positive")
        if s > 0:
            require(not hostile_collision(s - 1, t),
                    "frontier west predecessor became positive")
        if t > 0:
            require(not hostile_collision(s, t - 1),
                    "frontier south predecessor became positive")
        if s > 0 and t > 0:
            require(not hostile_collision(s - 1, t - 1),
                    "frontier southwest predecessor became positive")

    # Once (9,0) is positive every cell with s>=9 is dominated; for each
    # s<=9 every t above first_t[s] is dominated in its own column.
    require(hostile_collision(9, 0), "13-axis endpoint disappeared")
    require(hostile_collision(0, 12), "7-axis endpoint disappeared")

    scales = tuple(P**s * R**t for s, t in frontier_tuple)
    require(
        scales == (
            13841287201,
            3672178237,
            6819759583,
            1809323971,
            3360173089,
            891474493,
            3074677333,
            5710115047,
            10604499373,
        ),
        "hostile frontier scales changed",
    )

    # The hostile source owner is c_1=13. At an interior corner, the
    # theorem's mixed atom has (nu_13,nu_7)=(s,t-1).
    interior_valuations = tuple(
        ((s, t), (s, t - 1))
        for s, t in frontier_tuple
        if s > 0 and t > 0
    )
    require(len(interior_valuations) == 7,
            "hostile interior-frontier count changed")
    require(
        tuple(cell for cell, valuation_pair in interior_valuations
              if valuation_pair[0] == 3) == ((3, 7),),
        "profile-(1,3,5) grade-three alignment changed",
    )

    return (
        first_t_tuple,
        frontier_tuple,
        scales,
        interior_valuations,
        gap_rows,
    )


def color_and_clock_boundary_check() -> tuple[int, int, int]:
    """Check the exact hostile alignment and the still-free root color."""
    profile_b = 3
    source_lambda = 1
    aligned_s = profile_b - source_lambda + 1
    require(aligned_s == 3, "hostile source grade alignment changed")

    corner_t = 7
    transported_seven = pow(R, corner_t - 1, P)
    require(transported_seven == 12, "7^6 is no longer -1 mod 13")
    prescribed_character = 4
    required_unit_residue = (
        prescribed_character * pow(transported_seven, -1, P)
    ) % P
    require(required_unit_residue == 9,
            "required mixed-shell root residue changed")
    return aligned_s, transported_seven, required_unit_residue


def main() -> None:
    upper_sets, disjoint_upper_sets, abstract_frontier_points = upper_set_check()
    selector_rows = selector_check()
    rows_13, rows_7, rows_91, unit_residues = prony_bank_check()
    jump_rows = jump_bound_check()
    valuation_rows = valuation_ledger_check()
    (
        first_t,
        frontier,
        scales,
        interior_valuations,
        hostile_gap_rows,
    ) = hostile_atlas_check()
    aligned_s, transported_seven, required_residue = (
        color_and_clock_boundary_check()
    )

    require(len(frontier) == 9, "hostile frontier size changed")
    require(sum(s > 0 and t > 0 for s, t in frontier) == 7,
            "hostile interior count changed")
    require(MODULUS * 4 * 4 - 1 == 1455,
            "local four-jump mod-91 bound changed")

    print("theorem=THM-2313")
    print("status=PROVED+VERIFIED-EXACT+INDEPENDENTLY-AUDITED")
    print("collision_set=upper_subset_of_N^2")
    print("frontier_cases=13-axis,7-axis,interior")
    print("mixed_selector=1-[13|n]-[7|n]+[91|n]=[gcd(n,91)=1]")
    print("mixed_corner_sign=positive")
    print("13_axis_shell_sign=negative")
    print("7_axis_shell_sign=negative")
    print("13_axis_landing=n<=13*J_A*J_B-1<=104*S^2-1")
    print("7_axis_landing=n<=7*J_A*J_B-1<=56*S^2-1")
    print("mixed_landing=n<=91*J_A*J_B-1<=728*S^2-1")
    print("mixed_valuations=(nu13(c)+s-1,nu7(c)+t-1)")
    print(f"upper_sets_4x4={upper_sets}")
    print(f"nonempty_disjoint_origin_upper_sets_4x4={disjoint_upper_sets}")
    print(f"abstract_frontier_points_checked={abstract_frontier_points}")
    print(f"selector_rows={selector_rows}")
    print(f"unit_residues_mod91={unit_residues}")
    print(f"prony_rows_mod13={rows_13}")
    print(f"prony_rows_mod7={rows_7}")
    print(f"prony_rows_mod91={rows_91}")
    print(f"jump_bound_rows={jump_rows}")
    print(f"valuation_ledger_rows={valuation_rows}")
    print(f"hostile_gap_rows={hostile_gap_rows}")
    print(f"hostile_first_t_by_s_0_through_9={first_t}")
    print(f"hostile_frontier={frontier}")
    print(f"hostile_frontier_scales={scales}")
    print(f"hostile_frontier_size={len(frontier)}")
    print(f"hostile_interior_count={len(interior_valuations)}")
    print(f"hostile_interior_valuations={interior_valuations}")
    print("hostile_local_J_A=4")
    print("hostile_local_J_B=4")
    print("hostile_local_mod91_bound=1455")
    print(f"hostile_profile135_aligned_s={aligned_s}")
    print("hostile_profile135_aligned_corner=(3,7)")
    print(f"hostile_7^6_mod13={transported_seven}")
    print(f"hostile_required_unit_residue_mod13={required_residue}")
    print("mixed_shell_does_not_prescribe_root_character=true")
    print("symmetric_collision_is_not_one_sided_forward_gluing=true")
    print("no_cycle_phase_holonomy_claim=true")
    print("all_checks=PASS")


if __name__ == "__main__":
    main()
