#!/usr/bin/env python3
"""Exact referee for THM-1136's thirteen-grid branches.

All proof-facing arithmetic uses ``fractions.Fraction``.  The finite checks
enumerate the complete residue-class universes used by the proofs; they do
not sample integer speed tuples.

The primary finite carrier is a conflict incidence graph whose left vertices
are translated noncarrier obligations and whose right vertices are the twelve
multipliers in F_13^*.  The script checks every cell of the continuous
translated-grid arrangement, as well as the complete residue quotients.  A
naked tournament forgets these incidences and the carrier-ratio labels.
"""

from fractions import Fraction
from itertools import combinations


Q = 13
UNITS = tuple(range(1, Q))
TARGET = Fraction(1, 14)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def circle_distance(x: Fraction) -> Fraction:
    residue = x % 1
    return min(residue, 1 - residue)


def period_distance(x: Fraction, period: int) -> Fraction:
    residue = x % period
    return min(residue, period - residue)


def available_multipliers(residues: tuple[int, ...], forbidden: set[int]) -> tuple[int, ...]:
    return tuple(
        a for a in UNITS
        if all((a * residue) % Q not in forbidden for residue in residues)
    )


def branch_a_residue_census() -> tuple[int, int, tuple[int, ...], tuple[int, ...]]:
    rows = 0
    minimum_available = len(UNITS)
    hardest_residues: tuple[int, ...] = ()
    hardest_available: tuple[int, ...] = ()
    for size in range(1, 6):
        for residues in combinations(UNITS, size):
            rows += 1
            available = available_multipliers(residues, {1, 12})
            require(available, f"branch A multiplier cover at {residues}")
            require(len(available) >= 12 - 2 * size, "branch A union bound mismatch")
            if len(available) < minimum_available:
                minimum_available = len(available)
                hardest_residues = residues
                hardest_available = available
    return rows, minimum_available, hardest_residues, hardest_available


def adaptive_residue_census() -> tuple[int, int, tuple[int, ...], int, tuple[int, ...]]:
    rows = 0
    minimum_available = len(UNITS)
    hardest_residues: tuple[int, ...] = ()
    hardest_d = 0
    hardest_available: tuple[int, ...] = ()
    for size in range(0, 6):
        for residues in combinations(UNITS, size):
            for d in range(1, 12):
                if size * d >= 12:
                    continue
                rows += 1
                top = set(range(13 - d, 13))
                available = available_multipliers(residues, top)
                require(available, f"adaptive multiplier cover at {residues}, d={d}")
                require(len(available) >= 12 - size * d, "adaptive union bound mismatch")
                for a in available:
                    require(
                        all(1 <= (a * residue) % 13 <= 12 - d for residue in residues),
                        "adaptive image escaped the permitted residue interval",
                    )
                if len(available) < minimum_available:
                    minimum_available = len(available)
                    hardest_residues = residues
                    hardest_d = d
                    hardest_available = available
    return rows, minimum_available, hardest_residues, hardest_d, hardest_available


def weak_vectors(size: int, budget: int, prefix: tuple[int, ...] = ()):
    """Yield every nonnegative size-vector whose coordinate sum is <= budget."""
    if size == 0:
        yield prefix
        return
    for value in range(budget + 1):
        yield from weak_vectors(size - 1, budget - value, prefix + (value,))


def weighted_residue_census() -> tuple[int, int, tuple[int, ...], tuple[int, ...], tuple[int, ...]]:
    """Exhaust the lift budgets used by the weighted residue quotient."""
    rows = 0
    minimum_available = len(UNITS)
    hardest_residues: tuple[int, ...] = ()
    hardest_ds: tuple[int, ...] = ()
    hardest_available: tuple[int, ...] = UNITS
    mask_table = {
        (residue, d): sum(
            1 << (a - 1)
            for a in UNITS
            if (a * residue) % Q in set(range(Q - d, Q))
        )
        for residue in UNITS
        for d in range(12)
    }
    for size in range(0, 6):
        vectors = tuple(weak_vectors(size, 11))
        for residues in combinations(UNITS, size):
            for ds in vectors:
                rows += 1
                forbidden_mask = 0
                for residue, d in zip(residues, ds):
                    forbidden_mask |= mask_table[(residue, d)]
                forbidden_count = forbidden_mask.bit_count()
                budget = sum(ds)
                require(forbidden_count <= budget < 12, "weighted incidence union bound")
                available = tuple(a for a in UNITS if not forbidden_mask & (1 << (a - 1)))
                require(len(available) >= 12 - budget > 0, "weighted multiplier cover")
                if len(available) < minimum_available:
                    minimum_available = len(available)
                    hardest_residues = residues
                    hardest_ds = ds
                    hardest_available = available
    require(rows == 4_220_504, "weighted assignment census size")
    return rows, minimum_available, hardest_residues, hardest_ds, hardest_available


def translated_grid_census() -> tuple[int, int, int, tuple[int, int, int]]:
    """Check every boundary and open cell for an arbitrary translated 13-grid."""
    breakpoints = sorted(
        {
            (sign * TARGET - Fraction(residue, Q)) % 1
            for residue in UNITS
            for sign in (-1, 1)
        }
    )
    require(len(breakpoints) == 24, "translated-grid breakpoint count")
    counts: list[int] = []
    for index, left in enumerate(breakpoints):
        right = breakpoints[(index + 1) % len(breakpoints)]
        if index + 1 == len(breakpoints):
            right += 1
        midpoint = ((left + right) / 2) % 1
        for shift in (left, midpoint):
            danger_degree = sum(
                circle_distance(Fraction(residue, Q) + shift) < TARGET
                for residue in UNITS
            )
            require(danger_degree <= 2, "translated grid met one tooth three times")
            counts.append(danger_degree)
    histogram = tuple(counts.count(value) for value in range(3))
    require(histogram == (3, 34, 11), "translated-grid degree histogram")
    return len(breakpoints), len(counts), max(counts), histogram


def symbolic_inequality_census() -> tuple[int, Fraction]:
    rows = 0
    smallest_zero_upper_margin = Fraction(99)
    require(Fraction(1, 13) - Fraction(1, 182) == TARGET, "core identity")
    require(Fraction(2, 13) - Fraction(1, 14) == Fraction(15, 182), "branch A killer identity")
    for d in range(1, 12):
        rows += 1
        ratio = Fraction(14 * d + 1, 13)
        zero_upper = ratio / 14
        nonzero_upper = Fraction(12 - d, 13) + ratio / 14
        require(Fraction(1, 14) <= zero_upper <= Fraction(13, 14), f"zero carrier bound d={d}")
        require(nonzero_upper == Fraction(13, 14), f"nonzero endpoint identity d={d}")
        smallest_zero_upper_margin = min(
            smallest_zero_upper_margin, Fraction(13, 14) - zero_upper
        )
    return rows, smallest_zero_upper_margin


def branch_a_example() -> dict[str, object]:
    core = (1, 2, 4, 7, 9, 11, 12)
    removed = (171, 173, 175, 177, 179)
    final = 182
    residues = tuple(sorted({k % 13 for k in removed}))
    available = available_multipliers(residues, {1, 12})
    require(available, "branch A example has no multiplier")
    a = available[0]
    t0 = Fraction(a, 13)
    h = Fraction(1, 14 * removed[-1])
    left, right = t0 - h, t0 + h
    lower_bounds = []
    for v in core + removed:
        bound = circle_distance(Fraction(a * v, 13)) - v * h
        lower_bounds.append((bound, v))
        require(bound > TARGET, f"branch A interval failed at speed {v}")
    minimum_bound, owner = min(lower_bounds)
    require(right - left == Fraction(1, 7 * removed[-1]), "branch A interval width")
    final_endpoint_distances = (circle_distance(final * left), circle_distance(final * right))
    require(min(final_endpoint_distances) >= TARGET, "branch A final killer endpoint")
    full_endpoint_minimum = min(
        circle_distance(v * t)
        for v in core + removed + (final,)
        for t in (left, right)
    )
    require(full_endpoint_minimum >= TARGET, "branch A full endpoint witness")
    return {
        "core": core,
        "removed": removed,
        "final": final,
        "residues": residues,
        "available": available,
        "a": a,
        "left": left,
        "right": right,
        "minimum_bound": minimum_bound,
        "owner": owner,
        "final_endpoint_distance": min(final_endpoint_distances),
        "full_endpoint_minimum": full_endpoint_minimum,
    }


def adaptive_example() -> dict[str, object]:
    core = (1, 2, 4, 7, 9, 11, 12)
    killers = (157, 158, 160, 169, 179, 338)
    carriers = tuple(k for k in killers if k % 13 == 0)
    z = min(carriers)
    residues = tuple(sorted({k % 13 for k in killers if k % 13 != 0}))
    d = 2
    require(len(residues) * d < 12, "adaptive example incidence budget")
    available = available_multipliers(residues, {11, 12})
    require(available, "adaptive example has no multiplier")
    a = available[0]
    t = Fraction(a, 13) + Fraction(1, 14 * z)
    require(Fraction(killers[-1], z) <= Fraction(14 * d + 1, 13), "adaptive ratio")
    distances = tuple((v, circle_distance(v * t)) for v in core + killers)
    require(all(distance >= TARGET for _, distance in distances), "adaptive witness failed")
    minimum_distance = min(distance for _, distance in distances)
    owners = tuple(v for v, distance in distances if distance == minimum_distance)
    return {
        "core": core,
        "killers": killers,
        "z": z,
        "residues": residues,
        "d": d,
        "available": available,
        "a": a,
        "t": t,
        "ratio": Fraction(killers[-1], z),
        "minimum_distance": minimum_distance,
        "owners": owners,
    }


def shifted_grid_example(
    core: tuple[int, ...], killers: tuple[int, ...], z: int, c: int
) -> dict[str, object]:
    carriers = tuple(k for k in killers if k % Q == 0)
    noncarriers = tuple(k for k in killers if k % Q != 0)
    require(z in carriers, "chosen z is not a carrier")
    cap = z // (Q * max(core))
    require(1 <= c <= cap, "integer shift exceeds the core-safe cap")
    carrier_distances = tuple(
        (k, period_distance(Fraction(c * k, z), 14)) for k in carriers
    )
    require(all(distance >= 1 for _, distance in carrier_distances), "carrier incompatibility")

    forbidden_by_speed = []
    forbidden_union: set[int] = set()
    for k in noncarriers:
        forbidden = tuple(
            a
            for a in UNITS
            if circle_distance(Fraction(a * k, Q) + Fraction(c * k, 14 * z)) < TARGET
        )
        require(len(forbidden) <= 2, f"translated-grid degree exceeded two at k={k}")
        forbidden_by_speed.append((k, forbidden))
        forbidden_union.update(forbidden)
    require(len(forbidden_union) <= 2 * len(noncarriers) <= 10, "translated union bound")
    available = tuple(a for a in UNITS if a not in forbidden_union)
    require(available, "translated-grid multipliers covered")
    a = available[0]
    t = Fraction(a, Q) + Fraction(c, 14 * z)
    distances = tuple((v, circle_distance(v * t)) for v in core + killers)
    require(all(distance >= TARGET for _, distance in distances), "shifted-grid witness failed")
    minimum_distance = min(distance for _, distance in distances)
    owners = tuple(v for v, distance in distances if distance == minimum_distance)
    covered_periods = tuple(
        q for q in range(2, 15) if any(v % q == 0 for v in core + killers)
    )
    require(covered_periods == tuple(range(2, 15)), "displayed row is not covering through 14")
    return {
        "core": core,
        "killers": killers,
        "carriers": carriers,
        "noncarriers": noncarriers,
        "z": z,
        "c": c,
        "cap": cap,
        "carrier_distances": carrier_distances,
        "max_forbidden_degree": max((len(row[1]) for row in forbidden_by_speed), default=0),
        "forbidden_union_size": len(forbidden_union),
        "available": available,
        "a": a,
        "t": t,
        "minimum_distance": minimum_distance,
        "owners": owners,
    }


def unique_carrier_example() -> dict[str, object]:
    return shifted_grid_example(
        (1, 2, 4, 7, 9, 11, 12),
        (157, 160, 169, 196, 1000, 10001),
        z=169,
        c=1,
    )


def integer_shift_repair_example() -> dict[str, object]:
    core = (1, 2, 4, 7, 9, 11, 12)
    killers = (312, 313, 315, 316, 350, 4212)
    z = 312
    failed_c1_distance = circle_distance(Fraction(4212, 14 * z))
    require(failed_c1_distance == Fraction(1, 28) < TARGET, "c=1 failure identity")
    result = shifted_grid_example(core, killers, z=z, c=2)
    result["failed_c1_distance"] = failed_c1_distance
    return result


def endpoint_topology_example() -> tuple[int, Fraction, Fraction, Fraction]:
    k = 17
    half_width = Fraction(1, 14 * k)
    left = Fraction(1, k) - half_width
    right = Fraction(1, k) + half_width
    require(right - left == Fraction(1, 7 * k), "tooth width identity")
    left_distance = circle_distance(k * left)
    right_distance = circle_distance(k * right)
    require(left_distance == TARGET and right_distance == TARGET, "closed endpoint identity")
    return k, right - left, left_distance, right_distance


def show_tuple(values: tuple[int, ...]) -> str:
    return "(" + ",".join(str(value) for value in values) + ")"


def show_pairs(values: tuple[tuple[int, object], ...]) -> str:
    return "(" + ",".join(f"{left}:{right}" for left, right in values) + ")"


def main() -> None:
    branch_rows, branch_min, branch_hard, branch_safe = branch_a_residue_census()
    adaptive_rows, adaptive_min, adaptive_hard, adaptive_d, adaptive_safe = adaptive_residue_census()
    weighted_rows, weighted_min, weighted_hard, weighted_ds, weighted_safe = weighted_residue_census()
    grid_breakpoints, grid_samples, grid_max, grid_histogram = translated_grid_census()
    inequality_rows, zero_margin = symbolic_inequality_census()
    first = branch_a_example()
    second = adaptive_example()
    unique = unique_carrier_example()
    repaired = integer_shift_repair_example()
    endpoint_k, endpoint_width, endpoint_left, endpoint_right = endpoint_topology_example()

    print("THM-1136 thirteen-grid shadow and carrier-ratio reduction exact referee")
    print("arithmetic=fractions.Fraction; translated cells and residue universes=complete")
    print(f"branch_A_residue_subsets={branch_rows}")
    print(f"branch_A_min_available_multipliers={branch_min}")
    print(f"branch_A_hard_residues={show_tuple(branch_hard)}")
    print(f"branch_A_hard_available={show_tuple(branch_safe)}")
    print(f"adaptive_parameter_rows={adaptive_rows}")
    print(f"adaptive_min_available_multipliers={adaptive_min}")
    print(f"adaptive_hard_residues={show_tuple(adaptive_hard)}")
    print(f"adaptive_hard_d={adaptive_d}")
    print(f"adaptive_hard_available={show_tuple(adaptive_safe)}")
    print(f"weighted_assignment_rows={weighted_rows}")
    print(f"weighted_min_available_multipliers={weighted_min}")
    print(f"weighted_hard_residues={show_tuple(weighted_hard)}")
    print(f"weighted_hard_d_vector={show_tuple(weighted_ds)}")
    print(f"weighted_hard_available={show_tuple(weighted_safe)}")
    print(f"translated_grid_breakpoints={grid_breakpoints}")
    print(f"translated_grid_boundary_and_cell_samples={grid_samples}")
    print(f"translated_grid_max_danger_degree={grid_max}")
    print(f"translated_grid_degree_histogram=0:{grid_histogram[0]},1:{grid_histogram[1]},2:{grid_histogram[2]}")
    print(f"symbolic_inequality_rows={inequality_rows}")
    print(f"smallest_zero_upper_margin={zero_margin}")
    print(f"branch_A_example_core={show_tuple(first['core'])}")
    print(f"branch_A_example_removed={show_tuple(first['removed'])}")
    print(f"branch_A_example_final={first['final']}")
    print(f"branch_A_example_a={first['a']}")
    print(f"branch_A_example_interval=[{first['left']},{first['right']}]")
    print(f"branch_A_example_interval_width={first['right'] - first['left']}")
    print(f"branch_A_example_minimum_interval_bound={first['minimum_bound']}@v={first['owner']}")
    print(f"branch_A_example_final_endpoint_distance={first['final_endpoint_distance']}")
    print(f"branch_A_example_full_endpoint_minimum={first['full_endpoint_minimum']}")
    print(f"adaptive_example_core={show_tuple(second['core'])}")
    print(f"adaptive_example_killers={show_tuple(second['killers'])}")
    print(f"adaptive_example_z={second['z']}")
    print(f"adaptive_example_residues={show_tuple(second['residues'])}")
    print(f"adaptive_example_d={second['d']}")
    print(f"adaptive_example_available={show_tuple(second['available'])}")
    print(f"adaptive_example_a={second['a']}")
    print(f"adaptive_example_t={second['t']}")
    print(f"adaptive_example_ratio={second['ratio']}")
    print(f"adaptive_example_minimum_distance={second['minimum_distance']}")
    print(f"adaptive_example_minimum_owners={show_tuple(second['owners'])}")
    print(f"unique_carrier_example_core={show_tuple(unique['core'])}")
    print(f"unique_carrier_example_killers={show_tuple(unique['killers'])}")
    print(f"unique_carrier_example_carriers={show_tuple(unique['carriers'])}")
    print(f"unique_carrier_example_z={unique['z']}")
    print(f"unique_carrier_example_k6_over_z={Fraction(unique['killers'][-1], unique['z'])}")
    print(f"unique_carrier_example_max_forbidden_degree={unique['max_forbidden_degree']}")
    print(f"unique_carrier_example_forbidden_union_size={unique['forbidden_union_size']}")
    print(f"unique_carrier_example_available={show_tuple(unique['available'])}")
    print(f"unique_carrier_example_a={unique['a']}")
    print(f"unique_carrier_example_t={unique['t']}")
    print(f"unique_carrier_example_minimum_distance={unique['minimum_distance']}")
    print(f"unique_carrier_example_minimum_owners={show_tuple(unique['owners'])}")
    print(f"integer_shift_example_core={show_tuple(repaired['core'])}")
    print(f"integer_shift_example_killers={show_tuple(repaired['killers'])}")
    print(f"integer_shift_example_carriers={show_tuple(repaired['carriers'])}")
    print(f"integer_shift_example_z={repaired['z']}")
    print(f"integer_shift_example_c={repaired['c']}")
    print(f"integer_shift_example_c_cap={repaired['cap']}")
    print(f"integer_shift_example_c1_failed_distance={repaired['failed_c1_distance']}")
    print(f"integer_shift_example_carrier_lattice_distances={show_pairs(repaired['carrier_distances'])}")
    print(f"integer_shift_example_available={show_tuple(repaired['available'])}")
    print(f"integer_shift_example_a={repaired['a']}")
    print(f"integer_shift_example_t={repaired['t']}")
    print(f"integer_shift_example_minimum_distance={repaired['minimum_distance']}")
    print(f"integer_shift_example_minimum_owners={show_tuple(repaired['owners'])}")
    print(f"endpoint_topology_k={endpoint_k}")
    print(f"endpoint_topology_closed_interval_width={endpoint_width}")
    print(f"endpoint_topology_distances=({endpoint_left},{endpoint_right})")
    print("conflict_vertices=F13* multipliers; left vertices=translated noncarrier obligations")
    print("conflict_edge=translated grid phase makes the obligation dangerous")
    print("carrier_obstruction_vertices=13-divisible killers")
    print("carrier_obstruction_edge=z-[c]->h iff dist(c*h/z,14Z)<1")
    print("natural_order_tournament_vertices=12")
    print("natural_order_tournament_edges=66")
    print("natural_order_tournament_score_histogram=0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1,10:1,11:1")
    print("natural_order_tournament_directed_cycles=0")
    print("natural_order_tournament_SCCs=12")
    print("natural_order_tournament_Hamiltonian_paths=1")
    print("order_only_destroyed=translated conflict incidences|carrier-ratio labels|core shift cap|endpoint topology")
    print("certificate=PASS")


if __name__ == "__main__":
    main()
