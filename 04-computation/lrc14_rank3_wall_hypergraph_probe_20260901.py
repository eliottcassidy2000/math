#!/usr/bin/env python3
"""Exact rank-at-most-three wall-hypergraph probe for the THM-4326 pool.

This is deliberately a scoped research probe, not an LRC(14) certificate.
It checks the rank-three inclusion--exclusion identity, reproduces the frozen
THM-4326 hostile rank-two value, tests selected third-outsider controls or an
optional complete typed scan through a requested finite endpoint, and audits
full-wall restoration whenever the rank-three truncation falls below target.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction
from itertools import combinations
from math import gcd


POOL = (
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290,
)
HOSTILE_BODY = frozenset((80, 85, 88, 95, 143, 145, 168, 193, 240))
HOSTILE_MASK = 0x031C7400
PAIR = (50, 70)
THIRD_OUTSIDERS = (
    1, 2, 3, 4, 5, 6, 7, 9, 11, 12, 13, 14, 17, 18, 19,
    45, 55, 65, 90, 110, 130, 220, 260, 590,
)
TARGET = Fraction(4, 63)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def lcm(left: int, right: int) -> int:
    return left // gcd(left, right) * right


def grid_for(outsiders: tuple[int, ...]) -> int:
    grid = 1
    for speed in POOL + outsiders:
        grid = lcm(grid, 14 * speed)
    return grid


def safe_midpoint(speed: int, grid: int, left: int, right: int) -> bool:
    residue = speed * (left + right) % (2 * grid)
    return grid <= 7 * residue <= 13 * grid


def mask_of(labels: frozenset[int]) -> int:
    return sum(1 << index for index, speed in enumerate(POOL) if speed in labels)


def walls_for(grid: int, speeds: tuple[int, ...]) -> list[int]:
    walls = {0, grid}
    for speed in speeds:
        unit, remainder = divmod(grid, 14 * speed)
        require(remainder == 0, "nonintegral wall unit")
        for tooth in range(speed):
            walls.add((14 * tooth + 1) * unit)
            walls.add((14 * tooth + 13) * unit)
    return sorted(walls)


def wall_cells(
    outsiders: tuple[int, ...],
) -> tuple[int, list[tuple[int, int, int]]]:
    require(len(set(outsiders)) == len(outsiders), "outsiders must be distinct")
    require(all(speed > 0 and speed not in POOL for speed in outsiders),
            "outsider typing failed")
    grid = grid_for(outsiders)
    walls = walls_for(grid, POOL + outsiders)
    cells: list[tuple[int, int, int]] = []
    for left, right in zip(walls, walls[1:]):
        if not all(safe_midpoint(speed, grid, left, right)
                   for speed in outsiders):
            continue
        failure = 0
        for index, speed in enumerate(POOL):
            if not safe_midpoint(speed, grid, left, right):
                failure |= 1 << index
        cells.append((left, right, failure))
    return grid, cells


def weights_up_to(
    cells: list[tuple[int, int, int]], rank: int,
) -> dict[int, int]:
    weights: dict[int, int] = defaultdict(int)
    for left, right, failure in cells:
        if failure.bit_count() <= rank:
            weights[failure] += right - left
    return dict(weights)


def retained_direct(weights: dict[int, int], body: int) -> int:
    return sum(width for failure, width in weights.items()
               if failure & body == 0)


def full_mass(cells: list[tuple[int, int, int]], body: int) -> int:
    return sum(right - left for left, right, failure in cells
               if failure & body == 0)


def rank3_statistics(
    weights: dict[int, int],
) -> tuple[int, list[int], list[list[int]], dict[int, int]]:
    total = sum(weights.values())
    degree = [0] * len(POOL)
    codegree = [[0] * len(POOL) for _ in POOL]
    triples: dict[int, int] = {}
    for failure, width in weights.items():
        vertices = tuple(index for index in range(len(POOL))
                         if failure >> index & 1)
        for vertex in vertices:
            degree[vertex] += width
        for left, right in combinations(vertices, 2):
            codegree[left][right] += width
            codegree[right][left] += width
        if len(vertices) == 3:
            triples[failure] = width
    return total, degree, codegree, triples


def retained_inclusion_exclusion(weights: dict[int, int], body: int) -> int:
    total, degree, codegree, triples = rank3_statistics(weights)
    vertices = tuple(index for index in range(len(POOL)) if body >> index & 1)
    value = total - sum(degree[index] for index in vertices)
    value += sum(codegree[left][right]
                 for left, right in combinations(vertices, 2))
    value -= sum(triples.get(sum(1 << vertex for vertex in triple), 0)
                 for triple in combinations(vertices, 3))
    return value


def degree_lower_bound(weights: dict[int, int], body_size: int) -> int:
    total, degree, _, _ = rank3_statistics(weights)
    return total - sum(sorted(degree, reverse=True)[:body_size])


def exact_maximum_coverage(
    weights: dict[int, int], body_size: int,
) -> tuple[int, int, int, int]:
    """Exact branch-and-bound for the weighted maximum-coverage dual.

    Returns (maximum covered width, least maximizing mask, nodes, prunes).
    The bound is the sum of the largest current marginals.  Every newly
    covered hyperedge occurs in at least one chosen marginal, so this may
    overcount overlaps but can never undercount a completion.
    """
    edges = tuple((failure, width) for failure, width in weights.items()
                  if failure != 0 and width != 0)

    def greedy() -> tuple[int, int]:
        candidates = (1 << len(POOL)) - 1
        uncovered = edges
        covered = 0
        body = 0
        for _ in range(body_size):
            marginal = [0] * len(POOL)
            for failure, width in uncovered:
                active = failure & candidates
                while active:
                    bit = active & -active
                    marginal[bit.bit_length() - 1] += width
                    active -= bit
            vertex = max((index for index in range(len(POOL))
                          if candidates >> index & 1),
                         key=lambda index: (marginal[index], -index))
            bit = 1 << vertex
            covered += marginal[vertex]
            body |= bit
            candidates ^= bit
            uncovered = tuple(edge for edge in uncovered if edge[0] & bit == 0)
        return covered, body

    incumbent, incumbent_body = greedy()
    nodes = 0
    prunes = 0

    def search(
        candidates: int,
        chosen: int,
        need: int,
        uncovered: tuple[tuple[int, int], ...],
        covered: int,
    ) -> None:
        nonlocal incumbent, incumbent_body, nodes, prunes
        nodes += 1
        available = candidates.bit_count()
        require(available >= need, "branch cardinality invariant")
        if need == 0:
            if covered > incumbent or (covered == incumbent and chosen < incumbent_body):
                incumbent, incumbent_body = covered, chosen
            return
        if available == need:
            added = sum(width for failure, width in uncovered
                        if failure & candidates)
            value = covered + added
            body = chosen | candidates
            if value > incumbent or (value == incumbent and body < incumbent_body):
                incumbent, incumbent_body = value, body
            return

        marginal = [0] * len(POOL)
        for failure, width in uncovered:
            active = failure & candidates
            while active:
                bit = active & -active
                marginal[bit.bit_length() - 1] += width
                active -= bit
        ordered = sorted((marginal[index], index)
                         for index in range(len(POOL)) if candidates >> index & 1)
        upper = covered + sum(value for value, _ in ordered[-need:])
        if upper < incumbent:
            prunes += 1
            return

        # Include the highest-marginal vertex first.  On ties take the lower
        # label index, only to keep the transcript deterministic.
        pivot = max((index for index in range(len(POOL))
                     if candidates >> index & 1),
                    key=lambda index: (marginal[index], -index))
        bit = 1 << pivot
        remaining = candidates ^ bit
        included_uncovered = tuple(edge for edge in uncovered
                                   if edge[0] & bit == 0)
        search(remaining, chosen | bit, need - 1, included_uncovered,
               covered + marginal[pivot])
        if remaining.bit_count() >= need:
            search(remaining, chosen, need, uncovered, covered)

    search((1 << len(POOL)) - 1, 0, body_size, edges, 0)
    return incumbent, incumbent_body, nodes, prunes


def subthreshold_bodies(
    weights: dict[int, int], grid: int, body_size: int,
) -> tuple[list[tuple[int, int]], int, int]:
    """Enumerate every body with retained rank-three mass below 4D/63.

    The returned pairs are (body mask, retained mass).  A branch is discarded
    only when the same current-marginal upper bound proves that even its best
    possible coverage cannot cross the strict target.
    """
    edges = tuple((failure, width) for failure, width in weights.items()
                  if failure != 0 and width != 0)
    total = sum(weights.values())
    hits: list[tuple[int, int]] = []
    nodes = 0
    prunes = 0

    def record(body: int, covered: int) -> None:
        retained = total - covered
        if 63 * retained < 4 * grid:
            hits.append((body, retained))

    def search(
        candidates: int,
        chosen: int,
        need: int,
        uncovered: tuple[tuple[int, int], ...],
        covered: int,
    ) -> None:
        nonlocal nodes, prunes
        nodes += 1
        available = candidates.bit_count()
        require(available >= need, "threshold branch cardinality invariant")
        if need == 0:
            record(chosen, covered)
            return
        if available == need:
            added = sum(width for failure, width in uncovered
                        if failure & candidates)
            record(chosen | candidates, covered + added)
            return

        marginal = [0] * len(POOL)
        for failure, width in uncovered:
            active = failure & candidates
            while active:
                bit = active & -active
                marginal[bit.bit_length() - 1] += width
                active -= bit
        ordered = sorted((marginal[index], index)
                         for index in range(len(POOL)) if candidates >> index & 1)
        upper = covered + sum(value for value, _ in ordered[-need:])
        if 63 * (total - upper) >= 4 * grid:
            prunes += 1
            return

        pivot = max((index for index in range(len(POOL))
                     if candidates >> index & 1),
                    key=lambda index: (marginal[index], -index))
        bit = 1 << pivot
        remaining = candidates ^ bit
        included_uncovered = tuple(edge for edge in uncovered
                                   if edge[0] & bit == 0)
        search(remaining, chosen | bit, need - 1, included_uncovered,
               covered + marginal[pivot])
        if remaining.bit_count() >= need:
            search(remaining, chosen, need, uncovered, covered)

    search((1 << len(POOL)) - 1, 0, body_size, edges, 0)
    hits.sort()
    return hits, nodes, prunes


def split_witness(
    base_cells: list[tuple[int, int, int]], base_grid: int, newcomer: int,
) -> tuple[int, tuple[int, int], tuple[int, int]]:
    """Find one aggregate atom split by a proposed next outsider.

    The newcomer need not have walls on the base grid: midpoint safety is
    constant only after refinement, so refine by its walls first and then
    look for a common rank<=3 mask with both safe and unsafe subcells.
    """
    refined_grid = lcm(base_grid, 14 * newcomer)
    scale = refined_grid // base_grid
    base_boundaries = {left * scale for left, _, _ in base_cells}
    base_boundaries.update(right * scale for _, right, _ in base_cells)
    refined_walls = sorted(base_boundaries | set(walls_for(refined_grid, (newcomer,))))

    # Reclassify against the base outsider condition and pool.  This avoids
    # assuming that a newcomer wall respects an old open-cell boundary.
    buckets: dict[int, dict[bool, tuple[int, int]]] = defaultdict(dict)
    for left, right in zip(refined_walls, refined_walls[1:]):
        if not all(safe_midpoint(speed, refined_grid, left, right)
                   for speed in PAIR):
            continue
        failure = 0
        for index, speed in enumerate(POOL):
            if not safe_midpoint(speed, refined_grid, left, right):
                failure |= 1 << index
        if failure.bit_count() > 3:
            continue
        is_safe = safe_midpoint(newcomer, refined_grid, left, right)
        buckets[failure].setdefault(is_safe, (left, right))
        if True in buckets[failure] and False in buckets[failure]:
            return failure, buckets[failure][True], buckets[failure][False]
    raise AssertionError("no adaptive split witness found")


def check_abstract_identity() -> None:
    # Includes every rank 0--3 atom and tests every body of every size.
    weights = {mask: (17 * mask + 11) for mask in range(1 << 5)
               if mask.bit_count() <= 3}
    for body in range(1 << 5):
        direct = retained_direct(weights, body)
        # The production formula has 30 vertices; pad the five-vertex model.
        embedded = dict(weights)
        require(retained_inclusion_exclusion(embedded, body) == direct,
                "abstract inclusion--exclusion identity failed")


def check_small_optimizers() -> None:
    # Independent flat enumeration controls both branch engines on a small
    # nontrivial weighted 3-hypergraph embedded in the first eight vertices.
    weights = {
        mask: 5 + 13 * mask + 7 * mask.bit_count()
        for mask in range(1 << 8) if mask.bit_count() <= 3
    }
    flat: list[tuple[int, int]] = []
    for vertices in combinations(range(len(POOL)), 3):
        body = sum(1 << vertex for vertex in vertices)
        flat.append((retained_direct(weights, body), body))
    flat_minimum = min(flat)
    coverage, body, _, _ = exact_maximum_coverage(weights, 3)
    require((sum(weights.values()) - coverage, body) == flat_minimum,
            "small flat/branch optimizer disagreement")

    ordered_values = sorted(value for value, _ in flat)
    cutoff = ordered_values[len(ordered_values) // 2]
    grid = (63 * cutoff + 3) // 4
    expected = sorted((body, value) for value, body in flat
                      if 63 * value < 4 * grid)
    actual, _, _ = subthreshold_bodies(weights, grid, 3)
    require(actual == expected and 0 < len(actual) < len(flat),
            "small flat/threshold branch disagreement")


def ratio_text(numerator: int, denominator: int) -> str:
    value = Fraction(numerator, denominator)
    return f"{value.numerator}/{value.denominator}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--scan-through",
        type=int,
        default=0,
        help="also scan every typed third outsider from 1 through this endpoint",
    )
    args = parser.parse_args()
    require(args.scan_through >= 0, "negative scan endpoint")
    require(mask_of(HOSTILE_BODY) == HOSTILE_MASK, "hostile body mask changed")
    check_abstract_identity()
    check_small_optimizers()

    pair_grid, pair_cells = wall_cells(PAIR)
    pair_rank2 = weights_up_to(pair_cells, 2)
    pair_l2 = retained_direct(pair_rank2, HOSTILE_MASK)
    require(pair_grid == 91_205_797_082_400, "THM-4326 hostile grid changed")
    require(pair_l2 == 5_794_739_949_188, "THM-4326 hostile L2 changed")
    pair_coverage, pair_body, pair_nodes, pair_prunes = exact_maximum_coverage(
        pair_rank2, 9
    )
    require(sum(pair_rank2.values()) - pair_coverage == pair_l2,
            "rank-two optimizer minimum changed")
    require(pair_body == HOSTILE_MASK, "rank-two optimizer body changed")

    pair_rank3 = weights_up_to(pair_cells, 3)
    pair_l3 = retained_direct(pair_rank3, HOSTILE_MASK)
    require(pair_l3 == retained_inclusion_exclusion(pair_rank3, HOSTILE_MASK),
            "pair rank-three identity failed")
    require(pair_l3 >= pair_l2, "rank monotonicity failed")
    pair_rank3_coverage, pair_rank3_body, pair_rank3_nodes, pair_rank3_prunes = (
        exact_maximum_coverage(pair_rank3, 9)
    )
    pair_rank3_minimum = sum(pair_rank3.values()) - pair_rank3_coverage
    require(pair_rank3_minimum == retained_direct(pair_rank3, pair_rank3_body),
            "pair rank-three minimum/direct disagreement")
    retained_rank3_cells = sum(
        failure.bit_count() <= 3 for _, _, failure in pair_cells
    )
    cofinal_margin = (
        Fraction(6, 7) * Fraction(pair_rank3_minimum, pair_grid) - TARGET
    )
    require(cofinal_margin > 0, "rank-three cofinal margin vanished")
    cutoff_fraction = Fraction(6 * retained_rank3_cells, 49) / cofinal_margin
    cofinal_cutoff = (
        cutoff_fraction.numerator + cutoff_fraction.denominator - 1
    ) // cutoff_fraction.denominator
    cofinal_lower = (
        Fraction(6, 7) * Fraction(pair_rank3_minimum, pair_grid)
        - Fraction(6 * retained_rank3_cells, 49 * cofinal_cutoff)
    )
    cofinal_previous = (
        Fraction(6, 7) * Fraction(pair_rank3_minimum, pair_grid)
        - Fraction(6 * retained_rank3_cells, 49 * (cofinal_cutoff - 1))
    )
    require(cofinal_cutoff == 18_406 and cofinal_lower > TARGET
            and cofinal_previous < TARGET,
            "cofinal component-count certificate changed")

    print("LRC14_RANK3_WALL_HYPERGRAPH_PROBE_20260901")
    print(
        "SCOPE fixed THM-4326 pool; selected controls plus complete bounded "
        "fixed-pair third scan; no arbitrary-pair or unbounded triple universe"
    )
    print("ABSTRACT_IDENTITY all rank<=3 atoms on five vertices, all 32 bodies PASS")
    print("SMALL_FLAT_CONTROL weighted rank<=3 atoms on eight vertices, all C(30,3) bodies PASS")
    print(
        "PAIR_CONTROL outsiders=50,70 body=031c7400 "
        f"grid={pair_grid} L2={pair_l2} L3={pair_l3} "
        f"L3_ratio={ratio_text(pair_l3, pair_grid)} "
        f"optimizer_nodes={pair_nodes} optimizer_prunes={pair_prunes}"
    )
    print(
        f"PAIR_RANK3_MIN L3={pair_rank3_minimum} body={pair_rank3_body:08x} "
        f"ratio={Fraction(pair_rank3_minimum, pair_grid)} "
        f"nodes={pair_rank3_nodes} prunes={pair_rank3_prunes}"
    )
    print(
        f"COFINAL_COMPONENT_BOUND rank3_cells={retained_rank3_cells} "
        f"margin={cofinal_margin} cutoff={cofinal_cutoff} "
        f"lower_at_cutoff={cofinal_lower} lower_at_previous={cofinal_previous}"
    )

    positive = 0
    below_target = 0
    for newcomer in THIRD_OUTSIDERS:
        outsiders = PAIR + (newcomer,)
        grid, cells = wall_cells(outsiders)
        weights = weights_up_to(cells, 3)
        hostile_retained = retained_direct(weights, HOSTILE_MASK)
        via_ie = retained_inclusion_exclusion(weights, HOSTILE_MASK)
        coverage, minimizing_body, nodes, prunes = exact_maximum_coverage(weights, 9)
        retained = sum(weights.values()) - coverage
        full = full_mass(cells, minimizing_body)
        coarse = degree_lower_bound(weights, 9)
        require(hostile_retained == via_ie,
                f"rank-three hostile identity failed for {outsiders}")
        require(retained == retained_direct(weights, minimizing_body),
                f"optimizer/direct disagreement for {outsiders}")
        require(retained == retained_inclusion_exclusion(weights, minimizing_body),
                f"optimizer/inclusion-exclusion disagreement for {outsiders}")
        require(minimizing_body.bit_count() == 9 and coarse <= retained <= full,
                f"bound chain failed for {outsiders}")
        positive += retained > 0
        below_target += Fraction(retained, grid) < TARGET
        sign = "ABOVE" if Fraction(retained, grid) >= TARGET else "BELOW"
        print(
            f"TRIPLE outsiders=50,70,{newcomer} grid={grid} "
            f"degree_bound={coarse} min_L3={retained} "
            f"min_body={minimizing_body:08x} full_on_min_body={full} "
            f"old_hostile_L3={hostile_retained} "
            f"min_ratio={ratio_text(retained, grid)} target={sign} "
            f"nodes={nodes} prunes={prunes}"
        )

    failure, safe_cell, unsafe_cell = split_witness(pair_cells, pair_grid, 1)
    print(
        f"ADAPTIVE_SPLIT newcomer=1 rank={failure.bit_count()} mask={failure:08x} "
        f"safe_cell={safe_cell[0]}/{lcm(pair_grid, 14)},{safe_cell[1]}/{lcm(pair_grid, 14)} "
        f"unsafe_cell={unsafe_cell[0]}/{lcm(pair_grid, 14)},{unsafe_cell[1]}/{lcm(pair_grid, 14)}"
    )
    print(
        f"SUMMARY triples={len(THIRD_OUTSIDERS)} positive_L3={positive} "
        f"below_4_over_63={below_target}"
    )

    if args.scan_through:
        scan_rows: list[tuple[Fraction, int, int, int, int, int]] = []
        scan_below: list[tuple[int, Fraction, int, Fraction]] = []
        equality = 0
        for newcomer in range(1, args.scan_through + 1):
            if newcomer in POOL or newcomer in PAIR:
                continue
            grid, cells = wall_cells(PAIR + (newcomer,))
            weights = weights_up_to(cells, 3)
            coverage, body, nodes, _ = exact_maximum_coverage(weights, 9)
            retained = sum(weights.values()) - coverage
            ratio = Fraction(retained, grid)
            scan_rows.append((ratio, newcomer, retained, grid, body, nodes))
            equality += ratio == TARGET
            if ratio < TARGET:
                full = Fraction(full_mass(cells, body), grid)
                scan_below.append((newcomer, ratio, body, full))
        require(scan_rows, "empty scan universe")
        minimum = min(scan_rows)
        print(
            f"SCAN endpoint={args.scan_through} universe={len(scan_rows)} "
            f"above={sum(row[0] > TARGET for row in scan_rows)} "
            f"equal={equality} below={len(scan_below)}"
        )
        print(
            f"SCAN_MIN third={minimum[1]} ratio={minimum[0]} "
            f"body={minimum[4]:08x} L3={minimum[2]} grid={minimum[3]} "
            f"nodes={minimum[5]}"
        )
        for newcomer, ratio, body, full in scan_below:
            print(
                f"SCAN_BELOW third={newcomer} body={body:08x} "
                f"L3_ratio={ratio} full_on_L3_min_body={full}"
            )

        total_subthreshold = 0
        global_full_minimum: tuple[Fraction, int, int, int, int] | None = None
        for newcomer, _, _, _ in scan_below:
            grid, cells = wall_cells(PAIR + (newcomer,))
            weights = weights_up_to(cells, 3)
            hits, nodes, prunes = subthreshold_bodies(weights, grid, 9)
            require(hits, f"missing subthreshold witnesses at third={newcomer}")
            literal_rows: list[tuple[int, int, int]] = []
            for body, retained in hits:
                literal = full_mass(cells, body)
                require(literal >= retained,
                        f"full/rank-three monotonicity at third={newcomer}")
                require(63 * literal > 4 * grid,
                        f"literal full-mass hostile at third={newcomer}")
                literal_rows.append((literal, body, retained))
            literal, body, retained = min(literal_rows)
            total_subthreshold += len(hits)
            candidate = (Fraction(literal, grid), newcomer, body, literal, retained)
            if global_full_minimum is None or candidate < global_full_minimum:
                global_full_minimum = candidate
            print(
                f"SUBTHRESHOLD_AUDIT third={newcomer} bodies={len(hits)} "
                f"min_full_body={body:08x} min_full={literal} "
                f"min_full_ratio={Fraction(literal, grid)} "
                f"L3_on_min_full={retained} full_ticks={63 * literal - 4 * grid} "
                f"nodes={nodes} prunes={prunes}"
            )
        require(global_full_minimum is not None, "missing global literal minimum")
        print(
            f"SUBTHRESHOLD_SUMMARY bodies={total_subthreshold} all_full_strict=1 "
            f"global_min_third={global_full_minimum[1]} "
            f"global_min_body={global_full_minimum[2]:08x} "
            f"global_min_ratio={global_full_minimum[0]} "
            f"global_min_full={global_full_minimum[3]} "
            f"L3_on_global_min={global_full_minimum[4]}"
        )
    print("LOSS rank>=4 widths, cyclic cell address, owners, endpoints, and within-mask newcomer phase")
    print("VERDICT PASS_EXPLORATORY_NO_LRC14_CLAIM")


if __name__ == "__main__":
    main()
