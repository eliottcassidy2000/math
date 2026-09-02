#!/usr/bin/env python3
"""Independent exact verifier for the fixed-pair rank-three scan.

This verifier intentionally does not import or reuse the primary probe.  It
reconstructs the wall decomposition from the rational endpoints, classifies a
wall-cell by the integer part of 14*v*t, and sends the rank-at-most-three
minimum-coverage problem to OR-Tools CP-SAT.  Thus its optimization path is
independent of the primary custom branch-and-bound.

For every typed third outsider 1 <= s <= 769 it proves the minimum rank-three
retained mass over all nine-subsets of the fixed 30-label pool.  Whenever that
lower bound is below 4/63, it enumerates every such body and evaluates its full
safe mass from all wall cells (including failure masks of rank at least four).

All arithmetic used to build coefficients and audit solutions is integral.
CP-SAT is required to return OPTIMAL for minima and OPTIMAL after exhaustive
feasibility enumeration.
"""

from __future__ import annotations

import hashlib
import platform
import sys
import time
from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import lcm

from ortools.sat.python import cp_model
import ortools


POOL = (
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290,
)
FIXED_OUTSIDERS = (50, 70)
ENDPOINT = 769
BODY_SIZE = 9


def demand(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def typed_thirds() -> tuple[int, ...]:
    forbidden = set(POOL) | set(FIXED_OUTSIDERS)
    return tuple(speed for speed in range(1, ENDPOINT + 1)
                 if speed not in forbidden)


def midpoint_is_safe(speed: int, midpoint_numerator: int, grid: int) -> bool:
    """Classify t=midpoint_numerator/(2*grid) without residues.

    floor(14*v*t) modulo 14 is 1,...,12 precisely on the open safe
    interval (1/14,13/14) modulo one.  Every safety boundary was inserted in
    the wall list, so no tested midpoint is one of those boundaries.
    """
    scaled_numerator = 7 * speed * midpoint_numerator
    tooth, remainder = divmod(scaled_numerator, grid)
    demand(not (remainder == 0 and tooth % 14 in (1, 13)),
           "a wall-cell midpoint was itself a safety wall")
    return 1 <= tooth % 14 <= 12


def exact_wall_weights(
    outsiders: tuple[int, int, int],
) -> tuple[int, dict[int, int], dict[int, int], int]:
    """Return grid, rank<=3 weights, all-mask weights, and cell count."""
    speeds = POOL + outsiders
    grid = lcm(*(14 * speed for speed in speeds))
    walls = {0, grid}
    for speed in speeds:
        denominator = 14 * speed
        demand(grid % denominator == 0, "wall denominator does not divide grid")
        unit = grid // denominator
        for period in range(speed):
            walls.add((14 * period + 1) * unit)
            walls.add((14 * period + 13) * unit)
    ordered = sorted(walls)

    full: dict[int, int] = defaultdict(int)
    kept_cells = 0
    for left, right in zip(ordered, ordered[1:]):
        midpoint = left + right
        if not all(midpoint_is_safe(speed, midpoint, grid)
                   for speed in outsiders):
            continue
        failure = 0
        for index, speed in enumerate(POOL):
            if not midpoint_is_safe(speed, midpoint, grid):
                failure |= 1 << index
        full[failure] += right - left
        kept_cells += 1
    rank3 = {mask: width for mask, width in full.items()
             if mask.bit_count() <= 3}
    demand(sum(full.values()) <= grid, "safe outsider mass exceeds circle")
    return grid, rank3, dict(full), kept_cells


@dataclass(frozen=True)
class CoverageModel:
    model: cp_model.CpModel
    body: tuple[cp_model.IntVar, ...]
    retained: cp_model.LinearExpr


def make_rank3_model(
    weights: dict[int, int],
    upper_bound: int | None = None,
) -> CoverageModel:
    """Encode retained hyperedges via z_e = AND_{i in e}(not x_i)."""
    model = cp_model.CpModel()
    body = tuple(model.new_bool_var(f"x_{index}")
                 for index in range(len(POOL)))
    model.add(sum(body) == BODY_SIZE)
    retained_terms: list[cp_model.IntVar] = []
    retained_coefficients: list[int] = []
    for mask, width in sorted(weights.items()):
        if mask == 0:
            continue
        vertices = tuple(index for index in range(len(POOL))
                         if mask >> index & 1)
        atom = model.new_bool_var(f"z_{mask:08x}")
        for index in vertices:
            model.add(atom + body[index] <= 1)
        model.add(atom + sum(body[index] for index in vertices) >= 1)
        retained_terms.append(atom)
        retained_coefficients.append(width)
    retained = (
        cp_model.LinearExpr.weighted_sum(retained_terms, retained_coefficients)
        + weights.get(0, 0)
    )
    if upper_bound is not None:
        model.add(retained <= upper_bound)
    return CoverageModel(model, body, retained)


def exact_rank3_minimum(weights: dict[int, int]) -> tuple[int, int, int]:
    encoded = make_rank3_model(weights)
    encoded.model.minimize(encoded.retained)
    solver = cp_model.CpSolver()
    # Optimization may use independent parallel subsolvers; OPTIMAL is still
    # required.  Exhaustive solution enumeration below remains single-threaded.
    solver.parameters.num_search_workers = 8
    solver.parameters.random_seed = 314159
    solver.parameters.cp_model_presolve = True
    status = solver.solve(encoded.model)
    demand(status == cp_model.OPTIMAL, "CP-SAT did not prove a rank-three optimum")
    optimum = solver.value(encoded.retained)
    body_mask = sum(1 << index for index, variable in enumerate(encoded.body)
                    if solver.value(variable))
    direct = sum(width for mask, width in weights.items()
                 if mask & body_mask == 0)
    demand(body_mask.bit_count() == BODY_SIZE and direct == optimum,
           "CP-SAT optimum failed direct exact audit")
    return optimum, body_mask, solver.num_branches


class SubthresholdCollector(cp_model.CpSolverSolutionCallback):
    def __init__(
        self,
        body_variables: tuple[cp_model.IntVar, ...],
        rank3_weights: dict[int, int],
        full_weights: dict[int, int],
        grid: int,
    ) -> None:
        super().__init__()
        self.body_variables = body_variables
        self.rank3_items = tuple(rank3_weights.items())
        self.full_items = tuple(full_weights.items())
        self.grid = grid
        self.rows: list[tuple[int, int, int]] = []

    def on_solution_callback(self) -> None:
        body = sum(1 << index for index, variable
                   in enumerate(self.body_variables) if self.value(variable))
        retained3 = sum(width for mask, width in self.rank3_items
                        if mask & body == 0)
        retained_full = sum(width for mask, width in self.full_items
                            if mask & body == 0)
        demand(body.bit_count() == BODY_SIZE, "enumerated body has wrong size")
        demand(63 * retained3 < 4 * self.grid,
               "enumerated body is not rank-three subthreshold")
        demand(retained3 <= retained_full,
               "rank-three retained mass exceeds full retained mass")
        demand(63 * retained_full > 4 * self.grid,
               "full wall mass did not strictly clear 4/63")
        self.rows.append((body, retained3, retained_full))


def enumerate_subthreshold(
    rank3: dict[int, int], full: dict[int, int], grid: int,
) -> tuple[list[tuple[int, int, int]], int]:
    strict_upper = (4 * grid - 1) // 63
    encoded = make_rank3_model(rank3, upper_bound=strict_upper)
    # A fixed body-first search makes the exhaustive transcript reproducible
    # and was materially faster than CP-SAT's default portfolio on the largest
    # exceptional row.  This is still a solver path distinct from the primary
    # probe's weighted-marginal branch-and-bound.
    encoded.model.add_decision_strategy(
        encoded.body, cp_model.CHOOSE_FIRST, cp_model.SELECT_MIN_VALUE
    )
    collector = SubthresholdCollector(encoded.body, rank3, full, grid)
    solver = cp_model.CpSolver()
    solver.parameters.enumerate_all_solutions = True
    solver.parameters.num_search_workers = 1
    solver.parameters.random_seed = 271828
    solver.parameters.search_branching = cp_model.FIXED_SEARCH
    status = solver.solve(encoded.model, collector)
    demand(status == cp_model.OPTIMAL,
           "CP-SAT did not exhaust the subthreshold feasibility model")
    collector.rows.sort()
    demand(len({body for body, _, _ in collector.rows}) == len(collector.rows),
           "subthreshold enumeration repeated a body")
    return collector.rows, solver.num_branches


def direct_rank3(weights: dict[int, int], body: int) -> int:
    return sum(width for mask, width in weights.items() if mask & body == 0)


def small_solver_control() -> None:
    active_weights = {
        mask: 19 + 23 * mask + 5 * mask.bit_count()
        for mask in range(1 << 8) if mask.bit_count() <= 3
    }
    flat = min(
        (direct_rank3(active_weights, sum(1 << index for index in chosen)),
         sum(1 << index for index in chosen))
        for chosen in combinations(range(8), 3)
    )
    # Embed the eight active vertices in the 30-variable production model.
    # Six enormous singleton atoms force six known padding vertices into every
    # optimum, leaving exactly three choices for the exhaustively audited part.
    weights = dict(active_weights)
    forcing_weight = sum(active_weights.values()) + 1
    for index in range(8, 14):
        weights[1 << index] = forcing_weight
    optimum, body, _ = exact_rank3_minimum(weights)
    demand(optimum == flat[0], "small exhaustive/CP-SAT optimum disagreement")
    demand((body & ((1 << 8) - 1)).bit_count() == 3
           and all(body >> index & 1 for index in range(8, 14)),
           "small control padding vertices were not forced")


def digest_update(hasher: "hashlib._Hash", fields: tuple[object, ...]) -> None:
    hasher.update(("|".join(str(field) for field in fields) + "\n").encode())


def main() -> None:
    started = time.perf_counter()
    demand(len(POOL) == len(set(POOL)) == 30, "pool typing changed")
    thirds = typed_thirds()
    demand(len(thirds) == 737, "typed third-outsider universe changed")
    demand(FIXED_OUTSIDERS[0] not in POOL and FIXED_OUTSIDERS[1] not in POOL,
           "fixed outsiders entered the pool")
    # One-speed geometry control: the safe arc occupies exactly 12/14=6/7.
    unit_grid = 14
    unit_walls = (0, 1, 13, 14)
    unit_mass = sum(
        right - left for left, right in zip(unit_walls, unit_walls[1:])
        if midpoint_is_safe(1, left + right, unit_grid)
    )
    demand(unit_mass == 12, "one-speed wall geometry control failed")
    small_solver_control()

    row_digest = hashlib.sha256()
    geometry_digest = hashlib.sha256()
    candidate_digest = hashlib.sha256()
    row_records: list[tuple[Fraction, int, int, int, int, int, int]] = []
    below_payload: dict[int, tuple[int, dict[int, int], dict[int, int]]] = {}
    total_cells = 0
    total_rank3_masks = 0
    total_full_masks = 0

    for row_index, third in enumerate(thirds, start=1):
        grid, rank3, full, cells = exact_wall_weights(
            FIXED_OUTSIDERS + (third,)
        )
        optimum, witness, branches = exact_rank3_minimum(rank3)
        relation = -1 if 63 * optimum < 4 * grid else (
            0 if 63 * optimum == 4 * grid else 1
        )
        row_records.append(
            (Fraction(optimum, grid), third, optimum, grid,
             witness, branches, relation)
        )
        digest_update(row_digest, (third, grid, optimum, relation))
        digest_update(geometry_digest, ("ROW", third, grid, cells))
        for mask, width in sorted(full.items()):
            digest_update(geometry_digest, (third, mask, width))
        if relation < 0:
            below_payload[third] = (grid, rank3, full)
        total_cells += cells
        total_rank3_masks += len(rank3)
        total_full_masks += len(full)
        if row_index % 25 == 0 or row_index == len(thirds):
            print(f"progress rows={row_index}/{len(thirds)} third={third}",
                  file=sys.stderr, flush=True)

    minimum = min(row_records)
    above = sum(row[-1] > 0 for row in row_records)
    equal = sum(row[-1] == 0 for row in row_records)
    below = sum(row[-1] < 0 for row in row_records)
    demand((above, equal, below) == (729, 0, 8),
           "scan trichotomy differs from audited claim")

    print("LRC14_RANK3_WALL_HYPERGRAPH_INDEPENDENT_20260901")
    print("METHOD exact rational walls; floor(14*v*t) classification; "
          "OR-Tools CP-SAT optimization and exhaustive feasibility enumeration")
    print(f"ENV python={platform.python_version()} ortools={ortools.__version__}")
    print("SCOPE pool=30 body=9 fixed_outsiders=50,70 endpoint=769 "
          f"typed_rows={len(thirds)}")
    print("CONTROLS one_speed_mass=6/7 small_flat_vs_cpsat=PASS")
    print(f"GEOMETRY cells={total_cells} full_masks={total_full_masks} "
          f"rank3_masks={total_rank3_masks} sha256={geometry_digest.hexdigest()}")
    print(f"SCAN above={above} equal={equal} below={below} "
          f"row_sha256={row_digest.hexdigest()}")
    print(f"SCAN_MIN third={minimum[1]} ratio={minimum[0]} "
          f"mass={minimum[2]} grid={minimum[3]}")

    total_candidates = 0
    global_full: tuple[Fraction, int, int, int, int, int] | None = None
    for third in sorted(below_payload):
        grid, rank3, full = below_payload[third]
        rows, branches = enumerate_subthreshold(rank3, full, grid)
        demand(rows, f"missing subthreshold body at third={third}")
        for body, retained3, retained_full in rows:
            digest_update(candidate_digest,
                          (third, body, retained3, retained_full))
        min_rank3 = min((retained3, body) for body, retained3, _ in rows)
        min_full = min((retained_full, body, retained3)
                       for body, retained3, retained_full in rows)
        candidate = (
            Fraction(min_full[0], grid), third, min_full[1],
            min_full[0], min_full[2], grid,
        )
        if global_full is None or candidate < global_full:
            global_full = candidate
        total_candidates += len(rows)
        print(
            f"BELOW third={third} bodies={len(rows)} "
            f"min_L3={min_rank3[0]} min_L3_body={min_rank3[1]:08x} "
            f"min_full={min_full[0]} min_full_body={min_full[1]:08x} "
            f"L3_on_min_full={min_full[2]} "
            f"strict_margin={63 * min_full[0] - 4 * grid} "
            f"enumeration_branches={branches}"
        )

    demand(total_candidates == 40_511,
           "subthreshold-body total differs from audited claim")
    demand(global_full is not None, "global full-mass minimum missing")
    print(f"CANDIDATES total={total_candidates} all_full_strict=PASS "
          f"sha256={candidate_digest.hexdigest()}")
    print(f"GLOBAL_FULL_MIN third={global_full[1]} body={global_full[2]:08x} "
          f"ratio={global_full[0]} mass={global_full[3]} "
          f"L3={global_full[4]} grid={global_full[5]}")
    print("CONCLUSION FINITE_EXACT all 737 typed thirds and all nine-subset "
          "bodies have full safe mass strictly greater than 4/63")
    print("LOSS fixed pair and pool only; no arbitrary third outsider; "
          "no physical-entry or LRC(14) implication")
    print(f"RUNTIME_SECONDS {time.perf_counter() - started:.3f}")


if __name__ == "__main__":
    main()
