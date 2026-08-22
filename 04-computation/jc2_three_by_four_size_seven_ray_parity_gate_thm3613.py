#!/usr/bin/env python3
"""Exact size-seven ray parity gate for THM-3613.

The seven size-seven additive words in the exponent-two 3 x 4 atlas are
one-dimensional positive rays.  This companion proves that, from scale four
onward, every scalar-arm placement and every singleton sign graph stabilizes.
The only remaining variation in the singleton UFD exponent is the parity of
the scale.  Scales one, two, and three are exhausted separately.

This is a necessary nonentry gate.  A placement is rejected only when every
eligible negative arm coefficient belongs to a singleton component in which
its UFD base exponent is two, and hence it cannot have a simple zero at a root
of the squarefree arm polynomial.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
import importlib.util
import json
from math import gcd
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
PARENT_THEOREM = ROOT / "01-canon/theorems/THM-3606-exponent-two-three-by-four-scalar-singleton-gate-atlas.md"
PARENT_SCRIPT = ROOT / "04-computation/jc2_three_by_four_scalar_singleton_gate_atlas_thm3606.py"
PARENT_OUTPUT = ROOT / "05-knowledge/results/jc2_three_by_four_scalar_singleton_gate_atlas_thm3606.out"
EXPECTED_PARENT_THEOREM_SHA256 = "c9b99569cafa1c2b14dbc065e1c66d24ce97d3e1bd2d388e616da2f5ad5510c4"
EXPECTED_PARENT_SCRIPT_SHA256 = "9d558f5637c5cc3573214fe00817bd6acb4c749b85a9046644e2390bd1d0ad91"
EXPECTED_PARENT_OUTPUT_SHA256 = "58b1160bd8831139fc5a6d4eb5102c876df4b47a823293080413c8f6c3f995b4"
TAIL_START = 4

EXPECTED_RAYS = {
    "W002": (1, 1, 1, 1, 2),
    "W003": (1, 2, 1, 1, 1),
    "W004": (1, 2, 1, 1, 2),
    "W005": (2, 1, 2, 1, 1),
    "W006": (1, 1, 1, 2, 1),
    "W007": (2, 1, 1, 1, 1),
    "W008": (1, 1, 2, 1, 1),
}

# (number of surviving placements before this gate, number rejected here).
EXPECTED_SMALL = {
    "W002": {1: (0, 0), 2: (3, 1), 3: (4, 1)},
    "W003": {1: (1, 1), 2: (6, 4), 3: (6, 0)},
    "W004": {1: (3, 1), 2: (12, 5), 3: (10, 1)},
    "W005": {1: (4, 3), 2: (7, 5), 3: (10, 3)},
    "W006": {1: (2, 1), 2: (7, 5), 3: (8, 0)},
    "W007": {1: (2, 2), 2: (4, 4), 3: (6, 3)},
    "W008": {1: (0, 0), 2: (2, 1), 3: (4, 2)},
}

EXPECTED_TAIL = {
    "W002": {0: (6, 4), 1: (6, 2)},
    "W003": {0: (8, 4), 1: (8, 1)},
    "W004": {0: (12, 5), 1: (12, 2)},
    "W005": {0: (12, 5), 1: (12, 3)},
    "W006": {0: (12, 5), 1: (12, 2)},
    "W007": {0: (8, 4), 1: (8, 4)},
    "W008": {0: (6, 4), 1: (6, 4)},
}


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_bytes(path: Path) -> bytes:
    return path.read_bytes().replace(b"\r\n", b"\n")


def lf_sha256(path: Path) -> str:
    return sha256(lf_bytes(path)).hexdigest()


require(lf_sha256(PARENT_THEOREM) == EXPECTED_PARENT_THEOREM_SHA256,
        "THM-3606 theorem hash")
require(lf_sha256(PARENT_SCRIPT) == EXPECTED_PARENT_SCRIPT_SHA256,
        "THM-3606 script hash")
require(lf_sha256(PARENT_OUTPUT) == EXPECTED_PARENT_OUTPUT_SHA256,
        "THM-3606 output hash")

spec = importlib.util.spec_from_file_location("thm3606_parent_for_thm3613", PARENT_SCRIPT)
require(spec is not None and spec.loader is not None, "parent import spec")
G = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = G
spec.loader.exec_module(G)
P = G.parent


@dataclass(frozen=True, order=True)
class Affine:
    slope: int
    intercept: int = 0

    def __add__(self, other: "Affine") -> "Affine":
        return Affine(self.slope + other.slope, self.intercept + other.intercept)

    def evaluate(self, scale: int) -> int:
        return self.slope * scale + self.intercept

    def canonical(self) -> tuple[int, int]:
        return self.slope, self.intercept


class DSU:
    def __init__(self, vertices: tuple[str, ...]):
        self.parent = {vertex: vertex for vertex in vertices}

    def find(self, vertex: str) -> str:
        if self.parent[vertex] != vertex:
            self.parent[vertex] = self.find(self.parent[vertex])
        return self.parent[vertex]

    def union(self, left: str, right: str) -> None:
        left_root, right_root = self.find(left), self.find(right)
        if left_root != right_root:
            self.parent[max(left_root, right_root)] = min(left_root, right_root)


def rational_rank(rows: tuple[tuple[int, ...], ...]) -> int:
    matrix = [[Fraction(value) for value in row] for row in rows]
    pivot_row = 0
    for column in range(5):
        source = next(
            (index for index in range(pivot_row, len(matrix)) if matrix[index][column]),
            None,
        )
        if source is None:
            continue
        matrix[pivot_row], matrix[source] = matrix[source], matrix[pivot_row]
        pivot = matrix[pivot_row][column]
        matrix[pivot_row] = [value / pivot for value in matrix[pivot_row]]
        for index in range(len(matrix)):
            if index == pivot_row:
                continue
            multiplier = matrix[index][column]
            if multiplier:
                matrix[index] = [
                    value - multiplier * pivot_value
                    for value, pivot_value in zip(matrix[index], matrix[pivot_row])
                ]
        pivot_row += 1
    return pivot_row


def word_rays() -> dict[str, tuple[int, ...]]:
    flats = P.enumerate_flat_lattice()
    cones = P.enumerate_positive_connected_cones(flats)
    word_sample: dict[tuple[tuple[tuple[int, int], ...], ...], tuple[int, ...]] = {}
    for cone in cones:
        for sample in cone.samples:
            word_sample[P.fibre_word(P.sum_fibres(sample))] = sample
    ordered = tuple(sorted(word_sample, key=lambda word: (len(P.sum_fibres(word_sample[word])), word)))
    ray_map = {
        f"W{index:03d}": word_sample[word]
        for index, word in enumerate(ordered, 1)
        if len(P.sum_fibres(word_sample[word])) == 7
    }
    require(ray_map == EXPECTED_RAYS, ("seven primitive rays", ray_map))
    for word_id, ray in ray_map.items():
        equalities, _inequalities = G.word_constraints(ray)
        rows = tuple(row for row, rhs in equalities if rhs == 0)
        require(rational_rank(rows) == 4, (word_id, "one-dimensional equality cone"))
        require(all(sum(a * x for a, x in zip(row, ray)) == 0 for row in rows),
                (word_id, "ray in equality kernel"))
        require(gcd(*ray) == 1 and 1 in ray, (word_id, "primitive integral generator"))
        require(P.fibre_word(P.sum_fibres(tuple(2 * value for value in ray)))
                == P.fibre_word(P.sum_fibres(ray)), (word_id, "positive scaling preserves word"))
    return ray_map


def components_and_exponents(
    weights: dict[str, int],
    singleton_rows: tuple[tuple[int, int, int, int, str], ...],
) -> tuple[dict[str, tuple[str, ...]], dict[str, int]]:
    vertices = tuple(sorted(weights))
    dsu = DSU(vertices)
    active: set[str] = set()
    for i, j, _p_weight, _q_weight, verdict in singleton_rows:
        if verdict not in {"++", "--"}:
            continue
        left, right = f"P{i}", f"Q{j}"
        dsu.union(left, right)
        active.update((left, right))
    by_root: dict[str, list[str]] = {}
    for vertex in sorted(active):
        by_root.setdefault(dsu.find(vertex), []).append(vertex)
    component_of: dict[str, tuple[str, ...]] = {}
    exponents: dict[str, int] = {}
    for part_list in by_root.values():
        part = tuple(sorted(part_list))
        divisor = 0
        for vertex in part:
            divisor = gcd(divisor, abs(weights[vertex]))
        require(divisor > 0, ("nonzero singleton component", part))
        for vertex in part:
            component_of[vertex] = part
            exponents[vertex] = abs(weights[vertex]) // divisor
    return component_of, exponents


def direct_record(placement: object) -> dict[str, object]:
    p_weights = placement.p_weights
    q_weights = placement.q_weights
    weights = {
        **{f"P{i}": value for i, value in enumerate(p_weights)},
        **{f"Q{j}": value for j, value in enumerate(q_weights)},
    }
    component_of, exponents = components_and_exponents(weights, placement.singleton_rows)
    arm_data = []
    for i, j, orientation in placement.eligible:
        negative = f"P{i}" if orientation == "-+" else f"Q{j}"
        arm_data.append((i, j, orientation, exponents.get(negative), component_of.get(negative)))
    rejected = all(exponent not in {None, 1} for _i, _j, _o, exponent, _part in arm_data)
    return {
        "key": (placement.scalar_index, p_weights[0], q_weights[0]),
        "eligible": tuple(arm_data),
        "rejected": rejected,
        "weights": tuple(p_weights) + tuple(q_weights),
    }


def direct_scale(ray: tuple[int, ...], scale: int) -> tuple[dict[str, object], ...]:
    gaps = tuple(scale * value for value in ray)
    fibres = P.sum_fibres(gaps)
    records = []
    for placement in G.placements_for_sample(gaps):
        if P.rectangle_exposed(gaps, fibres[placement.scalar_index]):
            continue
        records.append(direct_record(placement))
    return tuple(sorted(records, key=lambda record: record["key"]))


def tail_sign(value: Affine) -> int:
    """Return the sign for every integer scale >=TAIL_START, proving stability."""
    if value.slope == 0:
        return (value.intercept > 0) - (value.intercept < 0)
    at_start = value.evaluate(TAIL_START)
    if value.slope > 0:
        require(at_start > 0, ("positive affine tail not stabilized", value))
        return 1
    require(at_start < 0, ("negative affine tail not stabilized", value))
    return -1


def affine_solution(value: Affine, target: int) -> str | int | None:
    """Solve value(n)=target over integers; 'all' denotes an identity."""
    if value.slope == 0:
        return "all" if value.intercept == target else None
    numerator = target - value.intercept
    if numerator % value.slope:
        return None
    return numerator // value.slope


def simultaneous_solution(
    left: Affine, left_target: int, right: Affine, right_target: int
) -> str | int | None:
    left_solution = affine_solution(left, left_target)
    right_solution = affine_solution(right, right_target)
    if left_solution == "all":
        return right_solution
    if right_solution == "all":
        return left_solution
    if left_solution is not None and left_solution == right_solution:
        return left_solution
    return None


def tail_records(ray: tuple[int, ...]) -> tuple[dict[str, object], ...]:
    a_support, b_support = P.supports(ray)
    fibres = P.sum_fibres(ray)
    candidates: dict[tuple[int, Affine, Affine], tuple[Affine, Affine]] = {}
    for scalar_index, scalar_fibre in enumerate(fibres):
        if len(scalar_fibre) < 2:
            continue
        for i, j in scalar_fibre:
            candidates[(scalar_index, Affine(-a_support[i], -2), Affine(-b_support[j], 1))] = (
                Affine(-a_support[i], -2), Affine(-b_support[j], 1)
            )
            candidates[(scalar_index, Affine(-a_support[i], 1), Affine(-b_support[j], -2))] = (
                Affine(-a_support[i], 1), Affine(-b_support[j], -2)
            )

    records: list[dict[str, object]] = []
    for (scalar_index, p0, q0), _ in sorted(candidates.items()):
        p_weights = tuple(p0 + Affine(value) for value in a_support)
        q_weights = tuple(q0 + Affine(value) for value in b_support)
        scalar_fibre = fibres[scalar_index]
        require(all((p_weights[i] + q_weights[j]).canonical() == (0, -1)
                    for i, j in scalar_fibre), "symbolic scalar fibre")
        eligible = []
        for i, j in scalar_fibre:
            pair = p_weights[i].canonical(), q_weights[j].canonical()
            if pair == ((0, -2), (0, 1)):
                eligible.append((i, j, "-+"))
            elif pair == ((0, 1), (0, -2)):
                eligible.append((i, j, "+-"))
            for orientation, p_target, q_target in (("-+", -2, 1), ("+-", 1, -2)):
                solution = simultaneous_solution(
                    p_weights[i], p_target, q_weights[j], q_target
                )
                if solution == "all":
                    require((i, j, orientation) in eligible,
                            ("identical eligibility recorded", i, j, orientation))
                elif solution is not None:
                    require(solution <= 3,
                            ("no exceptional eligible address in tail", i, j, orientation, solution))
        require(eligible, ("tail candidate has anchor", scalar_index, p0, q0))

        dsu = DSU(tuple([f"P{i}" for i in range(3)] + [f"Q{j}" for j in range(4)]))
        active: set[str] = set()
        singleton_ledger = []
        survives = True
        for index, fibre in enumerate(fibres):
            if index == scalar_index or len(fibre) != 1:
                continue
            i, j = fibre[0]
            left_sign, right_sign = tail_sign(p_weights[i]), tail_sign(q_weights[j])
            if left_sign == right_sign and left_sign in {-1, 1}:
                verdict = "--" if left_sign < 0 else "++"
                left, right = f"P{i}", f"Q{j}"
                dsu.union(left, right)
                active.update((left, right))
            elif left_sign == right_sign == 0:
                verdict = "00"
            else:
                survives = False
                verdict = "FAIL"
            singleton_ledger.append((i, j, verdict))
        if not survives or P.rectangle_exposed(ray, scalar_fibre):
            continue

        by_root: dict[str, list[str]] = {}
        for vertex in sorted(active):
            by_root.setdefault(dsu.find(vertex), []).append(vertex)
        component_of = {
            vertex: tuple(sorted(part))
            for part in by_root.values()
            for vertex in part
        }
        require(len(eligible) == 1, ("tail eligible address is unique", eligible))
        i, j, orientation = eligible[0]
        negative = f"P{i}" if orientation == "-+" else f"Q{j}"
        exponents = {}
        for parity in (0, 1):
            part = component_of.get(negative)
            if part is None:
                exponents[parity] = None
                continue
            component_weights = tuple(
                p_weights[int(vertex[1:])] if vertex.startswith("P")
                else q_weights[int(vertex[1:])]
                for vertex in part
            )
            # The component contains the negative arm weight -2, so its gcd is
            # either one or two.  It is two exactly when every weight is even.
            all_even = all((value.slope * parity + value.intercept) % 2 == 0
                           for value in component_weights)
            exponents[parity] = 1 if all_even else 2
        records.append({
            "affine_key": (scalar_index, p0.canonical(), q0.canonical()),
            "eligible": tuple(eligible),
            "component": component_of.get(negative),
            "arm_exponents": (exponents[0], exponents[1]),
            "singleton_ledger": tuple(singleton_ledger),
            "weights": tuple(value.canonical() for value in p_weights + q_weights),
        })
    return tuple(sorted(records, key=lambda record: record["affine_key"]))


def evaluated_tail_record(record: dict[str, object], scale: int) -> dict[str, object]:
    scalar_index, p0, q0 = record["affine_key"]
    weights = tuple(slope * scale + intercept for slope, intercept in record["weights"])
    i, j, orientation = record["eligible"][0]
    exponent = record["arm_exponents"][scale % 2]
    return {
        "key": (scalar_index, p0[0] * scale + p0[1], q0[0] * scale + q0[1]),
        "eligible": ((i, j, orientation, exponent, record["component"]),),
        "rejected": exponent == 2,
        "weights": weights,
    }


def count(records: tuple[dict[str, object], ...]) -> tuple[int, int]:
    return len(records), sum(bool(record["rejected"]) for record in records)


def compact_counts(table: dict[str, dict[int, tuple[int, int]]]) -> str:
    return ";".join(
        word_id + ":" + ",".join(
            f"n{scale}={values[0]}/{values[1]}"
            for scale, values in sorted(scales.items())
        )
        for word_id, scales in sorted(table.items())
    )


def main() -> None:
    rays = word_rays()
    small: dict[str, dict[int, tuple[int, int]]] = {}
    tails: dict[str, dict[int, tuple[int, int]]] = {}
    tail_ledgers: dict[str, tuple[dict[str, object], ...]] = {}
    dead_keys: dict[str, dict[int, tuple[object, ...]]] = {}

    for word_id, ray in rays.items():
        small[word_id] = {
            scale: count(direct_scale(ray, scale))
            for scale in range(1, TAIL_START)
        }
        symbolic = tail_records(ray)
        tail_ledgers[word_id] = symbolic
        # Independent bridge: symbolic tail evaluation agrees record-for-record
        # with the original integer implementation in both parity classes.
        for scale in (TAIL_START, TAIL_START + 1):
            evaluated = tuple(sorted(
                (evaluated_tail_record(record, scale) for record in symbolic),
                key=lambda record: record["key"],
            ))
            require(evaluated == direct_scale(ray, scale),
                    (word_id, "symbolic/direct bridge", scale))
        tails[word_id] = {}
        dead_keys[word_id] = {}
        for parity, scale in ((0, TAIL_START), (1, TAIL_START + 1)):
            evaluated = tuple(evaluated_tail_record(record, scale) for record in symbolic)
            tails[word_id][parity] = count(evaluated)
            dead_keys[word_id][parity] = tuple(
                record["affine_key"] for record in symbolic
                if record["arm_exponents"][parity] == 2
            )

    require(small == EXPECTED_SMALL, ("exceptional-scale table", small))
    require(tails == EXPECTED_TAIL, ("tail parity table", tails))
    require(sum(values[parity][0] for values in tails.values() for parity in (0, 1)) == 128,
            "two parity classes each have 64 placements")
    require(sum(values[parity][1] for values in tails.values() for parity in (0, 1)) == 49,
            "tail parity rejections 31+18")

    semantic_object = {
        "parents": (
            EXPECTED_PARENT_THEOREM_SHA256,
            EXPECTED_PARENT_SCRIPT_SHA256,
            EXPECTED_PARENT_OUTPUT_SHA256,
        ),
        "tail_start": TAIL_START,
        "eligibility_exception_bound": 3,
        "ray_equality_ranks": {word_id: 4 for word_id in rays},
        "rays": rays,
        "small": small,
        "tail": tails,
        "ledgers": tail_ledgers,
        "dead_keys": dead_keys,
    }
    semantic_sha = sha256(json.dumps(
        semantic_object, sort_keys=True, separators=(",", ":")
    ).encode("ascii")).hexdigest()

    even_total = tuple(map(sum, zip(*(tails[word_id][0] for word_id in sorted(tails)))))
    odd_total = tuple(map(sum, zip(*(tails[word_id][1] for word_id in sorted(tails)))))
    print("THM-3613 THREE-BY-FOUR SIZE-SEVEN RAY PARITY GATE")
    print(f"thm3606_theorem_sha256_lf={EXPECTED_PARENT_THEOREM_SHA256}")
    print(f"thm3606_script_sha256_lf={EXPECTED_PARENT_SCRIPT_SHA256}")
    print(f"thm3606_output_sha256_lf={EXPECTED_PARENT_OUTPUT_SHA256}")
    print("rays=" + ";".join(f"{word_id}:{','.join(map(str, ray))}:rank4:primitive" for word_id, ray in rays.items()))
    print("exceptional_counts_total/rejected=" + compact_counts(small))
    print("tail_even_counts_total/rejected=" + ";".join(
        f"{word_id}:{tails[word_id][0][0]}/{tails[word_id][0][1]}" for word_id in sorted(tails)
    ))
    print("tail_odd_counts_total/rejected=" + ";".join(
        f"{word_id}:{tails[word_id][1][0]}/{tails[word_id][1][1]}" for word_id in sorted(tails)
    ))
    print(f"tail_totals=even:{even_total[0]}/{even_total[1]};odd:{odd_total[0]}/{odd_total[1]}")
    print("mechanism=tail_sign_graph_constant_for_n>=4;no_new_eligible_after3;unique_arm;R=2/gcd(2,component_weights);parity_only")
    print(f"semantic_sha256={semantic_sha}")
    print(f"script_sha256_lf={lf_sha256(Path(__file__).resolve())}")
    print("boundary=necessary_scalar-arm/singleton gate only;surviving placements are not Darboux pairs")
    print("status=PROVED+EXACT+OPTIMIZATION-SAFE;size-seven ray residue classification")
    print("PASS")


if __name__ == "__main__":
    main()
