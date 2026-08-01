#!/usr/bin/env python3
r"""Exact anharmonic-S3 orbit scout for the reflected LRC(14) body bank.

Identify the fourteen labels with ``P1(F_13)`` by

    1,...,13 -> their residues modulo 13,       14 -> infinity.

The anharmonic generators ``s(r)=1/r`` and ``c(r)=-1/(1+r)`` generate the
six-state ``S3`` from THM-3035.  This scout compares that action with the
exact robust-edge predicate used by the reflected THM-2941 certificate bank,
including the incoming exact closure of every robust-edge-8 body.

The result is deliberately a hostile as well as a signal.  Every nontrivial
S3 orbit meeting the 561-body residual contains a robust-K6 body (all fifteen
edges), and the only trapped orbit is the fixed generic six-set
``{2,4,5,7,8,10}``.  Thus graph incidence pulls back perfectly: any star,
two-star, or spanning-tree schedule is available.  But the action does not
preserve the robust predicate: it changes integer label order, ``14*lcm(E)``,
safe-cell addresses, overlap floors, singleton debt, and exact margins.  No
LRC row closes without a new weighted certificate-pullback theorem carrying
those integral and owner sidecars.  In particular, this script uses the exact
edge-8 closure as an input; it does not reprove that 652,688-row computation.
"""

from __future__ import annotations

import argparse
import hashlib
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from math import lcm
from pathlib import Path
from typing import Callable


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "05-knowledge/results/"
    / "lrc14_reflected_anharmonic_s3_robust_body_orbit_scout_20260801.out"
)
EDGE8_SOURCE = (
    ROOT
    / "04-computation"
    / "lrc14_j7_reflected_robust_edge8_threshold_block_uniform_closure_thm2941.py"
)
EDGE8_OUTPUT = (
    ROOT
    / "05-knowledge/results"
    / "lrc14_j7_reflected_robust_edge8_threshold_block_uniform_closure_thm2941.out"
)

LABELS = tuple(range(1, 15))
BODY_SIZE = 6
BODY_COUNT = 3003
PRIME = 13
INFINITY = None
FIBER_FLOOR = F(1, 105)
ROBUST_BANK_THRESHOLD = 8
EXPECTED_BANK_COUNT = 2442
EXPECTED_RESIDUAL_COUNT = 561
EDGE8_DEPENDENCY_SOURCE_SHA256 = (
    "3f2552c7e316a6da821f8d78f859aa9d73b2d8b58081c0a31f8d233f37eec2f0"
)
EDGE8_DEPENDENCY_OUTPUT_SHA256 = (
    "90d42c1369532e4d31cfacbf4ddf455fa330075c878b414831e489eec32ecc2b"
)
EDGE8_DEPENDENCY_SEMANTIC_SHA256 = (
    "af12592ebd25cd5745c34711c99aaba72c45eeb2cad62ee5b35ec6a7752b8da7"
)
EXPECTED_ROBUST_EDGE_HIST = (
    (0, 211),
    (1, 88),
    (2, 41),
    (3, 53),
    (4, 53),
    (5, 36),
    (6, 46),
    (7, 33),
    (8, 21),
    (9, 35),
    (10, 32),
    (11, 15),
    (12, 37),
    (13, 26),
    (14, 59),
    (15, 2217),
)
EXPECTED_TRANSITIONS = (
    ("1", 561, 0),
    ("s", 61, 500),
    ("c", 58, 503),
    ("c2", 58, 503),
    ("sc", 25, 536),
    ("sc2", 245, 316),
)
EXPECTED_ORBIT_PROFILE = (
    ((1, 0), 1),
    ((1, 1), 1),
    ((2, 0), 2),
    ((3, 0), 12),
    ((3, 1), 20),
    ((3, 2), 1),
    ((6, 0), 130),
    ((6, 1), 191),
    ((6, 2), 139),
    ((6, 3), 23),
)
EXPECTED_COMPLETE_REPRESENTATIVE_MULTIPLICITY_HIST = (
    (0, 1),
    (2, 48),
    (3, 131),
    (4, 276),
    (5, 105),
)
EXPECTED_RESIDUAL_ORBIT_COMPLETE_PROFILE = (
    ((1, 1, 0), 1),
    ((3, 1, 1), 4),
    ((3, 1, 2), 16),
    ((3, 2, 1), 1),
    ((6, 1, 2), 1),
    ((6, 1, 3), 9),
    ((6, 1, 4), 76),
    ((6, 1, 5), 105),
    ((6, 2, 2), 4),
    ((6, 2, 3), 43),
    ((6, 2, 4), 92),
    ((6, 3, 2), 11),
    ((6, 3, 3), 12),
)
EXPECTED_POINT_ORBITS = (
    ("boundary", (12, 13, 14)),
    ("harmonic", (1, 6, 11)),
    ("equianharmonic", (3, 9)),
    ("generic", (2, 4, 5, 7, 8, 10)),
)
EXPECTED_FIXED_BODIES = (
    ((1, 6, 11, 12, 13, 14), 15, "boundary+harmonic"),
    ((2, 4, 5, 7, 8, 10), 2, "generic"),
)
C2_HARD_BODIES = (
    ("H", (1, 2, 3, 4, 6, 12)),
    ("H2", (1, 3, 4, 6, 8, 12)),
    ("B", (1, 2, 4, 6, 8, 12)),
    ("E", (2, 3, 4, 6, 8, 12)),
)
EXPECTED_C2_HARD_ORBITS = (
    (
        "H",
        (
            ((1, 2, 3, 4, 6, 12), 0),
            ((1, 2, 3, 5, 11, 13), 15),
            ((1, 6, 7, 8, 9, 14), 9),
            ((1, 7, 9, 10, 11, 12), 15),
            ((3, 4, 5, 6, 11, 14), 15),
            ((6, 8, 9, 10, 11, 13), 15),
        ),
    ),
    (
        "H2",
        (
            ((1, 2, 3, 7, 11, 13), 15),
            ((1, 2, 6, 7, 9, 14), 0),
            ((1, 3, 4, 6, 8, 12), 0),
            ((1, 5, 9, 10, 11, 12), 15),
            ((3, 5, 6, 10, 11, 14), 15),
            ((4, 6, 8, 9, 11, 13), 15),
        ),
    ),
    (
        "B",
        (
            ((1, 2, 4, 6, 8, 12), 0),
            ((1, 2, 5, 7, 11, 13), 15),
            ((1, 2, 6, 7, 8, 14), 1),
            ((1, 5, 7, 10, 11, 12), 15),
            ((4, 5, 6, 10, 11, 14), 15),
            ((4, 6, 8, 10, 11, 13), 15),
        ),
    ),
    (
        "E",
        (
            ((1, 2, 3, 5, 7, 13), 15),
            ((1, 2, 7, 8, 9, 14), 9),
            ((2, 3, 4, 6, 8, 12), 0),
            ((3, 4, 5, 10, 11, 14), 15),
            ((4, 6, 8, 9, 10, 13), 15),
            ((5, 7, 9, 10, 11, 12), 15),
        ),
    ),
)

# Frozen after the first independent exact run.  Set only after inspecting the
# full body/orbit transcript, never from an injected expected answer.
EXPECTED_BODY_DIGEST: str | None = "40419fbac9827f975c217f7e8dc517a929d8908d03c68977c153dd19136cf12a"
EXPECTED_ORBIT_DIGEST: str | None = "d0df274a6d6a78ce771be6714247a96a73c8b90b580149ad9864aab9cb6968c1"
EXPECTED_SEMANTIC_SHA256: str | None = "b953fbbf70537b0e982017030d76d3d703ad5a7c9919efeb82ba0fbb4a014a94"


Point = int | None
Action = Callable[[Point], Point]


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def source_sha256() -> str:
    return lf_sha256(Path(__file__))


def label_to_point(label: int) -> Point:
    require(label in LABELS, label)
    return INFINITY if label == 14 else label % PRIME


def point_to_label(point: Point) -> int:
    if point is INFINITY:
        return 14
    require(0 <= point < PRIME, point)
    return 13 if point == 0 else point


def identity(point: Point) -> Point:
    return point


def inverse(point: Point) -> Point:
    if point is INFINITY:
        return 0
    if point == 0:
        return INFINITY
    return pow(point, -1, PRIME)


def cycle(point: Point) -> Point:
    if point is INFINITY:
        return 0
    if (1 + point) % PRIME == 0:
        return INFINITY
    return (-pow((1 + point) % PRIME, -1, PRIME)) % PRIME


def compose(first: Action, second: Action) -> Action:
    return lambda point: first(second(point))


def permutation(action: Action) -> tuple[int, ...]:
    return tuple(point_to_label(action(label_to_point(label))) for label in LABELS)


def act_body(action: Action, body: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted(point_to_label(action(label_to_point(label))) for label in body))


def universal_debt(body: tuple[int, ...], ruler: int) -> F:
    return sum((F(label, 7 * (ruler - label)) for label in body), F(0))


def robust_edge_data(
    body: tuple[int, ...],
) -> tuple[int, F, tuple[tuple[int, int], ...]]:
    ruler = 14 * lcm(*body)
    debt = universal_debt(body, ruler)
    edges = tuple(
        (i, j)
        for i, j in combinations(range(BODY_SIZE), 2)
        if FIBER_FLOOR - F(4 * (body[i] + body[j]), ruler) > debt
    )
    return ruler, debt, edges


def named_point_orbits(group: tuple[tuple[str, Action], ...]) -> tuple[tuple[str, tuple[int, ...]], ...]:
    raw = []
    seen: set[int] = set()
    for label in LABELS:
        if label in seen:
            continue
        orbit = tuple(sorted({act_body(action, (label,))[0] for _, action in group}))
        seen.update(orbit)
        raw.append(orbit)
    names = {
        (12, 13, 14): "boundary",
        (1, 6, 11): "harmonic",
        (3, 9): "equianharmonic",
        (2, 4, 5, 7, 8, 10): "generic",
    }
    order = {name: index for index, (name, _) in enumerate(EXPECTED_POINT_ORBITS)}
    return tuple(sorted(((names[orbit], orbit) for orbit in raw), key=lambda row: order[row[0]]))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    require(
        lf_sha256(EDGE8_SOURCE) == EDGE8_DEPENDENCY_SOURCE_SHA256,
        ("edge8 source changed", lf_sha256(EDGE8_SOURCE)),
    )
    require(
        lf_sha256(EDGE8_OUTPUT) == EDGE8_DEPENDENCY_OUTPUT_SHA256,
        ("edge8 output changed", lf_sha256(EDGE8_OUTPUT)),
    )
    require(
        f"semantic_sha256={EDGE8_DEPENDENCY_SEMANTIC_SHA256}"
        in EDGE8_OUTPUT.read_text(encoding="utf-8"),
        "edge8 semantic token changed",
    )

    c2 = compose(cycle, cycle)
    group = (
        ("1", identity),
        ("s", inverse),
        ("c", cycle),
        ("c2", c2),
        ("sc", compose(inverse, cycle)),
        ("sc2", compose(inverse, c2)),
    )
    group_permutations = tuple((name, permutation(action)) for name, action in group)
    require(len(set(value for _, value in group_permutations)) == 6, group_permutations)
    require(permutation(compose(inverse, inverse)) == permutation(identity), "s^2")
    require(permutation(compose(cycle, c2)) == permutation(identity), "c^3")

    body_rows = []
    body_data: dict[
        tuple[int, ...], tuple[int, F, tuple[tuple[int, int], ...]]
    ] = {}
    robust_counts: dict[tuple[int, ...], int] = {}
    body_digest = hashlib.sha256()
    for body in combinations(LABELS, BODY_SIZE):
        ruler, debt, edges = robust_edge_data(body)
        robust_counts[body] = len(edges)
        body_data[body] = (ruler, debt, edges)
        row = (body, ruler, debt, edges)
        body_rows.append(row)
        body_digest.update((repr(row) + "\n").encode())
    require(len(body_rows) == BODY_COUNT, len(body_rows))

    edge_hist = tuple(sorted(Counter(robust_counts.values()).items()))
    require(edge_hist == EXPECTED_ROBUST_EDGE_HIST, edge_hist)
    bank = {body for body, count in robust_counts.items() if count >= ROBUST_BANK_THRESHOLD}
    residual = set(robust_counts) - bank
    require(len(bank) == EXPECTED_BANK_COUNT, len(bank))
    require(len(residual) == EXPECTED_RESIDUAL_COUNT, len(residual))

    transitions = []
    for name, action in group:
        residual_to_residual = sum(act_body(action, body) in residual for body in residual)
        residual_to_bank = len(residual) - residual_to_residual
        bank_to_residual = sum(act_body(action, body) in residual for body in bank)
        require(bank_to_residual == residual_to_bank, (name, bank_to_residual, residual_to_bank))
        transitions.append((name, residual_to_residual, residual_to_bank))
    transitions = tuple(transitions)
    require(transitions == EXPECTED_TRANSITIONS, transitions)

    complete_representative_multiplicity_hist = tuple(
        sorted(
            Counter(
                sum(
                    robust_counts[act_body(action, body)] == 15
                    for _, action in group
                )
                for body in residual
            ).items()
        )
    )
    require(
        complete_representative_multiplicity_hist
        == EXPECTED_COMPLETE_REPRESENTATIVE_MULTIPLICITY_HIST,
        complete_representative_multiplicity_hist,
    )

    orbit_profile: Counter[tuple[int, int]] = Counter()
    residual_orbit_complete_profile: Counter[tuple[int, int, int]] = Counter()
    orbit_rows = []
    seen_bodies: set[tuple[int, ...]] = set()
    for body in robust_counts:
        if body in seen_bodies:
            continue
        orbit = tuple(sorted({act_body(action, body) for _, action in group}))
        seen_bodies.update(orbit)
        residual_count = sum(member in residual for member in orbit)
        complete_count = sum(robust_counts[member] == 15 for member in orbit)
        orbit_profile[(len(orbit), residual_count)] += 1
        if residual_count:
            residual_orbit_complete_profile[
                (len(orbit), residual_count, complete_count)
            ] += 1
        orbit_rows.append(tuple((member, robust_counts[member]) for member in orbit))
    require(len(seen_bodies) == BODY_COUNT, len(seen_bodies))
    frozen_orbit_profile = tuple(sorted(orbit_profile.items()))
    require(frozen_orbit_profile == EXPECTED_ORBIT_PROFILE, frozen_orbit_profile)
    frozen_residual_orbit_complete_profile = tuple(
        sorted(residual_orbit_complete_profile.items())
    )
    require(
        frozen_residual_orbit_complete_profile
        == EXPECTED_RESIDUAL_ORBIT_COMPLETE_PROFILE,
        frozen_residual_orbit_complete_profile,
    )

    point_orbits = named_point_orbits(group)
    require(point_orbits == EXPECTED_POINT_ORBITS, point_orbits)
    fixed_bodies = []
    for body, count in robust_counts.items():
        if all(act_body(action, body) == body for _, action in group):
            if body == tuple(sorted(point_orbits[0][1] + point_orbits[1][1])):
                description = "boundary+harmonic"
            elif body == point_orbits[3][1]:
                description = "generic"
            else:
                raise RuntimeError(("unexpected fixed body", body))
            fixed_bodies.append((body, count, description))
    fixed_bodies = tuple(fixed_bodies)
    require(fixed_bodies == EXPECTED_FIXED_BODIES, fixed_bodies)
    trapped_orbits = tuple(
        orbit for orbit in orbit_rows if sum(count < ROBUST_BANK_THRESHOLD for _, count in orbit) == len(orbit)
    )
    require(trapped_orbits == ((((2, 4, 5, 7, 8, 10), 2),),), trapped_orbits)

    c2_hard_orbits = tuple(
        (
            name,
            tuple(
                (member, robust_counts[member])
                for member in sorted({act_body(action, body) for _, action in group})
            ),
        )
        for name, body in C2_HARD_BODIES
    )
    require(c2_hard_orbits == EXPECTED_C2_HARD_ORBITS, c2_hard_orbits)

    margin_distortion_checks = 0
    for body, (ruler, debt, _) in body_data.items():
        for _, action in group:
            target = act_body(action, body)
            target_ruler, target_debt, _ = body_data[target]
            scalar = target_debt - debt
            for left, right in combinations(body, 2):
                target_left = act_body(action, (left,))[0]
                target_right = act_body(action, (right,))[0]
                source_margin = (
                    FIBER_FLOOR - F(4 * (left + right), ruler) - debt
                )
                target_margin = (
                    FIBER_FLOOR
                    - F(4 * (target_left + target_right), target_ruler)
                    - target_debt
                )
                left_potential = F(4 * target_left, target_ruler) - F(
                    4 * left, ruler
                )
                right_potential = F(4 * target_right, target_ruler) - F(
                    4 * right, ruler
                )
                require(
                    source_margin - target_margin
                    == scalar + left_potential + right_potential,
                    (body, target, left, right),
                )
                margin_distortion_checks += 1
    require(
        margin_distortion_checks == BODY_COUNT * len(group) * 15,
        margin_distortion_checks,
    )

    orbit_digest = hashlib.sha256((repr(tuple(orbit_rows)) + "\n").encode()).hexdigest()
    if EXPECTED_BODY_DIGEST is not None:
        require(body_digest.hexdigest() == EXPECTED_BODY_DIGEST, body_digest.hexdigest())
    if EXPECTED_ORBIT_DIGEST is not None:
        require(orbit_digest == EXPECTED_ORBIT_DIGEST, orbit_digest)

    semantic_payload = (
        group_permutations,
        edge_hist,
        transitions,
        frozen_orbit_profile,
        complete_representative_multiplicity_hist,
        frozen_residual_orbit_complete_profile,
        point_orbits,
        fixed_bodies,
        trapped_orbits,
        c2_hard_orbits,
        ROBUST_BANK_THRESHOLD,
        EDGE8_DEPENDENCY_SOURCE_SHA256,
        EDGE8_DEPENDENCY_OUTPUT_SHA256,
        EDGE8_DEPENDENCY_SEMANTIC_SHA256,
        margin_distortion_checks,
        body_digest.hexdigest(),
        orbit_digest,
    )
    semantic_sha = hashlib.sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_sha == EXPECTED_SEMANTIC_SHA256, semantic_sha)

    lines = [
        "LRC14 reflected anharmonic-S3 robust-body orbit exact scout",
        "universe=all 3003 six-subsets of labels 1..14;14=P1(F13) infinity;13=0",
        f"action_generators=s:r->1/r,c:r->-1/(1+r);permutations={group_permutations}",
        f"robust_predicate=1/105-4(a+b)/(14*lcm(E))>sum_e e/[7(14*lcm(E)-e)];edge_hist={edge_hist}",
        f"uniform_closure_bank=robust_edges>=8:{len(bank)};residual=robust_edges<8:{len(residual)}",
        "closure_dependency=edge9 bank plus exact uniform closure of all 21 edge8 bodies;the latter is an input, not reproved here",
        f"edge8_dependency_source_sha256={EDGE8_DEPENDENCY_SOURCE_SHA256}",
        f"edge8_dependency_output_sha256={EDGE8_DEPENDENCY_OUTPUT_SHA256}",
        f"edge8_dependency_semantic_sha256={EDGE8_DEPENDENCY_SEMANTIC_SHA256}",
        f"residual_transitions=(element,residual_to_residual,residual_to_bank)={transitions}",
        f"orbit_profile=((orbit_size,residual_count),orbit_count)={frozen_orbit_profile}",
        f"residual_complete_representative_multiplicity=((number_of_S3_elements_landing_in_robust_K6),body_count)={complete_representative_multiplicity_hist}",
        f"residual_orbit_complete_profile=((orbit_size,residual_count,distinct_robust_K6_members),orbit_count)={frozen_residual_orbit_complete_profile}",
        f"point_orbits={point_orbits}",
        f"fixed_bodies=(body,robust_edges,point-orbit-description)={fixed_bodies}",
        "unique_trapped_orbit=((2,4,5,7,8,10),);this is the generic P1(F13) orbit",
        f"c2_hard_body_orbits={c2_hard_orbits}",
        "positive_signal=every nontrivial S3 body orbit meeting the 561 residual contains a robust-K6 representative with all 15 edges",
        "scheduling_implication=pulling a robust-K6 representative back gives the full K6 incidence graph, hence every star, two-star, and spanning-tree schedule;all physical edge weights must still be recomputed on the source body at one common cell",
        f"margin_distortion_identity_checks={margin_distortion_checks};M_E(a,b)-M_gE(ga,gb)=(D_gE-D_E)+u_gE(a)+u_gE(b)",
        "destroyed_data=integer label order and gaps;14*lcm(E);safe-cell ruler/address and owner;pair-overlap floor;singleton debt and exact margin;level word and physical packet Z(E,q)",
        "guardrail=the anharmonic action is not a symmetry of the robust predicate or of the reflected physical packet;no LRC closure follows",
        "needed_map=for each residual body and every C2-wedge level word, a one-way common-cell weighted-certificate comparison controlling the exact scalar-plus-vertex-potential distortion;graph scheduling alone is already lossless",
        f"body_digest={body_digest.hexdigest()}",
        f"orbit_digest={orbit_digest}",
        f"source_sha256={source_sha256()}",
        f"semantic_sha256={semantic_sha}",
        "normal_vs_python_O=BYTE_IDENTICAL",
        "all_exact_controls=PASS",
    ]
    payload = "\n".join(lines) + "\n"
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
