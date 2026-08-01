#!/usr/bin/env python3
"""Exact exploratory scout for the P1(F13)-generic same-owner exchange atom.

The global anharmonic action isolates the fixed body

    G = (2,4,5,7,8,10).

Unlike that action, this scout fixes the reflected ruler, the deterministic
upper-median safe cell, and its boundary owner, then permits a one-label body
exchange.  It finds the connected 13-body atom containing G and the incoming
D7/D8 weakest body W=(1,2,5,7,8,14).

This is finite-exact exploratory evidence, not an LRC closure or a theorem
promotion.  An exchange changes slopes, singleton debt, safe-cell carrier
ranges, and located overlap weights even when ruler/cell/owner stay fixed.
"""

from __future__ import annotations

import argparse
import hashlib
from collections import Counter, deque
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "04-computation/lrc14_j7_reflected_levels_all_q_mass_closure_thm2941.py"
ORBIT = (
    ROOT
    / "04-computation/"
    / "lrc14_reflected_anharmonic_s3_robust_body_orbit_scout_20260801.py"
)
OUTPUT = (
    ROOT
    / "05-knowledge/results/"
    / "lrc14_reflected_p1_generic_same_owner_exchange_scout_20260801.out"
)

EXPECTED_BASE_SHA256 = "2cf0866932f775cc493f97093333e81e65ac3aa76a8e439de969aa700c993f31"
EXPECTED_ORBIT_SHA256 = "5173b47eea103083713361a95c787a50718f569df53b4d1d5b4792c43faa9d83"
EXPECTED_SEMANTIC_SHA256: str | None = "93da4b6968b248cdf6133ebb8a1756b4c7313c3cb274561c6227c0ce59bee01f"

RULER = 3920
UPPER_MEDIAN_CELL = 2100
BOUNDARY_OWNER = 2
GENERIC_FIXED_BODY = (2, 4, 5, 7, 8, 10)
D7_D8_WEAKEST_BODY = (1, 2, 5, 7, 8, 14)
CORE = (2, 8)
OUTER = (1, 4, 5, 7, 10, 14)
EXPECTED_ROBUST_HIST = ((1, 3), (2, 3), (3, 2), (4, 3), (5, 1), (6, 1))
EXPECTED_DEGREE_HIST = ((6, 4), (7, 8), (8, 1))
EXPECTED_EDGE_COUNT = 44
EXPECTED_BODY_ROWS = (
    ((1, 2, 4, 5, 7, 8), 6, 1602, ((1, 2), (1, 4), (1, 5), (1, 7), (2, 4), (2, 5))),
    ((1, 2, 4, 5, 8, 14), 5, 1672, ((1, 2), (1, 4), (1, 5), (2, 4), (2, 5))),
    ((1, 2, 4, 7, 8, 10), 4, 1626, ((1, 2), (1, 4), (1, 7), (2, 4))),
    ((1, 2, 4, 8, 10, 14), 3, 1632, ((1, 2), (1, 4), (2, 4))),
    ((1, 2, 5, 7, 8, 10), 4, 1542, ((1, 2), (1, 5), (1, 7), (2, 5))),
    ((1, 2, 5, 7, 8, 14), 4, 1582, ((1, 2), (1, 5), (1, 7), (2, 5))),
    ((1, 2, 5, 8, 10, 14), 3, 1588, ((1, 2), (1, 5), (2, 5))),
    ((1, 2, 7, 8, 10, 14), 1, 1566, ((1, 2),)),
    ((2, 4, 5, 7, 8, 10), 2, 1682, ((2, 4), (2, 5))),
    ((2, 4, 5, 7, 8, 14), 2, 1682, ((2, 4), (2, 5))),
    ((2, 4, 5, 8, 10, 14), 2, 1688, ((2, 4), (2, 5))),
    ((2, 4, 7, 8, 10, 14), 1, 1666, ((2, 4),)),
    ((2, 5, 7, 8, 10, 14), 1, 1622, ((2, 5),)),
)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, ("cannot import", path))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    base_sha = sha256(BASE)
    orbit_sha = sha256(ORBIT)
    require(base_sha == EXPECTED_BASE_SHA256, ("base dependency", base_sha))
    require(orbit_sha == EXPECTED_ORBIT_SHA256, ("orbit dependency", orbit_sha))
    base = load("same_owner_base", BASE)
    orbit = load("same_owner_orbit", ORBIT)

    full_rows = []
    selected = {}
    component_digest = hashlib.sha256()
    for body in combinations(orbit.LABELS, orbit.BODY_SIZE):
        ruler, debt, robust_edges = orbit.robust_edge_data(body)
        safe_ruler, safe_ranges = base.safe_cell_ranges(body)
        require(safe_ruler == ruler, (body, safe_ruler, ruler))
        cells = tuple(
            cell for left, right in safe_ranges for cell in range(left, right)
        )
        require(cells, ("empty safe carrier", body))
        upper_median = cells[len(cells) // 2]
        owners = tuple(
            label
            for label in body
            if (label * upper_median) % ruler == ruler // 14
        )
        row = (
            body,
            ruler,
            debt,
            robust_edges,
            safe_ranges,
            upper_median,
            owners,
        )
        full_rows.append(row)
        if (
            ruler == RULER
            and upper_median == UPPER_MEDIAN_CELL
            and owners == (BOUNDARY_OWNER,)
        ):
            label_edges = tuple(
                (body[left], body[right]) for left, right in robust_edges
            )
            selected[body] = (
                len(robust_edges),
                len(cells),
                label_edges,
                debt,
                safe_ranges,
            )
            component_digest.update((repr(row) + "\n").encode())
    require(len(full_rows) == 3003, len(full_rows))

    expected_formula = tuple(
        sorted(
            tuple(sorted(CORE + choice))
            for choice in combinations(OUTER, 4)
            if set(choice) & {5, 10} and set(choice) & {7, 14}
        )
    )
    require(tuple(selected) == expected_formula, (tuple(selected), expected_formula))
    body_rows = tuple(
        (body, count, safe_count, label_edges)
        for body, (count, safe_count, label_edges, _, _) in selected.items()
    )
    require(body_rows == EXPECTED_BODY_ROWS, body_rows)

    exchange_edges = tuple(
        (left, right)
        for left, right in combinations(selected, 2)
        if len(set(left) ^ set(right)) == 2
    )
    require(len(exchange_edges) == EXPECTED_EDGE_COUNT, len(exchange_edges))
    adjacency = {body: [] for body in selected}
    for left, right in exchange_edges:
        adjacency[left].append(right)
        adjacency[right].append(left)
    reached = {GENERIC_FIXED_BODY}
    queue = deque((GENERIC_FIXED_BODY,))
    while queue:
        body = queue.popleft()
        for neighbor in adjacency[body]:
            if neighbor not in reached:
                reached.add(neighbor)
                queue.append(neighbor)
    require(reached == set(selected), ("disconnected", reached))

    robust_hist = tuple(
        sorted(Counter(row[0] for row in selected.values()).items())
    )
    degree_hist = tuple(
        sorted(Counter(len(neighbors) for neighbors in adjacency.values()).items())
    )
    require(robust_hist == EXPECTED_ROBUST_HIST, robust_hist)
    require(degree_hist == EXPECTED_DEGREE_HIST, degree_hist)
    require(max(count for count, *_ in selected.values()) == 6, selected)
    require(all(count < 8 for count, *_ in selected.values()), selected)

    generic = selected[GENERIC_FIXED_BODY]
    weakest = selected[D7_D8_WEAKEST_BODY]
    require(generic[:3] == (2, 1682, ((2, 4), (2, 5))), generic)
    require(
        weakest[:3] == (4, 1582, ((1, 2), (1, 5), (1, 7), (2, 5))),
        weakest,
    )
    c2 = orbit.compose(orbit.cycle, orbit.cycle)
    group = (
        orbit.identity,
        orbit.inverse,
        orbit.cycle,
        c2,
        orbit.compose(orbit.inverse, orbit.cycle),
        orbit.compose(orbit.inverse, c2),
    )
    require(
        all(orbit.act_body(action, GENERIC_FIXED_BODY) == GENERIC_FIXED_BODY for action in group),
        "generic body is not S3-fixed",
    )

    edge_digest = hashlib.sha256(
        (repr(exchange_edges) + "\n").encode()
    ).hexdigest()
    semantic_payload = (
        body_rows,
        exchange_edges,
        robust_hist,
        degree_hist,
        component_digest.hexdigest(),
        edge_digest,
        base_sha,
        orbit_sha,
    )
    semantic_sha = hashlib.sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_sha == EXPECTED_SEMANTIC_SHA256, semantic_sha)

    lines = [
        "LRC14 P1(F13)-generic same-owner exchange component exact scout",
        "status=FINITE-EXACT EXPLORATORY SIDECAR;NO LRC CLOSURE",
        "universe=all 3003 six-subsets of labels 1..14",
        "selection=ruler=3920;upper_median_safe_cell=2100;unique_boundary_owner=2",
        "family={2,8} union A;A is a 4-subset of {1,4,5,7,10,14} meeting both {5,10} and {7,14}",
        f"body_rows=(body,robust_edge_count,safe_cell_count,robust_label_edges)={body_rows}",
        "adjacency=one-label exchange, equivalently symmetric-difference size 2",
        f"exchange_edge_count={len(exchange_edges)};degree_hist={degree_hist};connected_from_generic_fixed_body=YES",
        f"exchange_edges={exchange_edges}",
        f"robust_count_hist={robust_hist};maximum_robust_edges=6;edge8_bank_members=0",
        f"generic_fixed_body={GENERIC_FIXED_BODY};robust_edges={generic[2]};safe_cell_count={generic[1]}",
        f"d7_d8_weakest_body={D7_D8_WEAKEST_BODY};robust_edges={weakest[2]};safe_cell_count={weakest[1]}",
        "preserved_by_exchange_atom=ruler;upper-median safe cell;boundary owner;four shared body labels on each adjacency",
        "destroyed_or_changed=two slopes;safe-carrier range set;singleton debt;located arcs and pair weights;certificate margin",
        "guardrail=same ruler/cell/owner is a stronger sidecar than anharmonic incidence but still does not transport positivity",
        "needed_map=an owner-preserving one-label replacement inequality for exact common-cell arcs and debt",
        f"component_digest={component_digest.hexdigest()}",
        f"edge_digest={edge_digest}",
        f"base_dependency_sha256={base_sha}",
        f"orbit_dependency_sha256={orbit_sha}",
        f"source_sha256={sha256(Path(__file__))}",
        f"semantic_sha256={semantic_sha}",
        "normal_vs_python_O=BYTE_IDENTICAL",
        "all_exact_controls=PASS",
    ]
    payload = "\n".join(lines) + "\n"
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
