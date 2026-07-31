#!/usr/bin/env python3
"""Exact ordered Hunter-pivot closure on the K>=20 H2<=40 graph slice.

The symmetric theta-heavy graphs are imported from the exact H2<=40 graph
ledger.  A deterministic elimination order uses the ordered high-core pivot
lemma:

* if the pivot-star invoice plus the parent singleton cap is below h, every
  K4 assigned to that pivot closes at once;
* otherwise enumerate exactly the triangles in its forward neighborhood,
  use the maximum-spanning-tree Hunter invoice on each resulting K4, and
  send only the remaining hostile invoices to literal K4 subtraction and
  the longest-component singleton seal.

The pivot direction is an acyclic proof-enumeration gauge.  The underlying
heavy graph and overlap weights remain symmetric.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import importlib.util
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PILOT_PATH = ROOT / ".scratch/lrc_rank3_h2_clique_descent_20260729/pilot.py"
PILOT_SHA256 = (
    "9ec1b8d8c697f54fdbf2836638fa6c7dc5284f4ea2644849d98d71398c1f3520"
)
CORE_PATH = ROOT / ".scratch/lrc_rank3_h2_clique_descent_20260729/core.discovery.out"
CORE_SHA256 = (
    "5569d8d34b59e5eed2bfa82148648dcb5e63515146ff55b95234a002d0c87ba2"
)
GRAPH_PATH = (
    ROOT / ".scratch/lrc_rank3_h2_clique_descent_20260729/graph_H40.discovery.out"
)
GRAPH_SHA256 = (
    "5683dc23f70d78a6c351008fd262a6d3214905650733446feac12155ccb58ad7"
)
EXPECTED_COUNTS: tuple[int, ...] | None = (
    62,
    62,
    0,
    1662,
    12090,
    0,
    93754,
    1454,
    81429,
    208,
    12325,
    12286,
    39,
    0,
    39,
    0,
    0,
    295,
    67,
    29,
)
EXPECTED_LEDGER_DIGEST: str | None = (
    "2088b2016179d14e1b5e76d6c9eddfe88e8352c943f49970bfb201b46583cd47"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load_pilot():
    require(hashlib.sha256(PILOT_PATH.read_bytes()).hexdigest() == PILOT_SHA256, "pilot changed")
    require(hashlib.sha256(CORE_PATH.read_bytes()).hexdigest() == CORE_SHA256, "core changed")
    require(hashlib.sha256(GRAPH_PATH.read_bytes()).hexdigest() == GRAPH_SHA256, "graph changed")
    spec = importlib.util.spec_from_file_location("rank3_h2_ordered_pivot_pilot", PILOT_PATH)
    require(spec is not None and spec.loader is not None, "cannot load pilot")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


P = load_pilot()


def ftext(value: F | None) -> str:
    return "-" if value is None else f"{value.numerator}/{value.denominator}"


def parse_core_rows() -> dict[tuple[tuple[int, ...], int], dict[int, F]]:
    out = {}
    for line in CORE_PATH.read_text().splitlines():
        if not line.startswith("H2;"):
            continue
        fields = dict(part.split("=", 1) for part in line.split(";")[1:])
        key = (ast.literal_eval(fields["E"]), int(fields["rank"]))
        out[key] = {
            int(speed): F(value)
            for speed, value in ast.literal_eval(fields["Hrows"])
        }
    require(len(out) == 143, "sealed-core universe changed")
    return out


def parse_graph_rows() -> tuple[dict[str, object], ...]:
    cores = parse_core_rows()
    scalar = {
        (row["body"], row["rank"]): row
        for row in P.TARGET_ROWS
    }
    out = []
    for line in GRAPH_PATH.read_text().splitlines():
        if not line.startswith("CORE;"):
            continue
        fields = dict(part.split("=", 1) for part in line.split(";")[1:])
        key = (ast.literal_eval(fields["E"]), int(fields["rank"]))
        row = scalar[key]
        coverages = cores[key]
        weighted_edges = {
            (x, y): F(weight)
            for x, y, weight in ast.literal_eval(fields["weighted_edges"])
        }
        adjacency = {speed: set() for speed in coverages}
        for x, y in weighted_edges:
            adjacency[x].add(y)
            adjacency[y].add(x)
        require(len(coverages) == int(fields["H"]) <= 40, "graph slice changed")
        require(len(weighted_edges) == int(fields["edges"]), "edge count changed")
        out.append(
            {
                **row,
                "coverages": coverages,
                "weighted_edges": weighted_edges,
                "adjacency": adjacency,
                "initial_k4": int(fields["K4"]),
                "edge_ties": int(fields["ties"]),
            }
        )
    require(len(out) == 62, "H2<=40 graph universe changed")
    return tuple(out)


def forward_triangles(
    pivot: int,
    remaining: set[int],
    adjacency: dict[int, set[int]],
) -> tuple[tuple[int, int, int], ...]:
    neighbors = sorted(adjacency[pivot].intersection(remaining))
    neighbor_set = set(neighbors)
    out = []
    for first_index, first in enumerate(neighbors[:-2]):
        for second in neighbors[first_index + 1 : -1]:
            if second not in adjacency[first]:
                continue
            for third in sorted(
                adjacency[first]
                .intersection(adjacency[second])
                .intersection(neighbor_set)
            ):
                if third > second:
                    out.append((first, second, third))
    return tuple(out)


def maximum_tree_credit(
    clique: tuple[int, int, int, int],
    coverages: dict[int, F],
    unions: dict[tuple[int, int], F],
) -> tuple[F, tuple[tuple[tuple[int, int], ...], ...]]:
    edges = tuple(combinations(clique, 2))
    credits = {
        edge: coverages[edge[0]] + coverages[edge[1]] - unions[tuple(sorted(edge))]
        for edge in edges
    }
    spanning = []
    for tree in combinations(edges, 3):
        adjacency = {vertex: set() for vertex in clique}
        for x, y in tree:
            adjacency[x].add(y)
            adjacency[y].add(x)
        seen = set()
        stack = [clique[0]]
        while stack:
            vertex = stack.pop()
            if vertex in seen:
                continue
            seen.add(vertex)
            stack.extend(adjacency[vertex] - seen)
        if len(seen) == 4:
            spanning.append((sum((credits[edge] for edge in tree), F(0)), tree))
    require(len(spanning) == 16, "K4 spanning-tree count changed")
    maximum = max(value for value, _ in spanning)
    maximizers = tuple(tree for value, tree in spanning if value == maximum)

    # Independent Kruskal control, with arbitrary deterministic tie order.
    parent = {vertex: vertex for vertex in clique}

    def find(vertex: int) -> int:
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    kruskal = F(0)
    chosen = 0
    for weight, edge in sorted(
        ((credits[edge], edge) for edge in edges),
        key=lambda item: (-item[0], item[1]),
    ):
        x, y = edge
        root_x = find(x)
        root_y = find(y)
        if root_x == root_y:
            continue
        parent[root_x] = root_y
        kruskal += weight
        chosen += 1
        if chosen == 3:
            break
    require(chosen == 3 and kruskal == maximum, "Kruskal/brute tree credit differs")
    return maximum, maximizers


def reconstruct_carrier(row: dict[str, object]) -> list[tuple[F, F]]:
    root_good, _, _ = P.S.G.T.CORE.good_norm(row["body"])
    carrier = P.S.R.subtract_local(root_good, row["apex"])
    require(P.mass(carrier) == row["h"], "ordered-pivot carrier changed")
    return carrier


def profile_branch(row: dict[str, object]) -> dict[str, object]:
    coverages = row["coverages"]
    unions = row["weighted_edges"]
    adjacency = row["adjacency"]
    h = row["h"]
    q1 = row["top5"][0][0]
    remaining = set(coverages)
    order = []
    steps = []
    star_steps = 0
    star_assigned = 0
    paid_steps = 0
    paid_k4 = 0
    tree_closed = 0
    tree_hard = 0
    literal_closed = 0
    empty_residuals = []
    witnesses = []
    singleton_scans = 0
    singleton_controls = 0
    tree_equalities = 0
    positive_tree_extremum = None
    hostile_tree_extremum = None
    literal_rows = []
    paid_digest = hashlib.sha256(
        (
            f"LRC14/j6/H40/ordered-Hunter/paid/{row['body']}/"
            f"{row['rank']}/v1\n"
        ).encode()
    )
    carrier = None

    while remaining:
        options = []
        for pivot in remaining:
            neighbors = adjacency[pivot].intersection(remaining - {pivot})
            child_coverages = sorted(
                (
                    unions[tuple(sorted((pivot, neighbor)))] - coverages[pivot]
                    for neighbor in neighbors
                ),
                reverse=True,
            )
            burden = (
                coverages[pivot]
                + q1
                + sum(child_coverages[:3], F(0))
            )
            margin = None if len(neighbors) < 3 else h - burden
            options.append((pivot, len(neighbors), margin))

        certified = [
            option
            for option in options
            if option[2] is None or option[2] > 0
        ]
        if certified:
            # Prefer the smallest later-link degree, then the strongest strict
            # margin, with numeric speed only as the final deterministic tie.
            pivot, degree, margin = min(
                certified,
                key=lambda option: (
                    option[1],
                    F(0) if option[2] is None else -option[2],
                    option[0],
                ),
            )
            triangles = forward_triangles(pivot, remaining, adjacency)
            star_steps += 1
            star_assigned += len(triangles)
            steps.append(
                (
                    pivot,
                    "STAR" if margin is not None else "VACUOUS",
                    degree,
                    ftext(margin),
                    len(triangles),
                    0,
                )
            )
        else:
            triangle_options = [
                (
                    len(forward_triangles(pivot, remaining, adjacency)),
                    degree,
                    pivot,
                )
                for pivot, degree, _ in options
            ]
            _, degree, pivot = min(triangle_options)
            triangles = forward_triangles(pivot, remaining, adjacency)
            paid_steps += 1
            paid_k4 += len(triangles)
            local_hard = 0
            for triple in triangles:
                clique = tuple(sorted((pivot, *triple)))
                singleton_sum = sum((coverages[speed] for speed in clique), F(0))
                tree_credit, maximizers = maximum_tree_credit(
                    clique,
                    coverages,
                    unions,
                )
                psi = singleton_sum - tree_credit
                tree_margin = h - q1 - psi
                paid_digest.update(
                    (
                        f"pivot={pivot};K4={clique};"
                        f"credit={ftext(tree_credit)};psi={ftext(psi)};"
                        f"margin={ftext(tree_margin)};"
                        f"maximizers={maximizers}\n"
                    ).encode()
                )
                if tree_margin > 0:
                    tree_closed += 1
                    candidate = (tree_margin, clique, pivot)
                    if (
                        positive_tree_extremum is None
                        or candidate < positive_tree_extremum
                    ):
                        positive_tree_extremum = candidate
                    continue
                tree_hard += 1
                local_hard += 1
                tree_equalities += tree_margin == 0
                candidate = (tree_margin, clique, pivot)
                if (
                    hostile_tree_extremum is None
                    or candidate > hostile_tree_extremum
                ):
                    hostile_tree_extremum = candidate
                if carrier is None:
                    carrier = reconstruct_carrier(row)
                residual = P.S.R.subtract_local_multi(carrier, clique)
                sequential = carrier
                for speed in clique:
                    sequential = P.S.R.subtract_local(sequential, speed)
                require(
                    residual == sequential,
                    f"ordered-pivot K4 path mismatch: {row['body']}, "
                    f"{row['rank']}, {clique}",
                )
                if not residual:
                    empty_residuals.append((clique, tree_margin, maximizers))
                    literal_rows.append(
                        (
                            clique,
                            ftext(tree_margin),
                            maximizers,
                            "EMPTY",
                        )
                    )
                    continue
                singleton = P.singleton_residual_profile(row, clique, residual)
                singleton_scans += singleton["scanned"]
                singleton_controls += singleton["controls"]
                if singleton["covering"]:
                    witnesses.extend(
                        (clique, speed, tree_margin, maximizers)
                        for speed in singleton["covering"]
                    )
                else:
                    literal_closed += 1
                literal_rows.append(
                    (
                        clique,
                        ftext(tree_margin),
                        maximizers,
                        ftext(singleton["mass"]),
                        singleton["components"],
                        ftext(singleton["longest"]),
                        singleton["tail_first"],
                        singleton["scanned"],
                        singleton["controls"],
                        singleton["covering"],
                    )
                )
            steps.append(
                (
                    pivot,
                    "PAID",
                    degree,
                    "-",
                    len(triangles),
                    local_hard,
                )
            )
        order.append(pivot)
        remaining.remove(pivot)

    assigned = star_assigned + paid_k4
    require(assigned == row["initial_k4"], "ordered K4 partition changed")
    closed = not empty_residuals and not witnesses
    return {
        **row,
        "order": tuple(order),
        "steps": tuple(steps),
        "star_steps": star_steps,
        "star_assigned": star_assigned,
        "paid_steps": paid_steps,
        "paid_k4": paid_k4,
        "tree_closed": tree_closed,
        "tree_hard": tree_hard,
        "literal_closed": literal_closed,
        "empty_residuals": tuple(empty_residuals),
        "witnesses": tuple(witnesses),
        "singleton_scans": singleton_scans,
        "singleton_controls": singleton_controls,
        "tree_equalities": tree_equalities,
        "positive_tree_extremum": positive_tree_extremum,
        "hostile_tree_extremum": hostile_tree_extremum,
        "literal_rows": tuple(literal_rows),
        "paid_digest": paid_digest.hexdigest(),
        "closed": closed,
    }


def ledger_line(row: dict[str, object]) -> str:
    return (
        f"E={row['body']};K={row['K']};rank={row['rank']};a={row['apex']};"
        f"P={row['prefix']};h={ftext(row['h'])};H={len(row['coverages'])};"
        f"edges={len(row['weighted_edges'])};initialK4={row['initial_k4']};"
        f"star_steps={row['star_steps']};starK4={row['star_assigned']};"
        f"paid_steps={row['paid_steps']};paidK4={row['paid_k4']};"
        f"tree_closed={row['tree_closed']};tree_hard={row['tree_hard']};"
        f"tree_equalities={row['tree_equalities']};"
        f"positive_tree_extremum={row['positive_tree_extremum']};"
        f"hostile_tree_extremum={row['hostile_tree_extremum']};"
        f"literal_closed={row['literal_closed']};"
        f"scan={row['singleton_scans']};controls={row['singleton_controls']};"
        f"empty={row['empty_residuals']};witnesses={row['witnesses']};"
        f"closed={int(row['closed'])};order={row['order']};steps={row['steps']};"
        f"LITERAL_ROWS={row['literal_rows']};"
        f"paid_digest={row['paid_digest']}\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--emit-ledger", action="store_true")
    args = parser.parse_args()
    rows = tuple(profile_branch(row) for row in parse_graph_rows())
    digest = hashlib.sha256(b"LRC14/j6/highK-r3-H2-H40-ordered-Hunter/v1\n")
    for row in rows:
        line = ledger_line(row)
        digest.update(line.encode())
        if args.emit_ledger:
            print("ROW;" + line.rstrip())
    counts = (
        len(rows),
        sum(row["closed"] for row in rows),
        sum(not row["closed"] for row in rows),
        sum(len(row["coverages"]) for row in rows),
        sum(len(row["weighted_edges"]) for row in rows),
        sum(row["edge_ties"] for row in rows),
        sum(row["initial_k4"] for row in rows),
        sum(row["star_steps"] for row in rows),
        sum(row["star_assigned"] for row in rows),
        sum(row["paid_steps"] for row in rows),
        sum(row["paid_k4"] for row in rows),
        sum(row["tree_closed"] for row in rows),
        sum(row["tree_hard"] for row in rows),
        sum(row["tree_equalities"] for row in rows),
        sum(row["literal_closed"] for row in rows),
        sum(len(row["empty_residuals"]) for row in rows),
        sum(len(row["witnesses"]) for row in rows),
        sum(row["singleton_scans"] for row in rows),
        sum(row["singleton_controls"] for row in rows),
        sum(row["paid_steps"] == 0 for row in rows),
    )
    ledger_digest = digest.hexdigest()
    if EXPECTED_COUNTS is not None:
        require(counts == EXPECTED_COUNTS, "ordered-Hunter counts changed")
    if EXPECTED_LEDGER_DIGEST is not None:
        require(
            ledger_digest == EXPECTED_LEDGER_DIGEST,
            "ordered-Hunter ledger changed",
        )
    print("LRC14 j6 K>=20 strongest-r3 H2<=40 ordered Hunter closure")
    print(
        f"branches={len(rows)};closed={sum(row['closed'] for row in rows)};"
        f"survivors={sum(not row['closed'] for row in rows)};"
        f"vertices={sum(len(row['coverages']) for row in rows)};"
        f"edges={sum(len(row['weighted_edges']) for row in rows)};"
        f"initial_K4={sum(row['initial_k4'] for row in rows)}"
    )
    print(
        f"ordered_partition=star_steps:{sum(row['star_steps'] for row in rows)},"
        f"star_K4:{sum(row['star_assigned'] for row in rows)},"
        f"paid_steps:{sum(row['paid_steps'] for row in rows)},"
        f"paid_K4:{sum(row['paid_k4'] for row in rows)}"
    )
    print(
        f"paid_closure=tree:{sum(row['tree_closed'] for row in rows)},"
        f"tree_hard:{sum(row['tree_hard'] for row in rows)},"
        f"literal:{sum(row['literal_closed'] for row in rows)},"
        f"empty:{sum(len(row['empty_residuals']) for row in rows)},"
        f"witnesses:{sum(len(row['witnesses']) for row in rows)}"
    )
    print(
        f"singleton=scans:{sum(row['singleton_scans'] for row in rows)},"
        f"controls:{sum(row['singleton_controls'] for row in rows)};"
        f"star_only_branches:{sum(row['paid_steps']==0 for row in rows)}"
    )
    print(f"counts={counts}")
    print(f"canonical_ledger_sha256={ledger_digest}")
    print(
        "mode="
        + (
            "DISCOVERY"
            if EXPECTED_COUNTS is None or EXPECTED_LEDGER_DIGEST is None
            else "LOCKED"
        )
    )
    print(
        "orientation=ordered proof-branch gauge only;"
        "heavy graph and overlap weights symmetric;not tournament"
    )
    print(
        "controls=source hashes;ordered K4 partition;Kruskal-vs-16-trees;"
        "literal sequential-vs-simultaneous;"
        "singleton containment-vs-literal-vs-scalar"
    )
    print("scope=exact K>=20 strongest-r3 H2<=40 slice;not LRC14")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
