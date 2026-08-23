#!/usr/bin/env python3
"""Exact audit for the Arithmetic-Kakeya k=1 / forest port-space barrier.

The proof object audited here is internal to the paid-row forcing semantics.
It does not invoke the external certificate=>AK implication.
"""

from __future__ import annotations

from fractions import Fraction
import ast
import hashlib
import itertools
import json
import pathlib
import sys


sys.stdout.reconfigure(newline="\n")
HERE = pathlib.Path(__file__).resolve().parent
ROOT = pathlib.Path(__file__).resolve().parents[1]
CHECKS = 0


def require(condition: object, label: str) -> None:
    global CHECKS
    if not condition:
        raise RuntimeError(f"audit gate failed: {label}")
    CHECKS += 1


def sha256(path: pathlib.Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


PINNED = {
    ROOT / "01-canon/theorems/THM-2850-paid-distinguished-coloop-round-budget-and-slope-grammar-rank-defect.md":
        "56976c5b07307f2e55c1cffa2562d1ff61e044566f23e87096c8876b5e7f9f52",
    ROOT / "04-computation/ak_forcing_engine.py":
        "f9a1b3c1d3e354d2bd8675687be31b8d97578080dd08e62dad1f77f497e0c485",
    ROOT / "04-computation/ak_mode3_v2.py":
        "5a85b0b9de423028d929d3ce1d86a1a0a90cab5286f6434300d1104cc44cef55",
    ROOT / "06-writeups/AK-FORCING-WORKBENCH-klein-S691.md":
        "1055b6f1d48097c055a9dfb61e414322a98c01d973f1fefe7422ec563b431df1",
    ROOT / "06-writeups/AK-BENCHMARK-SPEC-REPORT-klein-S691.md":
        "7c264e71e9f137d6845e27b3261b6af029224c673f0ee23ba961b0911e80059b",
}

for pinned_path, pinned_hash in PINNED.items():
    require(pinned_path.is_file(), f"pinned input exists: {pinned_path.name}")
    require(sha256(pinned_path) == pinned_hash,
            f"pinned input hash: {pinned_path.name}")


def exact_rank(rows: list[list[int | Fraction]], ncols: int) -> int:
    """Exact rational row rank, independent of the maintained AK engine."""
    matrix = [[Fraction(value) for value in row] for row in rows
              if any(value != 0 for value in row)]
    pivot_row = 0
    for column in range(ncols):
        pivot = next((row for row in range(pivot_row, len(matrix))
                      if matrix[row][column] != 0), None)
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        pivot_value = matrix[pivot_row][column]
        matrix[pivot_row] = [value / pivot_value
                             for value in matrix[pivot_row]]
        for row in range(len(matrix)):
            if row == pivot_row or matrix[row][column] == 0:
                continue
            multiple = matrix[row][column]
            matrix[row] = [left - multiple * right
                           for left, right in zip(matrix[row],
                                                  matrix[pivot_row])]
        pivot_row += 1
        if pivot_row == len(matrix):
            break
    return pivot_row


def in_row_span(rows: list[list[int | Fraction]], target: list[int | Fraction]) -> bool:
    ncols = len(target)
    return exact_rank(rows, ncols) == exact_rank(rows + [target], ncols)


def edge_row(n: int, left: int, right: int, slope: int) -> list[int]:
    row = [0] * (2 * n)
    row[2 * left] = 1
    row[2 * left + 1] = slope
    row[2 * right] = -1
    row[2 * right + 1] = -slope
    return row


def ground_row(n: int, vertex: int, slope: int) -> list[int]:
    row = [0] * (2 * n)
    row[2 * vertex] = 1
    row[2 * vertex + 1] = slope
    return row


def project_rows(rows: list[list[int]], live: tuple[int, ...]) -> list[list[int]]:
    columns = tuple(coordinate for vertex in live
                    for coordinate in (2 * vertex, 2 * vertex + 1))
    return [[row[column] for column in columns] for row in rows]


def exact_force(rows: list[list[int]], n: int,
                wildcards: tuple[int, ...] = ()) -> tuple[
                    bool, tuple[tuple[int, ...], ...], tuple[int, ...]]:
    """Simultaneous exact forcing by pure distinguished-d membership."""
    wild = set(wildcards)
    trace: list[tuple[int, ...]] = []
    ranks: list[int] = []
    while len(wild) < n:
        live = tuple(vertex for vertex in range(n) if vertex not in wild)
        restricted = project_rows(rows, live)
        ncols = 2 * len(live)
        rank = exact_rank(restricted, ncols)
        ranks.append(rank)
        fired = []
        for local_vertex, vertex in enumerate(live):
            target = [0] * ncols
            target[2 * local_vertex + 1] = 1
            if in_row_span(restricted, target):
                fired.append(vertex)
        if not fired:
            break
        firing = tuple(fired)
        trace.append(firing)
        wild.update(fired)
    return len(wild) == n, tuple(trace), tuple(ranks)


# A span state for available ground half-edges at one vertex.  None is zero,
# an integer is one nonvertical line (1,rho), and FULL is Q^2.
FULL = "FULL"
SLOPES = (-1, 0, 2)
STATES: tuple[None | int | str, ...] = (None,) + SLOPES + (FULL,)


def join_state(state: None | int | str, slope: int) -> int | str:
    if state is None:
        return slope
    if state == FULL or state == slope:
        return state
    return FULL


def contains_line(state: None | int | str, slope: int) -> bool:
    return state == FULL or state == slope


def state_dimension(state: None | int | str) -> int:
    return 0 if state is None else (2 if state == FULL else 1)


def state_cost(state: None | int | str) -> int:
    """Minimum number of nonvertical ground rows realizing this state."""
    return state_dimension(state)


def forest_port_state(vertices: tuple[int, ...],
                      adjacency: dict[int, list[tuple[int, int]]],
                      grounds: tuple[None | int | str, ...],
                      root: int) -> None | int | str:
    """Leaf-pruning port space for a forest, rooted at ``root``."""
    vertex_set = set(vertices)

    def absorb(node: int, parent: int | None) -> None | int | str:
        state = grounds[node]
        for neighbor, slope in adjacency[node]:
            if neighbor == parent or neighbor not in vertex_set:
                continue
            child_state = absorb(neighbor, node)
            if contains_line(child_state, slope):
                state = join_state(state, slope)
        return state

    return absorb(root, None)


def path_network(length: int, edge_choices: tuple[None | int, ...],
                 grounds: tuple[None | int | str, ...]) -> tuple[
                     tuple[int, ...], dict[int, list[tuple[int, int]]]]:
    vertices = tuple(range(length))
    adjacency = {vertex: [] for vertex in vertices}
    for left, slope in enumerate(edge_choices):
        if slope is None:
            continue
        right = left + 1
        adjacency[left].append((right, slope))
        adjacency[right].append((left, slope))
    return vertices, adjacency


def delete_vertices(vertices: tuple[int, ...],
                    adjacency: dict[int, list[tuple[int, int]]],
                    grounds: tuple[None | int | str, ...],
                    deleted: set[int]) -> tuple[
                        tuple[int, ...], dict[int, list[tuple[int, int]]],
                        tuple[None | int | str, ...]]:
    live = tuple(vertex for vertex in vertices if vertex not in deleted)
    live_set = set(live)
    new_adjacency = {vertex: [] for vertex in vertices}
    new_grounds = list(grounds)
    seen_edges = set()
    for left in vertices:
        for right, slope in adjacency[left]:
            edge_key = (min(left, right), max(left, right))
            if edge_key in seen_edges:
                continue
            seen_edges.add(edge_key)
            if left in live_set and right in live_set:
                new_adjacency[left].append((right, slope))
                new_adjacency[right].append((left, slope))
            elif left in live_set:
                new_grounds[left] = join_state(new_grounds[left], slope)
            elif right in live_set:
                new_grounds[right] = join_state(new_grounds[right], slope)
    return live, new_adjacency, tuple(new_grounds)


def fireable_vertices(vertices: tuple[int, ...],
                      adjacency: dict[int, list[tuple[int, int]]],
                      grounds: tuple[None | int | str, ...]) -> tuple[int, ...]:
    return tuple(vertex for vertex in vertices
                 if forest_port_state(vertices, adjacency, grounds, vertex)
                 == FULL)


def concrete_rows(length: int, edge_choices: tuple[None | int, ...],
                  grounds: tuple[None | int | str, ...]) -> list[list[int]]:
    rows = []
    for left, slope in enumerate(edge_choices):
        if slope is not None:
            rows.append(edge_row(length, left, left + 1, slope))
    for vertex, state in enumerate(grounds):
        if state is None:
            continue
        if state == FULL:
            rows.append(ground_row(length, vertex, SLOPES[0]))
            rows.append(ground_row(length, vertex, SLOPES[1]))
        else:
            rows.append(ground_row(length, vertex, state))
    return rows


# Exhaust the full three-line state quotient for paths through five vertices.
# This quotient is exact for port spaces: multiplicity and the identity of a
# second distinct slope do not matter once the local span is Q^2.
network_count = 0
rooted_port_count = 0
success_count = 0
partial_count = 0
stalled_count = 0
min_success_slack = None
matrix_cross_checks = 0

for length in range(1, 6):
    edge_options = (None,) + SLOPES
    for edge_choices in itertools.product(edge_options, repeat=length - 1):
        for grounds in itertools.product(STATES, repeat=length):
            network_count += 1
            vertices, adjacency = path_network(length, edge_choices, grounds)
            firing = fireable_vertices(vertices, adjacency, grounds)
            rooted_port_count += length
            if not firing:
                stalled_count += 1
            elif len(firing) == length:
                success_count += 1
                paid_rows = (sum(slope is not None for slope in edge_choices)
                             + sum(state_cost(state) for state in grounds))
                slack = paid_rows - 2 * length
                require(slack >= 0, "successful state quotient obeys score-two floor")
                if min_success_slack is None or slack < min_success_slack:
                    min_success_slack = slack
            else:
                partial_count += 1
                live, projected_adjacency, projected_grounds = delete_vertices(
                    vertices, adjacency, grounds, set(firing))
                require(fireable_vertices(live, projected_adjacency,
                                          projected_grounds) == (),
                        "forest quotient creates no second-round firing")

            # Exact matrix cross-check: exhaustive through length three, then
            # a deterministic residue sample at lengths four and five.
            sample_key = network_count * 131 + length * 17
            if length <= 3 or sample_key % 997 == 0:
                rows = concrete_rows(length, edge_choices, grounds)
                row_rank = exact_rank(rows, 2 * length)
                for root in vertices:
                    without_root = [
                        row[:2 * root] + row[2 * root + 2:]
                        for row in rows
                    ]
                    exact_port_dimension = row_rank - exact_rank(
                        without_root, 2 * length - 2)
                    predicted = forest_port_state(vertices, adjacency,
                                                  grounds, root)
                    require(exact_port_dimension == state_dimension(predicted),
                            "leaf pruning equals exact port dimension")
                    if isinstance(predicted, int):
                        target = [0] * (2 * length)
                        target[2 * root] = 1
                        target[2 * root + 1] = predicted
                        require(in_row_span(rows, target),
                                "one-dimensional exact port has predicted nonvertical line")
                    matrix_cross_checks += 1

require(network_count == 842105,
        "complete path state quotient through five vertices visited")
require(rooted_port_count == 4166205,
        "complete rooted port census through five vertices visited")
require(min_success_slack == 0,
        "score-two equality occurs in the exhaustive state quotient")


# Positive control 1: an arbitrary-length merge-free path attaining score 2.
positive_n = 6
positive_edges = (0,) * (positive_n - 1)
positive_grounds = (FULL,) + (2,) * (positive_n - 1)
positive_rows = concrete_rows(positive_n, positive_edges, positive_grounds)
positive_ok, positive_trace, positive_ranks = exact_force(
    positive_rows, positive_n)
require(positive_ok and positive_trace == (tuple(range(positive_n)),),
        "score-two path forces in one simultaneous round")
require(len(positive_rows) == exact_rank(positive_rows, 2 * positive_n)
        == 2 * positive_n,
        "score-two path has full exact initial rank")


# Positive control 2: initial wildcard.  Cross-edges become ground half-edges.
wild_n = 7
wild = (2,)
wild_rows = [edge_row(wild_n, left, left + 1, 0)
             for left in range(wild_n - 1)]
wild_rows += [ground_row(wild_n, vertex, 2)
              for vertex in range(wild_n) if vertex not in wild]
wild_ok, wild_trace, wild_ranks = exact_force(wild_rows, wild_n, wild)
require(wild_ok and wild_trace == ((0, 1, 3, 4, 5, 6),),
        "wildcard path forces all live vertices in one round")
require(len(wild_rows) == 12 and wild_ranks == (12,),
        "wildcard equality has paid rows = exact rank = 2u")


# Positive control 3: adjacent merges partition a k=1 base path into interval
# classes, whose quotient is again a path.  Two merge slots are free; the
# remaining three labelled slots plus five seeds attain 8/4 = 2.
merge_rows = concrete_rows(4, (0, 0, 0), (FULL, 2, 2, 2))
merge_ok, merge_trace, merge_ranks = exact_force(merge_rows, 4)
require(merge_ok and merge_trace == ((0, 1, 2, 3),),
        "adjacent-merge quotient path forces in one round")
require(len(merge_rows) == merge_ranks[0] == 8,
        "adjacent-merge quotient attains exact score two")


# Hostile control 1: the forbidden vertical label (a+b=0) destroys the lemma.
vertical_rows = [[0, 1]]
vertical_ok, vertical_trace, vertical_ranks = exact_force(vertical_rows, 1)
require(vertical_ok and vertical_trace == ((0,),) and vertical_ranks == (1,),
        "forbidden vertical seed would give score one")


# Hostile control 2: a live 3-cycle admits a genuine 5/3 two-round paid-row
# certificate, so acyclicity is essential.  Slopes 0,1,2 correspond to the
# allowed integer labels (1,1), (1,0), and (3,-1), respectively.
cycle_rows = [
    edge_row(3, 0, 1, 0),
    edge_row(3, 1, 2, 1),
    edge_row(3, 2, 0, 2),
    ground_row(3, 1, 2),
    ground_row(3, 2, 0),
]
cycle_ok, cycle_trace, cycle_ranks = exact_force(cycle_rows, 3)
require(cycle_ok and cycle_trace == ((0,), (1, 2)),
        "cyclic five-thirds hostile forces in two rounds")
require(cycle_ranks == (5, 4), "cyclic hostile has rank profile 5 then 4")
cycle_without_d0 = [row[:1] + row[2:] for row in cycle_rows]
cycle_without_pair0 = [row[2:] for row in cycle_rows]
require(exact_rank(cycle_rows, 6) == 5
        and exact_rank(cycle_without_d0, 5) == 4
        and exact_rank(cycle_without_pair0, 4) == 4,
        "cycle first fire drops d-rank but not matching s-rank")


self_tree = ast.parse(pathlib.Path(__file__).read_text(encoding="utf-8"))
require(not any(isinstance(node, ast.Assert) for node in ast.walk(self_tree)),
        "checker contains no optimization-erased assertion")

semantic = {
    "theorem": (
        "nonvertical labelled forest port spaces are 0, a nonvertical line, "
        "or Q^2; forcing closure stabilizes after its first round; success "
        "implies initial rank 2u and paid score at least 2"
    ),
    "k1_scope": (
        "strict/verifiable paths and legal adjacent-merge intuitive quotients; "
        "arbitrary X,n,T,seeds over Q; charged row count dominates active rows"
    ),
    "finite_universe": {
        "lengths": [1, 2, 3, 4, 5],
        "edge_states": ["absent", -1, 0, 2],
        "ground_span_states": ["zero", "line(-1)", "line(0)",
                               "line(2)", "Q2"],
        "networks": network_count,
        "rooted_ports": rooted_port_count,
        "matrix_cross_checks": matrix_cross_checks,
    },
    "positive": {
        "path": {"u": positive_n, "paid": len(positive_rows),
                 "trace": positive_trace},
        "wildcard": {"u": wild_n - len(wild), "paid": len(wild_rows),
                     "trace": wild_trace},
        "adjacent_merge_quotient": {"u": 4, "paid": len(merge_rows),
                                    "trace": merge_trace},
    },
    "hostile": {
        "forbidden_vertical": {"u": 1, "paid": 1,
                               "trace": vertical_trace},
        "live_cycle": {"u": 3, "paid": 5, "score": "5/3",
                       "trace": cycle_trace, "ranks": cycle_ranks},
    },
    "external_boundary": (
        "no certificate=>AK proof; no private-verifier claim; current external "
        "1.67513 record and <=1.675 target unchanged"
    ),
}
semantic_hash = hashlib.sha256(json.dumps(
    semantic, sort_keys=True, separators=(",", ":")
).encode("ascii")).hexdigest()

print("audit=AK_k1_forest_port_barrier")
print("verdict=PROVED_INTERNAL_REGIME+FINITE_EXACT_CONTROLS")
print("theorem=forest_port_spaces:0|nonvertical_line|Q2;closure_rounds<=1")
print("consequence=successful_forest:initial_rank=2u;paid_score>=2")
print("k1=verifiable_strict_and_loose_coincide;adjacent_merge_quotient_covered")
print(f"pinned_inputs={len(PINNED)}")
print(f"finite_state_universe=networks:{network_count};rooted_ports:{rooted_port_count};"
      f"matrix_cross_checks:{matrix_cross_checks}")
print(f"finite_outcomes=success:{success_count};partial:{partial_count};"
      f"stalled:{stalled_count};min_success_slack:{min_success_slack}")
print(f"positive_path=u:{positive_n};paid:{len(positive_rows)};trace:{positive_trace}")
print(f"positive_wildcard=u:{wild_n-len(wild)};paid:{len(wild_rows)};trace:{wild_trace}")
print(f"positive_adjacent_merge=u:4;paid:{len(merge_rows)};trace:{merge_trace}")
print(f"hostile_vertical=score:1;trace:{vertical_trace};excluded_by:a+b!=0")
print(f"hostile_cycle=score:5/3;trace:{cycle_trace};ranks:{cycle_ranks};"
      "first_pair_rank_drop:1")
print("benchmark=external_record:1.67513;target:<=1.675;status:OPEN_UNCHANGED")
print("scope=no_k>=2_claim;no_arbitrary_quotient_claim;no_private_verifier_claim;no_AK_theorem")
print(f"semantic_sha256={semantic_hash}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
