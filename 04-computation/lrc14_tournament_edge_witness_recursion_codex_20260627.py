#!/usr/bin/env python3
"""HYP-3124/S271: tournament edge-witness class-deck stress audit.

This is a finite scout and synthesis ledger, not a proof of LRC14.

Core finite object:
  edge_witness(tail -> tip) =
    outside four-sector deck
    + tail-deletion child
    + tip-deletion child
    + optional cross-sector orientation word.

LRC14 translation:
  The open covering branch from HYP-3121 is the decorrelation floor
  Rprime = mu(R-safe and Q-lonely)/(mu(R-safe) mu(Q-lonely)) >= c
  for the few-apex/multi-far cases r=2..6.  This scout treats that floor as
  an information-retention stress test.  HYP-3125 is the canonical measured
  Rprime floor scout; this script checks which finite tournament edge
  coordinates must be retained before such a quotient is trusted.

Tournament Analysis declaration:
  vertices: proof carriers and two-ended witness packets, not runners;
  alternate vertex sets considered: runners, gaps, fixed circle sections,
    section boundaries, wall-crossing events, residues, cover arcs, Fourier
    modes, matroid topes/circuits, directed edges, and proof obligations;
  pairwise observable: retained LRC predicate axes for the multi-far floor;
  switch/gauge: carrier A beats B if it wins more axes; ties use the fixed
    Hamiltonian path listed in TIE_ORDER.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from itertools import combinations
from typing import Dict, Iterable, List, Sequence, Tuple

import a000568_edge_perspective_extension_codex_s213 as edge_lift
import tournament_rigidity_cascade_s589 as tour


SECTORS: Tuple[Tuple[str, Tuple[int, int]], ...] = (
    ("++", (1, 1)),  # tail -> x and tip -> x
    ("+-", (1, 0)),  # tail -> x -> tip
    ("-+", (0, 1)),  # tip -> x -> tail
    ("--", (0, 0)),  # x beats both endpoints
)
SECTOR_LABELS = tuple(label for label, _ in SECTORS)


def classes_by_n(n: int) -> Tuple[int, ...]:
    if n == 6:
        return edge_lift.generated_classes6()
    return tour.classes(n)


def oriented_edges(mask: int, n: int) -> List[Tuple[int, int]]:
    out: List[Tuple[int, int]] = []
    for a, b in combinations(range(n), 2):
        out.append((a, b) if tour.edge(mask, n, a, b) else (b, a))
    return out


def induced_canonical(mask: int, n: int, vertices: Sequence[int]) -> int:
    if len(vertices) <= 1:
        return 0
    return tour.canonical(tour.induced(mask, n, tuple(vertices)), len(vertices))


def opposite(mask: int, n: int) -> int:
    out = 0
    for i, j in combinations(range(n), 2):
        out = tour.set_edge(out, n, i, j, not tour.edge(mask, n, i, j))
    return out


def cyclic_triples(mask: int, n: int) -> int:
    total = 0
    for a, b, c in combinations(range(n), 3):
        total += int(
            (tour.edge(mask, n, a, b) and tour.edge(mask, n, b, c) and tour.edge(mask, n, c, a))
            or (
                tour.edge(mask, n, a, c)
                and tour.edge(mask, n, c, b)
                and tour.edge(mask, n, b, a)
            )
        )
    return total


def hamiltonian_paths(mask: int, n: int) -> int:
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for seen in range(1 << n):
        for v in range(n):
            if not dp[seen][v]:
                continue
            for u in range(n):
                if (seen & (1 << u)) or not tour.edge(mask, n, v, u):
                    continue
                dp[seen | (1 << u)][u] += dp[seen][v]
    return sum(dp[-1])


def sector_groups(mask: int, n: int, tail: int, tip: int) -> Dict[str, Tuple[int, ...]]:
    groups: Dict[str, List[int]] = {label: [] for label in SECTOR_LABELS}
    value_to_label = {value: label for label, value in SECTORS}
    for x in range(n):
        if x in (tail, tip):
            continue
        value = (int(tour.edge(mask, n, tail, x)), int(tour.edge(mask, n, tip, x)))
        groups[value_to_label[value]].append(x)
    return {label: tuple(groups[label]) for label in SECTOR_LABELS}


def sector_counts(groups: Dict[str, Tuple[int, ...]]) -> Tuple[Tuple[str, int], ...]:
    return tuple((label, len(groups[label])) for label in SECTOR_LABELS)


def sector_internal(mask: int, n: int, groups: Dict[str, Tuple[int, ...]]) -> Tuple[Tuple[str, int], ...]:
    return tuple((label, induced_canonical(mask, n, groups[label])) for label in SECTOR_LABELS)


def cross_sector_word(mask: int, n: int, groups: Dict[str, Tuple[int, ...]]) -> Tuple[Tuple[str, str, int, int], ...]:
    word: List[Tuple[str, str, int, int]] = []
    for i, left in enumerate(SECTOR_LABELS):
        for right in SECTOR_LABELS[i + 1 :]:
            wins = sum(1 for u in groups[left] for v in groups[right] if tour.edge(mask, n, u, v))
            word.append((left, right, wins, len(groups[left]) * len(groups[right])))
    return tuple(word)


def child(mask: int, n: int, victim: int) -> int:
    if n <= 1:
        return 0
    return tour.canonical(tour.delete_vertex(mask, n, victim), n - 1)


def edge_signature(mask: int, n: int, tail: int, tip: int, mode: str) -> Tuple[object, ...]:
    groups = sector_groups(mask, n, tail, tip)
    counts = sector_counts(groups)
    internal = sector_internal(mask, n, groups)
    cross = cross_sector_word(mask, n, groups)
    tail_child = child(mask, n, tail)
    tip_child = child(mask, n, tip)

    if mode == "sector_counts":
        return (counts,)
    if mode == "sector_internal":
        return counts, internal
    if mode == "recursive_children":
        return counts, internal, ("tail_child", tail_child), ("tip_child", tip_child)
    if mode == "roleless_children":
        return counts, internal, ("children", tuple(sorted((tail_child, tip_child))))
    if mode == "full_edge_witness":
        return counts, internal, cross, ("tail_child", tail_child), ("tip_child", tip_child)
    raise ValueError(mode)


def class_deck(mask: int, n: int, mode: str) -> Tuple[Tuple[object, ...], ...]:
    return tuple(sorted(edge_signature(mask, n, tail, tip, mode) for tail, tip in oriented_edges(mask, n)))


def class_deck_stats(n: int) -> List[Tuple[str, int, int, int, Dict[int, int]]]:
    reps = classes_by_n(n)
    rows: List[Tuple[str, int, int, int, Dict[int, int]]] = []
    for mode in MODES:
        fibers = Counter(class_deck(rep, n, mode) for rep in reps)
        hist = dict(sorted(Counter(fibers.values()).items()))
        rows.append((mode, len(fibers), len(reps) - len(fibers), max(fibers.values()), hist))
    return rows


def edge_type_counts(n: int) -> List[Tuple[str, int]]:
    reps = classes_by_n(n)
    rows: List[Tuple[str, int]] = []
    for mode in MODES:
        types = {
            edge_signature(rep, n, tail, tip, mode)
            for rep in reps
            for tail, tip in oriented_edges(rep, n)
        }
        rows.append((mode, len(types)))
    return rows


def collision_examples(n: int, mode: str, limit: int = 4) -> List[str]:
    reps = classes_by_n(n)
    fibers: Dict[Tuple[Tuple[object, ...], ...], List[int]] = defaultdict(list)
    for rep in reps:
        fibers[class_deck(rep, n, mode)].append(rep)
    examples: List[str] = []
    for group in fibers.values():
        if len(group) <= 1:
            continue
        group = sorted(group)
        converse_flags = []
        for a, b in combinations(group, 2):
            converse_flags.append(tour.canonical(opposite(a, n), n) == b)
        summary = []
        for rep in group[:4]:
            summary.append(
                f"{rep}:score={tour.score_sequence(rep, n)},c3={cyclic_triples(rep, n)},H={hamiltonian_paths(rep, n)}"
            )
        examples.append(
            f"fiber_size={len(group)} converse_pair_seen={any(converse_flags)} reps=["
            + "; ".join(summary)
            + "]"
        )
        if len(examples) >= limit:
            break
    return examples


def recursion_repair_audit(n: int) -> List[str]:
    reps = classes_by_n(n)
    sector_fibers: Dict[Tuple[object, ...], List[Tuple[int, int, int]]] = defaultdict(list)
    recursive_splits = 0
    cross_splits_after_recursion = 0
    max_sector_fiber = 0
    max_recursive_fiber = 0

    for rep in reps:
        for tail, tip in oriented_edges(rep, n):
            key = edge_signature(rep, n, tail, tip, "sector_internal")
            sector_fibers[key].append((rep, tail, tip))

    for instances in sector_fibers.values():
        if len(instances) <= 1:
            continue
        max_sector_fiber = max(max_sector_fiber, len(instances))
        rec_keys = {
            edge_signature(rep, n, tail, tip, "recursive_children")
            for rep, tail, tip in instances
        }
        full_keys = {
            edge_signature(rep, n, tail, tip, "full_edge_witness")
            for rep, tail, tip in instances
        }
        recursive_splits += int(len(rec_keys) > 1)
        cross_splits_after_recursion += int(len(full_keys) > len(rec_keys))
        rec_fibers = Counter(
            edge_signature(rep, n, tail, tip, "recursive_children")
            for rep, tail, tip in instances
        )
        max_recursive_fiber = max(max_recursive_fiber, max(rec_fibers.values()))

    return [
        f"n={n} sector_internal_edge_fibers={len(sector_fibers)}",
        f"n={n} nontrivial_sector_fibers={sum(1 for v in sector_fibers.values() if len(v) > 1)}",
        f"n={n} sector_fibers_split_by_tail_tip_children={recursive_splits}",
        f"n={n} recursive_fibers_split_by_cross_sector_word={cross_splits_after_recursion}",
        f"n={n} max_sector_fiber={max_sector_fiber} max_recursive_fiber_after_child_pair={max_recursive_fiber}",
    ]


AXES = (
    "lrc_predicate_retention",
    "multi_far_floor_adjacency",
    "two_ended_recursion",
    "destroyed_coordinate_control",
    "average_distribution_power",
    "zero_free_or_positive_correlation",
    "finite_checkability",
    "desmoothing_or_reconstruction",
    "formalization_readiness",
)


@dataclass(frozen=True)
class Carrier:
    name: str
    kind: str
    scores: Dict[str, int]
    preserves: str
    destroys: str
    next_test: str


def carrier(
    name: str,
    kind: str,
    values: Sequence[int],
    preserves: str,
    destroys: str,
    next_test: str,
) -> Carrier:
    assert len(values) == len(AXES)
    return Carrier(name, kind, dict(zip(AXES, values)), preserves, destroys, next_test)


CARRIERS: Tuple[Carrier, ...] = (
    carrier(
        "edge_tail_tip_recursive_witness",
        "two_ended_tournament_packet",
        (5, 5, 5, 5, 3, 3, 5, 4, 4),
        "the overlap predicate as a two-ended witness with both marginal deletion children",
        "exact time coordinate and high Fourier phase unless cross-sector or analytic sidecar is retained",
        "add edge_witness_id, tail_child, tip_child, and cross_sector_orientation_word to few-apex rows",
    ),
    carrier(
        "cross_sector_orientation_word",
        "observer_cut_sidecar",
        (4, 4, 4, 5, 2, 2, 5, 4, 4),
        "the first concrete orientation coordinate that repairs edge-sector quotient collapse",
        "the analytic size of the overlap if used without measure or Fourier data",
        "test whether any multi-far packet collision survives sector counts plus orientation",
    ),
    carrier(
        "l2_cauchy_schwarz_tail_floor",
        "existing_spectral_floor",
        (5, 5, 2, 4, 4, 4, 4, 4, 4),
        "a rigorous-style low-frequency exact sum plus Parseval tail bound for Rprime",
        "local endpoint ownership and the reason a bad frequency appears",
        "turn lrc14_spectrum_L2tail_synthesis_kpswf12 into a finite HYP-3121 row certificate",
    ),
    carrier(
        "elliott_halberstam_level_packet",
        "average_distribution_schema",
        (4, 5, 2, 3, 5, 3, 3, 3, 3),
        "average residue-channel discrepancy control, analogous to level of distribution one",
        "worst individual residue/channel behavior and exact witness time",
        "define finite channel errors over q, residue, edge-witness, and Fourier shells; prove an averaged bound",
    ),
    carrier(
        "gaussian_heat_kernel_smoothing",
        "smooth_correlation_floor",
        (4, 5, 2, 3, 4, 4, 4, 5, 3),
        "a smoothed positive correlation floor with exponentially damped Fourier tails",
        "sharp boundary atoms unless a finite-ruler desmoothing margin is retained",
        "prove smoothed Rprime_t >= c_t and desmooth through finite denominator thresholds",
    ),
    carrier(
        "asano_endpoint_contraction",
        "multiaffine_zero_free_contraction",
        (4, 4, 5, 4, 3, 5, 3, 3, 3),
        "endpoint/tail-tip interaction variables while preserving a zero-free sector under contraction",
        "metric size of the safe arcs unless coefficients are tied back to Phi/P activation",
        "build a two-end partition polynomial and check whether contraction preserves nonvanishing",
    ),
    carrier(
        "endpoint_phi_activation_gate",
        "exact_gap_circuit_gate",
        (5, 4, 2, 4, 3, 4, 4, 4, 4),
        "the exact endpoint-cover gap sum and max activation from HYP-2108/HYP-2112",
        "distribution across residue channels if endpoint activations are scalarized",
        "join Phi_gap, P_max_activation, and edge_witness_id in one packet schema",
    ),
    carrier(
        "normal_fan_cech_finite_ruler",
        "component_topology_floor",
        (5, 5, 1, 4, 2, 2, 4, 4, 4),
        "normalized component and finite-ruler payload for the few-apex floor",
        "tail/tip covariance and Fourier phases if topology is isolated",
        "combine Cech beta, owner currents, and Gaussian/L2 floor certificates",
    ),
    carrier(
        "raw_bonferroni_scalar",
        "negative_control",
        (1, 2, 0, 1, 1, 1, 5, 1, 5),
        "only the measure-sum obstruction that already fails on few-apex examples",
        "all correlation structure; it cannot see Rprime when mu(A)+mu(B)<1",
        "keep only as a failure alarm, never as a closure route",
    ),
    carrier(
        "raw_edge_count_shadow",
        "negative_control",
        (1, 1, 1, 1, 0, 0, 5, 1, 5),
        "small tournament size telemetry",
        "tail/tip role, sector orientation, deletion recursion, and LRC predicate",
        "use only to sanity-check finite enumeration",
    ),
)


TIE_ORDER = tuple(c.name for c in CARRIERS)
MODES = (
    "sector_counts",
    "sector_internal",
    "roleless_children",
    "recursive_children",
    "full_edge_witness",
)


def carrier_adj() -> List[List[int]]:
    n = len(CARRIERS)
    adj = [[0] * n for _ in range(n)]
    order = {name: i for i, name in enumerate(TIE_ORDER)}
    for i, a in enumerate(CARRIERS):
        for j, b in enumerate(CARRIERS):
            if i == j:
                continue
            wins_a = sum(a.scores[axis] > b.scores[axis] for axis in AXES)
            wins_b = sum(b.scores[axis] > a.scores[axis] for axis in AXES)
            if wins_a > wins_b or (wins_a == wins_b and order[a.name] < order[b.name]):
                adj[i][j] = 1
    return adj


def directed_three_cycles_adj(adj: List[List[int]]) -> int:
    total = 0
    for a, b, c in combinations(range(len(adj)), 3):
        total += int(adj[a][b] and adj[b][c] and adj[c][a])
        total += int(adj[a][c] and adj[c][b] and adj[b][a])
    return total


def scc_sizes(adj: List[List[int]]) -> List[int]:
    def reach(start: int, graph: List[List[int]]) -> set[int]:
        seen = {start}
        q = deque([start])
        while q:
            v = q.popleft()
            for u, bit in enumerate(graph[v]):
                if bit and u not in seen:
                    seen.add(u)
                    q.append(u)
        return seen

    reverse = [[adj[j][i] for j in range(len(adj))] for i in range(len(adj))]
    remaining = set(range(len(adj)))
    sizes = []
    while remaining:
        v = min(remaining)
        comp = reach(v, adj) & reach(v, reverse)
        sizes.append(len(comp))
        remaining -= comp
    return sorted(sizes, reverse=True)


def hamiltonian_paths_adj(adj: List[List[int]]) -> int:
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for seen in range(1 << n):
        for v in range(n):
            if not dp[seen][v]:
                continue
            for u in range(n):
                if (seen & (1 << u)) or not adj[v][u]:
                    continue
                dp[seen | (1 << u)][u] += dp[seen][v]
    return sum(dp[-1])


def selected_path(adj: List[List[int]]) -> List[str]:
    remaining = set(range(len(adj)))
    path: List[int] = []
    while remaining:
        best = max(
            remaining,
            key=lambda i: (
                sum(adj[i][j] for j in remaining if j != i),
                -TIE_ORDER.index(CARRIERS[i].name),
            ),
        )
        path.append(best)
        remaining.remove(best)
    return [CARRIERS[i].name for i in path]


def print_edge_tables() -> None:
    print("1. EDGE-WITNESS RECURSION ON SMALL TOURNAMENT CLASSES")
    print("n  U(n)  mode                  unique_class_decks  collisions  max_fiber  fiber_hist")
    for n in (3, 4, 5, 6):
        reps = classes_by_n(n)
        for mode, unique, collisions, max_fiber, hist in class_deck_stats(n):
            print(f"{n}  {len(reps):4d}  {mode:21s} {unique:18d} {collisions:11d} {max_fiber:10d} {hist}")
        print("")

    print("2. INDIVIDUAL DIRECTED-EDGE TYPE COUNTS")
    print("n  mode                  unique_edge_types")
    for n in (3, 4, 5, 6):
        for mode, count in edge_type_counts(n):
            print(f"{n}  {mode:21s} {count:17d}")
        print("")

    print("3. WHERE RECURSION AND CROSS-SECTOR WORDS REPAIR EDGE QUOTIENTS")
    for n in (5, 6):
        for line in recursion_repair_audit(n):
            print(line)
    print("")

    print("4. CLASS-DECK COLLISION EXAMPLES")
    for mode in ("sector_internal", "recursive_children", "full_edge_witness"):
        examples = collision_examples(6, mode)
        print(f"mode={mode} collision_examples={len(examples)}")
        for ex in examples:
            print("  " + ex)
    print("")


def print_carrier_tournament() -> None:
    adj = carrier_adj()
    scores = {CARRIERS[i].name: sum(adj[i]) for i in range(len(adj))}
    print("5. HYP-3125 HANDOFF CARRIER TOURNAMENT")
    print(f"vertices={len(CARRIERS)} axes={list(AXES)}")
    print(f"score_hist={dict(sorted(Counter(scores.values()).items()))}")
    print(f"directed_3cycles={directed_three_cycles_adj(adj)}")
    print(f"scc_sizes={scc_sizes(adj)}")
    print(f"hamiltonian_path_count={hamiltonian_paths_adj(adj)}")
    print("selected_path=" + " > ".join(selected_path(adj)))
    print("")
    print("Carrier ledger:")
    for c in sorted(CARRIERS, key=lambda x: (-scores[x.name], TIE_ORDER.index(x.name))):
        print(f"- {c.name} [{c.kind}] score={scores[c.name]}")
        print(f"  preserves: {c.preserves}")
        print(f"  destroys/audits: {c.destroys}")
        print(f"  next_test: {c.next_test}")
    print("")


def print_floor_program() -> None:
    print("6. HYP-3125 HANDOFF ROUTES FOR THE UNIFORM Rprime >= c FLOOR")
    print(
        "A. Edge-witness route: encode R-safe and Q-lonely as the two ends of a "
        "witness edge.  Tail and tip deletion children are the marginals; the "
        "cross-sector word is the covariance payload.  A collapse to empty "
        "intersection is legal only if a child or cross-sector coordinate names "
        "the debt."
    )
    print(
        "B. Elliott-Halberstam-style route: do not ask for every residue/fiber "
        "to be uniform.  Ask for a level-of-distribution average over channels "
        "(q, residue, edge_witness_id, Fourier shell).  This is an analogy, not "
        "a prime theorem import: the target is a finite discrepancy inequality "
        "whose averaged error is below the product mass."
    )
    print(
        "C. Gaussian route: heat-smooth the two indicator events, prove a "
        "smoothed positive correlation using damped Fourier tails, then "
        "desmooth through a finite-ruler margin.  This is the analytic sibling "
        "of the existing L2-Cauchy-Schwarz tail certificate."
    )
    print(
        "D. Asano route: make endpoint/tail-tip interactions multiaffine, then "
        "contract the two endpoint variables while preserving a zero-free sector. "
        "If the contracted polynomial remains nonvanishing in the legal sector "
        "and coefficients still read back to Phi/P activation, empty overlap is "
        "forbidden or routed to named debt."
    )
    print("")
    print("Assumption challenge:")
    print(
        "Chosen vertices here are edge witnesses, residue/Fourier channels, and "
        "proof carriers.  The preserved predicate for the handoff is the LRC14 "
        "covering witness R-safe intersect Q-lonely nonempty, equivalently a "
        "positive Rprime floor in HYP-3125.  The quotient destroys exact time, "
        "endpoint owner, phase, and sometimes tail/tip role.  Those coordinates "
        "must be retained, averaged with a finite error budget, resurrected by "
        "Gaussian/Asano/endpoint sidecars, or emitted as named HYP-2968/HYP-3121 "
        "debt before the HYP-3125 edge_floor_packet compresses them."
    )


def main() -> None:
    print("HYP-3124/S271 EDGE-WITNESS CLASS-DECK STRESS AUDIT")
    print("script=04-computation/lrc14_tournament_edge_witness_recursion_codex_20260627.py")
    print("status=finite scout + synthesis; not a proof")
    print("")
    print_edge_tables()
    print_carrier_tournament()
    print_floor_program()


if __name__ == "__main__":
    main()
