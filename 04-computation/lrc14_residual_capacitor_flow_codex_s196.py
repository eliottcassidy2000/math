#!/usr/bin/env python3
"""Residual capacitor flow scout for LRC14.

The recent HYP-3027..HYP-3031 stack leaves two instructive open-route
collisions after cheap proof gates:

* exact M+q still identifies one petal packet with one covering packet;
* coarse boundary topology still identifies one K33 packet with one covering
  packet.

This script treats each collision as a residual capacitor: two open packets
hold the same charge under a coarse quotient, and the next sidecar is a cut
that discharges the route ambiguity.  The computation is deliberately tiny and
theorem-facing.  It recomputes only the four packets involved, attaches the
existing topology, stalk, and fusion sidecars, and reports the first usable cut.

Tournament Analysis declaration:
  vertices: residual-cut carriers, not runners or raw arcs;
  pairwise observable: how many residual capacitors are separated, route/status
    purity, retained magnitude, topology, stalk data, packet labels,
    non-circularity, and proof cost;
  switch/gauge: lexicographic comparison of those observables;
  tie Hamiltonian path: the STAGES order below.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from pathlib import Path
import argparse
import sys


REPO = Path(__file__).resolve().parents[1]


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    assert spec and spec.loader
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


lp = load_module(
    "lrc14_labelled_packet_classifier_s196",
    REPO / "04-computation" / "lrc14_labelled_packet_counterexample_classifier_codex_20260624.py",
)
ac = load_module(
    "lrc14_arc_cech_s196",
    REPO / "04-computation" / "lrc14_arc_cech_nerve_carrier_codex_s174.py",
)
fusion = load_module(
    "lrc14_fusion_s196",
    REPO / "04-computation" / "lrc14_carrier_fusion_switchboard_codex_s189.py",
)
stalk_mod = load_module(
    "lrc14_stalk_s196",
    REPO / "04-computation" / "lrc14_safe_component_stalk_codex_s193.py",
)


TARGET_PAIRS = (
    (
        "M_q_petal_covering_capacitor",
        "shares word_plus_M_q; petal route versus covering route",
        (
            "two drop(10, 13)->add(20, 26)",
            "two drop(8, 12)->add(16, 24)",
        ),
    ),
    (
        "boundary_topology_k33_covering_capacitor",
        "shares word_plus_boundary_topology; K33 state lift versus covering route",
        (
            "two drop(12, 13)->add(26, 36)",
            "single swap 12->72",
        ),
    ),
)

TARGET_NAMES = frozenset(name for _pid, _note, names in TARGET_PAIRS for name in names)
PAIR_BY_NAME = {name: pair_id for pair_id, _note, names in TARGET_PAIRS for name in names}


def fmt(fr: Fraction) -> str:
    return str(fr.numerator) if fr.denominator == 1 else f"{fr.numerator}/{fr.denominator}"


def mu_bucket(mu: Fraction) -> str:
    if mu == 0:
        return "0"
    if mu <= Fraction(1, 1000):
        return "(0,1e-3]"
    if mu <= Fraction(1, 100):
        return "(1e-3,1e-2]"
    return ">1e-2"


@dataclass(frozen=True)
class PacketView:
    packet: object
    topology: object
    stalk: object
    fusion_features: object


@dataclass(frozen=True)
class Stage:
    name: str
    key_name: str
    magnitude: int
    topology: int
    stalk: int
    packet_labels: int
    noncircular: int
    cost: int
    note: str


@dataclass(frozen=True)
class StageStats:
    stage: Stage
    fibers: int
    pair_splits: int
    mixed_route: int
    mixed_status: int
    max_mixed: int
    largest_key: object
    largest_routes: tuple[tuple[str, int], ...]


STAGES = (
    Stage("raw_residual_pair", "pair_id", 0, 0, 0, 0, 6, 1, "negative control: pair identity only"),
    Stage("automatic_word", "automatic_word", 0, 0, 0, 0, 6, 1, "finite automaton shadow"),
    Stage("word_plus_M_q", "word_M_q", 5, 0, 0, 0, 6, 3, "exact scale plus q threshold"),
    Stage("word_plus_boundary_topology", "word_boundary_topology", 2, 4, 0, 0, 5, 4, "coarse safe topology"),
    Stage("closed_arc_topology", "closed_arc_topology", 1, 5, 0, 0, 5, 4, "HYP-3030 closed/open Cech data"),
    Stage("safe_component_owner_stalk", "stalk_owner", 1, 3, 4, 0, 5, 4, "endpoint/peak owner stalk"),
    Stage("safe_component_exact_stalk", "stalk_exact", 2, 4, 5, 0, 5, 5, "exact largest-component stalk"),
    Stage("fusion_signature", "fusion_signature", 4, 4, 3, 3, 5, 6, "HYP-3026 fused packet sidecar"),
    Stage("packet_label_sink", "packet_labels", 2, 2, 1, 6, 5, 5, "C27/K33/transfer packet labels"),
    Stage("route_label_sink", "route", 1, 1, 1, 1, 0, 1, "circular target label; sanity check only"),
)


def topology_key(record) -> tuple[object, ...]:
    return (
        record.closed_arc_betti.beta0,
        record.closed_arc_betti.beta1,
        record.arc_betti.beta0,
        record.arc_betti.beta1,
        record.audit.safe_topes,
        record.runner_quotient_defect,
        record.boundary_pair_sums,
    )


def packet_label_key(p) -> tuple[object, ...]:
    return (
        p.packet_route,
        p.packet_rank,
        p.state_lift,
        p.transfer,
        tuple(p.unknown_pairs),
    )


def stage_key(view: PacketView, stage: Stage) -> object:
    p = view.packet
    if stage.key_name == "pair_id":
        return PAIR_BY_NAME[p.name]
    if stage.key_name == "automatic_word":
        return p.automatic_word
    if stage.key_name == "word_M_q":
        return p.automatic_word, p.M, p.q_threshold
    if stage.key_name == "word_boundary_topology":
        return (
            p.automatic_word,
            p.q_threshold,
            p.strict_components,
            p.boundary_count,
            mu_bucket(p.strict_safe_mu),
        )
    if stage.key_name == "closed_arc_topology":
        return p.automatic_word, topology_key(view.topology)
    if stage.key_name == "stalk_owner":
        return p.automatic_word, view.stalk.owner_only_key
    if stage.key_name == "stalk_exact":
        return p.automatic_word, view.stalk.exact_key
    if stage.key_name == "fusion_signature":
        return view.fusion_features.fusion_signature
    if stage.key_name == "packet_labels":
        return p.automatic_word, packet_label_key(p)
    if stage.key_name == "route":
        return p.automatic_word, p.route
    raise ValueError(stage.key_name)


def status(packet) -> str:
    return "boundary" if packet.strict_safe_mu == 0 else "open"


def target_rows(single_limit: int, two_swap_limit: int, alias_depth: int, lcm_tail_max: int):
    rows = lp.build_bank(single_limit, two_swap_limit, alias_depth, lcm_tail_max)
    selected = [row for row in rows if row[0] in TARGET_NAMES]
    found = {row[0] for row in selected}
    missing = sorted(TARGET_NAMES - found)
    if missing:
        raise RuntimeError(f"missing target rows: {missing}")
    return selected


def build_views(rows, workers: int) -> list[PacketView]:
    packets = lp.compute_packets(rows, workers)
    by_name = {p.name: p for p in packets}
    ordered = [by_name[name] for _pid, _note, names in TARGET_PAIRS for name in names]
    views: list[PacketView] = []
    for packet in ordered:
        topo = ac.row_record(packet.name, tuple(packet.speeds))
        stalk = stalk_mod.largest_component_stalk(packet)
        row = fusion.Row(packet.name, tuple(packet.speeds), packet.route)
        fused = fusion.row_features(row)
        views.append(PacketView(packet, topo, stalk, fused))
    return views


def stats_for_stage(views: list[PacketView], stage: Stage) -> StageStats:
    groups: dict[object, list[PacketView]] = defaultdict(list)
    for view in views:
        groups[stage_key(view, stage)].append(view)

    mixed_route = [rows for rows in groups.values() if len({v.packet.route for v in rows}) > 1]
    mixed_status = [rows for rows in groups.values() if len({status(v.packet) for v in rows}) > 1]
    pair_splits = 0
    for _pid, _note, names in TARGET_PAIRS:
        pair_views = [v for v in views if v.packet.name in names]
        if len({stage_key(v, stage) for v in pair_views}) > 1:
            pair_splits += 1
    if mixed_route:
        largest = max(mixed_route, key=lambda rows: (len(rows), sorted(v.packet.name for v in rows)))
        largest_key = stage_key(largest[0], stage)
        largest_routes = tuple(sorted(Counter(v.packet.route for v in largest).items()))
        max_mixed = len(largest)
    else:
        largest_key = None
        largest_routes = ()
        max_mixed = 0
    return StageStats(
        stage=stage,
        fibers=len(groups),
        pair_splits=pair_splits,
        mixed_route=len(mixed_route),
        mixed_status=len(mixed_status),
        max_mixed=max_mixed,
        largest_key=largest_key,
        largest_routes=largest_routes,
    )


def first_split_stage(pair_views: list[PacketView]) -> str:
    for stage in STAGES:
        if len({stage_key(v, stage) for v in pair_views}) > 1:
            return stage.name
    return "none"


def zeta_exit_class(pair_views: list[PacketView], first_cut: str) -> str:
    routes = {v.packet.route for v in pair_views}
    if "K33-STATE-LIFT" in routes:
        return "cross_handoff"
    if "BOUNDARY-PETAL-SPORADIC" in routes and "COVERING-MOMENT" in routes:
        return "nested_refinement"
    if first_cut in {"closed_arc_topology", "safe_component_owner_stalk", "safe_component_exact_stalk"}:
        return "owner_strip"
    if "BOUNDARY-AP-GW" in routes:
        return "same_tile_boundary"
    return "residual"


def key_state(pair_views: list[PacketView], stage: Stage) -> str:
    keys = [stage_key(v, stage) for v in pair_views]
    return "split" if len(set(map(repr, keys))) > 1 else "same"


def preview(obj: object, limit: int = 92) -> str:
    text = repr(obj)
    return text if len(text) <= limit else text[: limit - 3] + "..."


def carrier_vector(stats: StageStats) -> tuple[int, ...]:
    stage = stats.stage
    return (
        stats.pair_splits,
        10 - stats.mixed_route,
        10 - stats.mixed_status,
        stage.packet_labels,
        stage.topology,
        stage.stalk,
        stage.magnitude,
        stage.noncircular,
        10 - stage.cost,
        stats.fibers,
    )


def tournament_edges(stats: list[StageStats]) -> set[tuple[str, str]]:
    edges: set[tuple[str, str]] = set()
    order = {st.stage.name: i for i, st in enumerate(stats)}
    vectors = {st.stage.name: carrier_vector(st) for st in stats}
    names = [st.stage.name for st in stats]
    for a, b in combinations(names, 2):
        va, vb = vectors[a], vectors[b]
        if va > vb or (va == vb and order[a] < order[b]):
            edges.add((a, b))
        else:
            edges.add((b, a))
    return edges


def score_histogram(names: list[str], edges: set[tuple[str, str]]) -> dict[int, int]:
    wins = Counter()
    for a, _b in edges:
        wins[a] += 1
    return dict(sorted(Counter(wins[name] for name in names).items()))


def directed_3cycles(names: list[str], edges: set[tuple[str, str]]) -> int:
    count = 0
    edge = set(edges)
    for a, b, c in combinations(names, 3):
        if (
            ((a, b) in edge and (b, c) in edge and (c, a) in edge)
            or ((a, c) in edge and (c, b) in edge and (b, a) in edge)
        ):
            count += 1
    return count


def scc_sizes(names: list[str], edges: set[tuple[str, str]]) -> tuple[int, ...]:
    adjacency = {name: [b for a, b in edges if a == name] for name in names}
    reverse = {name: [a for a, b in edges if b == name] for name in names}

    def reach(start: str, graph: dict[str, list[str]]) -> set[str]:
        seen = {start}
        stack = [start]
        while stack:
            cur = stack.pop()
            for nxt in graph[cur]:
                if nxt not in seen:
                    seen.add(nxt)
                    stack.append(nxt)
        return seen

    unseen = set(names)
    sizes: list[int] = []
    while unseen:
        node = next(iter(unseen))
        comp = reach(node, adjacency) & reach(node, reverse)
        sizes.append(len(comp))
        unseen -= comp
    return tuple(sorted(sizes, reverse=True))


def hamiltonian_path_count(names: list[str], edges: set[tuple[str, str]]) -> int:
    index = {name: i for i, name in enumerate(names)}
    n = len(names)
    adj = [[False] * n for _ in range(n)]
    for a, b in edges:
        adj[index[a]][index[b]] = True
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1 << n):
        for last in range(n):
            val = dp[mask][last]
            if not val:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += val
    return sum(dp[(1 << n) - 1])


def print_assumption_challenge() -> None:
    print("[0] Assumption challenge")
    print("  considered vertices:")
    print("    runners, gaps, endpoint walls, wall-crossing events, residue words,")
    print("    cover arcs, safe-component bars, dyadic Haar tiles, packet rows,")
    print("    route labels, residual pairs, and proof obligations.")
    print("  chosen vertices:")
    print("    residual-cut carriers that may discharge a route capacitor.")
    print("  preserved LRC predicate:")
    print("    all four target rows are open; the tested predicate is route-purity")
    print("    after status has already been protected by HYP-3024/HYP-3030 gates.")
    print("  destroyed information:")
    print("    raw runner identity and exact route labels are destroyed until a sidecar")
    print("    separates the capacitor or declares named residual debt.")
    print()


def print_target_packets(views: list[PacketView]) -> None:
    print("[1] Residual capacitors")
    for pair_id, note, names in TARGET_PAIRS:
        pair_views = [v for v in views if v.packet.name in names]
        routes = dict(Counter(v.packet.route for v in pair_views))
        print(f"  {pair_id}: {note}; routes={routes}")
        for view in pair_views:
            p = view.packet
            topo = view.topology
            print(
                f"    {p.name:<34} route={p.route:<24} status={status(p):<4} "
                f"M={fmt(p.M):>6} q={p.q_threshold:<2d} mu={fmt(p.strict_safe_mu):>10} "
                f"components={p.strict_components:<2d} boundary={p.boundary_count:<2d} "
                f"closed=({topo.closed_arc_betti.beta0},{topo.closed_arc_betti.beta1}) "
                f"open=({topo.arc_betti.beta0},{topo.arc_betti.beta1}) "
                f"stalk={view.stalk.short()}"
            )
    print()


def print_stage_table(stats: list[StageStats]) -> None:
    print("[2] Cut-stage table on four target packets")
    print(
        "  stage                         fibers pair_splits mixed_route mixed_status "
        "max_mixed note"
    )
    for st in stats:
        print(
            f"  {st.stage.name:<30} {st.fibers:6d} {st.pair_splits:11d} "
            f"{st.mixed_route:11d} {st.mixed_status:12d} {st.max_mixed:9d} "
            f"{st.stage.note}"
        )
    print()


def print_pair_cuts(views: list[PacketView]) -> None:
    print("[3] First-cut readout")
    for pair_id, note, names in TARGET_PAIRS:
        pair_views = [v for v in views if v.packet.name in names]
        first = first_split_stage(pair_views)
        exit_class = zeta_exit_class(pair_views, first)
        print(f"  {pair_id}")
        print(f"    note={note}")
        print(f"    first_cut_stage={first}")
        print(f"    zeta_exit_class={exit_class}")
        print("    stage states:")
        for stage in STAGES:
            print(f"      {stage.name:<30} {key_state(pair_views, stage)}")
        print("    separating key previews:")
        for stage in STAGES:
            keys = [stage_key(view, stage) for view in pair_views]
            if len(set(map(repr, keys))) > 1:
                for view, key in zip(pair_views, keys):
                    print(f"      {stage.name}:{view.packet.name}: {preview(key)}")
                break
    print()


def print_tournament(stats: list[StageStats]) -> None:
    names = [st.stage.name for st in stats]
    edges = tournament_edges(stats)
    ranking = sorted(stats, key=lambda st: (carrier_vector(st), -STAGES.index(st.stage)), reverse=True)
    print("[4] Tournament Analysis")
    print("  vertices: residual-cut carriers, not runners or arcs.")
    print("  pairwise observable:")
    print("    capacitor separation count, route/status purity, retained packet labels,")
    print("    topology, stalk, magnitude, non-circularity, proof cost, fiber count.")
    print("  switch/gauge:")
    print("    A -> B iff A has lexicographically larger observable vector.")
    print("  tie Hamiltonian path:")
    print("    " + " -> ".join(names))
    print(
        f"  fingerprint: score_hist={score_histogram(names, edges)} "
        f"c3={directed_3cycles(names, edges)} "
        f"scc={scc_sizes(names, edges)} hp={hamiltonian_path_count(names, edges)}"
    )
    print("  high-retention path:")
    print("    " + " > ".join(st.stage.name for st in ranking))
    print()


def print_theorem_target() -> None:
    print("[5] Flow/cut theorem target")
    print("  Residual capacitor lemma:")
    print("    In any automatic/residue fiber after status is protected, a mixed open")
    print("    route pair is a two-plate capacitor.  The first nonzero cut among exact")
    print("    scale, coarse topology, closed-arc topology, safe-component stalk,")
    print("    fusion sidecar, and packet labels must either classify the pair as")
    print("    owner_strip, nested_refinement, cross_handoff, or named residual debt.")
    print("  Current four-packet audit:")
    print("    M_q_petal_covering_capacitor exits by nested_refinement.")
    print("    boundary_topology_k33_covering_capacitor exits by cross_handoff.")
    print("  New angle:")
    print("    Treat the proof ladder as max-flow/min-cut on residual debt.  Status")
    print("    gates remove boundary/open capacity; route capacitors are then finite")
    print("    cut obligations, not new counterexample families.")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--single-limit", type=int, default=180)
    parser.add_argument("--two-swap-limit", type=int, default=36)
    parser.add_argument("--alias-depth", type=int, default=4)
    parser.add_argument("--lcm-tail-max", type=int, default=5)
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="exact-M workers; default 1 because this script dynamically loads helper modules",
    )
    args = parser.parse_args()

    rows = target_rows(args.single_limit, args.two_swap_limit, args.alias_depth, args.lcm_tail_max)
    views = build_views(rows, args.workers)
    stats = [stats_for_stage(views, stage) for stage in STAGES]

    print("LRC14 residual capacitor flow scout (codex-2026-06-26-S196b)")
    print("=" * 78)
    print(
        f"parameters: single_limit={args.single_limit} two_swap_limit={args.two_swap_limit} "
        f"alias_depth={args.alias_depth} lcm_tail_max={args.lcm_tail_max} workers={args.workers}"
    )
    print(f"target_packets={len(views)} target_pairs={len(TARGET_PAIRS)}")
    print()
    print_assumption_challenge()
    print_target_packets(views)
    print_stage_table(stats)
    print_pair_cuts(views)
    print_tournament(stats)
    print_theorem_target()


if __name__ == "__main__":
    main()
