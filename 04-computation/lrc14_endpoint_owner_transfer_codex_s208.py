#!/usr/bin/env python3
"""Endpoint-owner transfer scout for LRC14 residual packets.

S201 showed that the q=23 petal/covering diagonal has the same coarse
endpoint-current word, `B18Z6`, but different external endpoint owners.
S196b showed that two residual capacitors are cut by different cheap sidecars:
one by coarse boundary topology and the other by exact M+q.

This pass checks a smaller hidden carrier: keep the owner names inside the
endpoint current.  The carrier is not the full route label and not a runner
tournament.  It is the external boundary-owner strip, optionally joined with
the largest-safe-component owner stalk.

Tournament Analysis declaration:
  vertices: owner-transfer proof carriers, not runners or raw arcs;
  pairwise observable: q=23 diagonal split, residual-capacitor split count,
    route purity, endpoint detail, stalk detail, topology/exact-scale retention,
    non-circularity, proof cost, and fiber count;
  switch/gauge: lexicographic comparison of those observables;
  tie Hamiltonian path: the CARRIERS order below.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from pathlib import Path
import sys


REPO = Path(__file__).resolve().parents[1]


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    assert spec and spec.loader
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


q23 = load_module(
    "lrc14_q23_square_s202",
    REPO / "04-computation" / "lrc14_q23_drop_add_haar_square_codex_s201.py",
)
cap = load_module(
    "lrc14_residual_capacitor_s202",
    REPO / "04-computation" / "lrc14_residual_capacitor_flow_codex_s196.py",
)


Q_TO_CAP_NAME = {
    "AA_petal_diag": "two drop(10, 13)->add(20, 26)",
    "BB_cover_diag": "two drop(8, 12)->add(16, 24)",
}


def fmt(fr: Fraction) -> str:
    return str(fr.numerator) if fr.denominator == 1 else f"{fr.numerator}/{fr.denominator}"


def external_owner_word(speeds: tuple[int, ...]) -> str:
    return q23.active_words(speeds)[0]


def parse_external_word(word: str) -> Counter[str]:
    if word == "none":
        return Counter()
    out: Counter[str] = Counter()
    for token in word.split(","):
        owner, count = token.split("x")
        out[owner] += int(count)
    return out


def parse_external_residues(word: str) -> Counter[str]:
    out: Counter[str] = Counter()
    for owner, count in parse_external_word(word).items():
        residue = owner.split(":", 1)[0]
        out[residue] += count
    return out


def delta_word(before: str, after: str, residues_only: bool = False) -> str:
    parser = parse_external_residues if residues_only else parse_external_word
    delta = parser(after)
    delta.subtract(parser(before))
    parts = []
    for key, val in sorted(delta.items(), key=lambda item: (item[0], item[1])):
        if val > 0:
            parts.append(f"+{key}x{val}")
        elif val < 0:
            parts.append(f"-{key}x{-val}")
    return ",".join(parts) if parts else "zero"


def stage(name: str):
    for item in cap.STAGES:
        if item.name == name:
            return item
    raise KeyError(name)


WORD_PLUS_BOUNDARY_TOPOLOGY = stage("word_plus_boundary_topology")


def stalk_word(view) -> str:
    stalk = view.stalk
    if stalk.status == "boundary":
        return f"boundary|B={len(stalk.boundary_owners)}|Z={stalk.boundary_zero_pairs}"
    return (
        f"{stalk.left_owners}->{stalk.peak_owners}->{stalk.right_owners}"
        f"|h={fmt(stalk.peak_height)}|len={fmt(stalk.length)}"
    )


@dataclass(frozen=True)
class QRecord:
    key: str
    route: str
    features: object
    external: str
    cap_view: object


@dataclass(frozen=True)
class Carrier:
    name: str
    q_key: object
    cap_key: object
    endpoint_detail: int
    stalk_detail: int
    topology_detail: int
    exact_detail: int
    noncircular: int
    cost: int
    destroys: str


def build_context():
    cells = q23.build_cells()
    rows = cap.target_rows(180, 36, 4, 5)
    views = cap.build_views(rows, 1)
    cap_by_name = {view.packet.name: view for view in views}
    q_records = []
    for key in ("AA_petal_diag", "BB_cover_diag"):
        cell = cells[key]
        q_records.append(
            QRecord(
                key=key,
                route=cell.spec.declared_route,
                features=cell.features,
                external=cell.external_owner_word,
                cap_view=cap_by_name[Q_TO_CAP_NAME[key]],
            )
        )
    return cells, views, q_records


def carrier_bank() -> tuple[Carrier, ...]:
    return (
        Carrier(
            "raw_residual_shadow",
            lambda q: "q23_diagonal",
            lambda view: cap.PAIR_BY_NAME[view.packet.name],
            0,
            0,
            0,
            0,
            7,
            1,
            "owner names, exact scale, topology, stalks, and routes",
        ),
        Carrier(
            "coarse_endpoint_count",
            lambda q: q.features.endpoint_current_word,
            lambda view: view.fusion_features.endpoint_current_word,
            1,
            0,
            0,
            0,
            7,
            1,
            "which endpoint owners realize the B/Z count",
        ),
        Carrier(
            "exact_M_q",
            lambda q: q.features.exact_m,
            lambda view: (view.packet.M, view.packet.q_threshold),
            0,
            0,
            0,
            4,
            7,
            3,
            "owner strip inside equal exact scale fibers",
        ),
        Carrier(
            "coarse_safe_body",
            lambda q: (
                q.features.bar_count,
                q.features.longest_bar,
                q.features.boundary_count,
                q.features.zero_sum_pairs,
            ),
            lambda view: cap.stage_key(view, WORD_PLUS_BOUNDARY_TOPOLOGY),
            1,
            0,
            4,
            1,
            6,
            4,
            "which owner caused the safe-body/topology split",
        ),
        Carrier(
            "safe_component_owner_stalk",
            lambda q: q.cap_view.stalk.owner_only_key,
            lambda view: view.stalk.owner_only_key,
            3,
            4,
            2,
            0,
            6,
            4,
            "external boundary speed token and nonlargest boundary facets",
        ),
        Carrier(
            "external_endpoint_owner_strip",
            lambda q: q.external,
            lambda view: external_owner_word(tuple(view.packet.speeds)),
            5,
            1,
            2,
            0,
            7,
            3,
            "safe-component interior stalk height/length",
        ),
        Carrier(
            "owner_transfer_carrier",
            lambda q: (q.external, q.cap_view.stalk.owner_only_key),
            lambda view: (external_owner_word(tuple(view.packet.speeds)), view.stalk.owner_only_key),
            5,
            4,
            2,
            0,
            7,
            4,
            "exact length/body data and full packet labels",
        ),
        Carrier(
            "route_label_sink",
            lambda q: q.route,
            lambda view: view.packet.route,
            1,
            1,
            1,
            1,
            0,
            1,
            "non-circular proof payload before route labels are known",
        ),
    )


def q_split(carrier: Carrier, q_records: list[QRecord]) -> bool:
    return len({repr(carrier.q_key(record)) for record in q_records}) > 1


def capacitor_splits(carrier: Carrier, views: list[object]) -> int:
    total = 0
    for _pair_id, _note, names in cap.TARGET_PAIRS:
        pair_views = [view for view in views if view.packet.name in names]
        if len({repr(carrier.cap_key(view)) for view in pair_views}) > 1:
            total += 1
    return total


def mixed_route_groups(carrier: Carrier, views: list[object]) -> int:
    groups: dict[str, list[object]] = defaultdict(list)
    for view in views:
        groups[repr(carrier.cap_key(view))].append(view)
    return sum(1 for group in groups.values() if len({view.packet.route for view in group}) > 1)


def carrier_vector(carrier: Carrier, views: list[object], q_records: list[QRecord]) -> tuple[int, ...]:
    qdiag = int(q_split(carrier, q_records))
    caps = capacitor_splits(carrier, views)
    groups: dict[str, list[object]] = defaultdict(list)
    for view in views:
        groups[repr(carrier.cap_key(view))].append(view)
    return (
        qdiag + caps,
        qdiag,
        caps,
        10 - mixed_route_groups(carrier, views),
        carrier.endpoint_detail,
        carrier.stalk_detail,
        carrier.topology_detail,
        carrier.exact_detail,
        carrier.noncircular,
        10 - carrier.cost,
        len(groups),
    )


def tournament_edges(carriers: tuple[Carrier, ...], views: list[object], q_records: list[QRecord]):
    order = {carrier.name: idx for idx, carrier in enumerate(carriers)}
    vectors = {carrier.name: carrier_vector(carrier, views, q_records) for carrier in carriers}
    edges = set()
    for a, b in combinations(carriers, 2):
        av, bv = vectors[a.name], vectors[b.name]
        if av > bv or (av == bv and order[a.name] < order[b.name]):
            edges.add((a.name, b.name))
        else:
            edges.add((b.name, a.name))
    return edges


def score_histogram(names: list[str], edges: set[tuple[str, str]]) -> dict[int, int]:
    wins = Counter()
    for a, _b in edges:
        wins[a] += 1
    return dict(sorted(Counter(wins[name] for name in names).items()))


def directed_3cycles(names: list[str], edges: set[tuple[str, str]]) -> int:
    edge = set(edges)
    count = 0
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
    print("    runners, gaps, fixed circle sections, section boundaries, wall-crossing")
    print("    events, residues, cover arcs, Fourier modes, Haar tiles, safe-component")
    print("    bars, endpoint owner tokens, matroid circuits, and proof obligations.")
    print("  chosen vertices:")
    print("    owner-transfer carriers: external boundary-owner strips and")
    print("    largest-safe-component owner stalks.")
    print("  preserved LRC predicate:")
    print("    status is already open on the tested rows; the preserved predicate is")
    print("    residual route schedulability without using route labels as a key.")
    print("  destroyed information:")
    print("    exact nonlargest safe components, full packet route labels, and global")
    print("    family proofs are destroyed until reattached downstream.")
    print()


def print_q23_diagonal(q_records: list[QRecord]) -> None:
    print("[1] q=23 diagonal owner transfer")
    print("  row                         M      endpoint ext_owner_strip          stalk_owner_word")
    for record in q_records:
        print(
            f"  {record.key:<27} {fmt(record.features.exact_m):<6} "
            f"{record.features.endpoint_current_word:<8} {record.external:<25} "
            f"{stalk_word(record.cap_view)}"
        )
    before, after = q_records
    print("  diagonal hidden split:")
    print(
        f"    coarse endpoint word stays {before.features.endpoint_current_word}; "
        f"token_delta({before.key}->{after.key})={delta_word(before.external, after.external)}"
    )
    print(
        f"    residue_delta({before.key}->{after.key})="
        f"{delta_word(before.external, after.external, residues_only=True)}"
    )
    print()


def print_capacitor_owner_table(views: list[object]) -> None:
    print("[2] Residual capacitor owner strips")
    for pair_id, note, names in cap.TARGET_PAIRS:
        pair_views = [view for view in views if view.packet.name in names]
        exact_state = "split" if len({(view.packet.M, view.packet.q_threshold) for view in pair_views}) > 1 else "same"
        topology_state = "split" if len({repr(cap.stage_key(view, WORD_PLUS_BOUNDARY_TOPOLOGY)) for view in pair_views}) > 1 else "same"
        external_state = "split" if len({external_owner_word(tuple(view.packet.speeds)) for view in pair_views}) > 1 else "same"
        stalk_state = "split" if len({repr(view.stalk.owner_only_key) for view in pair_views}) > 1 else "same"
        print(f"  {pair_id}: {note}")
        print(
            f"    exact_M_q={exact_state} boundary_topology={topology_state} "
            f"external_owner_strip={external_state} owner_stalk={stalk_state}"
        )
        for view in pair_views:
            p = view.packet
            print(
                f"    {p.name:<34} route={p.route:<24} M={fmt(p.M):<5} "
                f"endpoint={view.fusion_features.endpoint_current_word:<6} "
                f"ext={external_owner_word(tuple(p.speeds)):<25} "
                f"stalk={stalk_word(view)}"
            )
    print()


def print_transfer_fibers(views: list[object]) -> None:
    print("[3] B18Z6 fiber refinement")
    groups: dict[str, list[object]] = defaultdict(list)
    for view in views:
        groups[view.fusion_features.endpoint_current_word].append(view)
    for endpoint, group in sorted(groups.items()):
        routes = sorted(Counter(view.packet.route for view in group).items())
        print(f"  coarse_endpoint={endpoint} rows={len(group)} routes={routes}")
        strips = defaultdict(list)
        for view in group:
            strips[external_owner_word(tuple(view.packet.speeds))].append(view)
        for strip, strip_views in sorted(strips.items()):
            labels = ",".join(view.packet.name for view in strip_views)
            strip_routes = sorted({view.packet.route for view in strip_views})
            print(f"    ext={strip:<25} routes={strip_routes} rows={labels}")
    print("  hidden_statement=the owner names inside B18Z6, not the B/Z counts alone,")
    print("    refine the q=23 diagonal and both displayed residual capacitors.")
    print()


def print_tournament(views: list[object], q_records: list[QRecord]) -> None:
    carriers = carrier_bank()
    names = [carrier.name for carrier in carriers]
    edges = tournament_edges(carriers, views, q_records)
    ranking = sorted(carriers, key=lambda carrier: carrier_vector(carrier, views, q_records), reverse=True)
    print("[4] Tournament Analysis")
    print("  vertices: owner-transfer proof carriers, not runners or raw arcs.")
    print("  pairwise observable:")
    print("    q=23 diagonal split, capacitor split count, route purity, endpoint")
    print("    detail, stalk detail, topology/exact-scale retention, non-circularity,")
    print("    proof cost, and fiber count.")
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
    print("    " + " > ".join(carrier.name for carrier in ranking))
    print("  carrier table:")
    print("    carrier                         qdiag cap_splits mixed_route fibers destroys")
    for carrier in ranking:
        groups = {repr(carrier.cap_key(view)) for view in views}
        print(
            f"    {carrier.name:<31} {int(q_split(carrier, q_records)):5d} "
            f"{capacitor_splits(carrier, views):10d} "
            f"{mixed_route_groups(carrier, views):11d} {len(groups):6d} "
            f"{carrier.destroys}"
        )
    print()


def print_theorem_target() -> None:
    print("[5] Candidate proof statement")
    print("  Endpoint-owner transfer lemma target:")
    print("    Inside the protected strict-open B18Z6 residual surface, the external")
    print("    endpoint-owner strip is a non-route local address.  It refines exact-M")
    print("    and coarse-topology capacitor cuts, splits the q=23 diagonal, and can")
    print("    be joined with the largest-safe-component owner stalk before falling")
    print("    back to full packet labels.")
    print("  Current finite evidence:")
    print("    all four residual capacitor packets have endpoint_current_word=B18Z6;")
    print("    external owner strips split both residual pairs and distinguish petal,")
    print("    K33, q=23 covering, and single-swap covering packets.")
    print("  Guardrail:")
    print("    coarse B/Z endpoint counts are diagnostics only.  A proof quotient must")
    print("    retain owner names, reconstruct them, annihilate them by a dual")
    print("    certificate, or route the lost coordinate to named residual debt.")


def main() -> None:
    _cells, views, q_records = build_context()
    print("LRC14 endpoint-owner transfer scout (codex-2026-06-26-S208)")
    print("source_threads=HYP-3044,HYP-3042,HYP-3041,HYP-3040,HYP-3039,HYP-3038,HYP-3037,HYP-3036,HYP-3035,HYP-3032,HYP-3031,HYP-3027,HYP-3026,HYP-3018")
    print()
    print_assumption_challenge()
    print_q23_diagonal(q_records)
    print_capacitor_owner_table(views)
    print_transfer_fibers(views)
    print_tournament(views, q_records)
    print_theorem_target()


if __name__ == "__main__":
    main()
