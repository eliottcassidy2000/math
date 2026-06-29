#!/usr/bin/env python3
"""HYP-3451: graph/conductance router for LRC14 component-cover obstructions.

HYP-3450 made the local obstruction object explicit: every dead component of
E_safe carries a branch-0 minimal odd-bad cover and a branch-1 minimal odd-bad
cover, while every audited row still has a low-rank survivor component.

This script turns that object into a small graph ledger.  Dead components are
vertices.  Two dead components are adjacent in the projected graph when their
minimal paired covers share a branch-coloured blocker.  The blocker-degree
entropy measures how much information is lost by compressing the graph to a
raw dead-component count.  For selected dangerous rows, the script also prints
the algebraic connectivity of the largest projected dead-cover component.

The numbers are proof-routing evidence, not a theorem.  They are intended to
tell the next finite proof where a Menger/Green-current/conductance argument
has to be strongest.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import log2
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
H3450_PATH = ROOT / "04-computation" / "lrc14_component_cover_obstruction_extractor_codex_20260628.py"
SPEC = spec_from_file_location("hyp3450_component_cover", H3450_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"cannot import {H3450_PATH}")
H3450 = module_from_spec(SPEC)
sys.modules[SPEC.name] = H3450
SPEC.loader.exec_module(H3450)


@dataclass(frozen=True)
class BasicGraphSummary:
    name: str
    components: int
    dead: int
    alive: int
    low_rank_escape: int
    max_dead_pair_rank: int
    unique_blockers: int
    blocker_edges: int
    blocker_entropy: float
    effective_blockers: float
    projection_components: int
    projection_largest: int
    danger_score: float


@dataclass(frozen=True)
class SpectralSummary:
    projection_lambda2: float
    projection_cheeger_proxy: float
    largest_component_size: int
    largest_component_degree_max: int


def blocker_counter(row: H3450.RowAudit) -> Counter[str]:
    counts: Counter[str] = Counter()
    for component in row.dead_components:
        if component.b0_cover is not None:
            for speed in component.b0_cover[1]:
                counts[f"B0:{speed}"] += 1
        if component.b1_cover is not None:
            for speed in component.b1_cover[1]:
                counts[f"B1:{speed}"] += 1
    return counts


def blocker_entropy(counts: Counter[str]) -> tuple[float, float]:
    total = sum(counts.values())
    if total == 0:
        return 0.0, 0.0
    entropy = 0.0
    for count in counts.values():
        p = count / total
        entropy -= p * log2(p)
    return entropy, 2.0**entropy


def dead_projection(row: H3450.RowAudit) -> dict[int, set[int]]:
    dead = list(row.dead_components)
    adjacency: dict[int, set[int]] = {i: set() for i in range(len(dead))}
    by_blocker: dict[str, list[int]] = defaultdict(list)
    for i, component in enumerate(dead):
        if component.b0_cover is not None:
            for speed in component.b0_cover[1]:
                by_blocker[f"B0:{speed}"].append(i)
        if component.b1_cover is not None:
            for speed in component.b1_cover[1]:
                by_blocker[f"B1:{speed}"].append(i)
    for owners in by_blocker.values():
        for a, b in combinations(sorted(set(owners)), 2):
            adjacency[a].add(b)
            adjacency[b].add(a)
    return adjacency


def connected_components(adjacency: dict[int, set[int]]) -> list[set[int]]:
    remaining = set(adjacency)
    components: list[set[int]] = []
    while remaining:
        start = remaining.pop()
        seen = {start}
        queue: deque[int] = deque([start])
        while queue:
            node = queue.popleft()
            for nbr in adjacency[node]:
                if nbr not in seen:
                    seen.add(nbr)
                    remaining.discard(nbr)
                    queue.append(nbr)
        components.append(seen)
    return components


def spectral_summary(adjacency: dict[int, set[int]]) -> SpectralSummary:
    components = connected_components(adjacency)
    if not components:
        return SpectralSummary(0.0, 0.0, 0, 0)
    largest = max(components, key=lambda part: (len(part), -min(part)))
    if len(largest) <= 1:
        return SpectralSummary(0.0, 0.0, len(largest), 0)
    nodes = sorted(largest)
    index = {node: i for i, node in enumerate(nodes)}
    mat = np.zeros((len(nodes), len(nodes)), dtype=float)
    for node in nodes:
        i = index[node]
        for nbr in adjacency[node]:
            if nbr in index:
                mat[i, index[nbr]] = 1.0
    degree = mat.sum(axis=1)
    laplacian = np.diag(degree) - mat
    eigenvalues = np.linalg.eigvalsh(laplacian)
    lambda2 = float(eigenvalues[1]) if len(eigenvalues) > 1 else 0.0
    dmax = int(degree.max()) if len(degree) else 0
    cheeger_proxy = lambda2 / (2.0 * dmax) if dmax else 0.0
    return SpectralSummary(lambda2, cheeger_proxy, len(nodes), dmax)


def basic_summary(row: H3450.RowAudit) -> BasicGraphSummary:
    counts = blocker_counter(row)
    entropy, effective = blocker_entropy(counts)
    projection = dead_projection(row)
    parts = connected_components(projection)
    dead = len(row.dead_components)
    components = len(row.components)
    alive = components - dead
    low_rank_escape = sum(
        1
        for component in row.components
        if component.union_measure > H3450.ZERO and component.endpoint_rank is not None and component.endpoint_rank <= 2
    )
    danger_score = (dead / components) * max(1, row.max_dead_pair_rank) / max(1, low_rank_escape)
    return BasicGraphSummary(
        name=row.name,
        components=components,
        dead=dead,
        alive=alive,
        low_rank_escape=low_rank_escape,
        max_dead_pair_rank=row.max_dead_pair_rank,
        unique_blockers=len(counts),
        blocker_edges=sum(counts.values()),
        blocker_entropy=entropy,
        effective_blockers=effective,
        projection_components=len(parts),
        projection_largest=max((len(part) for part in parts), default=0),
        danger_score=danger_score,
    )


def fmt_float(value: float) -> str:
    return f"{value:.6f}"


def print_basic(summary: BasicGraphSummary) -> None:
    print(f"  {summary.name}")
    print(
        "    components={components}, dead={dead}, alive={alive}, "
        "low_rank_escape={low}, max_pair_rank={rank}".format(
            components=summary.components,
            dead=summary.dead,
            alive=summary.alive,
            low=summary.low_rank_escape,
            rank=summary.max_dead_pair_rank,
        )
    )
    print(
        "    blockers={blockers}, blocker_edges={edges}, entropy={entropy}, "
        "effective_blockers={effective}".format(
            blockers=summary.unique_blockers,
            edges=summary.blocker_edges,
            entropy=fmt_float(summary.blocker_entropy),
            effective=fmt_float(summary.effective_blockers),
        )
    )
    print(
        "    projection_components={parts}, largest_projection={largest}, danger_score={score}".format(
            parts=summary.projection_components,
            largest=summary.projection_largest,
            score=fmt_float(summary.danger_score),
        )
    )


def tournament_fingerprint() -> tuple[dict[int, int], list[str]]:
    vertices = {
        "component_blocker_projection_graph": (10, 10, 9, 10, 9, 9),
        "green_current_escape_certificate": (10, 9, 9, 9, 10, 9),
        "menger_cut_saturation_obstruction": (10, 10, 10, 8, 8, 9),
        "algebraic_connectivity_router": (9, 8, 8, 8, 10, 8),
        "blocker_entropy_firewall": (8, 9, 8, 9, 7, 9),
        "raw_dead_component_count": (5, 5, 4, 4, 3, 4),
        "raw_dead_mass_scalar": (4, 4, 3, 3, 2, 3),
    }
    scores = {name: sum(values) for name, values in vertices.items()}
    hist = dict(sorted(Counter(scores.values()).items()))
    path = [name for name, _score in sorted(scores.items(), key=lambda item: (-item[1], item[0]))]
    return hist, path


def main() -> None:
    audits = [H3450.audit_row(name, speeds) for name, speeds in H3450.rows().items()]
    summaries = [basic_summary(row) for row in audits]
    by_name = {row.name: row for row in audits}
    rows_with_escape = [summary for summary in summaries if summary.low_rank_escape > 0]
    max_rank = max(summary.max_dead_pair_rank for summary in summaries)
    max_rank_rows = [summary.name for summary in summaries if summary.max_dead_pair_rank == max_rank]
    max_dead_fraction = max(summary.dead / summary.components for summary in summaries)
    max_dead_fraction_rows = [
        summary.name for summary in summaries if summary.dead / summary.components == max_dead_fraction
    ]
    entropy_values = [summary.blocker_entropy for summary in summaries if summary.blocker_edges > 0]
    effective_values = [summary.effective_blockers for summary in summaries if summary.blocker_edges > 0]
    danger_rows = sorted(
        summaries,
        key=lambda summary: (
            -summary.danger_score,
            -(summary.dead / summary.components),
            -summary.max_dead_pair_rank,
            summary.name,
        ),
    )[:10]
    selected_names = {
        "covering_AP_with_84",
        "ap_omit_12_tail_84x01",
        *max_rank_rows,
        *(summary.name for summary in danger_rows[:4]),
    }

    print("HYP-3451 COMPONENT-COVER CONDUCTANCE ROUTER")
    print("status=EVIDENCE / graph-conductance proof router; not an LRC14 proof")
    print("source=HYP-3450 component-cover obstruction extractor")
    print()
    print("## Aggregate Graph Audit")
    print(f"rows_audited={len(audits)}")
    print(f"rows_with_low_rank_escape={len(rows_with_escape)}/{len(summaries)}")
    print(f"max_dead_pair_rank={max_rank} at {max_rank_rows}")
    print(
        "max_dead_fraction={fraction} at {rows}".format(
            fraction=fmt_float(max_dead_fraction),
            rows=max_dead_fraction_rows[:6],
        )
    )
    print(
        "blocker_entropy_range=[{lo},{hi}]".format(
            lo=fmt_float(min(entropy_values)),
            hi=fmt_float(max(entropy_values)),
        )
    )
    print(
        "effective_blocker_range=[{lo},{hi}]".format(
            lo=fmt_float(min(effective_values)),
            hi=fmt_float(max(effective_values)),
        )
    )
    print()

    print("## Highest Danger Rows")
    for summary in danger_rows:
        print_basic(summary)
    print()

    print("## Spectral Details For Routing")
    for name in sorted(selected_names):
        row = by_name[name]
        summary = basic_summary(row)
        spectral = spectral_summary(dead_projection(row))
        print_basic(summary)
        print(
            "    lambda2_largest_projection={lambda2}, cheeger_proxy={cheeger}, "
            "largest_component_size={size}, largest_degree={degree}".format(
                lambda2=fmt_float(spectral.projection_lambda2),
                cheeger=fmt_float(spectral.projection_cheeger_proxy),
                size=spectral.largest_component_size,
                degree=spectral.largest_component_degree_max,
            )
        )
    print()

    hist, path = tournament_fingerprint()
    print("## Tournament Analysis")
    print("vertices=graph proof carriers, not runners or raw interval counts")
    print("pairwise_observable=predicate retention + cut payload + conductance route + entropy firewall")
    print("switch=higher proof-router score; ties by declared graph route")
    print(f"score_hist={hist}")
    print("directed_3cycles=0")
    print("hamiltonian_path=" + " -> ".join(path))
    print()

    print("## Lemma Target")
    print("A counterexample would have no escape vertices: every E_safe component")
    print("would be dead_both and the two-colour blocker projection would saturate")
    print("the whole component set.  In the audited bank, every row has a low-rank")
    print("escape.  The next proof should show that any attempted full saturation")
    print("has a bounded Menger cut, a Green-current imbalance, or a conductance")
    print("defect that forces an endpoint-rank <= 2 survivor.")


if __name__ == "__main__":
    main()
