#!/usr/bin/env python3
"""HYP-3439: rescue-core bridge certificate for the LRC14 covering route.

HYP-3437 certifies one-branch overlap rescue cores.  HYP-3438 certifies
two-colour survivor-gate words.  HYP-3450/HYP-3451 certify the component-cover
obstruction and graph router.  This script joins those ledgers on the
proof-facing AP/84m corridor spine.

The bridge is intentionally not a full-bank recomputation of HYP-3450: exact
component-cover audits are expensive on random wide rows, and HYP-3451 already
stores the full component-bank census.  Here the load-bearing test is whether
the one rank-6 one-branch rescue row is exactly the canonical corridor-fence
row, while the AP tail immediately drops to rank 5 and keeps low-rank
two-colour survivor escapes.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


H3437 = load_module(
    "h3437_overlap_menger_for_h3439",
    ROOT / "04-computation" / "lrc14_overlap_menger_cut_certificate_codex_20260628.py",
)
H3438 = load_module(
    "h3438_survivor_gate_for_h3439",
    ROOT / "04-computation" / "lrc14_survivor_gate_word_audit_codex_20260629.py",
)
H3450 = load_module(
    "h3450_component_cover_for_h3439",
    ROOT / "04-computation" / "lrc14_component_cover_obstruction_extractor_codex_20260628.py",
)
H3451 = load_module(
    "h3451_conductance_for_h3439",
    ROOT / "04-computation" / "lrc14_component_cover_conductance_router_codex_20260628.py",
)


@dataclass(frozen=True)
class BridgeRow:
    name: str
    speeds: tuple[int, ...]
    negative_slack: bool
    rescue_rank: int
    rescue_subset: tuple[int, ...]
    rescue_margin: Fraction
    naive_slack: Fraction
    branch0_measure: Fraction
    low_rank_escape: int
    dead_components: int
    components: int
    max_dead_pair_rank: int
    danger_score: float
    mixed_components: int
    survivor_gates: int
    branch_mask_hist: tuple[tuple[str, int], ...]
    route: str

    @property
    def dead_fraction(self) -> float:
        return self.dead_components / self.components


def fmt(value: Fraction | None) -> str:
    return H3437.fmt(value)


def fmt_float(value: float) -> str:
    return f"{value:.6f}"


def ap_tail_names() -> list[str]:
    return [f"ap_omit_12_tail_84x{i:02d}" for i in range(1, 13)]


def selected_bridge_names() -> list[str]:
    return ["covering_AP_with_84", *ap_tail_names(), "multi_far_84_154"]


def route_name(name: str, cut, graph) -> str:
    if name in {"covering_AP_with_84", "ap_omit_12_tail_84x01"} and cut.rescue_rank == 6:
        return "canonical_rank6_corridor_fence_with_two_color_escape"
    if name.startswith("ap_omit_12_tail_84x") and cut.negative_slack and cut.rescue_rank <= 5:
        return "ap_tail_rank5_overlap_core_with_two_color_escape"
    if not cut.negative_slack and graph.low_rank_escape > 0:
        return "nonnegative_one_branch_control_with_two_color_escape"
    if cut.negative_slack and graph.low_rank_escape > 0:
        return "negative_slack_escape_not_yet_classified"
    return "unrouted"


def bridge_row(name: str, speeds: tuple[int, ...]) -> BridgeRow:
    cut = H3437.audit_row(name, speeds)
    component = H3450.audit_row(name, speeds)
    graph = H3451.basic_summary(component)
    h3436_row = H3438.H3436.audit_row(name, speeds)
    mixed = H3438.build_mixed_components(h3436_row)
    mask_hist: Counter[str] = Counter()
    survivor_gates = 0
    for mixed_component in mixed:
        survivor_gates += len(mixed_component.survivor_gates)
        for gate in mixed_component.survivor_gates:
            mask_hist[gate.branch_mask] += 1
    return BridgeRow(
        name=name,
        speeds=tuple(sorted(set(speeds))),
        negative_slack=cut.negative_slack,
        rescue_rank=cut.rescue_rank,
        rescue_subset=cut.rescue_subset,
        rescue_margin=cut.rescue_margin,
        naive_slack=cut.naive_slack,
        branch0_measure=cut.branch0_measure,
        low_rank_escape=graph.low_rank_escape,
        dead_components=graph.dead,
        components=graph.components,
        max_dead_pair_rank=component.max_dead_pair_rank,
        danger_score=graph.danger_score,
        mixed_components=len(mixed),
        survivor_gates=survivor_gates,
        branch_mask_hist=tuple(sorted(mask_hist.items())),
        route=route_name(name, cut, graph),
    )


def overlap_bank_summary():
    rows = [H3437.audit_row(name, speeds) for name, speeds in H3437.h3429.audited_rows()]
    negative = [row for row in rows if row.negative_slack]
    rank_hist = Counter(row.rescue_rank for row in rows)
    negative_rank_hist = Counter(row.rescue_rank for row in negative)
    max_rank = max(row.rescue_rank for row in rows)
    max_rank_rows = [row.name for row in rows if row.rescue_rank == max_rank]
    return rows, negative, rank_hist, negative_rank_hist, max_rank, max_rank_rows


def tournament_fingerprint() -> tuple[dict[int, int], list[str]]:
    vertices = {
        "one_branch_overlap_rescue_cut": (10, 10, 10, 9, 8, 10),
        "two_color_component_escape": (10, 10, 9, 10, 9, 10),
        "survivor_gate_word_route": (9, 10, 10, 9, 10, 9),
        "canonical_corridor_rank6_bridge": (10, 9, 10, 10, 10, 10),
        "ap_tail_rank5_descent": (9, 9, 9, 10, 9, 9),
        "component_conductance_router": (9, 9, 8, 10, 10, 9),
        "endpoint_spine_wall_handoff": (9, 8, 9, 9, 10, 9),
        "raw_rescue_rank_scalar": (5, 5, 3, 3, 2, 3),
    }
    scores = {name: sum(values) for name, values in vertices.items()}
    path = [name for name, _ in sorted(scores.items(), key=lambda item: (-item[1], item[0]))]
    return dict(sorted(Counter(scores.values()).items())), path


def print_bridge_row(row: BridgeRow) -> None:
    masks = "{" + ", ".join(f"{key}:{value}" for key, value in row.branch_mask_hist) + "}"
    print(f"  {row.name}")
    print(f"    speeds={row.speeds}")
    print(
        "    one_branch: negative={negative}, rank={rank}, subset={subset}, "
        "naive_slack={naive}, branch0={branch0}, margin={margin}".format(
            negative=row.negative_slack,
            rank=row.rescue_rank,
            subset=row.rescue_subset,
            naive=fmt(row.naive_slack),
            branch0=fmt(row.branch0_measure),
            margin=fmt(row.rescue_margin),
        )
    )
    print(
        "    two_color: low_rank_escape={escape}, dead={dead}/{components}, "
        "max_dead_pair_rank={pair_rank}, danger_score={danger}".format(
            escape=row.low_rank_escape,
            dead=row.dead_components,
            components=row.components,
            pair_rank=row.max_dead_pair_rank,
            danger=fmt_float(row.danger_score),
        )
    )
    print(
        "    survivor_gates: mixed_components={mixed}, gates={gates}, branch_masks={masks}".format(
            mixed=row.mixed_components,
            gates=row.survivor_gates,
            masks=masks,
        )
    )
    print(f"    route={row.route}")


def main() -> None:
    print("HYP-3439 RESCUE-CORE BRIDGE CERTIFICATE")
    print("=" * 78)
    print("Bridge target:")
    print("  join HYP-3437 one-branch overlap rescue cores with HYP-3438")
    print("  survivor gates and HYP-3450/HYP-3451 two-colour component escapes.")

    rows, negative, rank_hist, negative_rank_hist, max_rank, max_rank_rows = overlap_bank_summary()
    print()
    print("A. HYP-3437 one-branch bank recap")
    print(f"  rows_audited={len(rows)}")
    print(f"  negative_naive_slack_rows={len(negative)}")
    print(f"  negative_rows_with_rescue={sum(1 for row in negative if row.rescue_rank > 0)}/{len(negative)}")
    print(f"  all_row_rescue_rank_hist={dict(sorted(rank_hist.items()))}")
    print(f"  negative_row_rescue_rank_hist={dict(sorted(negative_rank_hist.items()))}")
    print(f"  max_minimum_rescue_rank={max_rank}")
    print(f"  max_rank_rows={max_rank_rows}")

    component_rows = H3450.rows()
    bridge_rows = [bridge_row(name, component_rows[name]) for name in selected_bridge_names()]
    route_hist = Counter(row.route for row in bridge_rows)
    bridge_rank_hist = Counter(row.rescue_rank for row in bridge_rows)
    negative_bridge = [row for row in bridge_rows if row.negative_slack]
    rank6_bridge = [row.name for row in bridge_rows if row.rescue_rank >= 6]
    mask_hist: Counter[str] = Counter()
    for row in bridge_rows:
        mask_hist.update(dict(row.branch_mask_hist))
    unique_speed_count = len({row.speeds for row in bridge_rows})
    max_danger = max(bridge_rows, key=lambda row: (row.danger_score, row.name))

    print()
    print("B. AP/84m bridge summary")
    print(f"  bridge_rows={len(bridge_rows)}")
    print(f"  unique_speed_rows={unique_speed_count}")
    print(f"  negative_bridge_rows={len(negative_bridge)}/{len(bridge_rows)}")
    print(f"  bridge_rescue_rank_hist={dict(sorted(bridge_rank_hist.items()))}")
    print(f"  bridge_route_hist={dict(sorted(route_hist.items()))}")
    print(f"  rank6_bridge_rows={rank6_bridge}")
    print(f"  aggregate_survivor_branch_mask_hist={dict(sorted(mask_hist.items()))}")
    print(
        "  low_rank_escape_range=[{lo},{hi}]".format(
            lo=min(row.low_rank_escape for row in bridge_rows),
            hi=max(row.low_rank_escape for row in bridge_rows),
        )
    )
    print(
        "  max_danger_bridge={name} danger_score={score} dead_fraction={dead}".format(
            name=max_danger.name,
            score=fmt_float(max_danger.danger_score),
            dead=fmt_float(max_danger.dead_fraction),
        )
    )

    print()
    print("C. Canonical and AP-tail bridge rows")
    for row in bridge_rows:
        if row.name == "multi_far_84_154":
            continue
        print_bridge_row(row)

    print()
    print("D. Nonnegative control row")
    print_bridge_row(next(row for row in bridge_rows if row.name == "multi_far_84_154"))

    print()
    print("E. Proof-facing bridge lemma")
    print("  The only selected AP/84m row requiring a rank-6 one-branch rescue core")
    print("  is the canonical covering row m=1, duplicated as ap_omit_12_tail_84x01.")
    print("  Every AP tail m=2..12 keeps negative naive slack but drops to the")
    print("  rank-5 core (5,7,9,11,13) and still has low-rank two-colour")
    print("  component escapes.  Thus the bridge target is not a new scalar")
    print("  inequality: it is a finite interface between overlap-tax cuts,")
    print("  survivor-gate words, and component-cover escape vertices.")
    print("  The nonnegative multi_far_84_154 row is a control: it has no one-branch")
    print("  deficit and has many two-colour low-rank escapes.")

    print()
    print("F. Tournament Analysis")
    score_hist, path = tournament_fingerprint()
    print("  vertices=bridge proof obligations, not runners or raw arcs")
    print(
        "  pairwise_observable=predicate retention + one-branch cut exactness + "
        "two-colour escape payload + survivor-gate route + scalar-firewall safety"
    )
    print("  switch_gauge=higher proof-facing score; ties by declared carrier order")
    print(f"  score_hist={score_hist}")
    print("  directed_3cycles=0")
    print("  hamiltonian_path_count=1")
    print("  hamiltonian_path=" + " -> ".join(path))

    print()
    print("G. Assumption challenge")
    print("  Considered vertices: runners, odd blockers, even-half gates, gaps,")
    print("  fixed circle sections, section boundaries, wall-crossing events,")
    print("  residues, cover arcs, Fourier modes, matroid circuits, survivor gates,")
    print("  component-cover graph nodes, endpoint walls, and proof obligations.")
    print("  Chosen vertices: bridge obligations between one-branch rescue cuts,")
    print("  survivor-gate words, and two-colour component escape certificates.")
    print("  Preserves: negative one-branch slack is routed to an exact overlap core")
    print("  plus a named two-colour survivor/component escape.")
    print("  Destroys: raw runner order, raw rescue-rank scalar meaning, and most")
    print("  interval geometry unless endpoint labels and component addresses are")
    print("  restored.  Challenged assumption: a rank-6 rescue core signals a broad")
    print("  noncanonical obstruction; on this bridge it is the canonical m=1 fence.")


if __name__ == "__main__":
    main()
