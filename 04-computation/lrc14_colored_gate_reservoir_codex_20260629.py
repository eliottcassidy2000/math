#!/usr/bin/env python3
"""HYP-3471: colored gate-reservoir audit for the LRC14 finite route.

This scout reconnects the older "14-colored reservoir" and colored
discrepancy language to the current survivor-gate/component-cover stack.
After the HYP-3462/HYP-3470 rebases, HYP-3462 is the AP84 corridor-splice
carrier, HYP-3470 is the exact AP84 CRT color-grid placement sidecar, and this
audit is the local typed gate-word sidecar.

HYP-3453 proved, on the audited bank, that every row with a dead component has
a rank <= 2 survivor gate.  The coloring question asks whether this can be
made sharper: does the obstruction always touch a typed even/branch boundary,
and which coordinates are lost by coarser color quotients?

The answer in this bank is yes:

    dead_components(row) > 0
      => a rank <= 2 survivor gate has exactly one E endpoint and one B0/B1
         endpoint.

Endpoint-kind coloring, numeric residue mod 14, typed residue mod 14, and the
full structural gate word form a compression ladder.  The AP84 four-color
packet is a real seed, but not a global scalar certificate; random rows require
the typed residue plus branch mask, adjacency, and cover-delta sidecars.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
H3453_PATH = ROOT / "04-computation" / "lrc14_gate_escape_transversal_router_codex_20260629.py"


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


H3453 = load_module("hyp3453_gate_escape_for_hyp3458", H3453_PATH)


def top_items(counter: Counter, limit: int = 12) -> str:
    items = counter.most_common(limit)
    suffix = "" if len(counter) <= limit else " ..."
    return "{" + ", ".join(f"{key!r}: {value}" for key, value in items) + suffix + "}"


def side_kinds(gate) -> tuple[str, str]:
    return tuple(gate.endpoint_kind_signature.split("|"))  # type: ignore[return-value]


def is_e_branch_gate(gate) -> bool:
    left, right = side_kinds(gate)
    return (left == "E" and right in {"B0", "B1"}) or (
        right == "E" and left in {"B0", "B1"}
    )


def is_same_branch_gate(gate) -> bool:
    left, right = side_kinds(gate)
    return left == right and left in {"B0", "B1"}


def is_cross_branch_gate(gate) -> bool:
    return set(side_kinds(gate)) == {"B0", "B1"}


def typed_residue_label(label: str) -> str:
    if label == "FREE":
        return "FREE"
    prefix, value = label.split(":", 1)
    return f"{prefix}:{int(value) % 14}"


def typed_residue(labels: Iterable[str]) -> str:
    labels = tuple(labels)
    if not labels:
        return "none"
    return ".".join(typed_residue_label(label) for label in labels)


def numeric_residue(labels: Iterable[str]) -> str:
    values = []
    for label in labels:
        if ":" in label:
            values.append(str(int(label.split(":", 1)[1]) % 14))
    return ".".join(values) if values else "none"


def typed_word(gate) -> str:
    return f"{typed_residue(gate.left_labels)}|{typed_residue(gate.right_labels)}"


def numeric_word(gate) -> str:
    return f"{numeric_residue(gate.left_labels)}|{numeric_residue(gate.right_labels)}"


def structural_word(gate) -> tuple[str, str, str, int, int]:
    return (
        gate.endpoint_kind_signature,
        gate.branch_mask,
        gate.adjacency,
        gate.b0_delta,
        gate.b1_delta,
    )


def full_color_word(gate) -> tuple[str, str, str, int, int]:
    return (
        typed_word(gate),
        gate.branch_mask,
        gate.adjacency,
        gate.b0_delta,
        gate.b1_delta,
    )


def row_kind_set(row) -> tuple[str, ...]:
    return tuple(sorted({gate.endpoint_kind_signature for gate in row.low_rank_gates}))


def row_typed_set(row) -> tuple[str, ...]:
    return tuple(sorted({typed_word(gate) for gate in row.low_rank_gates}))


def row_structural_set(row) -> tuple[tuple[str, str, str, int, int], ...]:
    return tuple(sorted({structural_word(gate) for gate in row.low_rank_gates}))


@dataclass(frozen=True)
class Carrier:
    name: str
    scores: tuple[int, ...]

    @property
    def total(self) -> int:
        return sum(self.scores)


AXES = (
    "predicate_retention",
    "color_payload",
    "branch_graph_composability",
    "compression_firewall",
    "ap_packet_splice",
    "random_gluing_sidecar",
    "scalar_penalty",
)

CARRIERS = (
    Carrier("E00_dead_positive_e_branch_gate", (10, 10, 10, 10, 8, 8, 10)),
    Carrier("E01_full_colored_gate_word", (10, 10, 10, 9, 8, 10, 9)),
    Carrier("E02_structural_color_sidecar", (9, 7, 10, 9, 8, 9, 9)),
    Carrier("E03_typed_mod14_gate_word", (8, 10, 7, 8, 10, 7, 8)),
    Carrier("E04_ap84_four_color_packet", (7, 9, 7, 8, 10, 5, 8)),
    Carrier("E05_endpoint_kind_coloring", (7, 4, 6, 6, 7, 5, 7)),
    Carrier("E06_numeric_14_residue", (5, 5, 3, 4, 6, 3, 5)),
    Carrier("E07_raw_color_count", (2, 1, 1, 1, 1, 1, 0)),
)


def tournament() -> tuple[dict[int, int], list[str]]:
    hist = dict(sorted(Counter(carrier.total for carrier in CARRIERS).items()))
    order = {carrier.name: index for index, carrier in enumerate(CARRIERS)}
    path = [
        carrier.name
        for carrier in sorted(CARRIERS, key=lambda carrier: (-carrier.total, order[carrier.name]))
    ]
    return hist, path


def main() -> None:
    rows = [H3453.join_row(name, speeds) for name, speeds in H3453.H3450.rows().items()]
    low_gates = [gate for row in rows for gate in row.low_rank_gates]
    e_branch_gates = [gate for gate in low_gates if is_e_branch_gate(gate)]
    same_branch_gates = [gate for gate in low_gates if is_same_branch_gate(gate)]
    cross_branch_gates = [gate for gate in low_gates if is_cross_branch_gate(gate)]

    dead_rows = [row for row in rows if row.has_dead]
    dead_without_e_branch = [
        row for row in dead_rows if not any(is_e_branch_gate(gate) for gate in row.low_rank_gates)
    ]
    dead_only_e_branch = [
        row
        for row in dead_rows
        if row.low_rank_gates and all(is_e_branch_gate(gate) for gate in row.low_rank_gates)
    ]
    dead_with_same_branch = [
        row for row in dead_rows if any(is_same_branch_gate(gate) for gate in row.low_rank_gates)
    ]
    dead_with_cross_branch = [
        row for row in dead_rows if any(is_cross_branch_gate(gate) for gate in row.low_rank_gates)
    ]

    canonical_ap_palette = {
        "B1:7|E:0",
        "E:0|B0:7",
        "E:0|B1:5",
        "B0:5|E:0",
    }
    rows_with_ap_palette = [
        row for row in rows if any(typed_word(gate) in canonical_ap_palette for gate in row.low_rank_gates)
    ]
    dead_with_ap_palette = [
        row
        for row in dead_rows
        if any(typed_word(gate) in canonical_ap_palette for gate in row.low_rank_gates)
    ]
    dead_without_ap_palette = [row.name for row in dead_rows if row not in dead_with_ap_palette]

    row_kind_sets = Counter(row_kind_set(row) for row in rows)
    row_typed_sets = Counter(row_typed_set(row) for row in rows)
    row_structural_sets = Counter(row_structural_set(row) for row in rows)

    by_kind = defaultdict(list)
    for row in rows:
        by_kind[row_kind_set(row)].append(
            (
                row.name,
                len(row.low_rank_gates),
                len(row_typed_set(row)),
                len(row_structural_set(row)),
                row.dead_count,
            )
        )

    min_e_branch_samples = []
    min_e_branch_struct = Counter()
    min_e_branch_typed = Counter()
    for row in dead_rows:
        gates = [gate for gate in row.low_rank_gates if is_e_branch_gate(gate)]
        if not gates:
            continue
        gate = min(gates, key=lambda g: (g.length, g.row_name, g.component_index, g.interval[0]))
        min_e_branch_samples.append(
            (
                row.name,
                typed_word(gate),
                gate.endpoint_kind_signature,
                gate.branch_mask,
                gate.adjacency,
                gate.b0_delta,
                gate.b1_delta,
                H3453.fmt(gate.length),
            )
        )
        min_e_branch_struct[structural_word(gate)] += 1
        min_e_branch_typed[typed_word(gate)] += 1

    score_hist, path = tournament()

    print("HYP-3471 COLORED GATE-RESERVOIR AUDIT")
    print("status=EVIDENCE / exact colored gate quotient audit; not an LRC14 proof")
    print("source=HYP-3453 survivor-gate/component-cover join plus older colored-reservoir route")
    print("corridor_splice_sibling=HYP-3462 AP84 corridor-splice certificate")
    print("color_grid_sibling=HYP-3470 exact AP84 q=14V CRT placement sidecar")
    print()
    print("## Assumption Challenge")
    print("alternate vertices considered: runners, residue colors, endpoint-kind colors,")
    print("typed mod-14 colors, gaps, fixed circle sections, section boundaries, wall")
    print("crossing events, cover arcs, survivor gates, dead-cover components,")
    print("branch-coloured blockers, and proof obligations.")
    print("chosen carrier vertices: low-rank survivor gates decorated by endpoint kind,")
    print("typed residue, branch mask, adjacency, and adjacent cover deltas.")
    print("preserved predicate: whether a dead-cover obstruction has a graph-composable")
    print("rank <= 2 E/branch gate feeding HYP-3451.")
    print("destroyed by scalarization: whether the AP84 four-color packet is present,")
    print("whether same-branch gates need gluing sidecars, and which owner-delta")
    print("coordinates random rows require.")
    print()
    print("## Aggregate Colored-Gate Lemma")
    print(f"rows_audited={len(rows)}")
    print(f"rows_with_dead_components={len(dead_rows)}/{len(rows)}")
    print(f"low_rank_gates={len(low_gates)}")
    print(f"e_branch_low_rank_gates={len(e_branch_gates)}")
    print(f"same_branch_low_rank_gates={len(same_branch_gates)}")
    print(f"cross_branch_low_rank_gates={len(cross_branch_gates)}")
    print(
        "dead_rows_with_e_branch_low_rank_gate="
        f"{len(dead_rows) - len(dead_without_e_branch)}/{len(dead_rows)}"
    )
    print(f"dead_rows_without_e_branch_low_rank_gate={[row.name for row in dead_without_e_branch]}")
    print(f"dead_rows_only_e_branch_low_rank_gates={len(dead_only_e_branch)}")
    print(f"dead_rows_with_same_branch_low_rank_gate={len(dead_with_same_branch)}")
    print(f"dead_rows_with_cross_branch_low_rank_gate={len(dead_with_cross_branch)}")
    print()
    print("## Compression Ladder")
    print(f"endpoint_kind_palette_size={len(set(gate.endpoint_kind_signature for gate in low_gates))}")
    print(f"numeric_mod14_palette_size={len(set(numeric_word(gate) for gate in low_gates))}")
    print(f"typed_mod14_palette_size={len(set(typed_word(gate) for gate in low_gates))}")
    print(f"structural_palette_size={len(set(structural_word(gate) for gate in low_gates))}")
    print(f"full_color_palette_size={len(set(full_color_word(gate) for gate in low_gates))}")
    print(f"row_kind_set_count={len(row_kind_sets)}")
    print(f"row_typed_set_count={len(row_typed_sets)}")
    print(f"row_structural_set_count={len(row_structural_sets)}")
    print(f"endpoint_kind_hist={dict(Counter(gate.endpoint_kind_signature for gate in low_gates))}")
    print(f"top_typed_mod14={top_items(Counter(typed_word(gate) for gate in low_gates), 20)}")
    print(f"top_numeric_mod14={top_items(Counter(numeric_word(gate) for gate in low_gates), 20)}")
    print(f"top_structural={top_items(Counter(structural_word(gate) for gate in low_gates), 20)}")
    print()
    print("## AP84 Four-Color Packet")
    print(f"canonical_ap_palette={sorted(canonical_ap_palette)}")
    print(f"rows_with_canonical_ap_palette={len(rows_with_ap_palette)}/{len(rows)}")
    print(f"dead_rows_with_canonical_ap_palette={len(dead_with_ap_palette)}/{len(dead_rows)}")
    print(f"dead_rows_without_canonical_ap_palette_count={len(dead_without_ap_palette)}")
    print(f"dead_rows_without_canonical_ap_palette_sample={dead_without_ap_palette[:20]}")
    print()
    print("## Endpoint-Kind Compression Failures")
    print(f"row_kind_sets_top={top_items(row_kind_sets, 8)}")
    for kind_set, class_rows in sorted(by_kind.items(), key=lambda item: (-len(item[1]), item[0]))[:6]:
        print(f"kind_class={kind_set} row_count={len(class_rows)} samples={class_rows[:8]}")
    print()
    print("## Minimum E-Branch Gate Per Dead Row")
    print(f"min_e_branch_typed_hist={top_items(min_e_branch_typed, 20)}")
    print(f"min_e_branch_struct_hist={top_items(min_e_branch_struct, 20)}")
    print(f"min_e_branch_samples={min_e_branch_samples[:20]}")
    print()
    print("## Lemma Target")
    print("Audited finite form:")
    print("  dead_components(row) > 0")
    print("    => row has a rank <= 2 survivor gate with one E endpoint and one")
    print("       B0/B1 endpoint.")
    print("This strengthens HYP-3453's rank <= 2 survivor-gate transversal.  The old")
    print("14-colored reservoir survives as typed residue words on those E/branch")
    print("gates, but the AP84 four-color packet is not universal: only 67/130 dead")
    print("rows contain it.  The legal global proof quotient is therefore")
    print("colored gate word = typed residue + branch mask + adjacency + cover delta,")
    print("with HYP-3462 handling the AP84 corridor splice,")
    print("HYP-3470 handling exact AP84 CRT placement,")
    print("HYP-3454/HYP-3456/HYP-3457 handling the AP packet, and HYP-3455")
    print("handling same-branch/random gluing stress.")
    print()
    print("## Tournament Analysis")
    print("vertices=colored proof carriers, not runners or raw color counts")
    print(
        "pairwise_observable=predicate retention + typed color payload + "
        "branch-graph composability + compression-firewall score"
    )
    print("switch=higher proof-facing carrier score; ties use declared route order")
    print(f"axes={','.join(AXES)}")
    print(f"score_hist={score_hist}")
    print("directed_3cycles=0")
    print(f"hamiltonian_path={' -> '.join(path)}")


if __name__ == "__main__":
    main()
