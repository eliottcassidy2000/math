#!/usr/bin/env python3
"""HYP-3438: survivor-gate word audit for the LRC14 bad-core route.

HYP-3436 made the local bad-core ledger exact:

    E_safe cap B0_odd cap B1_odd

with minimal B0/B1 odd-owner covers and endpoint labels for each bad-core
component.  This script inspects the complementary mixed components.  Inside
each even-safe component that contains both bad-core blocks and survivor gaps,
it emits a gate word:

    bad block / survivor gap / bad block
    + endpoint labels
    + adjacent B0/B1 cover deltas
    + branch mask of the survivor gap
    + parent even wall.

The goal is to expose the finite local-to-global obstruction needed after
HYP-3436 and HYP-3437: scalar survivor measure is not a certificate; the proof
object is the exact gate word and its legal route to corridor-fence,
endpoint-spine, owner-current, two-adic loss, or signed-SPEC sidecars.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
H3436_PATH = ROOT / "04-computation" / "lrc14_minimal_bad_core_cover_extractor_codex_20260629.py"


def load_h3436():
    spec = spec_from_file_location("h3436_bad_core_for_h3438", H3436_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {H3436_PATH}")
    module = module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


H3436 = load_h3436()
BASE = H3436.BASE

Interval = tuple[Fraction, Fraction]
Label = tuple[str, ...]


def fmt(x: Fraction | None) -> str:
    if x is None:
        return "None"
    return H3436.fmt(x)


def interval_text(interval: Interval | None) -> str:
    if interval is None:
        return "edge"
    return f"[{fmt(interval[0])},{fmt(interval[1])}]"


def point_in(intervals: list[Interval], x: Fraction) -> bool:
    return any(lo < x < hi for lo, hi in intervals)


def label_kinds(labels: Label) -> str:
    if not labels:
        return "none"
    kinds = sorted({label.split(":", 1)[0] for label in labels})
    return "+".join(kinds)


def label_values(labels: Label) -> str:
    if not labels:
        return "none"
    return ".".join(labels)


def cover_pair(cover: tuple[int, ...] | None) -> str:
    if cover is None:
        return "edge"
    return "(" + ",".join(map(str, cover)) + ")"


def cover_delta(left: tuple[int, ...] | None, right: tuple[int, ...] | None) -> int:
    if left is None and right is None:
        return 0
    if left is None:
        return len(right or ())
    if right is None:
        return len(left)
    return len(set(left) ^ set(right))


def branch_mask(mid: Fraction, b0_union: list[Interval], b1_union: list[Interval]) -> str:
    b0_good = not point_in(b0_union, mid)
    b1_good = not point_in(b1_union, mid)
    if b0_good and b1_good:
        return "both"
    if b0_good:
        return "branch0"
    if b1_good:
        return "branch1"
    raise AssertionError("survivor midpoint is still in both bad unions")


@dataclass(frozen=True)
class Segment:
    kind: str
    interval: Interval
    left_labels: Label
    right_labels: Label
    branch_mask: str = ""
    b0_cover: tuple[int, ...] | None = None
    b1_cover: tuple[int, ...] | None = None

    @property
    def length(self) -> Fraction:
        return self.interval[1] - self.interval[0]

    @property
    def cover_total(self) -> int:
        return len(self.b0_cover or ()) + len(self.b1_cover or ())

    @property
    def cover_pair(self) -> tuple[int, int] | None:
        if self.kind != "B":
            return None
        return (len(self.b0_cover or ()), len(self.b1_cover or ()))


@dataclass(frozen=True)
class SurvivorGate:
    row_name: str
    speeds: tuple[int, ...]
    component_index: int
    parent: Interval
    interval: Interval
    gate_word: str
    branch_mask: str
    left_labels: Label
    right_labels: Label
    left_bad: Segment | None
    right_bad: Segment | None
    parent_left_labels: Label
    parent_right_labels: Label

    @property
    def length(self) -> Fraction:
        return self.interval[1] - self.interval[0]

    @property
    def adjacency(self) -> str:
        if self.left_bad is not None and self.right_bad is not None:
            return "two_sided"
        if self.left_bad is not None:
            return "left_bad_edge"
        if self.right_bad is not None:
            return "right_bad_edge"
        return "isolated"

    @property
    def b0_delta(self) -> int:
        return cover_delta(
            self.left_bad.b0_cover if self.left_bad else None,
            self.right_bad.b0_cover if self.right_bad else None,
        )

    @property
    def b1_delta(self) -> int:
        return cover_delta(
            self.left_bad.b1_cover if self.left_bad else None,
            self.right_bad.b1_cover if self.right_bad else None,
        )

    @property
    def total_delta(self) -> int:
        return self.b0_delta + self.b1_delta

    @property
    def endpoint_kind_signature(self) -> str:
        return f"{label_kinds(self.left_labels)}|{label_kinds(self.right_labels)}"

    @property
    def parent_kind_signature(self) -> str:
        return f"{label_kinds(self.parent_left_labels)}|{label_kinds(self.parent_right_labels)}"

    @property
    def cover_signature(self) -> str:
        left = self.left_bad
        right = self.right_bad
        left_sig = "edge" if left is None else f"B0{cover_pair(left.b0_cover)}:B1{cover_pair(left.b1_cover)}"
        right_sig = "edge" if right is None else f"B0{cover_pair(right.b0_cover)}:B1{cover_pair(right.b1_cover)}"
        return f"{left_sig}->{right_sig}"


@dataclass(frozen=True)
class MixedComponent:
    row_name: str
    speeds: tuple[int, ...]
    component_index: int
    parent: Interval
    segments: tuple[Segment, ...]
    survivor_gates: tuple[SurvivorGate, ...]
    parent_left_labels: Label
    parent_right_labels: Label

    @property
    def word(self) -> str:
        return "-".join(segment.kind for segment in self.segments)

    @property
    def length(self) -> Fraction:
        return self.parent[1] - self.parent[0]

    @property
    def survivor_measure(self) -> Fraction:
        return sum((gate.length for gate in self.survivor_gates), Fraction(0))

    @property
    def max_cover_total(self) -> int:
        bad_segments = [segment.cover_total for segment in self.segments if segment.kind == "B"]
        return max(bad_segments) if bad_segments else 0


def build_mixed_components(row: H3436.RowAudit) -> tuple[MixedComponent, ...]:
    b0_union = H3436.odd_bad_union(row.odd, 0)
    b1_union = H3436.odd_bad_union(row.odd, 1)
    labels = H3436.endpoint_label_map(row.odd, row.even_half)
    core_by_interval = {component.interval: component for component in row.core_components}

    mixed: list[MixedComponent] = []
    for component_index, parent in enumerate(row.even_safe):
        local_bad = BASE.intersect_two([parent], row.bad_core)
        local_survivors = H3436.subtract_intervals([parent], local_bad)
        if not local_bad or not local_survivors:
            continue

        pieces: list[tuple[str, Interval]] = [("B", interval) for interval in local_bad]
        pieces.extend(("S", interval) for interval in local_survivors)
        pieces.sort(key=lambda item: (item[1][0], item[1][1], item[0]))

        segments: list[Segment] = []
        for kind, interval in pieces:
            left = labels.get(interval[0], ())
            right = labels.get(interval[1], ())
            if kind == "B":
                core = core_by_interval[interval]
                segments.append(
                    Segment(
                        kind="B",
                        interval=interval,
                        left_labels=left,
                        right_labels=right,
                        b0_cover=core.b0_cover,
                        b1_cover=core.b1_cover,
                    )
                )
            else:
                mid = (interval[0] + interval[1]) / 2
                segments.append(
                    Segment(
                        kind="S",
                        interval=interval,
                        left_labels=left,
                        right_labels=right,
                        branch_mask=branch_mask(mid, b0_union, b1_union),
                    )
                )

        word = "-".join(segment.kind for segment in segments)
        parent_left = labels.get(parent[0], ())
        parent_right = labels.get(parent[1], ())
        gates: list[SurvivorGate] = []
        for idx, segment in enumerate(segments):
            if segment.kind != "S":
                continue
            left_bad = segments[idx - 1] if idx > 0 and segments[idx - 1].kind == "B" else None
            right_bad = segments[idx + 1] if idx + 1 < len(segments) and segments[idx + 1].kind == "B" else None
            gates.append(
                SurvivorGate(
                    row_name=row.name,
                    speeds=row.speeds,
                    component_index=component_index,
                    parent=parent,
                    interval=segment.interval,
                    gate_word=word,
                    branch_mask=segment.branch_mask,
                    left_labels=segment.left_labels,
                    right_labels=segment.right_labels,
                    left_bad=left_bad,
                    right_bad=right_bad,
                    parent_left_labels=parent_left,
                    parent_right_labels=parent_right,
                )
            )

        mixed.append(
            MixedComponent(
                row_name=row.name,
                speeds=row.speeds,
                component_index=component_index,
                parent=parent,
                segments=tuple(segments),
                survivor_gates=tuple(gates),
                parent_left_labels=parent_left,
                parent_right_labels=parent_right,
            )
        )
    return tuple(mixed)


@dataclass(frozen=True)
class Carrier:
    name: str
    scores: tuple[int, ...]

    @property
    def total(self) -> int:
        return sum(self.scores)


AXES = (
    "predicate_retention",
    "finite_gate_exactness",
    "endpoint_payload",
    "cover_delta_payload",
    "corridor_fence_bridge",
    "owner_current_bridge",
    "scalar_firewall_safety",
)

CARRIERS = (
    Carrier("G00_exact_survivor_gate_word", (10, 10, 10, 10, 9, 8, 10)),
    Carrier("G01_adjacent_cover_delta_ledger", (10, 10, 9, 10, 8, 9, 10)),
    Carrier("G02_endpoint_wall_alternation", (9, 10, 10, 8, 9, 8, 10)),
    Carrier("G03_branch_mask_relocation_witness", (10, 9, 8, 8, 9, 8, 9)),
    Carrier("G04_corridor_fence_recognizer", (9, 8, 8, 7, 10, 7, 9)),
    Carrier("G05_overlap_cut_bridge", (8, 8, 8, 8, 8, 8, 9)),
    Carrier("G06_owner_current_exception_router", (8, 7, 8, 8, 7, 10, 9)),
    Carrier("G07_signed_SPEC_route", (8, 7, 7, 7, 7, 7, 9)),
    Carrier("G08_raw_survivor_measure", (4, 3, 2, 1, 1, 1, 0)),
)


def tournament() -> tuple[dict[int, int], list[str]]:
    hist = dict(sorted(Counter(carrier.total for carrier in CARRIERS).items()))
    order = {carrier.name: i for i, carrier in enumerate(CARRIERS)}
    path = [
        carrier.name
        for carrier in sorted(CARRIERS, key=lambda carrier: (-carrier.total, order[carrier.name]))
    ]
    return hist, path


def print_gate(gate: SurvivorGate, indent: str = "  ") -> None:
    print(f"{indent}{gate.row_name} component={gate.component_index} word={gate.gate_word}")
    print(
        f"{indent}  parent={interval_text(gate.parent)} parent_labels="
        f"{label_values(gate.parent_left_labels)}|{label_values(gate.parent_right_labels)}"
    )
    print(
        f"{indent}  survivor={interval_text(gate.interval)} len={fmt(gate.length)} "
        f"branch_mask={gate.branch_mask} endpoint_kinds={gate.endpoint_kind_signature}"
    )
    print(
        f"{indent}  labels={label_values(gate.left_labels)}|{label_values(gate.right_labels)} "
        f"adjacency={gate.adjacency} delta=({gate.b0_delta},{gate.b1_delta})"
    )
    print(f"{indent}  covers={gate.cover_signature}")


def print_component(component: MixedComponent, indent: str = "  ") -> None:
    print(
        f"{indent}{component.row_name} component={component.component_index} "
        f"parent={interval_text(component.parent)} word={component.word}"
    )
    print(
        f"{indent}  survivor_measure={fmt(component.survivor_measure)} "
        f"survivors={len(component.survivor_gates)} max_cover_total={component.max_cover_total}"
    )
    for gate in component.survivor_gates[:4]:
        print_gate(gate, indent + "  ")


def top_items(counter: Counter, limit: int = 12) -> str:
    items = counter.most_common(limit)
    suffix = "" if len(counter) <= limit else " ..."
    return "{" + ", ".join(f"{key!r}: {value}" for key, value in items) + suffix + "}"


def main() -> None:
    rows_by_name = H3436.row_bank()
    audits = [H3436.audit_row(name, speeds) for name, speeds in rows_by_name.items()]
    mixed_components = [component for row in audits for component in build_mixed_components(row)]
    gates = [gate for component in mixed_components for gate in component.survivor_gates]
    total_safe_components = sum(len(row.even_safe) for row in audits)
    full_bad_components = sum(row.fully_bad_safe_components for row in audits)
    clean_components = total_safe_components - full_bad_components - len(mixed_components)

    word_hist = Counter(component.word for component in mixed_components)
    branch_hist = Counter(gate.branch_mask for gate in gates)
    adjacency_hist = Counter(gate.adjacency for gate in gates)
    endpoint_hist = Counter(gate.endpoint_kind_signature for gate in gates)
    parent_hist = Counter(gate.parent_kind_signature for gate in gates)
    delta_hist = Counter((gate.b0_delta, gate.b1_delta) for gate in gates)
    cover_sig_hist = Counter(gate.cover_signature for gate in gates)
    row_gate_hist = Counter(gate.row_name for gate in gates)
    mixed_by_row = Counter(component.row_name for component in mixed_components)
    max_survivors_component = max(mixed_components, key=lambda component: len(component.survivor_gates))
    max_segments_component = max(mixed_components, key=lambda component: len(component.segments))
    smallest_gates = sorted(gates, key=lambda gate: (gate.length, gate.row_name, gate.component_index))[:8]
    hardest_delta = sorted(gates, key=lambda gate: (-gate.total_delta, gate.length, gate.row_name))[:8]
    largest_gate_rows = row_gate_hist.most_common(8)
    tight = [component for component in mixed_components if component.row_name == "covering_AP_with_84"]
    smallest_survivor_row = min(audits, key=lambda row: (row.survivor_measure, row.name))
    hist, path = tournament()

    print("HYP-3438 SURVIVOR-GATE WORD AUDIT")
    print("status=EVIDENCE / exact survivor-gap gate-word classification; not a proof")
    print(
        "source=HYP-3436 mixed bad-core components plus HYP-3437 overlap-cut context; "
        "same 135-row primitive covering bank as HYP-3435/HYP-3436"
    )
    print()

    print("## Assumption Challenge")
    print("  alternate vertices considered: runners, gaps, fixed sections, section boundaries,")
    print("  wall-crossing events, residues, cover arcs, endpoint walls, branch-bad")
    print("  owners, bad-core blocks, survivor gaps, cover-pair deltas, mixed E")
    print("  components, owner-current labels, and proof obligations.")
    print("  chosen carrier vertices: exact survivor gates inside mixed E_safe components.")
    print("  preserved predicate: whether a two-adic branch relocation witness survives.")
    print("  destroyed by scalarization: adjacent B0/B1 cover deltas, endpoint labels,")
    print("  branch mask, parent even wall, and legal sidecar route.")
    print()

    print("## Aggregate Gate Audit")
    print(f"rows_audited={len(audits)}")
    print(f"mixed_E_components={len(mixed_components)}")
    print(f"clean_E_components={clean_components}")
    print(f"fully_bad_E_components={full_bad_components}")
    print(f"survivor_gates_inside_mixed_components={len(gates)}")
    print(f"gate_word_hist={dict(sorted(word_hist.items()))}")
    print(f"survivor_branch_mask_hist={dict(sorted(branch_hist.items()))}")
    print(f"gate_adjacency_hist={dict(sorted(adjacency_hist.items()))}")
    print(f"survivor_endpoint_kind_hist={top_items(endpoint_hist)}")
    print(f"parent_endpoint_kind_hist={top_items(parent_hist)}")
    print(f"cover_delta_hist={dict(sorted(delta_hist.items()))}")
    print(f"top_cover_signature_hist={top_items(cover_sig_hist, 10)}")
    print(f"rows_with_most_gates={largest_gate_rows}")
    print(
        "max_survivors_in_one_mixed_component="
        f"{len(max_survivors_component.survivor_gates)} ({max_survivors_component.row_name}, "
        f"component {max_survivors_component.component_index}, word {max_survivors_component.word})"
    )
    print(
        "max_segments_in_one_mixed_component="
        f"{len(max_segments_component.segments)} ({max_segments_component.row_name}, "
        f"component {max_segments_component.component_index}, word {max_segments_component.word})"
    )
    print(
        "smallest_total_survivor_row="
        f"{smallest_survivor_row.name} survivor={fmt(smallest_survivor_row.survivor_measure)} "
        f"mixed_components={mixed_by_row[smallest_survivor_row.name]}"
    )
    print()

    print("## Tight Canonical Gate Words")
    for component in tight:
        print_component(component)
    print()

    print("## Smallest Survivor Gates")
    for gate in smallest_gates:
        print_gate(gate)
    print()

    print("## Hardest Adjacent Cover Deltas")
    for gate in hardest_delta:
        print_gate(gate)
    print()

    print("## Extracted Lemma Shape")
    print("  Mixed even-safe components have exact gate words.  A proof may use")
    print("  survivor mass only after retaining the local word:")
    print("    adjacent bad-core covers + endpoint labels + branch mask + parent E wall.")
    print("  Legal exits are now finite: a gate is a direct branch relocation witness,")
    print("  a corridor-fence instance, an endpoint-spine/wall certificate, an")
    print("  owner-current exception, a two-adic loss/debt case, an overlap-cut bridge,")
    print("  or a signed-SPEC route.  The raw scalar survivor count is a negative")
    print("  control unless it reconstructs one of those sidecars.")
    print()

    print("## Tournament Analysis")
    print("vertices=survivor-gate proof carriers, not runners or raw gaps")
    print("pairwise_observable=predicate retention + gate exactness + endpoint/cover-delta payload + scalar-firewall safety")
    print("switch=higher proof-facing weighted score; ties use declared carrier order")
    print(f"axes={','.join(AXES)}")
    print(f"score_hist={hist}")
    print("directed_3cycles=0")
    print("hamiltonian_path=" + " -> ".join(path))


if __name__ == "__main__":
    main()
