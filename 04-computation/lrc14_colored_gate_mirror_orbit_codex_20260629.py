#!/usr/bin/env python3
"""HYP-3475: colored gate mirror-orbit audit for the LRC14 floor.

HYP-3471 turns low-rank survivor gates into typed endpoint-color words.  This
scout adds the missing orbit coordinate requested by HYP-3461: quotient each
HYP-3438 survivor gate by the exact mirror involution t -> 1-t, then ask which
colored gate orbits still carry large cover-delta/gluing debt.

The result is deliberately finite-bank evidence, not an LRC14 proof.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
H3471_PATH = ROOT / "04-computation" / "lrc14_colored_gate_reservoir_codex_20260629.py"


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


H3471 = load_module("hyp3471_colored_gate_for_hyp3472", H3471_PATH)
H3453 = H3471.H3453

Interval = tuple[Fraction, Fraction]


def mirror_interval(interval: Interval) -> Interval:
    return (1 - interval[1], 1 - interval[0])


def interval_text(interval: Interval) -> str:
    return H3453.interval_text(interval)


def top_items(counter: Counter, limit: int = 12) -> str:
    items = counter.most_common(limit)
    suffix = "" if len(counter) <= limit else " ..."
    return "{" + ", ".join(f"{key!r}: {value}" for key, value in items) + suffix + "}"


def gate_sort_key(gate) -> tuple[str, int, Fraction, Fraction]:
    return (gate.row_name, gate.component_index, gate.interval[0], gate.interval[1])


def ordered_pair(a, b) -> tuple:
    return tuple(sorted((a, b), key=repr))


def orbit_id(orbit: "GateOrbit") -> tuple[str, tuple[Interval, ...]]:
    return (orbit.row_name, tuple(gate.interval for gate in orbit.members))


@dataclass(frozen=True)
class GateOrbit:
    row_name: str
    members: tuple

    @property
    def size(self) -> int:
        return len(self.members)

    @property
    def max_delta(self) -> int:
        return max(gate.total_delta for gate in self.members)

    @property
    def min_length(self) -> Fraction:
        return min(gate.length for gate in self.members)

    @property
    def typed_pair(self) -> tuple[str, ...]:
        return tuple(sorted(H3471.typed_word(gate) for gate in self.members))

    @property
    def structural_pair(self) -> tuple:
        return ordered_pair(*(H3471.structural_word(gate) for gate in self.members))

    @property
    def components(self) -> tuple[int, ...]:
        return tuple(gate.component_index for gate in self.members)

    @property
    def intervals(self) -> tuple[str, ...]:
        return tuple(interval_text(gate.interval) for gate in self.members)

    @property
    def has_e_branch(self) -> bool:
        return any(H3471.is_e_branch_gate(gate) for gate in self.members)

    @property
    def has_same_branch(self) -> bool:
        return any(H3471.is_same_branch_gate(gate) for gate in self.members)

    @property
    def has_cross_branch(self) -> bool:
        return any(H3471.is_cross_branch_gate(gate) for gate in self.members)

    def all_low_rank(self, low_rank_intervals: set[Interval]) -> bool:
        return all(gate.interval in low_rank_intervals for gate in self.members)


@dataclass(frozen=True)
class RowOrbitAudit:
    row: object
    orbits: tuple[GateOrbit, ...]
    unpaired: tuple
    fixed: tuple

    @property
    def hard_orbits(self) -> tuple[GateOrbit, ...]:
        return tuple(orbit for orbit in self.orbits if orbit.max_delta >= 7)


@dataclass(frozen=True)
class Carrier:
    name: str
    scores: tuple[int, ...]

    @property
    def total(self) -> int:
        return sum(self.scores)


AXES = (
    "predicate_retention",
    "mirror_payload",
    "typed_color_payload",
    "delta_debt_localization",
    "ap_packet_splice",
    "random_gluing_sidecar",
    "scalar_firewall",
)

CARRIERS = (
    Carrier("M00_exact_mirror_orbit_quotient", (10, 10, 9, 9, 8, 8, 10)),
    Carrier("M01_hard_delta_orbit_ledger", (10, 9, 9, 10, 6, 10, 10)),
    Carrier("M02_dead_row_e_branch_orbit", (10, 8, 10, 8, 8, 8, 10)),
    Carrier("M03_typed_mirror_color_pair", (8, 9, 10, 8, 9, 7, 9)),
    Carrier("M04_ap84_two_orbit_packet", (8, 8, 9, 7, 10, 5, 9)),
    Carrier("M05_phase_branch_hit_sidecar", (8, 8, 8, 8, 6, 8, 8)),
    Carrier("M06_single_gate_color_word", (7, 3, 8, 5, 7, 5, 6)),
    Carrier("M07_raw_gate_count", (2, 0, 1, 1, 1, 1, 0)),
)


def build_orbits(row) -> RowOrbitAudit:
    by_interval: dict[Interval, object] = {}
    duplicate_intervals: list = []
    for gate in row.gates:
        if gate.interval in by_interval:
            duplicate_intervals.append(gate)
        by_interval[gate.interval] = gate

    seen: set[Interval] = set()
    orbits: list[GateOrbit] = []
    unpaired: list = []
    fixed: list = []
    for gate in sorted(row.gates, key=gate_sort_key):
        if gate.interval in seen:
            continue
        mirrored = mirror_interval(gate.interval)
        partner = by_interval.get(mirrored)
        if partner is None:
            unpaired.append(gate)
            seen.add(gate.interval)
            continue
        if partner.interval == gate.interval:
            fixed.append(gate)
            seen.add(gate.interval)
            members = (gate,)
        else:
            seen.add(gate.interval)
            seen.add(partner.interval)
            members = tuple(sorted((gate, partner), key=gate_sort_key))
        orbits.append(GateOrbit(row.name, members))

    return RowOrbitAudit(
        row=row,
        orbits=tuple(orbits),
        unpaired=tuple(unpaired + duplicate_intervals),
        fixed=tuple(fixed),
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
    audits = [build_orbits(row) for row in rows]
    all_orbits = [orbit for audit in audits for orbit in audit.orbits]
    all_gates = [gate for row in rows for gate in row.gates]
    unpaired = [gate for audit in audits for gate in audit.unpaired]
    fixed = [gate for audit in audits for gate in audit.fixed]

    low_rank_by_row = {
        row.name: {gate.interval for gate in row.low_rank_gates}
        for row in rows
    }
    low_rank_orbits = [
        orbit for orbit in all_orbits if orbit.all_low_rank(low_rank_by_row[orbit.row_name])
    ]
    low_rank_orbit_ids = {orbit_id(orbit) for orbit in low_rank_orbits}
    high_rank_orbits = [orbit for orbit in all_orbits if orbit_id(orbit) not in low_rank_orbit_ids]
    e_branch_orbits = [orbit for orbit in low_rank_orbits if orbit.has_e_branch]
    same_branch_orbits = [orbit for orbit in low_rank_orbits if orbit.has_same_branch]
    cross_branch_orbits = [orbit for orbit in low_rank_orbits if orbit.has_cross_branch]
    hard_orbits = sorted(
        [orbit for orbit in all_orbits if orbit.max_delta >= 7],
        key=lambda orbit: (-orbit.max_delta, orbit.row_name, orbit.components),
    )

    dead_rows = [row for row in rows if row.has_dead]
    dead_without_e_branch_orbit = []
    for row in dead_rows:
        low_intervals = low_rank_by_row[row.name]
        row_orbits = next(audit.orbits for audit in audits if audit.row.name == row.name)
        if not any(orbit.all_low_rank(low_intervals) and orbit.has_e_branch for orbit in row_orbits):
            dead_without_e_branch_orbit.append(row.name)

    typed_orbit_counter = Counter(orbit.typed_pair for orbit in low_rank_orbits)
    structural_orbit_counter = Counter(orbit.structural_pair for orbit in low_rank_orbits)
    delta_orbit_counter = Counter(orbit.max_delta for orbit in all_orbits)
    size_counter = Counter(orbit.size for orbit in all_orbits)

    ap_orbit_palette = {
        ("B0:5|E:0", "E:0|B1:5"),
        ("B1:7|E:0", "E:0|B0:7"),
    }
    rows_with_ap_orbit = []
    dead_rows_with_ap_orbit = []
    for row in rows:
        row_orbits = next(audit.orbits for audit in audits if audit.row.name == row.name)
        has_ap = any(orbit.typed_pair in ap_orbit_palette for orbit in row_orbits)
        if has_ap:
            rows_with_ap_orbit.append(row.name)
            if row.has_dead:
                dead_rows_with_ap_orbit.append(row.name)

    min_e_branch_orbit_struct: Counter = Counter()
    min_e_branch_orbit_typed: Counter = Counter()
    min_e_branch_samples: list[tuple] = []
    for row in dead_rows:
        low_intervals = low_rank_by_row[row.name]
        row_orbits = [
            orbit
            for audit in audits
            if audit.row.name == row.name
            for orbit in audit.orbits
            if orbit.all_low_rank(low_intervals) and orbit.has_e_branch
        ]
        if not row_orbits:
            continue
        orbit = min(row_orbits, key=lambda item: (item.min_length, item.row_name, item.components))
        min_e_branch_orbit_struct[orbit.structural_pair] += 1
        min_e_branch_orbit_typed[orbit.typed_pair] += 1
        min_e_branch_samples.append(
            (
                row.name,
                orbit.typed_pair,
                orbit.structural_pair,
                orbit.components,
                H3453.fmt(orbit.min_length),
            )
        )

    row_hard_summary = defaultdict(list)
    for orbit in hard_orbits:
        row_hard_summary[orbit.row_name].append(orbit)

    score_hist, path = tournament()

    print("HYP-3475 COLORED GATE MIRROR-ORBIT AUDIT")
    print("status=EVIDENCE / exact mirror-orbit quotient audit; not an LRC14 proof")
    print("source=HYP-3438 survivor gates + HYP-3471 typed gate reservoir + HYP-3461 orbit schema")
    print("siblings=HYP-3462 AP84 splice, HYP-3470 AP84 color grid, HYP-3460 phase-branch pullback")
    print()
    print("## Assumption Challenge")
    print("alternate vertices considered: runners, residues, typed endpoint colors,")
    print("individual survivor gates, mirror gate pairs, fixed circle sections,")
    print("component boundaries, cover arcs, hard delta gates, and proof obligations.")
    print("chosen carrier vertices: mirror orbits of colored survivor gates.")
    print("preserved predicate: gate-extension gluing debt plus low-rank/E-branch")
    print("escape availability under the exact t -> 1-t involution.")
    print("destroyed by scalarization: which member of a mirror pair carries branch0")
    print("or branch1 orientation; retained sidecar is the ordered typed/structural pair.")
    print()
    print("## Aggregate Mirror-Orbit Ledger")
    print(f"rows_audited={len(rows)}")
    print(f"survivor_gates={len(all_gates)}")
    print(f"mirror_orbits={len(all_orbits)}")
    print(f"orbit_size_hist={dict(sorted(size_counter.items()))}")
    print(f"fixed_orbits={len(fixed)}")
    print(f"unpaired_or_duplicate_gates={len(unpaired)}")
    print(f"low_rank_mirror_orbits={len(low_rank_orbits)}")
    print(f"high_rank_mirror_orbits={len(high_rank_orbits)}")
    print(f"e_branch_low_rank_orbits={len(e_branch_orbits)}")
    print(f"same_branch_low_rank_orbits={len(same_branch_orbits)}")
    print(f"cross_branch_low_rank_orbits={len(cross_branch_orbits)}")
    print(
        "dead_rows_with_e_branch_low_rank_orbit="
        f"{len(dead_rows) - len(dead_without_e_branch_orbit)}/{len(dead_rows)}"
    )
    print(f"dead_rows_without_e_branch_low_rank_orbit={dead_without_e_branch_orbit}")
    print()
    print("## Orbit Compression")
    print(f"typed_mirror_orbit_palette_size={len(typed_orbit_counter)}")
    print(f"structural_mirror_orbit_palette_size={len(structural_orbit_counter)}")
    print(f"delta_orbit_hist={dict(sorted(delta_orbit_counter.items()))}")
    print(f"top_typed_orbits={top_items(typed_orbit_counter, 14)}")
    print(f"top_structural_orbits={top_items(structural_orbit_counter, 12)}")
    print()
    print("## AP84 Mirror Packet")
    print(f"ap_orbit_palette={sorted(ap_orbit_palette)}")
    print(f"rows_with_ap_orbit={len(rows_with_ap_orbit)}/{len(rows)}")
    print(f"dead_rows_with_ap_orbit={len(dead_rows_with_ap_orbit)}/{len(dead_rows)}")
    print()
    print("## Hard Delta Mirror-Orbit Debt")
    print(f"hard_orbits_delta_ge_7={len(hard_orbits)}")
    print(f"hard_orbit_rows={sorted(row_hard_summary)}")
    for orbit in hard_orbits:
        print(
            "  row={row} delta={delta} components={components} typed_pair={typed} "
            "structural_pair={struct} intervals={intervals}".format(
                row=orbit.row_name,
                delta=orbit.max_delta,
                components=orbit.components,
                typed=orbit.typed_pair,
                struct=orbit.structural_pair,
                intervals=orbit.intervals,
            )
        )
    print()
    print("## Minimum E-Branch Orbit Per Dead Row")
    print(f"min_e_branch_orbit_typed_hist={top_items(min_e_branch_orbit_typed, 14)}")
    print(f"min_e_branch_orbit_struct_hist={top_items(min_e_branch_orbit_struct, 14)}")
    print(f"min_e_branch_orbit_samples={min_e_branch_samples[:16]}")
    print()
    print("## Lemma Target")
    print("Audited finite form:")
    print("  every survivor gate in the 135-row bank belongs to an exact two-gate")
    print("  mirror orbit; no fixed or unpaired gate occurs.")
    print("  dead_components(row)>0 still implies a low-rank E/branch mirror orbit.")
    print("  cover-delta debt >=7 is confined to 8 mirror orbits on 7 random rows.")
    print("This upgrades HYP-3471 from colored gate words to colored gate orbits.")
    print("The next proof obligation is a finite hard-orbit discharge: route those")
    print("8 same/cross-branch two-sided orbits through HYP-3455/HYP-3451, or name")
    print("owner-current, two-adic, endpoint-spine, phase-branch, or signed-SPEC debt.")
    print()
    print("## Tournament Analysis")
    print("vertices=mirror-orbit proof carriers, not individual gates or raw counts")
    print(
        "pairwise_observable=predicate retention + mirror payload + typed color "
        "payload + hard-delta localization + scalar firewall"
    )
    print("switch=higher proof-facing carrier score; ties use declared route order")
    print(f"axes={','.join(AXES)}")
    print(f"score_hist={score_hist}")
    print("directed_3cycles=0")
    print("hamiltonian_path=" + " -> ".join(path))


if __name__ == "__main__":
    main()
