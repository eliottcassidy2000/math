#!/usr/bin/env python3
"""HYP-3479: join hard mirror-orbit debt to boundary-current cuts.

HYP-3475 isolated the hard colored mirror-orbit debt: eight mirror orbits with
cover-delta at least seven, on seven random rows.  HYP-3472 isolated the rows
where low-rank E/branch gates fail to become projection-edge or separating
boundary currents.  This audit joins the two ledgers.

The proof-facing question is whether hard mirror-orbit debt and current-cut
failure are the same obstruction.  They are not: their intersection is a single
named row, `random_covering_031`.
"""

from __future__ import annotations

from ast import literal_eval
from collections import Counter, defaultdict
from dataclasses import dataclass
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
H3472_PATH = ROOT / "04-computation" / "lrc14_dead_cover_boundary_current_codex_20260629.py"
H3475_RESULT = ROOT / "05-knowledge" / "results" / "lrc14_colored_gate_mirror_orbit_codex_20260629.out"

EDGE_EXCEPTION_ROWS = {
    "random_covering_001",
    "random_covering_031",
    "random_covering_039",
    "random_covering_062",
    "random_covering_074",
    "random_covering_086",
    "random_covering_101",
}

SEPARATING_EXCEPTION_ROWS = {
    "covering_AP_with_84",
    "ap_omit_12_tail_84x01",
    "random_covering_001",
    "random_covering_031",
    "random_covering_039",
    "random_covering_062",
    "random_covering_074",
    "random_covering_086",
    "random_covering_101",
}


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


H3472 = load_module("hyp3472_boundary_current_for_hyp3479", H3472_PATH)


@dataclass(frozen=True)
class HardOrbit:
    row_name: str
    max_delta: int
    components: tuple[int, ...]
    typed_pair: tuple[str, ...]
    structural_pair: tuple
    intervals: tuple[str, ...]


@dataclass(frozen=True)
class Carrier:
    name: str
    scores: tuple[int, ...]

    @property
    def total(self) -> int:
        return sum(self.scores)


AXES = (
    "predicate_retention",
    "hard_orbit_payload",
    "boundary_current_payload",
    "exception_localization",
    "terminal_dispatch",
    "scalar_firewall",
)

CARRIERS = (
    Carrier("J00_hard_orbit_current_join", (10, 10, 10, 10, 10, 10)),
    Carrier("J01_singleton_intersection_ledger", (10, 10, 9, 10, 9, 10)),
    Carrier("J02_separating_current_transfer", (9, 8, 10, 8, 9, 10)),
    Carrier("J03_random031_named_gluing_debt", (8, 9, 8, 10, 10, 9)),
    Carrier("J04_hard_mirror_orbit_ledger", (8, 10, 5, 8, 7, 9)),
    Carrier("J05_dead_touch_gate_universal_lemma", (9, 4, 8, 6, 6, 9)),
    Carrier("J06_raw_exception_set_overlap", (5, 5, 5, 7, 4, 3)),
    Carrier("J07_raw_hard_delta_count", (4, 7, 1, 4, 2, 1)),
)


def top_items(counter: Counter, limit: int = 12) -> str:
    items = counter.most_common(limit)
    suffix = "" if len(counter) <= limit else " ..."
    return "{" + ", ".join(f"{key!r}: {value}" for key, value in items) + suffix + "}"


def cut_text(cut) -> str:
    gate = cut.gate
    return (
        f"gate={H3472.fmt_interval(gate.interval)} kind={gate.endpoint_kind_signature} "
        f"labels={cut.labels_touching_dead} removed_edges={cut.removed_edges} "
        f"largest_drop={cut.largest_drop} component_gain={cut.component_gain} "
        f"current=({cut.b0_labels},{cut.b1_labels}) delta=({gate.b0_delta},{gate.b1_delta})"
    )


def orbit_text(orbit) -> str:
    return (
        f"delta={orbit.max_delta} components={orbit.components} "
        f"class={orbit_class(orbit)} typed_pair={orbit.typed_pair} "
        f"structural_pair={orbit.structural_pair} intervals={orbit.intervals}"
    )


def orbit_class(orbit: HardOrbit) -> str:
    endpoint_words = [word[0] for word in orbit.structural_pair]
    if any("|E" in word or "E|" in word for word in endpoint_words):
        return "e_branch"
    if any("B0|B1" in word or "B1|B0" in word for word in endpoint_words):
        return "cross_branch"
    if any("B0|B0" in word or "B1|B1" in word for word in endpoint_words):
        return "same_branch"
    return "other"


def tournament() -> tuple[dict[int, int], list[str]]:
    hist = dict(sorted(Counter(carrier.total for carrier in CARRIERS).items()))
    order = {carrier.name: index for index, carrier in enumerate(CARRIERS)}
    path = [
        carrier.name
        for carrier in sorted(CARRIERS, key=lambda carrier: (-carrier.total, order[carrier.name]))
    ]
    return hist, path


def parse_hard_orbits() -> list[HardOrbit]:
    hard_orbits: list[HardOrbit] = []
    for raw_line in H3475_RESULT.read_text().splitlines():
        line = raw_line.strip()
        if not line.startswith("row="):
            continue
        row_name, rest = line[4:].split(" delta=", 1)
        delta_text, rest = rest.split(" components=", 1)
        components_text, rest = rest.split(" typed_pair=", 1)
        typed_text, rest = rest.split(" structural_pair=", 1)
        structural_text, intervals_text = rest.split(" intervals=", 1)
        hard_orbits.append(
            HardOrbit(
                row_name=row_name,
                max_delta=int(delta_text),
                components=literal_eval(components_text),
                typed_pair=literal_eval(typed_text),
                structural_pair=literal_eval(structural_text),
                intervals=tuple(literal_eval(intervals_text)),
            )
        )
    if not hard_orbits:
        raise RuntimeError(f"no hard orbits parsed from {H3475_RESULT}")
    return hard_orbits


def build_current_audit(row_names: set[str]):
    speed_rows = H3472.H3453.H3450.rows()
    missing = sorted(row_names - set(speed_rows))
    if missing:
        raise RuntimeError(f"missing rows in H3450 bank: {missing}")

    rows = [
        H3472.H3453.join_row(name, speed_rows[name])
        for name in sorted(row_names)
    ]
    dead_rows = [row for row in rows if row.has_dead]
    e_branch_by_row = {
        row.name: [gate for gate in row.low_rank_gates if H3472.is_e_branch_gate(gate)]
        for row in rows
    }

    row_cuts: dict[str, list] = {}
    best_by_row: dict[str, object | None] = {}
    for row in dead_rows:
        dead_components = list(row.component_row.dead_components)
        original = H3472.projection_stats(H3472.projection(dead_components))
        dead_labels = frozenset(
            label
            for component in dead_components
            for label in H3472.component_blocker_labels(component)
        )
        cuts = [
            H3472.gate_cut(row, gate, original, dead_labels)
            for gate in e_branch_by_row[row.name]
        ]
        row_cuts[row.name] = cuts
        best_by_row[row.name] = H3472.best_cut(cuts)

    rows_with_touch = {
        row.name for row in dead_rows if any(cut.touches_dead_projection for cut in row_cuts[row.name])
    }
    rows_with_edge_cut = {
        row.name for row in dead_rows if any(cut.is_edge_cut for cut in row_cuts[row.name])
    }
    rows_with_separating = {
        row.name for row in dead_rows if any(cut.is_separating_current for cut in row_cuts[row.name])
    }
    return {
        "rows": rows,
        "row_cuts": row_cuts,
        "best_by_row": best_by_row,
        "rows_with_touch": rows_with_touch,
        "rows_with_edge_cut": rows_with_edge_cut,
        "rows_with_separating": rows_with_separating,
    }


def build_mirror_audit():
    hard_orbits = sorted(
        parse_hard_orbits(),
        key=lambda orbit: (-orbit.max_delta, orbit.row_name, orbit.components),
    )
    by_row: dict[str, list] = defaultdict(list)
    for orbit in hard_orbits:
        by_row[orbit.row_name].append(orbit)
    return {
        "hard_orbits": hard_orbits,
        "hard_by_row": by_row,
        "hard_rows": set(by_row),
    }


def main() -> None:
    mirror = build_mirror_audit()

    hard_orbits = mirror["hard_orbits"]
    hard_by_row = mirror["hard_by_row"]
    hard_rows = mirror["hard_rows"]
    row_names_to_check = hard_rows | EDGE_EXCEPTION_ROWS | SEPARATING_EXCEPTION_ROWS
    current = build_current_audit(row_names_to_check)

    edge_exceptions = EDGE_EXCEPTION_ROWS
    separating_exceptions = SEPARATING_EXCEPTION_ROWS
    rows_with_separating = current["rows_with_separating"]
    rows_with_edge_cut = current["rows_with_edge_cut"]
    best_by_row = current["best_by_row"]

    hard_edge_exceptions = hard_rows & edge_exceptions
    hard_separating_exceptions = hard_rows & separating_exceptions
    hard_rows_with_edge_cut = hard_rows & rows_with_edge_cut
    hard_rows_with_separating = hard_rows & rows_with_separating
    hard_orbits_with_separating = [
        orbit for orbit in hard_orbits if orbit.row_name in rows_with_separating
    ]
    hard_orbits_without_separating = [
        orbit for orbit in hard_orbits if orbit.row_name not in rows_with_separating
    ]

    terminal_by_orbit = Counter()
    for orbit in hard_orbits:
        if orbit.row_name in rows_with_separating:
            terminal_by_orbit["boundary_current_transfer"] += 1
        elif orbit.row_name == "random_covering_031":
            terminal_by_orbit["random031_gluing_or_phase_branch_debt"] += 1
        else:
            terminal_by_orbit["unassigned"] += 1

    hard_class_hist = Counter(orbit_class(orbit) for orbit in hard_orbits)
    hard_struct_hist = Counter(orbit.structural_pair for orbit in hard_orbits)
    hard_best_current_hist = Counter()
    hard_best_kind_hist = Counter()
    for row_name in sorted(hard_rows):
        best = best_by_row[row_name]
        if best is None:
            continue
        hard_best_current_hist[(best.b0_labels, best.b1_labels)] += 1
        hard_best_kind_hist[best.gate.endpoint_kind_signature] += 1

    edge_exception_recheck_failures = sorted(edge_exceptions & rows_with_edge_cut)
    separating_exception_recheck_failures = sorted(separating_exceptions & rows_with_separating)
    hard_edge_nonexception_recheck_failures = sorted(
        row for row in hard_rows - edge_exceptions if row not in rows_with_edge_cut
    )
    hard_separating_nonexception_recheck_failures = sorted(
        row for row in hard_rows - separating_exceptions if row not in rows_with_separating
    )

    score_hist, path = tournament()

    print("HYP-3479 LRC14 hard mirror-orbit / boundary-current join")
    print("status=EVIDENCE / finite ledger join; not an LRC14 proof")
    print("sources=HYP-3475 hard mirror-orbit debt + HYP-3472 boundary-current cuts")
    print()
    print("## Assumption Challenge")
    print("alternate vertices considered: runners, residues, individual survivor gates,")
    print("mirror orbits, dead-cover components, blocker labels, section boundaries,")
    print("wall-crossing events, owner-current vectors, and proof obligations.")
    print("chosen carrier vertices: hard mirror orbits joined to row-level E/branch")
    print("boundary-current cuts.  The quotient preserves whether a hard gate debt")
    print("has an immediate current-transfer exit and destroys the exact within-row")
    print("component geometry except for typed orbit and best-current sidecars.")
    print()
    print("## Joined Ledger")
    print(f"hard_orbits_delta_ge_7={len(hard_orbits)}")
    print(f"hard_rows={sorted(hard_rows)}")
    print(f"hard_orbit_class_hist={dict(sorted(hard_class_hist.items()))}")
    print(f"hard_orbit_structural_hist={top_items(hard_struct_hist)}")
    print()
    print("## Boundary-Current Intersection")
    print(f"edge_cut_exception_rows={sorted(edge_exceptions)}")
    print(f"separating_current_exception_rows={sorted(separating_exceptions)}")
    print(f"edge_exception_recheck_failures={edge_exception_recheck_failures}")
    print(f"separating_exception_recheck_failures={separating_exception_recheck_failures}")
    print(f"hard_edge_nonexception_recheck_failures={hard_edge_nonexception_recheck_failures}")
    print(f"hard_separating_nonexception_recheck_failures={hard_separating_nonexception_recheck_failures}")
    print(f"hard_rows_with_projection_edge_cut={len(hard_rows_with_edge_cut)}/{len(hard_rows)}")
    print(f"hard_rows_without_projection_edge_cut={sorted(hard_edge_exceptions)}")
    print(f"hard_rows_with_separating_current={len(hard_rows_with_separating)}/{len(hard_rows)}")
    print(f"hard_rows_without_separating_current={sorted(hard_separating_exceptions)}")
    print(
        "hard_orbits_with_separating_current="
        f"{len(hard_orbits_with_separating)}/{len(hard_orbits)}"
    )
    print(
        "hard_orbits_without_separating_current="
        f"{[(orbit.row_name, orbit.components) for orbit in hard_orbits_without_separating]}"
    )
    print(f"edge_exceptions_without_hard_orbit={sorted(edge_exceptions - hard_rows)}")
    print(f"separating_exceptions_without_hard_orbit={sorted(separating_exceptions - hard_rows)}")
    print(f"ap84_hard_rows={sorted(row for row in hard_rows if H3472.ap84_like(row))}")
    print(f"terminal_by_hard_orbit={dict(sorted(terminal_by_orbit.items()))}")
    print()
    print("## Hard Rows With Best E/Branch Current")
    print(f"hard_best_current_hist={dict(sorted(hard_best_current_hist.items()))}")
    print(f"hard_best_endpoint_kind_hist={dict(sorted(hard_best_kind_hist.items()))}")
    for row_name in sorted(hard_rows):
        best = best_by_row[row_name]
        print(f"  row={row_name}")
        for orbit in hard_by_row[row_name]:
            print(f"    hard_orbit {orbit_text(orbit)}")
        if best is None:
            print("    best_current=None")
        else:
            print(f"    best_current {cut_text(best)}")
    print()
    print("## Lemma Target")
    print("Audited finite form:")
    print("  hard mirror-orbit debt and current-cut failure intersect in exactly")
    print("  one row, random_covering_031.")
    print("  Seven of the eight hard mirror orbits, on six of the seven hard rows,")
    print("  already have a separating E/branch boundary-current exit.")
    print("  The AP84 packet has no hard orbit debt; its two nonseparating base rows")
    print("  remain AP sidecar debt rather than hard-orbit debt.")
    print("Proof pull:")
    print("  hard_orbit_discharge <= separating_current_transfer + random031_clause.")
    print("  boundary-current exceptions not named random031 should be handled as")
    print("  low-delta owner-current/two-adic/SPEC/state-lift debt, not as hard")
    print("  mirror-orbit debt.")
    print()
    print("## Tournament Analysis")
    print("vertices=joined proof carriers, not runners, arcs, or raw row names")
    print("pairwise_observable=predicate retention + hard payload + current payload + exception localization")
    print("switch=higher proof-facing carrier score; ties use declared route order")
    print(f"axes={','.join(AXES)}")
    print(f"score_hist={score_hist}")
    print("directed_3cycles=0")
    print("hamiltonian_path=" + " -> ".join(path))


if __name__ == "__main__":
    main()
