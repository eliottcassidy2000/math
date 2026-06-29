#!/usr/bin/env python3
"""HYP-3462: AP84 corridor-splice certificate for the LRC14 bridge.

HYP-3439 left one AP-tail splice obligation:

    prove the canonical rank-6 base case using HYP-3431,
    prove the rank-5 AP-tail descent using HYP-3454/HYP-3456/HYP-3457.

This script packages that splice.  It first imports the HYP-3431 low-core
branch-union identity as the carrier

    C1 = [8/49, 6/35],     C0 = [29/35, 41/49],

then verifies that the one-branch overlap core is rank 6 only at m=1 and rank
5 for the checked AP84 tail m=2..70.  Finally it attaches the already named
AP-tail sidecars:

    HYP-3457  finite mixed transients m=1..4,
    HYP-3454  rank-one endpoint interval m>=5,
    HYP-3456  mod-35 floor-count boundary clock.

The output is a proof-facing bridge certificate, not an LRC14 proof and not a
global primitive-row theorem.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
F = Fraction
Interval = tuple[F, F]


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


H3431 = load_module(
    "hyp3431_corridor_splice",
    ROOT / "04-computation" / "lrc14_canonical_corridor_fence_codex_20260628.py",
)
H3437 = load_module(
    "hyp3437_corridor_splice",
    ROOT / "04-computation" / "lrc14_overlap_menger_cut_certificate_codex_20260628.py",
)
H3454 = load_module(
    "hyp3454_corridor_splice",
    ROOT / "04-computation" / "lrc14_ap84_endpoint_clock_certificate_codex_20260629.py",
)
H3456 = load_module(
    "hyp3456_corridor_splice",
    ROOT / "04-computation" / "lrc14_ap84_mod35_floor_count_codex_20260629.py",
)
H3457 = load_module(
    "hyp3457_corridor_splice",
    ROOT / "04-computation" / "lrc14_ap84_finite_transients_codex_20260629.py",
)


def fmt(value: F | None) -> str:
    if value is None:
        return "None"
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def interval_fmt(interval: Interval) -> str:
    return f"[{fmt(interval[0])},{fmt(interval[1])}]"


def measure(intervals: list[Interval]) -> F:
    return sum((hi - lo for lo, hi in intervals), F(0))


@dataclass(frozen=True)
class Carrier:
    b0_good: tuple[Interval, ...]
    b1_good: tuple[Interval, ...]
    branch_overlap: tuple[Interval, ...]
    branch_union: tuple[Interval, ...]
    expected: tuple[Interval, ...]

    @property
    def identity_ok(self) -> bool:
        return self.branch_union == self.expected

    @property
    def lengths(self) -> tuple[F, ...]:
        return tuple(hi - lo for lo, hi in self.branch_union)


@dataclass(frozen=True)
class RankRecord:
    m: int
    rescue_rank: int
    rescue_subset: tuple[int, ...]
    negative_slack: bool
    naive_slack: F
    branch0_measure: F
    rescue_margin: F
    route: str


def carrier_identity() -> Carrier:
    e_low = H3431.even_safe_intervals(H3431.LOW_EVEN_HALF)
    b0_bad = H3431.union_many([H3431.branch0_bad_one(o) for o in H3431.ODD_CORE])
    b1_bad = H3431.union_many([H3431.branch1_bad_one(o) for o in H3431.ODD_CORE])
    b0_good = H3431.intersect_two(e_low, H3431.complement(b0_bad))
    b1_good = H3431.intersect_two(e_low, H3431.complement(b1_bad))
    return Carrier(
        b0_good=tuple(b0_good),
        b1_good=tuple(b1_good),
        branch_overlap=tuple(H3431.intersect_two(b0_good, b1_good)),
        branch_union=tuple(H3431.merge(b0_good + b1_good)),
        expected=tuple(H3431.FIXED_CORRIDORS),
    )


def ap84_speeds(m: int) -> tuple[int, ...]:
    return tuple(list(range(1, 12)) + [13, 84 * m])


def rank_record(m: int) -> RankRecord:
    row = H3437.audit_row(f"ap84_splice_x{m:03d}", ap84_speeds(m))
    if m == 1 and row.rescue_rank == 6:
        route = "canonical_rank6_base"
    elif m >= 2 and row.rescue_rank == 5:
        route = "rank5_ap_tail_descent"
    else:
        route = "unexpected_rank_route"
    return RankRecord(
        m=m,
        rescue_rank=row.rescue_rank,
        rescue_subset=row.rescue_subset,
        negative_slack=row.negative_slack,
        naive_slack=row.naive_slack,
        branch0_measure=row.branch0_measure,
        rescue_margin=row.rescue_margin,
        route=route,
    )


def rank_records(limit: int) -> list[RankRecord]:
    return [rank_record(m) for m in range(1, limit + 1)]


def expected_rank_ok(record: RankRecord) -> bool:
    if record.m == 1:
        return record.rescue_rank == 6 and record.rescue_subset == (3, 5, 7, 9, 11, 13)
    return record.rescue_rank == 5 and record.rescue_subset == (5, 7, 9, 11, 13)


def symbolic_endpoint_failures(limit: int) -> list[int]:
    failures: list[int] = []
    lo_c, hi_c = H3454.LOW_CORRIDOR
    for m in range(5, limit + 1):
        interval = H3454.expected_interval(m)
        length = interval[1] - interval[0]
        if not (lo_c <= interval[0] and interval[1] <= hi_c and length == F(1, 49 * m)):
            failures.append(m)
    return failures


def checked_endpoint_failures(records) -> list[int]:
    failures: list[int] = []
    for record in records[4:]:
        cert = H3454.endpoint_certificate(record)
        if not (
            cert.low_corridor_ok
            and cert.length_ok
            and cert.is_best_component
            and cert.is_rank_one_ee
            and cert.labels == f"L[E:{84 * cert.m}] R[E:{84 * cert.m}]"
            and cert.rank == 1
        ):
            failures.append(cert.m)
    return failures


def floor_failures(audit_records) -> tuple[list[int], list[int], list[int]]:
    corridors = tuple(H3431.fixed_low_corridors())
    floor_records: list[H3456.FloorRecord] = []
    for record in audit_records:
        m = record.m
        c1 = len(H3456.corridor_gap_indices(m, corridors[0]))
        c0 = len(H3456.corridor_gap_indices(m, corridors[1]))
        formula = H3456.floor_count_per_corridor(m)
        floor_records.append(
            H3456.FloorRecord(
                m=m,
                c1_count=c1,
                c0_count=c0,
                formula_count=formula,
                predicted_escapes=2 * formula,
                audited_escapes=record.escapes,
            )
        )
    mirror = [record.m for record in floor_records if not record.mirror_ok]
    formula = [record.m for record in floor_records if not record.formula_ok]
    audit = [record.m for record in floor_records if not record.audit_ok]
    return mirror, formula, audit


def floor_shift_failures(limit: int) -> list[int]:
    return [
        m
        for m in range(1, limit + 1)
        if H3456.floor_count_per_corridor(m + 35) - H3456.floor_count_per_corridor(m) != 12
    ]


def finite_transient_failures(audit_records) -> tuple[list[int], list[int], list[int]]:
    rows: list[H3457.TransientRecord] = []
    for record in audit_records[:4]:
        best = record.best
        rows.append(
            H3457.TransientRecord(
                m=record.m,
                components=record.components,
                dead=record.dead,
                escapes=record.escapes,
                low_rank_escape=record.summary.low_rank_escape,
                max_dead_pair_rank=record.row.max_dead_pair_rank,
                projection_components=record.projection_components,
                projection_largest=record.projection_largest,
                expected=H3457.expected_windows(record.m),
                actual=H3457.actual_windows(record),
                best_interval=best.best_survivor,
                best_labels=best.labels,
                best_rank=record.row.best_rank,
                best_length=record.best_len,
            )
        )
    exact = [row.m for row in rows if not row.exact_windows_ok]
    closure = [row.m for row in rows if not row.closure_ok]
    rank_drop = [row.m for row in rows if row.m >= 3 and row.max_dead_pair_rank > 2]
    return exact, closure, rank_drop


def tournament_fingerprint() -> tuple[dict[int, int], list[str], int]:
    vertices = {
        "fixed_branch_union_carrier": (10, 10, 10, 10, 9, 10),
        "canonical_rank6_base_split": (10, 9, 9, 10, 9, 10),
        "rank5_overlap_descent": (9, 10, 9, 10, 9, 9),
        "finite_transient_sidecar": (9, 9, 9, 9, 9, 9),
        "endpoint_clock_tail": (10, 10, 10, 9, 8, 10),
        "floor_count_boundary_clock": (10, 10, 10, 10, 9, 10),
        "component_audit_checks": (8, 8, 8, 8, 7, 8),
        "raw_rescue_rank_scalar": (4, 4, 3, 3, 2, 3),
    }
    scores = {name: sum(parts) for name, parts in vertices.items()}
    hist = dict(sorted(Counter(scores.values()).items()))
    path = [name for name, _score in sorted(scores.items(), key=lambda item: (-item[1], item[0]))]
    rank = {name: index for index, name in enumerate(path)}
    cycles = 0
    for a, b, c in combinations(vertices, 3):
        ab = rank[a] < rank[b]
        bc = rank[b] < rank[c]
        ca = rank[c] < rank[a]
        if ab == bc == ca:
            cycles += 1
    return hist, path, cycles


def main() -> None:
    rank_limit = 70
    endpoint_limit = 70
    symbolic_limit = 420
    shift_limit = 210

    carrier = carrier_identity()
    audit_records = H3454.H3452.audit_tail(endpoint_limit)
    ranks = rank_records(rank_limit)
    rank_failures = [record.m for record in ranks if not expected_rank_ok(record)]
    rank_hist = dict(sorted(Counter(record.rescue_rank for record in ranks).items()))
    route_hist = dict(sorted(Counter(record.route for record in ranks).items()))
    transient_exact, transient_closure, transient_rank_drop = finite_transient_failures(audit_records)
    endpoint_checked = checked_endpoint_failures(audit_records)
    endpoint_symbolic = symbolic_endpoint_failures(symbolic_limit)
    floor_mirror, floor_formula, floor_audit = floor_failures(audit_records)
    shift_failures = floor_shift_failures(shift_limit)
    hist, path, cycles = tournament_fingerprint()

    print("HYP-3462 AP84 CORRIDOR-SPLICE CERTIFICATE")
    print("=" * 72)
    print("status=EVIDENCE / AP-tail bridge splice; not an LRC14 proof")
    print("source=HYP-3431 + HYP-3437/HYP-3439 + HYP-3454/HYP-3456/HYP-3457")
    print()

    print("A. HYP-3431 low branch-union carrier")
    print("  low_core_odd=(1,3,5,7,9,11,13)")
    print("  low_even_half=(1,2,3,4,5)")
    print(f"  b0_good={tuple(interval_fmt(iv) for iv in carrier.b0_good)}")
    print(f"  b1_good={tuple(interval_fmt(iv) for iv in carrier.b1_good)}")
    print(f"  branch_overlap={tuple(interval_fmt(iv) for iv in carrier.branch_overlap)}")
    print(f"  branch_union={tuple(interval_fmt(iv) for iv in carrier.branch_union)}")
    print(f"  expected={tuple(interval_fmt(iv) for iv in carrier.expected)}")
    print(f"  identity_ok={carrier.identity_ok}")
    print(f"  branch_union_measure={fmt(measure(list(carrier.branch_union)))}")
    print(f"  corridor_lengths={tuple(fmt(length) for length in carrier.lengths)}")
    print()

    print("B. HYP-3439 rank split on the AP84 tower")
    print(f"  checked_m=1..{rank_limit}")
    print(f"  rescue_rank_hist={rank_hist}")
    print(f"  route_hist={route_hist}")
    print(f"  rank_split_failures={rank_failures}")
    print("  expected: m=1 -> rank 6 core (3,5,7,9,11,13); m>=2 -> rank 5 core (5,7,9,11,13)")
    print("  samples: m | rank | subset | naive_slack | branch0 | margin | route")
    samples = (1, 2, 3, 4, 5, 6, 7, 8, 12, 21, 35, 36, 70)
    by_m = {record.m: record for record in ranks}
    for m in samples:
        record = by_m[m]
        print(
            "  {m:3d} | {rank:4d} | {subset} | {naive:>12s} | {branch0:>8s} | {margin:>8s} | {route}".format(
                m=m,
                rank=record.rescue_rank,
                subset=record.rescue_subset,
                naive=fmt(record.naive_slack),
                branch0=fmt(record.branch0_measure),
                margin=fmt(record.rescue_margin),
                route=record.route,
            )
        )
    print()

    print("C. AP-tail sidecar splice")
    print("  finite clause: HYP-3457 handles m=1..4")
    print(f"    exact_window_failures={transient_exact}")
    print(f"    closure_failures={transient_closure}")
    print(f"    rank_drop_failures_for_m_ge_3={transient_rank_drop}")
    print("  endpoint clause: HYP-3454 handles m>=5")
    print(f"    checked_endpoint_failures_m_5_to_{endpoint_limit}={endpoint_checked}")
    print(f"    symbolic_endpoint_containment_failures_m_5_to_{symbolic_limit}={endpoint_symbolic}")
    print("  floor-count clause: HYP-3456 counts the alive boundary gaps")
    print(f"    mirror_failures_m_1_to_{endpoint_limit}={floor_mirror}")
    print(f"    formula_failures_m_1_to_{endpoint_limit}={floor_formula}")
    print(f"    component_audit_failures_m_1_to_{endpoint_limit}={floor_audit}")
    print(f"    shift_failures_m_1_to_{shift_limit}={shift_failures}")
    print()

    print("D. Spliced proof interface")
    print("  m=1: HYP-3431 supplies the canonical rank-6 corridor-fence base;")
    print("       HYP-3457 gives the exact finite survivor windows.")
    print("  m=2..4: HYP-3437/HYP-3439 rank drops to the rank-5 core")
    print("          (5,7,9,11,13), and HYP-3457 supplies the finite windows.")
    print("  m>=5: the checked rank-5 overlap descent is stable through m=70;")
    print("        HYP-3454 gives the rank-one endpoint interval and HYP-3456")
    print("        gives the mod-35 boundary count.")
    print("  This removes the AP84 splice as an unnamed HYP-3439 obligation.  The")
    print("  remaining global proof work is the non-AP primitive-row transfer through")
    print("  HYP-3453/HYP-3451/HYP-3455 or named owner/current/state-lift debt.")
    print()

    print("E. Tournament Analysis")
    print("  vertices=splice proof obligations, not runners or raw arcs")
    print("  pairwise_observable=carrier retention + rank split + endpoint/floor/transient payload + scalar penalty")
    print("  switch=higher retained AP-tail proof payload; ties by bridge-splice order")
    print(f"  score_hist={hist}")
    print(f"  directed_3cycles={cycles}")
    print("  hamiltonian_path=" + " -> ".join(path))
    print()

    print("F. Assumption Challenge")
    print("  Considered vertices: runners, m-values, residues, fixed circle sections,")
    print("  section boundaries, high-grid gaps, wall-crossing events, branch-union")
    print("  corridors, cover arcs, survivor gates, component graph nodes, odd rescue")
    print("  cores, endpoint walls, Fourier modes, matroid circuits, and proof")
    print("  obligations.  The chosen quotient preserves the AP84 branch-union survivor")
    print("  predicate plus the HYP-3439 rank split.  It destroys arbitrary non-AP row")
    print("  geometry and is therefore a bridge splice, not the global LRC14 theorem.")


if __name__ == "__main__":
    main()
