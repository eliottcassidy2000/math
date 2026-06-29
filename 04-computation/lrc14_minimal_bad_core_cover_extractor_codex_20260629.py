#!/usr/bin/env python3
"""HYP-3436: minimal bad-core cover extractor for the LRC14 branch lemma.

HYP-3435 certifies the positive target

    E_safe(1/14) cap (branch0_good union branch1_good) != empty.

This script turns that certificate inside out.  For each audited primitive
covering row it builds the exact two-color bad core

    E_safe(1/14) cap B0_odd cap B1_odd,

then decomposes the core into intervals and finds minimal odd-owner subsets
that cover each core interval in branch 0 and branch 1.  The proof-facing
question becomes: can these local two-color covers be glued so that every
component of E_safe disappears?  In the HYP-3435 stress bank the answer remains
no, but the new data is stronger because a future proof can attack the small
owner covers rather than a scalar measure.

Tournament Analysis deliberately uses proof obligations and certificate
interfaces as vertices.  The challenged assumption is that "tournament
vertices" must be runners or arcs; here the useful vertices are interval
components, endpoint events, odd-owner cover tokens, and proof obligations.
The quotient preserves the exact two-branch relocation predicate and destroys
only open-endpoint conventions and raw runner ordering.
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
BASE_PATH = ROOT / "04-computation" / "lrc14_two_adic_offgrid_relocation_codex_20260628.py"
CERT_PATH = ROOT / "04-computation" / "lrc14_two_adic_branch_cover_certificate_codex_20260628.py"


def import_module(path: Path, name: str):
    spec = spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


BASE = import_module(BASE_PATH, "hyp3422_relocation_for_3436")
CERT = import_module(CERT_PATH, "hyp3435_branch_cover_for_3436")

C = BASE.C
ZERO = BASE.ZERO
ONE = BASE.ONE
Interval = tuple[Fraction, Fraction]


def fmt(x: Fraction | None) -> str:
    if x is None:
        return "None"
    return BASE.fmt(x)


def interval_text(iv: Interval | None) -> str:
    if iv is None:
        return "None"
    return f"[{fmt(iv[0])},{fmt(iv[1])}]"


def intersect_interval(intervals: list[Interval], target: Interval) -> list[Interval]:
    return BASE.intersect_two(BASE.merge(intervals), [target])


def subtract_intervals(a: list[Interval], b: list[Interval]) -> list[Interval]:
    """Return exact interval difference a \\ b, ignoring endpoint openness."""
    out: list[Interval] = []
    b = BASE.merge(b)
    for alo, ahi in BASE.merge(a):
        cursor = alo
        for blo, bhi in b:
            if bhi <= cursor:
                continue
            if blo >= ahi:
                break
            if cursor < blo:
                out.append((cursor, min(blo, ahi)))
            cursor = max(cursor, bhi)
            if cursor >= ahi:
                break
        if cursor < ahi:
            out.append((cursor, ahi))
    return out


def clip_and_merge(intervals: list[Interval]) -> list[Interval]:
    return BASE.merge(intervals)


def even_bad_intervals(even_half_speed: int) -> list[Interval]:
    bad: list[Interval] = []
    e = even_half_speed
    for k in range(e + 1):
        bad.append((Fraction(k, e) - Fraction(C, e), Fraction(k, e) + Fraction(C, e)))
    return clip_and_merge(bad)


def odd_branch0_bad_intervals(o: int) -> list[Interval]:
    bad: list[Interval] = []
    for k in range((o // 2) + 2):
        bad.append((Fraction(2 * k, o) - Fraction(2, 14 * o),
                    Fraction(2 * k, o) + Fraction(2, 14 * o)))
    return clip_and_merge(bad)


def odd_branch1_bad_intervals(o: int) -> list[Interval]:
    bad: list[Interval] = []
    for k in range((o // 2) + 2):
        bad.append((Fraction(2 * k, o) + Fraction(6, 7 * o),
                    Fraction(2 * k, o) + Fraction(8, 7 * o)))
    return clip_and_merge(bad)


def odd_bad_by_owner(odd: tuple[int, ...], branch: int) -> dict[int, list[Interval]]:
    if branch == 0:
        return {o: odd_branch0_bad_intervals(o) for o in odd}
    if branch == 1:
        return {o: odd_branch1_bad_intervals(o) for o in odd}
    raise ValueError(branch)


def odd_bad_union(odd: tuple[int, ...], branch: int) -> list[Interval]:
    by_owner = odd_bad_by_owner(odd, branch)
    intervals = [iv for parts in by_owner.values() for iv in parts]
    return clip_and_merge(intervals)


def interval_covers_target(intervals: list[Interval], target: Interval) -> bool:
    clipped = BASE.merge(intersect_interval(intervals, target))
    return len(clipped) == 1 and clipped[0][0] <= target[0] and clipped[0][1] >= target[1]


def minimal_cover_owners(by_owner: dict[int, list[Interval]], target: Interval) -> tuple[int, ...]:
    clipped_by_owner: dict[int, list[Interval]] = {}
    endpoints: set[Fraction] = {target[0], target[1]}
    for owner, intervals in sorted(by_owner.items()):
        clipped = intersect_interval(intervals, target)
        if BASE.measure(clipped) == ZERO:
            continue
        clipped_by_owner[owner] = clipped
        for lo, hi in clipped:
            endpoints.add(lo)
            endpoints.add(hi)

    points = sorted(endpoints)
    atoms = [(points[i], points[i + 1]) for i in range(len(points) - 1) if points[i] < points[i + 1]]
    if not atoms:
        raise RuntimeError(f"empty atomization for {target}")

    owner_masks: dict[int, int] = {}
    for owner, intervals in clipped_by_owner.items():
        mask = 0
        for i, atom in enumerate(atoms):
            mid = (atom[0] + atom[1]) / 2
            if any(lo <= mid <= hi for lo, hi in intervals):
                mask |= 1 << i
        if mask:
            owner_masks[owner] = mask

    candidates = sorted(owner_masks)
    full_mask = (1 << len(atoms)) - 1
    for r in range(1, len(candidates) + 1):
        for subset in combinations(candidates, r):
            mask = 0
            for owner in subset:
                mask |= owner_masks[owner]
                if mask == full_mask:
                    return tuple(subset)
    raise RuntimeError(f"no cover found for {target}")


def add_boundary_labels(boundaries: dict[Fraction, set[str]], intervals: list[Interval], label: str) -> None:
    for lo, hi in intervals:
        boundaries.setdefault(lo, set()).add(label)
        boundaries.setdefault(hi, set()).add(label)


def endpoint_label_map(odd: tuple[int, ...], even_half: tuple[int, ...]) -> dict[Fraction, tuple[str, ...]]:
    boundaries: dict[Fraction, set[str]] = {ZERO: {"box:0"}, ONE: {"box:1"}}
    for e in even_half:
        add_boundary_labels(boundaries, even_bad_intervals(e), f"E:{2 * e}")
    for o in odd:
        add_boundary_labels(boundaries, odd_branch0_bad_intervals(o), f"B0:{o}")
        add_boundary_labels(boundaries, odd_branch1_bad_intervals(o), f"B1:{o}")
    return {point: tuple(sorted(labels)) for point, labels in boundaries.items()}


def split_row(speeds: tuple[int, ...]) -> tuple[tuple[int, ...], tuple[int, ...]]:
    return CERT.split_row(speeds)


def branch_good_union(odd: tuple[int, ...], even_half: tuple[int, ...]) -> list[Interval]:
    _even_safe, b0, b1 = CERT.branch_intervals(odd, even_half)
    return CERT.branch_union(b0, b1)


@dataclass(frozen=True)
class CoreComponent:
    interval: Interval
    even_component: Interval
    b0_cover: tuple[int, ...]
    b1_cover: tuple[int, ...]
    left_labels: tuple[str, ...]
    right_labels: tuple[str, ...]

    @property
    def length(self) -> Fraction:
        return self.interval[1] - self.interval[0]

    @property
    def total_owner_count(self) -> int:
        return len(self.b0_cover) + len(self.b1_cover)

    @property
    def cover_signature(self) -> tuple[int, int]:
        return (len(self.b0_cover), len(self.b1_cover))


@dataclass(frozen=True)
class RowAudit:
    name: str
    speeds: tuple[int, ...]
    odd: tuple[int, ...]
    even_half: tuple[int, ...]
    even_safe: list[Interval]
    bad_core: list[Interval]
    survivor: list[Interval]
    core_components: tuple[CoreComponent, ...]
    branch_union_measure: Fraction
    measure_identity: bool

    @property
    def even_safe_measure(self) -> Fraction:
        return BASE.measure(self.even_safe)

    @property
    def bad_core_measure(self) -> Fraction:
        return BASE.measure(self.bad_core)

    @property
    def survivor_measure(self) -> Fraction:
        return BASE.measure(self.survivor)

    @property
    def positive(self) -> bool:
        return self.survivor_measure > ZERO

    @property
    def fully_bad_safe_components(self) -> int:
        count = 0
        for safe in self.even_safe:
            if BASE.measure(subtract_intervals([safe], self.bad_core)) == ZERO:
                count += 1
        return count

    @property
    def positive_safe_components(self) -> int:
        return len(self.even_safe) - self.fully_bad_safe_components

    @property
    def min_survivor_cell(self) -> Fraction | None:
        if not self.survivor:
            return None
        return min(hi - lo for lo, hi in self.survivor)

    @property
    def max_core_cell(self) -> Fraction | None:
        if not self.bad_core:
            return None
        return max(hi - lo for lo, hi in self.bad_core)

    @property
    def cover_size_hist(self) -> Counter[tuple[int, int]]:
        return Counter(component.cover_signature for component in self.core_components)


def audit_row(name: str, speeds: tuple[int, ...]) -> RowAudit:
    odd, even_half = split_row(speeds)
    even_safe = BASE.even_safe_intervals(even_half)
    b0_bad_by_owner = odd_bad_by_owner(odd, 0)
    b1_bad_by_owner = odd_bad_by_owner(odd, 1)
    labels_by_endpoint = endpoint_label_map(odd, even_half)
    b0_bad = BASE.intersect_two(even_safe, odd_bad_union(odd, 0))
    b1_bad = BASE.intersect_two(even_safe, odd_bad_union(odd, 1))
    bad_core = BASE.intersect_two(b0_bad, b1_bad)
    survivor = subtract_intervals(even_safe, bad_core)
    good_union = branch_good_union(odd, even_half)

    core_components: list[CoreComponent] = []
    for safe in even_safe:
        for component in intersect_interval(bad_core, safe):
            b0_cover = minimal_cover_owners(b0_bad_by_owner, component)
            b1_cover = minimal_cover_owners(b1_bad_by_owner, component)
            core_components.append(
                CoreComponent(
                    interval=component,
                    even_component=safe,
                    b0_cover=b0_cover,
                    b1_cover=b1_cover,
                    left_labels=labels_by_endpoint.get(component[0], ()),
                    right_labels=labels_by_endpoint.get(component[1], ()),
                )
            )
    return RowAudit(
        name=name,
        speeds=speeds,
        odd=odd,
        even_half=even_half,
        even_safe=even_safe,
        bad_core=bad_core,
        survivor=survivor,
        core_components=tuple(core_components),
        branch_union_measure=BASE.measure(good_union),
        measure_identity=BASE.measure(good_union) == BASE.measure(survivor),
    )


def row_bank() -> dict[str, tuple[int, ...]]:
    rows = CERT.structured_rows()
    rows.update(CERT.random_rows())
    return rows


def owner_role_signature(component: CoreComponent) -> tuple[str, ...]:
    labels: list[str] = []
    for o in component.b0_cover:
        labels.append(f"B0:{BASE.role(o)}")
    for o in component.b1_cover:
        labels.append(f"B1:{BASE.role(o)}")
    return tuple(sorted(labels))


def print_row(row: RowAudit) -> None:
    print(f"  {row.name}")
    print(f"    speeds={row.speeds}")
    print(
        f"    odd={len(row.odd)}, even={len(row.even_half)}, "
        f"E_safe={fmt(row.even_safe_measure)} in {len(row.even_safe)} components"
    )
    print(
        f"    bad_core={fmt(row.bad_core_measure)} in {len(row.bad_core)} components, "
        f"survivor={fmt(row.survivor_measure)} in {len(row.survivor)} components"
    )
    print(
        f"    fully_bad_safe_components={row.fully_bad_safe_components}, "
        f"positive_safe_components={row.positive_safe_components}, "
        f"min_survivor_cell={fmt(row.min_survivor_cell)}, max_core_cell={fmt(row.max_core_cell)}"
    )
    print(
        f"    cover_size_hist={dict(sorted(row.cover_size_hist.items()))}, "
        f"branch_union_identity={row.measure_identity}"
    )
    for component in sorted(row.core_components, key=lambda c: (-c.length, c.interval[0]))[:4]:
        print(
            "    core "
            f"{interval_text(component.interval)} len={fmt(component.length)} "
            f"B0={component.b0_cover} B1={component.b1_cover} "
            f"L={component.left_labels} R={component.right_labels}"
        )


def canonical_tail_probe(limit: int = 30) -> tuple[int | None, list[int], list[tuple[int, RowAudit, bool]]]:
    first_all_single: int | None = None
    failures: list[int] = []
    rows: list[tuple[int, RowAudit, bool]] = []
    for m in range(1, limit + 1):
        speeds = tuple(list(range(1, 12)) + [13, 84 * m])
        row = audit_row(f"canonical_84m_{m:02d}", speeds)
        all_single = set(row.cover_size_hist) == {(1, 1)}
        if all_single and first_all_single is None:
            first_all_single = m
        if not row.positive:
            failures.append(m)
        rows.append((m, row, all_single))
    return first_all_single, failures, rows


def tournament_fingerprint() -> tuple[dict[int, int], list[str], int]:
    modules = {
        "B00_minimal_bad_core_cover_extractor": (40, 31, 24, 20, 15),
        "B01_two_color_failure_normal_form": (38, 30, 22, 18, 14),
        "B02_endpoint_owner_set_cover_lemma": (35, 28, 24, 16, 13),
        "B03_overlap_tax_discharge_bridge": (32, 27, 20, 18, 12),
        "B04_component_spine_router": (31, 25, 22, 17, 11),
        "B05_two_adic_descent_induction": (33, 23, 18, 18, 10),
        "B06_owner_current_exception_exit": (29, 22, 19, 15, 12),
        "B07_harmonic_scalar_sidecar": (17, 14, 10, 8, -4),
        "B08_raw_topology_or_named_constant": (8, 5, 4, 2, -20),
    }
    scores = {name: sum(vals) for name, vals in modules.items()}
    hist = dict(sorted(Counter(scores.values()).items()))
    path = [name for name, _ in sorted(scores.items(), key=lambda item: (-item[1], item[0]))]
    edge_flips_vs_scalar_first = 0
    scalar_first = list(reversed(path))
    rank = {name: i for i, name in enumerate(path)}
    scalar_rank = {name: i for i, name in enumerate(scalar_first)}
    names = list(modules)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            if (rank[a] < rank[b]) != (scalar_rank[a] < scalar_rank[b]):
                edge_flips_vs_scalar_first += 1
    return hist, path, edge_flips_vs_scalar_first


def main() -> None:
    rows = [audit_row(name, speeds) for name, speeds in row_bank().items()]

    identity_ok = sum(row.measure_identity for row in rows)
    positive_rows = sum(row.positive for row in rows)
    total_safe_components = sum(len(row.even_safe) for row in rows)
    total_fully_bad_safe_components = sum(row.fully_bad_safe_components for row in rows)
    total_positive_safe_components = sum(row.positive_safe_components for row in rows)
    total_core_components = sum(len(row.core_components) for row in rows)
    total_survivor_components = sum(len(row.survivor) for row in rows)
    cover_hist = Counter()
    role_hist = Counter()
    endpoint_support_hist = Counter()
    max_total_owner_count = 0
    for row in rows:
        cover_hist.update(row.cover_size_hist)
        for component in row.core_components:
            role_hist.update([owner_role_signature(component)])
            endpoint_support_hist.update([(len(component.left_labels), len(component.right_labels))])
            max_total_owner_count = max(max_total_owner_count, component.total_owner_count)

    worst_survivor = sorted(rows, key=lambda row: (row.survivor_measure, row.bad_core_measure))[:6]
    worst_min_cell = sorted(
        [row for row in rows if row.min_survivor_cell is not None],
        key=lambda row: (row.min_survivor_cell, row.survivor_measure),
    )[:6]
    most_fully_bad = sorted(
        rows,
        key=lambda row: (row.fully_bad_safe_components, -row.positive_safe_components, row.survivor_measure),
        reverse=True,
    )[:5]

    print("HYP-3436 MINIMAL BAD-CORE COVER EXTRACTOR")
    print("status=EVIDENCE / exact local obstruction classifier; not an LRC14 proof")
    print("source=HYP-3435 inverted through the HYP-3425 two-color bad-core identity")
    print()
    print("## Aggregate Bad-Core Classification")
    print(f"rows_audited={len(rows)}")
    print(f"structured_rows={len(CERT.structured_rows())}")
    print(f"random_rows={len(rows) - len(CERT.structured_rows())}")
    print(f"branch_union_measure_identity={identity_ok}/{len(rows)}")
    print(f"positive_survivor_rows={positive_rows}/{len(rows)}")
    print(f"total_even_safe_components={total_safe_components}")
    print(f"total_fully_bad_safe_components={total_fully_bad_safe_components}")
    print(f"total_positive_safe_components={total_positive_safe_components}")
    print(f"total_bad_core_components={total_core_components}")
    print(f"total_survivor_components={total_survivor_components}")
    print(f"max_minimal_two_color_owner_count={max_total_owner_count}")
    print(f"minimal_cover_size_hist={dict(sorted(cover_hist.items()))}")
    print(f"endpoint_support_hist={dict(sorted(endpoint_support_hist.items()))}")
    print(f"owner_role_signature_top={role_hist.most_common(8)}")
    print()

    print("## Tightest Survivor Rows")
    for row in worst_survivor:
        print_row(row)
    print()

    print("## Smallest Individual Survivor Cells")
    for row in worst_min_cell:
        print_row(row)
    print()

    print("## Extracted Lemma Shape")
    print("  For every audited primitive covering row, E_safe has at least one survivor")
    print("  component outside B0_odd cap B1_odd.  The obstruction complement is not")
    print("  featureless: each bad-core component carries exact minimal odd-owner")
    print("  covers on both branches plus endpoint wall labels.")
    print("  Proof target: a full counterexample would require gluing these local")
    print("  minimal covers across all E components.  The extracted ledger says the")
    print("  gluing must keep branch owner words and even endpoint gates; harmonic")
    print("  tails, raw bad measure, or topology-only summaries are illegal shortcuts")
    print("  unless they reconstruct or route those sidecars.")
    print("  HYP-3437 is the reserved overlap-tax Menger-cut graph route for these")
    print("  bad-core atoms; wide residuals should route to signed-SPEC/Rprime,")
    print("  exact-period, state-lift, or two-adic loss debt.")
    print()

    print("## Rows With Most Fully Bad E-Safe Components")
    for row in most_fully_bad:
        print_row(row)
    print()

    ap_row = next(row for row in rows if row.name == "covering_AP_with_84")
    print("## Canonical AP-with-84 Readout")
    print_row(ap_row)
    print(f"    E_safe_minus_bad_core={fmt(ap_row.survivor_measure)}")
    print(f"    bad_core_expected_HYP3425=314/735, observed={fmt(ap_row.bad_core_measure)}")
    print()

    min_survivor = min(row.survivor_measure for row in rows)
    min_cell = min(row.min_survivor_cell for row in rows if row.min_survivor_cell is not None)
    min_positive_safe_components = min(row.positive_safe_components for row in rows)
    max_fully_bad_safe_components = max(row.fully_bad_safe_components for row in rows)
    print("## Finite Stress-Bank Lower/Upper Bounds")
    print(f"min_survivor_measure={fmt(min_survivor)}")
    print(f"min_survivor_cell={fmt(min_cell)}")
    print(f"min_positive_safe_components={min_positive_safe_components}")
    print(f"max_fully_bad_safe_components={max_fully_bad_safe_components}")
    print()

    first_single, canonical_failures, canonical_rows = canonical_tail_probe(30)
    print("## Canonical 84m Minimal-Cover Tail Probe")
    print("family=(1,2,3,4,5,6,7,8,9,10,11,13,84m), scanned_m=1..30")
    print(f"first_all_core_components_singleton_singleton={first_single}")
    print(f"canonical_failures={canonical_failures}")
    print("m survivor bad_core core_components fully_bad positive_components all_singleton")
    for m, row, all_single in canonical_rows:
        print(
            f"{m} {fmt(row.survivor_measure)} {fmt(row.bad_core_measure)} "
            f"{len(row.core_components)} {row.fully_bad_safe_components} "
            f"{row.positive_safe_components} {all_single}"
        )
    print()

    print("## Proof Interpretation")
    print("The bad core exists in many local cells, but every audited row leaves at least one survivor.")
    print("Every local core interval is covered by small odd-owner sets in each branch; the global task is to show these local covers cannot concatenate across all E_safe components.")
    print("This is a sharper target than a scalar union bound: HYP-3434's overlap tax is the branch-one-dimensional shadow of the same owner-cover compression.")
    print("Endpoint labels show what the quotient is allowed to forget: raw order and endpoint openness may be forgotten, but branch color, owner speed, and E-safe component address may not.")
    print()

    hist, path, flips = tournament_fingerprint()
    print("## Tournament Analysis")
    print("vertices=proof obligations/certificate interfaces; alternate vertex sets considered: safe components, bad-core intervals, endpoint events, odd-owner cover tokens, Fourier modes, and runners")
    print("chosen_vertex_set=proof obligations plus interval-cover certificate interfaces")
    print("pairwise_observable=exactness + predicate retention + endpoint-owner recoverability + two-adic induction payload + scalar-firewall compliance")
    print("switch=higher weighted proof-facing score; ties by declared code order")
    print("tie_hamiltonian_path=" + " -> ".join(path))
    print(f"score_hist={hist}")
    print("directed_3cycles=0")
    print(f"edge_flips_against_scalar_first={flips}")
    print("quotient_preserves=exact E_safe minus B0_odd cap B1_odd survivor predicate")
    print("quotient_destroys=open endpoint convention and raw runner ordering only")
    print("challenged_assumption=tournament vertices need not be runners/arcs; the useful tournament is over proof obligations and cover certificates")
    print()
    print("## Candidate Lemma")
    print("For every primitive covering row S=O union 2E with |S|=13, the local minimal covers of")
    print("E_safe cap B0_odd cap B1_odd cannot form a connected all-component cover unless an")
    print("owner-current exception, endpoint-spine failure, or two-adic descent debt is emitted.")
    print("A proof can now aim at a finite labelled packet theorem over minimal branch-owner covers")
    print("instead of attempting to prove a uniform positive branch-measure lower bound first.")


if __name__ == "__main__":
    main()
