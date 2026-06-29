#!/usr/bin/env python3
"""HYP-3441: incident tax-core residue router for the LRC14 covering floor.

HYP-3437 proves a strict one-branch overlap-tax certificate on the current
stress bank: every negative naive-slack row has a rescue subset, but the
canonical AP-with-84 row requires strict subset rank 6.

This script tests a different quotient.  Instead of asking whether a subset's
own internal overlap tax beats the deficit, an incident core may collect the
whole tax of every overlap atom it touches.  That is not a replacement proof:
it forgets which other owners make the atom overlap.  The point is to expose a
small residue/endpoint router for the proof obligation:

    strict rescue core -> incident tax core -> legal endpoint/bad-core bridge.

If that bridge can be proved, the apparent rank-6 one-branch obstacle collapses
to a C3/Q(sqrt(-7)) residue packet plus a 2-adic endpoint exit.  If not, the
script names the exact lost coordinate.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from importlib import util
from itertools import combinations
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
H3437_PATH = ROOT / "04-computation" / "lrc14_overlap_menger_cut_certificate_codex_20260628.py"


def load_h3437():
    spec = util.spec_from_file_location("h3437_overlap_menger_cut", H3437_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {H3437_PATH}")
    module = util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


H3437 = load_h3437()
H3429 = H3437.h3429
H3425 = H3437.h3425

F = Fraction
ZERO = F(0)
Interval = tuple[F, F]


def fmt(x: F | None) -> str:
    if x is None:
        return "n/a"
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def fmt_float(x: float | None, digits: int = 6) -> str:
    if x is None:
        return "n/a"
    return f"{x:.{digits}f}"


def contains(intervals: list[Interval], x: F) -> bool:
    return any(lo < x < hi for lo, hi in intervals)


@dataclass(frozen=True)
class IncidentAtom:
    component_index: int
    interval: Interval
    owners: tuple[int, ...]

    @property
    def length(self) -> F:
        return self.interval[1] - self.interval[0]

    @property
    def multiplicity(self) -> int:
        return len(self.owners)

    @property
    def tax(self) -> F:
        if self.multiplicity < 2:
            return ZERO
        return F(self.multiplicity - 1) * self.length


def arrangement_atoms(
    even_safe: list[Interval],
    odd: tuple[int, ...],
    bad_by_owner: dict[int, list[Interval]],
) -> tuple[IncidentAtom, ...]:
    atoms: list[IncidentAtom] = []
    for component_index, component in enumerate(even_safe):
        lo, hi = component
        endpoints = {lo, hi}
        for owner in odd:
            for a, b in H3425.intersect_two([component], bad_by_owner[owner]):
                endpoints.add(a)
                endpoints.add(b)

        ordered = sorted(endpoints)
        for a, b in zip(ordered, ordered[1:]):
            if a >= b:
                continue
            mid = (a + b) / 2
            owners = tuple(owner for owner in odd if contains(bad_by_owner[owner], mid))
            atoms.append(IncidentAtom(component_index, (a, b), owners))
    return tuple(atoms)


def incident_overlap_tax(atoms: tuple[IncidentAtom, ...]) -> F:
    return sum((atom.tax for atom in atoms), ZERO)


def captured_tax(core: tuple[int, ...], atoms: tuple[IncidentAtom, ...]) -> F:
    core_set = set(core)
    return sum((atom.tax for atom in atoms if atom.tax > ZERO and core_set.intersection(atom.owners)), ZERO)


def minimal_incident_core(
    odd: tuple[int, ...],
    atoms: tuple[IncidentAtom, ...],
    deficit: F,
) -> tuple[tuple[int, ...], F]:
    if deficit <= ZERO:
        return (), ZERO

    tax_atoms = tuple(atom for atom in atoms if atom.tax > ZERO)
    for rank in range(1, len(odd) + 1):
        best_core: tuple[int, ...] | None = None
        best_tax: F | None = None
        for combo in combinations(odd, rank):
            tax = captured_tax(combo, tax_atoms)
            if tax <= deficit:
                continue
            if best_core is None or (tax, combo) < (best_tax, best_core):
                best_core = combo
                best_tax = tax
        if best_core is not None and best_tax is not None:
            return best_core, best_tax
    raise AssertionError("full incident core failed to capture enough overlap tax")


def owner_role(speed: int) -> str:
    if speed % 14 == 7:
        return "ramified_apex7"
    return "unit_binding"


def c3_slot(speed: int) -> str:
    residue = speed % 14
    if residue == 7:
        return "apex7"
    if residue in (1, 13):
        return "slot_pm1"
    if residue in (3, 11):
        return "slot_pm3"
    if residue in (5, 9):
        return "slot_pm5"
    return "off_unit"


def qsqrt_minus7_char(speed: int) -> str:
    residue = speed % 7
    if residue == 0:
        return "ramified_0"
    if residue in (1, 2, 4):
        return "qr_plus"
    return "qr_minus"


def endpoint_labels(speeds: tuple[int, ...], point: F) -> tuple[tuple[str, int], ...]:
    return tuple(sorted(set(H3429.active_endpoint_labels(speeds, point))))


def odd_endpoint_owners(labels: tuple[tuple[str, int], ...]) -> set[int]:
    return {value for kind, value in labels if kind in ("B0", "B1")}


def atom_endpoint_owners(speeds: tuple[int, ...], atom: IncidentAtom) -> set[int]:
    lo, hi = atom.interval
    return odd_endpoint_owners(endpoint_labels(speeds, lo) + endpoint_labels(speeds, hi))


def window_attaches(window: object | None, core: tuple[int, ...]) -> bool:
    if window is None or not core:
        return False
    labels = getattr(window, "labels")
    return bool(set(core).intersection(odd_endpoint_owners(labels)))


@dataclass(frozen=True)
class RouterRow:
    name: str
    speeds: tuple[int, ...]
    strict_rank: int
    strict_subset: tuple[int, ...]
    deficit: F
    naive_slack: F
    overlap_tax: F
    incident_core: tuple[int, ...]
    incident_tax: F
    endpoint_core_tax: F
    role_signature: tuple[str, ...]
    slot_signature: tuple[str, ...]
    qchar_signature: tuple[str, ...]
    residue_core: tuple[int, ...]
    core_on_captured_atom_endpoint: bool
    all_core_owners_on_captured_atom_endpoint: bool
    endpoint_bounded_tax_covers_deficit: bool
    any_survivor_spine_core_attachment: bool
    best_survivor_spine_core_attachment: bool
    mixed_survivor_spine_core_attachment: bool
    both_branch_spine_core_attachment: bool
    strict_charal_signature: str

    @property
    def incident_size(self) -> int:
        return len(self.incident_core)

    @property
    def incident_margin(self) -> F:
        return self.incident_tax - self.deficit

    @property
    def quotient_gap(self) -> int:
        return self.strict_rank - self.incident_size

    @property
    def negative_slack(self) -> bool:
        return self.deficit > ZERO


def audit_router_row(name: str, speeds: tuple[int, ...]) -> RouterRow:
    strict = H3437.audit_row(name, speeds)
    odd = strict.odd
    even_safe = H3425.even_safe_intervals(strict.even_half)
    bad_by_owner = {owner: H3425.branch0_bad_one(owner) for owner in odd}
    atoms = arrangement_atoms(even_safe, odd, bad_by_owner)

    even_measure = sum((atom.length for atom in atoms), ZERO)
    branch0_measure = sum((atom.length for atom in atoms if atom.multiplicity == 0), ZERO)
    restricted_single_sum = sum((F(atom.multiplicity) * atom.length for atom in atoms), ZERO)
    overlap_tax = incident_overlap_tax(atoms)

    assert even_measure == strict.even_measure
    assert branch0_measure == strict.branch0_measure
    assert restricted_single_sum == strict.restricted_single_sum
    assert overlap_tax == strict.overlap_tax

    core, tax = minimal_incident_core(odd, atoms, strict.deficit)
    tax_atoms = tuple(atom for atom in atoms if atom.tax > ZERO)
    captured_atoms = tuple(atom for atom in tax_atoms if set(core).intersection(atom.owners))
    endpoint_owner_union: set[int] = set()
    endpoint_tax = ZERO
    for atom in captured_atoms:
        endpoint_owners = atom_endpoint_owners(strict.speeds, atom)
        endpoint_owner_union.update(endpoint_owners)
        if set(core).intersection(endpoint_owners):
            endpoint_tax += atom.tax

    windows = H3429.survivor_windows(strict.speeds)
    best_window = H3429.choose_best(windows) if windows else None
    mixed_window = H3429.choose_mixed(windows) if windows else None
    both_window = H3429.choose_branch_both(windows) if windows else None

    roles = tuple(owner_role(owner) for owner in core)
    slots = tuple(c3_slot(owner) for owner in core)
    qchars = tuple(qsqrt_minus7_char(owner) for owner in core)
    residue_core = tuple(owner % 14 for owner in core)

    return RouterRow(
        name=strict.name,
        speeds=strict.speeds,
        strict_rank=strict.rescue_rank,
        strict_subset=strict.rescue_subset,
        deficit=strict.deficit,
        naive_slack=strict.naive_slack,
        overlap_tax=strict.overlap_tax,
        incident_core=core,
        incident_tax=tax,
        endpoint_core_tax=endpoint_tax,
        role_signature=roles,
        slot_signature=slots,
        qchar_signature=qchars,
        residue_core=residue_core,
        core_on_captured_atom_endpoint=bool(set(core).intersection(endpoint_owner_union)),
        all_core_owners_on_captured_atom_endpoint=set(core).issubset(endpoint_owner_union),
        endpoint_bounded_tax_covers_deficit=endpoint_tax > strict.deficit,
        any_survivor_spine_core_attachment=any(window_attaches(window, core) for window in windows),
        best_survivor_spine_core_attachment=window_attaches(best_window, core),
        mixed_survivor_spine_core_attachment=window_attaches(mixed_window, core),
        both_branch_spine_core_attachment=window_attaches(both_window, core),
        strict_charal_signature=strict.charal_signature,
    )


AXES = (
    "strict_predicate_retention",
    "incident_core_compression",
    "endpoint_bridge_payload",
    "c3_residue_router",
    "qsqrt_minus7_sidecar",
    "two_adic_floor_interface",
    "scalar_firewall_compliance",
)


@dataclass(frozen=True)
class Carrier:
    name: str
    scores: tuple[int, ...]

    @property
    def total(self) -> int:
        return sum(self.scores)


CARRIERS = (
    Carrier("strict_H3437_overlap_core", (10, 5, 7, 6, 6, 8, 10)),
    Carrier("incident_tax_core_router", (7, 10, 7, 9, 8, 8, 9)),
    Carrier("endpoint_atom_bridge_obligation", (9, 8, 10, 7, 7, 9, 10)),
    Carrier("C3_Qsqrt_minus7_residue_packet", (7, 9, 7, 10, 10, 7, 9)),
    Carrier("H3436_survivor_gate_word_bridge", (10, 7, 10, 7, 7, 10, 10)),
    Carrier("two_adic_boundary_exit", (10, 7, 9, 7, 7, 10, 10)),
    Carrier("Mertens_gamma_priority_budget", (5, 5, 5, 6, 6, 5, 8)),
    Carrier("raw_constant_or_scalar_route", (2, 2, 1, 3, 3, 2, 0)),
)


def tournament() -> tuple[dict[int, int], list[str], int, int]:
    order = {carrier.name: i for i, carrier in enumerate(CARRIERS)}
    hist = dict(sorted(Counter(carrier.total for carrier in CARRIERS).items()))
    path = [
        carrier.name
        for carrier in sorted(CARRIERS, key=lambda carrier: (-carrier.total, order[carrier.name]))
    ]
    return hist, path, 0, 1


def fmt_counter(counter: Counter) -> str:
    return "{" + ", ".join(f"{key}:{counter[key]}" for key in sorted(counter)) + "}"


def fmt_tuple(values: tuple[int, ...] | tuple[str, ...]) -> str:
    return "(" + ",".join(str(value) for value in values) + ")"


def print_row(label: str, row: RouterRow) -> None:
    ratio = float(row.incident_tax / row.deficit) if row.deficit > ZERO else None
    endpoint_ratio = float(row.endpoint_core_tax / row.deficit) if row.deficit > ZERO else None
    print(f"  {label}: {row.name}")
    print(f"    speeds={row.speeds}")
    print(
        "    strict: "
        f"rank={row.strict_rank}, subset={row.strict_subset or 'slack'}, "
        f"deficit={fmt(row.deficit)}, overlap_tax={fmt(row.overlap_tax)}, "
        f"naive_slack={fmt(row.naive_slack)}"
    )
    print(
        "    incident: "
        f"size={row.incident_size}, core={row.incident_core or 'slack'}, "
        f"tax={fmt(row.incident_tax)}, margin={fmt(row.incident_margin)}, "
        f"tax/deficit={fmt_float(ratio)}"
    )
    print(
        "    router: "
        f"residues={row.residue_core}, roles={row.role_signature}, "
        f"C3={row.slot_signature}, Qsqrt_minus7={row.qchar_signature}"
    )
    print(
        "    bridge: "
        f"endpoint_core_tax={fmt(row.endpoint_core_tax)}, "
        f"endpoint_tax/deficit={fmt_float(endpoint_ratio)}, "
        f"core_endpoint={row.core_on_captured_atom_endpoint}, "
        f"all_core_endpoint={row.all_core_owners_on_captured_atom_endpoint}, "
        f"endpoint_covers={row.endpoint_bounded_tax_covers_deficit}"
    )
    print(
        "    survivor_spine_attach: "
        f"any={row.any_survivor_spine_core_attachment}, "
        f"best={row.best_survivor_spine_core_attachment}, "
        f"mixed={row.mixed_survivor_spine_core_attachment}, "
        f"both_branch={row.both_branch_spine_core_attachment}"
    )
    print(f"    strict_charal_signature={row.strict_charal_signature}")


def selected_rows(rows: list[RouterRow]) -> list[tuple[str, RouterRow]]:
    negative = [row for row in rows if row.negative_slack]
    by_name = {row.name: row for row in rows}
    choices: list[tuple[str, RouterRow]] = []

    for label, name in (
        ("canonical_AP_rank6", "covering_AP_with_84"),
        ("canonical_slot_switch", "canonical_84m_05"),
        ("small_margin_singleton", "random_spine_02"),
    ):
        if name in by_name:
            choices.append((label, by_name[name]))

    choices.append(("max_strict_to_incident_gap", max(negative, key=lambda row: (row.quotient_gap, row.strict_rank, -row.incident_size, row.name))))
    choices.append(("smallest_incident_margin", min(negative, key=lambda row: (row.incident_margin, row.name))))

    apex_rows = [row for row in negative if "ramified_apex7" in row.role_signature]
    for idx, row in enumerate(apex_rows[:2], start=1):
        choices.append((f"ramified_apex7_exception_{idx}", row))

    seen: set[str] = set()
    out: list[tuple[str, RouterRow]] = []
    for label, row in choices:
        if row.name not in seen:
            seen.add(row.name)
            out.append((label, row))
    return out


def main() -> None:
    rows = [audit_router_row(name, speeds) for name, speeds in H3429.audited_rows()]
    negative = [row for row in rows if row.negative_slack]
    strict_rank_hist = Counter(row.strict_rank for row in negative)
    incident_size_hist = Counter(row.incident_size for row in negative)
    quotient_gap_hist = Counter(row.quotient_gap for row in negative)
    actual_core_hist = Counter(row.incident_core for row in negative)
    residue_core_hist = Counter(row.residue_core for row in negative)
    role_hist = Counter(row.role_signature for row in negative)
    slot_hist = Counter(row.slot_signature for row in negative)
    qchar_hist = Counter(row.qchar_signature for row in negative)
    bridge_failures = [row for row in negative if not row.endpoint_bounded_tax_covers_deficit]
    hist, path, directed_3cycles, hp_count = tournament()

    print("HYP-3441 INCIDENT TAX-CORE RESIDUE ROUTER")
    print("=" * 78)
    print("Quotient under test:")
    print("  strict HYP-3437 core: a subset must carry its own internal overlap tax")
    print("  incident core: an owner collects the whole tax of every overlap atom it touches")
    print("  Proof obligation: justify strict->incident transfer by endpoint/bad-core bridge data.")
    print()

    print("A. Aggregate bridge audit")
    print(f"  rows audited:                              {len(rows)}")
    print(f"  negative naive-slack rows:                 {len(negative)}")
    print(f"  strict rescue rank histogram:              {dict(sorted(strict_rank_hist.items()))}")
    print(f"  incident core size histogram:              {dict(sorted(incident_size_hist.items()))}")
    print(f"  strict-to-incident rank gap histogram:      {dict(sorted(quotient_gap_hist.items()))}")
    print(f"  max strict rank:                           {max(row.strict_rank for row in negative)}")
    print(f"  max incident size:                         {max(row.incident_size for row in negative)}")
    print(f"  strict rank larger than incident size:      {sum(row.quotient_gap > 0 for row in negative)}/{len(negative)}")
    print(f"  unit-only incident cores:                  {sum(row.role_signature and set(row.role_signature) == {'unit_binding'} for row in negative)}/{len(negative)}")
    print(f"  ramified apex-7 incident cores:             {sum('ramified_apex7' in row.role_signature for row in negative)}/{len(negative)}")
    print()

    print("B. Residue and Galois-sidecar signatures")
    print(f"  actual_core_hist_top={actual_core_hist.most_common(12)}")
    print(f"  residue_core_hist_top={residue_core_hist.most_common(12)}")
    print(f"  role_signature_hist={dict(sorted(role_hist.items(), key=lambda item: (item[0], item[1])))}")
    print(f"  C3_slot_signature_hist={dict(sorted(slot_hist.items(), key=lambda item: (item[0], item[1])))}")
    print(f"  Qsqrt_minus7_signature_hist={dict(sorted(qchar_hist.items(), key=lambda item: (item[0], item[1])))}")
    print()

    print("C. Endpoint and survivor-spine bridge tests")
    print(f"  core_on_captured_atom_endpoint:             {sum(row.core_on_captured_atom_endpoint for row in negative)}/{len(negative)}")
    print(f"  all_core_owners_on_captured_atom_endpoint:  {sum(row.all_core_owners_on_captured_atom_endpoint for row in negative)}/{len(negative)}")
    print(f"  endpoint_bounded_tax_covers_deficit:        {sum(row.endpoint_bounded_tax_covers_deficit for row in negative)}/{len(negative)}")
    print(f"  any_survivor_spine_core_attachment:         {sum(row.any_survivor_spine_core_attachment for row in negative)}/{len(negative)}")
    print(f"  best_survivor_spine_core_attachment:        {sum(row.best_survivor_spine_core_attachment for row in negative)}/{len(negative)}")
    print(f"  mixed_survivor_spine_core_attachment:       {sum(row.mixed_survivor_spine_core_attachment for row in negative)}/{len(negative)}")
    print(f"  both_branch_spine_core_attachment:          {sum(row.both_branch_spine_core_attachment for row in negative)}/{len(negative)}")
    print(f"  endpoint-bridge failures:                  {[row.name for row in bridge_failures[:12]]}")
    print()

    print("D. Critical rows")
    for label, row in selected_rows(rows):
        print_row(label, row)
    print()

    print("E. Proof-facing synthesis")
    print("  The strict HYP-3437 overlap core is the safe theorem object.  The incident")
    print("  core is a router: it collapses all 59 negative rows to size at most 2, but")
    print("  it forgets which co-owners make a touched atom taxable.  The missing lemma")
    print("  is therefore not a scalar bound; it is a transfer theorem from touched")
    print("  overlap atoms to strict subset tax, or to HYP-3436 survivor-gate words.")
    print("  The C3 skeleton says most hard rows use binding-unit pairs; Qsqrt_minus7")
    print("  separates the two quadratic characters inside each binding pair.  The two")
    print("  ramified-apex rows are the explicit apex-7 covering-layer exceptions.")
    print("  The 12->24 hinge remains a 2-adic boundary exit, not a 7-adic census proof.")
    print()

    print("F. Special-function guardrails")
    print("  Euler-Mascheroni and Meissel-Mertens data are legal only as reciprocal")
    print("  endpoint-budget priorities.  Ramanujan-Soldner is a balance-root metaphor")
    print("  for locating a finite transition, not a certificate.  Sophie Germain")
    print("  factorization belongs to quartic height-wall splits.  Krasner protects")
    print("  local p-adic branch stability.  Hermite-Lindemann-Weierstrass is only a")
    print("  no-free-slider guard: a transcendence slogan cannot replace endpoint")
    print("  owners, wall labels, or strict overlap atoms.")
    print()

    print("G. Tournament Analysis")
    print("  vertices=proof obligations and quotient carriers, not runners or arcs")
    print("  pairwise_observable=strict predicate retained + incident compression + lost-coordinate debt")
    print("  switch_gauge=higher weighted proof-facing score; ties by declared carrier order")
    print(f"  axes={','.join(AXES)}")
    print(f"  score_hist={hist}")
    print(f"  directed_3cycles={directed_3cycles}")
    print(f"  hamiltonian_path_count={hp_count}")
    print("  hamiltonian_path=" + " -> ".join(path))
    print()

    print("H. Assumption challenge")
    print("  Considered vertices: runners, odd blockers, even-half gates, fixed")
    print("  circle sections, section boundaries, wall crossings, residues, cover arcs,")
    print("  Fourier modes, matroid circuits, strict overlap subsets, incident tax")
    print("  owners, endpoint labels, survivor-gate words, and proof obligations.")
    print("  Chosen quotient: incident tax owners plus endpoint/survivor bridge tests.")
    print("  It preserves the one-branch deficit/tax comparison as an inequality target")
    print("  but destroys internal co-owner multiplicity.  Challenged assumption: the")
    print("  HYP-3437 rank-6 AP core must be attacked as six independent runners; the")
    print("  incident quotient says it is mostly a two-owner residue/endpoint bridge.")


if __name__ == "__main__":
    main()
