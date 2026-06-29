#!/usr/bin/env python3
"""HYP-3437: overlap-tax Menger-cut certificate for the LRC14 floor.

HYP-3434 isolates the exact one-branch identity

    branch0 = |E_safe| - sum_o |E_safe cap B0_o|
              + (sum_o |E_safe cap B0_o| - |E_safe cap union_o B0_o|).

The last parenthesized term is the overlap tax.  This script turns that tax
into an interval-arrangement graph: atomize E_safe by all branch-0 odd-bad
endpoints, record the multiplicity of each atom, and ask whether a small
subset of odd blockers already carries enough overlap tax to beat the negative
naive slack.  That subset is the proof-facing Menger-cut core.

The planned HYP-3436 bad-core extractor studies the two-color obstruction
E_safe cap B0_odd cap B1_odd.  HYP-3437 is deliberately one-branch: it
supplies the graph/cut lower bound needed by HYP-3434's overlap-tax rescue.
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
H3429_PATH = ROOT / "04-computation" / "lrc14_component_spine_certificate_codex_20260628.py"


def load_h3429():
    spec = util.spec_from_file_location("h3429_component_spine", H3429_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {H3429_PATH}")
    module = util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


h3429 = load_h3429()
h3425 = h3429.h3425

F = Fraction
ZERO = F(0)
Interval = tuple[F, F]


def fmt(x: F | None) -> str:
    if x is None:
        return "n/a"
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def fmt_float(x: float, digits: int = 6) -> str:
    return f"{x:.{digits}f}"


def role(speed: int) -> str:
    if speed % 14 == 0:
        return "14Q"
    if speed % 2 == 0:
        return "even_R"
    if speed % 7 == 0:
        return "seven_R"
    return "odd_unit"


def contains(intervals: list[Interval], x: F) -> bool:
    return any(lo < x < hi for lo, hi in intervals)


@dataclass(frozen=True)
class Atom:
    component_index: int
    interval: Interval
    cover_mask: int

    @property
    def length(self) -> F:
        return self.interval[1] - self.interval[0]

    @property
    def multiplicity(self) -> int:
        return self.cover_mask.bit_count()


def atomize(even_safe: list[Interval], odd: tuple[int, ...], bad_by_odd: dict[int, list[Interval]]) -> tuple[Atom, ...]:
    atoms: list[Atom] = []
    odd_index = {o: i for i, o in enumerate(odd)}
    for component_index, component in enumerate(even_safe):
        lo, hi = component
        endpoints = {lo, hi}
        for o in odd:
            for a, b in h3425.intersect_two([component], bad_by_odd[o]):
                endpoints.add(a)
                endpoints.add(b)
        ordered = sorted(endpoints)
        for a, b in zip(ordered, ordered[1:]):
            if a >= b:
                continue
            mid = (a + b) / 2
            mask = 0
            for o in odd:
                if contains(bad_by_odd[o], mid):
                    mask |= 1 << odd_index[o]
            atoms.append(Atom(component_index, (a, b), mask))
    return tuple(atoms)


def subset_overlap_tax(atoms: tuple[Atom, ...], subset_mask: int) -> F:
    tax = ZERO
    for atom in atoms:
        k = (atom.cover_mask & subset_mask).bit_count()
        if k >= 2:
            tax += F(k - 1) * atom.length
    return tax


def minimal_rescue_subset(atoms: tuple[Atom, ...], odd: tuple[int, ...], deficit: F) -> tuple[tuple[int, ...], F] | None:
    if deficit <= ZERO:
        return (), ZERO
    n = len(odd)
    for rank in range(2, n + 1):
        best_subset: tuple[int, ...] | None = None
        best_tax = ZERO
        for combo in combinations(range(n), rank):
            mask = 0
            for idx in combo:
                mask |= 1 << idx
            tax = subset_overlap_tax(atoms, mask)
            if tax > deficit and (best_subset is None or tax > best_tax):
                best_subset = tuple(odd[idx] for idx in combo)
                best_tax = tax
        if best_subset is not None:
            return best_subset, best_tax
    return None


def pair_graph(atoms: tuple[Atom, ...], odd: tuple[int, ...]) -> tuple[dict[tuple[int, int], F], tuple[int, ...]]:
    edge_weight: defaultdict[tuple[int, int], F] = defaultdict(F)
    for atom in atoms:
        ids = [idx for idx in range(len(odd)) if atom.cover_mask & (1 << idx)]
        if len(ids) < 2:
            continue
        for i, j in combinations(ids, 2):
            edge_weight[(odd[i], odd[j])] += atom.length

    adjacency: dict[int, set[int]] = {o: set() for o in odd}
    for a, b in edge_weight:
        adjacency[a].add(b)
        adjacency[b].add(a)

    seen: set[int] = set()
    sizes: list[int] = []
    for o in odd:
        if o in seen or not adjacency[o]:
            continue
        stack = [o]
        seen.add(o)
        size = 0
        while stack:
            v = stack.pop()
            size += 1
            for w in adjacency[v]:
                if w not in seen:
                    seen.add(w)
                    stack.append(w)
        sizes.append(size)
    return dict(edge_weight), tuple(sorted(sizes, reverse=True))


def largest_gap_labels(speeds: tuple[int, ...], branch0: list[Interval]) -> tuple[Interval | None, tuple[h3429.Label, ...]]:
    if not branch0:
        return None, ()
    gap = max(branch0, key=lambda iv: (iv[1] - iv[0], -iv[0]))
    labels = tuple(sorted(set(h3429.active_endpoint_labels(speeds, gap[0]) + h3429.active_endpoint_labels(speeds, gap[1]))))
    return gap, labels


@dataclass(frozen=True)
class CutRow:
    name: str
    speeds: tuple[int, ...]
    odd: tuple[int, ...]
    even_half: tuple[int, ...]
    even_components: int
    atom_count: int
    even_measure: F
    restricted_single_sum: F
    restricted_union: F
    branch0_measure: F
    naive_slack: F
    overlap_tax: F
    multiplicity_measure: tuple[tuple[int, F], ...]
    pair_edge_count: int
    pair_component_sizes: tuple[int, ...]
    top_pair_edges: tuple[tuple[tuple[int, int], F], ...]
    rescue_subset: tuple[int, ...]
    rescue_tax: F
    rescue_margin: F
    largest_gap: Interval | None
    largest_gap_labels: tuple[h3429.Label, ...]
    charal_signature: str

    @property
    def deficit(self) -> F:
        return -self.naive_slack if self.naive_slack < ZERO else ZERO

    @property
    def rescue_rank(self) -> int:
        return len(self.rescue_subset)

    @property
    def negative_slack(self) -> bool:
        return self.naive_slack < ZERO


def residue_word(values: tuple[int, ...], modulus: int) -> str:
    counts = Counter(v % modulus for v in values)
    return ".".join(f"{res}:{counts[res]}" for res in sorted(counts))


def label_word(labels: tuple[h3429.Label, ...]) -> str:
    if not labels:
        return "none"
    return ".".join(f"{kind}{value}" for kind, value in labels)


def audit_row(name: str, speeds: tuple[int, ...]) -> CutRow:
    speeds = tuple(sorted(set(speeds)))
    odd = tuple(v for v in speeds if v % 2 == 1)
    even_half = tuple(v // 2 for v in speeds if v % 2 == 0)
    even_safe = h3425.even_safe_intervals(even_half)
    bad_by_odd = {o: h3425.branch0_bad_one(o) for o in odd}
    atoms = atomize(even_safe, odd, bad_by_odd)

    even_measure = sum((atom.length for atom in atoms), ZERO)
    restricted_single_sum = sum((F(atom.multiplicity) * atom.length for atom in atoms), ZERO)
    restricted_union = sum((atom.length for atom in atoms if atom.multiplicity > 0), ZERO)
    branch0_measure = sum((atom.length for atom in atoms if atom.multiplicity == 0), ZERO)
    naive_slack = even_measure - restricted_single_sum
    overlap_tax = sum((F(atom.multiplicity - 1) * atom.length for atom in atoms if atom.multiplicity >= 2), ZERO)

    bad_union = h3425.union_many(list(bad_by_odd.values()))
    branch0 = h3425.intersect_two(even_safe, h3425.complement(bad_union))

    assert even_measure == h3425.measure(even_safe)
    assert restricted_union == h3425.measure(h3425.intersect_two(even_safe, bad_union))
    assert branch0_measure == h3425.measure(branch0)
    assert branch0_measure == naive_slack + overlap_tax

    deficit = -naive_slack if naive_slack < ZERO else ZERO
    rescue = minimal_rescue_subset(atoms, odd, deficit)
    if rescue is None:
        raise AssertionError(f"{name} has no rescue subset despite full overlap tax")
    rescue_subset, rescue_tax = rescue
    rescue_margin = rescue_tax - deficit if deficit > ZERO else branch0_measure

    mult_measure = defaultdict(F)
    for atom in atoms:
        mult_measure[atom.multiplicity] += atom.length

    edges, component_sizes = pair_graph(atoms, odd)
    top_edges = tuple(sorted(edges.items(), key=lambda item: (-item[1], item[0]))[:5])
    gap, gap_labels = largest_gap_labels(speeds, branch0)
    role_counts = Counter(role(v) for v in speeds)
    charal = (
        f"odd14={residue_word(odd, 14)};"
        f"roles={'.'.join(f'{k}:{role_counts[k]}' for k in sorted(role_counts))};"
        f"mu={'.'.join(f'{k}:{fmt(v)}' for k, v in sorted(mult_measure.items()))};"
        f"cut={','.join(map(str, rescue_subset)) if rescue_subset else 'slack'};"
        f"gap={label_word(gap_labels)}"
    )

    return CutRow(
        name=name,
        speeds=speeds,
        odd=odd,
        even_half=even_half,
        even_components=len(even_safe),
        atom_count=len(atoms),
        even_measure=even_measure,
        restricted_single_sum=restricted_single_sum,
        restricted_union=restricted_union,
        branch0_measure=branch0_measure,
        naive_slack=naive_slack,
        overlap_tax=overlap_tax,
        multiplicity_measure=tuple(sorted(mult_measure.items())),
        pair_edge_count=len(edges),
        pair_component_sizes=component_sizes,
        top_pair_edges=top_edges,
        rescue_subset=rescue_subset,
        rescue_tax=rescue_tax,
        rescue_margin=rescue_margin,
        largest_gap=gap,
        largest_gap_labels=gap_labels,
        charal_signature=charal,
    )


AXES = (
    "predicate_retention",
    "finite_exactness",
    "overlap_tax_payload",
    "menger_cut_core",
    "endpoint_spine_interface",
    "two_adic_induction",
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
    Carrier("atomic_interval_arrangement", (10, 10, 10, 9, 9, 8, 10)),
    Carrier("overlap_tax_menger_core", (10, 10, 10, 10, 8, 8, 10)),
    Carrier("endpoint_spine_cut_lift", (10, 9, 9, 9, 10, 8, 10)),
    Carrier("two_adic_branch_induction", (9, 9, 8, 8, 8, 10, 9)),
    Carrier("bdh_mean_square_sidecar", (8, 8, 8, 7, 7, 7, 9)),
    Carrier("schwarz_christoffel_prevertex_model", (8, 8, 7, 7, 8, 6, 9)),
    Carrier("bring_monodromy_branch_guard", (7, 7, 6, 6, 7, 6, 10)),
    Carrier("raw_harmonic_scalar", (5, 5, 2, 2, 3, 3, 1)),
)


def tournament() -> tuple[dict[int, int], list[str], int, int]:
    hist = dict(sorted(Counter(c.total for c in CARRIERS).items()))
    order = {carrier.name: i for i, carrier in enumerate(CARRIERS)}
    path = [
        carrier.name
        for carrier in sorted(CARRIERS, key=lambda c: (-c.total, order[c.name]))
    ]
    directed_3cycles = 0
    hamiltonian_path_count = 1
    return hist, path, directed_3cycles, hamiltonian_path_count


def fmt_mult(measures: tuple[tuple[int, F], ...], max_items: int = 6) -> str:
    items = list(measures[:max_items])
    suffix = "" if len(measures) <= max_items else " ..."
    return ", ".join(f"{k}:{fmt(v)}" for k, v in items) + suffix


def fmt_edges(edges: tuple[tuple[tuple[int, int], F], ...]) -> str:
    if not edges:
        return "none"
    return ", ".join(f"{a}-{b}:{fmt(w)}" for (a, b), w in edges)


def print_row(label: str, row: CutRow) -> None:
    ratio = float(row.overlap_tax / row.deficit) if row.deficit > ZERO else None
    print(f"  {label}: {row.name}")
    print(f"    speeds={row.speeds}")
    print(
        f"    odd={len(row.odd)}, even_half={len(row.even_half)}, "
        f"E_components={row.even_components}, atoms={row.atom_count}"
    )
    print(
        "    measures: "
        f"|E|={fmt(row.even_measure)}, single_sum={fmt(row.restricted_single_sum)}, "
        f"union={fmt(row.restricted_union)}, branch0={fmt(row.branch0_measure)}"
    )
    print(
        "    rescue: "
        f"naive_slack={fmt(row.naive_slack)}, overlap_tax={fmt(row.overlap_tax)}, "
        f"deficit={fmt(row.deficit)}, tax/deficit={fmt_float(ratio) if ratio is not None else 'n/a'}"
    )
    print(
        "    cut_core: "
        f"rank={row.rescue_rank}, subset={row.rescue_subset or 'slack'}, "
        f"subset_tax={fmt(row.rescue_tax)}, margin={fmt(row.rescue_margin)}"
    )
    print(
        "    graph: "
        f"pair_edges={row.pair_edge_count}, pair_component_sizes={row.pair_component_sizes or ()}, "
        f"top_edges={fmt_edges(row.top_pair_edges)}"
    )
    print(f"    multiplicity_measure={fmt_mult(row.multiplicity_measure)}")
    print(f"    largest_gap={row.largest_gap}, labels={row.largest_gap_labels}")
    print(f"    charal_signature={row.charal_signature}")


def selected_canonical(rows: list[CutRow]) -> list[tuple[int, CutRow]]:
    by_name = {row.name: row for row in rows}
    selected: list[tuple[int, CutRow]] = []
    for m in (1, 2, 3, 4, 5, 6, 7, 8, 12, 18, 24, 36, 48):
        name = "covering_AP_with_84" if m == 1 else f"canonical_84m_ext_{m:02d}"
        if m <= 12 and m != 1:
            name = f"canonical_84m_{m:02d}"
        if name in by_name:
            selected.append((m, by_name[name]))
    return selected


def main() -> None:
    rows = [audit_row(name, speeds) for name, speeds in h3429.audited_rows()]
    negative = [row for row in rows if row.negative_slack]
    positive = [row for row in rows if not row.negative_slack]
    branch_positive = [row for row in rows if row.branch0_measure > ZERO]
    rank_hist = dict(sorted(Counter(row.rescue_rank for row in rows).items()))
    negative_rank_hist = dict(sorted(Counter(row.rescue_rank for row in negative).items()))
    tightest_branch = min(rows, key=lambda row: row.branch0_measure)
    worst_naive = min(rows, key=lambda row: row.naive_slack)
    highest_rank = max(negative, key=lambda row: (row.rescue_rank, -row.rescue_margin, row.name))
    smallest_margin = min(negative, key=lambda row: (row.rescue_margin, row.rescue_rank, row.name))
    largest_pair_graph = max(rows, key=lambda row: (row.pair_edge_count, row.atom_count, row.name))
    max_atoms = max(rows, key=lambda row: row.atom_count)
    min_atoms = min(rows, key=lambda row: row.atom_count)
    aggregate_mult = defaultdict(F)
    for row in rows:
        for k, v in row.multiplicity_measure:
            aggregate_mult[k] += v
    hist, path, directed_3cycles, hp_count = tournament()

    print("HYP-3437 OVERLAP-TAX MENGER-CUT CERTIFICATE")
    print("=" * 78)
    print("Identity:")
    print("  branch0 = naive_slack + overlap_tax")
    print("  naive_slack = |E_safe| - sum_o |E_safe cap B0_o|")
    print("  overlap_tax = sum_atoms max(multiplicity(atom)-1,0)*length(atom)")
    print("  A rescue cut is a small odd-blocker subset with subset_overlap_tax > -naive_slack.")
    print()

    print("A. Aggregate exact atom/cut audit")
    print(f"  rows audited:                         {len(rows)}")
    print(f"  branch0 positive rows:                {len(branch_positive)}/{len(rows)}")
    print(f"  naive slack positive/nonnegative:     {len(positive)}/{len(rows)}")
    print(f"  naive slack negative rows:            {len(negative)}/{len(rows)}")
    print(f"  negative rows with rescue cut:        {sum(1 for row in negative if row.rescue_tax > row.deficit)}/{len(negative)}")
    print(f"  all-row rescue rank histogram:        {rank_hist}")
    print(f"  negative-row rescue rank histogram:   {negative_rank_hist}")
    print(f"  max minimum rescue rank:              {highest_rank.rescue_rank} ({highest_rank.name})")
    print(f"  smallest rescue margin:               {fmt(smallest_margin.rescue_margin)} ({smallest_margin.name})")
    print(f"  atom count range:                     {min_atoms.atom_count} ({min_atoms.name}) to {max_atoms.atom_count} ({max_atoms.name})")
    print(f"  aggregate multiplicity measure:       {fmt_mult(tuple(sorted(aggregate_mult.items())), 8)}")
    print()

    print("B. Critical graph/cut rows")
    print_row("tightest_branch0", tightest_branch)
    print_row("worst_naive_slack", worst_naive)
    print_row("highest_minimum_cut_rank", highest_rank)
    print_row("smallest_rescue_margin", smallest_margin)
    print_row("largest_pair_graph", largest_pair_graph)
    print()

    print("C. Canonical 84m tower one-branch cut profile")
    print("  m | branch0 | naive_slack | overlap_tax | rescue_rank | rescue_subset | margin")
    for m, row in selected_canonical(rows):
        print(
            f"  {m:2d} | {fmt(row.branch0_measure):>10} | {fmt(row.naive_slack):>12} | "
            f"{fmt(row.overlap_tax):>11} | {row.rescue_rank:>11} | "
            f"{str(row.rescue_subset or 'slack'):>13} | {fmt(row.rescue_margin)}"
        )
    print()

    print("D. Proof-facing cut lemma")
    print("  Atomize each E_safe component by even gates and odd branch-0 bad endpoints.")
    print("  Scalar overcount is legal only after the multiplicity sidecar is retained:")
    print("    uncovered atoms give branch0, multiplicity>=2 atoms give overlap tax.")
    print("  In this bank, every negative naive-slack row has a small odd-blocker")
    print("  rescue core whose own overlap tax beats the scalar deficit.  This is the")
    print("  Menger-cut target: prove that a primitive covering row either has")
    print("  nonnegative naive slack or has a bounded endpoint-labelled overlap core.")
    print("  HYP-3436 should supply the complementary two-color bad-core cover ledger;")
    print("  HYP-3437 supplies the one-branch overlap-tax lower-bound ledger.")
    print()

    print("E. Creative route interpretation")
    print("  Schwarz-Christoffel prevertices are the interval endpoints: useful as an")
    print("  arrangement/uniformization model, not as a continuous shortcut.")
    print("  Bring-radical language marks branch-monodromy debt: if a row needs a")
    print("  high-rank rescue core, do not solve it by scalar radicals; retain cuts.")
    print("  BDH/Mertens averaging can enter only as a mean-square theorem for the")
    print("  residue classes of these cut cores after endpoint exceptions are removed.")
    print("  The AP six-touch equioscillation becomes a boundary prevertex model for")
    print("  equality; the covering floor here is moved by the 2-adic branch atoms.")
    print()

    print("F. Tournament Analysis")
    print("  vertices=proof carriers and cut interfaces, not runners/arcs/raw rows")
    print("  pairwise_observable=retained LRC predicate + overlap-tax payload + scalar-forgetting debt")
    print("  switch_gauge=higher weighted proof-facing score; ties by declared carrier order")
    print(f"  axes={','.join(AXES)}")
    print(f"  score_hist={hist}")
    print(f"  directed_3cycles={directed_3cycles}")
    print(f"  hamiltonian_path_count={hp_count}")
    print("  hamiltonian_path=" + " -> ".join(path))
    print()

    print("G. Assumption challenge")
    print("  Considered vertices: runners, odd blockers, even-half gates, E_safe")
    print("  components, atomic cells, endpoint labels, wall-crossing events, residues,")
    print("  cover arcs, Fourier modes, matroid circuits, graph cuts, and proof")
    print("  obligations.  Chosen vertices: atomic cells plus proof-carrier cuts.")
    print("  This preserves the one-branch LRC predicate exactly because branch0 and")
    print("  overlap_tax are reconstructed from atom multiplicities.  It destroys raw")
    print("  runner order and most continuous geometry.  Challenged assumption: the")
    print("  overlap correction is diffuse harmonic noise; the audit says it has a")
    print("  small labelled cut core in every hard row scanned.")


if __name__ == "__main__":
    main()
