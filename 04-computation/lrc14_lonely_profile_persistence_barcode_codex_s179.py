#!/usr/bin/env python3
"""Exact lonely-profile persistence barcode scout for LRC14.

This is a proof-interface computation, not a proof of LRC14.

The usual exact safe-component check asks whether

    min_v ||v t|| > 1/14

on an open interval.  This script keeps more geometry: each connected safe
component at threshold 1/14 is treated as a persistence bar of the superlevel
profile m_S(t)=min_v ||v t||, with exact length, peak time, peak height, and
height margin above 1/14.

Tournament Analysis declaration:
  vertices: proof carriers / geometric certificates, not runners;
  pairwise observable: retained LRC predicate, topology, exact height, endpoint
    geometry, stability margin, certificate handoff, and anti-scalar discipline;
  switch/gauge: majority comparison of the observable vector;
  tie Hamiltonian path: the listed carrier order.

Assumption challenge:
  Alternate vertices considered: runners, gaps, fixed circle sections, section
  boundaries, wall-crossing events, residues, cover arcs, Fourier modes,
  matroid circuits, persistence bars, and proof obligations.  This script uses
  persistence bars because they preserve the LRC predicate and safe-component
  topology while exposing exactly what a raw maximin scalar destroys.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations, permutations


N = 14
THRESHOLD = Fraction(1, N)


def fmt(fr: Fraction) -> str:
    return str(fr.numerator) if fr.denominator == 1 else f"{fr.numerator}/{fr.denominator}"


def frac_floor(x: Fraction) -> int:
    return x.numerator // x.denominator


def frac_part(x: Fraction) -> Fraction:
    return x - frac_floor(x)


def dist_to_integer(x: Fraction) -> Fraction:
    y = frac_part(x)
    return min(y, Fraction(1) - y)


def exact_score(row: tuple[int, ...], t: Fraction) -> Fraction:
    return min(dist_to_integer(Fraction(v) * t) for v in row)


def breakpoints(row: tuple[int, ...]) -> list[Fraction]:
    pts = {Fraction(0), Fraction(1)}
    for v in row:
        for k in range(2 * v + 1):
            pts.add(Fraction(k, 2 * v))
    return sorted(pts)


@dataclass(frozen=True)
class Affine:
    slope: int
    intercept: int

    def value(self, t: Fraction) -> Fraction:
        return Fraction(self.slope) * t + self.intercept


def cell_affines(row: tuple[int, ...], left: Fraction, right: Fraction) -> list[Affine]:
    mid = (left + right) / 2
    affines: list[Affine] = []
    for v in row:
        x = Fraction(v) * mid
        base = frac_floor(x)
        residue = x - base
        if residue < Fraction(1, 2):
            affines.append(Affine(v, -base))
        else:
            affines.append(Affine(-v, base + 1))
    return affines


def min_affine_value(affines: list[Affine], t: Fraction) -> Fraction:
    return min(a.value(t) for a in affines)


def max_min_on_segment(affines: list[Affine], left: Fraction, right: Fraction) -> tuple[Fraction, Fraction]:
    candidates = {left, right}
    for a, b in combinations(affines, 2):
        if a.slope == b.slope:
            continue
        t = Fraction(b.intercept - a.intercept, a.slope - b.slope)
        if left <= t <= right:
            candidates.add(t)

    best_t = left
    best_val = Fraction(-1)
    for t in sorted(candidates):
        val = min_affine_value(affines, t)
        if val > best_val or (val == best_val and t < best_t):
            best_val = val
            best_t = t
    return best_val, best_t


def safe_interval_in_cell(
    affines: list[Affine], left: Fraction, right: Fraction, threshold: Fraction
) -> tuple[Fraction, Fraction] | None:
    lo, hi = left, right
    for a in affines:
        # slope is never zero for these distance pieces.
        bound = Fraction(threshold - a.intercept, a.slope)
        if a.slope > 0:
            if bound > lo:
                lo = bound
        else:
            if bound < hi:
                hi = bound
        if lo > hi:
            return None
    if lo == hi:
        return None
    return lo, hi


@dataclass(frozen=True)
class CellPiece:
    left: Fraction
    right: Fraction
    affines: tuple[Affine, ...]


@dataclass
class Component:
    left: Fraction
    right: Fraction
    pieces: list[CellPiece]


@dataclass(frozen=True)
class Bar:
    left: Fraction
    right: Fraction
    length: Fraction
    peak_time: Fraction
    height: Fraction
    persistence: Fraction


def merge_safe_pieces(row: tuple[int, ...], pieces: list[CellPiece]) -> list[Component]:
    if not pieces:
        return []
    pieces = sorted(pieces, key=lambda p: (p.left, p.right))
    components: list[Component] = [
        Component(pieces[0].left, pieces[0].right, [pieces[0]])
    ]
    for piece in pieces[1:]:
        cur = components[-1]
        if piece.left < cur.right:
            cur.right = max(cur.right, piece.right)
            cur.pieces.append(piece)
        elif piece.left == cur.right and exact_score(row, piece.left) > THRESHOLD:
            cur.right = piece.right
            cur.pieces.append(piece)
        else:
            components.append(Component(piece.left, piece.right, [piece]))
    return components


def barcode(row: tuple[int, ...], threshold: Fraction = THRESHOLD) -> tuple[Fraction, Fraction, list[Bar]]:
    row = tuple(sorted(row))
    pts = breakpoints(row)
    safe_pieces: list[CellPiece] = []
    global_height = Fraction(-1)
    global_peak = Fraction(0)

    for left, right in zip(pts, pts[1:]):
        if left == right:
            continue
        affines = cell_affines(row, left, right)
        h, t = max_min_on_segment(affines, left, right)
        if h > global_height or (h == global_height and t < global_peak):
            global_height, global_peak = h, t

        interval = safe_interval_in_cell(affines, left, right, threshold)
        if interval is not None:
            safe_left, safe_right = interval
            # Discard threshold-only slivers; strict topology is handled during
            # merge by checking whether the touching point is strictly safe.
            if safe_left < safe_right:
                safe_pieces.append(CellPiece(safe_left, safe_right, tuple(affines)))

    bars: list[Bar] = []
    for comp in merge_safe_pieces(row, safe_pieces):
        best_h = Fraction(-1)
        best_t = comp.left
        for piece in comp.pieces:
            left = max(comp.left, piece.left)
            right = min(comp.right, piece.right)
            if left > right:
                continue
            h, t = max_min_on_segment(list(piece.affines), left, right)
            if h > best_h or (h == best_h and t < best_t):
                best_h, best_t = h, t
        bars.append(
            Bar(
                left=comp.left,
                right=comp.right,
                length=comp.right - comp.left,
                peak_time=best_t,
                height=best_h,
                persistence=best_h - threshold,
            )
        )

    return global_height, global_peak, bars


def is_fibbinary(x: int) -> bool:
    return x > 0 and (x & (x << 1)) == 0


def is_moser_de_bruijn(x: int) -> bool:
    if x <= 0:
        return False
    y = x
    while y:
        if y % 4 not in (0, 1):
            return False
        y //= 4
    return True


def first_terms(predicate, count: int) -> tuple[int, ...]:
    out: list[int] = []
    x = 1
    while len(out) < count:
        if predicate(x):
            out.append(x)
        x += 1
    return tuple(out)


def replace(row: tuple[int, ...], old: int, new: int) -> tuple[int, ...]:
    return tuple(sorted(new if x == old else x for x in row))


AP13 = tuple(range(1, 14))


ROWS: list[tuple[str, tuple[int, ...]]] = [
    ("AP13", AP13),
    ("GW_12_to_24", replace(AP13, 12, 24)),
    ("K33_12_to_36", replace(AP13, 12, 36)),
    ("petal_10_to_20", replace(AP13, 10, 20)),
    ("petal_13_to_26", replace(AP13, 13, 26)),
    ("covering_12_to_84", replace(AP13, 12, 84)),
    ("fibbinary_first13", first_terms(is_fibbinary, 13)),
    ("moser_de_bruijn_first13", first_terms(is_moser_de_bruijn, 13)),
]


@dataclass(frozen=True)
class RowBarcode:
    name: str
    row: tuple[int, ...]
    M: Fraction
    peak_time: Fraction
    bars: tuple[Bar, ...]

    @property
    def margin(self) -> Fraction:
        return self.M - THRESHOLD

    @property
    def safe_mass(self) -> Fraction:
        return sum((b.length for b in self.bars), Fraction(0))

    @property
    def longest_bar(self) -> Fraction:
        return max((b.length for b in self.bars), default=Fraction(0))

    @property
    def bar_count(self) -> int:
        return len(self.bars)

    @property
    def min_positive_persistence(self) -> Fraction:
        return min((b.persistence for b in self.bars), default=Fraction(0))


def row_barcode(name: str, row: tuple[int, ...]) -> RowBarcode:
    M, peak_time, bars = barcode(row)
    bars = tuple(sorted(bars, key=lambda b: (-b.persistence, -b.length, b.left)))
    return RowBarcode(name, tuple(sorted(row)), M, peak_time, bars)


@dataclass(frozen=True)
class Carrier:
    name: str
    predicate: int
    topology: int
    exact_height: int
    endpoint_geometry: int
    stability: int
    certificate_handoff: int
    anti_scalar: int
    note: str

    @property
    def vector(self) -> tuple[int, ...]:
        return (
            self.predicate,
            self.topology,
            self.exact_height,
            self.endpoint_geometry,
            self.stability,
            self.certificate_handoff,
            self.anti_scalar,
        )


CARRIERS = [
    Carrier(
        "lonely_profile_persistence_barcode",
        3,
        3,
        3,
        2,
        3,
        3,
        3,
        "keeps exact superlevel components, peak heights, and margins",
    ),
    Carrier(
        "exact_safe_component_front",
        3,
        3,
        2,
        3,
        2,
        3,
        3,
        "current Haar/Baire component certificate with endpoint owners",
    ),
    Carrier(
        "fejer_interval_certificate",
        3,
        1,
        2,
        2,
        3,
        3,
        3,
        "dual certificate over a labelled packet and chosen interval",
    ),
    Carrier(
        "ramanujan_exact_period_projector",
        2,
        1,
        1,
        2,
        2,
        3,
        3,
        "period/phase pre-split that still needs endpoint geometry",
    ),
    Carrier(
        "automaton_gap_state",
        2,
        1,
        1,
        1,
        2,
        2,
        3,
        "finite-state side channel for carry/gap packets",
    ),
    Carrier(
        "divisor_abundancy_side_channel",
        2,
        0,
        1,
        1,
        1,
        2,
        3,
        "factorization and unit-excess product guardrail",
    ),
    Carrier(
        "raw_maximin_scalar",
        1,
        0,
        2,
        0,
        1,
        0,
        1,
        "knows M but forgets where and how the witness survives",
    ),
    Carrier(
        "raw_sequence_name",
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        "names a row family but preserves no LRC certificate",
    ),
]


def carrier_edge(a: Carrier, b: Carrier) -> bool:
    va, vb = a.vector, b.vector
    wins = sum(x > y for x, y in zip(va, vb))
    losses = sum(x < y for x, y in zip(va, vb))
    if wins != losses:
        return wins > losses
    return CARRIERS.index(a) < CARRIERS.index(b)


def tournament_matrix(carriers: list[Carrier]) -> list[list[bool]]:
    n = len(carriers)
    mat = [[False] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        if carrier_edge(carriers[i], carriers[j]):
            mat[i][j] = True
        else:
            mat[j][i] = True
    return mat


def directed_3cycles(mat: list[list[bool]]) -> int:
    count = 0
    for a, b, c in combinations(range(len(mat)), 3):
        edges = int(mat[a][b]) + int(mat[b][c]) + int(mat[c][a])
        if edges in (0, 3):
            count += 1
    return count


def scc_sizes(mat: list[list[bool]]) -> list[int]:
    n = len(mat)

    def reach(start: int, graph: list[list[bool]]) -> set[int]:
        seen = {start}
        stack = [start]
        while stack:
            u = stack.pop()
            for v, ok in enumerate(graph[u]):
                if ok and v not in seen:
                    seen.add(v)
                    stack.append(v)
        return seen

    reverse = [[mat[j][i] for j in range(n)] for i in range(n)]
    remaining = set(range(n))
    sizes: list[int] = []
    while remaining:
        seed = min(remaining)
        comp = reach(seed, mat) & reach(seed, reverse)
        sizes.append(len(comp))
        remaining -= comp
    return sizes


def hamiltonian_path_count(mat: list[list[bool]]) -> int:
    count = 0
    for order in permutations(range(len(mat))):
        if all(mat[order[i]][order[i + 1]] for i in range(len(order) - 1)):
            count += 1
    return count


def print_row_report(records: list[RowBarcode]) -> None:
    print("LRC14 lonely-profile persistence barcode scout")
    print("threshold=1/14")
    print()
    print("[ROW BARCODE SUMMARY]")
    for rec in records:
        print(f"- {rec.name}")
        print(f"  row={rec.row}")
        print(
            "  M={M} peak_t={t} margin={margin} bar_count={bc} "
            "safe_mass={mass} longest_bar={longest} min_positive_persistence={minp}".format(
                M=fmt(rec.M),
                t=fmt(rec.peak_time),
                margin=fmt(rec.margin),
                bc=rec.bar_count,
                mass=fmt(rec.safe_mass),
                longest=fmt(rec.longest_bar),
                minp=fmt(rec.min_positive_persistence),
            )
        )
        for idx, bar in enumerate(rec.bars[:3], 1):
            print(
                "    top_bar_{idx}: interval=({left},{right}) length={length} "
                "peak_t={peak_t} height={height} persistence={persistence}".format(
                    idx=idx,
                    left=fmt(bar.left),
                    right=fmt(bar.right),
                    length=fmt(bar.length),
                    peak_t=fmt(bar.peak_time),
                    height=fmt(bar.height),
                    persistence=fmt(bar.persistence),
                )
            )
    print()


def print_barcode_comparisons(records: list[RowBarcode]) -> None:
    zero_bar = [r.name for r in records if r.bar_count == 0]
    positive = [r for r in records if r.bar_count > 0]
    print("[PERSISTENCE READOUT]")
    print(f"zero_bar_rows={zero_bar}")
    print(
        "positive_rows_by_margin="
        + str([(r.name, fmt(r.margin), r.bar_count) for r in sorted(positive, key=lambda x: (-x.margin, x.name))])
    )
    print(
        "positive_rows_by_total_safe_mass="
        + str([(r.name, fmt(r.safe_mass), r.bar_count) for r in sorted(positive, key=lambda x: (-x.safe_mass, x.name))])
    )
    print(
        "smallest_positive_bar="
        + str(
            [
                (r.name, fmt(r.min_positive_persistence), fmt(r.longest_bar))
                for r in sorted(positive, key=lambda x: (x.min_positive_persistence, x.name))[:4]
            ]
        )
    )
    print()
    print("[INTERPRETATION]")
    print(
        "AP/GW are the zero-bar atoms at threshold 1/14.  Positive rows do not merely "
        "have M>1/14; they have exact bars whose height-margin is a stability radius "
        "for local threshold or packet perturbations."
    )
    print(
        "A proof quotient that keeps only M forgets bar count, component length, peak "
        "location, and whether the witness is a wide robust interval or a thin "
        "near-boundary sliver."
    )
    print()


def print_tournament() -> None:
    mat = tournament_matrix(CARRIERS)
    scores = [sum(row) for row in mat]
    score_hist = dict(sorted(Counter(scores).items()))
    print("[TOURNAMENT ANALYSIS]")
    print("vertices=" + str([c.name for c in CARRIERS]))
    print(
        "pairwise_observable=(predicate, topology, exact_height, endpoint_geometry, "
        "stability, certificate_handoff, anti_scalar)"
    )
    print("switch_gauge=majority comparison of retained-coordinate vector")
    print("tie_hamiltonian_path=" + " > ".join(c.name for c in CARRIERS))
    print(f"score_hist={score_hist}")
    print(f"directed_3cycles={directed_3cycles(mat)}")
    print(f"scc_sizes={scc_sizes(mat)}")
    print(f"hamiltonian_path_count={hamiltonian_path_count(mat)}")
    print("score_order=" + " > ".join(CARRIERS[i].name for i in sorted(range(len(CARRIERS)), key=lambda i: -scores[i])))
    print()


def main() -> None:
    records = [row_barcode(name, row) for name, row in ROWS]
    print_row_report(records)
    print_barcode_comparisons(records)
    print_tournament()
    print("[PACKET FIELD PROPOSAL]")
    print("lonely_profile_bar_count")
    print("lonely_profile_total_length")
    print("lonely_profile_longest_bar")
    print("lonely_profile_peak_height")
    print("lonely_profile_peak_time")
    print("lonely_profile_persistence_margin")
    print("lonely_profile_component_signature")
    print()
    print("[ASSUMPTION CHALLENGE]")
    print(
        "Chosen vertices are persistence/proof carriers.  This preserves the LRC "
        "predicate 'there is an open safe bar above 1/14' and destroys runner identity "
        "only after exact scale, topology, and peak geometry are retained."
    )


if __name__ == "__main__":
    main()
