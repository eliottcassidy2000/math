#!/usr/bin/env python3
"""Exact active-bottleneck normal-fan scout for LRC14.

This is a proof-interface computation, not a proof of LRC14.

HYP-3015 kept the lonely-profile persistence barcode of

    m_S(t)=min_v ||v t||.

This script keeps the local normal-fan data behind each bar: which runners
actually constrain the lower envelope at the bar endpoints and at the peak.
Those active bottleneck owners are closer to a checkable certificate than a bar
length or a raw maximin scalar, because they name the local inequalities a
Fejer, endpoint-owner, tope/cocircuit, or state-lift route must discharge.

Tournament Analysis declaration:
  vertices: proof carriers and local certificate schemas, not runners;
  pairwise observable: retained predicate, local bottleneck support, endpoint
    geometry, quotient-purity repair, stability, handoff quality, anti-scalar
    discipline;
  switch/gauge: majority comparison of the observable vector;
  tie Hamiltonian path: the carrier order listed below.

Assumption challenge:
  Alternate vertices considered: runners, gaps, fixed circle sections, section
  boundaries, wall-crossing events, residues, cover arcs, Fourier modes,
  matroid circuits, persistence bars, active bottleneck sets, and proof
  obligations.  This script uses active bottleneck / normal-fan records because
  they preserve the LRC open-bar predicate while adding the local owner support
  destroyed by raw barcodes and automatic-word quotients.
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
    speed: int
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
            affines.append(Affine(v, v, -base))
        else:
            affines.append(Affine(v, -v, base + 1))
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


def candidate_times(row: tuple[int, ...]) -> set[Fraction]:
    row = tuple(sorted(row))
    pts = breakpoints(row)
    candidates = set(pts)
    for left, right in zip(pts, pts[1:]):
        if left == right:
            continue
        affines = cell_affines(row, left, right)
        for a, b in combinations(affines, 2):
            if a.slope == b.slope:
                continue
            t = Fraction(b.intercept - a.intercept, a.slope - b.slope)
            if left <= t <= right:
                candidates.add(t)
    return candidates


def safe_interval_in_cell(
    affines: list[Affine], left: Fraction, right: Fraction, threshold: Fraction
) -> tuple[Fraction, Fraction] | None:
    lo, hi = left, right
    for a in affines:
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
    components: list[Component] = [Component(pieces[0].left, pieces[0].right, [pieces[0]])]
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


def barcode(row: tuple[int, ...], threshold: Fraction = THRESHOLD) -> tuple[Fraction, list[Bar]]:
    row = tuple(sorted(row))
    pts = breakpoints(row)
    safe_pieces: list[CellPiece] = []
    global_height = Fraction(-1)

    for left, right in zip(pts, pts[1:]):
        if left == right:
            continue
        affines = cell_affines(row, left, right)
        h, _ = max_min_on_segment(affines, left, right)
        global_height = max(global_height, h)
        interval = safe_interval_in_cell(affines, left, right, threshold)
        if interval is not None:
            safe_left, safe_right = interval
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
                best_h = h
                best_t = t
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
    return global_height, bars


def owner_side(v: int, t: Fraction) -> str:
    residue = frac_part(Fraction(v) * t)
    if residue < Fraction(1, 2):
        return "+"
    if residue > Fraction(1, 2):
        return "-"
    return "pm"


@dataclass(frozen=True)
class OwnerSet:
    labels: tuple[str, ...]
    speeds: tuple[int, ...]

    @property
    def size(self) -> int:
        return len(self.speeds)

    @property
    def residue_sum(self) -> int:
        return sum(self.speeds) % N

    def short(self) -> str:
        if not self.labels:
            return "-"
        return "{" + ",".join(self.labels) + f"}}[sum={self.residue_sum}]"


def active_owners(row: tuple[int, ...], t: Fraction, score: Fraction | None = None) -> OwnerSet:
    if score is None:
        score = exact_score(row, t)
    labels: list[str] = []
    speeds: list[int] = []
    for v in sorted(row):
        if dist_to_integer(Fraction(v) * t) == score:
            speeds.append(v)
            labels.append(f"{v}{owner_side(v, t)}")
    return OwnerSet(tuple(labels), tuple(speeds))


@dataclass(frozen=True)
class FanBar:
    bar: Bar
    left_owners: OwnerSet
    peak_owners: OwnerSet
    right_owners: OwnerSet

    @property
    def support(self) -> tuple[int, ...]:
        return tuple(sorted(set(self.left_owners.speeds + self.peak_owners.speeds + self.right_owners.speeds)))

    @property
    def support_size(self) -> int:
        return len(self.support)

    @property
    def endpoint_zero_sum(self) -> bool:
        return self.left_owners.residue_sum == 0 or self.right_owners.residue_sum == 0

    def signature(self) -> str:
        return (
            f"L={self.left_owners.short()} "
            f"P={self.peak_owners.short()} "
            f"R={self.right_owners.short()}"
        )


@dataclass(frozen=True)
class RowFan:
    name: str
    row: tuple[int, ...]
    M: Fraction
    global_peak_times: tuple[Fraction, ...]
    global_peak_owner_sets: tuple[OwnerSet, ...]
    bars: tuple[FanBar, ...]

    @property
    def bar_count(self) -> int:
        return len(self.bars)

    @property
    def safe_mass(self) -> Fraction:
        return sum((b.bar.length for b in self.bars), Fraction(0))

    @property
    def min_persistence(self) -> Fraction:
        return min((b.bar.persistence for b in self.bars), default=Fraction(0))

    @property
    def peak_owner_hist(self) -> Counter[int]:
        return Counter(b.peak_owners.size for b in self.bars)

    @property
    def endpoint_owner_hist(self) -> Counter[int]:
        hist: Counter[int] = Counter()
        for b in self.bars:
            hist[b.left_owners.size] += 1
            hist[b.right_owners.size] += 1
        return hist


def row_fan(name: str, row: tuple[int, ...]) -> RowFan:
    row = tuple(sorted(row))
    M, bars = barcode(row)
    peaks = tuple(sorted(t for t in candidate_times(row) if exact_score(row, t) == M))
    peak_owners = tuple(active_owners(row, t, M) for t in peaks)
    fan_bars = []
    for bar in bars:
        fan_bars.append(
            FanBar(
                bar=bar,
                left_owners=active_owners(row, bar.left, THRESHOLD),
                peak_owners=active_owners(row, bar.peak_time, bar.height),
                right_owners=active_owners(row, bar.right, THRESHOLD),
            )
        )
    fan_bars = sorted(
        fan_bars,
        key=lambda fb: (-fb.bar.persistence, -fb.bar.length, fb.peak_owners.labels, fb.bar.left),
    )
    return RowFan(name, row, M, peaks, peak_owners, tuple(fan_bars))


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


def replace_two(row: tuple[int, ...], old_a: int, new_a: int, old_b: int, new_b: int) -> tuple[int, ...]:
    out: list[int] = []
    for x in row:
        if x == old_a:
            out.append(new_a)
        elif x == old_b:
            out.append(new_b)
        else:
            out.append(x)
    return tuple(sorted(out))


AP13 = tuple(range(1, 14))

ROWS: list[tuple[str, tuple[int, ...]]] = [
    ("AP13", AP13),
    ("GW_12_to_24", replace(AP13, 12, 24)),
    ("K33_12_to_36", replace(AP13, 12, 36)),
    ("petal_10_to_20", replace(AP13, 10, 20)),
    ("petal_13_to_26", replace(AP13, 13, 26)),
    ("P10_plus_GW", replace_two(AP13, 10, 20, 12, 24)),
    ("covering_12_to_84", replace(AP13, 12, 84)),
    ("mixed_AP_fiber_12_to_26", replace(AP13, 12, 26)),
    ("mixed_AP_fiber_12_to_96", replace(AP13, 12, 96)),
    ("mixed_GW_fiber_12_to_38", replace(AP13, 12, 38)),
    ("mixed_GW_fiber_12_to_52", replace(AP13, 12, 52)),
    ("fibbinary_first13", first_terms(is_fibbinary, 13)),
    ("moser_de_bruijn_first13", first_terms(is_moser_de_bruijn, 13)),
]


@dataclass(frozen=True)
class Carrier:
    name: str
    predicate: int
    local_bottleneck_support: int
    endpoint_geometry: int
    quotient_purity_repair: int
    stability: int
    handoff_quality: int
    anti_scalar: int
    note: str

    @property
    def vector(self) -> tuple[int, ...]:
        return (
            self.predicate,
            self.local_bottleneck_support,
            self.endpoint_geometry,
            self.quotient_purity_repair,
            self.stability,
            self.handoff_quality,
            self.anti_scalar,
        )


CARRIERS = [
    Carrier(
        "active_bottleneck_normal_fan",
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        "keeps the local active inequalities at endpoints and peaks",
    ),
    Carrier(
        "lonely_profile_persistence_barcode",
        3,
        2,
        2,
        2,
        3,
        3,
        3,
        "keeps open bars, heights, lengths, and persistence margins",
    ),
    Carrier(
        "tope_cocircuit_wall_language",
        3,
        2,
        3,
        2,
        2,
        3,
        3,
        "keeps wall/topology and endpoint cocircuits but not peak support",
    ),
    Carrier(
        "fejer_interval_certificate",
        3,
        2,
        2,
        2,
        3,
        3,
        3,
        "dual certificate after an interval and packet are chosen",
    ),
    Carrier(
        "exact_safe_component_front",
        3,
        1,
        2,
        1,
        2,
        3,
        3,
        "knows strict components but not the bottleneck normal fan",
    ),
    Carrier(
        "automaton_fiber_magnitude_cocycle",
        2,
        1,
        1,
        3,
        2,
        2,
        3,
        "repairs mixed automatic fibers by retaining Farey/magnitude data",
    ),
    Carrier(
        "raw_maximin_scalar",
        1,
        0,
        0,
        1,
        1,
        0,
        1,
        "knows M but forgets local owners and topology",
    ),
    Carrier(
        "raw_automaton_word",
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        "can mix AP/GW boundary atoms with open rows",
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


def print_fan_report(records: list[RowFan]) -> None:
    print("LRC14 active-bottleneck normal-fan scout")
    print("threshold=1/14")
    print()
    print("[ROW NORMAL-FAN SUMMARY]")
    for rec in records:
        peak_owner_preview = "; ".join(os.short() for os in rec.global_peak_owner_sets[:4])
        if len(rec.global_peak_owner_sets) > 4:
            peak_owner_preview += "; ..."
        print(f"- {rec.name}")
        print(f"  row={rec.row}")
        print(
            "  M={M} global_peak_count={gpc} bar_count={bc} safe_mass={mass} "
            "min_persistence={mp}".format(
                M=fmt(rec.M),
                gpc=len(rec.global_peak_times),
                bc=rec.bar_count,
                mass=fmt(rec.safe_mass),
                mp=fmt(rec.min_persistence),
            )
        )
        print(f"  global_peak_owners={peak_owner_preview}")
        print(
            "  bar_peak_owner_hist={ph} endpoint_owner_hist={eh}".format(
                ph=dict(sorted(rec.peak_owner_hist.items())),
                eh=dict(sorted(rec.endpoint_owner_hist.items())),
            )
        )
        for idx, fb in enumerate(rec.bars[:3], 1):
            print(
                "    fan_bar_{idx}: interval=({left},{right}) length={length} "
                "peak_t={peak_t} height={height} persistence={pers} support={support} {sig}".format(
                    idx=idx,
                    left=fmt(fb.bar.left),
                    right=fmt(fb.bar.right),
                    length=fmt(fb.bar.length),
                    peak_t=fmt(fb.bar.peak_time),
                    height=fmt(fb.bar.height),
                    pers=fmt(fb.bar.persistence),
                    support=fb.support,
                    sig=fb.signature(),
                )
            )

    print()
    print("[AGGREGATE NORMAL-FAN FINGERPRINTS]")
    all_peak_hist: Counter[int] = Counter()
    all_endpoint_hist: Counter[int] = Counter()
    support_hist: Counter[int] = Counter()
    endpoint_zero_sum = 0
    total_bars = 0
    for rec in records:
        all_peak_hist.update(rec.peak_owner_hist)
        all_endpoint_hist.update(rec.endpoint_owner_hist)
        for fb in rec.bars:
            total_bars += 1
            support_hist[fb.support_size] += 1
            endpoint_zero_sum += int(fb.endpoint_zero_sum)
    print(f"total_positive_bars={total_bars}")
    print(f"peak_owner_count_hist={dict(sorted(all_peak_hist.items()))}")
    print(f"endpoint_owner_count_hist={dict(sorted(all_endpoint_hist.items()))}")
    print(f"bar_support_size_hist={dict(sorted(support_hist.items()))}")
    print(f"bars_with_endpoint_zero_sum={endpoint_zero_sum}")

    print()
    print("[LOW-PERSISTENCE DISCHARGE QUEUE]")
    queue = sorted(
        (
            (fb.bar.persistence, fb.bar.length, rec.name, fb)
            for rec in records
            for fb in rec.bars
        ),
        key=lambda item: (item[0], item[1], item[2]),
    )[:10]
    for pers, length, name, fb in queue:
        print(
            "  {name}: persistence={pers} length={length} peak={peak} support={support} {sig}".format(
                name=name,
                pers=fmt(pers),
                length=fmt(length),
                peak=fmt(fb.bar.peak_time),
                support=fb.support,
                sig=fb.signature(),
            )
        )

    print()
    print("[MIXED AUTOMATON-FIBER REPAIR READOUT]")
    for name in (
        "AP13",
        "mixed_AP_fiber_12_to_26",
        "mixed_AP_fiber_12_to_96",
        "GW_12_to_24",
        "mixed_GW_fiber_12_to_38",
        "mixed_GW_fiber_12_to_52",
    ):
        rec = next(r for r in records if r.name == name)
        if rec.bars:
            fb = rec.bars[0]
            sig = fb.signature()
            cert = f"open_bar peak={fmt(fb.bar.peak_time)} support={fb.support} {sig}"
        else:
            cert = "zero_bar boundary_peak_owners=" + "; ".join(os.short() for os in rec.global_peak_owner_sets[:4])
        print(f"  {name}: M={fmt(rec.M)} {cert}")


def print_tournament() -> None:
    mat = tournament_matrix(CARRIERS)
    scores = {CARRIERS[i].name: sum(mat[i]) for i in range(len(CARRIERS))}
    print()
    print("[TOURNAMENT ANALYSIS]")
    print("vertices=proof carriers / local certificate schemas, not runners")
    print(
        "observable=(predicate, local_bottleneck_support, endpoint_geometry, "
        "quotient_purity_repair, stability, handoff_quality, anti_scalar)"
    )
    print("switch=majority comparison; tie path=listed carrier order")
    print("score_hist=" + str(dict(sorted(Counter(scores.values()).items()))))
    print(f"directed_3cycles={directed_3cycles(mat)}")
    print(f"scc_sizes={scc_sizes(mat)}")
    print(f"hamiltonian_path_count={hamiltonian_path_count(mat)}")
    score_order = sorted(CARRIERS, key=lambda c: (-scores[c.name], CARRIERS.index(c)))
    print("tie_path=" + " > ".join(c.name for c in CARRIERS))
    print("score_order=" + " > ".join(c.name for c in score_order))
    for carrier in CARRIERS:
        print(f"  {carrier.name}: score={scores[carrier.name]} vector={carrier.vector} note={carrier.note}")


def main() -> None:
    records = [row_fan(name, row) for name, row in ROWS]
    print_fan_report(records)
    print_tournament()


if __name__ == "__main__":
    main()
