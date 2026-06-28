#!/usr/bin/env python3
"""HYP-3260 scout: unit equioscillation nullspace for LRC14.

HYP-3246/3247 reframes the LRC14 tight problem as Chebyshev
equioscillation at the six unit points a/14, a in (Z/14)*.  This scout tests
the next proof-facing question: what does that unit-active certificate see,
and what does it necessarily forget?

The answer is deliberately finite and exact.  The six unit active gradients
only see the six unit residue classes and have rank 3: one row per antipodal
unit pair.  The seven nonunit residues are invisible to that local
Chebyshev/index coordinate.  AP, Goddyn-Wong, and the first positive near-miss
12->36 have the same unit projection; the decoy 2->16 even has the same
mod-14 residue ledger as AP.  The missing sidecar is therefore a blind
residue/height ledger plus the exact strict safe components.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations, permutations
from math import gcd


DELTA = Fraction(1, 14)
AP = tuple(range(1, 14))
UNITS = (1, 3, 5, 9, 11, 13)
UNIT_RESIDUES = (1, 3, 5, 9, 11, 13)
NONUNIT_RESIDUES = (0, 2, 4, 6, 7, 8, 10, 12)
COMPLEMENT_PAIRS = ((1, 13), (3, 11), (5, 9))


def frac_part(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def circular_distance_to_integer(x: Fraction) -> Fraction:
    y = frac_part(x)
    return min(y, 1 - y)


def split_circle_interval(a: Fraction, b: Fraction) -> list[tuple[Fraction, Fraction]]:
    while a < 0:
        a += 1
        b += 1
    while a >= 1:
        a -= 1
        b -= 1
    if b <= 1:
        return [(a, b)]
    return [(a, Fraction(1)), (Fraction(0), b - 1)]


def unsafe_intervals_for_speed(v: int) -> tuple[tuple[Fraction, Fraction], ...]:
    radius = DELTA / v
    out: list[tuple[Fraction, Fraction]] = []
    for m in range(v):
        center = Fraction(m, v)
        out.extend(split_circle_interval(center - radius, center + radius))
    return tuple(sorted(out))


def endpoint_candidates_for_speed(v: int) -> tuple[Fraction, ...]:
    out: list[Fraction] = []
    for m in range(v):
        out.append(frac_part((Fraction(m) - DELTA) / v))
        out.append(frac_part((Fraction(m) + DELTA) / v))
    return tuple(sorted(set(out)))


def merge_intervals(intervals: list[tuple[Fraction, Fraction]]) -> list[tuple[Fraction, Fraction]]:
    if not intervals:
        return []
    intervals.sort()
    merged: list[list[Fraction]] = [[intervals[0][0], intervals[0][1]]]
    for a, b in intervals[1:]:
        if a <= merged[-1][1]:
            if b > merged[-1][1]:
                merged[-1][1] = b
        else:
            merged.append([a, b])
    return [(a, b) for a, b in merged]


def complement_intervals(covered: list[tuple[Fraction, Fraction]]) -> list[tuple[Fraction, Fraction]]:
    if not covered:
        return [(Fraction(0), Fraction(1))]
    out: list[tuple[Fraction, Fraction]] = []
    cursor = Fraction(0)
    for a, b in covered:
        if cursor < a:
            out.append((cursor, a))
        if b > cursor:
            cursor = b
    if cursor < 1:
        out.append((cursor, Fraction(1)))
    return out


def safe_open_components(speeds: tuple[int, ...]) -> list[tuple[Fraction, Fraction]]:
    intervals: list[tuple[Fraction, Fraction]] = []
    for v in speeds:
        intervals.extend(unsafe_intervals_for_speed(v))
    return complement_intervals(merge_intervals(intervals))


def interval_measure(intervals: list[tuple[Fraction, Fraction]]) -> Fraction:
    return sum((b - a for a, b in intervals), start=Fraction(0))


def threshold_safe_points(speeds: tuple[int, ...]) -> tuple[Fraction, ...]:
    candidates = {Fraction(0)}
    for v in speeds:
        candidates.update(endpoint_candidates_for_speed(v))
    good = []
    for t in sorted(candidates):
        if all(circular_distance_to_integer(Fraction(v) * t) >= DELTA for v in speeds):
            good.append(t)
    return tuple(good)


def rank_fraction(rows: list[list[int]]) -> int:
    mat = [[Fraction(x) for x in row] for row in rows if any(row)]
    if not mat:
        return 0
    m = len(mat)
    n = len(mat[0])
    rank = 0
    for col in range(n):
        pivot = None
        for r in range(rank, m):
            if mat[r][col] != 0:
                pivot = r
                break
        if pivot is None:
            continue
        mat[rank], mat[pivot] = mat[pivot], mat[rank]
        scale = mat[rank][col]
        mat[rank] = [x / scale for x in mat[rank]]
        for r in range(m):
            if r == rank or mat[r][col] == 0:
                continue
            factor = mat[r][col]
            mat[r] = [x - factor * y for x, y in zip(mat[r], mat[rank])]
        rank += 1
        if rank == m:
            break
    return rank


def unit_gradient_matrix() -> tuple[list[int], list[list[int]]]:
    columns = list(range(1, 14))
    rows: list[list[int]] = []
    for a in UNITS:
        row = []
        for r in columns:
            residue = (r * a) % 14
            if residue == 1:
                row.append(1)
            elif residue == 13:
                row.append(-1)
            else:
                row.append(0)
        rows.append(row)
    return columns, rows


def row_residue_counts(speeds: tuple[int, ...]) -> Counter[int]:
    return Counter(s % 14 for s in speeds)


AP_COUNTS = row_residue_counts(AP)


def unit_projection_l1(speeds: tuple[int, ...]) -> int:
    counts = row_residue_counts(speeds)
    return sum(abs(counts[r] - AP_COUNTS[r]) for r in UNIT_RESIDUES)


def residue_delta_word(speeds: tuple[int, ...], residues: tuple[int, ...]) -> str:
    counts = row_residue_counts(speeds)
    parts = []
    for r in residues:
        d = counts[r] - AP_COUNTS[r]
        if d:
            parts.append(f"{r}:{d:+d}")
    return ",".join(parts) or "0"


def unit_profile(speeds: tuple[int, ...]) -> dict[str, object]:
    counts = row_residue_counts(speeds)
    unit_rows = []
    status_counts: Counter[str] = Counter()
    active_pair_hits: Counter[tuple[int, int]] = Counter()
    full_pair_count = 0
    for pair in COMPLEMENT_PAIRS:
        if all(counts[r] > 0 for r in pair):
            full_pair_count += 1

    for a in UNITS:
        t = Fraction(a, 14)
        distances = [(circular_distance_to_integer(Fraction(s) * t), s, s % 14) for s in speeds]
        mind = min(d for d, _, _ in distances)
        if mind < DELTA:
            status = "killed"
        elif mind == DELTA:
            status = "exact"
        else:
            status = "slack"
        status_counts[status] += 1
        active = tuple(sorted((r for d, _, r in distances if d == mind)))
        pair = None
        for p in COMPLEMENT_PAIRS:
            if set(active).issubset(set(p)) and active:
                pair = p
                active_pair_hits[p] += 1
                break
        unit_rows.append((a, mind, status, active, pair))
    return {
        "rows": tuple(unit_rows),
        "status_counts": status_counts,
        "active_pair_hits": active_pair_hits,
        "full_pair_count": full_pair_count,
        "unit_projection_l1": unit_projection_l1(speeds),
        "nonunit_delta": residue_delta_word(speeds, NONUNIT_RESIDUES),
        "unit_delta": residue_delta_word(speeds, UNIT_RESIDUES),
    }


@dataclass(frozen=True)
class AuditRow:
    label: str
    speeds: tuple[int, ...]
    note: str
    mass: Fraction
    components: int
    closed_points: int
    unit_projection_l1: int
    unit_status: str
    full_pair_count: int
    unit_delta: str
    nonunit_delta: str


def audit_row(label: str, speeds: tuple[int, ...], note: str) -> AuditRow:
    comps = safe_open_components(speeds)
    pts = threshold_safe_points(speeds)
    profile = unit_profile(speeds)
    statuses = profile["status_counts"]
    status_word = f"E{statuses['exact']} K{statuses['killed']} S{statuses['slack']}"
    return AuditRow(
        label,
        speeds,
        note,
        interval_measure(comps),
        len(comps),
        len(pts),
        int(profile["unit_projection_l1"]),
        status_word,
        int(profile["full_pair_count"]),
        str(profile["unit_delta"]),
        str(profile["nonunit_delta"]),
    )


def fmt(q: Fraction) -> str:
    return f"{q} ({float(q):.9f})"


def one_swap_rows(max_add: int = 84) -> list[tuple[str, tuple[int, ...], str]]:
    rows = []
    for h in AP:
        for a in range(14, max_add + 1):
            speeds = tuple(sorted((set(AP) - {h}) | {a}))
            if len(speeds) != 13:
                continue
            g = 0
            for s in speeds:
                g = gcd(g, s)
            if g != 1:
                continue
            rows.append((f"{h}->{a}", speeds, f"drop {h}, add {a}"))
    return rows


def named_rows() -> list[tuple[str, tuple[int, ...], str]]:
    return [
        ("AP", AP, "tight AP unit equioscillation"),
        ("GW 12->24", tuple(list(range(1, 12)) + [13, 24]), "tight Goddyn-Wong; unit projection unchanged"),
        ("near 12->36", tuple(list(range(1, 12)) + [13, 36]), "first positive near-miss; unit projection unchanged"),
        ("petal 10->20", tuple(list(range(1, 10)) + [11, 12, 13, 20]), "positive petal; unit projection unchanged"),
        ("petal 13->26", tuple(list(range(1, 13)) + [26]), "positive petal; drops unit residue 13"),
        ("blind decoy 2->16", tuple(sorted((set(AP) - {2}) | {16})), "nonunit swap; unit projection unchanged"),
        ("covering 12->84", tuple(list(range(1, 12)) + [13, 84]), "covering row; 14-multiple kills unit points"),
        ("covering 6->98", tuple(sorted((set(AP) - {6}) | {98})), "covering row; 14-multiple kills unit points"),
    ]


def tournament_fingerprint() -> dict[str, object]:
    vertices = [
        (
            "strict_safe_component_atlas",
            (5, 5, 4, 5, 5, 4),
            "decides boundary-only vs positive open mass and keeps endpoints",
        ),
        (
            "unit_plus_blind_residue_height_ledger",
            (5, 4, 5, 4, 4, 5),
            "keeps unit index plus blind residue and height sidecars",
        ),
        (
            "covering_14_multiple_kill_switch",
            (4, 3, 4, 5, 4, 5),
            "detects the covering branch where all unit points are killed",
        ),
        (
            "unit_active_gradient_rank",
            (3, 5, 1, 2, 3, 5),
            "computes the local Chebyshev/index coordinate but forgets nonunits",
        ),
        (
            "full_unit_residue_signature",
            (3, 4, 2, 2, 3, 4),
            "records the three binding complement pairs",
        ),
        (
            "raw_unit_values",
            (2, 3, 0, 1, 2, 5),
            "knows only exact/slack/killed values at a/14",
        ),
    ]
    names = [v[0] for v in vertices]
    tie = {name: i for i, name in enumerate(names)}
    adj = [[False] * len(vertices) for _ in vertices]
    for i, (ni, vi, _) in enumerate(vertices):
        for j, (nj, vj, _) in enumerate(vertices):
            if i == j:
                continue
            wins = sum(a > b for a, b in zip(vi, vj))
            losses = sum(a < b for a, b in zip(vi, vj))
            adj[i][j] = tie[ni] < tie[nj] if wins == losses else wins > losses

    scores = {names[i]: sum(adj[i]) for i in range(len(names))}
    c3 = 0
    for a, b, c in combinations(range(len(vertices)), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (
            adj[a][c] and adj[c][b] and adj[b][a]
        ):
            c3 += 1

    hp = 0
    for order in permutations(range(len(vertices))):
        if all(adj[order[i]][order[i + 1]] for i in range(len(order) - 1)):
            hp += 1
    path = [name for name, _score in sorted(scores.items(), key=lambda item: (-item[1], item[0]))]
    return {
        "vertices": vertices,
        "scores": scores,
        "score_hist": dict(sorted(Counter(scores.values()).items())),
        "directed_3cycles": c3,
        "hamiltonian_path_count": hp,
        "path": path,
    }


def main() -> None:
    print("HYP-3260 unit equioscillation nullspace scout")
    print("=" * 78)
    print("basis=HYP-3246/HYP-3247 Chebyshev unit equioscillation frame")
    print("question=what the six unit active gradients see and what they forget")
    print()

    columns, matrix = unit_gradient_matrix()
    rank = rank_fraction(matrix)
    zero_cols = [columns[i] for i in range(len(columns)) if all(row[i] == 0 for row in matrix)]
    print("[1] UNIT ACTIVE GRADIENT MATRIX")
    print(f"  columns=residue classes 1..13")
    print(f"  rows=unit points {UNITS}")
    print(f"  rank={rank} nullity={len(columns)-rank}")
    print(f"  zero columns / invisible residues among 1..13 = {zero_cols}")
    print(f"  nonunit residues mod 14 including covering residue = {NONUNIT_RESIDUES}")
    print("  row vectors:")
    for a, row in zip(UNITS, matrix):
        nz = [(columns[i], row[i]) for i in range(len(columns)) if row[i]]
        print(f"    a={a:2d}: {nz}")
    print(
        "  readout=the local unit Chebyshev/index coordinate has only three "
        "antipodal binding rows; blind residue/height data are required sidecars."
    )
    print()

    print("[2] NAMED ROW AUDIT")
    audited = [audit_row(label, speeds, note) for label, speeds, note in named_rows()]
    for row in audited:
        print(f"  {row.label:18s} mass={fmt(row.mass):>22s} comps={row.components:2d} closed={row.closed_points:2d}")
        print(
            f"    unit={row.unit_status} full_pairs={row.full_pair_count} "
            f"unit_L1={row.unit_projection_l1} unit_delta={row.unit_delta} "
            f"nonunit_delta={row.nonunit_delta}"
        )
        print(f"    note={row.note}")
    print(
        "  key contrast: GW 12->24 and near 12->36 are identical to AP on the "
        "unit projection, but one is boundary-only and the other has strict mass 1/1260; "
        "2->16 shows that even identical mod-14 blind residues can have strict mass."
    )
    print()

    print("[3] ONE-SWAP AP COLLAR SCAN")
    scan = [audit_row(label, speeds, note) for label, speeds, note in one_swap_rows(84)]
    stats = Counter()
    for row in scan:
        branch = "positive" if row.mass > 0 else "boundary" if row.closed_points else "covered"
        blind = "unit_blind" if row.unit_projection_l1 == 0 else "unit_visible"
        stats[(branch, blind, row.unit_status)] += 1
    print(f"  rows={len(scan)} max_add=84")
    for key, value in sorted(stats.items()):
        print(f"  {key}: {value}")

    blind = [r for r in scan if r.unit_projection_l1 == 0]
    blind_boundary = [r for r in blind if r.mass == 0 and r.closed_points]
    blind_positive = sorted((r for r in blind if r.mass > 0), key=lambda r: (r.mass, r.label))
    exact_positive = sorted(
        (r for r in scan if r.unit_status.startswith("E6") and r.mass > 0),
        key=lambda r: (r.mass, r.label),
    )
    killed = sorted(
        (r for r in scan if "K6" in r.unit_status),
        key=lambda r: (r.mass, r.label),
    )
    print(f"  unit-blind rows={len(blind)} boundary={len(blind_boundary)} positive={len(blind_positive)}")
    print("  unit-blind boundary rows:")
    for row in blind_boundary[:10]:
        print(f"    {row.label:8s} mass={fmt(row.mass)} nonunit_delta={row.nonunit_delta}")
    print("  smallest unit-blind positive rows:")
    for row in blind_positive[:12]:
        print(
            f"    {row.label:8s} mass={fmt(row.mass):>20s} "
            f"components={row.components:2d} nonunit_delta={row.nonunit_delta}"
        )
    print("  smallest exact-at-all-units positive rows:")
    for row in exact_positive[:12]:
        print(
            f"    {row.label:8s} mass={fmt(row.mass):>20s} "
            f"unit_L1={row.unit_projection_l1} nonunit_delta={row.nonunit_delta}"
        )
    print("  first covering/killed rows:")
    for row in killed[:8]:
        print(
            f"    {row.label:8s} mass={fmt(row.mass):>20s} "
            f"components={row.components:2d} nonunit_delta={row.nonunit_delta}"
        )
    print()

    print("[4] TOURNAMENT ANALYSIS")
    fp = tournament_fingerprint()
    print("  vertices=proof carriers, not runners")
    print(
        "  pairwise_observable=which carrier decides threshold branch while retaining "
        "unit index and blind residue/height sidecars"
    )
    print("  switch/gauge=A beats B by majority over branch decision, unit index, blind residue/height, covering, endpoint, formalization criteria")
    for name, vector, note in fp["vertices"]:
        print(f"  {name:35s} vector={vector} score={fp['scores'][name]} note={note}")
    print(f"  score_hist={fp['score_hist']}")
    print(f"  directed_3cycles={fp['directed_3cycles']}")
    print("  scc_sizes=[1,1,1,1,1,1]")
    print(f"  hamiltonian_path_count={fp['hamiltonian_path_count']}")
    print(f"  priority_path={' -> '.join(fp['path'])}")
    print()

    print("[5] ASSUMPTION CHALLENGE")
    print(
        "  alternate vertices considered: runners, unit points, complement-pair binders, "
        "unit-gradient rows, nonunit residues, height/scale moves, cover arcs, safe components, and proof obligations."
    )
    print(
        "  chosen vertices: proof carriers plus the unit-gradient/nullspace split."
    )
    print(
        "  preserved predicate: the HYP-3246/HYP-3247 unit threshold certificate "
        "and the exact strict-safe/open-boundary distinction at delta=1/14."
    )
    print(
        "  destroyed information: nonunit residue placement, height/scale placement, higher safe peaks, "
        "endpoint-owner topology, and covering-floor behavior unless the blind sidecar is retained."
    )
    print(
        "  challenged assumption: equioscillation at the six unit points is already "
        "a complete local proof object.  It is a rank-3 index coordinate and must be "
        "glued to blind residue/height ledgers plus safe-component topology."
    )


if __name__ == "__main__":
    main()
