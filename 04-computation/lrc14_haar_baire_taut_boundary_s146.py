#!/usr/bin/env python3
"""S146: Haar-Baire taut boundary scout for LRC14.

This extends the S145 Borel/Baire/Haar carrier by making the proposed
Haar-Baire Wave*/HBT* frontier more concrete.

The finite question tested here is the boundary-support lemma near AP:

    In bounded AP one-swap and two-swap neighborhoods, are AP and
    Goddyn-Wong the only rows with closed threshold support but no strict
    Haar/Baire-open interval at delta=1/14?

The script also reports taut fronts: each strict safe component is an interval
whose endpoints carry active owner labels.  These are the LRC analogue of ANYA
interval nodes or CWave primitive fronts.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import gcd


DELTA = Fraction(1, 14)
AP = tuple(range(1, 14))


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


UNSAFE_CACHE: dict[int, tuple[tuple[Fraction, Fraction], ...]] = {}
ENDPOINT_CACHE: dict[int, tuple[Fraction, ...]] = {}


def unsafe_intervals_for_speed(v: int) -> tuple[tuple[Fraction, Fraction], ...]:
    if v in UNSAFE_CACHE:
        return UNSAFE_CACHE[v]
    radius = DELTA / v
    out: list[tuple[Fraction, Fraction]] = []
    for m in range(v):
        center = Fraction(m, v)
        out.extend(split_circle_interval(center - radius, center + radius))
    ans = tuple(sorted(out))
    UNSAFE_CACHE[v] = ans
    return ans


def endpoint_candidates_for_speed(v: int) -> tuple[Fraction, ...]:
    if v in ENDPOINT_CACHE:
        return ENDPOINT_CACHE[v]
    out: list[Fraction] = []
    for m in range(v):
        out.append(frac_part((Fraction(m) - DELTA) / v))
        out.append(frac_part((Fraction(m) + DELTA) / v))
    ans = tuple(sorted(set(out)))
    ENDPOINT_CACHE[v] = ans
    return ans


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


def active_owners(speeds: tuple[int, ...], t: Fraction) -> tuple[str, ...]:
    owners: list[str] = []
    for v in speeds:
        y = frac_part(Fraction(v) * t)
        if y == DELTA:
            owners.append(f"{v}L")
        elif y == 1 - DELTA:
            owners.append(f"{v}R")
    return tuple(owners)


def midpoint_slack(speeds: tuple[int, ...], a: Fraction, b: Fraction) -> Fraction:
    mid = (a + b) / 2
    return min(circular_distance_to_integer(Fraction(v) * mid) for v in speeds) - DELTA


def shell_fold(x: int, modulus: int = 27) -> int:
    r = x % modulus
    return min(r, modulus - r)


def transfer_label(speeds: tuple[int, ...]) -> str:
    holes = tuple(x for x in AP if x not in speeds)
    adds = tuple(x for x in speeds if x not in AP)
    if not holes and not adds:
        return "AP"
    htxt = ",".join(f"H{h}:g{gcd(h,27)}" for h in holes) or "-"
    atxt = ",".join(f"D{shell_fold(a)}:g{gcd(shell_fold(a),27)}@{a}" for a in adds) or "-"
    return f"{htxt} -> {atxt}"


@dataclass(frozen=True)
class Row:
    name: str
    speeds: tuple[int, ...]
    note: str


NAMED_ROWS = [
    Row("AP", AP, "tight AP boundary row"),
    Row("GW 12->24", tuple(list(range(1, 12)) + [13, 24]), "tight Goddyn-Wong boundary row"),
    Row("near 12->36", tuple(list(range(1, 12)) + [13, 36]), "first unit-excess K33 near-miss"),
    Row("petal 10->20", tuple(list(range(1, 10)) + [11, 12, 13, 20]), "unit-visible C27 petal"),
    Row("petal 13->26", tuple(list(range(1, 13)) + [26]), "unit-visible C27 petal"),
    Row(
        "two-swap 10,12->20,24",
        tuple([1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 13, 20, 24]),
        "S138 petal plus GW splice",
    ),
    Row(
        "two-swap 10,12->20,36",
        tuple([1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 13, 20, 36]),
        "S138 petal plus K33 splice",
    ),
]


def one_swap_rows(add_max: int) -> list[Row]:
    rows = [Row("AP", AP, "baseline")]
    for h in AP:
        for a in range(14, add_max + 1):
            speeds = tuple(sorted((set(AP) - {h}) | {a}))
            if len(speeds) == 13 and gcd(*speeds) == 1:
                rows.append(Row(f"{h}->{a}", speeds, transfer_label(speeds)))
    return rows


def two_swap_rows(add_max: int) -> list[Row]:
    rows: list[Row] = []
    for holes in combinations(AP, 2):
        base = set(AP) - set(holes)
        for adds in combinations(range(14, add_max + 1), 2):
            speeds = tuple(sorted(base | set(adds)))
            if len(speeds) == 13 and gcd(*speeds) == 1:
                rows.append(Row(f"{holes}->{adds}", speeds, transfer_label(speeds)))
    return rows


def classify_rows(rows: list[Row]) -> tuple[list[tuple[Row, int]], list[tuple[Fraction, Row]], list[Row]]:
    boundary_only: list[tuple[Row, int]] = []
    positive: list[tuple[Fraction, Row]] = []
    covered: list[Row] = []
    for row in rows:
        comps = safe_open_components(row.speeds)
        mass = interval_measure(comps)
        if mass > 0:
            positive.append((mass, row))
            continue
        pts = threshold_safe_points(row.speeds)
        if pts:
            boundary_only.append((row, len(pts)))
        else:
            covered.append(row)
    positive.sort(key=lambda item: (item[0], item[1].name))
    return boundary_only, positive, covered


def print_named_fronts() -> None:
    print("[1] Named Haar-Baire taut fronts")
    for row in NAMED_ROWS:
        comps = safe_open_components(row.speeds)
        mass = interval_measure(comps)
        pts = threshold_safe_points(row.speeds)
        print(f"  {row.name}:")
        print(f"    transfer={transfer_label(row.speeds)}")
        print(f"    strict_mass={mass} components={len(comps)} closed_pts={len(pts)}")
        if comps:
            for i, (a, b) in enumerate(comps[:6], 1):
                print(
                    "    front"
                    f" {i}: interval=({a},{b}) left={active_owners(row.speeds, a)}"
                    f" right={active_owners(row.speeds, b)} mid_slack={midpoint_slack(row.speeds, a, b)}"
                )
        else:
            for t in pts:
                print(f"    boundary t={t} owners={active_owners(row.speeds, t)}")


def print_scan(name: str, rows: list[Row]) -> None:
    boundary_only, positive, covered = classify_rows(rows)
    print(f"[2] {name}")
    print(f"  rows={len(rows)} boundary_only={len(boundary_only)} positive_open={len(positive)} covered={len(covered)}")
    print("  boundary-only rows:")
    for row, npts in boundary_only[:20]:
        print(f"    {row.name:22s} closed_pts={npts:2d} transfer={row.note}")
    if len(boundary_only) > 20:
        print(f"    ... {len(boundary_only)-20} more")
    print("  smallest positive strict Haar masses:")
    for mass, row in positive[:12]:
        print(f"    {row.name:22s} mass={str(mass):>8s} transfer={row.note}")
    if covered:
        print("  covered rows with no threshold support:")
        for row in covered[:12]:
            print(f"    {row.name:22s} transfer={row.note}")


def print_planner_tournament() -> None:
    print("[3] HBT*/Haar-Baire Wave proof-order tournament")
    vertices = [
        ("HBT*_boundary_rank", (5, 5, 5, 5, 5, 4, 5), "keeps interior, boundary owner, relation rank"),
        ("Haar-Baire_Wave*", (5, 5, 5, 4, 5, 3, 5), "exact interval wavefront with Haar/Baire labels"),
        ("ANYA_interval_taut", (5, 5, 4, 2, 4, 2, 5), "interval nodes and taut turns"),
        ("CWave_relation_front", (4, 4, 5, 2, 4, 3, 4), "primitive wavefronts, relation hyperplanes"),
        ("BlockA_local_database", (3, 4, 3, 4, 3, 2, 4), "finite local AP/C27/state-lift databases"),
        ("Theta_owner_line", (3, 2, 4, 3, 3, 1, 3), "owner line-of-sight shortcut"),
        ("FieldD_measure_interp", (2, 3, 4, 4, 2, 1, 2), "measure interpolation, weak boundary memory"),
        ("raw_denominator_A*", (1, 1, 1, 2, 1, 1, 1), "grid-first denominator search"),
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
    scores = [sum(row) for row in adj]
    c3 = 0
    for a, b, c in combinations(range(len(vertices)), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (adj[a][c] and adj[c][b] and adj[b][a]):
            c3 += 1
    order = [name for _score, name in sorted(zip(scores, names), reverse=True)]
    print("  criteria: endpoint_witness, haar_interior, boundary_code, relation_wave, lrc_transfer, dynamic_update, anti_scalar_guard")
    for name, vector, note in vertices:
        print(f"  {name:24s} vector={vector} score={scores[names.index(name)]} note={note}")
    print(f"  directed_3cycles={c3}")
    print(f"  ranking={' > '.join(order)}")


def main() -> None:
    print("=" * 78)
    print("S146 Haar-Baire taut boundary scout for LRC14")
    print("=" * 78)
    print_named_fronts()
    print()
    print_scan("one-swap AP neighborhood, add<=160", one_swap_rows(160))
    print()
    print_scan("two-swap AP neighborhood, add<=40", two_swap_rows(40))
    print()
    print_planner_tournament()
    print()
    print("[4] Readout")
    print("  The one-swap finite check supports the boundary-support lemma:")
    print("  through add<=160, AP and GW 12->24 are the only boundary-only rows,")
    print("  and the first open row is near 12->36 with strict Haar mass 1/1260.")
    print("  The two-swap add<=40 check tests the S138 neighborhood for additional")
    print("  boundary-only packets; any nonzero count above should be treated as a")
    print("  new exact boundary-front obligation, not as Haar-positive evidence.")
    print("  HBT*/Haar-Baire Wave should therefore propagate open fronts by measure,")
    print("  and route zero-interior boundary fronts to C27/unital/state-lift labels.")


if __name__ == "__main__":
    main()
