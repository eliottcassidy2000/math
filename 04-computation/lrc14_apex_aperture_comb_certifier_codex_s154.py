#!/usr/bin/env python3
"""S154: apex-aperture comb certifier for the remaining LRC14 branch.

The current proof frontier is no longer the apex-majority case
`|M14|>=7`; THM-571 closes it.  The remaining covering-strictness branch has
at most six multiples of 14.  This script tests a local theorem-shaped
certificate for that branch.

Write S = C union 14Q, where C has no multiples of 14.  At a unit apex a/14,
each c in C is closed-safe.  If all denominator-14 boundary contacts move the
same way, one side of a/14 contains a strict C-safe aperture (0,U).  A row is
then strictly lonely if the danger combs of Q do not cover that aperture:

    exists h in (0,U) such that ||q h|| > 1/14 for every q in Q,
    hence t = a/14 +/- h/14 is strict-safe for S.

Rows not certified by this test are not counterexamples.  They are sharper
obligations: every unit apex is first-order balanced/blocked, or all one-sided
apertures are saturated by the few 14-multiple combs.

Tournament Analysis:
  For each row, vertices are the six denominator-14 unit apexes.  The pairwise
  observable is the best uncovered aperture measure after the 14-multiple combs
  are removed.  Continuous aperture data is discretized by exact comparison of
  rational measures; ties follow the unit Hamiltonian path.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import gcd
from pathlib import Path
import argparse
import sys


REPO = Path(__file__).resolve().parents[1]
AP = tuple(range(1, 14))
U14 = (1, 3, 5, 9, 11, 13)
DELTA = Fraction(1, 14)


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    assert spec and spec.loader
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


s150 = load_module(
    "s154_packet_migration",
    REPO / "04-computation" / "lrc14_packet_migration_gauntlet_codex_s150.py",
)


@dataclass(frozen=True)
class Aperture:
    unit: int
    side: str
    upper: Fraction
    uncovered_measure: Fraction
    first_gap: tuple[Fraction, Fraction] | None
    blocker_count: int


@dataclass(frozen=True)
class Cert:
    status: str
    row_name: str
    speeds: tuple[int, ...]
    qdiv: int
    m14_count: int
    unit: int | None
    side: str | None
    aperture_upper: Fraction
    uncovered_measure: Fraction
    witness: Fraction | None
    min_margin: Fraction | None
    tournament_class: str
    tournament_c3: int
    tournament_scores: tuple[int, ...]


def frac_part(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def dist_to_integer(x: Fraction) -> Fraction:
    y = frac_part(x)
    return min(y, 1 - y)


def min_margin(speeds: tuple[int, ...], t: Fraction) -> Fraction:
    return min(dist_to_integer(Fraction(v) * t) for v in speeds)


def merge_intervals(intervals: list[tuple[Fraction, Fraction]]) -> list[tuple[Fraction, Fraction]]:
    if not intervals:
        return []
    intervals = sorted(intervals)
    out = [intervals[0]]
    for lo, hi in intervals[1:]:
        plo, phi = out[-1]
        if lo <= phi:
            if hi > phi:
                out[-1] = (plo, hi)
        else:
            out.append((lo, hi))
    return out


def danger_intervals(qs: tuple[int, ...], upper: Fraction) -> list[tuple[Fraction, Fraction]]:
    intervals: list[tuple[Fraction, Fraction]] = []
    if upper <= 0:
        return intervals
    for q in qs:
        # ||14*q*h|| <= 1/14 iff h is near k/(14q) with radius 1/(196q).
        max_k = int(Fraction(14 * q) * upper) + 2
        radius = Fraction(1, 196 * q)
        denom = 14 * q
        for k in range(max_k + 1):
            center = Fraction(k, denom)
            lo = max(Fraction(0), center - radius)
            hi = min(upper, center + radius)
            if hi > 0 and lo < upper and lo <= hi:
                intervals.append((lo, hi))
    return merge_intervals(intervals)


def complement_gaps(
    intervals: list[tuple[Fraction, Fraction]], upper: Fraction
) -> list[tuple[Fraction, Fraction]]:
    gaps: list[tuple[Fraction, Fraction]] = []
    cur = Fraction(0)
    for lo, hi in intervals:
        if lo > cur:
            gaps.append((cur, lo))
        if hi > cur:
            cur = hi
    if cur < upper:
        gaps.append((cur, upper))
    return [(lo, hi) for lo, hi in gaps if hi > lo]


def aperture_upper(core: tuple[int, ...], unit: int, side: int) -> tuple[Fraction, int]:
    """Return the one-sided strict core aperture in h>0.

    side=+1 means t=a/14+h.  side=-1 means t=a/14-h.
    """
    uppers: list[Fraction] = []
    blockers = 0
    for c in core:
        r = (unit * c) % 14
        if r == 0:
            # core excludes multiples of 14, but keep the function total.
            blockers += 1
            uppers.append(Fraction(0))
        elif side > 0:
            u = Fraction(13 - r, 14 * c)
            if u == 0:
                blockers += 1
            uppers.append(u)
        else:
            u = Fraction(r - 1, 14 * c)
            if u == 0:
                blockers += 1
            uppers.append(u)
    if not uppers:
        return Fraction(0), blockers
    return min(uppers), blockers


def apertures(speeds: tuple[int, ...]) -> list[Aperture]:
    core = tuple(v for v in speeds if v % 14 != 0)
    qs = tuple(v // 14 for v in speeds if v % 14 == 0)
    out: list[Aperture] = []
    for unit in U14:
        for side_int, side_name in ((1, "+"), (-1, "-")):
            upper, blockers = aperture_upper(core, unit, side_int)
            if upper <= 0:
                out.append(Aperture(unit, side_name, upper, Fraction(0), None, blockers))
                continue
            gaps = complement_gaps(danger_intervals(qs, upper), upper)
            uncovered = sum((hi - lo for lo, hi in gaps), Fraction(0))
            out.append(
                Aperture(
                    unit,
                    side_name,
                    upper,
                    uncovered,
                    gaps[0] if gaps else None,
                    blockers,
                )
            )
    return out


def row_name(holes: tuple[int, ...], adds: tuple[int, ...]) -> str:
    return s150.row_name(holes, adds)


def unit_tournament(aps: list[Aperture]) -> tuple[str, int, tuple[int, ...]]:
    by_unit: dict[int, Fraction] = {}
    for unit in U14:
        by_unit[unit] = max(a.uncovered_measure for a in aps if a.unit == unit)
    n = len(U14)
    adj = [[False] * n for _ in range(n)]
    tie = {u: i for i, u in enumerate(U14)}
    for i, u in enumerate(U14):
        for j, v in enumerate(U14):
            if i == j:
                continue
            if by_unit[u] > by_unit[v]:
                adj[i][j] = True
            elif by_unit[u] == by_unit[v] and tie[u] < tie[v]:
                adj[i][j] = True
    # Exact comparison of scalar aperture masses always produces a transitive
    # total preorder, with the fixed unit path breaking ties.  Hence every row
    # lands in the same six-vertex tournament isomorphism class; retaining the
    # class name makes that quotient explicit without spending 6! work per row.
    label = "transitive6:000000000000000"
    c3 = directed_triangles(adj)
    scores = tuple(sorted((sum(1 for j in range(n) if adj[i][j]) for i in range(n)), reverse=True))
    return label, c3, scores


def directed_triangles(adj: list[list[bool]]) -> int:
    n = len(adj)
    total = 0
    for i, j, k in combinations(range(n), 3):
        e = [adj[i][j], adj[j][k], adj[k][i]]
        f = [adj[j][i], adj[k][j], adj[i][k]]
        if all(e) or all(f):
            total += 1
    return total


def certify(
    name: str, speeds: tuple[int, ...], holes: tuple[int, ...] = (), adds: tuple[int, ...] = ()
) -> Cert:
    q = s150.qdiv(speeds)
    m14 = sum(1 for v in speeds if v % 14 == 0)
    aps = apertures(speeds)
    tclass, c3, scores = unit_tournament(aps)

    if q <= 14:
        return Cert("qdiv-direct", name, speeds, q, m14, None, None, Fraction(0), Fraction(0), None, None, tclass, c3, scores)
    if m14 == 0:
        return Cert("no-14-multiple", name, speeds, q, m14, None, None, Fraction(0), Fraction(0), None, None, tclass, c3, scores)
    if m14 >= 7:
        return Cert("THM-571-apex-majority", name, speeds, q, m14, None, None, Fraction(0), Fraction(0), None, None, tclass, c3, scores)

    best = max(aps, key=lambda a: (a.uncovered_measure, a.upper))
    if best.uncovered_measure > 0 and best.first_gap is not None:
        lo, hi = best.first_gap
        h = (lo + hi) / 2
        sign = 1 if best.side == "+" else -1
        t = Fraction(best.unit, 14) + sign * h
        margin = min_margin(speeds, t)
        if margin <= DELTA:
            raise AssertionError((name, best, h, t, margin, speeds))
        return Cert(
            "aperture-comb-certified",
            name,
            speeds,
            q,
            m14,
            best.unit,
            best.side,
            best.upper,
            best.uncovered_measure,
            t,
            margin,
            tclass,
            c3,
            scores,
        )

    positive_apertures = [a for a in aps if a.upper > 0]
    if not positive_apertures:
        status = "all-apertures-first-order-blocked"
    else:
        status = "all-apertures-comb-saturated"
    return Cert(
        status,
        name,
        speeds,
        q,
        m14,
        best.unit if best.upper > 0 else None,
        best.side if best.upper > 0 else None,
        best.upper,
        best.uncovered_measure,
        None,
        None,
        tclass,
        c3,
        scores,
    )


def named_rows() -> list[tuple[str, tuple[int, ...], tuple[int, ...], tuple[int, ...]]]:
    rows: list[tuple[str, tuple[int, ...], tuple[int, ...], tuple[int, ...]]] = [
        ("AP", AP, (), ()),
        ("GW 12->24", tuple(sorted((set(AP) - {12}) | {24})), (12,), (24,)),
        ("covering comb 12->84", tuple(sorted((set(AP) - {12}) | {84})), (12,), (84,)),
        ("covering comb 12->168", tuple(sorted((set(AP) - {12}) | {168})), (12,), (168,)),
        ("unlabelled repair drop(4,6)->add(19,42)", tuple(sorted((set(AP) - {4, 6}) | {19, 42})), (4, 6), (19, 42)),
        ("unlabelled repair drop(2,6)->add(17,42)", tuple(sorted((set(AP) - {2, 6}) | {17, 42})), (2, 6), (17, 42)),
        ("P10+GW", tuple(sorted((set(AP) - {10, 12}) | {20, 24})), (10, 12), (20, 24)),
        ("P10+K33", tuple(sorted((set(AP) - {10, 12}) | {20, 36})), (10, 12), (20, 36)),
    ]
    return rows


def generated_rows(k: int, add_max: int):
    rows, _ = s150.generate_bank(k, add_max)
    for row in rows:
        yield row.name, row.speeds, row.holes, row.adds


def summarize(label: str, certs: list[Cert], sample: int) -> None:
    counts = Counter(c.status for c in certs)
    classes = Counter(c.tournament_class for c in certs)
    score_hist = Counter(c.tournament_scores for c in certs)
    live = [c for c in certs if c.qdiv > 14 and 1 <= c.m14_count <= 6]
    live_counts = Counter(c.status for c in live)
    print(f"[{label}]")
    print(f"  rows_seen={sample}")
    print(f"  status_counts={dict(sorted(counts.items()))}")
    print(f"  live_low_m14_rows={len(live)}")
    print(f"  live_low_m14_status_counts={dict(sorted(live_counts.items()))}")
    print(f"  unit_tournament_iso_classes={len(classes)}")
    print(f"  unit_tournament_score_histograms={dict(score_hist.most_common(6))}")
    residues = [c for c in live if c.status not in {"aperture-comb-certified", "qdiv-direct"}]
    if residues:
        print("  uncertified live rows:")
        for c in residues[:8]:
            print(
                f"    {c.row_name:42s} qdiv={c.qdiv:2d} m14={c.m14_count} "
                f"status={c.status} best={c.unit}{c.side or ''} U={c.aperture_upper}"
            )
    else:
        print("  uncertified live rows: none")
    examples = [c for c in live if c.status == "aperture-comb-certified"]
    if examples:
        print("  certificate examples:")
        for c in examples[:6]:
            print(
                f"    {c.row_name:42s} apex={c.unit}{c.side} "
                f"U={c.aperture_upper} free={c.uncovered_measure} "
                f"t={c.witness} margin={c.min_margin}"
            )
    print()


def print_assumption_challenge() -> None:
    print("S154 LRC14 APEX-APERTURE COMB CERTIFIER")
    print("=" * 78)
    print("[assumption challenge]")
    print("  considered vertices:")
    print("    runners, 14-multiple speeds, residual core speeds, apex units,")
    print("    one-sided apertures, danger-comb teeth, packet families, and proof gates.")
    print("  chosen vertices:")
    print("    six denominator-14 unit apexes for row tournaments; proof gates for")
    print("    the global certificate path.")
    print("  preserved predicate:")
    print("    strict LRC14 witness existence near a denominator-14 apex after")
    print("    THM-571 removes the |M14|>=7 branch.")
    print("  destroyed data:")
    print("    off-apex witnesses far from the denominator-14 shell; those remain")
    print("    for exact Haar/Baire classifiers and state-lift routes.")
    print()


def print_tournament_summary(all_certs: list[Cert]) -> None:
    vertices = [
        "qdiv-direct",
        "THM-571-apex-majority",
        "aperture-comb-certified",
        "all-apertures-first-order-blocked",
        "all-apertures-comb-saturated",
        "state-lift-or-bounded-core-residual",
    ]
    scores = {v: len(vertices) - 1 - i for i, v in enumerate(vertices)}
    print("[Tournament Analysis: proof-gate tournament]")
    print("  vertices=proof gates, not runners")
    print("  observable=which gate preserves strict-counterexample status while")
    print("    adding a constructive witness or sharper residual obligation")
    print("  tie_hamiltonian_path=" + " > ".join(vertices))
    print(f"  score_hist={dict(Counter(scores.values()))}")
    print("  directed_3_cycles=0")
    print("  hamiltonian_paths=1")
    unit_classes = Counter(c.tournament_class for c in all_certs)
    print("[Tournament Analysis: six-unit aperture classes]")
    print(f"  achieved_isomorphism_classes={len(unit_classes)}")
    for cls, count in unit_classes.most_common(6):
        print(f"    class={cls} count={count}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--one-max", type=int, default=420)
    parser.add_argument("--two-max", type=int, default=60)
    parser.add_argument("--three-max", type=int, default=30)
    args = parser.parse_args()

    print_assumption_challenge()
    all_certs: list[Cert] = []

    named = [certify(name, speeds, holes, adds) for name, speeds, holes, adds in named_rows()]
    all_certs.extend(named)
    summarize("named rows", named, len(named))

    for label, k, add_max in [
        ("one-swap AP bank", 1, args.one_max),
        ("two-swap AP bank", 2, args.two_max),
        ("three-swap AP bank", 3, args.three_max),
    ]:
        certs = [certify(name, speeds, holes, adds) for name, speeds, holes, adds in generated_rows(k, add_max)]
        all_certs.extend(certs)
        summarize(f"{label} add<={add_max}", certs, len(certs))

    live = [c for c in all_certs if c.qdiv > 14 and 1 <= c.m14_count <= 6]
    print("[global readout]")
    print(f"  total_rows={len(all_certs)}")
    print(f"  live_low_m14_rows={len(live)}")
    print(f"  live_low_m14_status_counts={dict(sorted(Counter(c.status for c in live).items()))}")
    print("  theorem-shaped necessary condition:")
    print("    A remaining strict counterexample with <=6 multiples of 14 must")
    print("    defeat every denominator-14 one-sided core aperture: each aperture")
    print("    is either first-order AP/GW-balanced or saturated by the few")
    print("    14-multiple danger combs.  Uncertified rows are therefore not")
    print("    generic counterexamples; they are exactly local-apex-obstructed")
    print("    rows needing off-apex Haar/source-kernel or state-lift machinery.")
    print()
    print_tournament_summary(all_certs)


if __name__ == "__main__":
    main()
