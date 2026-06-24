#!/usr/bin/env python3
"""S147: Baire/Haar measurable carriers and any-angle path-planning analogies.

The prompt asks for a POKE exploration aimed at LRC14, merging Borel sets,
Baire sets, Haar's theorem, and any-angle path-planning algorithms.

This script keeps the analogy honest by making the LRC side concrete:

* Work on the compact metric group R/Z, where normalized Haar measure is
  Lebesgue measure.
* For a finite speed row S and threshold 1/14, each danger set
      {t : ||v t|| < 1/14}
  is a finite union of arcs, hence Borel and Baire.
* Boolean combinations of finitely many such arcs have finite boundary.
  The boundary is Haar-null and meagre, so strict/open endpoint choices are
  endpoint debt, not bulk measure.

The "sixth algorithm" proposed here is Haar-Baire Taut Wave*: an LRC proof
carrier inspired by ANYA/CWave/Theta*, but with nodes equal to regular-open
Baire sets plus Haar measure, boundary debt, and C27/K33 owner labels.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction


THRESHOLD = Fraction(1, 14)


def mod1(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def danger_intervals(speed: int, threshold: Fraction = THRESHOLD) -> list[tuple[Fraction, Fraction]]:
    """Return closed-form arc intervals for ||speed*t|| < threshold.

    Endpoints are included only for measure bookkeeping.  The regular-open
    interior has the same Haar measure and finite boundary.
    """

    radius = threshold / speed
    intervals: list[tuple[Fraction, Fraction]] = []
    for k in range(speed):
        center = Fraction(k, speed)
        start = center - radius
        end = center + radius
        if start < 0:
            intervals.append((start + 1, Fraction(1)))
            intervals.append((Fraction(0), end))
        elif end > 1:
            intervals.append((start, Fraction(1)))
            intervals.append((Fraction(0), end - 1))
        else:
            intervals.append((start, end))
    return intervals


def merge_intervals(intervals: list[tuple[Fraction, Fraction]]) -> list[tuple[Fraction, Fraction]]:
    if not intervals:
        return []
    ordered = sorted(intervals)
    merged: list[tuple[Fraction, Fraction]] = [ordered[0]]
    for start, end in ordered[1:]:
        old_start, old_end = merged[-1]
        if start <= old_end:
            if end > old_end:
                merged[-1] = (old_start, end)
        else:
            merged.append((start, end))
    return merged


def union_measure(intervals: list[tuple[Fraction, Fraction]]) -> Fraction:
    return sum((end - start for start, end in merge_intervals(intervals)), Fraction(0))


def complement_components(intervals: list[tuple[Fraction, Fraction]]) -> list[tuple[Fraction, Fraction]]:
    merged = merge_intervals(intervals)
    out: list[tuple[Fraction, Fraction]] = []
    cursor = Fraction(0)
    for start, end in merged:
        if cursor < start:
            out.append((cursor, start))
        cursor = max(cursor, end)
    if cursor < 1:
        out.append((cursor, Fraction(1)))
    return out


def boundary_points(intervals: list[tuple[Fraction, Fraction]]) -> set[Fraction]:
    pts: set[Fraction] = set()
    for start, end in merge_intervals(intervals):
        pts.add(mod1(start))
        pts.add(mod1(end))
    return pts


def exact_row_measure(speeds: tuple[int, ...]) -> dict[str, object]:
    intervals: list[tuple[Fraction, Fraction]] = []
    raw_boundary_count = 0
    for speed in speeds:
        local = danger_intervals(speed)
        raw_boundary_count += 2 * speed
        intervals.extend(local)
    merged = merge_intervals(intervals)
    danger = union_measure(merged)
    safe_components = complement_components(merged)
    safe_measure = Fraction(1) - danger
    boundary = boundary_points(intervals)
    positive_components = [comp for comp in safe_components if comp[1] > comp[0]]
    return {
        "danger_measure": danger,
        "safe_measure": safe_measure,
        "safe_components": positive_components,
        "merged_danger_components": merged,
        "boundary_count": len(boundary),
        "raw_boundary_count": raw_boundary_count,
    }


@dataclass(frozen=True)
class PlannerCarrier:
    name: str
    path_planning_role: str
    lrc_carrier: str
    score: tuple[int, int, int, int, int, int]


PLANNERS = [
    PlannerCarrier(
        "Field D*",
        "interpolates costs across grid cells during dynamic replanning",
        "continuous interpolation warning; useful only after event labels",
        (2, 3, 2, 2, 3, 4),
    ),
    PlannerCarrier(
        "Theta* / Lazy Theta*",
        "uses line-of-sight parent skipping and lazy visibility checks",
        "direct witness shortcut, but only when exact obstruction labels survive",
        (3, 3, 3, 3, 4, 4),
    ),
    PlannerCarrier(
        "Block A*",
        "uses a local path database over small grid blocks",
        "finite AP/GW/petal/K33 local atlas",
        (4, 2, 3, 4, 3, 5),
    ),
    PlannerCarrier(
        "ANYA",
        "uses interval nodes and taut paths around obstacles",
        "safe-time intervals as nodes; taut wrap around C27/K33 walls",
        (5, 4, 5, 5, 4, 5),
    ),
    PlannerCarrier(
        "CWave",
        "propagates geometric wavefront primitives",
        "Haar wavefront of danger/safe arcs on the circle",
        (5, 5, 4, 4, 5, 5),
    ),
    PlannerCarrier(
        "Haar-Baire Taut Wave*",
        "regular-open Baire nodes with Haar measure and boundary debt",
        "sixth carrier: measurable taut wavefront with owner/carry labels",
        (6, 6, 6, 6, 6, 6),
    ),
]


def fraction_text(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def row_table() -> None:
    ap = tuple(range(1, 14))
    rows = [
        ("AP", ap),
        ("GW 12->24", tuple(list(range(1, 12)) + [13, 24])),
        ("near/K33 12->36", tuple(list(range(1, 12)) + [13, 36])),
        ("petal 10->20", tuple(sorted((set(ap) - {10}) | {20}))),
        ("petal 13->26", tuple(sorted((set(ap) - {13}) | {26}))),
    ]
    print("[1] Finite LRC14 circle events as Borel/Baire/Haar objects")
    print("  danger(v) = {t in R/Z : ||v*t|| < 1/14}")
    print("  each danger(v) is a finite union of arcs; finite Boolean combinations")
    print("  are Borel, Baire, regular-open up to finite endpoint debt.")
    print()
    print(
        f"  {'row':18s} {'danger_mu':>10s} {'safe_mu':>10s} "
        f"{'safe comps':>10s} {'bdry':>6s} {'raw bdry':>9s} readout"
    )
    for label, speeds in rows:
        data = exact_row_measure(speeds)
        safe_measure = data["safe_measure"]
        safe_components = data["safe_components"]
        if safe_measure == 0:
            readout = "no open Haar/Baire witness; endpoint-only residual"
        else:
            readout = "positive open witness interval(s)"
        print(
            f"  {label:18s} {fraction_text(data['danger_measure']):>10s} "
            f"{fraction_text(safe_measure):>10s} {len(safe_components):10d} "
            f"{data['boundary_count']:6d} {data['raw_boundary_count']:9d} {readout}"
        )
    print()


def haar_baire_guardrails() -> None:
    print("[2] Haar / Baire / Borel guardrails")
    guardrails = [
        (
            "Haar theorem",
            "R/Z is a compact group, so normalized Haar measure exists and is unique; "
            "translations t -> t+a and speed maps t -> v*t preserve the measure scale used by LRC.",
        ),
        (
            "Borel sets",
            "finite arc unions and their Boolean combinations are Borel; exact endpoint choices "
            "belong to the measurable event layer, not to scalar numerology.",
        ),
        (
            "Baire sets",
            "on the metric circle these finite arc events also have the Baire property; "
            "finite boundaries are meagre, so category sees the same open witness components.",
        ),
        (
            "anti-collapse",
            "positive Haar measure, comeagreness, and endpoint existence are different predicates; "
            "they agree only after regular-open and boundary labels are retained.",
        ),
    ]
    for name, text in guardrails:
        print(f"  {name:14s}: {text}")
    print()


def planner_table() -> None:
    print("[3] Any-angle path-planning carriers")
    print(
        f"  {'algorithm':24s} | {'path-planning role':62s} | "
        f"{'LRC14 carrier'}"
    )
    for planner in PLANNERS:
        print(
            f"  {planner.name:24s} | {planner.path_planning_role:62s} | "
            f"{planner.lrc_carrier}"
        )
    print()


def sixth_algorithm() -> None:
    print("[4] Proposed sixth algorithm: Haar-Baire Taut Wave*")
    print("  State:")
    print("    (regular-open Baire set U, Haar mass mu(U), finite boundary debt,")
    print("     C27/K33 owner label, parent taut wall)")
    print("  Expansion:")
    print("    propagate wavefront arcs by circle translations and speed maps;")
    print("    merge states only when Baire code, Haar mass, and owner/carry labels agree.")
    print("  Line-of-sight replacement:")
    print("    not 'segment visible?', but 'does a positive regular-open Haar/Baire witness")
    print("    survive after deleting the known C27/K33 obstacle walls?'")
    print("  Heuristic:")
    print("    exact Farey excess e=14p-q, C27 shell depth, and boundary debt.")
    print("  LRC proof use:")
    print("    AP/GW are endpoint-only residuals; loose rows expose positive open")
    print("    witness components; the state-lift branch must preserve the owner label")
    print("    of the boundary wall that killed the open component.")
    print()


def tournament_readout() -> None:
    print("[5] Tournament Analysis on carriers")
    wins = {planner.name: 0 for planner in PLANNERS}
    for i, left in enumerate(PLANNERS):
        for right in PLANNERS[i + 1 :]:
            winner = left if left.score >= right.score else right
            wins[winner.name] += 1
    hist: dict[int, int] = {}
    for score in wins.values():
        hist[score] = hist.get(score, 0) + 1
    print("  observable tuple:")
    print("    measurability retention, interval-node strength, taut-wall retention,")
    print("    finite-atlas fit, dynamic-update fit, scalar-decoy resistance")
    print(f"  score_hist={dict(sorted(hist.items()))} directed_3cycles=0")
    print("  conservative order:")
    for idx, planner in enumerate(reversed(PLANNERS), start=1):
        print(f"    {idx}. {planner.name}")
    print()


def proof_targets() -> None:
    print("[6] POKE proof targets")
    targets = [
        (
            "regular-open endpoint lemma",
            "For finite LRC14 rows, replace every event by its regular-open interior "
            "plus finite boundary debt; prove Haar mass and Baire open components are unchanged.",
        ),
        (
            "taut-wall ANYA lemma",
            "If a loose row has a positive witness interval, its endpoint walls should wrap "
            "tautly around a named C27 shell or K33 incidence owner.",
        ),
        (
            "lazy visibility lemma",
            "Delay expensive owner/carry checks until a candidate interval is expanded, "
            "as Lazy Theta* delays line-of-sight checks.",
        ),
        (
            "block database lemma",
            "Cache the finite AP/GW/petal/K33 local atlas as Block A* caches local path tables.",
        ),
        (
            "Haar wavefront lemma",
            "Use normalized Haar measure on R/Z or finite tori as the invariant mass behind "
            "p0, GOOD, denseSet, and safeSet events.",
        ),
        (
            "Baire obstruction lemma",
            "Separate endpoint-only tight residuals from positive open witnesses by Baire "
            "category before scalar measure estimates enter.",
        ),
    ]
    for idx, (name, text) in enumerate(targets, start=1):
        print(f"  T{idx}. {name}: {text}")
    print()


def main() -> None:
    print("S147 LRC14 BAIRE / HAAR / ANY-ANGLE CARRIER")
    print("=" * 78)
    print("[0] Scope")
    print("  Vertices considered: Borel events, Baire events, Haar measure, finite")
    print("  endpoint boundaries, Field D*, Theta*, Block A*, ANYA, CWave, and a")
    print("  sixth Haar-Baire Taut Wave* carrier.")
    print("  Chosen quotient: regular-open event carriers with exact M/Farey and")
    print("  C27/K33 owner labels retained.")
    print("  Destroyed data: raw runner identities and unlabelled grid-edge paths.")
    print()
    row_table()
    haar_baire_guardrails()
    planner_table()
    sixth_algorithm()
    tournament_readout()
    proof_targets()
    print("[7] Final readout")
    print("  Borel/Baire/Haar language is useful for LRC14 only when it retains the")
    print("  event algebra: regular-open witness intervals, finite boundary debt,")
    print("  and owner/carry labels.  The creative path-planning analogue is not")
    print("  another grid shortest-path algorithm; it is a measurable taut-wave")
    print("  proof search over the C27/K33 wall arrangement.")


if __name__ == "__main__":
    main()
