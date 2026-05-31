#!/usr/bin/env python3
"""
lrc_analogy_atlas_s450.py

codex-2026-05-31 S450

Build a small exact atlas of Lonely Runner analogues.

The point is not another search for a counterexample.  It is to keep the
equivalent formulations honest by pinning them to the same finite boundary
object: a union of forbidden arcs on the multiplier circle and its protected
endpoint incidence.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import gcd, lcm


def mod1(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def norm_dist(x: Fraction) -> Fraction:
    y = mod1(x)
    return min(y, 1 - y)


def fmt(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def forbidden_intervals(speeds: tuple[int, ...]) -> list[tuple[Fraction, Fraction]]:
    n = len(speeds) + 1
    pieces: list[tuple[Fraction, Fraction]] = []
    for v in speeds:
        half = Fraction(1, n * v)
        for m in range(v):
            center = Fraction(m, v)
            lo = center - half
            hi = center + half
            if lo < 0:
                pieces.append((lo + 1, Fraction(1)))
                pieces.append((Fraction(0), hi))
            elif hi > 1:
                pieces.append((lo, Fraction(1)))
                pieces.append((Fraction(0), hi - 1))
            else:
                pieces.append((lo, hi))
    return pieces


def merge_intervals(pieces: list[tuple[Fraction, Fraction]]) -> list[tuple[Fraction, Fraction]]:
    if not pieces:
        return []
    pieces = sorted(pieces)
    merged: list[list[Fraction]] = [[pieces[0][0], pieces[0][1]]]
    for lo, hi in pieces[1:]:
        cur = merged[-1]
        if lo <= cur[1]:
            if hi > cur[1]:
                cur[1] = hi
        else:
            merged.append([lo, hi])
    return [(lo, hi) for lo, hi in merged]


def endpoints(speeds: tuple[int, ...]) -> set[Fraction]:
    n = len(speeds) + 1
    out: set[Fraction] = set()
    for v in speeds:
        den = n * v
        for m in range(v):
            for eps in (-1, 1):
                out.add(mod1(Fraction(n * m + eps, den)))
    return out


def protected_endpoint(e: Fraction, speeds: tuple[int, ...]) -> bool:
    n = len(speeds) + 1
    threshold = Fraction(1, n)
    return any(norm_dist(v * e) < threshold for v in speeds)


def complement_gaps(merged: list[tuple[Fraction, Fraction]]) -> list[tuple[Fraction, Fraction]]:
    if not merged:
        return [(Fraction(0), Fraction(1))]
    gaps: list[tuple[Fraction, Fraction]] = []
    for (_, hi), (lo2, _) in zip(merged, merged[1:]):
        if lo2 > hi:
            gaps.append((hi, lo2))
    first_lo = merged[0][0]
    last_hi = merged[-1][1]
    if first_lo > 0 or last_hi < 1:
        gaps.append((last_hi, first_lo + 1))
    return gaps


@dataclass(frozen=True)
class LrcBoundarySummary:
    label: str
    n: int
    speeds: tuple[int, ...]
    quotient: int
    forbidden_measure: Fraction
    status: str
    endpoint_count: int
    unprotected_count: int
    max_gap: Fraction
    first_unprotected: Fraction | None


def summarize(label: str, speeds: tuple[int, ...]) -> LrcBoundarySummary:
    n = len(speeds) + 1
    speed_lcm = 1
    for v in speeds:
        speed_lcm = lcm(speed_lcm, v)
    q = n * speed_lcm

    merged = merge_intervals(forbidden_intervals(speeds))
    forbidden_measure = sum(hi - lo for lo, hi in merged)
    eps = endpoints(speeds)
    unprotected = sorted(e for e in eps if not protected_endpoint(e, speeds))
    gaps = complement_gaps(merged)
    max_gap = max((hi - lo for lo, hi in gaps), default=Fraction(0))
    if forbidden_measure < 1:
        status = "positive_gap"
    elif unprotected:
        status = "boundary_only"
    else:
        status = "full_open_cover_candidate"
    return LrcBoundarySummary(
        label=label,
        n=n,
        speeds=tuple(sorted(speeds)),
        quotient=q,
        forbidden_measure=forbidden_measure,
        status=status,
        endpoint_count=len(eps),
        unprotected_count=len(unprotected),
        max_gap=max_gap,
        first_unprotected=unprotected[0] if unprotected else None,
    )


def initial(n: int) -> tuple[int, ...]:
    return tuple(range(1, n))


def print_boundary_table() -> None:
    examples = [
        ("initial n=8", initial(8)),
        ("initial n=14", initial(14)),
        ("n14 seven-ladder", (1, 7, 14, 21, 28, 35, 49, 56, 63, 70, 77, 84, 91)),
        ("n14 S380 gate ladder", (1, 14, 28, 42, 56, 70, 98, 112, 126, 140, 154, 168, 182)),
        ("initial n=16", initial(16)),
    ]
    print("Exact finite-boundary checks")
    print("=" * 88)
    print(
        f"{'example':<22} {'n':>3} {'status':<26} {'forbidden':>12} "
        f"{'max_gap/th':>12} {'unprot':>8} {'first_unprot':>14} {'Q':>14}"
    )
    for label, speeds in examples:
        row = summarize(label, speeds)
        threshold = Fraction(1, row.n)
        first = "-" if row.first_unprotected is None else fmt(row.first_unprotected)
        print(
            f"{row.label:<22} {row.n:>3} {row.status:<26} "
            f"{fmt(row.forbidden_measure):>12} {fmt(row.max_gap / threshold):>12} "
            f"{row.unprotected_count:>8} {first:>14} {row.quotient:>14}"
        )


ANALOGUES = [
    (
        "endpoint arc cover",
        "exact equivalence",
        "forbidden intervals in R/Z",
        "unprotected endpoint or positive complement",
        "THM-357 finite protected-boundary hypergraph",
    ),
    (
        "distance graph coloring",
        "exact equivalence",
        "multiplier circular coloring of G(Z,D)",
        "regular (n)-coloring by x -> t*x",
        "same endpoints are failed tight colorings",
    ),
    (
        "subtorus/cube avoidance",
        "exact equivalence",
        "one-dimensional subgroup t*v in (R/Z)^k",
        "hit the central safe cube",
        "same orbit, different projection",
    ),
    (
        "view obstruction",
        "exact equivalence family",
        "ray through periodic obstacle field",
        "line of sight avoiding obstacles",
        "forbidden arcs are shadow intervals of obstacles",
    ),
    (
        "zonotope covering radius",
        "exact equivalence family",
        "projected cube as lattice zonotope",
        "covering-radius inequality",
        "dualizes endpoint debt into lattice covering debt",
    ),
    (
        "lonely runner spectra",
        "same moduli space",
        "1D subtori and L-infinity distance to center",
        "spectrum bound at dimension k",
        "recursion n -> n-1 mirrors endpoint deletion/runner deletion",
    ),
    (
        "nowhere-zero flow compression",
        "application, not isomorphism",
        "flow values on graph edges",
        "compress values to {1,...,k}",
        "shares modular avoidance but not the same boundary object",
    ),
    (
        "Danzer/dense forest",
        "analogy",
        "line/ray hits periodic or sparse obstacles",
        "unobstructed line of sight",
        "same visual metaphor, different quantifiers and dimension",
    ),
]


def print_analogue_table() -> None:
    print()
    print("Analogue atlas")
    print("=" * 88)
    print(f"{'analogue':<30} {'relation':<24} {'object':<34} core witness")
    for name, relation, obj, witness, _ in ANALOGUES:
        print(f"{name:<30} {relation:<24} {obj:<34} {witness}")
    print()
    print("Repo-relevant structure")
    print("=" * 88)
    for name, _, _, _, lesson in ANALOGUES:
        print(f"- {name}: {lesson}")


def print_synthesis() -> None:
    print()
    print("S450 synthesis")
    print("=" * 88)
    print(
        "The invariant shared by the exact formulations is not the speed set "
        "itself.  It is a one-parameter subgroup plus a finite protected "
        "boundary of forbidden neighborhoods."
    )
    print(
        "A counterexample has the same abstract shape in every exact model: "
        "an open cover with no exposed boundary.  The repo's endpoint debt, "
        "zonotope covering debt, and IP row debt are three names for this "
        "same missing boundary witness."
    )
    print(
        "The genuinely useful analogues are those that preserve boundary "
        "incidence.  Analogues that preserve only density, volume, or a "
        "visual obstruction are suggestive but lose the thing LRC seems to "
        "care about."
    )


def main() -> None:
    print_boundary_table()
    print_analogue_table()
    print_synthesis()


if __name__ == "__main__":
    main()
