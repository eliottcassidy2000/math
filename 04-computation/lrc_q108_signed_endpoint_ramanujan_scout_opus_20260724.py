#!/usr/bin/env python3
"""Exact signed-endpoint scout for the reduced OPEN-Q-108 target.

For an 11-speed core C and a new speed W, theta=1/14, write G_C as a
disjoint union of closed intervals [l,r].  If D is the centered danger arc
of density 1/7 and H is the continuous periodic primitive of 1_D-1/7, then

    eps_W(G_C) = |G_C intersect D_W|-|G_C|/7
               = (1/W) sum_[l,r] (H(Wr)-H(Wl)).

This script evaluates that identity in Fraction arithmetic, retains every
active endpoint owner, and proves the finite part of THM-2162's exact
one-swap stability theorem.  An optional scout mode ranks a bounded bank of
configurations.  Neither mode proves the universal OPEN-Q-108 inequality.
"""

from __future__ import annotations

import argparse
import itertools
from dataclasses import dataclass
from fractions import Fraction as F
from math import gcd
from functools import reduce

THETA = F(1, 14)
P = F(1, 7)
TARGET = F(7, 858)
LOCAL_NEXT = F(3859, 420420)
EXTREMIZER = tuple(v for v in range(1, 14) if v != 6)


@dataclass(frozen=True)
class EndpointTerm:
    endpoint: F
    orientation: int
    owners: tuple[int, ...]
    phase: F
    value: F


def safe_intervals(v: int) -> list[tuple[F, F]]:
    return [
        (F(14 * k + 1, 14 * v), F(14 * k + 13, 14 * v))
        for k in range(v)
    ]


def intersect(a: list[tuple[F, F]], b: list[tuple[F, F]]) -> list[tuple[F, F]]:
    out: list[tuple[F, F]] = []
    i = j = 0
    while i < len(a) and j < len(b):
        lo = max(a[i][0], b[j][0])
        hi = min(a[i][1], b[j][1])
        if lo <= hi:
            if out and lo <= out[-1][1]:
                out[-1] = (out[-1][0], max(out[-1][1], hi))
            else:
                out.append((lo, hi))
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return out


def good_intervals(speeds: tuple[int, ...]) -> list[tuple[F, F]]:
    cur = [(F(0), F(1))]
    for v in speeds:
        cur = intersect(cur, safe_intervals(v))
        if not cur:
            break
    return cur


def measure(intervals: list[tuple[F, F]]) -> F:
    return sum((r - l for l, r in intervals), F(0))


def frac(x: F) -> F:
    return x - (x.numerator // x.denominator)


def primitive_h(x: F) -> F:
    """Periodic primitive H of 1_[||x||<1/14]-1/7, with H(0)=0."""
    y = frac(x)
    if y <= THETA:
        return F(6, 7) * y
    if y <= 1 - THETA:
        return F(1, 14) - y / 7
    return -F(3, 49) + F(6, 7) * (y - (1 - THETA))


def endpoint_owners(core: tuple[int, ...], x: F) -> tuple[int, ...]:
    owners = []
    for v in core:
        y = frac(v * x)
        if y == THETA or y == 1 - THETA:
            owners.append(v)
    return tuple(owners)


def signed_terms(
    core: tuple[int, ...], intervals: list[tuple[F, F]], w: int
) -> list[EndpointTerm]:
    terms = []
    for l, r in intervals:
        for x, orientation in ((l, -1), (r, 1)):
            phase = frac(w * x)
            terms.append(
                EndpointTerm(
                    endpoint=x,
                    orientation=orientation,
                    owners=endpoint_owners(core, x),
                    phase=phase,
                    value=orientation * primitive_h(phase),
                )
            )
    return terms


def direct_final_measure(core: tuple[int, ...], w: int) -> F:
    return measure(good_intervals(tuple(sorted((*core, w)))))


def independent_danger_union_measure(speeds: tuple[int, ...]) -> F:
    """Second exact path: merge clipped danger arcs and subtract from one."""
    danger: list[tuple[F, F]] = []
    for v in speeds:
        for k in range(v):
            lo = F(k, v) - THETA / v
            hi = F(k, v) + THETA / v
            for shift in (-1, 0, 1):
                a = max(F(0), lo + shift)
                b = min(F(1), hi + shift)
                if a < b:
                    danger.append((a, b))
    danger.sort()
    merged: list[list[F]] = []
    for a, b in danger:
        if not merged or a > merged[-1][1]:
            merged.append([a, b])
        elif b > merged[-1][1]:
            merged[-1][1] = b
    return F(1) - sum((b - a for a, b in merged), F(0))


def analyze(core: tuple[int, ...], w: int) -> dict[str, object]:
    intervals = good_intervals(core)
    mu = measure(intervals)
    terms = signed_terms(core, intervals, w)
    endpoint_sum = sum((t.value for t in terms), F(0))
    epsilon = endpoint_sum / w
    final_from_identity = F(6, 7) * mu - epsilon
    final_direct = direct_final_measure(core, w)
    if final_from_identity != final_direct:
        raise AssertionError((core, w, final_from_identity, final_direct))
    n_raw = len(intervals)
    n_positive = sum(l < r for l, r in intervals)
    isolated = n_raw - n_positive
    # A zero-length component contributes one left and one right endpoint at
    # the same phase, hence exactly zero to the cocycle.  Only positive-length
    # components consume periodic-primitive variation.
    sharp_budget = F(6 * n_positive, 49 * w)
    owner_mult = sum(len(t.owners) for t in terms)
    return {
        "core": core,
        "w": w,
        "mu": mu,
        "n_raw": n_raw,
        "n_positive": n_positive,
        "isolated": isolated,
        "epsilon": epsilon,
        "endpoint_sum": endpoint_sum,
        "final": final_direct,
        "budget": sharp_budget,
        "budget_ratio": epsilon / sharp_budget if sharp_budget else F(0),
        "owner_mult": owner_mult,
        "terms": terms,
    }


def theorem_companion() -> None:
    print("THM-2162 SIGNED ENDPOINT / ONE-SWAP STABILITY -- EXACT COMPANION")
    print(f"extremizer={EXTREMIZER}")
    print(f"sharp_target={TARGET}")
    print(f"proper_neighbor_floor={LOCAL_NEXT}=target+{LOCAL_NEXT-TARGET}")
    print()
    print(
        "drop  core_mu        N+ iso  analytic_cutoff_for_neighbor_floor  "
        "finite_W  finite_min"
    )

    checked = 0
    equal_rows: list[tuple[int, int, F]] = []
    for drop in EXTREMIZER:
        core = tuple(v for v in EXTREMIZER if v != drop)
        intervals = good_intervals(core)
        mu = measure(intervals)
        n_positive = sum(l < r for l, r in intervals)
        isolated = len(intervals) - n_positive
        slack = F(6, 7) * mu - LOCAL_NEXT
        if slack <= 0:
            raise AssertionError(("nonpositive analytic slack", drop, slack))
        cutoff = F(6 * n_positive, 49) / slack
        finite_top = cutoff.numerator // cutoff.denominator

        finite_min: tuple[F, int] | None = None
        for w in range(14, finite_top + 1):
            row = analyze(core, w)
            checked += 1
            if row["final"] < LOCAL_NEXT:
                raise AssertionError(("neighbor below floor", drop, w, row["final"]))
            if row["final"] == LOCAL_NEXT:
                equal_rows.append((drop, w, row["final"]))
            if finite_min is None or row["final"] < finite_min[0]:
                finite_min = (row["final"], w)

        finite_text = (
            "none"
            if finite_min is None
            else f"{finite_min[0]}@W={finite_min[1]}"
        )
        finite_range = "empty" if finite_top < 14 else f"14..{finite_top}"
        print(
            f"{drop:>4}  {str(mu):<14} {n_positive:>2} {isolated:>3}  "
            f"{str(cutoff):<35} {finite_range:<9} {finite_text}"
        )

    expected_equal = [(10, 20, LOCAL_NEXT)]
    if equal_rows != expected_equal:
        raise AssertionError(("wrong equality classifier", equal_rows))

    # The only proper small replacement is the formerly missing speed 6.
    # These are the other AP-drop rows.  Their exact minimum is the drop-12
    # value, strictly above the newly found far-neighbor floor.
    small_rows = []
    for drop in EXTREMIZER:
        core = tuple(v for v in EXTREMIZER if v != drop)
        small_rows.append((direct_final_measure(core, 6), drop))
    small_rows.sort()
    if small_rows[0] != (F(426, 35035), 12):
        raise AssertionError(("unexpected small-neighbor minimum", small_rows[0]))
    if not LOCAL_NEXT < small_rows[0][0]:
        raise AssertionError("far-neighbor floor should beat every small neighbor")

    witness = tuple(sorted((*tuple(v for v in EXTREMIZER if v != 10), 20)))
    primary = measure(good_intervals(witness))
    independent = independent_danger_union_measure(witness)
    if primary != LOCAL_NEXT or independent != LOCAL_NEXT:
        raise AssertionError(("independent replay mismatch", primary, independent))

    print()
    print(f"finite_far_rows_checked={checked}")
    print(
        "far_equality_classifier="
        + ",".join(
            f"drop={drop}:W={w}:measure={value}"
            for drop, w, value in equal_rows
        )
    )
    print(f"small_neighbor_min={small_rows[0][0]} at drop={small_rows[0][1]}")
    print(f"independent_replay_row={witness}")
    print(f"independent_replay_measure={independent}")
    print("ALL EXACT CHECKS PASS")


def bounded_scout(bank: int, w_max: int, top: int) -> None:
    smallest: list[dict[str, object]] = []
    adverse: list[dict[str, object]] = []
    checked = 0
    for core in itertools.combinations(range(1, bank + 1), 11):
        if reduce(gcd, core) != 1:
            continue
        for w in range(14, w_max + 1):
            if w in core:
                continue
            row = analyze(core, w)
            checked += 1
            smallest.append(row)
            adverse.append(row)

    smallest.sort(key=lambda r: (r["final"], r["core"], r["w"]))
    adverse.sort(
        key=lambda r: (r["budget_ratio"], r["epsilon"], -r["final"]),
        reverse=True,
    )
    smallest = smallest[:top]
    adverse = adverse[:top]

    print("OPEN-Q-108 SIGNED-ENDPOINT / RAMANUJAN SCOUT -- EXACT")
    print(
        f"universe: primitive 11-subsets of [1,{bank}], "
        f"14<=W<={w_max}, W not in core; rows={checked}"
    )
    print(f"target={TARGET}")
    print()
    print("SMALLEST FINAL SAFE MEASURES")
    for row in smallest:
        print(
            f"final={row['final']} mu={row['mu']} eps={row['epsilon']} "
            f"N+={row['n_positive']} iso={row['isolated']} "
            f"W={row['w']} core={row['core']}"
        )
    print()
    print("LARGEST FRACTION OF THE SHARP PER-INTERVAL BUDGET")
    for row in adverse:
        print(
            f"ratio={row['budget_ratio']} eps={row['epsilon']} "
            f"sum={row['endpoint_sum']} N+={row['n_positive']} "
            f"iso={row['isolated']} W={row['w']} "
            f"owners={row['owner_mult']} core={row['core']}"
        )
    print()
    print(
        "identity verified exactly on every row; "
        f"below_target={sum(r['final'] < TARGET for r in smallest)} "
        "(count only within displayed smallest rows)"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("theorem", "scout"), default="theorem")
    parser.add_argument("--bank", type=int, default=14)
    parser.add_argument("--w-max", type=int, default=28)
    parser.add_argument("--top", type=int, default=12)
    args = parser.parse_args()
    if args.mode == "theorem":
        theorem_companion()
    else:
        bounded_scout(args.bank, args.w_max, args.top)


if __name__ == "__main__":
    main()
