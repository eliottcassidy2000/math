#!/usr/bin/env python3
"""HYP-3114 scout: irrational/transcendental approximation sidecars for LRC14.

The scout computes exact LRC14 witness intervals for several named rows, then
places algebraic, transcendental, and Liouville-like targets inside the widest
positive component.  It measures continued-fraction denominators needed to
enter that component and compares them with the exact interval grid bound.

This is route selection, not a proof of LRC14.
"""

from __future__ import annotations

import itertools
import math
from collections import Counter
from dataclasses import dataclass
from decimal import Decimal, getcontext
from fractions import Fraction as F
from math import lcm

getcontext().prec = 120

N = 14


def fstr(x: F) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def dfrac(x: F) -> Decimal:
    return Decimal(x.numerator) / Decimal(x.denominator)


def ceil_frac(x: F) -> int:
    return -(-x.numerator // x.denominator)


def merge_intervals(intervals: list[tuple[F, F]]) -> list[tuple[F, F]]:
    intervals = sorted((lo, hi) for lo, hi in intervals if hi > lo)
    out: list[tuple[F, F]] = []
    for lo, hi in intervals:
        if out and lo <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], hi))
        else:
            out.append((lo, hi))
    return out


def forbidden_intervals(speed: int) -> list[tuple[F, F]]:
    radius = F(1, N * speed)
    intervals: list[tuple[F, F]] = []
    for m in range(speed):
        center = F(m, speed)
        lo = center - radius
        hi = center + radius
        if lo < 0:
            intervals.append((lo + 1, F(1)))
            intervals.append((F(0), hi))
        elif hi > 1:
            intervals.append((lo, F(1)))
            intervals.append((F(0), hi - 1))
        else:
            intervals.append((lo, hi))
    return intervals


def witness_intervals(row: tuple[int, ...]) -> list[tuple[F, F]]:
    bad: list[tuple[F, F]] = []
    for speed in row:
        bad.extend(forbidden_intervals(speed))
    merged = merge_intervals(bad)
    safe: list[tuple[F, F]] = []
    cur = F(0)
    for lo, hi in merged:
        if lo > cur:
            safe.append((cur, lo))
        cur = max(cur, hi)
    if cur < 1:
        safe.append((cur, F(1)))
    return safe


def circular_distance_to_integer(x: Decimal) -> Decimal:
    frac = x - Decimal(math.floor(x))
    return min(frac, Decimal(1) - frac)


def witness_margin_decimal(row: tuple[int, ...], t: Decimal) -> Decimal:
    threshold = Decimal(1) / Decimal(N)
    margins = [circular_distance_to_integer(Decimal(s) * t) - threshold for s in row]
    return min(margins)


def decimal_floor(x: Decimal) -> int:
    return int(x.to_integral_value(rounding="ROUND_FLOOR"))


def continued_fraction_convergents(x: Decimal, max_terms: int = 80):
    a_terms: list[int] = []
    y = +x
    p_nm2, p_nm1 = 0, 1
    q_nm2, q_nm1 = 1, 0
    for _ in range(max_terms):
        a = decimal_floor(y)
        a_terms.append(a)
        p = a * p_nm1 + p_nm2
        q = a * q_nm1 + q_nm2
        yield a_terms[:], p, q, Decimal(p) / Decimal(q)
        frac = y - Decimal(a)
        if abs(frac) < Decimal("1e-110"):
            break
        y = Decimal(1) / frac
        p_nm2, p_nm1 = p_nm1, p
        q_nm2, q_nm1 = q_nm1, q


def beta_constants() -> list[tuple[str, str, Decimal]]:
    sqrt2 = Decimal(2).sqrt()
    sqrt5 = Decimal(5).sqrt()
    e_minus_2 = Decimal(
        "0.718281828459045235360287471352662497757247093699959574966967627724076630353"
    )
    pi_minus_3 = Decimal(
        "0.141592653589793238462643383279502884197169399375105820974944592307816406286"
    )
    liouville = sum(Decimal(10) ** Decimal(-math.factorial(k)) for k in range(1, 7))
    return [
        ("phi_minus_1", "algebraic_badly_approximable", (sqrt5 - Decimal(1)) / Decimal(2)),
        ("sqrt2_minus_1", "quadratic_algebraic", sqrt2 - Decimal(1)),
        ("e_minus_2", "transcendental_explicit_cf", e_minus_2),
        ("pi_minus_3", "transcendental_unknown_measure", pi_minus_3),
        ("liouville_10", "transcendental_liouville_spike", liouville),
    ]


@dataclass(frozen=True)
class Row:
    name: str
    speeds: tuple[int, ...]
    note: str


ROWS = [
    Row("AP13_tight", tuple(range(1, 14)), "tight AP row; no positive interval"),
    Row("AP12_tail84", tuple(list(range(1, 12)) + [13, 84]), "small positive tail row"),
    Row(
        "loaded_B6",
        tuple(list(range(1, 12)) + [13, 84 * lcm(*range(1, 7))]),
        "THM-575 divisor-loaded raw-time obstruction",
    ),
    Row("single_tail168", tuple(list(range(1, 12)) + [13, 168]), "wider comparison tail"),
]


def first_cf_hit(target: Decimal, lo: Decimal, hi: Decimal) -> dict[str, object]:
    best = None
    for terms, p, q, value in continued_fraction_convergents(target):
        err = abs(target - value)
        inside = lo < value < hi
        if inside:
            exponent = float(-err.ln() / Decimal(q).ln()) if err > 0 and q > 1 else float("inf")
            best = {
                "p": p,
                "q": q,
                "err": err,
                "terms": terms,
                "max_partial": max(terms[1:] or [terms[0]]),
                "local_exponent": exponent,
            }
            break
    return best or {
        "p": None,
        "q": None,
        "err": None,
        "terms": [],
        "max_partial": None,
        "local_exponent": None,
    }


def interval_report(row: Row) -> dict[str, object]:
    intervals = witness_intervals(row.speeds)
    lengths = [hi - lo for lo, hi in intervals]
    if not intervals:
        return {
            "row": row,
            "intervals": intervals,
            "measure": F(0),
            "widest": None,
            "widest_len": F(0),
            "grid_bound": None,
        }
    widest_i = max(range(len(intervals)), key=lambda i: lengths[i])
    widest = intervals[widest_i]
    widest_len = lengths[widest_i]
    grid_bound = (widest_len.denominator // widest_len.numerator) + 1
    return {
        "row": row,
        "intervals": intervals,
        "measure": sum(lengths, F(0)),
        "widest": widest,
        "widest_len": widest_len,
        "grid_bound": grid_bound,
    }


AXES = {
    "preserves_lrc_predicate": {
        "finite_interval_margin": 5,
        "observer_gluing_packet": 5,
        "continued_fraction_packet": 4,
        "algebraic_roth_height_fence": 4,
        "transcendental_measure_sidecar": 3,
        "liouville_spike_schedule": 3,
        "raw_named_constant": 0,
    },
    "controls_destroyed_coordinate": {
        "observer_gluing_packet": 5,
        "finite_interval_margin": 5,
        "continued_fraction_packet": 4,
        "liouville_spike_schedule": 4,
        "algebraic_roth_height_fence": 4,
        "transcendental_measure_sidecar": 3,
        "raw_named_constant": 0,
    },
    "finite_checkable": {
        "finite_interval_margin": 5,
        "continued_fraction_packet": 5,
        "observer_gluing_packet": 4,
        "liouville_spike_schedule": 4,
        "algebraic_roth_height_fence": 3,
        "transcendental_measure_sidecar": 2,
        "raw_named_constant": 5,
    },
    "handles_exceptional_approximants": {
        "algebraic_roth_height_fence": 5,
        "liouville_spike_schedule": 5,
        "continued_fraction_packet": 4,
        "transcendental_measure_sidecar": 4,
        "finite_interval_margin": 3,
        "observer_gluing_packet": 3,
        "raw_named_constant": 0,
    },
    "feeds_denominator_grid": {
        "finite_interval_margin": 5,
        "continued_fraction_packet": 5,
        "observer_gluing_packet": 4,
        "algebraic_roth_height_fence": 3,
        "transcendental_measure_sidecar": 3,
        "liouville_spike_schedule": 2,
        "raw_named_constant": 1,
    },
}

VERTICES = [
    "finite_interval_margin",
    "observer_gluing_packet",
    "continued_fraction_packet",
    "algebraic_roth_height_fence",
    "transcendental_measure_sidecar",
    "liouville_spike_schedule",
    "raw_named_constant",
]


def tournament_report() -> dict[str, object]:
    tie_order = {v: i for i, v in enumerate(VERTICES)}
    edges = {v: set() for v in VERTICES}
    margins = []
    for a, b in itertools.combinations(VERTICES, 2):
        aw = bw = 0
        for scores in AXES.values():
            if scores[a] > scores[b]:
                aw += 1
            elif scores[b] > scores[a]:
                bw += 1
        if aw > bw or (aw == bw and tie_order[a] < tie_order[b]):
            edges[a].add(b)
            margins.append((aw - bw, a, b, aw, bw))
        else:
            edges[b].add(a)
            margins.append((bw - aw, b, a, bw, aw))
    score_hist = Counter(len(edges[v]) for v in VERTICES)
    cycles = 0
    for a, b, c in itertools.combinations(VERTICES, 3):
        if (
            (b in edges[a] and c in edges[b] and a in edges[c])
            or (c in edges[a] and b in edges[c] and a in edges[b])
        ):
            cycles += 1
    return {
        "scores": {v: len(edges[v]) for v in VERTICES},
        "score_hist": dict(sorted(score_hist.items())),
        "directed_3cycles": cycles,
        "priority_path": sorted(VERTICES, key=lambda v: -len(edges[v])),
        "edge_flip_risks": sorted(margins)[:10],
    }


def main() -> None:
    print("HYP-3114 irrational/transcendental approximation sidecar scout -- codex S265")
    print("Approximation enters LRC14 through positive witness intervals and margin retention.")
    print()
    print("=" * 96)
    print("MAP 1: EXACT LRC14 WITNESS INTERVALS")
    print("=" * 96)
    reports = [interval_report(row) for row in ROWS]
    for rep in reports:
        row: Row = rep["row"]  # type: ignore[assignment]
        print(f"{row.name:16s} max_speed={max(row.speeds):5d} components={len(rep['intervals']):4d} measure={fstr(rep['measure'])}")  # type: ignore[arg-type]
        if rep["widest"] is None:
            print("  no positive component")
        else:
            lo, hi = rep["widest"]  # type: ignore[misc]
            print(
                f"  widest=({fstr(lo)}, {fstr(hi)}) len={fstr(rep['widest_len'])} "
                f"all-denominator grid bound q>{rep['grid_bound'] - 1}"
            )
    print()
    print("=" * 96)
    print("MAP 2: CONTINUED-FRACTION ENTRY INTO POSITIVE COMPONENTS")
    print("=" * 96)
    for rep in reports:
        row: Row = rep["row"]  # type: ignore[assignment]
        if rep["widest"] is None:
            continue
        lo_f, hi_f = rep["widest"]  # type: ignore[misc]
        lo = dfrac(lo_f)
        hi = dfrac(hi_f)
        length = hi - lo
        print(f"\n{row.name}: widest_len={fstr(rep['widest_len'])}, grid_bound={rep['grid_bound']}")
        print("constant                 type                              q_hit      max_a  exponent  robust_radius")
        for name, typ, beta in beta_constants():
            target = lo + beta * length
            margin = witness_margin_decimal(row.speeds, target)
            robust_radius = margin / Decimal(max(row.speeds))
            hit = first_cf_hit(target, lo, hi)
            q = hit["q"]
            q_text = str(q) if q is not None else "none"
            max_a = hit["max_partial"]
            exp = hit["local_exponent"]
            exp_text = "inf" if exp == float("inf") else ("none" if exp is None else f"{exp:.3f}")
            print(
                f"{name:24s} {typ:31s} {q_text:>9s} {str(max_a):>10s} "
                f"{exp_text:>8s}  {robust_radius:.3E}"
            )
    print()
    print("=" * 96)
    print("MAP 3: TOURNAMENT ANALYSIS ON APPROXIMATION CARRIERS")
    print("=" * 96)
    tr = tournament_report()
    print("vertices are proof carriers, not constants.")
    print("axes=" + ", ".join(AXES))
    print(f"score_hist={tr['score_hist']}")
    print(f"scores={tr['scores']}")
    print(f"directed_3cycles={tr['directed_3cycles']}")
    print("priority_path=" + " -> ".join(tr["priority_path"]))  # type: ignore[arg-type]
    print("low-margin edge-flip risks:")
    for margin, winner, loser, aw, bw in tr["edge_flip_risks"]:  # type: ignore[assignment]
        print(f"  margin={margin}: {winner} -> {loser} ({aw}-{bw} axes)")
    print()
    print("=" * 96)
    print("SYNTHESIS")
    print("=" * 96)
    print("strengthened:")
    print("  AP13 has no positive interval; approximation cannot repair a boundary-only row.")
    print("  The loaded_B6 row has widest component 1/5880, matching the raw-time denominator wall.")
    print("  Algebraic/transcendental/Liouville labels matter only through retained CF/measure sidecars.")
    print("still open:")
    print("  Replace widest-component scouts by the normalized THM-565 slow/ruler-coordinate intervals.")
    print("  Attach interval-margin fields to HYP-3112 ear payload and HYP-3098 observer-gluing rows.")
    print("challenged assumption:")
    print("  The useful vertex is not an irrational constant.  It is the interval-plus-margin packet")
    print("  that lets approximation produce a legal finite-grid witness.")
    print("DONE.")


if __name__ == "__main__":
    main()
