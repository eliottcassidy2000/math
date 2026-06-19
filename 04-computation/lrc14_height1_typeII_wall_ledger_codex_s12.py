#!/usr/bin/env python3
"""
LRC(14) height-1 type-II support-six wall ledger.

This is the next finite slice after HYP-2612/HYP-2613.  The failed naive
Minkowski hope was "large span forces large coefficients"; HYP-2612 found
height-1 counterexamples such as

    21 = 1 + 2 + 3 + 4 + 5 + 6.

Here we isolate that obstruction instead of treating it as tail noise:

  * one large offset M above the bounded-spread finite-check ceiling B(k);
  * the remaining k-2 offsets lie in the bounded core {1,...,B(k)};
  * there is a height-1 type-II relation touching M, i.e.
        +/- M + sum_{i in S} +/- e_i = 0
    with at least five core terms, so THM-538's support-six floor is active.

For k=8,9,10, enumerate every primitive row of this form and check the exact
seven-sector cover measure against the canonical cap.  This does not prove
LRC(14): it only says the height-1 one-large resonance walls are harmless and
may be moved into the finite ledger before the true signed theta/Minkowski tail.

Tournament Analysis uses offsets as vertices with a wall-incidence observable;
this is a finite-obligation quotient, not a witness-time quotient.
"""
from __future__ import annotations

import itertools
from collections import defaultdict
from fractions import Fraction as F
from functools import reduce
from math import gcd
from time import perf_counter

try:
    import sys

    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)
except Exception:
    pass


CAP = {8: F(2243, 5880), 9: F(1979, 4004), 10: F(55, 91)}
BOUND = {8: 16, 9: 15, 10: 13}


def banner(title: str) -> None:
    print("\n" + "=" * 88)
    print(title)
    print("=" * 88)


def lcm(a: int, b: int) -> int:
    return a * b // gcd(a, b)


def gcd_all(xs: tuple[int, ...]) -> int:
    return reduce(gcd, [abs(x) for x in xs if x], 0)


def meas_s7(E: tuple[int, ...]) -> F:
    """Exact measure of hitting all seven fixed sectors."""
    nonzero = tuple(e for e in sorted(set(E)) if e)
    D = 7 * reduce(lcm, nonzero, 1)
    cuts = {0, D}
    for e in nonzero:
        step = D // (7 * e)
        for x in range(0, D + 1, step):
            cuts.add(x)
    total = F(0)
    ordered = sorted(cuts)
    for a, b in zip(ordered, ordered[1:]):
        if b <= a:
            continue
        num = a + b
        den = 2 * D
        sectors = {0}
        for e in nonzero:
            sectors.add((7 * e * num) // den % 7)
        if len(sectors) == 7:
            total += F(b - a, D)
    return total


def height1_typeII_rows(k: int) -> dict[tuple[tuple[int, ...], int], int]:
    """
    Enumerate primitive one-large rows with a height-1 type-II relation.

    Return {(core, M): support_size_of_shortest_height1_relation}.  The core
    has k-2 elements from [1,B(k)] and M>B(k).
    """
    B = BOUND[k]
    d = k - 2
    values = tuple(range(1, B + 1))
    rows: dict[tuple[tuple[int, ...], int], int] = {}
    signs = (-1, 1)

    for support_size_core in range(5, d + 1):
        for support in itertools.combinations(values, support_size_core):
            rest = tuple(x for x in values if x not in support)
            extras = list(itertools.combinations(rest, d - support_size_core)) or [()]
            for coeffs in itertools.product(signs, repeat=support_size_core):
                signed_sum = sum(c * v for c, v in zip(coeffs, support))
                for coeff_M in signs:
                    M = -signed_sum // coeff_M
                    if M <= B or -signed_sum != coeff_M * M:
                        continue
                    for extra in extras:
                        core = tuple(sorted(support + extra))
                        if gcd_all(core + (M,)) != 1:
                            continue
                        key = (core, M)
                        support_size = support_size_core + 1
                        old = rows.get(key)
                        if old is None or support_size < old:
                            rows[key] = support_size
    return rows


def row_report(k: int) -> dict:
    B = BOUND[k]
    cap = CAP[k]
    started = perf_counter()
    rows = height1_typeII_rows(k)
    generated = perf_counter() - started

    by_M: dict[int, int] = defaultdict(int)
    by_support: dict[int, int] = defaultdict(int)
    participation: dict[int, int] = defaultdict(int)
    worst = (F(-1), (), 0, 0)
    top: list[tuple[F, tuple[int, ...], int, int]] = []
    over_cap = 0

    measured_started = perf_counter()
    for (core, M), support_size in rows.items():
        by_M[M] += 1
        by_support[support_size] += 1
        for v in core + (M,):
            participation[v] += 1
        value = meas_s7((0,) + core + (M,))
        if value > cap:
            over_cap += 1
        item = (value, core, M, support_size)
        if value > worst[0]:
            worst = item
        top.append(item)

    measured = perf_counter() - measured_started
    top.sort(reverse=True)

    return {
        "k": k,
        "B": B,
        "cap": cap,
        "rows": rows,
        "row_count": len(rows),
        "generated_seconds": generated,
        "measured_seconds": measured,
        "distinct_M": len(by_M),
        "M_min": min(by_M) if by_M else None,
        "M_max": max(by_M) if by_M else None,
        "by_M": dict(by_M),
        "by_support": dict(sorted(by_support.items())),
        "participation": dict(participation),
        "over_cap": over_cap,
        "worst": worst,
        "top": top[:8],
    }


def tournament_fingerprint(report: dict) -> dict:
    """
    Tournament Analysis: vertices are offsets appearing in the wall ledger.

    Observable p(v) = number of rows in which offset v participates.  Orient
    v -> w when p(v) >= p(w), with numeric order breaking ties.  This quotient
    preserves which offsets carry height-1 wall obligations and destroys
    phase-location data.
    """
    scores = report["participation"]
    vertices = sorted(scores)
    outdeg = {v: 0 for v in vertices}
    flips = 0
    cycles = 0
    for a, b in itertools.combinations(vertices, 2):
        if scores[a] > scores[b] or (scores[a] == scores[b] and a < b):
            outdeg[a] += 1
        else:
            outdeg[b] += 1
            if a < b:
                flips += 1
    for a, b, c in itertools.combinations(vertices, 3):
        ab = scores[a] > scores[b] or (scores[a] == scores[b] and a < b)
        bc = scores[b] > scores[c] or (scores[b] == scores[c] and b < c)
        ac = scores[a] > scores[c] or (scores[a] == scores[c] and a < c)
        if (ab and bc and not ac) or ((not ab) and (not bc) and ac):
            cycles += 1
    hist: dict[int, int] = defaultdict(int)
    for deg in outdeg.values():
        hist[deg] += 1
    ordered = sorted(vertices, key=lambda v: (-scores[v], v))
    return {
        "vertices": len(vertices),
        "top_path": ordered[:12],
        "score_hist": dict(sorted(hist.items())),
        "edge_flips": flips,
        "directed_3_cycles": cycles,
    }


def main() -> None:
    banner("LRC(14) height-1 type-II wall ledger")
    print("Finite slice: one large offset, height-1 support-six relation, exact S7 measure.")
    print(f"Caps: {CAP}")
    print(f"Bounded core ceilings: {BOUND}")

    all_reports = [row_report(k) for k in (8, 9, 10)]

    for report in all_reports:
        k = report["k"]
        cap = report["cap"]
        worst_value, worst_core, worst_M, worst_support = report["worst"]
        print()
        print(f"k={k}, B={report['B']}: {report['row_count']} primitive rows")
        print(
            f"  distinct large offsets M={report['distinct_M']} "
            f"range=[{report['M_min']},{report['M_max']}], support histogram={report['by_support']}"
        )
        print(
            f"  generation {report['generated_seconds']:.3f}s, exact-measure scan "
            f"{report['measured_seconds']:.3f}s"
        )
        print(
            f"  cap={cap}={float(cap):.6f}; over_cap={report['over_cap']}; "
            f"worst={worst_value}={float(worst_value):.6f}; "
            f"margin={cap - worst_value}={float(cap - worst_value):.6f}"
        )
        print(f"  worst row: E=(0,{','.join(map(str, worst_core))},{worst_M}), support={worst_support}")
        print("  top rows:")
        for value, core, M, support_size in report["top"][:5]:
            print(
                f"    meas={value}={float(value):.6f}, "
                f"E=(0,{','.join(map(str, core))},{M}), support={support_size}"
            )
        tour = tournament_fingerprint(report)
        print("  Tournament Analysis (vertices=offsets, gauge=wall incidence):")
        print(f"    Hamiltonian tie path / incidence order: {tour['top_path']}")
        print(f"    score_hist={tour['score_hist']}, edge_flips={tour['edge_flips']}, "
              f"directed_3_cycles={tour['directed_3_cycles']}")

    banner("Reading")
    print(
        "1. The exact height-1 one-large support-six walls are finite and harmless: "
        "0 rows exceed cap for k=8,9,10."
    )
    print(
        "2. This directly discharges the HYP-2612 examples 21=sum(1..6) and "
        "the k=10 signed height-1 wall as proof-threatening cases; they are "
        "resonant but safely below cap."
    )
    print(
        "3. The result does not close LRC(14).  It covers one-large, height-1 "
        "type-II walls only.  Height >=2 walls, multi-large walls, and signed "
        "hyperplane tails remain the open HYP-2608/HYP-2613/HYP-2614 work."
    )
    print(
        "4. The proof split is now sharper: bounded finite check + finite "
        "height-1 wall ledger + relative signed support-six theta tail + "
        "signed-mass sequence spine."
    )
    print(
        "5. Assumption challenged: tournament vertices need not be runners. "
        "For this slice, wall obligations/offset incidences preserve the LRC "
        "predicate relevant to the support-six tail and discard witness-time geometry."
    )


if __name__ == "__main__":
    main()
