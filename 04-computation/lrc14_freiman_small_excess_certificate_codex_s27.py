#!/usr/bin/env python3
"""
HYP-2638 / T886: Freiman 3k-4 finite pocket for the LRC(14)
additive-energy route.

For a k-point integer set E, write

    exc(E) = |E+E| - (2k-1).

Freiman 3k-4 says that if exc(E) <= k-3, then E is contained in an arithmetic
progression of length k+exc(E).  Translation and dilation preserve the sector
functional used in HYP-2607/HYP-2635, so this whole small-excess pocket has a
finite primitive normal-form certificate.

This script enumerates those normal forms for k=8,9,10, computes exact L_y and
p0=meas(S7), and records maxima by exact excess.

Tournament Analysis note:
  vertices are proof obligations, not runners/arcs:
    Freiman_3k4_pocket > exact_excess_layer > AP_invariance >
    relation_fiber_cover > higher_rank_GAP_tail > raw_pair_sum_energy >
    raw_runner_vertices.
  The quotient preserves the small-excess LRC sector predicate only because
  every normal form is evaluated exactly; it destroys reciprocal-tail sign data.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction
from functools import reduce
from itertools import combinations
from math import gcd


CAPS = {
    8: Fraction(2243, 5880),
    9: Fraction(1979, 4004),
    10: Fraction(55, 91),
}


def gcd_all(values: tuple[int, ...]) -> int:
    return reduce(gcd, values, 0)


def primitive(E: tuple[int, ...]) -> bool:
    return gcd_all(E) == 1


def N_at(E: tuple[int, ...], x: Fraction) -> int:
    hit = set()
    for e in E:
        v = e * x
        v -= v.numerator // v.denominator
        hit.add((v.numerator * 7) // v.denominator)
    return sum(1 for j in range(1, 7) if j not in hit)


def dist_p(E: tuple[int, ...]) -> list[Fraction]:
    E = tuple(sorted(set(E)))
    bps = {Fraction(0), Fraction(1)}
    for e in E:
        if e == 0:
            continue
        for a in range(0, 7 * e + 1):
            bps.add(Fraction(a, 7 * e))
    ordered = sorted(b for b in bps if 0 <= b <= 1)
    p = [Fraction(0)] * 7
    for lo, hi in zip(ordered, ordered[1:]):
        if hi == lo:
            continue
        p[N_at(E, (lo + hi) / 2)] += hi - lo
    return p


def g_poly(k: int, t: int) -> Fraction:
    if k == 8:
        return Fraction((t - 1) * (t - 2) * (t - 4) * (t - 5), 40)
    if k in (9, 10):
        return Fraction(-(t - 2) * (t - 3) * (t - 6), 36)
    return Fraction((t - 3) * (t - 4), 12)


def L_y(E: tuple[int, ...], k: int) -> Fraction:
    p = dist_p(E)
    return sum(p[t] * g_poly(k, t) for t in range(7))


def sumset_excess(E: tuple[int, ...]) -> int:
    sums = {a + b for a in E for b in E}
    return len(sums) - (2 * len(E) - 1)


def normal_forms(k: int) -> list[tuple[int, ...]]:
    """Primitive affine normal forms in the largest 3k-4 AP hull."""
    max_coord = 2 * k - 4
    forms = []
    for tail in combinations(range(1, max_coord + 1), k - 1):
        E = (0,) + tail
        if primitive(E):
            forms.append(E)
    return forms


def fmt(fr: Fraction) -> str:
    return f"{fr} = {float(fr):.9f}"


def summarize_k(k: int) -> dict[str, object]:
    ap = tuple(range(k))
    ap_L = L_y(ap, k)
    ap_p = dist_p(ap)
    ap_p0 = ap_p[0]
    cap = CAPS[k]
    max_allowed_excess = k - 3

    by_excess: dict[int, dict[str, object]] = {}
    hull_failures: list[tuple[int, tuple[int, ...]]] = []
    count_by_excess = Counter()
    total_small_excess = 0
    total_forms = 0

    for E in normal_forms(k):
        total_forms += 1
        exc = sumset_excess(E)
        if exc > max_allowed_excess:
            continue
        total_small_excess += 1
        count_by_excess[exc] += 1
        hull_length = max(E) + 1
        if hull_length > k + exc:
            hull_failures.append((exc, E))

        p = dist_p(E)
        Ly = sum(p[t] * g_poly(k, t) for t in range(7))
        row = by_excess.setdefault(
            exc,
            {
                "max_L": (Fraction(-1), ()),
                "max_p0": (Fraction(-1), ()),
                "max_span": (0, ()),
                "min_L": (Fraction(10), ()),
            },
        )
        if Ly > row["max_L"][0]:
            row["max_L"] = (Ly, E)
        if p[0] > row["max_p0"][0]:
            row["max_p0"] = (p[0], E)
        if hull_length > row["max_span"][0]:
            row["max_span"] = (hull_length, E)
        if Ly < row["min_L"][0]:
            row["min_L"] = (Ly, E)

    positive_L = [
        (row["max_L"][0], exc, row["max_L"][1])
        for exc, row in by_excess.items()
        if exc > 0
    ]
    positive_p0 = [
        (row["max_p0"][0], exc, row["max_p0"][1])
        for exc, row in by_excess.items()
        if exc > 0
    ]
    best_pos_L = max(positive_L) if positive_L else (Fraction(0), None, ())
    best_pos_p0 = max(positive_p0) if positive_p0 else (Fraction(0), None, ())

    return {
        "k": k,
        "cap": cap,
        "ap_L": ap_L,
        "ap_p0": ap_p0,
        "total_forms": total_forms,
        "total_small_excess": total_small_excess,
        "count_by_excess": count_by_excess,
        "by_excess": by_excess,
        "hull_failures": hull_failures,
        "best_pos_L": best_pos_L,
        "best_pos_p0": best_pos_p0,
    }


def tournament_fingerprint() -> dict[str, object]:
    vertices = [
        "Freiman_3k4_pocket",
        "exact_excess_layer",
        "AP_invariance",
        "relation_fiber_cover",
        "higher_rank_GAP_tail",
        "raw_pair_sum_energy",
        "raw_runner_vertices",
    ]
    n = len(vertices)
    scores = {v: n - i - 1 for i, v in enumerate(vertices)}
    score_hist = Counter(scores.values())
    directed_3_cycles = 0
    # The declared order orients every pair forward, so SCCs are singletons and
    # the only Hamiltonian path is the order itself.
    return {
        "vertices": vertices,
        "score_hist": dict(sorted(score_hist.items())),
        "directed_3_cycles": directed_3_cycles,
        "scc_sizes": [1] * n,
        "hamiltonian_paths": 1,
    }


def main() -> None:
    print("=== HYP-2638 / T886: Freiman 3k-4 small-excess finite certificate ===")
    print("Method: primitive affine normal forms after translation/dilation.")
    print("Freiman pocket: exc(E)<=k-3, so expected AP hull length <= k+exc(E).\n")

    all_ok = True
    summaries = [summarize_k(k) for k in (8, 9, 10)]
    for summary in summaries:
        k = summary["k"]
        cap = summary["cap"]
        ap_L = summary["ap_L"]
        ap_p0 = summary["ap_p0"]
        best_pos_L, best_pos_L_exc, best_pos_L_E = summary["best_pos_L"]
        best_pos_p0, best_pos_p0_exc, best_pos_p0_E = summary["best_pos_p0"]
        if summary["hull_failures"] or best_pos_L >= ap_L or best_pos_L >= cap:
            all_ok = False

        print(f"--- k={k} ---")
        print(f"primitive normal forms in [0,{2*k-4}]: {summary['total_forms']}")
        print(
            f"small-excess forms exc<=k-3={k-3}: "
            f"{summary['total_small_excess']}"
        )
        print(f"AP L_y: {fmt(ap_L)}")
        print(f"AP p0=meas(S7): {fmt(ap_p0)}")
        print(f"cap_k: {fmt(cap)}")
        print(
            f"best positive-excess L_y: {fmt(best_pos_L)} "
            f"at exc={best_pos_L_exc}, E={best_pos_L_E}"
        )
        print(f"AP - best_positive_L_y: {fmt(ap_L - best_pos_L)}")
        print(f"cap - best_positive_L_y: {fmt(cap - best_pos_L)}")
        print(
            f"best positive-excess p0: {fmt(best_pos_p0)} "
            f"at exc={best_pos_p0_exc}, E={best_pos_p0_E}"
        )
        print(f"AP p0 - best_positive_p0: {fmt(ap_p0 - best_pos_p0)}")
        print("maxima by exact excess:")
        for exc in range(0, k - 2):
            row = summary["by_excess"].get(exc)
            if not row:
                continue
            max_L, max_L_E = row["max_L"]
            max_p0, max_p0_E = row["max_p0"]
            max_span, max_span_E = row["max_span"]
            print(
                f"  exc={exc:2d} count={summary['count_by_excess'][exc]:5d} "
                f"maxL={float(max_L):.9f} E={max_L_E} "
                f"maxp0={float(max_p0):.9f} p0E={max_p0_E} "
                f"max_hull={max_span} hullE={max_span_E}"
            )
        if summary["hull_failures"]:
            print("HULL FAILURES:")
            for exc, E in summary["hull_failures"][:10]:
                print(f"  exc={exc} hull={max(E)+1} bound={k+exc} E={E}")
        else:
            print("Freiman hull check: 0 failures in enumerated small-excess forms.")
        print()

    print("=== Tournament Analysis ===")
    fp = tournament_fingerprint()
    print("Hamiltonian path:")
    print("  " + " > ".join(fp["vertices"]))
    print(f"score_hist = {fp['score_hist']}")
    print(f"directed_3_cycles = {fp['directed_3_cycles']}")
    print(f"SCC_sizes = {fp['scc_sizes']}")
    print(f"hamiltonian_paths = {fp['hamiltonian_paths']}")
    print()

    if all_ok:
        print("CERTIFICATE STATUS: PASS for k=8,9,10 finite Freiman 3k-4 pocket.")
        print("Interpretation: excess 0 is the AP row; every positive-excess")
        print("small-doubling normal form has a strict L_y margin below AP and cap.")
    else:
        print("CERTIFICATE STATUS: ATTENTION REQUIRED.")
        print("A hull failure or positive-excess AP/cap challenge appeared.")


if __name__ == "__main__":
    main()
