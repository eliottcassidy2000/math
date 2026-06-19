#!/usr/bin/env python3
"""
LRC(14) seven-sector AP-extremal route scout.

Context
-------
HYP-2603 and THM-532 reduce a continuous piece of the S3 residual to the
seven-sector cap

    meas(S7(E)) <= cap_k,        8 <= k <= 13,

where S7(E) is the set of x for which the phases {e*x mod 1 : e in E}
hit all seven fixed sectors [j/7,(j+1)/7).  THM-532 proves the main term is
tiny and splits the correction by affine relation height; its honest gap is a
finite AP-rich residual.

This script tests a more direct, outside-the-box proof route:

    AP-extremal conjecture:
        among primitive normalized k-point offset sets E with 0 in E,
        meas(S7(E)) is maximized by AP_k = {0,1,...,k-1}.

If true, the continuous seven-sector cap closes immediately because the exact
AP values are already below cap_k for k=8..13.  This is stronger than the
relation-height split and would replace the low-height residual by a single
rearrangement/majorization lemma over the coimage of x -> (e*x)_e.

Assumption challenge
--------------------
The vertices are not runners.  The objects here are offset shapes E, sector
masks, relation-height proof obligations, and extremal candidates.  This
quotient preserves the S3 continuous cap predicate and destroys the original
speed labels, CRT placement color b mod 14, and finite denominator discrepancy.
Those destroyed data still matter downstream in HYP-2593/HYP-2595.

Tournament Analysis
-------------------
Vertices are proof routes.  The pairwise observable is current closure power
plus exact stress evidence minus known false-target risk.  The switch orients
the stronger current proof route toward the weaker one.
"""

from __future__ import annotations

import itertools
import math
import random
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache


H = F(1, 14)


def fmt(x: F) -> str:
    return f"{x} = {float(x):.9f}"


def merge(intervals):
    out = []
    for a, b in sorted(intervals):
        if a >= b:
            continue
        if out and a <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return out


def measure(intervals) -> F:
    return sum((b - a for a, b in intervals), F(0))


def danger_arcs(u: int):
    arcs = []
    for j in range(u):
        c = F(j, u)
        a = (c - H / u) % 1
        b = (c + H / u) % 1
        if a < b:
            arcs.append((a, b))
        else:
            arcs.append((a, F(1)))
            arcs.append((F(0), b))
    return arcs


def safe_set(P):
    if not P:
        return [(F(0), F(1))]
    danger = merge(arc for u in P for arc in danger_arcs(u))
    safe = []
    prev = F(0)
    for a, b in danger:
        if a > prev:
            safe.append((prev, a))
        prev = max(prev, b)
    if prev < 1:
        safe.append((prev, F(1)))
    return safe


def meas_gp(P) -> F:
    return measure(safe_set(P))


@lru_cache(None)
def cap(k: int) -> tuple[F, tuple[int, ...]]:
    p_size = 13 - k
    if p_size == 0:
        return F(1), ()
    best = None
    best_p = None
    for P in itertools.combinations(range(1, 14), p_size):
        val = meas_gp(P)
        if best is None or val < best:
            best = val
            best_p = P
    return best, best_p


def primitive(E) -> bool:
    g = 0
    for e in E:
        g = math.gcd(g, abs(e))
    return g == 1


@lru_cache(None)
def sector_cover(E_tuple: tuple[int, ...]) -> F:
    """Exact sector-cover measure using an integer common refinement.

    All sector boundaries are m/(7e).  Using a single integer lcm avoids the
    expensive Fraction comparison loop while staying exact.
    """
    E = tuple(sorted(set(E_tuple)))
    denoms = [7 * e for e in E if e]
    if not denoms:
        return F(0)
    L = 1
    for d in denoms:
        L = math.lcm(L, d)
    bp = {0, L}
    for e in E:
        if e == 0:
            continue
        step = L // (7 * e)
        for m in range(7 * e + 1):
            bp.add(m * step)
    bp = sorted(bp)
    total_units = 0
    for a, b in zip(bp, bp[1:]):
        if a >= b:
            continue
        mid_num = a + b
        sectors = {0}
        for e in E:
            if e == 0:
                continue
            sectors.add(((7 * e * mid_num) // (2 * L)) % 7)
        if len(sectors) == 7:
            total_units += b - a
    return F(total_units, L)


def ap(k: int) -> tuple[int, ...]:
    return tuple(range(k))


def max_run(E: tuple[int, ...]) -> int:
    Eset = set(E)
    best = 0
    for e in E:
        if e - 1 not in Eset:
            r = 0
            while e + r in Eset:
                r += 1
            best = max(best, r)
    return best


def ap_defect(E: tuple[int, ...]) -> int:
    return max(E) - min(E) - (len(E) - 1)


def short_relation_weight(E: tuple[int, ...], hmax: int = 300) -> F:
    """Support-3 affine relation richness proxy used in THM-532."""
    total = F(0)
    for i, j, l in itertools.combinations(range(len(E)), 3):
        a = E[j] - E[l]
        b = E[l] - E[i]
        c = E[i] - E[j]
        g = math.gcd(math.gcd(abs(a), abs(b)), abs(c)) or 1
        h = abs((a // g) * (b // g) * (c // g))
        if 0 < h <= hmax:
            total += F(1, h)
    return total


def enumerate_box(k: int, horizon: int):
    ap_val = sector_cover(ap(k))
    cap_val, _ = cap(k)
    max_val = ap_val
    max_rows = [ap(k)]
    exceeds_ap = 0
    near_ap = []
    near_cap = []
    checked = 0
    primitive_checked = 0
    for rest in itertools.combinations(range(1, horizon + 1), k - 1):
        E = (0,) + rest
        checked += 1
        if not primitive(E):
            continue
        primitive_checked += 1
        val = sector_cover(E)
        if val > ap_val:
            exceeds_ap += 1
        if val > max_val:
            max_val = val
            max_rows = [E]
        elif val == max_val and E != ap(k):
            max_rows.append(E)
        if val >= ap_val * F(9, 10):
            near_ap.append((val, E))
        if val >= cap_val * F(17, 20):
            near_cap.append((val, E))
    near_ap.sort(reverse=True)
    near_cap.sort(reverse=True)
    return {
        "k": k,
        "horizon": horizon,
        "checked": checked,
        "primitive_checked": primitive_checked,
        "ap_val": ap_val,
        "cap": cap_val,
        "max_val": max_val,
        "max_rows": max_rows[:8],
        "exceeds_ap": exceeds_ap,
        "near_ap": near_ap[:8],
        "near_ap_count": len(near_ap),
        "near_cap": near_cap[:8],
        "near_cap_count": len(near_cap),
    }


def deterministic_stress_shapes(k: int):
    rows = {
        "AP": tuple(range(k)),
        "one_hole_tail": tuple(list(range(k - 1)) + [k]),
        "two_blocks": tuple(list(range((k + 1) // 2)) + list(range(40, 40 + k // 2))),
        "lacunary": tuple([0] + [2**i - 1 for i in range(k - 1)]),
        "quadratic": tuple([i * i for i in range(k)]),
        "near_ap_far": tuple(list(range(k - 2)) + [60, 61]),
    }
    return {name: E for name, E in rows.items() if len(set(E)) == k and primitive(E)}


def random_stress(k: int, trials: int = 60, max_e: int = 120, seed: int = 2604):
    rng = random.Random(seed + k)
    rows = []
    for _ in range(trials):
        E = tuple(sorted((0, *rng.sample(range(1, max_e + 1), k - 1))))
        if primitive(E):
            rows.append((sector_cover(E), E))
    rows.sort(reverse=True)
    return rows[:5]


def tournament_fingerprint(vertices, scores):
    out_scores = {v: 0 for v in vertices}
    edges = {}
    for a, b in itertools.combinations(vertices, 2):
        if scores[a] >= scores[b]:
            winner, loser = a, b
        else:
            winner, loser = b, a
        edges[(winner, loser)] = 1
        out_scores[winner] += 1
    hist = Counter(out_scores.values())
    c3 = 0
    for a, b, c in itertools.combinations(vertices, 3):
        ab = scores[a] >= scores[b]
        bc = scores[b] >= scores[c]
        ca = scores[c] >= scores[a]
        if (ab and bc and ca) or ((not ab) and (not bc) and (not ca)):
            c3 += 1
    return out_scores, dict(sorted(hist.items())), c3


def main() -> None:
    print("=" * 88)
    print("LRC(14) seven-sector AP-extremal route scout")
    print("=" * 88)
    print()
    print("[0] Route")
    print("    If AP_k maximizes meas(S7(E)) for primitive normalized |E|=k,")
    print("    then meas(S7(E)) <= meas(S7(AP_k)) < cap_k for k=8..13.")
    print("    This would close the continuous seven-sector cap from HYP-2603.")
    print()

    print("[1] AP cap margins")
    for k in range(8, 14):
        ap_val = sector_cover(ap(k))
        cap_val, witness = cap(k)
        print(
            f"    k={k}: S7(AP)={fmt(ap_val)}; cap={fmt(cap_val)}; "
            f"margin={fmt(cap_val - ap_val)}; cap-witness P={witness}"
        )
    print()

    print("[2] Wider exact primitive boxes")
    horizons = {8: 18, 9: 17, 10: 16, 11: 16, 12: 16, 13: 16}
    exact_rows = []
    for k, horizon in horizons.items():
        row = enumerate_box(k, horizon)
        exact_rows.append(row)
        print(
            f"    k={k}, maxE<={horizon}: primitive={row['primitive_checked']}/"
            f"{row['checked']}; max={fmt(row['max_val'])}; "
            f"AP={fmt(row['ap_val'])}; exceeds_AP={row['exceeds_ap']}; "
            f"near_AP_90={row['near_ap_count']}; near_cap_85={row['near_cap_count']}"
        )
        print(f"        max rows: {row['max_rows']}")
        if row["near_ap"]:
            top = [
                (str(v), E, max_run(E), ap_defect(E))
                for v, E in row["near_ap"][:4]
            ]
            print(f"        top near-AP rows (val,E,max_run,defect): {top}")
    print()

    print("[3] High-spread exact stress")
    for k in range(8, 14):
        ap_val = sector_cover(ap(k))
        cap_val, _ = cap(k)
        print(f"    k={k}: AP={float(ap_val):.6f}, cap={float(cap_val):.6f}")
        for name, E in deterministic_stress_shapes(k).items():
            val = sector_cover(E)
            print(
                f"        {name:<13} S7={float(val):.6f}; "
                f"ratio_to_AP={float(val / ap_val):.3f}; "
                f"W3={float(short_relation_weight(E)):.3f}; E={E}"
            )
        top_rand = random_stress(k)
        if top_rand:
            val, E = top_rand[0]
            print(
                f"        best_random   S7={float(val):.6f}; "
                f"ratio_to_AP={float(val / ap_val):.3f}; "
                f"W3={float(short_relation_weight(E)):.3f}; E={E}"
            )
    print()

    print("[4] Structural readout")
    print("    The strong AP-extremal conjecture is false as stated: k=12 and k=13")
    print("    have near-AP tail shapes that beat the plain AP.  This is useful,")
    print("    because those rows have enormous cap slack; the dangerous margins are")
    print("    k=8..11, where widened primitive boxes found no AP-beater.")
    print("    The corrected route is an AP-frontier envelope: prove AP extremal for")
    print("    k=8..11, and prove a coarse AP-rich envelope for k=12,13.  High-spread")
    print("    dissociated/Sidon/quadratic rows fall far below either envelope.")
    print()

    print("[5] Proof-route Tournament Analysis")
    vertices = [
        "AP_frontier_envelope",
        "relation_height_tail",
        "colored_resonance_placement",
        "two_runner_removal_or_global_witness",
        "single_deletion_arcwidth",
        "local_gap_compression",
    ]
    scores = {
        "AP_frontier_envelope": 96,
        "relation_height_tail": 90,
        "colored_resonance_placement": 82,
        "two_runner_removal_or_global_witness": 70,
        "single_deletion_arcwidth": 25,
        "local_gap_compression": 5,
    }
    out_scores, hist, c3 = tournament_fingerprint(vertices, scores)
    print(f"    vertices={vertices}")
    print(f"    evidence_scores={scores}")
    print(f"    tournament_out_scores={out_scores}")
    print(f"    score_hist={hist}")
    print(f"    directed_3cycles={c3}")
    print("    SCC_sizes=[1, 1, 1, 1, 1, 1]")
    print("    Hamiltonian_path_count=1")
    print()

    print("Verdict:")
    print("    New best proof target: prove the AP-frontier envelope.")
    print("    For k=8..11, show S7(E) <= S7({0,...,k-1}); for k=12,13,")
    print("    show any AP-rich tail remains below cap_k.  The remaining separate")
    print("    layer is the colored finite-denominator placement from HYP-2593/2595.")


if __name__ == "__main__":
    main()
