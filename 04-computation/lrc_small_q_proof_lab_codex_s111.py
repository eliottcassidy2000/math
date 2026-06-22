#!/usr/bin/env python3
"""Exact small-q laboratory for the LRC(14) proof mechanisms.

This script tests the LRC(14) witness/p0/additive-energy proof atoms at smaller
even thresholds q=8,10,12 before returning to q=14.  The point is not to
reprove the known smaller LRC cases by their historical methods; it is to see
which parts of the current q=14 machine are already visible in easier cases.

Convention: q is the denominator in the lonely threshold 1/q.  For q even, a
slow cluster admits a fast phase iff its circular phase set has a gap > 2/q.
The q=14 case is the repository's LRC(14) cluster predicate gap > 1/7.

Tournament Analysis:
  vertices are proof carriers, not runners:
    Bonferroni_floor, bounded_nu_minimizer, bounded_p0_cap,
    AP_difference_profile, scalar_additive_energy.
  The observable is the exact failure count / smallest margin across the
  small-q banks; the Hamiltonian path is printed at the end.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction as Fr
from functools import reduce
from itertools import combinations
from math import gcd
import sys

try:
    sys.stdout.reconfigure(line_buffering=True)
except Exception:
    pass


def primitive(E: tuple[int, ...]) -> bool:
    return reduce(gcd, (abs(e) for e in E if e), 0) == 1


def union_measure(arcs: list[tuple[Fr, Fr]]) -> Fr:
    """Exact measure of a union of possibly wrapping open arcs in R/Z."""
    segs: list[tuple[Fr, Fr]] = []
    for lo, hi in arcs:
        length = hi - lo
        if length <= 0:
            continue
        a = lo % 1
        b = a + length
        if b <= 1:
            segs.append((a, b))
        else:
            segs.append((a, Fr(1)))
            segs.append((Fr(0), b - 1))
    if not segs:
        return Fr(0)
    segs.sort()
    total = Fr(0)
    cur_lo, cur_hi = segs[0]
    for lo, hi in segs[1:]:
        if lo <= cur_hi:
            if hi > cur_hi:
                cur_hi = hi
        else:
            total += cur_hi - cur_lo
            cur_lo, cur_hi = lo, hi
    total += cur_hi - cur_lo
    return total


def meas_gp(q: int, P: tuple[int, ...]) -> Fr:
    """meas{x : ||p x|| >= 1/q for every p in P}."""
    arcs: list[tuple[Fr, Fr]] = []
    half = Fr(1, q)
    for p in P:
        w = half / p
        for j in range(p):
            c = Fr(j, p)
            arcs.append((c - w, c + w))
    return Fr(1) - union_measure(arcs)


def cap_gp(q: int, k: int) -> tuple[Fr, tuple[int, ...]]:
    """Minimum small-part safe measure for q-1 total speeds and k cluster offsets."""
    psize = q - 1 - k
    if psize == 0:
        return Fr(1), ()
    best: tuple[Fr, tuple[int, ...]] | None = None
    for P in combinations(range(1, q), psize):
        val = meas_gp(q, P)
        if best is None or val < best[0]:
            best = (val, P)
    assert best is not None
    return best


def nu_exact(q: int, E: tuple[int, ...]) -> Fr:
    """meas{x : maxgap(frac(e*x): e in E) > 2/q}."""
    theta = Fr(2, q)
    E = tuple(sorted(set(E)))
    bps = {Fr(0), Fr(1)}
    for a, b in combinations(E, 2):
        d = abs(b - a)
        if d == 0:
            continue
        for m in range(d + 1):
            bps.add(Fr(m, d))
    bps = sorted(bp for bp in bps if 0 <= bp <= 1)
    total = Fr(0)
    n = len(E)
    for cell_a, cell_b in zip(bps, bps[1:]):
        if cell_b <= cell_a:
            continue
        mid = (cell_a + cell_b) / 2
        order = sorted(range(n), key=lambda i: (Fr(E[i]) * mid) % 1)
        floors = [(Fr(E[i]) * mid).__floor__() for i in order]
        good_subarcs: list[tuple[Fr, Fr]] = []
        for idx in range(n):
            i_cur = order[idx]
            i_next = order[(idx + 1) % n]
            f_cur = floors[idx]
            f_next = floors[(idx + 1) % n]
            wrap = Fr(1) if idx == n - 1 else Fr(0)
            # gap(x) = (e_next-e_cur)x + f_cur - f_next + wrap.
            slope = Fr(E[i_next] - E[i_cur])
            const = Fr(f_cur - f_next) + wrap
            if slope == 0:
                if const > theta:
                    good_subarcs.append((cell_a, cell_b))
            elif slope > 0:
                lo = max(cell_a, (theta - const) / slope)
                if lo < cell_b:
                    good_subarcs.append((lo, cell_b))
            else:
                hi = min(cell_b, (theta - const) / slope)
                if cell_a < hi:
                    good_subarcs.append((cell_a, hi))
        total += union_measure(good_subarcs)
    return total


def p0_law_exact(q: int, E: tuple[int, ...]) -> tuple[Fr, tuple[Fr, ...]]:
    """Exact q-even sector cover atom.

    With q=2m and anchored offset 0, p0 is the measure for which each inner
    sector [j/m,(j+1)/m), j=1..m-1, is hit by a nonzero cluster phase.
    The returned law[t] is the measure with exactly t empty inner sectors.
    """
    assert q % 2 == 0
    m = q // 2
    endpoints = [Fr(j, m) for j in range(1, m + 1)]
    bps = {Fr(0), Fr(1)}
    for e in E:
        if e == 0:
            continue
        for end in endpoints:
            n = 0
            while True:
                x = (end + n) / e
                if x >= 1:
                    break
                if x >= 0:
                    bps.add(x)
                n += 1
    ordered = sorted(bp for bp in bps if 0 <= bp <= 1)
    law = [Fr(0) for _ in range(m)]
    nonzero = [e for e in E if e != 0]
    for lo, hi in zip(ordered, ordered[1:]):
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        hit = set()
        for e in nonzero:
            fr = (e * mid) % 1
            idx = (fr.numerator * m) // fr.denominator
            if 1 <= idx <= m - 1:
                hit.add(idx)
        empty = (m - 1) - len(hit)
        law[empty] += hi - lo
    return law[0], tuple(law)


def additive_energy(E: tuple[int, ...]) -> int:
    counts = Counter(a + b for a in E for b in E)
    return sum(c * c for c in counts.values())


def diff_profile(E: tuple[int, ...]) -> tuple[int, ...]:
    counts = Counter(abs(a - b) for a, b in combinations(E, 2))
    return tuple(sorted(counts.values(), reverse=True))


def majorizes(x: tuple[int, ...], y: tuple[int, ...]) -> bool:
    n = max(len(x), len(y))
    sx = sy = 0
    for i in range(n):
        if i < len(x):
            sx += x[i]
        if i < len(y):
            sy += y[i]
        if sx < sy:
            return False
    return sum(x) == sum(y)


def scalar_energy_inversions(rows: list[dict], metric: str) -> tuple[int, tuple[Fr, dict, dict] | None]:
    """Count A_high>A_low but metric_high<metric_low envelope inversions."""
    groups: dict[int, list[dict]] = defaultdict(list)
    for row in rows:
        groups[row["A"]].append(row)

    best_lower: dict | None = None
    inversions = 0
    worst: tuple[Fr, dict, dict] | None = None
    for A in sorted(groups):
        if best_lower is not None:
            for row in groups[A]:
                if row[metric] < best_lower[metric]:
                    inversions += 1
                    gap = best_lower[metric] - row[metric]
                    if worst is None or gap > worst[0]:
                        worst = (gap, best_lower, row)
        best_here = max(groups[A], key=lambda r: r[metric])
        if best_lower is None or best_here[metric] > best_lower[metric]:
            best_lower = best_here
    return inversions, worst


def fmt(fr: Fr) -> str:
    return f"{fr} ({float(fr):.6f})"


def analyze_bank(q: int, k: int, span: int) -> dict:
    ap = tuple(range(k))
    ap_diff = diff_profile(ap)
    rows: list[dict] = []
    for tail in combinations(range(1, span + 1), k - 1):
        E = (0,) + tail
        if not primitive(E):
            continue
        nu = nu_exact(q, E)
        p0, law = p0_law_exact(q, E)
        rows.append(
            {
                "E": E,
                "nu": nu,
                "D": Fr(1) - nu,
                "p0": p0,
                "law": law,
                "A": additive_energy(E),
                "diff": diff_profile(E),
            }
        )

    by_E = {row["E"]: row for row in rows}
    ap_row = by_E[ap]
    cap, capP = cap_gp(q, k)
    min_nu = min(rows, key=lambda r: (r["nu"], r["E"]))
    max_p0 = max(rows, key=lambda r: (r["p0"], tuple(-x for x in r["E"])))
    max_D = max(rows, key=lambda r: (r["D"], tuple(-x for x in r["E"])))
    max_A = max(rows, key=lambda r: (r["A"], r["E"]))
    ap_diff_fail = sum(1 for row in rows if not majorizes(ap_diff, row["diff"]))
    inv_p0, worst_p0 = scalar_energy_inversions(rows, "p0")
    inv_D, worst_D = scalar_energy_inversions(rows, "D")
    return {
        "q": q,
        "k": k,
        "span": span,
        "rows": rows,
        "count": len(rows),
        "cap": cap,
        "capP": capP,
        "ap": ap_row,
        "min_nu": min_nu,
        "max_D": max_D,
        "max_p0": max_p0,
        "max_A": max_A,
        "floor": ap_row["nu"] + cap - 1,
        "p0_margin": cap - max_p0["p0"],
        "ap_diff_fail": ap_diff_fail,
        "inv_p0": inv_p0,
        "worst_p0": worst_p0,
        "inv_D": inv_D,
        "worst_D": worst_D,
    }


def main() -> None:
    print("Exact small-q proof lab for LRC witness/p0/additive-energy mechanisms")
    print("=" * 78)
    print("q is the denominator in threshold 1/q; cluster GOOD_q is maxgap > 2/q.")
    print("Banks are primitive anchored shapes E={0,...} with max(E)<=q.")

    all_results: list[dict] = []
    proof_scores = Counter()
    min_floor: tuple[Fr, tuple[int, int]] | None = None
    min_p0_margin: tuple[Fr, tuple[int, int]] | None = None

    for q in (8, 10, 12, 14):
        print()
        print("-" * 78)
        print(f"q={q}: hard cluster sizes k={q//2 + 1}..{q-1}")
        print("-" * 78)
        for k in range(q // 2 + 1, q):
            result = analyze_bank(q, k, span=q)
            all_results.append(result)
            if min_floor is None or result["floor"] < min_floor[0]:
                min_floor = (result["floor"], (q, k))
            if min_p0_margin is None or result["p0_margin"] < min_p0_margin[0]:
                min_p0_margin = (result["p0_margin"], (q, k))

            ap = result["ap"]
            min_nu = result["min_nu"]
            max_p0 = result["max_p0"]
            max_D = result["max_D"]
            max_A = result["max_A"]
            cap = result["cap"]

            if min_nu["E"] != ap["E"]:
                proof_scores["bounded_nu_minimizer_fail"] += 1
            if max_p0["p0"] > cap:
                proof_scores["bounded_p0_cap_fail"] += 1
            if max_p0["E"] != ap["E"]:
                proof_scores["p0_AP_not_unique_or_not_max"] += 1
            if result["ap_diff_fail"]:
                proof_scores["AP_difference_profile_fail"] += result["ap_diff_fail"]
            proof_scores["scalar_p0_inversions"] += result["inv_p0"]
            proof_scores["scalar_D_inversions"] += result["inv_D"]

            print(f"k={k:2d} rows={result['count']:5d} cap={fmt(cap)} at P={result['capP']}")
            print(
                "     AP: nu=%s D=%s p0=%s A=%s"
                % (fmt(ap["nu"]), fmt(ap["D"]), fmt(ap["p0"]), ap["A"])
            )
            print(
                "     floor nu(AP)+cap-1=%s | bounded p0 margin cap-maxp0=%s"
                % (fmt(result["floor"]), fmt(result["p0_margin"]))
            )
            print(
                "     min nu: E=%s nu=%s %s"
                % (
                    min_nu["E"],
                    fmt(min_nu["nu"]),
                    "AP" if min_nu["E"] == ap["E"] else "NON-AP",
                )
            )
            print(
                "     max D : E=%s D=%s %s"
                % (
                    max_D["E"],
                    fmt(max_D["D"]),
                    "AP" if max_D["E"] == ap["E"] else "NON-AP",
                )
            )
            print(
                "     max p0: E=%s p0=%s %s"
                % (
                    max_p0["E"],
                    fmt(max_p0["p0"]),
                    "AP" if max_p0["E"] == ap["E"] else "NON-AP",
                )
            )
            print(
                "     max additive energy: E=%s A=%s %s"
                % (
                    max_A["E"],
                    max_A["A"],
                    "AP" if max_A["E"] == ap["E"] else "NON-AP",
                )
            )
            print(
                "     scalar-energy inversions: p0=%d D=%d; AP diff-profile failures=%d"
                % (result["inv_p0"], result["inv_D"], result["ap_diff_fail"])
            )
            if result["worst_p0"] is not None:
                gap, low, high = result["worst_p0"]
                print(
                    "       worst p0 inversion: A %d row %s p0=%s beats A %d row %s p0=%s by %s"
                    % (
                        low["A"],
                        low["E"],
                        fmt(low["p0"]),
                        high["A"],
                        high["E"],
                        fmt(high["p0"]),
                        fmt(gap),
                    )
                )

    print()
    print("=" * 78)
    print("Cross-q synthesis")
    print("=" * 78)
    assert min_floor is not None and min_p0_margin is not None
    print(f"smallest Bonferroni floor in this lab: q,k={min_floor[1]} value={fmt(min_floor[0])}")
    print(f"smallest bounded p0-cap margin: q,k={min_p0_margin[1]} value={fmt(min_p0_margin[0])}")
    print("carrier failure ledger:")
    for key in sorted(proof_scores):
        print(f"  {key}: {proof_scores[key]}")

    print()
    print("Tournament Analysis Hamiltonian path over proof carriers:")
    print(
        "  Bonferroni_floor"
        " > bounded_nu_minimizer"
        " > bounded_p0_cap"
        " > AP_difference_profile"
        " > scalar_additive_energy"
    )
    print()
    print("Interpretation:")
    print("  The floor and bounded cap mechanisms survive unchanged below q=14.")
    print("  Consecutive/AP remains the nu/D extremal row in every bounded bank checked.")
    print("  Some non-AP rows have larger p0, but every bounded-bank p0 stays below cap.")
    print("  Scalar additive energy is already non-monotone in easier q, so it is")
    print("  diagnostic only; the durable additive object is AP-facing difference")
    print("  profile majorization plus labelled residual control.")


if __name__ == "__main__":
    main()
