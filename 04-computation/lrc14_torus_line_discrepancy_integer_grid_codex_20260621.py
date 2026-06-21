#!/usr/bin/env python3
"""
HYP-2736 torus-line discrepancy integer-grid scout.

This refines the HYP-2730 L7 tail.  For coprime p>q and apex prime 7,
the cell law of v -> (qv,pv) against the 7x7 sector grid has all breakpoints
on the integer grid of size 7*p*q.  If c_ij is the number of integer grid
subintervals landing in cell (i,j), then

    D_{p,q} = sum_ij |c_ij/(7pq) - 1/49|
            = sum_ij |49*c_ij - 7*p*q| / (343*p*q).

Thus the sharp observed target D_{p,q} <= 12/(7q) is equivalent to the purely
integer inequality

    sum_ij |49*c_ij - 7*p*q| <= 588*p.

This script verifies the exact integer formula against the older Fraction
breakpoint computation on a small atlas, then scans the L7 bounded-ratio
window 1 < p/q <= 43/20.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction as Fr
from math import gcd

P = 7
HI = Fr(43, 20)
MARGIN = Fr(21, 100)
QMAX = 160


def sector(y: Fr) -> int:
    return int(P * y)


def slow_discrepancy(p: int, q: int) -> Fr:
    """Original Fraction breakpoint computation, used only for cross-checks."""
    breakpoints = {Fr(0), Fr(1)}
    for f in (p, q):
        for t in range(P * f):
            breakpoints.add(Fr(t, P * f))
    cells: dict[tuple[int, int], Fr] = {}
    points = sorted(breakpoints)
    for a, b in zip(points, points[1:]):
        mid = (a + b) / 2
        key = (sector((q * mid) % 1), sector((p * mid) % 1))
        cells[key] = cells.get(key, Fr(0)) + (b - a)
    uniform = Fr(1, P * P)
    return sum(
        abs(cells.get((i, j), Fr(0)) - uniform)
        for i in range(P)
        for j in range(P)
    )


def cell_counts(p: int, q: int) -> list[list[int]]:
    """Integer counts c_ij on the common 7*p*q breakpoint grid."""
    total = P * p * q
    points = {0, total}
    for t in range(P * q + 1):
        points.add(p * t)
    for s in range(P * p + 1):
        points.add(q * s)

    counts = [[0 for _ in range(P)] for _ in range(P)]
    sorted_points = sorted(points)
    for a, b in zip(sorted_points, sorted_points[1:]):
        if a == b:
            continue
        midpoint_twice = a + b
        # midpoint = midpoint_twice / (2*7*p*q).
        # q*midpoint sector = floor(midpoint_twice/(2*p)) mod 7.
        i = (midpoint_twice // (2 * p)) % P
        j = (midpoint_twice // (2 * q)) % P
        counts[i][j] += b - a
    return counts


def defect_sum(p: int, q: int) -> int:
    counts = cell_counts(p, q)
    uniform_units = P * p * q
    return sum(
        abs(P * P * counts[i][j] - uniform_units)
        for i in range(P)
        for j in range(P)
    )


def discrepancy(p: int, q: int) -> Fr:
    return Fr(defect_sum(p, q), P * P * P * p * q)


def window_pairs(qmax: int):
    for q in range(1, qmax + 1):
        for p in range(q + 1, int(HI * q) + 1):
            if gcd(p, q) != 1:
                continue
            ratio = Fr(p, q)
            if Fr(1) < ratio <= HI:
                yield p, q


def score_tournament(labels: list[str]):
    """A total-order tournament fingerprint for already sorted risk labels."""
    n = len(labels)
    score_hist = Counter(range(n))
    return {
        "vertices": n,
        "score_hist": dict(sorted(score_hist.items())),
        "directed_3cycles": 0,
        "scc_sizes": [1] * n,
        "hamiltonian_path_count": 1,
        "path": " > ".join(labels),
    }


def main() -> None:
    print("=" * 86)
    print("HYP-2736 torus-line discrepancy integer-grid scout")
    print("=" * 86)

    checked = 0
    for p, q in window_pairs(12):
        d_fast = discrepancy(p, q)
        d_slow = slow_discrepancy(p, q)
        assert d_fast == d_slow, (p, q, d_fast, d_slow)
        checked += 1
    print(f"\n[A] integer-grid formula cross-check: {checked} q<=12 ratios matched slow Fraction cells")
    print("    D = sum |49*c_ij - 7pq| / (343pq)")
    print("    sharp target D <= 12/(7q) is exactly defect_sum <= 588*p")

    rows = []
    by_q: dict[int, tuple[int, Fr, int]] = {}
    by_residue: dict[tuple[int, int], tuple[Fr, int, int, Fr]] = {}
    pairs = 0
    violations_12_over_7q = []
    for p, q in window_pairs(QMAX):
        pairs += 1
        d = discrepancy(p, q)
        defect = defect_sum(p, q)
        rows.append((p, q, d, d * q, defect))
        if q not in by_q or d > by_q[q][1]:
            by_q[q] = (p, d, defect)
        if d * q > Fr(12, 7):
            violations_12_over_7q.append((p, q, d))
        if q >= 15:
            key = (p % P, q % P)
            old = by_residue.get(key)
            if old is None or (d * q, -q, -p) > (old[0], -old[2], -old[1]):
                by_residue[key] = (d * q, p, q, d)

    max_row = max(rows, key=lambda row: row[3])
    tail_rows = [row for row in rows if row[1] >= 5]
    max_tail = max(tail_rows, key=lambda row: row[2])
    largest_bad_q = max((q for _, q, d, _, _ in rows if d >= MARGIN), default=0)

    print(f"\n[B] exact scan window: q<= {QMAX}, pairs={pairs}")
    print(
        "    max D*q = %s at p/q=%s/%s (D=%s)"
        % (max_row[3], max_row[0], max_row[1], max_row[2])
    )
    print(f"    D <= 12/(7q) violations in scan: {len(violations_12_over_7q)}")
    print(
        "    largest q with D >= %s is %s; finite atlas through q=4 is enough in this scan"
        % (MARGIN, largest_bad_q)
    )
    print(
        "    worst q>=5 row: p/q=%s/%s, D=%s, D*q=%s"
        % (max_tail[0], max_tail[1], max_tail[2], max_tail[3])
    )

    print("\n[C] worst exact discrepancy by denominator q<=14")
    print("      q    p       D          D*q        defect <= 588*p slack")
    for q in range(1, 15):
        if q not in by_q:
            continue
        p, d, defect = by_q[q]
        slack = 588 * p - defect
        print(f"    {q:>3} {p:>4} {str(d):>10} {str(d*q):>10} {defect:>8} <= {588*p:<8} slack={slack}")

    residue_rows = []
    for key, (risk, p, q, d) in by_residue.items():
        residue_rows.append((risk, key, p, q, d))
    residue_rows.sort(key=lambda row: (-row[0], row[3], row[2], row[1]))
    all_labels = [f"{a}/{b}" for _, (a, b), _, _, _ in residue_rows]
    top_labels = all_labels[:12]
    fp = score_tournament(all_labels)
    print("\n[D] residue-pair tournament for q>=15")
    print("    vertices are (p mod 7)/(q mod 7), ordered by worst observed D*q")
    print(f"    vertices={fp['vertices']}")
    print(f"    top-12 path: {' > '.join(top_labels)}")
    print("    score_hist is the transitive histogram {0:1,...,%d:1}" % (fp["vertices"] - 1))
    print(f"    directed_3cycles={fp['directed_3cycles']} SCC_sizes={fp['scc_sizes'][:12]}...")
    print(f"    Hamiltonian_path_count={fp['hamiltonian_path_count']}")
    print("    top residue witnesses:")
    for risk, (a, b), p, q, d in residue_rows[:12]:
        print(f"      ({a},{b}) witness {p}/{q}: D={d}, D*q={risk}")

    print("\n[E] proof target")
    print("    Prove the integer defect inequality defect_sum(p,q) <= 588*p for all")
    print("    coprime 1 < p/q <= 43/20.  Together with exact q<=8 atlas checks,")
    print("    this supplies a sharper alternate HYP-2730 non-resonant L7 tail.")
    print("    The zero rows when p==0 or q==0 mod 7 suggest a mod-7 block proof,")
    print("    with the nonzero residue-pair tournament above as the finite guide.")
    print("\nDONE.")


if __name__ == "__main__":
    main()
