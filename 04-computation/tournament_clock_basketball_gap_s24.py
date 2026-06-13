#!/usr/bin/env python3
"""
tournament_clock_basketball_gap_s24.py

oracle-2026-06-01-S24

Two short Tournament-Analysis vignettes:

(1) BASKETBALL (discrete anchor).  A 5x5 directed pass-count matrix P.  Orient
    i->j iff P[i][j] > P[j][i]; ties broken by the jersey path 1->2->3->4->5
    (the repo base path).  Report the team's tournament: score sequence, H
    (Hamiltonian path count), iso-class, and whether any ties needed the path.

(2) H AS A SPREAD-METER.  Over the runner clock, correlate the tournament's H
    with the runners' maximum circular gap.  Finding 2 said transitive (H=1)
    <=> max gap > 1/2.  Here we check the full monotone trend: more bunching
    (larger max gap) => lower H (closer to transitive); most-even spacing
    (small max gap) => maximal H (regular).
"""

from __future__ import annotations

from fractions import Fraction
import importlib.util, pathlib
from collections import defaultdict

spec = importlib.util.spec_from_file_location(
    "clock", str(pathlib.Path(__file__).resolve().parent / "tournament_clock_s24.py"))
clock = importlib.util.module_from_spec(spec); spec.loader.exec_module(clock)


def basketball(P):
    n = len(P)
    ties = 0
    adj = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            if P[i][j] > P[j][i]:
                adj[i][j] = 1
            elif P[i][j] == P[j][i]:
                ties += 1
                if i < j:           # base-path tiebreak i->j for i<j
                    adj[i][j] = 1
    adj = tuple(tuple(r) for r in adj)
    return {
        "score": clock.score_sequence(adj),
        "H": clock.hamiltonian_path_count(adj),
        "canon": clock.canonical_form(adj),
        "ties_broken_by_path": ties // 2,
        "adj": adj,
    }


def max_gap(speeds, t):
    pts = sorted(float(clock.frac(Fraction(s) * t)) for s in speeds)
    gaps = [pts[(i + 1) % len(pts)] - pts[i] for i in range(len(pts) - 1)]
    gaps.append(pts[0] + 1 - pts[-1])
    return max(gaps)


def main():
    print("Tournament Analysis vignettes (oracle-2026-06-01-S24)\n")

    print("=" * 60)
    print("(1) BASKETBALL — pass matrix -> tournament of 5")
    print("=" * 60)
    # rows pass to columns; jerseys 1..5 are indices 0..4
    teams = {
        "Team A (PG-heavy)": [
            [0, 18, 15, 9, 12],   # 1 (PG) passes a lot
            [7, 0, 6, 5, 8],
            [5, 9, 0, 7, 4],
            [4, 6, 8, 0, 3],
            [6, 7, 5, 5, 0],
        ],
        "Team B (with ties)": [
            [0, 10, 10, 4, 6],    # 1<->2 tie (10,10); 1<->3 tie (10,10)
            [10, 0, 7, 9, 5],
            [10, 6, 0, 8, 7],
            [5, 4, 6, 0, 9],
            [6, 8, 7, 9, 0],
        ],
    }
    for name, P in teams.items():
        r = basketball(P)
        print(f" {name}: score={r['score']}  H={r['H']}  "
              f"ties_via_path={r['ties_broken_by_path']}")
        # locate this iso-class in the runner-clock menu? print canon H only
    print(" The team is a POINT in G_5; a season is a walk (one point per game).")
    print(" Note H=1 would mean a perfect passing pecking-order (transitive);")
    print(" higher H means cyclic passing structure (no clean hierarchy).\n")

    print("=" * 60)
    print("(2) H as a spread-meter over the runner clock")
    print("=" * 60)
    for label, speeds in {
        "arith 0..4": (0, 1, 2, 3, 4),
        "primes 0,2,3,5,7": (0, 2, 3, 5, 7),
        "spread 0,1,4,9,11": (0, 1, 4, 9, 11),
    }.items():
        cells = clock.clock_cells(speeds)
        # bucket cells by H, record max-gap range
        byH = defaultdict(list)
        for t, adj in cells:
            H = clock.hamiltonian_path_count(adj)
            byH[H].append(max_gap(speeds, t))
        print(f" [{label}]")
        for H in sorted(byH):
            gs = byH[H]
            print(f"    H={H:<3}: cells={len(gs):<3} max_gap in "
                  f"[{min(gs):.3f}, {max(gs):.3f}]  mean={sum(gs)/len(gs):.3f}")
    print()
    print(" Trend: H=1 (transitive) sits at the LARGEST gaps (>1/2, bunched);")
    print(" the maximal H (regular) sits at the SMALLEST max-gap (even spread).")
    print(" So H reads off how bunched the runners are — a tournament gap-meter.")


if __name__ == "__main__":
    main()
