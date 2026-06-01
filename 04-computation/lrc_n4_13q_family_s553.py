#!/usr/bin/env python3
"""
lrc_n4_13q_family_s553.py    codex-2026-06-01-S553

Exact formula for the n=4 LRC family (1,3,q).

For threshold 1/4, define
  M(q) = |{t in [0,1): ||t||, ||3t||, ||q t|| >= 1/4}|.

The first two speeds force t into [5/12,7/12] up to endpoints.  Writing
t=1/2+u, |u|<=1/12, reduces the formula to a one-runner count over
[-q/12,q/12].  With q=12m+r, the closed form can be written as M(q)=N_r(q)/(12q):

  r:      0    1      2      3      4      5    6    7      8      9      10     11
  N_r:    q    q+1    q-2    q+3    q-2    q+1  q    q-1    q+2    q-3    q+2    q-1

Consequences for distinct triples:
  q=2 is the AP (1,2,3) after sorting, and M(2)=0.
  q=4 is the adjacent/additive row (1,3,4), and M(4)=1/24.
  For q not in {1,2,3,4}, M(q)>=1/18, equality first at q=9.

Tournament Analysis declaration:
  Vertices: selected q rows representing triples (1,3,q).
  Pairwise observable: exact safe measure M(q).
  Switch/gauge: q_i -> q_j iff M(q_i) < M(q_j), so tighter rows win.
  Tie Hamiltonian path: natural q order in the selected list.
  Fingerprints: score histogram, directed 3-cycles, SCCs, edge flips, and
    Hamiltonian-path count are reported.

Assumption challenge:
  Vertices could be q rows, the forced interval [5/12,7/12], residual q mod 12
  classes, safe subintervals, endpoints, or proof obligations.  Here q rows
  preserve the exact family measure and the next-obstruction ordering, but
  destroy all triples not containing the pair (1,3).
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as F

THR = F(1, 4)


def fnorm(x: F) -> F:
    r = x - (x.numerator // x.denominator)
    return min(r, 1 - r)


def exact_measure_and_intervals(q: int) -> tuple[F, list[tuple[F, F]]]:
    speeds = (1, 3, q)
    bps = {F(0), F(1)}
    for s in speeds:
        for num in range(1, 4 * s):
            if num % 4 in (1, 3):
                bps.add(F(num, 4 * s))
    pts = sorted(bps)
    measure = F(0)
    intervals: list[tuple[F, F]] = []
    for lo, hi in zip(pts, pts[1:]):
        mid = (lo + hi) / 2
        if min(fnorm(F(s) * mid) for s in speeds) >= THR:
            intervals.append((lo, hi))
            measure += hi - lo
    return measure, intervals


def numerator(q: int) -> int:
    r = q % 12
    if r in (0, 6):
        return q
    if r in (1, 5):
        return q + 1
    if r in (2, 4):
        return q - 2
    if r == 3:
        return q + 3
    if r in (7, 11):
        return q - 1
    if r in (8, 10):
        return q + 2
    if r == 9:
        return q - 3
    raise AssertionError("unreachable residue")


def formula(q: int) -> F:
    return F(numerator(q), 12 * q)


def residual_length(q: int) -> F:
    """Length in x-space of the residual piece after full periods."""
    r = q % 12
    if q % 2 == 0:
        table = {
            0: F(0),
            2: F(0),
            4: F(1, 6),
            6: F(1, 2),
            8: F(5, 6),
            10: F(1),
        }
    else:
        table = {
            1: F(1, 6),
            3: F(1, 2),
            5: F(1, 2),
            7: F(1, 2),
            9: F(1, 2),
            11: F(5, 6),
        }
    return table[r]


def interval_count_formula(q: int) -> F:
    m = q // 12
    return (F(m) + residual_length(q)) / q


def h_count(adj: list[list[int]]) -> int:
    n = len(adj)
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for v in range(n):
            c = dp[mask][v]
            if not c:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v][u]:
                    dp[mask | (1 << u)][u] += c
    return sum(dp[full])


def sccs(adj: list[list[int]]) -> list[list[int]]:
    n = len(adj)
    seen = [False] * n
    order: list[int] = []

    def dfs(v: int) -> None:
        seen[v] = True
        for u, edge in enumerate(adj[v]):
            if edge and not seen[u]:
                dfs(u)
        order.append(v)

    for v in range(n):
        if not seen[v]:
            dfs(v)

    radj = [[adj[j][i] for j in range(n)] for i in range(n)]
    seen = [False] * n
    comps: list[list[int]] = []

    def rdfs(v: int, comp: list[int]) -> None:
        seen[v] = True
        comp.append(v)
        for u, edge in enumerate(radj[v]):
            if edge and not seen[u]:
                rdfs(u, comp)

    for v in reversed(order):
        if not seen[v]:
            comp: list[int] = []
            rdfs(v, comp)
            comps.append(sorted(comp))
    return comps


def directed_triangles(adj: list[list[int]]) -> int:
    n = len(adj)
    out = 0
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                cyclic = adj[i][j] + adj[j][k] + adj[k][i]
                if cyclic in (0, 3):
                    out += 1
    return out


def tournament(qs: list[int], measures: dict[int, F]) -> list[list[int]]:
    n = len(qs)
    adj = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            qi, qj = qs[i], qs[j]
            if measures[qi] == measures[qj]:
                winner = i
            else:
                winner = i if measures[qi] < measures[qj] else j
            loser = j if winner == i else i
            adj[winner][loser] = 1
    return adj


def main() -> None:
    print("LRC n=4 (1,3,q) family formula (codex-S553)")
    print("Family: speeds (1,3,q), threshold 1/4\n")

    print("Closed form: q=12m+r and M(q)=N_r(q)/(12q)")
    print("  r= 0: N=q      r= 1: N=q+1    r= 2: N=q-2    r= 3: N=q+3")
    print("  r= 4: N=q-2    r= 5: N=q+1    r= 6: N=q      r= 7: N=q-1")
    print("  r= 8: N=q+2    r= 9: N=q-3    r=10: N=q+2    r=11: N=q-1\n")

    bad_formula = 0
    bad_interval = 0
    measures: dict[int, F] = {}
    print("Verification table")
    for q in range(1, 81):
        exact, intervals = exact_measure_and_intervals(q)
        pred = formula(q)
        interval_pred = interval_count_formula(q)
        measures[q] = exact
        if exact != pred:
            bad_formula += 1
        if interval_pred != pred:
            bad_interval += 1
        if q <= 40 or exact <= F(1, 12):
            print(
                f"  q={q:2d} r={q%12:2d} exact={str(exact):>6s} "
                f"formula={str(pred):>6s} interval_count={str(interval_pred):>6s} "
                f"intervals={len(intervals):2d}"
            )
    print(f"\nformula mismatches for q=1..80: {bad_formula}")
    print(f"interval-count mismatches for q=1..80: {bad_interval}")

    distinct = [(m, q) for q, m in measures.items() if q not in (1, 3)]
    print("smallest (1,3,q) measures with distinct triples")
    for m, q in sorted(distinct)[:12]:
        print(f"  q={q:2d}, speeds={tuple(sorted((1,3,q)))}, M={m}")
    residual = [(m, q) for q, m in measures.items() if q not in (1, 2, 3, 4)]
    min_m, min_q = min(residual)
    print(f"\nminimum after excluding duplicate/AP/adjacent-additive rows q in {{1,2,3,4}}: q={min_q}, M={min_m}")

    qs = [2, 4, 9, 7, 14, 11, 16, 6, 13, 5, 8, 10]
    selected = {q: measures[q] for q in qs}
    adj = tournament(qs, selected)
    scores = [sum(row) for row in adj]
    flips = []
    for i, qi in enumerate(qs):
        for j, qj in enumerate(qs):
            if i < j and adj[j][i]:
                flips.append(f"q={qj}->q={qi}")

    print("\nTournament Analysis")
    print("  vertices: selected q rows for triples (1,3,q)")
    print("  observable: exact safe measure M(q); lower means tighter")
    print("  switch: q_i -> q_j iff M(q_i) < M(q_j); ties use selected q order")
    print(f"  tie path: {' -> '.join('q='+str(q) for q in qs)}")
    print(f"  scores: {dict(zip(qs, scores))}")
    print(f"  score_hist: {dict(sorted(Counter(scores).items()))}")
    print(f"  directed_3_cycles: {directed_triangles(adj)}")
    print(f"  SCCs: {[[qs[i] for i in comp] for comp in sccs(adj)]}")
    print(f"  Hamiltonian_path_count: {h_count(adj)}")
    print(f"  edge_flips_vs_tie_path: {flips[:24]}{' ...' if len(flips) > 24 else ''}")


if __name__ == "__main__":
    main()
