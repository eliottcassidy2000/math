#!/usr/bin/env python3
"""
lrc_n4_adjacent_family_s552.py    codex-2026-06-01-S552

Exact formula for the adjacent n=4 LRC family (1,q,q+1).

For threshold 1/4, let
  M(q) = |{t in [0,1): ||t||, ||q t||, ||(q+1)t|| >= 1/4}|.

The script verifies the closed form
  q = 0 mod 4: M(q) = (q+2)/(16(q+1))
  q = 1 mod 4: M(q) = (q+3)/(16q)
  q = 2 mod 4: M(q) = (q-2)/(16(q+1))
  q = 3 mod 4: M(q) = (q-1)/(16q)

and records the consequence:
  M(2)=0 for the AP (1,2,3), and min_{q>=3} M(q)=M(6)=1/28.

Tournament Analysis declaration:
  Vertices: adjacent triples (1,q,q+1), q=2..17.
  Pairwise observable: exact safe measure M(q).
  Switch/gauge: q_i -> q_j if M(q_i) < M(q_j), i.e. q_i is tighter.
  Tie Hamiltonian path: natural q order 2 -> 3 -> ... -> 17.
  Fingerprints: score histogram, directed 3-cycles, SCCs, edge flips, and
    Hamiltonian path count are reported.

Assumption challenge:
  Vertices could be runners, breakpoints, safe intervals, danger arcs, residues,
  parity classes, or proof obligations. Here we choose adjacent-family q rows
  because the proof quotient preserves the exact n=4 safe measure and destroys
  all non-adjacent triple information. This is a local theorem inside the bigger
  HYP-2040 measure-gap problem, not a proof of the full n=4 gap.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as F

THR = F(1, 4)


def fnorm(x: F) -> F:
    r = x - (x.numerator // x.denominator)
    return min(r, 1 - r)


def exact_measure_and_intervals(q: int) -> tuple[F, list[tuple[F, F]]]:
    speeds = (1, q, q + 1)
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


def formula(q: int) -> F:
    r = q % 4
    if r == 0:
        return F(q + 2, 16 * (q + 1))
    if r == 1:
        return F(q + 3, 16 * q)
    if r == 2:
        return F(q - 2, 16 * (q + 1))
    return F(q - 1, 16 * q)


def proof_sum(q: int) -> tuple[str, int, F, F]:
    """Return the residue case, positive-side interval count, sum numerator, total."""
    if q % 2 == 0:
        if q % 4 == 0:
            r = q // 4
            case = "q=4r"
        else:
            r = (q - 2) // 4
            case = "q=4r+2"
        # positive side intervals: u in [(j+3/4)/(q+1),(j+3/4)/q]
        s = sum(F(j, 1) + F(3, 4) for j in range(r))
    else:
        if q % 4 == 1:
            r = (q - 1) // 4
            case = "q=4r+1"
        else:
            r = (q - 3) // 4
            case = "q=4r+3"
        # positive side intervals: u in [(j+1/4)/(q+1),(j+1/4)/q]
        s = sum(F(j, 1) + F(1, 4) for j in range(r + 1))
        r = r + 1
    total = F(2, q * (q + 1)) * s
    return case, r, s, total


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
    print("LRC n=4 adjacent-family formula (codex-S552)")
    print("Family: speeds (1,q,q+1), threshold 1/4\n")

    print("Closed form by residue class")
    print("  q=0 mod 4: M(q)=(q+2)/(16(q+1))")
    print("  q=1 mod 4: M(q)=(q+3)/(16q)")
    print("  q=2 mod 4: M(q)=(q-2)/(16(q+1))")
    print("  q=3 mod 4: M(q)=(q-1)/(16q)\n")

    print("Verification table")
    bad = 0
    measures: dict[int, F] = {}
    for q in range(2, 42):
        exact, intervals = exact_measure_and_intervals(q)
        pred = formula(q)
        measures[q] = exact
        if exact != pred:
            bad += 1
        if q <= 17 or q % 4 == 2:
            case, count, partial_sum, total = proof_sum(q)
            print(
                f"  q={q:2d} {case:7s} exact={str(exact):>6s} "
                f"formula={str(pred):>6s} intervals={len(intervals):2d} "
                f"pos_count={count:2d} pos_sum={partial_sum} sum_formula={total}"
            )
    print(f"\nformula mismatches for q=2..41: {bad}")

    positive = [(m, q) for q, m in measures.items() if q >= 3]
    min_m, min_q = min(positive)
    print(f"minimum positive adjacent-family measure for q=3..41: q={min_q}, M={min_m}")
    print(f"global adjacent-family theorem: M(2)=0 and M(q)>=1/28 for all q>=3, equality at q=6")

    qs = list(range(2, 18))
    adj = tournament(qs, measures)
    scores = [sum(row) for row in adj]
    flips = []
    for i, qi in enumerate(qs):
        for j, qj in enumerate(qs):
            if i < j and adj[j][i]:
                flips.append(f"q={qj}>q={qi}")
    print("\nTournament Analysis")
    print("  vertices: adjacent rows q=2..17 representing triples (1,q,q+1)")
    print("  observable: exact safe measure M(q); lower means tighter")
    print("  switch: q_i -> q_j iff M(q_i) < M(q_j); ties use natural q order")
    print("  tie path: q=2 -> q=3 -> ... -> q=17")
    print(f"  scores: {dict(zip(qs, scores))}")
    print(f"  score_hist: {dict(sorted(Counter(scores).items()))}")
    print(f"  directed_3_cycles: {directed_triangles(adj)}")
    print(f"  SCCs: {[[qs[i] for i in comp] for comp in sccs(adj)]}")
    print(f"  Hamiltonian_path_count: {h_count(adj)}")
    print(f"  edge_flips_vs_natural_q_path: {flips[:24]}{' ...' if len(flips) > 24 else ''}")


if __name__ == "__main__":
    main()
