#!/usr/bin/env python3
"""
lrc_n4_additive_return_s553.py    codex-2026-06-01-S553

Exact formula for primitive additive-return n=4 LRC triples (a,b,a+b).

For threshold 1/4 and coprime positive integers a,b, define
  M(a,b) = |{t in [0,1): ||a t||, ||b t||, ||(a+b)t|| >= 1/4}|.

The script verifies the closed form used in THM-393:

  If a,b are odd:
    a+b = 0 mod 4: M(a,b) = (ab - 1)/(16ab)
    a+b = 2 mod 4: M(a,b) = (ab + 3)/(16ab)

  If exactly one of a,b is even, let e be the even generator, o the odd
  generator, and c=a+b:
    e = 0 mod 4: M(a,b) = (oc + 1)/(16oc)
    e = 2 mod 4: M(a,b) = (oc - 3)/(16oc)

Consequences:
  M(a,b)=0 only for {a,b}={1,2}, i.e. speeds (1,2,3).
  Every other primitive additive-return triple has M(a,b) >= 1/28.
  Equality occurs only for {a,b}={1,6}, i.e. speeds (1,6,7).

Tournament Analysis declaration:
  Vertices: selected primitive generator pairs (a,b), representing triples
    (a,b,a+b).
  Pairwise observable: exact safe measure M(a,b).
  Switch/gauge: v_i -> v_j iff M(v_i) < M(v_j), so the tighter additive
    return wins.
  Tie Hamiltonian path: lexicographic generator-pair order.
  Fingerprints: score histogram, directed 3-cycles, SCCs, edge flips, and
    Hamiltonian-path count are reported.

Assumption challenge:
  Vertices could be runners, generator pairs, sum triples, branch centers,
  endpoint obligations, residues, parity classes, or proof intervals.  Here
  vertices are generator pairs because the quotient preserves the additive
  relation c=a+b and exact safe measure, but destroys non-additive triples.
  The theorem is a local HYP-2040 slice, not a proof of the full n=4 gap.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from math import gcd

THR = F(1, 4)


def fnorm(x: F) -> F:
    r = x - (x.numerator // x.denominator)
    return min(r, 1 - r)


def exact_measure_and_intervals(speeds: tuple[int, int, int]) -> tuple[F, list[tuple[F, F]]]:
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


def formula(a: int, b: int) -> F:
    assert gcd(a, b) == 1
    if a % 2 == 1 and b % 2 == 1:
        prod = a * b
        if (a + b) % 4 == 0:
            return F(prod - 1, 16 * prod)
        return F(prod + 3, 16 * prod)

    e = a if a % 2 == 0 else b
    o = b if a % 2 == 0 else a
    c = a + b
    prod = o * c
    if e % 4 == 0:
        return F(prod + 1, 16 * prod)
    return F(prod - 3, 16 * prod)


def branch_sum_measure(a: int, b: int) -> F:
    """Exact moving-window sum from the speed-a safe arcs.

    On the h-th speed-a safe arc, write t=(2h+1)/(2a)+u with
    |u|<=1/(4a).  If u>=0 then {bt} must lie in [3/4-au,3/4]; if u<=0
    and w=-u then {bt} must lie in [1/4,1/4+aw].  This function sums those
    affine intervals exactly.  It is intentionally independent of the closed
    form, and acts as a second proof-shape verifier.
    """
    c = a + b
    radius = F(1, 4 * a)
    total = F(0)

    # Enough j values for all affine windows that can meet [0,radius].
    # The largest useful numerator is < c/(4a), so this small exact bound is
    # harmless even for the asymmetric rows in the verification range.
    max_j = c // a + 4

    for h in range(a):
        beta = F((b * (2 * h + 1)) % (2 * a), 2 * a)
        for j in range(-2, max_j + 1):
            alpha = F(j, 1) + F(3, 4) - beta
            lo = alpha / c
            hi = alpha / b
            overlap = min(hi, radius) - max(lo, F(0))
            if overlap > 0:
                total += overlap

            gamma = beta + F(j, 1) - F(1, 4)
            lo = gamma / c
            hi = gamma / b
            overlap = min(hi, radius) - max(lo, F(0))
            if overlap > 0:
                total += overlap

    return total


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


def tournament(vertices: list[tuple[int, int]], measures: dict[tuple[int, int], F]) -> list[list[int]]:
    n = len(vertices)
    adj = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            vi, vj = vertices[i], vertices[j]
            if measures[vi] == measures[vj]:
                winner = i
            else:
                winner = i if measures[vi] < measures[vj] else j
            loser = j if winner == i else i
            adj[winner][loser] = 1
    return adj


def case_name(a: int, b: int) -> str:
    if a % 2 == 1 and b % 2 == 1:
        return "odd/odd c=0mod4" if (a + b) % 4 == 0 else "odd/odd c=2mod4"
    e = a if a % 2 == 0 else b
    return "even e=0mod4" if e % 4 == 0 else "even e=2mod4"


def main() -> None:
    print("LRC n=4 additive-return formula (codex-S553)")
    print("Family: primitive generator pairs (a,b), speeds (a,b,a+b), threshold 1/4\n")

    print("Closed form")
    print("  a,b odd and a+b=0 mod 4: M=(ab-1)/(16ab)")
    print("  a,b odd and a+b=2 mod 4: M=(ab+3)/(16ab)")
    print("  one even generator e, odd generator o, c=a+b, e=0 mod 4: M=(oc+1)/(16oc)")
    print("  one even generator e, odd generator o, c=a+b, e=2 mod 4: M=(oc-3)/(16oc)\n")

    samples = [
        (1, 2),
        (1, 6),
        (1, 3),
        (2, 3),
        (2, 5),
        (3, 5),
        (3, 10),
        (1, 10),
        (4, 1),
        (4, 5),
        (5, 9),
        (7, 11),
    ]
    print("Sample table")
    for a, b in samples:
        exact, intervals = exact_measure_and_intervals((a, b, a + b))
        branch = branch_sum_measure(a, b)
        pred = formula(a, b)
        print(
            f"  (a,b)=({a:2d},{b:2d}) speeds={(a,b,a+b)!s:>12s} "
            f"{case_name(a,b):18s} exact={str(exact):>6s} "
            f"formula={str(pred):>6s} branch={str(branch):>6s} intervals={len(intervals):2d}"
        )

    bad_exact = 0
    bad_branch = 0
    minima: list[tuple[F, tuple[int, int]]] = []
    for a in range(1, 61):
        for b in range(1, 61):
            if gcd(a, b) != 1:
                continue
            exact, _ = exact_measure_and_intervals((a, b, a + b))
            pred = formula(a, b)
            branch = branch_sum_measure(a, b)
            if exact != pred:
                bad_exact += 1
            if branch != pred:
                bad_branch += 1
            minima.append((exact, tuple(sorted((a, b)))))
    minima = sorted(set(minima))
    print(f"\nformula mismatches for ordered coprime a,b<=60: {bad_exact}")
    print(f"branch-sum mismatches for ordered coprime a,b<=60: {bad_branch}")
    print("smallest additive-return measures")
    for m, pair in minima[:12]:
        print(f"  generators={pair}, speeds={tuple(sorted((pair[0], pair[1], pair[0]+pair[1])))} M={m}")

    positive = [(m, pair) for m, pair in minima if m > 0]
    min_m, min_pair = min(positive)
    print(f"\nminimum positive additive-return measure in box: generators={min_pair}, M={min_m}")
    print("theorem consequence: only zero is generators {1,2}; positive minimum is {1,6}, M=1/28")

    vertices = [
        (1, 2),
        (1, 6),
        (1, 3),
        (2, 3),
        (2, 5),
        (3, 5),
        (3, 10),
        (1, 10),
        (4, 1),
        (4, 5),
        (5, 9),
        (7, 11),
    ]
    measures = {v: formula(*v) for v in vertices}
    adj = tournament(vertices, measures)
    scores = [sum(row) for row in adj]
    flips = []
    for i, vi in enumerate(vertices):
        for j, vj in enumerate(vertices):
            if i < j and adj[j][i]:
                flips.append(f"{vj}->{vi}")
    print("\nTournament Analysis")
    print("  vertices: selected additive-return generator pairs (a,b)")
    print("  observable: exact safe measure M(a,b); lower means tighter")
    print("  switch: v_i -> v_j iff M(v_i) < M(v_j); ties use lexicographic path")
    print(f"  tie path: {' -> '.join(str(v) for v in vertices)}")
    print(f"  scores: {dict(zip(vertices, scores))}")
    print(f"  score_hist: {dict(sorted(Counter(scores).items()))}")
    print(f"  directed_3_cycles: {directed_triangles(adj)}")
    print(f"  SCCs: {[[vertices[i] for i in comp] for comp in sccs(adj)]}")
    print(f"  Hamiltonian_path_count: {h_count(adj)}")
    print(f"  edge_flips_vs_tie_path: {flips[:24]}{' ...' if len(flips) > 24 else ''}")


if __name__ == "__main__":
    main()
