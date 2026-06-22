#!/usr/bin/env python3
"""
S114 lcm-committed denominator wall.

Family:
    S_X = {1,2,...,11,13,lcm(2,...,X)}.

For every D <= X, the committed speed is 0 mod D, so no reduced fraction a/D
can be a 1/14 lonely witness.  This is a rigorous obstruction to any fixed
finite denominator basis.  The first witness above the wall is not always
nextprime(X): the observed control is "first denominator whose prime-power
content is not fully swallowed by lcm(2..X)" plus compatibility with the base
row {1..11,13}.

Tournament Analysis:
    vertices are proof routes/lemmas rather than runners.  The pairwise
    observable is the number of S114 criteria each route satisfies:
    unbounded-denominator wall, compatibility with the lcm family, preservation
    of the finite/analytic node split, formalizability, and fit to the first-D
    data.  Ties use the declared Hamiltonian order.
"""

from __future__ import annotations

from collections import Counter, deque
from math import gcd, lcm


BASE = list(range(1, 12)) + [13]


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def next_prime(n: int) -> int:
    m = n + 1
    while not is_prime(m):
        m += 1
    return m


def factor(n: int) -> str:
    pieces = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            e = 0
            while n % d == 0:
                n //= d
                e += 1
            pieces.append(f"{d}^{e}" if e > 1 else str(d))
        d += 1 if d == 2 else 2
    if n > 1:
        pieces.append(str(n))
    return "*".join(pieces) if pieces else "1"


def lcm_upto(x: int) -> int:
    out = 1
    for d in range(2, x + 1):
        out = lcm(out, d)
    return out


def lonely_at_denominator(speeds: list[int], d: int, a: int) -> bool:
    """Inclusive 1/14 threshold, using exact integer arithmetic."""
    return all(14 * min((s * a) % d, d - ((s * a) % d)) >= d for s in speeds)


def witnesses(speeds: list[int], d: int) -> list[int]:
    return [
        a
        for a in range(1, d)
        if gcd(a, d) == 1 and lonely_at_denominator(speeds, d, a)
    ]


def first_witness_denominator(x: int, scan_extra: int = 500) -> tuple[int | None, list[int], int]:
    committed = lcm_upto(x)
    speeds = BASE + [committed]
    for d in range(2, x + scan_extra + 1):
        ws = witnesses(speeds, d)
        if ws:
            return d, ws, committed % d
    return None, [], -1


def divisibility_wall_holds(x: int) -> bool:
    committed = lcm_upto(x)
    speeds = BASE + [committed]
    for d in range(2, x + 1):
        if committed % d != 0:
            return False
        if witnesses(speeds, d):
            return False
    return True


def tournament_analysis(rows: list[tuple[int, int | None, list[int], int]]) -> None:
    vertices = [
        "fixed_finite_basis",
        "nextprime_rule",
        "prime_power_wall",
        "radical_plus_residue_control",
        "torus_equidistribution_tail",
        "legendre_venn_finite_node",
        "three_gap_ap_hull_node",
    ]
    # Criteria in the order used for the pairwise observable.
    criteria = {
        "fixed_finite_basis": (0, 0, 0, 1, 0),
        "nextprime_rule": (1, 1, 0, 1, 0),
        "prime_power_wall": (1, 1, 1, 1, 1),
        "radical_plus_residue_control": (1, 1, 1, 0, 1),
        "torus_equidistribution_tail": (1, 1, 1, 0, 1),
        "legendre_venn_finite_node": (0, 0, 1, 1, 1),
        "three_gap_ap_hull_node": (0, 0, 1, 0, 1),
    }
    tie_path = {v: i for i, v in enumerate(vertices)}

    def score(v: str) -> int:
        return sum(criteria[v])

    edges: dict[str, list[str]] = {v: [] for v in vertices}
    indeg = {v: 0 for v in vertices}
    for i, u in enumerate(vertices):
        for v in vertices[i + 1 :]:
            su, sv = score(u), score(v)
            if su > sv or (su == sv and tie_path[u] < tie_path[v]):
                winner, loser = u, v
            else:
                winner, loser = v, u
            edges[winner].append(loser)
            indeg[loser] += 1

    score_hist = Counter(indeg.values())

    def has_triangle() -> bool:
        edge_set = {(u, v) for u, outs in edges.items() for v in outs}
        for a in vertices:
            for b in vertices:
                for c in vertices:
                    if a != b and b != c and a != c:
                        if (a, b) in edge_set and (b, c) in edge_set and (c, a) in edge_set:
                            return True
        return False

    # In a transitive tournament, Kahn's algorithm returns the Hamiltonian path.
    q = deque([v for v in vertices if indeg[v] == 0])
    topo = []
    indeg2 = dict(indeg)
    while q:
        v = q.popleft()
        topo.append(v)
        for w in edges[v]:
            indeg2[w] -= 1
            if indeg2[w] == 0:
                q.append(w)

    print("\nTournament Analysis")
    print("  vertices:", ", ".join(vertices))
    print("  observable: route score on (unbounded wall, lcm-compatible, node-split preserving, formalizable, data-fit)")
    print("  tie Hamiltonian path:", " -> ".join(vertices))
    print("  score histogram by indegree:", dict(sorted(score_hist.items())))
    print("  directed 3-cycle present:", has_triangle())
    print("  Hamiltonian/topological path:", " -> ".join(topo))
    print("  strongest route:", topo[0])


def main() -> None:
    print("S114 lcm committed-speed denominator wall")
    print("Base row B={1,...,11,13}; threshold is inclusive >= 1/14.\n")

    base_first = None
    for d in range(2, 120):
        ws = witnesses(BASE, d)
        if ws:
            base_first = (d, ws)
            break
    print(f"Base B first witness denominator: D={base_first[0]}, witnesses={base_first[1]}")
    print()

    sample_x = [14, 15, 16, 20, 24, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 150, 180, 200]
    rows = []
    for x in sample_x:
        wall = divisibility_wall_holds(x)
        d, ws, residue = first_witness_denominator(x)
        rows.append((x, d, ws, residue))
        np = next_prime(x)
        relation = "nextprime" if d == np else "not-nextprime"
        print(
            f"X={x:3d} wall_D<=X={wall} firstD={d:3d} "
            f"factor={factor(d or 1):>5s} nextprime={np:3d} {relation:13s} "
            f"LmodD={residue:3d} witnesses={ws[:6]}"
        )

    mismatches = []
    for x in range(14, 141):
        d, ws, residue = first_witness_denominator(x)
        if d != next_prime(x):
            mismatches.append((x, d, next_prime(x), factor(d or 1), ws[:4], residue))

    print("\nExact conclusions")
    print("  1. For every sampled X, no denominator D<=X can witness: lcm(2..X) is 0 mod D.")
    print("  2. Therefore any universal finite denominator basis is impossible.")
    print("  3. firstD=nextprime(X) is not a theorem; mismatches occur often.")
    print(f"  4. In X=14..140, mismatches with nextprime(X): {len(mismatches)} / {140-14+1}.")
    print("  5. The first opening can be composite; X=110,120 open at D=121=11^2.")

    print("\nFirst ten nextprime-rule mismatches")
    for x, d, np, fac, ws, residue in mismatches[:10]:
        print(f"  X={x:3d}: firstD={d:3d} ({fac}), nextprime={np:3d}, LmodD={residue:3d}, witnesses={ws}")

    print("\nProof interpretation")
    print("  The rigorous wall is a divisor-lattice statement, not a primality statement.")
    print("  Node 3 cannot be replaced by a finite residue certificate: the committed speed")
    print("  forces witness denominators beyond every fixed bound.  The useful replacement")
    print("  is an effective equidistribution/prime-power-residue lemma for denominators")
    print("  just beyond the lcm wall, composed with the finite Node 2 AP/three-gap cap.")

    tournament_analysis(rows)


if __name__ == "__main__":
    main()
