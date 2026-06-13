#!/usr/bin/env python3
"""Faulhaber odd-moment anchors and the OCF odd-cycle analogy.

This script is intentionally small and reproducible.  It checks the user's
power-balance identity, tests the proposed asymptotic anchor, and records a
Tournament Analysis carrier whose vertices are proof interfaces rather than
runners or arcs.
"""

from __future__ import annotations

from collections import Counter
from decimal import Decimal, getcontext
from functools import lru_cache
from itertools import combinations, permutations
from math import comb


def s_power(n: int, r: int) -> int:
    return sum(j**r for j in range(1, n + 1))


def balance_direct(c: int, p: int, n: int) -> int:
    left = c**p + sum((c - j) ** p for j in range(1, n + 1))
    right = sum((c + j) ** p for j in range(1, n + 1))
    return left - right


def balance_odd_formula(c: int, p: int, n: int) -> int:
    return c**p - 2 * sum(
        comb(p, r) * c ** (p - r) * s_power(n, r)
        for r in range(1, p + 1, 2)
    )


def decimal_balance(c: Decimal, p: int, n: int) -> Decimal:
    total = c**p
    for r in range(1, p + 1, 2):
        total -= Decimal(2 * comb(p, r) * s_power(n, r)) * (c ** (p - r))
    return total


def asymptotic_a(p: int, n: int) -> Decimal:
    pD = Decimal(p)
    nD = Decimal(n)
    u = nD * (nD + 1)
    a0 = pD * nD * nD + Decimal(p - 1) * nD
    a1 = Decimal((p - 1) * (p - 2)) / (Decimal(12) * pD)
    b_num = Decimal(-((p - 1) * (p - 2) * (2 * p * p - 4 * p - 1)))
    b_den = Decimal(180) * (pD**3)
    return a0 + a1 + b_num / b_den / u


def root_a(p: int, n: int) -> Decimal:
    """Bisection for the real anchor a_p(n), bracketed in base+[0,1]."""
    pD = Decimal(p)
    nD = Decimal(n)
    base = pD * nD * nD + Decimal(p - 1) * nD
    lo = base
    hi = base + Decimal(1)
    # Balance is written in c=a+n.
    flo = decimal_balance(lo + nD, p, n)
    fhi = decimal_balance(hi + nD, p, n)
    if flo == 0:
        return lo
    if fhi == 0:
        return hi
    if flo * fhi > 0:
        raise ValueError(f"root not bracketed for p={p}, n={n}: {flo}, {fhi}")
    for _ in range(160):
        mid = (lo + hi) / 2
        fm = decimal_balance(mid + nD, p, n)
        if fm == 0:
            return mid
        if flo * fm <= 0:
            hi = mid
            fhi = fm
        else:
            lo = mid
            flo = fm
    return (lo + hi) / 2


def fmt_decimal(x: Decimal, places: int = 12) -> str:
    q = Decimal(10) ** -places
    return str(x.quantize(q))


def verify_odd_formula() -> tuple[int, int]:
    checks = 0
    failures = 0
    for p in range(1, 10):
        for n in range(1, 9):
            for c in range(n + 1, n + 12):
                checks += 1
                if balance_direct(c, p, n) != balance_odd_formula(c, p, n):
                    failures += 1
    return checks, failures


def count_hamiltonian_paths(adj: list[list[bool]]) -> int:
    n = len(adj)

    @lru_cache(None)
    def dp(last: int, mask: int) -> int:
        if mask == (1 << last):
            return 1
        prev_mask = mask ^ (1 << last)
        total = 0
        for prev in range(n):
            if (prev_mask >> prev) & 1 and adj[prev][last]:
                total += dp(prev, prev_mask)
        return total

    full = (1 << n) - 1
    return sum(dp(last, full) for last in range(n))


def canonical_cycle(cyc: tuple[int, ...]) -> tuple[int, ...]:
    variants = []
    m = len(cyc)
    for i in range(m):
        rot = cyc[i:] + cyc[:i]
        variants.append(rot)
    return min(variants)


def directed_odd_cycles(adj: list[list[bool]]) -> list[tuple[int, ...]]:
    n = len(adj)
    cycles: set[tuple[int, ...]] = set()
    for k in range(3, n + 1, 2):
        for subset in combinations(range(n), k):
            start = min(subset)
            rest = [v for v in subset if v != start]
            for perm in permutations(rest):
                cyc = (start,) + perm
                if all(adj[cyc[i]][cyc[(i + 1) % k]] for i in range(k)):
                    cycles.add(canonical_cycle(cyc))
    return sorted(cycles)


def ocf_alpha_illustration() -> tuple[int, int, int, int]:
    """Two disjoint cyclic triangles with all cross edges one way."""
    n = 6
    adj = [[False] * n for _ in range(n)]
    tri1 = [(0, 1), (1, 2), (2, 0)]
    tri2 = [(3, 4), (4, 5), (5, 3)]
    for u, v in tri1 + tri2:
        adj[u][v] = True
        adj[v][u] = False
    for u in range(3):
        for v in range(3, 6):
            adj[u][v] = True
            adj[v][u] = False
    cycles = directed_odd_cycles(adj)
    disjoint_pairs = 0
    for a, b in combinations(cycles, 2):
        if set(a).isdisjoint(b):
            disjoint_pairs += 1
    ocf_value = 1 + 2 * len(cycles) + 4 * disjoint_pairs
    return len(cycles), disjoint_pairs, count_hamiltonian_paths(adj), ocf_value


def tournament_fingerprints(vertices: list[str], edges: dict[str, set[str]]) -> dict[str, object]:
    score_hist = Counter(len(edges[v]) for v in vertices)
    c3 = 0
    for tri in combinations(vertices, 3):
        outdegrees = [sum(1 for u in tri if u in edges[v]) for v in tri]
        if sorted(outdegrees) == [1, 1, 1]:
            c3 += 1

    def reach(start: str) -> set[str]:
        seen = {start}
        stack = [start]
        while stack:
            v = stack.pop()
            for u in edges[v]:
                if u not in seen:
                    seen.add(u)
                    stack.append(u)
        return seen

    unseen = set(vertices)
    scc_sizes = []
    while unseen:
        v = next(iter(unseen))
        rv = reach(v)
        comp = {u for u in unseen if v in reach(u) and u in rv}
        scc_sizes.append(len(comp))
        unseen -= comp
    scc_sizes.sort(reverse=True)

    n = len(vertices)

    @lru_cache(None)
    def dp(last: int, mask: int) -> int:
        if mask == (1 << last):
            return 1
        total = 0
        prev_mask = mask ^ (1 << last)
        for prev in range(n):
            if (prev_mask >> prev) & 1 and vertices[last] in edges[vertices[prev]]:
                total += dp(prev, prev_mask)
        return total

    full = (1 << n) - 1
    ham_paths = sum(dp(i, full) for i in range(n))
    return {
        "score_hist": dict(sorted(score_hist.items())),
        "directed_3cycles": c3,
        "scc_sizes": scc_sizes,
        "hamiltonian_paths": ham_paths,
        "scores": {v: len(edges[v]) for v in vertices},
    }


def carrier_tournament() -> tuple[list[str], dict[str, set[str]], dict[str, object]]:
    vertices = [
        "faulhaber_odd_moments",
        "ocf_odd_cycles",
        "omega_alpha_packets",
        "beatty_pell_addresses",
        "convolution_hidden_lifts",
        "lrc14_q27_support",
        "code72_support_design",
        "unit_distance_nonproduct_fibers",
    ]
    # Criteria order:
    # odd atom, compatibility layer, address/carry layer, support transfer,
    # exactness, computation-readiness, hidden-lift retention.
    profiles = {
        "faulhaber_odd_moments": [3, 1, 2, 1, 3, 3, 2],
        "ocf_odd_cycles": [3, 3, 1, 1, 3, 2, 3],
        "omega_alpha_packets": [2, 3, 1, 2, 3, 3, 2],
        "beatty_pell_addresses": [1, 1, 3, 3, 3, 3, 2],
        "convolution_hidden_lifts": [0, 2, 2, 3, 3, 3, 3],
        "lrc14_q27_support": [1, 3, 3, 3, 2, 2, 3],
        "code72_support_design": [1, 3, 2, 3, 2, 1, 2],
        "unit_distance_nonproduct_fibers": [0, 2, 3, 2, 2, 2, 2],
    }
    tie_path = {v: i for i, v in enumerate(vertices)}
    edges = {v: set() for v in vertices}
    for a, b in combinations(vertices, 2):
        pa = profiles[a]
        pb = profiles[b]
        aw = sum(x > y for x, y in zip(pa, pb))
        bw = sum(y > x for x, y in zip(pa, pb))
        if aw > bw or (aw == bw and tie_path[a] < tie_path[b]):
            edges[a].add(b)
        else:
            edges[b].add(a)
    return vertices, edges, tournament_fingerprints(vertices, edges)


def main() -> None:
    getcontext().prec = 90

    print("FAULHABER ODD-MOMENT / OCF BRIDGE")
    print("=" * 44)
    print()

    checks, failures = verify_odd_formula()
    print("1. Odd-moment balance identity")
    print(f"   exact checks: {checks}")
    print(f"   failures: {failures}")
    print(
        "   identity: D_p(c,n)=c^p-2*sum_{r odd<=p} binom(p,r)c^(p-r)S_r(n)"
    )
    print("   reason: paired terms (c-j)^p-(c+j)^p cancel all even r")
    print()

    print("2. Exact integer anchors")
    for p in (1, 2):
        rows = []
        for n in range(1, 8):
            c = p * n * (n + 1)
            a = c - n
            rows.append(f"n={n}:a={a},c={c},D={balance_odd_formula(c,p,n)}")
        print(f"   p={p}: " + "; ".join(rows))
    print(
        "   p=2 geometry side note: 6*S_2(n)=n(n+1)(2n+1), the cuboid"
    )
    print(
        "   packing volume.  The balance equation itself uses only odd S_1(n),"
    )
    print("   while the square-pyramid volume is the even shell package around it.")
    print()

    print("3. Asymptotic anchor check")
    print(
        "   a_asym=p*n^2+(p-1)*n+(p-1)(p-2)/(12p)"
        "-(p-1)(p-2)(2p^2-4p-1)/(180*p^3*n(n+1))"
    )
    print("   rows show true-a minus asymptotic-a; n^4 scaling should stabilize")
    for p in range(3, 9):
        cells = []
        for n in (10, 30, 100):
            ar = root_a(p, n)
            aa = asymptotic_a(p, n)
            err = ar - aa
            scaled = err * (Decimal(n) ** 4)
            cells.append(
                f"n={n}:err={fmt_decimal(err, 14)},n4err={fmt_decimal(scaled, 8)}"
            )
        print(f"   p={p}: " + "; ".join(cells))
    print()

    cycles, pairs, hpaths, ocf_value = ocf_alpha_illustration()
    print("4. OCF compatibility illustration")
    print(
        "   two disjoint directed triangles with all cross edges one way:"
        f" odd_cycles={cycles}, disjoint_pairs={pairs}, H={hpaths}, OCF={ocf_value}"
    )
    print(
        "   reading: Faulhaber balance is linear in odd moments S_r; OCF adds"
    )
    print(
        "   compatibility packets alpha_k, i.e. independent collections in Omega(T)."
    )
    print()

    vertices, _edges, fp = carrier_tournament()
    print("5. Tournament Analysis carrier")
    print("   vertices are proof carriers, not runners/arcs:")
    print("   " + ", ".join(vertices))
    print("   pairwise observable: majority over 7 retained-channel criteria")
    print("   switch/gauge: higher retained-channel score wins; ties use listed path")
    print("   preserves for LRC: blocker/support address, carry class, witness target")
    print("   destroys for LRC: raw time geometry and exact interval endpoints")
    print(
        "   challenged assumption: the useful tournament may live on proof obligations,"
    )
    print("   odd atoms, or hidden lifts rather than on runner pairs.")
    print(f"   score_hist={fp['score_hist']}")
    print(f"   directed_3cycles={fp['directed_3cycles']}")
    print(f"   scc_sizes={fp['scc_sizes']}")
    print(f"   hamiltonian_paths={fp['hamiltonian_paths']}")
    print(f"   scores={fp['scores']}")
    print()

    print("6. Working hypotheses generated")
    print("   H1: c=p*u+corrections with u=n(n+1) is the common triangular carrier;")
    print("       p=1,2 are integral because the Bernoulli correction vanishes.")
    print("   H2: odd Faulhaber moments are the scalar shadow of odd-cycle atoms;")
    print("       OCF suggests the missing higher structure is compatibility of atoms.")
    print("   H3: LRC14 Q27 and code72 support ledgers should retain odd atoms plus")
    print("       compatibility packets, mirroring alpha_k rather than only S_r.")
    print("   H4: HYP-2456's Beatty-Pell address is the endpoint/carry layer of the")
    print("       same carrier; the Faulhaber expansion is its moment layer.")


if __name__ == "__main__":
    main()
