#!/usr/bin/env python3
"""
lrc_multitier_crt_sieve_s562.py

codex-2026-06-02 S562

Abstract/exact probe for recursive multi-tier CRT sieves.

The S561 two-tier CRT prototype filters k+1=2q by a cheap mod-q tier and sends
only a generated residual to the expensive tier.  This script looks one level
above that: when a residual gate packet is lifted in one local factor, does the
visible Archimedean gap merely trade places with endpoint/frontier debt?

Vertices are whole gate packets, not runners.  A packet is

    {breaker} union {scale*q : 1 <= q < n, q != skip}.

For the known n=14 and n=18 hard ladders, multiplying scale by 2 translates the
same odd-payload packet in the dyadic tier.  The exact audit also includes the
n=17 skip-8 packet and its dyadic translates as a rank-one control.
"""

from __future__ import annotations

from collections import Counter, deque
from dataclasses import dataclass
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from itertools import combinations, permutations
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
S356 = SourceFileLoader(
    "lonely_runner_residue_probe_s356",
    str(ROOT / "04-computation" / "lonely_runner_residue_probe_s356.py"),
).load_module()


@dataclass(frozen=True)
class Packet:
    n: int
    scale: int
    skip: int
    breaker: int
    speeds: tuple[int, ...]
    factor_n: tuple[tuple[int, int], ...]
    v_scale: tuple[tuple[int, int], ...]
    v_skip: tuple[tuple[int, int], ...]
    gap_ratio: Fraction
    boundary: int
    frontier_product: Fraction
    missing_zero_branches: tuple[int, ...]
    classification: str

    @property
    def label(self) -> str:
        return f"n{self.n}:scale{self.scale}:skip{self.skip}"


def factor(value: int) -> tuple[tuple[int, int], ...]:
    out: list[tuple[int, int]] = []
    p = 2
    x = value
    while p * p <= x:
        if x % p == 0:
            e = 0
            while x % p == 0:
                x //= p
                e += 1
            out.append((p, e))
        p += 1 if p == 2 else 2
    if x > 1:
        out.append((x, 1))
    return tuple(out)


def vp(value: int, p: int) -> int:
    e = 0
    x = value
    while x and x % p == 0:
        x //= p
        e += 1
    return e


def valuation_vector(value: int, primes: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    return tuple((p, vp(value, p)) for p in primes)


def fmt(value: Fraction | None) -> str:
    return S356.fmt_frac(value)


def gcd_all(values: tuple[int, ...]) -> int:
    g = 0
    for value in values:
        g = gcd(g, value)
    return g


def packet_speeds(n: int, scale: int, skip: int, breaker: int = 1) -> tuple[int, ...]:
    raw = (breaker,) + tuple(scale * q for q in range(1, n) if q != skip)
    normalized = S356.normalize_speed_set(list(raw))
    if gcd_all(normalized) != 1:
        raise AssertionError(f"non-primitive packet after normalization: n={n}, scale={scale}")
    return normalized


def missing_zero_branches(n: int, speeds: tuple[int, ...]) -> tuple[int, ...]:
    # THM-369/390 semantics: if no speed is divisible by q<n, then t=1/q is
    # an open sieve witness; q=n is the compactified wall node.
    return tuple(q for q in range(2, n + 1) if all(speed % q for speed in speeds))


def make_packet(n: int, scale: int, skip: int, breaker: int = 1) -> Packet:
    speeds = packet_speeds(n, scale, skip, breaker)
    report = S356.report(f"n{n}_scale{scale}_skip{skip}", list(speeds))
    gap_ratio = report.max_gap / report.threshold
    boundary = report.boundary_witness_count
    missing = missing_zero_branches(n, speeds)
    primes = tuple(sorted({p for p, _ in factor(n)} | {p for p, _ in factor(scale)} | {p for p, _ in factor(skip)}))
    classification = "positive_gap" if report.max_gap else ("boundary_only" if boundary else "open_cover_candidate")
    return Packet(
        n=n,
        scale=scale,
        skip=skip,
        breaker=breaker,
        speeds=speeds,
        factor_n=factor(n),
        v_scale=valuation_vector(scale, primes),
        v_skip=valuation_vector(skip, primes),
        gap_ratio=gap_ratio,
        boundary=boundary,
        frontier_product=gap_ratio * boundary,
        missing_zero_branches=missing,
        classification=classification,
    )


def lift_edges(packets: list[Packet]) -> list[tuple[Packet, Packet, int]]:
    edges: list[tuple[Packet, Packet, int]] = []
    for left in packets:
        for right in packets:
            if left.n != right.n or left.skip != right.skip or right.scale <= left.scale:
                continue
            if right.scale % left.scale:
                continue
            quotient = right.scale // left.scale
            factors = factor(quotient)
            if len(factors) == 1 and factors[0][1] == 1:
                edges.append((left, right, factors[0][0]))
    return edges


def tournament_fingerprint(rows: list[Packet]) -> dict[str, object]:
    # Vertices are recursive CRT packet states.  The switch ranks the most
    # counterexample-like row higher: smaller visible gap, then same/lower
    # product debt, then deeper scale translation.
    def key(row: Packet) -> tuple[Fraction, Fraction, int, int]:
        depth = sum(e for _, e in row.v_scale)
        return (-row.gap_ratio, -row.frontier_product, depth, -row.boundary)

    n = len(rows)
    adj = [[False] * n for _ in range(n)]
    for i, left in enumerate(rows):
        for j, right in enumerate(rows):
            if i != j:
                adj[i][j] = key(left) > key(right)

    scores = [sum(adj[i]) for i in range(n)]
    c3 = 0
    for i, j, k in combinations(range(n), 3):
        cyc = (
            adj[i][j] and adj[j][k] and adj[k][i]
        ) or (
            adj[i][k] and adj[k][j] and adj[j][i]
        )
        c3 += int(cyc)

    def reaches(start: int) -> set[int]:
        seen = {start}
        todo = deque([start])
        while todo:
            u = todo.popleft()
            for v in range(n):
                if adj[u][v] and v not in seen:
                    seen.add(v)
                    todo.append(v)
        return seen

    remaining = set(range(n))
    sccs: list[int] = []
    while remaining:
        u = next(iter(remaining))
        ru = reaches(u)
        comp = {v for v in remaining if v in ru and u in reaches(v)}
        sccs.append(len(comp))
        remaining -= comp

    hamiltonian_paths: int | str = 0
    if n <= 9:
        for perm in permutations(range(n)):
            if all(adj[perm[i]][perm[i + 1]] for i in range(n - 1)):
                hamiltonian_paths += 1
    else:
        hamiltonian_paths = "skipped(n>9)"

    return {
        "vertices": [row.label for row in rows],
        "score_hist": dict(sorted(Counter(scores).items())),
        "directed_3_cycles": c3,
        "sccs": sorted(sccs, reverse=True),
        "hamiltonian_paths": hamiltonian_paths,
    }


def main() -> None:
    packets = [
        make_packet(14, 7, 6),
        make_packet(14, 14, 6),
        make_packet(14, 28, 6),
        make_packet(17, 17, 8),
        make_packet(17, 34, 8),
        make_packet(17, 68, 8),
        make_packet(18, 9, 8),
        make_packet(18, 18, 8),
        make_packet(18, 36, 8),
    ]

    print("S562 recursive multi-tier CRT sieve packet audit")
    print("=" * 72)
    print("1. Packet states")
    print("  label                n-factors       v(scale)              v(skip)        gap/th    boundary  gap*boundary missing")
    for row in packets:
        missing = ",".join(str(q) for q in row.missing_zero_branches) or "-"
        print(
            f"  {row.label:<20} {str(dict(row.factor_n)):<15} "
            f"{str(dict(row.v_scale)):<21} {str(dict(row.v_skip)):<14} "
            f"{fmt(row.gap_ratio):>8} {row.boundary:>9} {fmt(row.frontier_product):>13} {missing}"
        )

    print()
    print("2. Recursive lift edges")
    for left, right, p in lift_edges(packets):
        gap_factor = right.gap_ratio / left.gap_ratio if left.gap_ratio else Fraction(0)
        boundary_factor = Fraction(right.boundary, left.boundary) if left.boundary else Fraction(0)
        product_ok = left.frontier_product == right.frontier_product
        print(
            f"  {left.label} --*{p}--> {right.label}: "
            f"gap_factor={fmt(gap_factor)} boundary_factor={fmt(boundary_factor)} "
            f"product_preserved={product_ok}"
        )

    print()
    print("3. Chain products")
    for n, skip in [(14, 6), (17, 8), (18, 8)]:
        chain = [row for row in packets if row.n == n and row.skip == skip]
        products = sorted({row.frontier_product for row in chain})
        vectors = [dict(row.v_scale) for row in chain]
        print(
            f"  n={n} skip={skip}: products={[fmt(x) for x in products]} "
            f"scale_vectors={vectors}"
        )

    print()
    print("4. Tournament Analysis")
    print("  vertices: recursive CRT packet states, not runners")
    print("  pairwise observable: (gap/th, gap*boundary, scale valuation depth, boundary)")
    print("  switch/gauge: smaller visible gap wins; ties prefer conserved lower debt and deeper translation")
    print(f"  fingerprints: {tournament_fingerprint(packets)}")

    print()
    print("5. Abstract synthesis")
    print("  A multi-tier CRT sieve should be recursive on product-tree nodes:")
    print("    node = (base denominator n, skipped quotient label, valuation vector of scale).")
    print("  A successful local tier either gives a THM-369/S561 witness or exports a residual packet.")
    print("  On these residual packets, a dyadic lift divides the visible gap by 2 and multiplies exposed boundary debt by 2.")
    print("  Thus gap*boundary is conserved while the packet translates in one local tree.")
    print("  THM-391 says a bare zero-branch star peels, so a proof must retain descendant endpoints, owner labels, or cross-prime coupling.")


if __name__ == "__main__":
    main()
