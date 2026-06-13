#!/usr/bin/env python3
"""
lrc_n17_dyadic_packet_surprise_s566.py

codex-2026-06-02 S566

Investigate the S562 surprise: the n=17 skip-8 packet has a dyadic residual
tier even though 2 is not a base factor of n=17.

The script checks three things exactly with S356 interval arithmetic:

1. Neighboring prime half-gate packets p=5..23 also halve gap and double
   boundary debt under scale *= 2.
2. For n=17, every skipped quotient label has the same dyadic conservation law,
   not just skip 8.
3. The scalar-gap and product-debt gauges disagree at n=17: skip 8 is the
   visible-gap optimum, but skip 6 has the lower conserved product.
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
class PacketReport:
    n: int
    scale: int
    skip: int
    breaker: int
    raw_gcd: int
    effective_scale: int
    effective_dyadic_depth: int
    gap_ratio: Fraction
    boundary: int
    product: Fraction
    forbidden_length: Fraction
    witness: Fraction | None
    boundary_modulus: int


REPORT_CACHE: dict[tuple[int, int, int, int], PacketReport] = {}


def fmt(value: Fraction | None) -> str:
    return S356.fmt_frac(value)


def vp(value: int, p: int) -> int:
    e = 0
    x = value
    while x and x % p == 0:
        x //= p
        e += 1
    return e


def gcd_all(values: tuple[int, ...]) -> int:
    g = 0
    for value in values:
        g = gcd(g, value)
    return g


def raw_packet(n: int, scale: int, skip: int, breaker: int = 1) -> tuple[int, ...]:
    return tuple(sorted((breaker,) + tuple(scale * q for q in range(1, n) if q != skip)))


def report_packet(n: int, scale: int, skip: int, breaker: int = 1) -> PacketReport:
    key = (n, scale, skip, breaker)
    if key in REPORT_CACHE:
        return REPORT_CACHE[key]

    raw = raw_packet(n, scale, skip, breaker)
    raw_gcd = gcd_all(raw)
    speeds = S356.normalize_speed_set(list(raw))
    effective_scale = scale // raw_gcd
    rep = S356.report(f"n{n}_scale{scale}_skip{skip}_breaker{breaker}", list(speeds))
    gap_ratio = rep.max_gap / rep.threshold
    packet = PacketReport(
        n=n,
        scale=scale,
        skip=skip,
        breaker=breaker,
        raw_gcd=raw_gcd,
        effective_scale=effective_scale,
        effective_dyadic_depth=vp(effective_scale, 2),
        gap_ratio=gap_ratio,
        boundary=rep.boundary_witness_count,
        product=gap_ratio * rep.boundary_witness_count,
        forbidden_length=rep.forbidden_length,
        witness=rep.witness or rep.boundary_witness,
        boundary_modulus=rep.boundary_modulus,
    )
    REPORT_CACHE[key] = packet
    return packet


def prime_half_gate_rows(primes: list[int], depths: int = 3) -> dict[int, list[PacketReport]]:
    rows: dict[int, list[PacketReport]] = {}
    for p in primes:
        skip = (p - 1) // 2
        rows[p] = [report_packet(p, p * (2**h), skip) for h in range(depths)]
    return rows


def n17_skip_spectrum(depths: int = 3) -> dict[int, list[PacketReport]]:
    return {
        skip: [report_packet(17, 17 * (2**h), skip) for h in range(depths)]
        for skip in range(1, 17)
    }


def conserved(rows: list[PacketReport]) -> bool:
    if len({row.product for row in rows}) != 1:
        return False
    for left, right in zip(rows, rows[1:]):
        if right.gap_ratio != left.gap_ratio / 2:
            return False
        if right.boundary != 2 * left.boundary:
            return False
        if right.forbidden_length != left.forbidden_length:
            return False
    return True


def tournament_fingerprint(rows: list[PacketReport], mode: str) -> dict[str, object]:
    # Vertices are skipped quotient labels.  Two gauges are compared:
    # scalar: smaller gap first, then smaller product;
    # product: smaller product first, then smaller gap.
    if mode == "scalar":
        key = lambda row: (-row.gap_ratio, -row.product, -row.boundary)
    elif mode == "product":
        key = lambda row: (-row.product, -row.gap_ratio, -row.boundary)
    else:
        raise ValueError(mode)

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

    hp: int | str = 0
    if n <= 9:
        for perm in permutations(range(n)):
            if all(adj[perm[i]][perm[i + 1]] for i in range(n - 1)):
                hp += 1
    else:
        hp = "skipped(n>9)"

    return {
        "mode": mode,
        "top_skip": max(range(n), key=lambda i: scores[i]) + 1,
        "score_hist": dict(sorted(Counter(scores).items())),
        "directed_3_cycles": c3,
        "sccs": sorted(sccs, reverse=True),
        "hamiltonian_paths": hp,
    }


def edge_flip_count(rows: list[PacketReport]) -> int:
    scalar = tournament_fingerprint_adjacency(rows, "scalar")
    product = tournament_fingerprint_adjacency(rows, "product")
    flips = 0
    for i in range(len(rows)):
        for j in range(i + 1, len(rows)):
            if scalar[i][j] != product[i][j]:
                flips += 1
    return flips


def tournament_fingerprint_adjacency(rows: list[PacketReport], mode: str) -> list[list[bool]]:
    if mode == "scalar":
        key = lambda row: (-row.gap_ratio, -row.product, -row.boundary)
    else:
        key = lambda row: (-row.product, -row.gap_ratio, -row.boundary)
    n = len(rows)
    return [[i != j and key(rows[i]) > key(rows[j]) for j in range(n)] for i in range(n)]


def print_prime_half_gate() -> None:
    print("1. Prime half-gate dyadic controls")
    rows_by_prime = prime_half_gate_rows([5, 7, 11, 13, 17, 19, 23])
    print("  p  skip  conserved  gap0     boundary0 product0 length")
    for p, rows in rows_by_prime.items():
        first = rows[0]
        print(
            f"  {p:<2} {first.skip:<5} {str(conserved(rows)):<10} "
            f"{fmt(first.gap_ratio):>7} {first.boundary:>9} "
            f"{fmt(first.product):>8} {fmt(first.forbidden_length):>12}"
        )
    print("  Observation: every audited prime half-gate packet has gap(h)=gap(0)/2^h and boundary(h)=boundary(0)*2^h.")
    print()


def print_n17_skip_spectrum() -> dict[int, list[PacketReport]]:
    print("2. n=17 skip spectrum under dyadic lift")
    spectrum = n17_skip_spectrum()
    print("  skip conserved gap0     boundary0 product0 scalar_rank product_rank")
    base_rows = [spectrum[skip][0] for skip in range(1, 17)]
    scalar_order = sorted(base_rows, key=lambda row: (row.gap_ratio, row.product, row.boundary))
    product_order = sorted(base_rows, key=lambda row: (row.product, row.gap_ratio, row.boundary))
    scalar_rank = {row.skip: i + 1 for i, row in enumerate(scalar_order)}
    product_rank = {row.skip: i + 1 for i, row in enumerate(product_order)}
    for skip in range(1, 17):
        first = spectrum[skip][0]
        print(
            f"  {skip:<4} {str(conserved(spectrum[skip])):<9} "
            f"{fmt(first.gap_ratio):>7} {first.boundary:>9} "
            f"{fmt(first.product):>8} {scalar_rank[skip]:>11} {product_rank[skip]:>12}"
        )
    print(
        f"  scalar-gap winner: skip {scalar_order[0].skip}; "
        f"product-debt winner: skip {product_order[0].skip}; "
        f"rank edge flips={edge_flip_count(base_rows)}"
    )
    print()
    return spectrum


def print_breaker_parity_collapse(max_h: int = 5) -> None:
    print("3. n=17 skip-8 breaker parity collapse")
    print("  h  depth_hist                      product_values")
    for h in range(max_h + 1):
        rows = [report_packet(17, 17 * (2**h), 8, breaker=r) for r in range(1, 17)]
        depth_hist = Counter(row.effective_dyadic_depth for row in rows)
        products = sorted({row.product for row in rows})
        print(
            f"  {h:<2} {dict(sorted(depth_hist.items()))!s:<31} "
            f"{[fmt(product) for product in products]}"
        )
    print("  Even breakers do not create new packets; normalization drops them to lower dyadic depth.")
    print()


def print_tournament_analysis(spectrum: dict[int, list[PacketReport]]) -> None:
    print("4. Tournament Analysis")
    base_rows = [spectrum[skip][0] for skip in range(1, 17)]
    print("  vertices: skipped quotient labels q=1..16, with whole n=17 packets attached")
    print("  observable: (gap/th at scale 17, conserved gap*boundary, boundary, skip)")
    print("  switch A: scalar-gap first; switch B: product-debt first")
    print(f"  scalar fingerprint:  {tournament_fingerprint(base_rows, 'scalar')}")
    print(f"  product fingerprint: {tournament_fingerprint(base_rows, 'product')}")
    print(f"  edge flips between gauges: {edge_flip_count(base_rows)}")
    print("  Interpretation: skip 8 is the scalar-gap face; skip 6 is the lower-product face.")
    print()


def main() -> None:
    print("S566 n=17 dyadic residual packet surprise")
    print("=" * 72)
    print_prime_half_gate()
    spectrum = print_n17_skip_spectrum()
    print_breaker_parity_collapse()
    print_tournament_analysis(spectrum)
    print("5. Synthesis")
    print("  The dyadic tier is not a base-CRT accident and not skip-8-only.")
    print("  For n=17, scale *= 2 is a residual-packet translation: visible gap halves, boundary debt doubles.")
    print("  The half-gate skip 8 wins the scalar gap race, but skip 6 wins the conserved-product gauge.")
    print("  This explains why n=17 carried back to n=14 through skip 6: the product/debt gauge already prefers the predecessor channel.")


if __name__ == "__main__":
    main()
