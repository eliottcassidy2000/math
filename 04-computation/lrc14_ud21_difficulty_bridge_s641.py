#!/usr/bin/env python3
"""S641: compare why LRC n=14 and unit-distance n=21 are difficult.

This is a bridge computation.  It does not prove LRC(14), nor an exact
unit-distance upper bound.  It records the shared carrier shape:

    27-quantum + retained section + bulk/lift side channel.

On the LRC side, 27 is the THM-401 shell modulus C=2n-1.  On the
unit-distance side, 27 is the edge increment per full Moser slab in THM-408.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import combinations
from math import gcd
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.append(str(SCRIPT_DIR))

import unit_distance_spine_ladder_s628 as s628  # noqa: E402


def factor(n: int) -> str:
    if n == 0:
        return "0"
    x = abs(n)
    pieces: list[str] = []
    p = 2
    while p * p <= x:
        if x % p == 0:
            e = 0
            while x % p == 0:
                x //= p
                e += 1
            pieces.append(str(p) if e == 1 else f"{p}^{e}")
        p += 1 if p == 2 else 2
    if x > 1:
        pieces.append(str(x))
    return " * ".join(pieces) if pieces else "1"


def residue_bits(residues: tuple[int, ...]) -> int:
    mask = 0
    for r in residues:
        mask |= 1 << r
    return mask


def lrc_full_mask(C: int) -> int:
    return residue_bits(tuple(range(1, C)))


def lrc_kill_mask(C: int, speed: int) -> int:
    mask = 0
    for j in range(1, C):
        if (speed * j) % C in (0, 1, C - 1):
            mask |= 1 << j
    return mask


def lrc_survivors(C: int, speeds: tuple[int, ...]) -> int:
    mask = lrc_full_mask(C)
    for speed in speeds:
        mask &= ~lrc_kill_mask(C, speed)
    return mask


def lrc_depths(C: int, speeds: tuple[int, ...]) -> dict[int, int]:
    depths = {j: 0 for j in range(1, C)}
    for speed in speeds:
        kill = lrc_kill_mask(C, speed)
        for j in range(1, C):
            if kill & (1 << j):
                depths[j] += 1
    return depths


def lrc_orbits(C: int) -> tuple[tuple[int, ...], ...]:
    unseen = set(range(1, C))
    orbits: list[tuple[int, ...]] = []
    while unseen:
        root = min(unseen)
        orbit: set[int] = set()
        stack = [root]
        while stack:
            x = stack.pop() % C
            if x == 0 or x in orbit:
                continue
            orbit.add(x)
            for y in ((2 * x) % C, (-x) % C):
                if y and y not in orbit:
                    stack.append(y)
        frozen = tuple(sorted(orbit))
        orbits.append(frozen)
        unseen -= orbit
    return tuple(orbits)


def ap_row(n: int) -> tuple[int, ...]:
    return tuple(range(1, n))


def vstar_row(n: int) -> tuple[int, ...]:
    return tuple(list(range(1, n - 2)) + [n - 1, 2 * n - 4])


@dataclass(frozen=True)
class LrcOrbitPacket:
    gcd_class: int
    orbit_size: int
    shell_reps: int
    released: int
    depth_hist: tuple[tuple[int, int], ...]
    redundancy_price: int


@dataclass(frozen=True)
class MoserPacket:
    name: str
    m: int
    vertices: int
    unit_edges: int
    spine_edges: int
    bulk_edges: int
    formula: str
    layer_quantum: int
    increment_from_previous: int | None


def lrc_packets(n: int, speeds: tuple[int, ...]) -> tuple[LrcOrbitPacket, ...]:
    C = 2 * n - 1
    survivors = lrc_survivors(C, speeds)
    depths = lrc_depths(C, speeds)
    packets: list[LrcOrbitPacket] = []
    for orbit in lrc_orbits(C):
        released = (residue_bits(orbit) & ~survivors).bit_count()
        hist = tuple(sorted(Counter(depths[j] for j in orbit).items()))
        price = sum(depth * count for depth, count in hist if depth > 0)
        packets.append(
            LrcOrbitPacket(
                gcd_class=gcd(orbit[0], C),
                orbit_size=len(orbit),
                shell_reps=len(orbit) // 2,
                released=released,
                depth_hist=hist,
                redundancy_price=price,
            )
        )
    return tuple(packets)


def moser_packet(name: str, m: int, builder, formula) -> MoserPacket:
    points = builder(m)
    vertices = len(points)
    edges = s628.unit_edges(tuple(sorted(points)))
    expected = formula(m)
    if edges != expected:
        raise AssertionError((name, m, edges, expected))
    spine_edges = vertices - 1
    prev_increment = None
    if m > 0:
        prev_edges = s628.unit_edges(tuple(sorted(builder(m - 1))))
        prev_increment = edges - prev_edges
    return MoserPacket(
        name=name,
        m=m,
        vertices=vertices,
        unit_edges=edges,
        spine_edges=spine_edges,
        bulk_edges=edges - spine_edges,
        formula=f"{expected}",
        layer_quantum=27,
        increment_from_previous=prev_increment,
    )


@dataclass(frozen=True)
class Route:
    name: str
    preserves_predicate: int
    cross_transfer: int
    computable_now: int
    proof_value: int
    scalar_risk: int


ROUTES = [
    Route("27_quantum_section_bulk_ledger", 5, 5, 5, 5, 1),
    Route("lift_ear_conservativity_transfer", 5, 5, 4, 5, 1),
    Route("side_channel_jackknife_damage_table", 5, 4, 5, 4, 1),
    Route("same_scalar_twin_search", 4, 4, 4, 4, 2),
    Route("centered_hex_vs_gcd9_microchannel", 4, 3, 4, 3, 2),
    Route("raw_14_21_number_match", 1, 1, 5, 1, 5),
]


def beats(a: Route, b: Route) -> bool:
    score_a = (
        3 * a.preserves_predicate
        + 2 * a.cross_transfer
        + 2 * a.proof_value
        + a.computable_now
        - 2 * a.scalar_risk
    )
    score_b = (
        3 * b.preserves_predicate
        + 2 * b.cross_transfer
        + 2 * b.proof_value
        + b.computable_now
        - 2 * b.scalar_risk
    )
    if score_a != score_b:
        return score_a > score_b
    return a.name < b.name


def hamiltonian_paths_tournament(routes: list[Route]) -> int:
    n = len(routes)
    dp: dict[tuple[int, int], int] = {}
    for i in range(n):
        dp[(1 << i, i)] = 1
    for mask in range(1 << n):
        for last in range(n):
            count = dp.get((mask, last), 0)
            if not count:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if beats(routes[last], routes[nxt]):
                    dp[(mask | (1 << nxt), nxt)] = dp.get((mask | (1 << nxt), nxt), 0) + count
    full = (1 << n) - 1
    return sum(dp.get((full, last), 0) for last in range(n))


def directed_triangles(routes: list[Route]) -> int:
    total = 0
    for a, b, c in combinations(range(len(routes)), 3):
        # A directed triangle occurs when no vertex beats both others.
        wins = []
        for x, y in ((a, b), (a, c), (b, c)):
            wins.append((x, y) if beats(routes[x], routes[y]) else (y, x))
        out = Counter(x for x, _ in wins)
        if sorted(out.values()) == [1, 1, 1]:
            total += 1
    return total


def print_lrc_section() -> None:
    n = 14
    C = 2 * n - 1
    print("LRC n=14 side")
    print("-------------")
    print(f"runner count minus observer: {n - 1}")
    print(f"C=2n-1={C}={factor(C)}")
    print("THM-407 action: residues folded by <2,-1>.")
    print("Canonical AP row packets:")
    for packet in lrc_packets(n, ap_row(n)):
        print(
            f"  gcd={packet.gcd_class:<2} orbit_size={packet.orbit_size:<2} "
            f"shell_reps={packet.shell_reps:<1} released={packet.released:<2} "
            f"depth_hist={dict(packet.depth_hist)} redundancy_price={packet.redundancy_price}"
        )
    print("Canonical V* row has the same 18/6/2 release profile in S624.")
    print("Open seam: arbitrary integer lifts need carry/CRT conservativity.")
    print("Carry identity: v=r+27k gives v == r-k (mod 14), since 27 == -1 (mod 14).")
    print()


def print_moser_section() -> None:
    print("Unit-distance Moser side")
    print("------------------------")
    rows = [
        moser_packet("P_1^+", 1, s628.path_plus, lambda m: 8 if m == 0 else 27 * m + 6),
        moser_packet("P_2^-", 2, s628.path_minus, lambda m: 6 if m == 0 else 27 * m + 3),
        moser_packet("P_2^+", 2, s628.path_plus, lambda m: 8 if m == 0 else 27 * m + 6),
    ]
    for packet in rows:
        print(
            f"{packet.name}: m={packet.m} vertices={packet.vertices} "
            f"unit_edges={packet.unit_edges} spine={packet.spine_edges} "
            f"bulk={packet.bulk_edges} formula={packet.formula} "
            f"increment_from_previous={packet.increment_from_previous}"
        )
    print("For m>=1, the closed formulas are:")
    print("  plus:  |P_m^+|=8m+6, E=27m+6, bulk=19m+1")
    print("  minus: |P_m^-|=8m+5, E=27m+3, bulk=19m-1")
    print("The exact n=21 row is P_2^-: E=2*27+3 and bulk=37=C_hex(3).")
    print("Open seam: endpoint-compatible ears and bulk/embedding labels for n=22.")
    print()


def print_bridge() -> None:
    print("Bridge thesis")
    print("-------------")
    print("The connection is not that 14 and 21 are numerically adjacent in some")
    print("mystical way.  It is that both frontiers are governed by a 27-quantum")
    print("carrier whose proof-relevant section is easy after it is named, while")
    print("the remaining difficulty is a small retained side channel.")
    print()
    print("Dictionary:")
    print("  LRC C=27 shell quotient  <->  Moser full-slab edge quantum 27")
    print("  LRC unit-shell section   <->  unit Hamiltonian spine section")
    print("  LRC nonunit/carry bulk   <->  unit-distance tile/bulk and ears")
    print("  lift/CRT conservativity  <->  endpoint-compatible ear conservativity")
    print("  gcd-9 tiny high-depth channel <-> direction/gain packets with small mass")
    print()
    print("Practical test:")
    print("  Build one channel-complete ledger that can delete a channel, measure")
    print("  false survivors or lost edges, and then ask whether the lost channel")
    print("  is exactly the lift/ear data needed for the proof.")
    print()


def print_external_cross_signal() -> None:
    print("External tournament-spectrum cross-signal")
    print("-----------------------------------------")
    print("The parallel monad S6 exhaustive n=9 H-spectrum gives:")
    print("  iso classes=191536, distinct H values=1520, H range=[1,3357]")
    print("  low odd gaps <=609 are still exactly [7, 21]")
    print("  the n=8 high gaps above 609 all fill in at n=9")
    print("This supports treating 21 as a durable obstruction value, while the")
    print("LRC/UD bridge itself still runs through the 27-section/bulk mechanism,")
    print("not through a literal unit-distance realization of H=21.")
    print()


def print_tournament() -> None:
    print("Tournament Analysis over bridge routes")
    print("--------------------------------------")
    routes = list(ROUTES)
    scores = {route.name: sum(1 for other in routes if route is not other and beats(route, other)) for route in routes}
    hist = Counter(scores.values())
    ranking = sorted(routes, key=lambda route: (-scores[route.name], route.name))
    print(f"vertices={len(routes)}")
    print(f"score_hist={dict(sorted(hist.items()))}")
    print(f"directed_3cycles={directed_triangles(routes)}")
    print(f"hamiltonian_paths={hamiltonian_paths_tournament(routes)}")
    print("tie Hamiltonian path:")
    for route in ranking:
        print(f"  score={scores[route.name]} {route.name}")
    print()


def main() -> None:
    print("S641 LRC n=14 / Unit-Distance n=21 Difficulty Bridge")
    print("====================================================")
    print()
    print_lrc_section()
    print_moser_section()
    print_bridge()
    print_external_cross_signal()
    print_tournament()
    print("Assumption challenge")
    print("--------------------")
    print("Vertices considered: LRC runners, LRC residues, gcd strata, carry")
    print("cocycles, Moser points, Moser directions, slabs, ears, exact cores,")
    print("and proof routes.  The script uses proof routes because the preserved")
    print("predicate is which transfer can explain both difficulties without")
    print("forgetting the side channel.")


if __name__ == "__main__":
    main()
