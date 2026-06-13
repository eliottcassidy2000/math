#!/usr/bin/env python3
"""
lrc_unit_distance_channel_ledger_s624.py

S624: bounce between LRC and the unit-distance frontier by damaging retained
channels one at a time.

The LRC side uses the S599w-x state-local witness mask over Z/(2n-1).  The
unit-distance side reuses the S623 Moser-carrier helpers but keeps the run
small.  The point is not to prove either frontier; it is to compare which
side channels become load-bearing when the raw frontier scalar is impaired.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import combinations
from math import gcd
import sys

import unit_distance_impairment_atlas_s623 as s623


def bit_count(mask: int) -> int:
    return mask.bit_count()


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


def lrc_survival_mask(C: int, speed: int) -> int:
    return lrc_full_mask(C) & ~lrc_kill_mask(C, speed)


def lrc_survivors(C: int, speeds: tuple[int, ...]) -> int:
    mask = lrc_full_mask(C)
    for speed in speeds:
        mask &= lrc_survival_mask(C, speed)
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
    # Matches the repo n=14 row (1,...,11,13,24): replace n-2 by 2n-4.
    return tuple(list(range(1, n - 2)) + [n - 1, 2 * n - 4])


def doubled_ap_row(n: int) -> tuple[int, ...]:
    return tuple(2 * i for i in range(1, n))


def fib_chain_row(n: int) -> tuple[int, ...]:
    values = [1, 3]
    while len(values) < n - 1:
        values.append(values[-1] + values[-2])
    return tuple(values)


@dataclass(frozen=True)
class LrcOrbitDamage:
    orbit: tuple[int, ...]
    orbit_gcd: int
    released_residues: int
    depth_hist: tuple[tuple[int, int], ...]
    min_depth: int
    max_depth: int


@dataclass(frozen=True)
class LrcRowDamage:
    n: int
    C: int
    row_name: str
    speeds: tuple[int, ...]
    survivor_count: int
    blocked_at_C_tick: bool
    orbit_damage: tuple[LrcOrbitDamage, ...]


def lrc_row_damage(n: int, row_name: str, speeds: tuple[int, ...]) -> LrcRowDamage:
    C = 2 * n - 1
    survivors = lrc_survivors(C, speeds)
    depths = lrc_depths(C, speeds)
    damages: list[LrcOrbitDamage] = []
    for orbit in lrc_orbits(C):
        omask = residue_bits(orbit)
        released = bit_count(omask & ~survivors)
        vals = [depths[j] for j in orbit]
        hist = tuple(sorted(Counter(vals).items()))
        damages.append(
            LrcOrbitDamage(
                orbit=orbit,
                orbit_gcd=gcd(orbit[0], C),
                released_residues=released,
                depth_hist=hist,
                min_depth=min(vals),
                max_depth=max(vals),
            )
        )
    return LrcRowDamage(
        n=n,
        C=C,
        row_name=row_name,
        speeds=speeds,
        survivor_count=bit_count(survivors),
        blocked_at_C_tick=survivors == 0,
        orbit_damage=tuple(damages),
    )


def print_lrc_jackknife() -> list[LrcRowDamage]:
    print("LRC WITNESS-ORBIT JACKKNIFE")
    print("---------------------------")
    print("Predicate: C-tick witness survival over C=2n-1.")
    print("Impairment: delete one <2,-1> witness orbit from the test.")
    print("Damage: released residues that would become false survivors in the impaired quotient.")
    print()

    jobs: list[tuple[int, str, tuple[int, ...]]] = []
    for n in (8, 11, 12, 14):
        jobs.append((n, "AP", ap_row(n)))
        jobs.append((n, "Vstar", vstar_row(n)))
    jobs.append((14, "doubled_AP", doubled_ap_row(14)))
    jobs.append((14, "fib_chain", fib_chain_row(14)))

    rows: list[LrcRowDamage] = []
    for n, name, speeds in jobs:
        row = lrc_row_damage(n, name, speeds)
        rows.append(row)
        status = "blocked" if row.blocked_at_C_tick else "survives"
        print(
            f"n={row.n:2d} C={row.C:2d} {row.row_name:10s} "
            f"survivors={row.survivor_count:2d} status={status} speeds={row.speeds}"
        )
        for damage in row.orbit_damage:
            depth_total = sum(depth * count for depth, count in damage.depth_hist)
            avg_depth = depth_total / len(damage.orbit)
            redundancy_price = damage.released_residues * avg_depth
            print(
                f"  orbit_gcd={damage.orbit_gcd:<2d} size={len(damage.orbit):2d} "
                f"released={damage.released_residues:2d} "
                f"depth=[{damage.min_depth},{damage.max_depth}] "
                f"hist={dict(damage.depth_hist)} "
                f"redundancy_price={redundancy_price:.1f}"
            )
        print()

    print("Reading:")
    print("- Composite C shells split into wide shallow channels and narrow redundant channels.")
    print("- At n=14, AP/Vstar/doubled_AP release 18 unit residues, 6 gcd-3 residues,")
    print("  and 2 gcd-9 residues when those channels are omitted; the gcd-9 residues")
    print("  are covered four times, so small mass can still be high-certainty structure.")
    print("- The prime C=23 control has one orbit: the raw bitmask has no shell-channel")
    print("  choice to damage, so side-channel scarcity appears only after lift/carry data.")
    print()
    return rows


@dataclass(frozen=True)
class UdDirectionLoss:
    rep: tuple[int, int, int, int]
    loss: int
    true_edges: int
    span: int


def print_unit_distance_fingerprint() -> tuple[list[UdDirectionLoss], dict[int, int]]:
    print("UNIT-DISTANCE LIGHTWEIGHT MOSER FINGERPRINT")
    print("-------------------------------------------")
    print("This recomputes a small target=9 Moser jackknife and cites S623 for the heavier atlas.")
    target = 9
    width = 30
    full = s623.beam_search(
        "moser",
        target,
        width,
        s623.MOSER_UNITS,
        s623.MOSER_UNITS,
        s623.add4,
        s623.canon4,
        s623.span4,
        "healthy",
    )
    losses: list[UdDirectionLoss] = []
    reps = s623.antipodal_reps(s623.MOSER_UNITS)
    full_profile = s623.shell_direction_profile(full.cluster, reps, s623.add4)
    for rep in reps:
        kept = s623.drop_direction(s623.MOSER_UNITS, rep)
        result = s623.beam_search(
            "moser",
            target,
            width,
            kept,
            s623.MOSER_UNITS,
            s623.add4,
            s623.canon4,
            s623.span4,
            "healthy",
        )
        losses.append(UdDirectionLoss(rep, full.true_edges - result.true_edges, result.true_edges, result.span))

    print(f"Moser target={target} width={width}: full true_edges={full.true_edges}, span={full.span}")
    hist = Counter(loss.loss for loss in losses)
    print(f"Direction-drop loss histogram: {dict(sorted(hist.items()))}")
    for loss, usage in zip(losses, full_profile):
        print(
            f"  drop {loss.rep}: usage_in_full={usage}, loss={loss.loss}, "
            f"usage_loss_price={usage * loss.loss}, true_edges={loss.true_edges}, span={loss.span}"
        )
    print()

    ceiling_results: dict[int, int] = {}
    print("Gain-ceiling comparison at the same target:")
    for ceiling in (1, 2, 3, 4):
        result = s623.beam_search(
            "moser",
            target,
            width,
            s623.MOSER_UNITS,
            s623.MOSER_UNITS,
            s623.add4,
            s623.canon4,
            s623.span4,
            "healthy",
            gain_ceiling=ceiling,
        )
        ceiling_results[ceiling] = result.true_edges
        print(
            f"  ceiling={ceiling}: true_edges={result.true_edges}, span={result.span}, "
            f"gain_hist={dict(result.gain_hist)}"
        )
    print()
    print("Stored S623 parent signal: target=10 loses one edge for six of nine")
    print("antipodal Moser directions; target=14 needs width 30 to reach 33 edges.")
    print()
    return losses, ceiling_results


@dataclass(frozen=True)
class Channel:
    name: str
    proof_relevance: int
    scaling_value: int
    small_damage: int
    predicate_fidelity: int
    side_info: int
    cost: int
    note: str


CHANNELS = (
    Channel("channel-complete ledger", 5, 4, 5, 5, 5, 3, "keeps both scalar frontier and side channels"),
    Channel("LRC witness-orbit jackknife", 4, 4, 5, 4, 4, 1, "tests which C=2n-1 residue channels are load-bearing"),
    Channel("UD direction-drop jackknife", 4, 4, 5, 4, 4, 2, "tests unit-shell indispensability by leave-one-direction-out"),
    Channel("UD gain-threshold solver", 5, 5, 4, 4, 4, 3, "asks for gain-4/5 extension packets before geometry checks"),
    Channel("LRC owner/carry/pinch fiber", 5, 3, 4, 5, 5, 4, "reattaches the rare labels needed by n=14 quotient proofs"),
    Channel("UD canonical orbit budget", 3, 5, 3, 3, 4, 2, "separates duplicate traffic from mathematical scarcity"),
    Channel("raw survival bitmask", 2, 5, 2, 2, 1, 1, "fast LRC scalar frontier but forgets lift/carry proof data"),
    Channel("raw wider Moser beam", 2, 4, 2, 2, 1, 5, "more compute without naming the retained invariant"),
    Channel("count-only quotient", 1, 3, 1, 1, 1, 1, "same scalar count can be non-equivalent across channels"),
)


def proof_gauge(channel: Channel) -> tuple[int, ...]:
    return (
        channel.proof_relevance,
        channel.predicate_fidelity,
        channel.side_info,
        channel.small_damage,
        channel.scaling_value,
        -channel.cost,
    )


def scaling_gauge(channel: Channel) -> tuple[int, ...]:
    return (
        channel.scaling_value,
        -channel.cost,
        channel.small_damage,
        channel.side_info,
        channel.predicate_fidelity,
        channel.proof_relevance,
    )


def make_tournament(gauge) -> list[list[int]]:
    n = len(CHANNELS)
    adj = [[0] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        ai = gauge(CHANNELS[i])
        aj = gauge(CHANNELS[j])
        if ai > aj or (ai == aj and i < j):
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return adj


def score_hist(adj: list[list[int]]) -> dict[int, int]:
    hist: dict[int, int] = {}
    for row in adj:
        score = sum(row)
        hist[score] = hist.get(score, 0) + 1
    return dict(sorted(hist.items()))


def directed_triangles(adj: list[list[int]]) -> int:
    count = 0
    for i, j, k in combinations(range(len(adj)), 3):
        if adj[i][j] and adj[j][k] and adj[k][i]:
            count += 1
        if adj[i][k] and adj[k][j] and adj[j][i]:
            count += 1
    return count


def hamiltonian_path_count(adj: list[list[int]]) -> int:
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for last in range(n):
            ways = dp[mask][last]
            if not ways:
                continue
            for nxt in range(n):
                if not (mask >> nxt) & 1 and adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += ways
    return sum(dp[-1])


def hamiltonian_path(adj: list[list[int]]) -> list[int]:
    path: list[int] = []
    for vertex in range(len(adj)):
        pos = 0
        while pos < len(path) and adj[path[pos]][vertex]:
            pos += 1
        path.insert(pos, vertex)
    return path


def scc_sizes(adj: list[list[int]]) -> list[int]:
    n = len(adj)

    def reach_from(start: int, transpose: bool) -> set[int]:
        seen = {start}
        stack = [start]
        while stack:
            v = stack.pop()
            for w in range(n):
                edge = adj[w][v] if transpose else adj[v][w]
                if edge and w not in seen:
                    seen.add(w)
                    stack.append(w)
        return seen

    remaining = set(range(n))
    sizes = []
    while remaining:
        v = min(remaining)
        comp = reach_from(v, False) & reach_from(v, True)
        sizes.append(len(comp))
        remaining -= comp
    return sorted(sizes, reverse=True)


def edge_flips(a: list[list[int]], b: list[list[int]]) -> list[tuple[str, str]]:
    flips: list[tuple[str, str]] = []
    for i, j in combinations(range(len(a)), 2):
        if a[i][j] != b[i][j]:
            flips.append((CHANNELS[i].name, CHANNELS[j].name))
    return flips


def print_tournament_analysis() -> None:
    print("TOURNAMENT ANALYSIS: CHANNELS AND PROOF OBLIGATIONS")
    print("---------------------------------------------------")
    print("Vertices: channels/proof obligations, not LRC runners or unit-distance points.")
    print("Pairwise observable: which channel better preserves frontier decisions under impairment.")
    print("Switches: proof-relevance gauge versus scaling-value gauge.")
    print()

    proof_adj = make_tournament(proof_gauge)
    scaling_adj = make_tournament(scaling_gauge)
    for label, adj in (("proof gauge", proof_adj), ("scaling gauge", scaling_adj)):
        scores = [sum(row) for row in adj]
        path = hamiltonian_path(adj)
        print(label)
        print(f"  score_hist={score_hist(adj)}")
        print(f"  directed_3cycles={directed_triangles(adj)}")
        print(f"  scc_sizes={scc_sizes(adj)}")
        print(f"  hamiltonian_path_count={hamiltonian_path_count(adj)}")
        print("  tie Hamiltonian path:")
        for idx in path:
            print(f"    score={scores[idx]} {CHANNELS[idx].name}: {CHANNELS[idx].note}")
        print()

    flips = edge_flips(proof_adj, scaling_adj)
    print(f"Edge flips between gauges: {len(flips)} / {len(CHANNELS) * (len(CHANNELS) - 1) // 2}")
    for left, right in flips:
        print(f"  flip: {left} <> {right}")
    print()


def print_assumption_challenge() -> None:
    print("ASSUMPTION CHALLENGE")
    print("--------------------")
    print("Alternate vertices considered: LRC runners, LRC witness residues, shell/gcd")
    print("strata, owner routes, carry cocycles, unit-distance points, unit directions,")
    print("gain packets, deletion cores, canonical orbit classes, obstruction filters,")
    print("and proof obligations.")
    print("Chosen vertices: channels/proof obligations, because the preserved predicate")
    print("is whether a damaged method changes the frontier decision.")
    print("Information preserved: C-tick survival response, witness-orbit release mass,")
    print("Moser direction loss, gain-ceiling loss, and channel ranking under two gauges.")
    print("Information destroyed: exact continuous LRC time geometry, arbitrary integer")
    print("lift provenance, full unit-distance embedding outside the finite carrier, and")
    print("proof-grade totally-unfaithful obstruction certificates.")
    print("Challenged assumption: widening a frontier-gain beam is the next unit of")
    print("progress.  The bridge table says widen only after naming the side channel that")
    print("the width is preserving.")
    print()


def main() -> None:
    sys.stdout.reconfigure(line_buffering=True)
    print("S624 LRC / UNIT-DISTANCE CHANNEL LEDGER")
    print("=======================================")
    print()
    print_lrc_jackknife()
    print_unit_distance_fingerprint()
    print_tournament_analysis()
    print_assumption_challenge()


if __name__ == "__main__":
    main()
