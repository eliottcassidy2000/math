#!/usr/bin/env python3
"""
S115 depth guard for the corrected-Venn Bonferroni-3 target.

The KPS S31t target says the far-runner Newton expansion should satisfy a
Bonferroni-3 upper bound

    p0(B union F) <= T_1 + T_2 + T_3.

This script isolates the exact pointwise object.  At slow time x, let M be the
set of inner sectors not already hit by the base B, with d=|M|.  If far runner
i hits color c_i, then the pointwise cover predicate is

    prod_{a in M} (OR_{i:c_i=a} z_i).

Therefore the r-th Newton packet is

    T_r(x) = (-1)^(r+d) * C_{d,r}(x),

where C_{d,r}(x) counts r-subsets of far runners whose colors all lie in M and
whose colors cover M.  Thus the high tail is controlled by a depth-parity
ledger.  Venn containment alone does not force the r>=4 tail to be nonpositive.
"""

from __future__ import annotations

from collections import Counter, deque
from fractions import Fraction
from itertools import combinations
from math import gcd
from random import Random


INNER = set(range(1, 7))


def sector(speed: int, x: Fraction) -> int:
    y = (speed * x) % 1
    return (7 * y.numerator) // y.denominator


def breakpoints(speeds: list[int]) -> list[Fraction]:
    out = {Fraction(0), Fraction(1)}
    for speed in speeds:
        if speed == 0:
            continue
        for j in range(0, 7 * speed + 1):
            out.add(Fraction(j, 7 * speed))
    return sorted(out)


def p0_exact(speeds: list[int]) -> Fraction:
    speeds = sorted(set(s for s in speeds if s != 0))
    if not speeds:
        return Fraction(0)
    total = Fraction(0)
    bps = breakpoints(speeds)
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        hit = {sector(speed, mid) for speed in speeds} & INNER
        if hit == INNER:
            total += hi - lo
    return total


def cover_counts_by_r(colors: list[int], missing: set[int], far_count: int) -> list[int]:
    """Counts r-subsets using only missing colors and covering all missing colors."""
    missing_list = sorted(missing)
    pos = {color: i for i, color in enumerate(missing_list)}
    full = (1 << len(missing_list)) - 1
    dp = [[0] * (far_count + 1) for _ in range(1 << len(missing_list))]
    dp[0][0] = 1
    for color in colors:
        if color not in pos:
            continue
        bit = 1 << pos[color]
        for mask in range((1 << len(missing_list)) - 1, -1, -1):
            for r in range(far_count - 1, -1, -1):
                if dp[mask][r]:
                    dp[mask | bit][r + 1] += dp[mask][r]
    return dp[full]


def packet_decomposition(base: list[int], far: list[int]):
    speeds = sorted(set(s for s in base + far if s != 0))
    far_count = len(far)
    bps = breakpoints(speeds)
    t = [Fraction(0) for _ in range(far_count + 1)]
    depth_packet = {
        d: [Fraction(0) for _ in range(far_count + 1)] for d in range(7)
    }
    tail_by_depth = {d: Fraction(0) for d in range(7)}
    p0_base = Fraction(0)
    p0_full = Fraction(0)

    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        length = hi - lo
        base_hit = {sector(speed, mid) for speed in base if speed != 0} & INNER
        missing = INNER - base_hit

        if not missing:
            p0_base += length
            p0_full += length
            continue

        far_colors = [sector(speed, mid) for speed in far]
        if base_hit | ({c for c in far_colors} & INNER) == INNER:
            p0_full += length

        d = len(missing)
        counts = cover_counts_by_r(far_colors, missing, far_count)
        for r, count in enumerate(counts):
            if r == 0 or count == 0:
                continue
            amount = count * length
            depth_packet[d][r] += amount
            signed = amount if (r + d) % 2 == 0 else -amount
            t[r] += signed
            if r >= 4:
                tail_by_depth[d] += signed

    t[0] = p0_base
    return p0_base, p0_full, t, depth_packet, tail_by_depth


def direct_newton_packets(base: list[int], far: list[int]) -> list[Fraction]:
    """Slow independent check from p0 on all far subsets."""
    far_count = len(far)
    values: dict[tuple[int, ...], Fraction] = {}
    for r in range(far_count + 1):
        for subset in combinations(range(far_count), r):
            speeds = base + [far[i] for i in subset]
            values[subset] = p0_exact(speeds)

    totals = [Fraction(0) for _ in range(far_count + 1)]
    totals[0] = values[()]
    for r in range(1, far_count + 1):
        total = Fraction(0)
        for subset in combinations(range(far_count), r):
            subtotal = Fraction(0)
            for j in range(r + 1):
                for part in combinations(subset, j):
                    subtotal += ((-1) ** (r - j)) * values[part]
            total += subtotal
        totals[r] = total
    return totals


def fmt(frac: Fraction) -> str:
    return f"{float(frac): .9f} ({frac})"


def report_case(label: str, base: list[int], far: list[int], check_direct: bool = False):
    p0_base, p0_full, packets, depth_packet, tail_by_depth = packet_decomposition(base, far)
    bonf3 = sum(packets[:4])
    tail = p0_full - bonf3
    print(f"\n{label}")
    print(f"  base={base} far={far} total_k={len(base) + len(far)}")
    print(f"  p0(base)={fmt(p0_base)}")
    print(f"  p0(full)={fmt(p0_full)}")
    print(f"  T packets={[fmt(x) for x in packets]}")
    print(f"  Bonferroni3=sum T0..T3={fmt(bonf3)}")
    print(f"  tail>=4=p0-Bonf3={fmt(tail)}  upper_bound_holds={tail <= 0}")
    print("  tail_by_missing_depth:")
    for d in sorted(tail_by_depth):
        if tail_by_depth[d]:
            print(f"    d={d}: {fmt(tail_by_depth[d])}")
    print("  nonzero unsigned depth packets A[d][r]:")
    for d, arr in depth_packet.items():
        if any(arr):
            print(f"    d={d}: {[fmt(x) for x in arr]}")
    if check_direct:
        direct = direct_newton_packets(base, far)
        print(f"  direct subset-p0 packets match formula: {direct == packets}")


def random_stress(samples: int = 80) -> None:
    rng = Random(2903)
    positive = []
    worst_tail = (Fraction(-10), None)
    worst_p0_positive = (Fraction(-1), None)

    for _ in range(samples):
        total = rng.randint(8, 13)
        far_count = rng.randint(4, min(6, total))
        base_count = total - far_count
        base_pool = list(range(1, 15))
        base = [0] + sorted(rng.sample(base_pool, base_count - 1))
        far = sorted(rng.sample(range(15, 38), far_count))
        # Avoid duplicate speed sets with a nonprimitive full gcd only as a diagnostic.
        g = 0
        for speed in base + far:
            g = gcd(g, speed)
        p0_base, p0_full, packets, _, tail_by_depth = packet_decomposition(base, far)
        tail = p0_full - sum(packets[:4])
        if tail > worst_tail[0]:
            worst_tail = (tail, (base, far, p0_full, tail_by_depth, g))
        if tail > 0:
            positive.append((base, far, p0_full, tail, tail_by_depth, g))
            if p0_full > worst_p0_positive[0]:
                worst_p0_positive = (p0_full, (base, far, tail, tail_by_depth, g))

    print("\nRandom exact stress")
    print(f"  samples={samples}")
    print(f"  positive Bonferroni3 tails={len(positive)}")
    print(
        "  max tail:",
        fmt(worst_tail[0]),
        "case=",
        worst_tail[1],
    )
    if worst_p0_positive[1] is not None:
        print(
            "  max p0 among positive-tail cases:",
            fmt(worst_p0_positive[0]),
            "case=",
            worst_p0_positive[1],
        )


def tournament_analysis() -> None:
    vertices = [
        "pointwise_depth_formula",
        "depth_parity_tail_inequality",
        "doublet_plus_triple_cap",
        "bonferroni3_unconditional",
        "raw_venn_containment",
        "raw_runner_vertices",
    ]
    # Criteria: exactness, falsifiability survived, proof usefulness, label retention.
    criteria = {
        "pointwise_depth_formula": (1, 1, 1, 1),
        "depth_parity_tail_inequality": (0, 1, 1, 1),
        "doublet_plus_triple_cap": (0, 1, 1, 1),
        "bonferroni3_unconditional": (0, 0, 1, 1),
        "raw_venn_containment": (0, 0, 0, 1),
        "raw_runner_vertices": (0, 0, 0, 0),
    }
    tie = {v: i for i, v in enumerate(vertices)}

    def score(v: str):
        return sum(criteria[v])

    edges = {v: [] for v in vertices}
    indeg = {v: 0 for v in vertices}
    for i, u in enumerate(vertices):
        for v in vertices[i + 1 :]:
            if score(u) > score(v) or (score(u) == score(v) and tie[u] < tie[v]):
                winner, loser = u, v
            else:
                winner, loser = v, u
            edges[winner].append(loser)
            indeg[loser] += 1

    edge_set = {(u, v) for u, outs in edges.items() for v in outs}
    cycles = 0
    for a, b, c in combinations(vertices, 3):
        if (a, b) in edge_set and (b, c) in edge_set and (c, a) in edge_set:
            cycles += 1
        if (a, c) in edge_set and (c, b) in edge_set and (b, a) in edge_set:
            cycles += 1

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
    print("  vertices: proof carriers, not runners")
    print("  observable: exactness + survived falsification + proof usefulness + label retention")
    print(f"  score_hist={dict(sorted(Counter(indeg.values()).items()))}")
    print(f"  directed_3cycles={cycles}")
    print("  Hamiltonian path=" + " > ".join(topo))


def main() -> None:
    print("S115: Bonferroni-3 needs a missing-depth guard")
    print("=" * 78)
    print("Pointwise formula:")
    print("  T_r(x)=(-1)^(r+d) C_{d,r}(x), d=# base-missed inner sectors.")
    print("  C_{d,r} counts r-subsets of far runners whose colors cover the")
    print("  missing sectors and use no already-covered/zero sector colors.")

    report_case(
        "KPS S31t binding-style sample: tail is negative",
        [0, 1, 2, 3, 4],
        [20, 21, 22, 23, 24],
        check_direct=False,
    )
    report_case(
        "Exact counterexample to unconditional Bonferroni-3, k=8",
        [0, 1, 2, 3],
        [16, 28, 29, 32],
        check_direct=True,
    )
    report_case(
        "Second positive-tail case with base size 5",
        [0, 1, 7, 10, 13],
        [15, 23, 24, 31],
        check_direct=False,
    )

    random_stress()

    print("\nProof correction")
    print("  The Newton expansion is exact, but the high-tail sign is a depth")
    print("  parity inequality, not a consequence of Venn containment alone.")
    print("  A viable proof target is:")
    print("    (positive even-depth high packets) <= (negative odd-depth high packets) + cap slack,")
    print("  with deep-missing rows routed to an easy small-p0/decorrelation bound.")

    tournament_analysis()


if __name__ == "__main__":
    main()
