#!/usr/bin/env python3
"""Exact atlas for HYP-6785 and verification of THM-760.

The exact proof object is the endogenous pair-sum blocker complex.  For every
pair (i,j), including i=j, put q=v_i+v_j and retain every multiplier p mod q
(NONUNITS INCLUDED).  If the protected endpoint runner(s) are safe at a target
c, the remaining unsafe runners form the blocker hyperedge.  THM-668 gives:
M(S)>=c iff some retained hyperedge is empty.

Tournament Analysis is deliberately secondary.  Its vertices are runners and
the pairwise observable is exclusive blocker responsibility over obligations
where neither competitor is protected.  The switch is the sign of that
observable; ties follow the increasing-speed Hamiltonian path.  It preserves
pairwise obstruction pressure but destroys simultaneous hyperedge coverage,
so it is not used as a proof quotient.

Alternate vertices considered: runners, distinct pair-sum rulers, witness
sheets, residues, wall events, and proof obligations.  Proof obligations are
selected for the exact object; witness sheets prove THM-760; runners are kept
only for the requested tournament fingerprints.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from math import gcd
import random


LRC = F(1, 14)
COVER_FLOOR = F(14, 183)
CORE_FLOOR = F(1, 13)


def dist_mod(x: int, q: int) -> int:
    r = x % q
    return min(r, q - r)


def dist_fraction(x: F) -> F:
    r = x.numerator % x.denominator
    return F(min(r, x.denominator - r), x.denominator)


def margin_at(speeds: tuple[int, ...], t: F) -> F:
    return min(dist_fraction(F(v) * t) for v in speeds)


def exact_M(speeds: tuple[int, ...]) -> tuple[F, F, tuple[int, int]]:
    """THM-668 exact maximum; all pair sums and all (possibly nonunit) p."""
    best = F(-1)
    best_t = F(0)
    best_pair = (0, 0)
    n = len(speeds)
    for i in range(n):
        for j in range(i, n):
            q = speeds[i] + speeds[j]
            for p in range(1, q):
                d = min(dist_mod(v * p, q) for v in speeds)
                m = F(d, q)
                if m > best:
                    best, best_t, best_pair = m, F(p, q), (i, j)
    return best, best_t, best_pair


def sheet_witness(core: tuple[int, ...], c: int, w: int) -> tuple[F, F, F, int]:
    """Return core M, chosen lifted time, exceptional clearance, and sheet."""
    assert c >= 2 and gcd(c, w) == 1
    core_M, t0, _ = exact_M(core)
    choices = []
    scaled = tuple(c * p for p in core)
    for k in range(c):
        t = (t0 + k) / c
        choices.append((dist_fraction(F(w) * t), t, k))
        assert margin_at(scaled, t) == margin_at(core, t0)
    exceptional, t, k = max(choices)
    assert exceptional >= F(1, 2) - F(1, 2 * c)
    assert margin_at(tuple(sorted(scaled + (w,))), t) >= min(core_M, exceptional)
    return core_M, t, exceptional, k


def blocker_obligations(speeds: tuple[int, ...], target: F):
    """Yield (protected pair, q, p, blocker frozenset) at endpoint-safe events."""
    n = len(speeds)
    A, B = target.numerator, target.denominator
    for i in range(n):
        for j in range(i, n):
            q = speeds[i] + speeds[j]
            protected = {i, j}
            for p in range(1, q):
                if any(B * dist_mod(speeds[h] * p, q) < A * q for h in protected):
                    continue
                blockers = frozenset(
                    h for h, v in enumerate(speeds)
                    if h not in protected and B * dist_mod(v * p, q) < A * q
                )
                yield (i, j, q, p, blockers)


def tournament(speeds: tuple[int, ...], target: F):
    obs = list(blocker_obligations(speeds, target))
    n = len(speeds)
    edge = [[False] * n for _ in range(n)]
    raw = {}
    for a in range(n):
        for b in range(a + 1, n):
            ab = ba = 0
            for i, j, _q, _p, blockers in obs:
                if a in (i, j) or b in (i, j):
                    continue
                if a in blockers and b not in blockers:
                    ab += 1
                elif b in blockers and a not in blockers:
                    ba += 1
            raw[(a, b)] = ab - ba
            if ab > ba or (ab == ba and speeds[a] < speeds[b]):
                edge[a][b] = True
            else:
                edge[b][a] = True
    scores = [sum(edge[i]) for i in range(n)]
    cycles = sum(
        edge[a][b] and edge[b][c] and edge[c][a]
        or edge[a][c] and edge[c][b] and edge[b][a]
        for a, b, c in combinations(range(n), 3)
    )
    reach = [row[:] for row in edge]
    for i in range(n):
        reach[i][i] = True
    for k in range(n):
        for i in range(n):
            if reach[i][k]:
                for j in range(n):
                    reach[i][j] = reach[i][j] or reach[k][j]
    scc_sizes = sorted(
        Counter(tuple(j for j in range(n) if reach[i][j] and reach[j][i]) for i in range(n)).values(),
        reverse=True,
    )
    # Count Hamiltonian paths by subset DP (n=13).
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1 << n):
        for last in range(n):
            if not dp[mask][last]:
                continue
            for nxt in range(n):
                if not (mask >> nxt) & 1 and edge[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += dp[mask][last]
    hpaths = sum(dp[-1])
    return edge, {
        "score_hist": dict(sorted(Counter(scores).items())),
        "cycles3": cycles,
        "scc_sizes": scc_sizes,
        "hamiltonian_paths": hpaths,
    }


def complex_summary(name: str, speeds0: list[int]):
    speeds = tuple(sorted(speeds0))
    M, t, pair = exact_M(speeds)
    mult = Counter(speeds[i] + speeds[j] for i in range(13) for j in range(i + 1, 13))
    print(f"{name}: M={M} ({float(M):.6f}), witness={t}, pair=({speeds[pair[0]]},{speeds[pair[1]]}), max_pair_mult={max(mult.values())}")
    tournaments = {}
    for target in (LRC, COVER_FLOOR, CORE_FLOOR):
        obs = list(blocker_obligations(speeds, target))
        empty = sum(not row[-1] for row in obs)
        blockers = Counter(len(row[-1]) for row in obs)
        edge, fp = tournament(speeds, target)
        tournaments[target] = edge
        print(f"  c={target}: obligations={len(obs)}, empty={empty}, blocker_sizes={dict(sorted(blockers.items()))}")
        print(f"    tournament={fp}")
    flips = 0
    e0, e1 = tournaments[COVER_FLOOR], tournaments[CORE_FLOOR]
    for i in range(13):
        for j in range(i + 1, 13):
            flips += e0[i][j] != e1[i][j]
    print(f"  edge_flips(14/183 -> 1/13)={flips}")


def main():
    print("THM-760 SHEET-DODGE VERIFICATION")
    rng = random.Random(20260714)
    trials = 0
    tightest = None
    for _ in range(80):
        core = tuple(sorted(rng.sample(range(1, 45), 12)))
        c = rng.randint(2, 24)
        choices = [w for w in range(1, 6 * c + 20) if gcd(c, w) == 1 and w not in {c * p for p in core}]
        w = rng.choice(choices)
        core_M, t, exceptional, _k = sheet_witness(core, c, w)
        S = tuple(sorted(tuple(c * p for p in core) + (w,)))
        got = margin_at(S, t)
        assert got >= min(core_M, F(1, 2) - F(1, 2 * c)) >= CORE_FLOOR
        trials += 1
        row = (got, core, c, w, t, exceptional)
        if tightest is None or row[0] < tightest[0]:
            tightest = row
    print(f"  exact random trials={trials}, failures=0")
    print(f"  smallest exhibited sheet margin={tightest[0]} at c={tightest[2]}, w={tightest[3]}, t={tightest[4]}")
    print("  theorem floor: min(M(core), 1/2-1/(2c)); LRC(<=13) corollary >= 1/13")

    # Exact scale-ray and safe-killer checks central to HYP-6780/THM-757.
    for c, w in ((26, 339), (52, 677), (104, 1353), (26, 15)):
        core = tuple(range(1, 13))
        core_M, t, exceptional, k = sheet_witness(core, c, w)
        got = margin_at(tuple(sorted(tuple(c * p for p in core) + (w,))), t)
        print(f"  scaled-block c={c}, w={w}: sheet={k}, t={t}, exceptional={exceptional}, margin={got}, coreM={core_M}")

    print("\nENDOGENOUS BLOCKER-COMPLEX ATLAS")
    families = [
        ("AP", list(range(1, 14))),
        ("GW", list(range(1, 12)) + [13, 24]),
        ("deep_well", list(range(1, 13)) + [182]),
        ("second_rung", list(range(1, 12)) + [13, 84]),
        ("compressed", [2 * k for k in range(1, 13)] + [13]),
        ("near_dilate_c26", [26 * k for k in range(1, 13)] + [339]),
        ("safe_killer_c26", [26 * k for k in range(1, 13)] + [15]),
        ("spread_blocker", [200, 496, 540, 656, 851, 921, 935, 1122, 1482, 1680, 1835, 1849, 1856]),
    ]
    for name, speeds in families:
        complex_summary(name, speeds)

    print("\nMETHODOLOGY")
    print("  exact predicate retained: an empty protected blocker edge <=> M(S)>=target (THM-668)")
    print("  destroyed by runner tournament: simultaneous cover intersections, exact margins, ruler labels")
    print("  challenged assumption: witness sheets/proof obligations, not runners, are the natural proof vertices")


if __name__ == "__main__":
    main()
