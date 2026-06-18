#!/usr/bin/env python3
"""
lrc14_integer_sequence_carriers_codex.py

Codex 2026-06-18: integer-sequence carriers for the LRC(14) S3 proof.

This builds on HYP-2597's denominator-cap reframe.  The universal good centers
for the pure cluster floor are the rationals with denominator b < 7/2:

    0, 1/2, 1/3, 2/3.

The new question is how those centers survive the small-speed safe set G_P.
That produces a simple binomial sequence:

  * center 1/2 survives exactly when P is all odd;
  * centers 1/3 and 2/3 survive exactly when P has no multiple of 3;
  * the asymptotic universal-center skeleton is dead exactly when P contains
    at least one even speed and at least one multiple of 3.

For |P|=s the survivor count is

    C(7,s) + C(9,s) - C(5,s),

because {1..13} has 7 odds, 9 nonmultiples of 3, and 5 odd nonmultiples of 3.
The complement is the mixed parity/triadic residual.  This is a proof-useful
integer sequence, not merely decoration: it says exactly which small parts are
handled by the fixed denominator-cap skeleton, and which require the genuinely
Diophantine recurring intervals / colored-resonance machinery.

Tournament Analysis declaration.
  Vertex set: proof carriers, not runners:
    triadic_centers, half_center, both_center_families, mixed_residual,
    AP_mu_sequence, colored_resonance_tail.
  Pairwise observable: (proved status, uniformity in spread, coefficient size,
    residual shrinkage).  The gauge is lexicographic, ties following the listed
    Hamiltonian path.
  Fingerprints: score histogram, directed 3-cycles, SCCs, Hamiltonian paths.

Assumption challenge.
  I considered runners, offsets E, center denominators, small-part divisibility
  types, residue colors, safe components, and proof carriers as vertices.  The
  chosen quotient preserves whether the universal good centers produce a
  rigorous 1/R reservoir inside G_P.  It destroys exact finite q=14V alignment,
  so it cannot by itself prove LRC(14); it isolates the bounded-spread skeleton
  from the large-spread colored discrepancy remainder.
"""

from __future__ import annotations

import itertools
import math
import sys
from collections import Counter
from fractions import Fraction as F

sys.path.insert(0, "04-computation")
import lrc14_global_threshold_ladder_codex as ladder  # noqa: E402


def binom(n: int, k: int) -> int:
    if k < 0 or k > n:
        return 0
    return math.comb(n, k)


def small_part_type(P: tuple[int, ...]) -> tuple[bool, bool]:
    all_odd = all(p % 2 for p in P)
    no_three = all(p % 3 for p in P)
    return all_odd, no_three


def survivor_count(s: int) -> int:
    return binom(7, s) + binom(9, s) - binom(5, s)


def universal_coeff(P: tuple[int, ...]) -> F:
    """
    For R>=4 and max(E)=R, the conservative universal center windows contribute
    at least coeff/R to meas(G_P cap Good_E).

    all odd P: center 1/2 contributes length 3/(14R).
    3-free P: centers 1/3 and 2/3 contribute total length 2/(21R).
    """
    all_odd, no_three = small_part_type(P)
    coeff = F(0)
    if all_odd:
        coeff += F(3, 14)
    if no_three:
        coeff += F(2, 21)
    return coeff


def wrap_intervals(intervals: list[tuple[F, F]]) -> list[tuple[F, F]]:
    out: list[tuple[F, F]] = []
    for a, b in intervals:
        if a < 0:
            out.append((F(0), b))
            out.append((a + 1, F(1)))
        elif b > 1:
            out.append((a, F(1)))
            out.append((F(0), b - 1))
        else:
            out.append((a, b))
    return ladder.merge(out)


def universal_windows(R: int) -> list[tuple[F, F]]:
    """
    Conservative HYP-2590/HYP-2597 windows for every offset shape E with maxE=R.
    """
    r = F(R, 1)
    return wrap_intervals(
        [
            (F(0), F(5, 7) / r),
            (F(1) - F(5, 7) / r, F(1)),
            (F(1, 2) - F(3, 28) / r, F(1, 2) + F(3, 28) / r),
            (F(1, 3) - F(1, 42) / r, F(1, 3) + F(1, 42) / r),
            (F(2, 3) - F(1, 42) / r, F(2, 3) + F(1, 42) / r),
        ]
    )


def skeleton_mass(P: tuple[int, ...], R: int) -> F:
    return ladder.measure(ladder.intersect(ladder.safe_set(P), universal_windows(R)))


def hamiltonian_paths_count(adj: list[list[bool]]) -> int:
    n = len(adj)
    dp = [Counter() for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for v, count in list(dp[mask].items()):
            if not count:
                continue
            for w in range(n):
                if not (mask & (1 << w)) and adj[v][w]:
                    dp[mask | (1 << w)][w] += count
    return sum(dp[-1].values())


def count_directed_3cycles(adj: list[list[bool]]) -> int:
    total = 0
    for a, b, c in itertools.combinations(range(len(adj)), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (adj[a][c] and adj[c][b] and adj[b][a]):
            total += 1
    return total


def scc_sizes(adj: list[list[bool]]) -> list[int]:
    n = len(adj)
    radj = [[adj[j][i] for j in range(n)] for i in range(n)]
    seen: set[int] = set()
    order: list[int] = []

    def dfs(v: int) -> None:
        seen.add(v)
        for w, ok in enumerate(adj[v]):
            if ok and w not in seen:
                dfs(w)
        order.append(v)

    for v in range(n):
        if v not in seen:
            dfs(v)
    seen.clear()
    out = []
    for start in reversed(order):
        if start in seen:
            continue
        stack = [start]
        seen.add(start)
        size = 0
        while stack:
            v = stack.pop()
            size += 1
            for w, ok in enumerate(radj[v]):
                if ok and w not in seen:
                    seen.add(w)
                    stack.append(w)
        out.append(size)
    return sorted(out, reverse=True)


def proof_carrier_tournament() -> dict:
    names = [
        "triadic_centers",
        "half_center",
        "both_center_families",
        "mixed_residual",
        "AP_mu_sequence",
        "colored_resonance_tail",
    ]
    # (proved, spread-uniform, coefficient/proof power, residual shrinkage)
    metrics = {
        "triadic_centers": (1, 0, 2, 3),          # coefficient 2/21
        "half_center": (1, 0, 3, 2),              # coefficient 3/14
        "both_center_families": (1, 0, 4, 4),     # coefficient 13/42
        "mixed_residual": (1, -1, 0, 0),          # classified but not solved
        "AP_mu_sequence": (1, 0, 3, 1),           # closed-form AP skeleton
        "colored_resonance_tail": (0, 1, 5, 5),   # not proved, but only route to uniform tail
    }
    n = len(names)
    adj = [[False] * n for _ in range(n)]
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            if i == j:
                continue
            adj[i][j] = metrics[a] > metrics[b] or (metrics[a] == metrics[b] and i < j)
    scores = {names[i]: sum(adj[i]) for i in range(n)}
    return {
        "scores": scores,
        "hist": dict(sorted(Counter(scores.values()).items())),
        "cycles3": count_directed_3cycles(adj),
        "sccs": scc_sizes(adj),
        "hamiltonian_paths": hamiltonian_paths_count(adj),
        "leader": max(scores, key=scores.get),
    }


def main() -> None:
    print("=" * 88)
    print("LRC14 integer-sequence carriers: universal centers versus mixed residual")
    print("=" * 88)

    print("\n1. Center-survivor sequence by small-part size s=|P|")
    print("   survivor(s)=C(7,s)+C(9,s)-C(5,s); mixed=C(13,s)-survivor")
    print("   s | survivor | mixed | total | cluster size k=13-s")
    for s in range(14):
        surv = survivor_count(s)
        total = binom(13, s)
        print(f"  {s:2d} | {surv:8d} | {total-surv:5d} | {total:5d} | {13-s:2d}")

    print("\n2. Asymptotic center coefficients")
    print("   For R=max(E)>=4, the universal-center reservoir has coefficient coeff/R.")
    reps = {
        "both odd and 3-free": (1, 5, 7, 11, 13),
        "odd but has multiple of 3": (1, 3, 5),
        "3-free but has even": (1, 2, 5),
        "mixed residual": (1, 2, 3),
    }
    for name, P in reps.items():
        coeff = universal_coeff(P)
        print(f"  {name:26s} P={P}: coeff={coeff} = {float(coeff):.6f}")

    print("\n3. Conservative-window stabilization check")
    print("   Positive skeleton counts converge to the survivor sequence as R grows.")
    sample_R = [4, 7, 10, 20, 40, 80, 160]
    selected_s = [0, 1, 2, 3, 4, 5, 8, 9, 10, 11, 12, 13]
    for R in sample_R:
        print(f"\n   R={R}")
        for s in selected_s:
            pos = 0
            mixed_pos = 0
            min_mass: F | None = None
            min_P: tuple[int, ...] | None = None
            for P in itertools.combinations(range(1, 14), s):
                P = tuple(P)
                mass = skeleton_mass(P, R)
                if mass > 0:
                    pos += 1
                    if min_mass is None or mass < min_mass:
                        min_mass = mass
                        min_P = P
                if universal_coeff(P) == 0 and mass > 0:
                    mixed_pos += 1
            print(
                f"    s={s:2d}: positive={pos:4d}, formula_survivors={survivor_count(s):4d}, "
                f"mixed_positive={mixed_pos:4d}, min={min_mass if min_mass is not None else 0} at {min_P}"
            )

    print("\n4. LRC proof readout by cluster size")
    print("   k=|E|, s=13-k.  The mixed column is the residual needing non-universal recurrence.")
    for k in range(3, 14):
        s = 13 - k
        total = binom(13, s)
        surv = survivor_count(s)
        print(f"  k={k:2d}: center-surviving P={surv:4d}, mixed residual P={total-surv:4d}, total P={total:4d}")

    print("\n5. Named hard rows classified")
    rows = [
        ("HYP2595_constant100_failure", (1, 2, 11)),
        ("quarter_min", (1, 2, 3)),
        ("near_via_min", (1, 2, 3, 11)),
        ("via_zero_k7", (1, 2, 3, 6, 12, 13)),
        ("hard_core_84m", tuple(list(range(1, 12)) + [13])),
        ("P_empty", ()),
    ]
    for label, P in rows:
        all_odd, no_three = small_part_type(P)
        print(
            f"  {label:28s}: |P|={len(P):2d}, all_odd={all_odd}, no_mult3={no_three}, "
            f"coeff={universal_coeff(P)}"
        )

    print("\n6. Proof-carrier tournament")
    tour = proof_carrier_tournament()
    print(f"  scores={tour['scores']}")
    print(
        f"  score_hist={tour['hist']}, cycles3={tour['cycles3']}, sccs={tour['sccs']}, "
        f"Hamiltonian_paths={tour['hamiltonian_paths']}, leader={tour['leader']}"
    )

    print("\n7. Takeaway")
    print("  The integer sequence survivor(s)=C(7,s)+C(9,s)-C(5,s) is the exact")
    print("  small-part split for the universal denominator-cap skeleton.")
    print("  It proves a c/R reservoir for all-odd or 3-free small parts, and")
    print("  identifies the mixed parity/triadic residual as the real Diophantine tail.")
    print("  This does not finish LRC(14); it makes the remaining proof obligation smaller")
    print("  and cleaner: solve the mixed residual by colored resonance / spread recurrence.")


if __name__ == "__main__":
    main()
