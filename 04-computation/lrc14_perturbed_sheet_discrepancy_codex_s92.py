#!/usr/bin/env python3
"""
lrc14_perturbed_sheet_discrepancy_codex_s92.py

Exact sheet-discrepancy bridge for clustered LRC14 rows.

The S91 pure-dilation proof used that, for S={b,2b,...,12b,V} with
gcd(b,V)=1, the b sheet offsets frac(V*n/b) are exactly equally spaced.
This script records the general replacement:

  For a parked runner W and denominator b, write r = W mod b and g=gcd(r,b).
  On the sheets tau=(n+u)/b, n=0..b-1, the dangerous n's sample one interval
  of length 1/7 on a b/g grid, each point repeated g times.  Hence the danger
  count is at most b/7 + g.

For h parked runners W_l added to an exactly dilated base bE, every u in the
base-good set G_E keeps at least

    b*(1 - h/7) - sum_l gcd(W_l,b)

sheets.  This is the discrete, lossless analogue of an Erdos-Turan/Koksma
bound, with the exact discrepancy term sum gcd(W_l,b).
"""

from __future__ import annotations

from fractions import Fraction as F
from math import ceil, gcd
from random import Random

DELTA = F(1, 14)
DANGER = F(1, 7)


def fmt(x: F) -> str:
    return f"{x.numerator}/{x.denominator}" if x.denominator != 1 else str(x.numerator)


def norm(x: F) -> F:
    x %= 1
    return min(x, 1 - x)


def in_interval(x: F, lo: F, hi: F) -> bool:
    x %= 1
    if lo <= hi:
        return lo <= x < hi
    return x >= lo or x < hi


def danger_count(b: int, r: int, theta: F) -> int:
    return sum(norm(F(r * n, b) + theta) < DELTA for n in range(b))


def interval_count(b: int, r: int, lo: F, hi: F, theta: F = F(0)) -> int:
    return sum(in_interval(F(r * n, b) + theta, lo, hi) for n in range(b))


def discrepancy_bound_check(trials: int = 2000) -> tuple[int, tuple[int, int, F, int, int] | None]:
    rng = Random(92)
    failures = 0
    worst: tuple[int, int, F, int, int] | None = None
    worst_slack = -10**9
    for _ in range(trials):
        b = rng.randint(2, 500)
        r = rng.randint(0, 4 * b)
        if r % b == 0:
            continue
        g = gcd(r, b)
        theta = F(rng.randrange(0, 5000), 5000)
        count = danger_count(b, r, theta)
        bound = F(b, 7) + g
        if F(count) > bound:
            failures += 1
            worst = (b, r, theta, count, g)
            break
        slack_num = count - ceil(bound)
        if slack_num > worst_slack:
            worst_slack = slack_num
            worst = (b, r, theta, count, g)
    return failures, worst


def union_survivors(b: int, residues: list[int], theta_seed: int = 0) -> tuple[int, int]:
    rng = Random(theta_seed)
    thetas = [F(rng.randrange(0, 10007), 10007) for _ in residues]
    killed = set()
    for r, theta in zip(residues, thetas):
        for n in range(b):
            if norm(F(r * n, b) + theta) < DELTA:
                killed.add(n)
    return b - len(killed), len(killed)


def positive_floor(b: int, residues: list[int]) -> F:
    h = len(residues)
    return F(1) - F(h, 7) - F(sum(gcd(r, b) for r in residues), b)


def quotient_gate_table() -> list[tuple[int, F]]:
    rows: list[tuple[int, F]] = []
    for h in range(1, 7):
        # If positivity fails, some gcd(r_i,b) >= b*(7-h)/(7h),
        # so the quotient b/gcd(r_i,b) is at most 7h/(7-h).
        rows.append((h, F(7 * h, 7 - h)))
    return rows


def main() -> None:
    print("LRC14 perturbed sheet-discrepancy bridge (Codex S92)")
    print()
    print("Exact one-runner sheet discrepancy")
    print("  danger count <= b/7 + gcd(r,b)")
    failures, worst = discrepancy_bound_check()
    print(f"  randomized exact checks: failures={failures}")
    if worst is not None:
        b, r, theta, count, g = worst
        print(
            "  tightest sampled row: "
            f"b={b}, r={r}, gcd={g}, theta={fmt(theta)}, "
            f"count={count}, bound={fmt(F(b,7)+g)}"
        )
    print()

    print("Multi-parked-runner floor over an exact base bE")
    print("  surviving sheet fraction >= 1 - h/7 - sum gcd(r_i,b)/b")
    examples = [
        (43, [1, 5, 11, 17, 23, 29]),
        (35, [1, 6, 11, 16, 21]),
        (30, [2, 4, 8, 14, 22]),
        (84, [6, 10, 14, 22, 26]),
        (210, [30, 42, 70, 105]),
    ]
    print("    b residues                      floor       actual_survive/killed")
    for b, residues in examples:
        floor = positive_floor(b, residues)
        survive, killed = union_survivors(b, residues, theta_seed=b + sum(residues))
        print(
            f"  {b:3d} {str(residues):28s} {fmt(floor):>10s} "
            f"{survive:4d}/{killed:<4d}"
        )
    print()

    print("Bounded-quotient fallback when the sheet floor is not positive")
    for h, qmax in quotient_gate_table():
        print(
            f"  h={h}: if 1-h/7-sum(g_i)/b <= 0, "
            f"some quotient b/g_i <= {fmt(qmax)} = {float(qmax):.3f}"
        )
    print()

    print("Proof interface")
    print("  Large quotient / small gcd: sheet discrepancy gives a positive finite-V floor.")
    print("  Small quotient / large gcd: route to a bounded residue atlas; if actual")
    print("  speed ratio is <=13, THM-405 closes immediately, otherwise use the")
    print("  existing bounded-ratio/finitary LRC<=13 handoff.")
    print()

    print("Tournament Analysis")
    print("  vertices: sheet_discrepancy > bounded_quotient_gate > THM405_first_window")
    print("            > finite_residue_atlas > raw_Erdos_Turan_bound > runner_vertices")
    print("  pairwise observable: surviving sheet fraction after parked-runner damage")
    print("  challenged assumption: finite-V requires a continuous equidistribution limit.")


if __name__ == "__main__":
    main()
