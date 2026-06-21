#!/usr/bin/env python3
"""HYP-2708 addendum: live-depth ladder for one/many far hits.

After the S23 OPEN-Q-108 localization, the wide p0 problem is mostly
single-far, while HYP-2708 handles the two-far survival deviation.  This
script puts those branches in one exact depth-kernel table.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from math import comb


def survival_currency(t: int) -> F:
    if 1 <= t <= 4:
        return F(1)
    if t == 6:
        return F(-4)
    return F(0)


def cover_currency(t: int) -> F:
    return F(1) if t == 0 else F(0)


def exact_hit_prob(t: int, h: int, r: int) -> F:
    """Probability r iid colors hit exactly h distinct colors from a t-set."""

    if not 0 <= h <= min(t, r):
        return F(0)
    total = F(0)
    # Choose the h missed colors that are hit.  Then force every chosen color
    # to appear at least once while avoiding the t-h missed colors not chosen.
    for j in range(h + 1):
        total += ((-1) ** j) * comb(h, j) * F(7 - t + h - j, 7) ** r
    return F(comb(t, h)) * total


def boundary_value(currency, t: int, r: int) -> F:
    return sum(
        exact_hit_prob(t, h, r) * currency(t - h)
        for h in range(min(t, r) + 1)
    )


def coeffs(currency, t: int, r: int) -> list[F]:
    boundary = boundary_value(currency, t, r)
    return [currency(t - h) - boundary for h in range(min(t, r) + 1)]


def live_depths(currency, r: int) -> tuple[int, ...]:
    live = []
    for t in range(7):
        if any(c != 0 for c in coeffs(currency, t, r)):
            live.append(t)
    return tuple(live)


def show_survival_ladder() -> None:
    print("SURVIVAL CURRENCY LIVE-DEPTH LADDER")
    print("  C(t)=p1+p2+p3+p4-4p6 coefficient by missed depth:")
    print("  " + ", ".join(f"C({t})={survival_currency(t)}" for t in range(7)))
    print()
    for r in range(1, 7):
        live = live_depths(survival_currency, r)
        silent = tuple(t for t in range(7) if t not in live)
        print(f"r={r} far hits: live_depths={live} silent_depths={silent}")
        for t in range(7):
            boundary = boundary_value(survival_currency, t, r)
            cs = coeffs(survival_currency, t, r)
            mark = "LIVE" if any(c != 0 for c in cs) else "silent"
            print(
                f"  t={t}: K{r}={boundary} "
                f"coeffs_h={[str(c) for c in cs]} {mark}"
            )
        print()


def show_cover_single_far() -> None:
    print("DIRECT p0 SINGLE-FAR KERNEL")
    print("  D(t)=1[t=0].  With one inserted far color, only t=1 can change p0.")
    for t in range(7):
        boundary = boundary_value(cover_currency, t, 1)
        cs = coeffs(cover_currency, t, 1)
        mark = "LIVE" if any(c != 0 for c in cs) else "silent"
        print(f"  t={t}: P1={boundary} coeffs_h={[str(c) for c in cs]} {mark}")
    print()


def tournament() -> None:
    branches = [
        ("S23 single-far p0 finite window", live_depths(cover_currency, 1), 1),
        ("one-far survival transfer", live_depths(survival_currency, 1), 1),
        ("HYP-2708 two-far survival", live_depths(survival_currency, 2), 2),
        ("three-far survival aggregate", live_depths(survival_currency, 3), 3),
        ("four-plus synchronizer branch", live_depths(survival_currency, 4), 4),
    ]
    scores = [0] * len(branches)
    edges = set()
    for i, (_, live_i, r_i) in enumerate(branches):
        for j, (_, live_j, r_j) in enumerate(branches):
            if i >= j:
                continue
            key_i = (len(live_i), r_i, i)
            key_j = (len(live_j), r_j, j)
            if key_i <= key_j:
                edges.add((i, j))
                scores[i] += 1
            else:
                edges.add((j, i))
                scores[j] += 1
    cycles = 0
    for a in range(len(branches)):
        for b in range(a + 1, len(branches)):
            for c in range(b + 1, len(branches)):
                if (a, b) in edges and (b, c) in edges and (c, a) in edges:
                    cycles += 1
                if (a, c) in edges and (c, b) in edges and (b, a) in edges:
                    cycles += 1

    print("TOURNAMENT ANALYSIS")
    print("  vertices: live-depth proof branches")
    print("  pairwise observable: fewer live depths, then smaller far-hit count")
    print("  switch/gauge: choose the currency before scalarizing row risk")
    print(f"  score_hist={dict(sorted(Counter(scores).items()))}")
    print(f"  directed_3cycles={cycles}")
    print("  Hamiltonian proof-pressure path:")
    for name, live, r in sorted(branches, key=lambda x: (len(x[1]), x[2], x[0])):
        print(f"    r={r}: {name}; live={live}")
    print()


def main() -> None:
    print("HYP-2708 addendum -- survival live-depth ladder")
    print("Exact arithmetic: closed-form iid hit law on seven colors.\n")
    show_survival_ladder()
    show_cover_single_far()
    tournament()
    print("SYNTHESIS")
    print("  S23's single-far p0 window is a pure t=1 closure problem.")
    print("  The one-far survival gate adds only high-tail t=5,6 debt.")
    print("  HYP-2708's two-far survival target is exactly t=1,2,5,6.")
    print("  For r=3 only t=4 remains silent; for r>=4 every positive depth is live.")


if __name__ == "__main__":
    main()
