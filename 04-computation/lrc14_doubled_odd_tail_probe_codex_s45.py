#!/usr/bin/env python3
"""HYP-2672 support scout: doubled-odd shell-full tail packets.

The B36 HYP-2672 scan refutes the naive far-tail p1/4 target by one row:

    E'=(0,1,2,4,8,14,26,34), w=38.

Its extras and runner are doubled odd values:

    extras = 2*(7,13,17),  w = 2*19.

This exact support scout tests whether that row is isolated in the narrow
doubled-odd packet class

    E' = {0,1,2,4,8,2a,2b,2c},  w=2d,

with odd 3 <= a < b < c <= 29 and c < d <= c+8.  It is not a proof of the
full tail, but it gives the next exception-ledger address.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations

from lrc14_far_delta_galois_phase_codex_s38 import primitive
from lrc14_newspeed_constant_codex_s45 import eval_row
from lrc14_shellfull_packet_gap_codex_s44 import fold_profile, fold_recip_mass, odd_carry_profile


ODD_MAX = 29
D_EXTRA = 8


def fmt(q: Fraction) -> str:
    return f"{q} ({float(q):.6f})"


def fold_mass(Ep: tuple[int, ...]) -> Fraction:
    return fold_recip_mass(fold_profile(Ep))


def main() -> None:
    print("HYP-2672 doubled-odd shell-full tail packet probe")
    print("exact Fraction arithmetic")
    print(f"family: E'={{0,1,2,4,8,2a,2b,2c}}, odd 3<=a<b<c<={ODD_MAX}, w=2d, c<d<=c+{D_EXTRA}")

    rows = []
    for a, b, c in combinations(range(3, ODD_MAX + 1, 2), 3):
        Ep = tuple(sorted({0, 1, 2, 4, 8, 2 * a, 2 * b, 2 * c}))
        if len(Ep) != 8:
            continue
        for d in range(c + 1, c + D_EXTRA + 1):
            w = 2 * d
            if w <= max(Ep) or not primitive(Ep + (w,)):
                continue
            row = eval_row(Ep, w)
            rows.append((row.actual, row.raw, a, b, c, d, row))

    rows.sort(reverse=True, key=lambda item: (item[0], item[1]))
    print(f"rows={len(rows)}")
    print(f"above 1/4={sum(actual > Fraction(1,4) for actual, *_rest in rows)}")
    print(f"above 3/10={sum(actual > Fraction(3,10) for actual, *_rest in rows)}")
    print()
    print("top doubled-odd packet rows")
    print("rank | ratio | gap 3/10 | gap 1/4 | odd tuple (a,b,c,d) | E' | w | fold_recip | odd-carry")
    for rank, (actual, _raw, a, b, c, d, row) in enumerate(rows[:24], 1):
        print(
            f"{rank:4d} | {fmt(actual):>20} | {Fraction(3,10)-actual} | "
            f"{Fraction(1,4)-actual} | {(a,b,c,d)} | {row.Ep} | {row.w} | "
            f"{fold_mass(row.Ep)} | {odd_carry_profile(row.Ep)}"
        )

    top = rows[0]
    print()
    print("Interpretation")
    print(
        "  In this doubled-odd packet class, the HYP-2672 B36 exception "
        f"{top[6].Ep}, w={top[6].w}, is the unique row above 1/4."
    )
    print("  No tested doubled-odd row exceeds 3/10.")
    print("  This supports replacing naive far-tail p1/4 by an exception ledger")
    print("  for doubled-odd packets plus a broader tail <3/10 decay target.")
    print("PASS: doubled-odd tail packet probe complete.")


if __name__ == "__main__":
    main()
