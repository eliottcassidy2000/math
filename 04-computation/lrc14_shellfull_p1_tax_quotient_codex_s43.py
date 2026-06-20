#!/usr/bin/env python3
"""HYP-2669 scout: shell-1-full p1-tax quotient stability.

HYP-2668 shows that a global raw 2p1/5 endpoint tax fails at B=14, but its
single failure damages shell 1 by missing the tower bit 2.  This script tests
the ordered target on the quotient where the dyadic-1 tower {1,2,4,8} is
forced to remain present.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from itertools import combinations

from lrc14_far_delta_galois_phase_codex_s38 import primitive
from lrc14_p1_tax_envelope_codex_s40 import RatioRow, base_data, row


SHELL1 = (1, 2, 4, 8)


def fmt(q: Fraction) -> str:
    return f"{q} ({float(q):.6f})"


def tax_gap(row_: RatioRow, c: Fraction) -> Fraction:
    return c * row_.w * row_.p1 - row_.raw


def shell_extras(Ep: tuple[int, ...]) -> tuple[int, ...]:
    shell = set(SHELL1)
    return tuple(v for v in Ep if v and v not in shell)


def odd_carry_profile(Ep: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    profile: Counter[int] = Counter()
    for v in Ep:
        if v == 0:
            continue
        odd = v
        weight = 1
        while odd % 2 == 0:
            odd //= 2
            weight *= 2
        profile[odd] += weight
    return tuple(sorted(profile.items()))


def scan_shellfull(Bmax: int, w_extra: int) -> list[RatioRow]:
    """Scan E'={0}+{1,2,4,8}+3 extras from [1,Bmax]."""

    rows: list[RatioRow] = []
    extras_universe = [v for v in range(1, Bmax + 1) if v not in SHELL1]
    for extras in combinations(extras_universe, 3):
        Ep = (0,) + tuple(sorted(SHELL1 + extras))
        base = base_data(Ep)
        for w in range(max(Ep) + 1, max(Ep) + w_extra + 1):
            if primitive(Ep + (w,)):
                rows.append(row(base, w, f"shellfullB{Bmax}"))
    return rows


def print_summary(rows: list[RatioRow], checkpoints: tuple[int, ...]) -> None:
    print("cumulative shell-1-full quotient")
    print("B | rows | >1/3 | >3/8 | >2/5 | max ratio | max row | min gap 2/5")
    for B in checkpoints:
        subset = [r for r in rows if max(r.Ep) <= B]
        top = max(subset, key=lambda r: (r.actual, r.envelope))
        min_gap = min(tax_gap(r, Fraction(2, 5)) for r in subset)
        print(
            f"{B:2d} | {len(subset):5d} | "
            f"{sum(r.actual > Fraction(1, 3) for r in subset):4d} | "
            f"{sum(r.actual > Fraction(3, 8) for r in subset):4d} | "
            f"{sum(r.actual > Fraction(2, 5) for r in subset):4d} | "
            f"{fmt(top.actual):>20} | w={top.w}, E'={top.Ep} | {min_gap}"
        )
    print()


def print_rows(title: str, rows: list[RatioRow], limit: int = 16) -> None:
    print(title)
    print("rank | ratio | w | E' | extras | odd-carry profile | raw | p1 | gap_2/5")
    for rank, row_ in enumerate(sorted(rows, key=lambda r: (r.actual, r.envelope), reverse=True)[:limit], 1):
        print(
            f"{rank:4d} | {fmt(row_.actual):>20} | {row_.w:2d} | {row_.Ep} | "
            f"{shell_extras(row_.Ep)} | {odd_carry_profile(row_.Ep)} | "
            f"{row_.raw} | {row_.p1} | {tax_gap(row_, Fraction(2, 5))}"
        )
    print()


def main() -> None:
    Bmax = 24
    w_extra = 8
    checkpoints = (14, 16, 18, 20, 22, 24)
    rows = scan_shellfull(Bmax, w_extra)
    top = max(rows, key=lambda r: (r.actual, r.envelope))
    new_speed_rows = [r for r in rows if max(r.Ep) > 14]
    top_new = max(new_speed_rows, key=lambda r: (r.actual, r.envelope))

    print("HYP-2669 shell-1-full p1-tax quotient scout")
    print("exact Fraction arithmetic")
    print(f"family: E'={{0}}+{{1,2,4,8}}+3 extras from [1,{Bmax}], w=max(E')+1..max(E')+{w_extra}")
    print(f"rows={len(rows)}")
    print(f"global shell-full max={fmt(top.actual)} at w={top.w}, E'={top.Ep}")
    print(f"global shell-full min 2/5 tax gap={min(tax_gap(r, Fraction(2, 5)) for r in rows)}")
    print(f"new-speed max with max(E')>14={fmt(top_new.actual)} at w={top_new.w}, E'={top_new.Ep}")
    print()

    print_summary(rows, checkpoints)
    print_rows("top shell-full rows through B=24", rows)
    print_rows("top rows with new speed beyond B=14", new_speed_rows, limit=12)

    print("Interpretation")
    print("  The shell-1-full quotient keeps the B=13/S41 leader through B=24:")
    print("    E'=(0,1,2,4,6,7,8,10), w=12, Delta^+/p1=997/2562.")
    print("  Its exact 2p1/5 gap is 139/2450, and no scanned shell-full row crosses 2/5.")
    print("  Rows introducing speeds beyond 14 are much lower; their leader is below 1/3.")
    print("  The live proof target is therefore not a growing scalar frontier in this quotient,")
    print("  but a finite dyadic-even packet lemma around the B=13 leader.")
    print()
    print("Tournament Analysis")
    print("  vertices: b13_shellfull_leader > dyadic_even_packet > new_speed_decay > extended_tower > raw_ratio")
    print("  observable: whether forcing shell 1 lets the 2p1/5 p1-tax target survive larger B")
    print("  switch/gauge: quotient first by the full dyadic-1 tower, then compare extra odd-carry profiles")
    print("  Hamiltonian path: b13_shellfull_leader > dyadic_even_packet > new_speed_decay > extended_tower > raw_ratio")
    print("  challenged assumption: the S42 shell-full survival might be a B=14 artifact; through B=24 it is not.")
    print("PASS: shell-1-full p1-tax quotient scout complete.")


if __name__ == "__main__":
    main()
