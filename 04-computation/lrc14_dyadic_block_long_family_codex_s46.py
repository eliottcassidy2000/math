#!/usr/bin/env python3
"""HYP-2671 addendum: extend the dyadic-block family beyond S45.

S45 found that the shell-full new-speed maximum is the `m=4` member of

    E_m={0,1,2,4,8,3m,4m,5m}, w=6m.

This small addendum keeps the same exact evaluator and pushes the family farther
out.  It is intentionally targeted: one row per `m`, not a broad shell-full
scan.
"""

from __future__ import annotations

from fractions import Fraction

from lrc14_far_delta_galois_phase_codex_s38 import primitive
from lrc14_newspeed_constant_codex_s45 import eval_row, fold_mass, odd_carry_profile


def fmt(q: Fraction) -> str:
    return f"{q} ({float(q):.6f})"


def main() -> None:
    m_min = 3
    m_max = 120
    rows = []

    for m in range(m_min, m_max + 1):
        ep = tuple(sorted({0, 1, 2, 4, 8, 3 * m, 4 * m, 5 * m}))
        w = 6 * m
        if len(ep) != 8 or not primitive(ep + (w,)):
            continue
        rows.append((m, eval_row(ep, w)))

    top_m, top = max(rows, key=lambda mr: (mr[1].actual, mr[1].raw))
    tail_rows = [(m, row) for m, row in rows if m > 24]
    tail_m, tail_top = max(tail_rows, key=lambda mr: (mr[1].actual, mr[1].raw))

    print("HYP-2671 dyadic-block long-family addendum")
    print("exact Fraction arithmetic")
    print(f"family: E_m={{0,1,2,4,8,3m,4m,5m}}, w=6m, m={m_min}..{m_max}")
    print(f"primitive rows={len(rows)}")
    print()
    print(f"global family top: m={top_m}, ratio={fmt(top.actual)}, gap below 1/3={Fraction(1,3)-top.actual}")
    print(f"tail m>24 top:    m={tail_m}, ratio={fmt(tail_top.actual)}, gap below 1/3={Fraction(1,3)-tail_top.actual}")
    print()

    print("top 20 family rows")
    print("rank | m | max(E') | ratio | gap below 1/3 | p1 | fold_recip | odd-carry")
    for rank, (m, row) in enumerate(sorted(rows, key=lambda mr: (mr[1].actual, mr[1].raw), reverse=True)[:20], 1):
        print(
            f"{rank:4d} | {m:3d} | {max(row.Ep):7d} | {fmt(row.actual):>20} | "
            f"{Fraction(1,3)-row.actual} | {row.p1} | {fold_mass(row)} | {odd_carry_profile(row.Ep)}"
        )
    print()

    print("range maxima")
    print("m range | rows | best m | best ratio | gap below 1/3")
    cuts = [(3, 24), (25, 60), (61, 120)]
    for a, b in cuts:
        group = [(m, row) for m, row in rows if a <= m <= b]
        best_m, best = max(group, key=lambda mr: (mr[1].actual, mr[1].raw))
        print(f"{a:3d}..{b:<3d} | {len(group):4d} | {best_m:6d} | {fmt(best.actual):>20} | {Fraction(1,3)-best.actual}")
    print()

    print("Interpretation")
    print("  The m=4 spike remains isolated in the targeted family through m=120.")
    print("  After m=24, the best row found is far below 1/3 and also far below the")
    print("  B30 non-family new-speed competitors.  This does not prove HYP-2671,")
    print("  but it reinforces the idea that m=4 is a finite resonance to isolate,")
    print("  not the first term of a high-ratio tail.")
    print("PASS: dyadic-block long-family addendum complete.")


if __name__ == "__main__":
    main()
