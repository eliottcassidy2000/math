#!/usr/bin/env python3
"""HYP-2665 scout: stress-test the p1 tax and its interval envelope.

HYP-2664 proposes that the positive far discrepancy satisfies

    Delta_w^+ <= p1(E')/3.

This script tests two things:

1. Broad exact evidence for the raw ratio Delta_w^+/p1(E').
2. Whether the elementary interval envelope for the centered sector indicator
   already proves the p1/3 tax.

For a single-missed cell [a,b] with missed sector s, the contribution is

    integral_a^b (1_{s/7 <= frac(w*x) < (s+1)/7} - 1/7) dx.

Ignoring phase and keeping only cell length L=b-a, the largest possible
positive contribution is

    envelope(w*L)/w,

where for r in [0,1) the max interval integral of 1_[0,1/7)-1/7 over length r
is r*6/7 for r<=1/7 and (1-r)/7 for r>=1/7.

If sum envelope(w*L)/w <= p1/3, the p1 tax has a direct coarea proof.  If not,
the excess identifies exactly where phase-packet cancellation is needed.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import reduce
from itertools import combinations
from math import gcd

from lrc14_far_delta_galois_phase_codex_s38 import G0, endpoint_terms, frac, missed_sector_cells, primitive
from lrc14_residual_plateau_packet_codex_s39 import CAP9, missed_distribution


@dataclass(frozen=True)
class BaseData:
    Ep: tuple[int, ...]
    cells: tuple[tuple[Fraction, Fraction, int], ...]
    terms: tuple[tuple[Fraction, int, int], ...]
    p1: Fraction
    phi: Fraction


@dataclass(frozen=True)
class RatioRow:
    actual: Fraction
    envelope: Fraction
    raw: Fraction
    p1: Fraction
    phi: Fraction
    cap_slack_if_third: Fraction
    Ep: tuple[int, ...]
    w: int
    cells: int
    label: str


def envelope_mass(t: Fraction) -> Fraction:
    r = frac(t)
    if r <= Fraction(1, 7):
        return Fraction(6, 7) * r
    return Fraction(1, 7) * (1 - r)


def cell_envelope(Ep: tuple[int, ...], w: int, cells=None) -> Fraction:
    if cells is None:
        cells = missed_sector_cells(Ep)
    total = Fraction(0)
    for lo, hi, _s in cells:
        total += envelope_mass(w * (hi - lo)) / w
    return total


def base_data(Ep: tuple[int, ...]) -> BaseData:
    cells = missed_sector_cells(Ep)
    dist = missed_distribution(Ep)
    p1 = dist[1]
    phi = dist[0] + p1 / 7
    return BaseData(Ep=Ep, cells=cells, terms=endpoint_terms(cells), p1=p1, phi=phi)


def raw_wdelta(base: BaseData, w: int) -> Fraction:
    total = Fraction(0)
    for x, s, coeff in base.terms:
        total += coeff * G0(w * x - Fraction(s, 7))
    return total


def row(base: BaseData, w: int, label: str) -> RatioRow:
    p1 = base.p1
    phi = base.phi
    raw = raw_wdelta(base, w)
    actual = max(raw, Fraction(0)) / (w * p1) if p1 else Fraction(0)
    env = cell_envelope(base.Ep, w, base.cells)
    env_ratio = env / p1 if p1 else Fraction(0)
    return RatioRow(
        actual=actual,
        envelope=env_ratio,
        raw=raw,
        p1=p1,
        phi=phi,
        cap_slack_if_third=CAP9 - phi - p1 / 3,
        Ep=base.Ep,
        w=w,
        cells=len(base.cells),
        label=label,
    )


def fmt(q: Fraction) -> str:
    return f"{q} ({float(q):.6f})"


def print_table(title: str, rows: list[RatioRow], limit: int = 12) -> None:
    print(title)
    print("rank | actual Delta+/p1 | envelope/p1 | raw wDelta | p1 | slack_third | w | cells | label | E'")
    for i, r in enumerate(rows[:limit], 1):
        print(
            f"{i:4d} | {fmt(r.actual):>22} | {fmt(r.envelope):>18} | {str(r.raw):>14} | "
            f"{str(r.p1):>12} | {str(r.cap_slack_if_third):>12} | {r.w:4d} | {r.cells:5d} | "
            f"{r.label} | {r.Ep}"
        )
    print()


def bounded_bank(B: int, w_extra: int, max_rows: int | None = None) -> list[RatioRow]:
    rows: list[RatioRow] = []
    count = 0
    for comb in combinations(range(1, B + 1), 7):
        Ep = (0,) + comb
        if reduce(gcd, Ep, 0) != 1:
            continue
        base = base_data(Ep)
        start = max(Ep) + 1
        stop = max(Ep) + w_extra
        for w in range(start, stop + 1):
            if not primitive(Ep + (w,)):
                continue
            rows.append(row(base, w, f"B{B}"))
        count += 1
        if max_rows is not None and count >= max_rows:
            break
    return rows


def structured_rows() -> list[RatioRow]:
    out: list[RatioRow] = []
    cases: list[tuple[str, tuple[int, ...], range]] = [
        ("consec8", tuple(range(8)), range(8, 220)),
        ("odd_struct", (0, 1, 3, 5, 7, 9, 10, 11), range(12, 220)),
        ("two_scale_20_40", (0, 1, 2, 20, 21, 22, 40), range(41, 180)),
        ("multiscale_30_60", (0, 1, 2, 30, 31, 32, 60, 61), range(62, 240)),
        ("cluster_40_80", (0, 1, 2, 40, 41, 42, 80, 81), range(82, 260)),
    ]
    for m in [10, 20, 30, 40, 50, 60, 70, 90, 110]:
        cases.append((f"scale_M={m}", (0, 1, 2, m, m + 1, m + 2, 2 * m, 2 * m + 1), range(2 * m + 2, 2 * m + 25)))
    for label, Ep, wrange in cases:
        base = base_data(Ep)
        for w in wrange:
            if primitive(Ep + (w,)):
                out.append(row(base, w, label))
    return out


def summarize(name: str, rows: list[RatioRow]) -> None:
    rows_actual = sorted(rows, key=lambda r: (r.actual, r.envelope), reverse=True)
    rows_env = sorted(rows, key=lambda r: (r.envelope, r.actual), reverse=True)
    violations = [r for r in rows if r.actual > Fraction(1, 3)]
    env_fail = [r for r in rows if r.envelope > Fraction(1, 3)]
    slack_fail = [r for r in rows if r.cap_slack_if_third < 0]
    print(f"=== {name} ===")
    print(f"rows={len(rows)}")
    print(f"actual Delta+/p1 violations of 1/3: {len(violations)}")
    print(f"envelope/p1 above 1/3: {len(env_fail)}")
    print(f"bounded-bank slack failures for p1/3 tax: {len(slack_fail)}")
    print(f"max actual={fmt(rows_actual[0].actual)} at {rows_actual[0].label}, w={rows_actual[0].w}, E'={rows_actual[0].Ep}")
    print(f"max envelope={fmt(rows_env[0].envelope)} at {rows_env[0].label}, w={rows_env[0].w}, E'={rows_env[0].Ep}")
    print()
    print_table("top actual ratios", rows_actual)
    print_table("top envelope ratios", rows_env)


def main() -> None:
    print("HYP-2665 p1-tax envelope scout")
    print("exact Fraction arithmetic")
    print()

    structured = structured_rows()
    summarize("structured resonant families", structured)

    bank13 = bounded_bank(13, 8, max_rows=250)
    summarize("targeted bounded bank B=13, first 250 primitive bases, w=max+1..max+8", bank13)

    bank16 = bounded_bank(16, 5, max_rows=150)
    summarize("targeted bounded bank B=16, first 150 primitive bases, w window 5", bank16)

    all_rows = structured + bank13 + bank16
    max_actual = max(r.actual for r in all_rows)
    max_env = max(r.envelope for r in all_rows)
    print("global summary")
    print(f"  total rows={len(all_rows)}")
    print(f"  max actual Delta+/p1 = {fmt(max_actual)}")
    print(f"  max envelope/p1      = {fmt(max_env)}")
    print("  direct interval envelope proves p1/3 exactly when envelope/p1<=1/3;")
    print("  rows with envelope excess but actual below 1/3 are the packet-cancellation target.")
    print()
    print("Tournament Analysis")
    print("  vertices: actual_p1_tax > interval_envelope > phase_packet_cancellation > endpoint_count > raw_speed")
    print("  pairwise observable: which bound certifies Delta_w^+ <= p1/3")
    print("  switch/gauge: replace endpoint count by interval-length envelope, then by phase packets")
    print("  challenged assumption: p1 tax is not automatically an interval-length fact;")
    print("    some rows require phase-packet cancellation beyond the coarse envelope.")
    print("PASS: p1-tax envelope scout complete.")


if __name__ == "__main__":
    main()
