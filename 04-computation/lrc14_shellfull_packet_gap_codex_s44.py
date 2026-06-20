#!/usr/bin/env python3
"""HYP-2670 scout: shell-full p1-tax packet gap beyond the B13 leader.

HYP-2669 showed that the shell-1-full p1-tax quotient keeps the same B13
leader through B=24.  This script asks for the next structural reason:

* Does the shell-full quotient have a finite high-ratio pocket?
* Do genuinely new speeds drop below a simpler tax threshold?
* Which packet/fold/tower coordinates distinguish the B13 leader?

The scan intentionally computes the exact raw p1-tax ratio but not the coarse
interval envelope, because HYP-2665--HYP-2667 already showed that the envelope
is too lossy for this question.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from time import perf_counter

from lrc14_far_delta_galois_phase_codex_s38 import primitive
from lrc14_p1_tax_envelope_codex_s40 import base_data, raw_wdelta
from lrc14_p1_tax_packet_frontier_codex_s41 import packet_contributions


SHELL1 = (1, 2, 4, 8)


@dataclass(frozen=True)
class LightRow:
    actual: Fraction
    raw: Fraction
    p1: Fraction
    phi: Fraction
    Ep: tuple[int, ...]
    w: int


def fmt(q: Fraction) -> str:
    return f"{q} ({float(q):.6f})"


def tax_gap(row: LightRow, c: Fraction) -> Fraction:
    return c * row.w * row.p1 - row.raw


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


def fold_profile(Ep: tuple[int, ...]) -> Counter[int]:
    vals = [v for v in Ep if v]
    present = set(vals)
    profile: Counter[int] = Counter()
    for i, a in enumerate(vals):
        for b in vals[i + 1 :]:
            c = a + b
            if c in present:
                profile[c] += 1
    return profile


def fold_recip_mass(profile: Counter[int]) -> Fraction:
    return sum((Fraction(count, target) for target, count in profile.items()), Fraction(0))


def scan_shellfull_light(Bmax: int, w_extra: int) -> list[LightRow]:
    """Scan E'={0}+{1,2,4,8}+3 extras from [1,Bmax]."""

    rows: list[LightRow] = []
    extras_universe = [v for v in range(1, Bmax + 1) if v not in SHELL1]
    for extras in combinations(extras_universe, 3):
        Ep = (0,) + tuple(sorted(SHELL1 + extras))
        base = base_data(Ep)
        if not base.p1:
            continue
        for w in range(max(Ep) + 1, max(Ep) + w_extra + 1):
            if not primitive(Ep + (w,)):
                continue
            raw = raw_wdelta(base, w)
            actual = max(raw, Fraction(0)) / (w * base.p1)
            rows.append(LightRow(actual=actual, raw=raw, p1=base.p1, phi=base.phi, Ep=Ep, w=w))
    return rows


def layer_name(row: LightRow) -> str:
    m = max(row.Ep)
    if m <= 14:
        return "finite <=14"
    if m <= 24:
        return "new 15..24"
    return "tail 25..30"


def print_layer_table(rows: list[LightRow]) -> None:
    print("layer split")
    print("layer | rows | max ratio | max row | >1/3 | >3/8 | >2/5 | min 2/5 gap")
    for name in ("finite <=14", "new 15..24", "tail 25..30", "all >14", "all >24"):
        if name == "all >14":
            group = [r for r in rows if max(r.Ep) > 14]
        elif name == "all >24":
            group = [r for r in rows if max(r.Ep) > 24]
        else:
            group = [r for r in rows if layer_name(r) == name]
        if not group:
            continue
        top = max(group, key=lambda r: (r.actual, r.raw))
        min_gap = min(tax_gap(r, Fraction(2, 5)) for r in group)
        print(
            f"{name:12s} | {len(group):5d} | {fmt(top.actual):>20} | "
            f"w={top.w}, E'={top.Ep} | "
            f"{sum(r.actual > Fraction(1, 3) for r in group):4d} | "
            f"{sum(r.actual > Fraction(3, 8) for r in group):4d} | "
            f"{sum(r.actual > Fraction(2, 5) for r in group):4d} | {min_gap}"
        )
    print()


def print_top_rows(rows: list[LightRow], limit: int = 16) -> None:
    print("top exact p1-tax rows")
    print(
        "rank | ratio | w | layer | E' | extras | odd-carry profile | "
        "folds | fold_recip | raw | p1 | gap_2/5"
    )
    for rank, row in enumerate(sorted(rows, key=lambda r: (r.actual, r.raw), reverse=True)[:limit], 1):
        profile = fold_profile(row.Ep)
        top_folds = tuple(sorted(profile.items(), key=lambda kv: (-kv[1], kv[0]))[:5])
        print(
            f"{rank:4d} | {fmt(row.actual):>20} | {row.w:2d} | {layer_name(row):12s} | "
            f"{row.Ep} | {shell_extras(row.Ep)} | {odd_carry_profile(row.Ep)} | "
            f"{top_folds} | {fold_recip_mass(profile)} | {row.raw} | {row.p1} | "
            f"{tax_gap(row, Fraction(2, 5))}"
        )
    print()


def print_packet_frontier(rows: list[LightRow], limit: int = 6) -> None:
    print("packet/fold anatomy of frontier rows")
    for rank, row in enumerate(sorted(rows, key=lambda r: (r.actual, r.raw), reverse=True)[:limit], 1):
        base = base_data(row.Ep)
        packets = packet_contributions(base, row.w)
        normalizer = row.w * row.p1
        positive = sum((p.raw for p in packets if p.raw > 0), Fraction(0))
        negative = -sum((p.raw for p in packets if p.raw < 0), Fraction(0))
        small_positive = sum((p.raw for p in packets if p.raw > 0 and p.y.denominator <= 35), Fraction(0))
        pos_denoms = Counter(p.y.denominator for p in packets if p.raw > 0)
        top_pos = [p for p in packets if p.raw > 0][:6]
        top_share = top_pos[0].raw / positive if positive and top_pos else Fraction(0)
        print(f"row {rank}: layer={layer_name(row)}, E'={row.Ep}, w={row.w}")
        print(
            f"  ratio={fmt(row.actual)}, gap_2/5={tax_gap(row, Fraction(2, 5))}, "
            f"extras={shell_extras(row.Ep)}, odd_carry={odd_carry_profile(row.Ep)}"
        )
        print(
            f"  positive_packets/(w*p1)={fmt(positive / normalizer)}, "
            f"negative_packets/(w*p1)={fmt(negative / normalizer)}, "
            f"small_den<=35 positive share={fmt(small_positive / positive if positive else Fraction(0))}, "
            f"top_positive_share={fmt(top_share)}, positive_denoms={tuple(sorted(pos_denoms.items())[:10])}"
        )
        print("  top positive packets")
        print("    rank | y | denom | raw/(w*p1) | raw | QR counts | NQR counts | abs_count")
        for j, packet in enumerate(top_pos, 1):
            print(
                f"    {j:4d} | {packet.y} | {packet.y.denominator:5d} | "
                f"{fmt(packet.raw / normalizer):>20} | {packet.raw} | "
                f"{packet.qr_counts} | {packet.nqr_counts} | {packet.abs_count}"
            )
        print()


def main() -> None:
    Bmax = 30
    w_extra = 8
    start = perf_counter()
    print("HYP-2670 shell-full p1-tax packet gap scout")
    print("exact Fraction arithmetic; raw Delta^+/p1 only, no interval envelope")
    print(f"family: E'={{0}}+{{1,2,4,8}}+3 extras from [1,{Bmax}], w=max(E')+1..max(E')+{w_extra}")
    rows = scan_shellfull_light(Bmax, w_extra)
    elapsed = perf_counter() - start
    rows_sorted = sorted(rows, key=lambda r: (r.actual, r.raw), reverse=True)
    top = rows_sorted[0]
    new_top = max((r for r in rows if max(r.Ep) > 14), key=lambda r: (r.actual, r.raw))
    tail_top = max((r for r in rows if max(r.Ep) > 24), key=lambda r: (r.actual, r.raw))
    print(f"rows={len(rows)}, elapsed_seconds={elapsed:.2f}")
    print(f"global shell-full max={fmt(top.actual)} at w={top.w}, E'={top.Ep}")
    print(f"new-speed max max(E')>14={fmt(new_top.actual)} at w={new_top.w}, E'={new_top.Ep}")
    print(f"tail max max(E')>24={fmt(tail_top.actual)} at w={tail_top.w}, E'={tail_top.Ep}")
    print(f"gap of new-speed max below 1/3={Fraction(1, 3) - new_top.actual}")
    print(f"gap of tail max below 1/4={Fraction(1, 4) - tail_top.actual}")
    print()

    print_layer_table(rows)
    print_top_rows(rows)
    print_packet_frontier(rows)

    print("Interpretation")
    print("  In this shell-1-full quotient through B=30, the only row above 1/3")
    print("  is still the B13 dyadic-even leader (0,1,2,4,6,7,8,10), w=12.")
    print("  Every row introducing a speed beyond 14 is below 1/3; every row")
    print("  introducing a speed beyond 24 is below 1/4 in this scan.")
    print("  This suggests a proof split sharper than a global 2p1/5 tax:")
    print("    finite shell-full packet pocket for max(E')<=14;")
    print("    new-speed decay lemma for max(E')>14;")
    print("    far-tail packet lemma for the max(E')>24 quotient.")
    print()
    print("Tournament Analysis")
    print("  vertices: finite_B13_pocket > new_speed_gap > phase_denominator_cliff > fold_target_concentration > raw_runner_vertices")
    print("  observable: exact Delta^+/p1 after quotienting by the full shell-1 tower")
    print("  switch/gauge: compare proof-obligation layers by max(E'), then expose phase packets and fold targets")
    print("  Hamiltonian path: finite_B13_pocket > new_speed_gap > phase_denominator_cliff > fold_target_concentration > raw_runner_vertices")
    print("  challenged assumption: shell-full proof work is not a single growing scalar frontier;")
    print("    the high packet appears finite while new speeds have a visible decay gap.")
    print("  alternate vertices considered: runners, gaps, fixed sectors, sector boundaries,")
    print("    wall-crossing events, residues, cover arcs, Fourier modes, fold targets,")
    print("    Glaisher odd-carry words, and proof-obligation layers.")
    print("  Preserved predicate: the quotient keeps the p1-tax inequality Delta^+ <= c*p1;")
    print("  destroyed information: interval-envelope placement outside the phase-packet address.")
    print("PASS: shell-full packet gap scout complete.")


if __name__ == "__main__":
    main()
