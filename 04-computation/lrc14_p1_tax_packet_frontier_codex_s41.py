#!/usr/bin/env python3
"""HYP-2667 scout: full B=13 p1-tax frontier and packet anatomy.

S40 refuted the raw HYP-2664 target

    Delta_w^+ <= p1(E')/3

on a targeted sample, but all sampled rows stayed below 3*p1/8.  This script
pushes the bounded B=13 bank without sampling and then dissects the rows above
1/3 into endpoint phase packets.  The purpose is to decide whether the next
constant target is still scalar, or whether the proof must immediately split
by dyadic/phase-packet motifs.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction

from lrc14_far_delta_galois_phase_codex_s38 import (
    NQR,
    QR,
    frac,
    phase_decomposition,
)
from lrc14_p1_tax_envelope_codex_s40 import BaseData, RatioRow, base_data, bounded_bank
from lrc14_residual_plateau_packet_codex_s39 import CAP9


SHELL1 = {1, 2, 4, 8}


@dataclass(frozen=True)
class Packet:
    y: Fraction
    raw: Fraction
    trace: Fraction
    quadratic: Fraction
    residual: Fraction
    counts: tuple[tuple[int, int], ...]

    @property
    def qr_counts(self) -> tuple[int, int, int]:
        data = dict(self.counts)
        return tuple(data.get(s, 0) for s in sorted(QR))

    @property
    def nqr_counts(self) -> tuple[int, int, int]:
        data = dict(self.counts)
        return tuple(data.get(s, 0) for s in sorted(NQR))

    @property
    def abs_count(self) -> int:
        return sum(abs(v) for _s, v in self.counts)


def fmt(q: Fraction) -> str:
    return f"{q} ({float(q):.6f})"


def ap_holes(Ep: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(v for v in range(1, 14) if v not in set(Ep))


def shell1_missing(Ep: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted(SHELL1 - set(Ep)))


def packet_contributions(base: BaseData, w: int) -> list[Packet]:
    grouped: dict[Fraction, Counter[int]] = defaultdict(Counter)
    for x, s, coeff in base.terms:
        grouped[frac(w * x)][s] += coeff

    packets: list[Packet] = []
    for y, counts in grouped.items():
        raw_total = Fraction(0)
        trace_total = Fraction(0)
        quadratic_total = Fraction(0)
        residual_total = Fraction(0)
        for s, coeff in counts.items():
            raw, trace, quadratic, residual = phase_decomposition(y, s)
            raw_total += coeff * raw
            trace_total += coeff * trace
            quadratic_total += coeff * quadratic
            residual_total += coeff * residual
        assert raw_total == trace_total + quadratic_total + residual_total
        packets.append(
            Packet(
                y=y,
                raw=raw_total,
                trace=trace_total,
                quadratic=quadratic_total,
                residual=residual_total,
                counts=tuple(sorted(counts.items())),
            )
        )
    return sorted(packets, key=lambda p: (p.raw, abs(p.raw), p.abs_count), reverse=True)


def tax_gap(row: RatioRow, c: Fraction) -> Fraction:
    """Return c*w*p1 - raw wDelta.  Positive means Delta^+ <= c*p1."""

    return c * row.w * row.p1 - row.raw


def cap_slack(row: RatioRow, c: Fraction) -> Fraction:
    return CAP9 - row.phi - c * row.p1


def print_frontier(rows: list[RatioRow]) -> None:
    sorted_rows = sorted(rows, key=lambda r: (r.actual, r.envelope), reverse=True)
    thresholds = (Fraction(1, 3), Fraction(3, 8), Fraction(2, 5))
    print("full B=13 p1-tax frontier")
    print("family: E'={0}+7-subsets of [1,13], w=max(E')+1..max(E')+8")
    print(f"rows={len(rows)}")
    for c in thresholds:
        count = sum(1 for row in rows if row.actual > c)
        print(f"  rows with Delta^+/p1 > {c}: {count}")
    print(f"  max actual Delta^+/p1 = {fmt(sorted_rows[0].actual)}")
    print(f"  min 2/5 tax gap c*w*p1-raw = {min(tax_gap(r, Fraction(2, 5)) for r in rows)}")
    print(f"  min cap slack for c=2/5 = {min(cap_slack(r, Fraction(2, 5)) for r in rows)}")
    print()
    print("top rows")
    print(
        "rank | ratio | w | E' | AP holes | shell1 missing | p1 | raw wDelta | "
        "tax_gap_2/5 | cap_slack_2/5 | envelope/p1"
    )
    for i, row in enumerate(sorted_rows[:24], 1):
        print(
            f"{i:4d} | {fmt(row.actual):>20} | {row.w:2d} | {row.Ep} | "
            f"{ap_holes(row.Ep)} | {shell1_missing(row.Ep)} | {row.p1} | {row.raw} | "
            f"{tax_gap(row, Fraction(2, 5))} | {cap_slack(row, Fraction(2, 5))} | "
            f"{fmt(row.envelope)}"
        )
    print()


def print_packet_atlas(rows: list[RatioRow], limit: int = 8) -> None:
    sorted_rows = sorted(rows, key=lambda r: (r.actual, r.envelope), reverse=True)
    print("packet anatomy of the top p1-tax rows")
    for i, row in enumerate(sorted_rows[:limit], 1):
        base = base_data(row.Ep)
        packets = packet_contributions(base, row.w)
        positive = sum((p.raw for p in packets if p.raw > 0), Fraction(0))
        negative = -sum((p.raw for p in packets if p.raw < 0), Fraction(0))
        normalizer = row.w * row.p1
        top_pos = [p for p in packets if p.raw > 0][:8]
        pos_l1 = positive / normalizer if normalizer else Fraction(0)
        neg_l1 = negative / normalizer if normalizer else Fraction(0)
        top_share = top_pos[0].raw / positive if positive and top_pos else Fraction(0)
        print(f"row {i}: E'={row.Ep}, w={row.w}")
        print(
            f"  ratio={fmt(row.actual)}, raw={row.raw}, p1={row.p1}, "
            f"holes={ap_holes(row.Ep)}, shell1_missing={shell1_missing(row.Ep)}"
        )
        print(
            f"  positive_packets/(w*p1)={fmt(pos_l1)}, negative_packets/(w*p1)={fmt(neg_l1)}, "
            f"top_positive_share={fmt(top_share)}, packet_count={len(packets)}"
        )
        print(
            f"  excess over 1/3={row.raw - Fraction(1, 3) * normalizer}, "
            f"excess over 3/8={row.raw - Fraction(3, 8) * normalizer}, "
            f"gap below 2/5={tax_gap(row, Fraction(2, 5))}"
        )
        print("  top positive phase packets")
        print("    rank | y | raw/(w*p1) | raw | trace | quad | residual | QR counts | NQR counts | abs_count")
        for j, packet in enumerate(top_pos, 1):
            print(
                f"    {j:4d} | {packet.y} | {fmt(packet.raw / normalizer):>20} | "
                f"{packet.raw} | {packet.trace} | {packet.quadratic} | {packet.residual} | "
                f"{packet.qr_counts} | {packet.nqr_counts} | {packet.abs_count}"
            )
        print()


def main() -> None:
    print("HYP-2667 full B=13 p1-tax packet frontier")
    print("exact Fraction arithmetic")
    print()

    rows = bounded_bank(13, 8, max_rows=None)
    print_frontier(rows)
    print_packet_atlas(rows)

    print("Interpretation")
    print("  S41 refutes the S40 provisional 3p1/8 scalar target on the full B=13 bank.")
    print("  The full bank has 17 rows above 1/3 and 2 rows above 3/8, but no row above 2/5.")
    print("  The two 3/8 failures are low bounded rows with w=12 and strong even/dyadic packet structure.")
    print("  This points to either a 2p1/5 raw tax, or a split theorem:")
    print("    generic phase packets <=3p1/8, dyadic-even packet frontier <=2p1/5.")
    print()
    print("Tournament Analysis")
    print("  vertices: dyadic_even_packet > phase_packet_concentration > p1_tax_constant > interval_envelope > endpoint_count")
    print("  observable: which carrier explains rows with Delta^+/p1 above 3/8")
    print("  switch/gauge: first sort by exact p1 ratio, then expose endpoint packets y=frac(w*x)")
    print("  Hamiltonian path: dyadic_even_packet > phase_packet_concentration > p1_tax_constant > interval_envelope > endpoint_count")
    print("  challenged assumption: after p1/3 failed, 3p1/8 looked viable in the sample; full B=13 shows")
    print("    a small dyadic-even frontier above 3/8, so constants must be stress-tested before canonization.")
    print("PASS: B=13 packet frontier complete.")


if __name__ == "__main__":
    main()
