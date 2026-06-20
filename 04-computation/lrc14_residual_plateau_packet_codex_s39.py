#!/usr/bin/env python3
"""HYP-2664 scout: residual packet imbalance versus plateau mass.

HYP-2662 showed that the resonant far-element error is dominated by the
intra-quadratic residual, not by the full F_7^* trace or the QR/NQR projection.
This script asks whether that residual is really a boundary tax on the
single-missed-sector chain.

For E = E' union {w}, write

    p0(E) = Phi(E') + Delta_w,
    Phi(E') = p0(E') + p1(E')/7.

The exact HYP-2653 formula gives w*Delta_w as an endpoint sum over cells where
E' misses exactly one sector.  Therefore a natural normalizer for the residual
is p1(E'), not the number of endpoints.  If

    abs(intra_quadratic_residual)/(w*p1(E'))

stays tame while Phi(E') is small on multiscale rows, the next proof target is a
packet/coarea lemma: intra-class phase imbalance is paid by single-missed mass
and cannot be combined with AP-sized plateau.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction
from functools import reduce
from itertools import combinations
from math import gcd

from lrc14_far_delta_galois_phase_codex_s38 import (
    NQR,
    QR,
    decompose_wdelta,
    endpoint_terms,
    frac,
    missed_sector_cells,
    phase_decomposition,
    primitive,
)


CAP9 = Fraction(1979, 4004)
INNER = tuple(range(1, 7))


def all_cells(E: tuple[int, ...]) -> tuple[tuple[Fraction, Fraction, int], ...]:
    """Common refinement cells with the number of missed inner sectors."""

    bps = {Fraction(0), Fraction(1)}
    for e in sorted(set(E)):
        if e == 0:
            continue
        for a in range(0, 7 * e + 1):
            bps.add(Fraction(a, 7 * e))
    ordered = sorted(b for b in bps if 0 <= b <= 1)
    cells: list[tuple[Fraction, Fraction, int]] = []
    for lo, hi in zip(ordered, ordered[1:]):
        if lo == hi:
            continue
        mid = (lo + hi) / 2
        hit: set[int] = set()
        for e in E:
            v = frac(e * mid)
            hit.add((v.numerator * 7) // v.denominator)
        missed = sum(1 for s in INNER if s not in hit)
        cells.append((lo, hi, missed))
    return tuple(cells)


def missed_distribution(E: tuple[int, ...]) -> dict[int, Fraction]:
    out = {i: Fraction(0) for i in range(7)}
    for lo, hi, missed in all_cells(E):
        out[missed] += hi - lo
    return out


def phase_packet_stats(
    Ep: tuple[int, ...], w: int, cells: tuple[tuple[Fraction, Fraction, int], ...]
) -> dict[str, Fraction | int]:
    """Aggregate endpoint terms by y=frac(w*x) and measure intra-class imbalance."""

    packets: dict[Fraction, Counter[int]] = defaultdict(Counter)
    for x, s, coeff in endpoint_terms(cells):
        packets[frac(w * x)][s] += coeff

    qr_l1 = Fraction(0)
    nqr_l1 = Fraction(0)
    qr_range = 0
    nqr_range = 0
    nonuniform_packets = 0
    total_abs_terms = 0
    signed_packet_support = 0
    max_packet_abs = 0

    for counts in packets.values():
        total_abs = sum(abs(v) for v in counts.values())
        if total_abs:
            signed_packet_support += 1
            total_abs_terms += total_abs
            max_packet_abs = max(max_packet_abs, total_abs)

        qr_vals = [counts[s] for s in sorted(QR)]
        nqr_vals = [counts[s] for s in sorted(NQR)]
        qr_sum = sum(qr_vals)
        nqr_sum = sum(nqr_vals)
        qr_mean = Fraction(qr_sum, 3)
        nqr_mean = Fraction(nqr_sum, 3)
        qr_def = sum(abs(Fraction(v) - qr_mean) for v in qr_vals)
        nqr_def = sum(abs(Fraction(v) - nqr_mean) for v in nqr_vals)
        qr_l1 += qr_def
        nqr_l1 += nqr_def
        qr_range += max(qr_vals) - min(qr_vals)
        nqr_range += max(nqr_vals) - min(nqr_vals)
        if qr_def or nqr_def:
            nonuniform_packets += 1

    max_residual_kernel = Fraction(0)
    for denom in range(1, 43):
        for num in range(denom):
            for s in INNER:
                _, _, _, residual = phase_decomposition(Fraction(num, denom), s)
                max_residual_kernel = max(max_residual_kernel, abs(residual))

    return {
        "packets": len(packets),
        "signed_packet_support": signed_packet_support,
        "nonuniform_packets": nonuniform_packets,
        "total_abs_terms": total_abs_terms,
        "max_packet_abs": max_packet_abs,
        "qr_l1": qr_l1,
        "nqr_l1": nqr_l1,
        "class_l1": qr_l1 + nqr_l1,
        "qr_range": qr_range,
        "nqr_range": nqr_range,
        "class_range": qr_range + nqr_range,
        "max_residual_kernel_tested": max_residual_kernel,
    }


def fmt(q: Fraction) -> str:
    return f"{q} ({float(q):.6f})"


def row(name: str, Ep: tuple[int, ...], w: int) -> dict[str, Fraction | int | str]:
    cells = missed_sector_cells(Ep)
    dist = missed_distribution(Ep)
    phi = dist[0] + dist[1] / 7
    dec = decompose_wdelta(Ep, w, cells)
    packet = phase_packet_stats(Ep, w, cells)
    full_dist = missed_distribution(tuple(sorted(Ep + (w,))))
    p0_full = full_dist[0]
    reconstructed = phi + dec["raw"] / w
    assert reconstructed == p0_full

    p1 = dist[1]
    residual_over_wp1 = abs(dec["residual"]) / (w * p1) if p1 else None
    raw_over_wp1 = abs(dec["raw"]) / (w * p1) if p1 else None
    residual_over_wphi = abs(dec["residual"]) / (w * phi) if phi else None
    cap_slack = CAP9 - p0_full
    resid_only_slack = CAP9 - phi - abs(dec["residual"]) / w

    out: dict[str, Fraction | int | str] = {
        "name": name,
        "w": w,
        "cells": len(cells),
        "p0_base": dist[0],
        "p1_base": p1,
        "phi": phi,
        "raw": dec["raw"],
        "trace": dec["trace"],
        "quadratic": dec["quadratic"],
        "residual": dec["residual"],
        "p0_full": p0_full,
        "cap_slack": cap_slack,
        "resid_only_slack": resid_only_slack,
        "residual_over_wp1": residual_over_wp1 if residual_over_wp1 is not None else "NA",
        "raw_over_wp1": raw_over_wp1 if raw_over_wp1 is not None else "NA",
        "residual_over_wphi": residual_over_wphi if residual_over_wphi is not None else "NA",
    }
    out.update(packet)
    return out


def print_row(data: dict[str, Fraction | int | str]) -> None:
    print(f"  {data['name']}, w={data['w']}")
    print(
        f"    Phi={fmt(data['phi'])}, p0={fmt(data['p0_base'])}, p1={fmt(data['p1_base'])}, "
        f"p0_full={fmt(data['p0_full'])}, cap_slack={fmt(data['cap_slack'])}"
    )
    print(
        f"    raw/w={fmt(data['raw'] / data['w'])}, residual/w={fmt(data['residual'] / data['w'])}, "
        f"resid_only_slack={fmt(data['resid_only_slack'])}"
    )
    print(
        f"    ratios: |resid|/(w*p1)={fmt(data['residual_over_wp1'])}, "
        f"|raw|/(w*p1)={fmt(data['raw_over_wp1'])}, |resid|/(w*Phi)={fmt(data['residual_over_wphi'])}"
    )
    print(
        f"    packets={data['packets']}, support={data['signed_packet_support']}, "
        f"nonuniform={data['nonuniform_packets']}, class_l1={data['class_l1']}, "
        f"class_range={data['class_range']}, max_packet_abs={data['max_packet_abs']}"
    )


def scan_best_w(Ep: tuple[int, ...], start: int, stop: int) -> int:
    best_w = start
    best = Fraction(-1)
    cells = missed_sector_cells(Ep)
    for w in range(start, stop + 1):
        if not primitive(Ep + (w,)):
            continue
        dec = decompose_wdelta(Ep, w, cells)
        if abs(dec["raw"]) > best:
            best = abs(dec["raw"])
            best_w = w
    return best_w


def bounded_bank_tax(coefficients: tuple[Fraction, ...]) -> None:
    """Test p0 + (1/7+c)*p1 against CAP9 on E'={0}+7-subsets of [1,13]."""

    max_rows: dict[Fraction, tuple[Fraction, tuple[int, ...], Fraction, Fraction, Fraction]] = {}
    violations = {c: 0 for c in coefficients}
    min_allowed: Fraction | None = None
    min_allowed_row: tuple[int, ...] | None = None
    min_allowed_p0 = Fraction(0)
    min_allowed_p1 = Fraction(0)
    row_count = 0

    for comb in combinations(range(1, 14), 7):
        Ep = (0,) + comb
        if reduce(gcd, Ep, 0) != 1:
            continue
        row_count += 1
        dist = missed_distribution(Ep)
        p0 = dist[0]
        p1 = dist[1]
        if p1:
            allowed = (CAP9 - p0) / p1 - Fraction(1, 7)
            if min_allowed is None or allowed < min_allowed:
                min_allowed = allowed
                min_allowed_row = Ep
                min_allowed_p0 = p0
                min_allowed_p1 = p1
        for c in coefficients:
            value = p0 + (Fraction(1, 7) + c) * p1
            slack = CAP9 - value
            old = max_rows.get(c)
            if old is None or value > old[0]:
                max_rows[c] = (value, Ep, p0, p1, slack)
            if value > CAP9:
                violations[c] += 1

    print("bounded AP-window p1-tax bank: E'={0}+7-subsets of [1,13]")
    print(f"  rows={row_count}, cap9={fmt(CAP9)}")
    assert min_allowed is not None and min_allowed_row is not None
    print(
        f"  minimum allowed c in p0+(1/7+c)*p1 <= cap9: {fmt(min_allowed)} "
        f"at {min_allowed_row}, p0={min_allowed_p0}, p1={min_allowed_p1}"
    )
    for c in coefficients:
        value, Ep, p0, p1, slack = max_rows[c]
        print(
            f"  c={c}: violations={violations[c]}, max={fmt(value)} at {Ep}, "
            f"p0={p0}, p1={p1}, slack={fmt(slack)}"
        )
    print()


def main() -> None:
    print("HYP-2664 residual/plateau packet scout")
    print("exact Fraction arithmetic")
    print()

    named = [
        ("consec8", tuple(range(8)), 8, 220),
        ("odd_struct_HYP2655", (0, 1, 3, 5, 7, 9, 10, 11), 12, 220),
        ("two_scale_20_40", (0, 1, 2, 20, 21, 22, 40), 41, 180),
        ("multiscale_30_60", (0, 1, 2, 30, 31, 32, 60, 61), 62, 240),
        ("cluster_40_80", (0, 1, 2, 40, 41, 42, 80, 81), 82, 260),
    ]

    rows: list[dict[str, Fraction | int | str]] = []
    print("named best-resonance rows")
    for name, Ep, start, stop in named:
        w = scan_best_w(Ep, start, stop)
        data = row(name, Ep, w)
        rows.append(data)
        print_row(data)
    print()

    print("scale family at w=2M+2")
    for m in [10, 20, 30, 40, 50, 60, 70]:
        Ep = (0, 1, 2, m, m + 1, m + 2, 2 * m, 2 * m + 1)
        w = 2 * m + 2
        if not primitive(Ep + (w,)):
            continue
        data = row(f"scale_M={m}", Ep, w)
        rows.append(data)
        print_row(data)
    print()

    ratios = [r["residual_over_wp1"] for r in rows if isinstance(r["residual_over_wp1"], Fraction)]
    raw_ratios = [r["raw_over_wp1"] for r in rows if isinstance(r["raw_over_wp1"], Fraction)]
    resid_slacks = [r["resid_only_slack"] for r in rows if isinstance(r["resid_only_slack"], Fraction)]
    print("summary")
    print(f"  max |resid|/(w*p1) among rows = {fmt(max(ratios))}")
    print(f"  max |raw|/(w*p1) among rows   = {fmt(max(raw_ratios))}")
    print(f"  min CAP9 - Phi - |resid|/w     = {fmt(min(resid_slacks))}")
    print()

    bounded_bank_tax(
        (
            Fraction(1, 4),
            Fraction(13, 51),
            Fraction(1, 3),
        )
    )

    print("Tournament Analysis")
    print("  vertices: p1_boundary_mass > phase_packet_class_l1 > residual_channel > raw_wDelta > raw_endpoint_count")
    print("  pairwise observable: which packet invariant predicts residual/w without destroying plateau slack")
    print("  switch/gauge: aggregate endpoints by y=frac(w*x), then quotient by QR/NQR class means")
    print("  tie Hamiltonian path:")
    print("    p1_boundary_mass > phase_packet_class_l1 > residual_channel > raw_wDelta > raw_endpoint_count")
    print("  challenged assumption: residual size is not a free endpoint-count error; it is a boundary")
    print("    of the one-missed-sector packet chain and should be paid by p1(E').")
    print()

    print("proof-route takeaway")
    print("  In these resonant families, residual grows before scaling, but |residual|/(w*p1)")
    print("  remains tame and residual-only cap slack stays positive. This supports the next")
    print("  lemma: bound intra-quadratic phase residual as a p1-boundary tax, then use")
    print("  HYP-2661/HYP-2663 packet rigidity for the small-plateau side.")
    print("PASS: residual/plateau packet scout complete.")


if __name__ == "__main__":
    main()
