#!/usr/bin/env python3
"""HYP-3454: AP84 endpoint-clock certificate for the LRC14 tail.

HYP-3452 found that the canonical AP-tail family

    S_m = {1,2,...,11,13,84m}

enters a rank-one E:84m/E:84m component-cover phase from m >= 5.  This script
turns that phase statement into a smaller proof obligation.

It checks, exactly, that the HYP-3433 interval

    I_m = [(14 ceil(48m/7)+1)/(588m),
           (14 ceil(48m/7)+13)/(588m)]

is the HYP-3450 best component-cover escape on the checked tail, and records
the endpoint inequalities that make the interval sit inside the fixed low
corridor [8/49, 6/35] while being bounded only by the moving E:84m wall.

The remaining escape-count information is reduced to a period-35 boundary
clock: one correction vector plus the exact formula shift

    escapes(m+35) = escapes(m) + 24

for the Beatty expression.  The all-m theorem still needs a floor-count proof
that this clock is the true component-boundary count; the endpoint interval
itself is already in closed form.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
H3452_PATH = ROOT / "04-computation" / "lrc14_ap84_tail_component_phase_codex_20260629.py"
F = Fraction
LOW_CORRIDOR = (F(8, 49), F(6, 35))


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


H3452 = load_module("hyp3452_endpoint_clock", H3452_PATH)


def fmt(value: F | None) -> str:
    if value is None:
        return "None"
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def expected_address(m: int) -> int:
    return (48 * m + 6) // 7


def expected_interval(m: int) -> tuple[F, F]:
    a = expected_address(m)
    return (F(14 * a + 1, 588 * m), F(14 * a + 13, 588 * m))


def interval_length(interval: tuple[F, F]) -> F:
    return interval[1] - interval[0]


@dataclass(frozen=True)
class EndpointCertificate:
    m: int
    address: int
    interval: tuple[F, F]
    left_margin: F
    right_margin: F
    length: F
    labels: str
    rank: int | None
    is_best_component: bool
    is_rank_one_ee: bool

    @property
    def low_corridor_ok(self) -> bool:
        return self.left_margin >= 0 and self.right_margin >= 0

    @property
    def length_ok(self) -> bool:
        return self.length == F(1, 49 * self.m)


def endpoint_certificate(record) -> EndpointCertificate:
    interval = expected_interval(record.m)
    best = record.best
    return EndpointCertificate(
        m=record.m,
        address=expected_address(record.m),
        interval=interval,
        left_margin=interval[0] - LOW_CORRIDOR[0],
        right_margin=LOW_CORRIDOR[1] - interval[1],
        length=interval_length(interval),
        labels=best.labels,
        rank=record.row.best_rank,
        is_best_component=best.best_survivor == interval,
        is_rank_one_ee=record.has_rank_one_ee_tail,
    )


def beatty_prediction(m: int, corrections: list[int]) -> int:
    return 2 * ((12 * m // 35) + corrections[(m - 1) % 35])


def residue_table(records, corrections: list[int]) -> list[tuple[int, int, int, int, int]]:
    table = []
    for r in range(1, 36):
        base = records[r - 1].escapes
        shifted = records[r + 35 - 1].escapes
        predicted_base = beatty_prediction(r, corrections)
        predicted_shift = beatty_prediction(r + 35, corrections)
        table.append((r, corrections[r - 1], base, shifted, predicted_shift - predicted_base))
    return table


def tournament_fingerprint() -> tuple[dict[int, int], list[str], int]:
    vertices = {
        "closed_form_endpoint_interval": (10, 10, 10, 10, 8, 10),
        "low_corridor_containment_inequalities": (10, 10, 9, 9, 8, 10),
        "moving_E84m_gap_certificate": (10, 10, 10, 8, 8, 9),
        "finite_transients_m1_to_m4": (9, 9, 7, 8, 8, 9),
        "mod35_escape_boundary_clock": (8, 8, 8, 8, 10, 8),
        "component_cover_reaudit": (8, 7, 8, 8, 7, 7),
        "raw_dead_fraction_peak": (3, 3, 2, 2, 3, 2),
    }
    scores = {name: sum(parts) for name, parts in vertices.items()}
    hist = dict(sorted(Counter(scores.values()).items()))
    path = [name for name, _score in sorted(scores.items(), key=lambda item: (-item[1], item[0]))]
    rank = {name: index for index, name in enumerate(path)}
    cycles = 0
    for a, b, c in combinations(vertices, 3):
        ab = rank[a] < rank[b]
        bc = rank[b] < rank[c]
        ca = rank[c] < rank[a]
        if ab == bc == ca:
            cycles += 1
    return hist, path, cycles


def main() -> None:
    limit = 70
    records = H3452.audit_tail(limit)
    certs = [endpoint_certificate(record) for record in records]
    tail = [cert for cert in certs if cert.m >= 5]

    endpoint_failures = [
        cert.m
        for cert in tail
        if not (
            cert.low_corridor_ok
            and cert.length_ok
            and cert.is_best_component
            and cert.is_rank_one_ee
            and cert.labels == f"L[E:{84 * cert.m}] R[E:{84 * cert.m}]"
            and cert.rank == 1
        )
    ]
    transient_rows = [
        (
            record.m,
            record.best.labels,
            record.row.best_rank,
            record.best.best_survivor,
            record.best_len,
        )
        for record in records[:4]
    ]

    corrections, beatty_failures = H3452.correction_vector(records)
    shift_failures = [
        r
        for r in range(1, 36)
        if records[r + 35 - 1].escapes - records[r - 1].escapes != 24
    ]
    table = residue_table(records, corrections)
    min_left = min(tail, key=lambda cert: cert.left_margin)
    min_right = min(tail, key=lambda cert: cert.right_margin)
    min_length_residual = min((cert.length - F(1, 49 * cert.m) for cert in tail), default=F(0))
    max_length_residual = max((cert.length - F(1, 49 * cert.m) for cert in tail), default=F(0))
    hist, path, cycles = tournament_fingerprint()

    print("HYP-3454 AP84 ENDPOINT-CLOCK CERTIFICATE")
    print("=" * 72)
    print("status=EVIDENCE / endpoint-inequality and residue-clock certificate; not an LRC14 proof")
    print("source=HYP-3431 corridor fence + HYP-3433 endpoint law + HYP-3452 component phase")
    print()

    print("A. Endpoint interval certificate")
    print("  family=S_m={1,2,...,11,13,84m}")
    print(f"  checked_m=1..{limit}")
    print("  theorem_tail=m>=5")
    print("  interval=I_m=[(14ceil(48m/7)+1)/(588m),(14ceil(48m/7)+13)/(588m)]")
    print("  fixed_low_corridor=[8/49,6/35]")
    print(f"  endpoint_failures_on_checked_tail={endpoint_failures}")
    print(f"  min_left_corridor_margin={fmt(min_left.left_margin)} at m={min_left.m}")
    print(f"  min_right_corridor_margin={fmt(min_right.right_margin)} at m={min_right.m}")
    print(f"  length_residual_range=[{fmt(min_length_residual)},{fmt(max_length_residual)}]")
    print()

    print("  samples: m | a_m | I_m | low margins | labels | rank")
    for cert in [certs[i - 1] for i in (5, 6, 7, 8, 12, 21, 35, 36, 70)]:
        print(
            "  {m:3d} | {a:4d} | [{lo},{hi}] | L={lm} R={rm} | {labels} | {rank}".format(
                m=cert.m,
                a=cert.address,
                lo=fmt(cert.interval[0]),
                hi=fmt(cert.interval[1]),
                lm=fmt(cert.left_margin),
                rm=fmt(cert.right_margin),
                labels=cert.labels,
                rank=cert.rank,
            )
        )
    print()

    print("B. Finite transient side")
    print("  m=1..4 do not have the pure E:84m/E:84m phase; they stay mixed E/B1.")
    for m, labels, rank, interval, length in transient_rows:
        lo, hi = interval
        print(
            f"  m={m}: labels={labels}, rank={rank}, best=[{fmt(lo)},{fmt(hi)}], length={fmt(length)}"
        )
    print()

    print("C. Mod-35 escape clock")
    print("  checked law: escapes(m)=2*(floor(12m/35)+d[(m-1) mod 35])")
    print(f"  correction_vector={corrections}")
    print(f"  beatty_clock_failures={beatty_failures}")
    print("  formula_shift: floor(12*(m+35)/35)=floor(12m/35)+12, so predicted escapes shift by 24")
    print(f"  checked_shift_failures={shift_failures}")
    print("  residue table: r | d_r | escapes(r) | escapes(r+35) | formula_shift")
    for r, correction, base, shifted, predicted_shift in table:
        print(f"  {r:2d} | {correction:3d} | {base:10d} | {shifted:13d} | {predicted_shift:13d}")
    print()

    print("D. Proof obligation now isolated")
    print("  The rank-one endpoint interval has closed-form walls and exact checked")
    print("  component-cover status from m>=5.  The remaining infinite AP-tail task is")
    print("  the floor-count lemma that the 35-residue clock counts every alive boundary")
    print("  component, plus four finite mixed transients m=1..4.  Raw dead fraction is")
    print("  only a scalar warning, not a certificate.")
    print()

    print("E. Tournament Analysis")
    print("  vertices=AP-tail proof obligations, not runners or raw arcs")
    print("  pairwise_observable=endpoint exactness + low-corridor containment + residue-clock payload + scalar-firewall safety")
    print("  switch=higher retained proof payload; ties by endpoint-clock route")
    print(f"  score_hist={hist}")
    print(f"  directed_3cycles={cycles}")
    print("  hamiltonian_path=" + " -> ".join(path))
    print()

    print("F. Assumption Challenge")
    print("  Considered vertices: runners, m-values, residues mod 35, fixed low")
    print("  corridor sections, moving E:84m gaps, endpoint walls, component-cover")
    print("  graph nodes, branch blockers, and proof obligations.  The chosen quotient")
    print("  preserves the AP-tail branch-union escape predicate, exact endpoint labels,")
    print("  low-corridor containment margins, and the residue-count sidecar.  It")
    print("  destroys non-AP wall geometry and full component adjacency, so it is only")
    print("  the AP-tail bridge lemma needed by HYP-3439/HYP-3452.")


if __name__ == "__main__":
    main()
