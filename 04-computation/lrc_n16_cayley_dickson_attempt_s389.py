#!/usr/bin/env python3
"""
lrc_n16_cayley_dickson_attempt_s389.py

codex-2026-05-31 S389

A proof attempt for the Lonely Runner denominator n=16, guided by the
Cayley-Dickson analogy.

The analogy used here is intentionally concrete:

* Cayley-Dickson doubling keeps a norm-like ledger through the octonions and
  then exposes zero-divisor behavior at dimension 16.
* The LRC n=16 gate protects the old unit endpoints, but the exact endpoint
  ledger shows the protection is paid for by descendant dyadic endpoints.

So the proof attempt is a three-way split:

1. no 16-gate: odd unit endpoints a/16 are unprotected;
2. dyadic gate/lpd ladders: visible gaps halve while endpoint debt doubles;
3. scalar quotient: non-scalar dyadic defects have a finite cell moat.

This is not claimed as a proof.  The output records the missing lemma precisely.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
S356 = SourceFileLoader(
    "lonely_runner_residue_probe_s356",
    str(ROOT / "04-computation" / "lonely_runner_residue_probe_s356.py"),
).load_module()
S360 = SourceFileLoader(
    "lonely_runner_endpoint_protection_s360",
    str(ROOT / "04-computation" / "lonely_runner_endpoint_protection_s360.py"),
).load_module()
S372 = SourceFileLoader(
    "lonely_runner_creative_multiroute_s372",
    str(ROOT / "04-computation" / "lonely_runner_creative_multiroute_s372.py"),
).load_module()


N = 16
K = N - 1


@dataclass(frozen=True)
class FamilyRow:
    label: str
    speeds: tuple[int, ...]
    classification: str
    gap_ratio: Fraction
    max_gap: Fraction
    unprotected: int
    first_unprotected: Fraction | None
    peel_depth: int
    core_endpoints: int
    mult16: int
    endpoint_q: int
    unprotected_layer_hist: tuple[tuple[tuple[int, int], int], ...]
    owner_v2_hist: tuple[tuple[int, int], ...]

    @property
    def debt_norm(self) -> Fraction:
        return self.gap_ratio * self.unprotected


def fmt_frac(value: Fraction | None) -> str:
    return S356.fmt_frac(value)


def fmt_decimal(value: Fraction | None) -> str:
    if value is None:
        return "-"
    return f"{float(value):.6f}"


def v2(value: int) -> int:
    if value == 0:
        return 99
    count = 0
    while value % 2 == 0:
        count += 1
        value //= 2
    return count


def primitive(speeds: tuple[int, ...]) -> bool:
    g = 0
    for speed in speeds:
        g = gcd(g, speed)
    return g == 1


def dyadic_ladder(d: int, skip: int = 14) -> tuple[int, ...]:
    speeds = tuple(sorted({1} | {d * q for q in range(1, N) if q != skip}))
    if len(speeds) != K or not primitive(speeds):
        raise ValueError(f"bad dyadic ladder d={d}, skip={skip}: {speeds}")
    return speeds


def best_dyadic_ladder(d: int) -> tuple[int, tuple[int, ...]]:
    best: tuple[Fraction, int, int, tuple[int, ...]] | None = None
    for skip in range(1, N):
        speeds = dyadic_ladder(d, skip)
        summary = S360.summarize(list(speeds))
        gap_ratio = summary.max_gap / summary.threshold
        key = (gap_ratio, summary.unprotected_count, skip, speeds)
        if best is None or key < best:
            best = key
    assert best is not None
    return best[2], best[3]


def endpoint_layers(speeds: tuple[int, ...]) -> tuple[
    tuple[tuple[tuple[int, int], int], ...],
    tuple[tuple[int, int], ...],
]:
    points_by_value: dict[Fraction, list[object]] = defaultdict(list)
    for endpoint in S360.endpoints(speeds):
        points_by_value[endpoint.value].append(endpoint)

    layer_hist: Counter[tuple[int, int]] = Counter()
    owner_v2_hist: Counter[int] = Counter()
    for value, labels in points_by_value.items():
        protectors = [
            speed for speed in speeds if S360.direct_protects(speeds, speed, value)
        ]
        if protectors:
            continue
        denominator = value.denominator
        layer = (v2(denominator), denominator >> v2(denominator))
        layer_hist[layer] += 1
        for label in labels:
            owner_v2_hist[v2(label.speed)] += 1

    return (
        tuple(sorted(layer_hist.items(), key=lambda item: (item[0][0], item[0][1]))),
        tuple(sorted(owner_v2_hist.items())),
    )


def summarize_family(label: str, speeds: tuple[int, ...]) -> FamilyRow:
    endpoint_summary = S360.summarize(list(speeds))
    descent = S372.S362.summarize(list(endpoint_summary.speeds))
    layer_hist, owner_v2_hist = endpoint_layers(tuple(endpoint_summary.speeds))
    return FamilyRow(
        label=label,
        speeds=tuple(endpoint_summary.speeds),
        classification=endpoint_summary.classification,
        gap_ratio=endpoint_summary.max_gap / endpoint_summary.threshold,
        max_gap=endpoint_summary.max_gap,
        unprotected=endpoint_summary.unprotected_count,
        first_unprotected=endpoint_summary.first_unprotected,
        peel_depth=len(descent.peel_layers),
        core_endpoints=descent.core_endpoint_count,
        mult16=sum(1 for speed in endpoint_summary.speeds if speed % N == 0),
        endpoint_q=endpoint_summary.boundary_modulus,
        unprotected_layer_hist=layer_hist,
        owner_v2_hist=owner_v2_hist,
    )


def family_rows() -> tuple[FamilyRow, ...]:
    raw: list[tuple[str, tuple[int, ...]]] = [
        ("initial 1..15", tuple(range(1, N))),
        ("single 16 gate", tuple(list(range(1, N - 1)) + [N])),
        ("drop 8 add 16", (1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 16)),
        ("drop 14 add 16", (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 16)),
    ]
    for d in (2, 4, 8, 16):
        skip, speeds = best_dyadic_ladder(d)
        raw.append((f"d={d} ladder skip {skip}", speeds))
    return tuple(summarize_family(label, speeds) for label, speeds in raw)


def gate_sample_rows() -> tuple[FamilyRow, ...]:
    samples = [
        (
            "fixed gate sample A",
            (13, 25, 48, 61, 63, 72, 75, 90, 94, 110, 118, 125, 126, 127, 128),
        ),
        (
            "fixed gate sample B",
            (9, 34, 47, 48, 57, 81, 85, 90, 95, 101, 103, 113, 117, 120, 128),
        ),
        (
            "fixed gate sample C",
            (2, 16, 23, 25, 47, 48, 55, 57, 67, 97, 99, 101, 105, 111, 112),
        ),
    ]
    return tuple(summarize_family(label, speeds) for label, speeds in samples)


def unit_gate_lemma_audit() -> None:
    print("UNIT-GATE LEMMA")
    print("=" * 78)
    print("For odd a, the unit endpoint t=a/16 can be strictly protected by")
    print("speed v only if v*a == 0 mod 16, hence only if 16 divides v.")
    print("This gives a real first branch of a proof attempt:")
    print("  no 16-multiple -> an odd a/16 endpoint is a boundary witness.")
    print()
    for v in range(1, 33):
        protected = [
            a
            for a in range(1, N, 2)
            if S360.circular_distance_to_integer(Fraction(v * a, N)) < Fraction(1, N)
        ]
        if protected:
            print(f"  v={v:2d} protects odd unit endpoints {protected}")
    print()


def print_family_table(rows: tuple[FamilyRow, ...]) -> None:
    print("DYADIC GATE AND DESCENDANT-DEBT AUDIT")
    print("=" * 78)
    print(
        "family                 mult16 class         gap/th      unprot "
        "debt_norm  peel core first"
    )
    print("-" * 92)
    for row in rows:
        print(
            f"{row.label:22s} {row.mult16:6d} {row.classification:12s} "
            f"{fmt_decimal(row.gap_ratio):>10s} {row.unprotected:7d} "
            f"{fmt_frac(row.debt_norm):>9s} {row.peel_depth:5d} "
            f"{row.core_endpoints:4d} {fmt_frac(row.first_unprotected):>8s}"
        )
    print()

    print("Unprotected endpoint layers for the dyadic ladders")
    print("layer=(v2(denominator), odd_part(denominator))")
    for row in rows:
        if "ladder" not in row.label:
            continue
        layer_text = ", ".join(
            f"{layer}:{count}" for layer, count in row.unprotected_layer_hist
        )
        owners = ", ".join(f"v2={layer}:{count}" for layer, count in row.owner_v2_hist)
        print(f"  {row.label:22s} layers {{{layer_text}}}; owners {{{owners}}}")
    print()


def print_norm_ledger(rows: tuple[FamilyRow, ...]) -> None:
    ladder_rows = [row for row in rows if "ladder" in row.label]
    print("CAYLEY-DICKSON-STYLE NORM LEDGER")
    print("=" * 78)
    print("Each doubling pushes the first unprotected endpoints one dyadic level")
    print("deeper.  The visible gap halves, but endpoint debt doubles, so the")
    print("product unprotected*(gap/th) stays essentially fixed.")
    print()
    for row in ladder_rows:
        print(
            f"  {row.label:22s} gap/th={fmt_frac(row.gap_ratio):>7s} "
            f"unprot={row.unprotected:3d} norm={fmt_frac(row.debt_norm)}"
        )
    print()
    print("Interpretation:")
    print("  The 16-gate is not closing a cover; it is acting like a zero-divisor")
    print("  repair.  It kills the old unit witnesses and creates a kernel one")
    print("  dyadic level lower.  A proof should preserve this debt norm until an")
    print("  endpoint leaf appears.")
    print()


def print_scalar_moat() -> None:
    print("NORMALIZED SCALAR-QUOTIENT MOAT AT n=16")
    print("=" * 78)
    system = S372.build_pattern_system(N)
    scalar_full, scalar_hist = S372.scalar_audit(system)
    print(
        f"cell_system: patterns={len(system.patterns)} "
        f"candidates={system.candidate_count} scalar_full={scalar_full} "
        f"scalar_gcd_hist={scalar_hist}"
    )
    for radius in (1, 2):
        summary = S372.puncture_summary(system, radius)
        print(
            f"radius={radius}: checked={summary.checked} "
            f"best_missed={summary.best_missed} best_count={summary.best_count} "
            f"positions={summary.positions} deltas={summary.deltas}"
        )
        print(f"  first_example={summary.examples[0]}")
    print()
    print("Interpretation:")
    print("  The best one-defect is a half-turn in the last coordinate, missing")
    print("  128 cells.  The best two-defect still misses 160 cells.  This is the")
    print("  finite scalar moat a full n=16 proof would need to lift from the")
    print("  normalized residue model back to arbitrary primitive speed sets.")
    print()


def print_gate_samples(rows: tuple[FamilyRow, ...]) -> None:
    print("FIXED GATE SAMPLE AUDIT")
    print("=" * 78)
    print("These are deterministic high-speed 16-gate samples retained from a")
    print("short random pressure pass; they are exact interval/endpoint audits.")
    for row in rows:
        print(
            f"  {row.label}: class={row.classification} mult16={row.mult16} "
            f"gap/th={fmt_decimal(row.gap_ratio)} unprot={row.unprotected} "
            f"peel={row.peel_depth} core_E={row.core_endpoints}"
        )
    print()


def print_attempt_summary() -> None:
    print("PROOF ATTEMPT STATUS")
    print("=" * 78)
    print("A plausible n=16 proof now has this shape:")
    print("  1. No 16-gate is impossible for a counterexample by the unit-gate lemma.")
    print("  2. A 16-gate closes the old unit skeleton only by exporting endpoint")
    print("     debt to v2-denominator layers 5,6,7,8.  The dyadic ladders preserve")
    print("     a positive debt norm instead of approaching an open cover.")
    print("  3. In the normalized scalar quotient, every one- or two-coordinate")
    print("     non-scalar dyadic defect has a large moat.")
    print()
    print("Missing lemma:")
    print("  Every primitive n=16 all-protected endpoint system with a 16-gate")
    print("  either descends to one of the dyadic debt ledgers above, or contains")
    print("  a labelled endpoint cycle with a private dyadic leaf.  Proving that")
    print("  would turn the Cayley-Dickson analogy into an actual n=16 proof.")


def main() -> None:
    print("LRC n=16 Cayley-Dickson proof attempt (codex-2026-05-31 S389)")
    print("n=16 means k=15 moving speeds and threshold 1/16.\n")
    print("Cayley-Dickson dictionary:")
    print("  doubling layer      -> dyadic quotient layer")
    print("  norm conservation   -> gap/debt ledger")
    print("  dimension-16 zero divisor -> 16-gate protecting old units but leaking")
    print("  lost alternativity  -> labelled endpoint cycle must be checked, not")
    print("                         just the projected support cycle")
    print()

    unit_gate_lemma_audit()
    rows = family_rows()
    print_family_table(rows)
    print_norm_ledger(rows)
    print_scalar_moat()
    print_gate_samples(gate_sample_rows())
    print_attempt_summary()


if __name__ == "__main__":
    main()
