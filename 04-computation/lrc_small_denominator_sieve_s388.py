#!/usr/bin/env python3
"""
lrc_small_denominator_sieve_s388.py

codex-2026-05-31 S388

Formal/computational support for the small-denominator divisibility sieve in
the Lonely Runner endpoint program.

The theorem-level fact is elementary: at level n, a primitive point a/m with
m <= n can be strictly protected only by speeds divisible by m.  This script
uses that filter as a microscope:

* verify the predicted small-denominator lonely witnesses for sampled sets;
* compare naive n-gate swaps with lcm(n,d) sieve-completing swaps;
* record the endpoint-core and endpoint-debt residue left after the full
  small-denominator sieve has been satisfied.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from math import gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
S356 = SourceFileLoader(
    "lonely_runner_residue_probe_s356",
    str(ROOT / "04-computation" / "lonely_runner_residue_probe_s356.py"),
).load_module()
S362 = SourceFileLoader(
    "lonely_runner_bohr_descent_s362",
    str(ROOT / "04-computation" / "lonely_runner_bohr_descent_s362.py"),
).load_module()

ONE = Fraction(1, 1)


@dataclass(frozen=True)
class SieveWitness:
    modulus: int
    residue: int
    point: Fraction
    open_margin: Fraction
    kind: str


@dataclass(frozen=True)
class AuditRow:
    label: str
    n: int
    speeds: tuple[int, ...]
    missing_moduli: tuple[int, ...]
    first_witness: SieveWitness | None
    forbidden_length: Fraction
    max_gap_ratio: Fraction
    boundary_witnesses: int
    components: int
    unprotected: int
    core_endpoints: int
    first_layer: int
    first_layer_modulus: int
    layer_hist: tuple[tuple[int, int], ...]


@dataclass(frozen=True)
class FastRow:
    label: str
    n: int
    speeds: tuple[int, ...]
    missing_moduli: tuple[int, ...]
    first_witness: SieveWitness | None
    forbidden_length: Fraction
    max_gap_ratio: Fraction
    boundary_witnesses: int
    components: int


def fmt(x: Fraction | None) -> str:
    return S356.fmt_frac(x)


def circle(x: Fraction) -> Fraction:
    return x % ONE


def dist_to_integer(x: Fraction) -> Fraction:
    y = circle(x)
    return min(y, ONE - y)


def normalized(raw_speeds: list[int]) -> tuple[int, ...]:
    return S356.normalize_speed_set(raw_speeds)


def missing_moduli(n: int, speeds: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(
        m for m in range(2, n + 1) if all(v % m != 0 for v in speeds)
    )


def primitive_residue(m: int) -> int:
    # Use the smallest unit.  Keeping this deterministic makes the output
    # suitable for session-to-session diffs.
    for a in range(1, m):
        if gcd(a, m) == 1:
            return a
    return 1


def small_denominator_witness(
    n: int, speeds: tuple[int, ...]
) -> SieveWitness | None:
    for m in missing_moduli(n, speeds):
        a = primitive_residue(m)
        t = Fraction(a, m)
        margin = min(dist_to_integer(v * t) for v in speeds) - Fraction(1, n)
        kind = "open" if m < n else "boundary"
        return SieveWitness(m, a, t, margin, kind)
    return None


def endpoint_layer_hist(values: set[Fraction]) -> tuple[tuple[int, int], ...]:
    hist: dict[int, int] = {}
    for value in values:
        hist[value.denominator] = hist.get(value.denominator, 0) + 1
    return tuple(sorted(hist.items()))


def audit(label: str, raw_speeds: list[int]) -> AuditRow:
    report = S356.report(label, raw_speeds)
    speeds = report.speeds
    n = len(speeds) + 1
    endpoints, intervals, _owners, protectors, boundary = S362.build_endpoint_system(
        speeds
    )
    unprotected = {endpoint for endpoint in endpoints if not protectors[endpoint]}
    layers, core_endpoints, _core_intervals = S362.peel_protection_core(
        report.boundary_modulus, endpoints, intervals, protectors, boundary
    )
    first_layer_values = unprotected
    first_layer_modulus = S362.subgroup_modulus(report.boundary_modulus, first_layer_values)
    first = small_denominator_witness(n, speeds)
    return AuditRow(
        label=label,
        n=n,
        speeds=speeds,
        missing_moduli=missing_moduli(n, speeds),
        first_witness=first,
        forbidden_length=report.forbidden_length,
        max_gap_ratio=report.max_gap / report.threshold,
        boundary_witnesses=report.boundary_witness_count,
        components=report.components,
        unprotected=len(unprotected),
        core_endpoints=len(core_endpoints),
        first_layer=layers[0].removed_endpoints if layers else 0,
        first_layer_modulus=first_layer_modulus,
        layer_hist=endpoint_layer_hist(unprotected),
    )


def fast_audit(label: str, raw_speeds: list[int]) -> FastRow:
    """Gap/sieve audit without the quadratic endpoint-incidence pass."""

    report = S356.report(label, raw_speeds)
    speeds = report.speeds
    n = len(speeds) + 1
    return FastRow(
        label=label,
        n=n,
        speeds=speeds,
        missing_moduli=missing_moduli(n, speeds),
        first_witness=small_denominator_witness(n, speeds),
        forbidden_length=report.forbidden_length,
        max_gap_ratio=report.max_gap / report.threshold,
        boundary_witnesses=report.boundary_witness_count,
        components=report.components,
    )


def print_row(row: AuditRow) -> None:
    first = row.first_witness
    if first is None:
        witness = "sieve-complete"
    else:
        witness = (
            f"{first.kind} t={fmt(first.point)} "
            f"(missing m={first.modulus}, margin={fmt(first.open_margin)})"
        )
    hist = " ".join(f"{d}:{c}" for d, c in row.layer_hist[:6])
    if len(row.layer_hist) > 6:
        hist += " ..."
    print(f"[{row.label}]")
    print(f"  speeds={row.speeds}")
    print(
        f"  missing={row.missing_moduli if row.missing_moduli else '-'}  "
        f"{witness}"
    )
    print(
        f"  length={fmt(row.forbidden_length)}  "
        f"gap/th={float(row.max_gap_ratio):.6f}  "
        f"boundary={row.boundary_witnesses}  components={row.components}"
    )
    print(
        f"  unprotected={row.unprotected}  coreE={row.core_endpoints}  "
        f"first_peel={row.first_layer}  first_layer_group={row.first_layer_modulus}"
    )
    print(f"  unprotected_denoms={hist or '-'}")
    print(flush=True)


def initial(n: int) -> list[int]:
    return list(range(1, n))


def drop_add(n: int, d: int, replacement: int) -> list[int]:
    speeds = set(initial(n))
    speeds.remove(d)
    speeds.add(replacement)
    return sorted(speeds)


def section(title: str) -> None:
    print("=" * 88)
    print(title)
    print("=" * 88)


def theorem_audit() -> None:
    section("SMALL-DENOMINATOR SIEVE WITNESS AUDIT")
    samples = [
        ("initial n=8", initial(8)),
        ("n8 sporadic tight", [1, 4, 5, 6, 7, 11, 13]),
        ("initial n=14", initial(14)),
        ("drop 13 add 14", drop_add(14, 13, 14)),
        ("drop 13 add 182", drop_add(14, 13, lcm(14, 13))),
        ("n14 seven-ladder", [1, 7, 14, 21, 28, 35, 49, 56, 63, 70, 77, 84, 91]),
        ("initial n=15", initial(15)),
        ("drop 14 add 15", drop_add(15, 14, 15)),
        ("drop 14 add 210", drop_add(15, 14, lcm(15, 14))),
        ("initial n=18", initial(18)),
        ("drop 17 add 18", drop_add(18, 17, 18)),
        ("drop 17 add 306", drop_add(18, 17, lcm(18, 17))),
    ]
    for label, speeds in samples:
        print_row(audit(label, speeds))


def swap_table(n: int) -> None:
    section(f"INITIAL-SEGMENT GATE VS LCM-SIEVE TRANSFERS AT n={n}")
    print("drop | add n: missing gap/th boundary | add lcm(n,drop): missing gap/th boundary")
    print("-" * 88)
    for d in range(2, n):
        naive = fast_audit(f"n{n} drop {d} add {n}", drop_add(n, d, n))
        repaired = fast_audit(
            f"n{n} drop {d} add lcm", drop_add(n, d, lcm(n, d))
        )
        naive_missing = ",".join(str(m) for m in naive.missing_moduli) or "-"
        repaired_missing = ",".join(str(m) for m in repaired.missing_moduli) or "-"
        print(
            f"{d:>4} | "
            f"{naive_missing:<9} {float(naive.max_gap_ratio):>7.4f} "
            f"{naive.boundary_witnesses:>8} | "
            f"{repaired_missing:<9} {float(repaired.max_gap_ratio):>7.4f} "
            f"{repaired.boundary_witnesses:>8}",
            flush=True,
        )
    print(flush=True)


def main() -> None:
    print("LRC small-denominator sieve and debt audit (codex-2026-05-31 S388)", flush=True)
    print(
        "At level n, primitive point a/m with m<=n can only be strictly "
        "protected by speeds divisible by m.",
        flush=True,
    )
    print(flush=True)
    theorem_audit()
    for n in (14, 15, 18):
        swap_table(n)
    section("SYNTHESIS")
    print(
        "The lcm repair rows satisfy the full small-denominator sieve, but in "
        "the n=14,15,18 transfer tables they still remain positive-gap.  The "
        "focused endpoint audits of the largest lcm repairs are also "
        "terminal-core empty."
    )
    print(
        "This separates two obstructions: missing a small denominator gives an "
        "immediate witness; covering it by a larger lcm gate exports the "
        "problem into endpoint debt rather than producing a cover."
    )


if __name__ == "__main__":
    main()
