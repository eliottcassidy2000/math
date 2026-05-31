#!/usr/bin/env python3
"""
lrc_above14_runner_favorites_s383.py

codex-2026-05-31 S383

Creative triage for Lonely Runner denominators above n=14.

The goal is not to search for a counterexample directly.  It is to rank the
first denominators above 14 by the structural signals this repo has been using:

* unit-skeleton density phi(n)/(n-1);
* divisor/quotient richness;
* exact largest-proper-divisor ladder pressure;
* product-sum factor-packing type;
* tractability as the next proof/disproof target.

The script prints several rankings, then makes one explicit mathematical pick.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from math import gcd, isqrt
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
OPS = SourceFileLoader(
    "natural_operation_modes_s366",
    str(ROOT / "04-computation" / "natural_operation_modes_s366.py"),
).load_module()


N_MIN = 15
N_MAX = 24
FAVORITE_N = 18


@dataclass(frozen=True)
class LadderRow:
    n: int
    factors: str
    phi: int
    tau: int
    unit_density: Fraction
    lpd: int | None
    skip: int | None
    classification: str
    gap_ratio: Fraction | None
    unprotected: int | None
    first_unprotected: Fraction | None
    pressure: Fraction | None
    pressure_density: Fraction | None
    ps_product: int
    ps_seed: tuple[int, ...]
    ps_ones: int


def fmt_frac(value: Fraction | None) -> str:
    if value is None:
        return "-"
    return S356.fmt_frac(value)


def fmt_float(value: Fraction | None) -> str:
    if value is None:
        return "-"
    return f"{float(value):.6f}"


def factorization(n: int) -> str:
    x = n
    parts: list[str] = []
    p = 2
    while p * p <= x:
        if x % p == 0:
            e = 0
            while x % p == 0:
                x //= p
                e += 1
            parts.append(f"{p}^{e}" if e > 1 else str(p))
        p += 1 if p == 2 else 2
    if x > 1:
        parts.append(str(x))
    return "*".join(parts)


def divisors(n: int) -> tuple[int, ...]:
    out: list[int] = []
    for d in range(1, isqrt(n) + 1):
        if n % d == 0:
            out.append(d)
            if d * d != n:
                out.append(n // d)
    return tuple(sorted(out))


def proper_divisors(n: int) -> tuple[int, ...]:
    return tuple(d for d in divisors(n) if 1 < d < n)


def phi(n: int) -> int:
    return sum(1 for a in range(1, n) if gcd(a, n) == 1)


def primitive(speeds: tuple[int, ...]) -> bool:
    g = 0
    for v in speeds:
        g = gcd(g, v)
    return g == 1


def lpd(n: int) -> int | None:
    ds = proper_divisors(n)
    return max(ds) if ds else None


def lpd_ladder(n: int, d: int) -> tuple[int, tuple[int, ...]]:
    best: tuple[Fraction, int, Fraction, int, tuple[int, ...]] | None = None
    for skip in range(1, n):
        speeds = tuple(sorted({1} | {d * q for q in range(1, n) if q != skip}))
        if len(speeds) != n - 1 or not primitive(speeds):
            continue
        report = S356.report("rank", list(speeds))
        rank = (
            report.max_gap / report.threshold,
            report.boundary_witness_count,
            ONE_MINUS(report.forbidden_length),
            skip,
            speeds,
        )
        if best is None or rank < best:
            best = rank
    if best is None:
        raise RuntimeError(f"no lpd ladder found for n={n}, d={d}")
    return best[3], best[4]


def ONE_MINUS(value: Fraction) -> Fraction:
    return Fraction(1, 1) - value


def build_rows() -> tuple[LadderRow, ...]:
    ps = OPS.enumerate_witnesses(N_MAX)
    rows: list[LadderRow] = []
    for n in range(N_MIN, N_MAX + 1):
        best_ps = sorted(ps[n], key=lambda w: (w.product, w.length, w.seed))[0]
        d = lpd(n)
        if d is None:
            rows.append(
                LadderRow(
                    n=n,
                    factors=factorization(n),
                    phi=phi(n),
                    tau=len(divisors(n)),
                    unit_density=Fraction(phi(n), n - 1),
                    lpd=None,
                    skip=None,
                    classification="-",
                    gap_ratio=None,
                    unprotected=None,
                    first_unprotected=None,
                    pressure=None,
                    pressure_density=None,
                    ps_product=best_ps.product,
                    ps_seed=best_ps.seed,
                    ps_ones=best_ps.ones,
                )
            )
            continue

        skip, speeds = lpd_ladder(n, d)
        summary = S360.summarize(list(speeds))
        gap_ratio = summary.max_gap / summary.threshold
        pressure = None if gap_ratio == 0 else Fraction(summary.unprotected_count, 1) / gap_ratio
        pressure_density = None if pressure is None else pressure / (n * n)
        rows.append(
            LadderRow(
                n=n,
                factors=factorization(n),
                phi=phi(n),
                tau=len(divisors(n)),
                unit_density=Fraction(phi(n), n - 1),
                lpd=d,
                skip=skip,
                classification=summary.classification,
                gap_ratio=gap_ratio,
                unprotected=summary.unprotected_count,
                first_unprotected=summary.first_unprotected,
                pressure=pressure,
                pressure_density=pressure_density,
                ps_product=best_ps.product,
                ps_seed=best_ps.seed,
                ps_ones=best_ps.ones,
            )
        )
    return tuple(rows)


def print_rows(rows: tuple[LadderRow, ...]) -> None:
    print("Exact lpd-ladder triage for denominators above 14")
    print("n is the forbidden denominator; moving speeds k=n-1.")
    print()
    print(
        "n factors phi/(n-1) tau lpd skip class gap/th unprot "
        "pressure/n^2 ps_min(seed)"
    )
    print("-" * 108)
    for row in rows:
        print(
            f"{row.n:2d} {row.factors:>7} {fmt_float(row.unit_density):>10} "
            f"{row.tau:3d} {str(row.lpd or '-'):>3} {str(row.skip or '-'):>4} "
            f"{row.classification:>12} {fmt_float(row.gap_ratio):>9} "
            f"{str(row.unprotected or '-'):>6} {fmt_float(row.pressure_density):>12} "
            f"{row.ps_product}:{row.ps_seed}+1^{row.ps_ones}"
        )
    print()


def print_rankings(rows: tuple[LadderRow, ...]) -> None:
    composites = [row for row in rows if row.pressure is not None]
    print("Rankings")
    print("  raw endpoint/gap pressure:")
    for row in sorted(composites, key=lambda r: r.pressure or 0, reverse=True)[:6]:
        print(
            f"    n={row.n:2d} pressure={fmt_float(row.pressure)} "
            f"gap/th={fmt_float(row.gap_ratio)} unprot={row.unprotected}"
        )
    print("  normalized pressure pressure/n^2:")
    for row in sorted(composites, key=lambda r: r.pressure_density or 0, reverse=True)[:6]:
        print(
            f"    n={row.n:2d} pressure/n^2={fmt_float(row.pressure_density)} "
            f"factors={row.factors} lpd={row.lpd}"
        )
    print("  lowest unit-skeleton density:")
    for row in sorted(rows, key=lambda r: r.unit_density)[:6]:
        print(
            f"    n={row.n:2d} phi/(n-1)={fmt_float(row.unit_density)} "
            f"factors={row.factors} tau={row.tau}"
        )
    print()


def print_pick(rows: tuple[LadderRow, ...]) -> None:
    row = next(r for r in rows if r.n == FAVORITE_N)
    print(f"My favorite next denominator: n={row.n}")
    print(
        f"  This is {row.n} total runners in the repo convention, "
        f"or {row.n - 1} moving speeds against a stationary runner."
    )
    print(
        "  Why: it is still computationally close to 14, but it is not just "
        "the same anomaly repeated."
    )
    print(
        f"  Its lpd ladder uses d={row.lpd}, has gap/th={fmt_float(row.gap_ratio)}, "
        f"unprotected={row.unprotected}, and pressure/n^2={fmt_float(row.pressure_density)}."
    )
    print(
        "  Structurally, 18=2*3^2 combines an even gate, a square 3-torsion "
        "layer, and a largest-proper-divisor layer 9.  That gives both a "
        "clean quotient ladder and a CRT/torsion interface."
    )
    print(
        f"  Its product-sum comparison is m(18)={row.ps_product} from seed "
        f"{row.ps_seed} with {row.ps_ones} ones, so the natural-operation side "
        "is already multi-factor rather than purely binary."
    )
    print(
        "  My read: n=16 is the clean proof laboratory, n=21 is the pretty "
        "3-by-7 transfer, and n=24 is the stress test.  But n=18 is the best "
        "next battlefield: small enough to certify, rich enough to expose the "
        "recursive mechanism."
    )


def main() -> None:
    print("Above-14 Lonely Runner candidate triage (codex-2026-05-31 S383)")
    print()
    rows = build_rows()
    print_rows(rows)
    print_rankings(rows)
    print_pick(rows)


if __name__ == "__main__":
    main()
