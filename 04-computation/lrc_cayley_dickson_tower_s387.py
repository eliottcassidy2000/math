#!/usr/bin/env python3
"""
lrc_cayley_dickson_tower_s387.py

codex-2026-05-31 S387

Explore a Cayley-Dickson-style reading of the Lonely Runner denominator tower.

The analogy used here:

* Cayley-Dickson doubles dimension and loses algebraic structure at each
  stage: order, commutativity, associativity, alternativity/division.
* LRC denominators split as n = 2^r * odd_core.  The 2^r part is the
  doubling row; the odd core is the torsion payload carried by that row.
* Inserting the mandatory n-gate closes the unit boundary layer, but creates
  descendant endpoint debt.  That debt plays the same heuristic role as
  zero-divisor leakage in late Cayley-Dickson levels.
"""

from __future__ import annotations

from collections import Counter
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


N_VALUES = tuple(range(14, 25))
ZOOM_N = 18
ONE = Fraction(1, 1)


CD_LEVELS = {
    0: ("R", "ordered division baseline"),
    1: ("C", "order lost; phase appears"),
    2: ("H", "commutativity lost; handedness appears"),
    3: ("O", "associativity lost; Fano-plane control"),
    4: ("S", "division lost; zero-divisor debt"),
    5: ("P", "higher zero-divisor jungle"),
}


@dataclass(frozen=True)
class TowerRow:
    n: int
    factors: str
    v2: int
    odd_core: int
    cd_symbol: str
    cd_meaning: str
    phi: int
    unit_density: Fraction
    tau: int
    lpd: int | None
    lpd_gap_ratio: Fraction | None
    lpd_unprotected: int | None
    debt_per_unit: Fraction | None
    layer_hist: tuple[tuple[int, int], ...]
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


def totient(n: int) -> int:
    return sum(1 for a in range(1, n) if gcd(a, n) == 1)


def v2_and_odd(n: int) -> tuple[int, int]:
    r = 0
    x = n
    while x % 2 == 0:
        r += 1
        x //= 2
    return r, x


def primitive(speeds: tuple[int, ...]) -> bool:
    g = 0
    for v in speeds:
        g = gcd(g, v)
    return g == 1


def lpd(n: int) -> int | None:
    ds = proper_divisors(n)
    return max(ds) if ds else None


def best_lpd_ladder(n: int, d: int) -> tuple[int, tuple[int, ...]]:
    best: tuple[Fraction, int, Fraction, int, tuple[int, ...]] | None = None
    for skip in range(1, n):
        speeds = tuple(sorted({1} | {d * q for q in range(1, n) if q != skip}))
        if len(speeds) != n - 1 or not primitive(speeds):
            continue
        report = S356.report("rank", list(speeds))
        rank = (
            report.max_gap / report.threshold,
            report.boundary_witness_count,
            ONE - report.forbidden_length,
            skip,
            speeds,
        )
        if best is None or rank < best:
            best = rank
    if best is None:
        raise RuntimeError(f"no lpd ladder for n={n}")
    return best[3], best[4]


def unprotected_points(speeds: tuple[int, ...]) -> list[Fraction]:
    points = {endpoint.value for endpoint in S360.endpoints(speeds)}
    return [
        point
        for point in sorted(points)
        if not any(S360.direct_protects(speeds, protector, point) for protector in speeds)
    ]


def debt_layer(n: int, point: Fraction) -> int:
    """Layer 1 is the unit boundary; lpd ladders often expose layer lpd."""

    return point.denominator // gcd(point.denominator, n)


def layer_histogram(n: int, speeds: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    return tuple(sorted(Counter(debt_layer(n, point) for point in unprotected_points(speeds)).items()))


def build_rows() -> tuple[TowerRow, ...]:
    ps = OPS.enumerate_witnesses(max(N_VALUES))
    rows: list[TowerRow] = []
    for n in N_VALUES:
        r, odd = v2_and_odd(n)
        symbol, meaning = CD_LEVELS.get(r, (f"2^{r}", "beyond named Cayley-Dickson levels"))
        best_ps = sorted(ps[n], key=lambda w: (w.product, w.length, w.seed))[0]
        d = lpd(n)
        if d is None:
            rows.append(
                TowerRow(
                    n=n,
                    factors=factorization(n),
                    v2=r,
                    odd_core=odd,
                    cd_symbol=symbol,
                    cd_meaning=meaning,
                    phi=totient(n),
                    unit_density=Fraction(totient(n), n - 1),
                    tau=len(divisors(n)),
                    lpd=None,
                    lpd_gap_ratio=None,
                    lpd_unprotected=None,
                    debt_per_unit=None,
                    layer_hist=(),
                    ps_product=best_ps.product,
                    ps_seed=best_ps.seed,
                    ps_ones=best_ps.ones,
                )
            )
            continue

        _skip, speeds = best_lpd_ladder(n, d)
        summary = S360.summarize(list(speeds))
        gap_ratio = summary.max_gap / summary.threshold
        rows.append(
            TowerRow(
                n=n,
                factors=factorization(n),
                v2=r,
                odd_core=odd,
                cd_symbol=symbol,
                cd_meaning=meaning,
                phi=totient(n),
                unit_density=Fraction(totient(n), n - 1),
                tau=len(divisors(n)),
                lpd=d,
                lpd_gap_ratio=gap_ratio,
                lpd_unprotected=summary.unprotected_count,
                debt_per_unit=Fraction(summary.unprotected_count, max(1, totient(n))),
                layer_hist=layer_histogram(n, speeds),
                ps_product=best_ps.product,
                ps_seed=best_ps.seed,
                ps_ones=best_ps.ones,
            )
        )
    return tuple(rows)


def print_cd_table() -> None:
    print("Cayley-Dickson reference tower")
    print("  level dim symbol meaning")
    for level in range(0, 6):
        symbol, meaning = CD_LEVELS[level]
        print(f"  {level:>5} {2**level:>3} {symbol:>6} {meaning}")
    print()


def print_lrc_tower(rows: tuple[TowerRow, ...]) -> None:
    print("LRC denominators as Cayley-Dickson rows with odd torsion payloads")
    print(
        "n factors 2^r odd CD phi/(n-1) tau lpd gap/th unprot debt/unit "
        "layers ps(seed)"
    )
    print("-" * 118)
    for row in rows:
        layers = ",".join(f"{layer}:{count}" for layer, count in row.layer_hist[:4])
        if len(row.layer_hist) > 4:
            layers += ",..."
        print(
            f"{row.n:2d} {row.factors:>7} {2**row.v2:>3} {row.odd_core:>3} "
            f"{row.cd_symbol:>2} {fmt_float(row.unit_density):>10} {row.tau:3d} "
            f"{str(row.lpd or '-'):>3} {fmt_float(row.lpd_gap_ratio):>8} "
            f"{str(row.lpd_unprotected or '-'):>6} {fmt_float(row.debt_per_unit):>9} "
            f"{layers or '-':>16} {row.ps_product}:{row.ps_seed}+1^{row.ps_ones}"
        )
    print()


def zoom_n18() -> None:
    n = ZOOM_N
    d = lpd(n)
    assert d is not None
    skip, speeds = best_lpd_ladder(n, d)
    summary = S360.summarize(list(speeds))
    leaks = unprotected_points(speeds)
    layer_hist = Counter(debt_layer(n, point) for point in leaks)
    denom_hist = Counter(point.denominator for point in leaks)
    residue_hist = Counter(
        int((point * n) % n) if (point * n).denominator == 1 else -1 for point in leaks
    )

    print("Zoom: n=18 lpd ladder as twisted Cayley-Dickson layer")
    print(f"  speeds={speeds}")
    print(
        f"  lpd={d} skip={skip} class={summary.classification} "
        f"gap/th={fmt_float(summary.max_gap / summary.threshold)} "
        f"unprotected={summary.unprotected_count} first={fmt_frac(summary.first_unprotected)}"
    )
    print(f"  layer_hist={tuple(sorted(layer_hist.items()))}")
    print(f"  denominator_hist_prefix={tuple(sorted(denom_hist.items())[:10])}")
    print(f"  residue_mod_18_hist={tuple(sorted(residue_hist.items()))}")
    print(f"  first_leaks={[fmt_frac(point) for point in leaks[:14]]}")
    print()


def print_synthesis(rows: tuple[TowerRow, ...]) -> None:
    favorite = next(row for row in rows if row.n == ZOOM_N)
    pure = next(row for row in rows if row.n == 16)
    stress = next(row for row in rows if row.n == 24)

    print("Synthesis")
    print(
        "  The Cayley-Dickson analogy is not that LRC literally is an algebra.  "
        "It is that every doubling row closes one old freedom while exposing a "
        "new defect layer."
    )
    print(
        f"  n=16 is the pure {pure.cd_symbol}-row test: 2-adic, clean, and close "
        "to the sedenion zero-divisor threshold."
    )
    print(
        f"  n=18 is the first attractive mixed row: CD row {2**favorite.v2} "
        f"({favorite.cd_symbol}) carrying odd core {favorite.odd_core}.  Its "
        f"unprotected endpoints concentrate at layer {favorite.lpd}, exactly "
        "the largest-proper-divisor payload."
    )
    print(
        f"  n=24 is the stress test: row {2**stress.v2} with odd core "
        f"{stress.odd_core}, lowest unit density in this window, and the loudest "
        "normalized pressure."
    )
    print(
        "  My favorite remains n=18 because it is where the tower stops being "
        "pure doubling but has not yet become a large divisor-rich mess."
    )


def main() -> None:
    print("LRC Cayley-Dickson tower analogy (codex-2026-05-31 S387)")
    print()
    rows = build_rows()
    print_cd_table()
    print_lrc_tower(rows)
    zoom_n18()
    print_synthesis(rows)


if __name__ == "__main__":
    main()
