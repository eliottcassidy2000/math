#!/usr/bin/env python3
"""
lrc_operation_shadow_bridge_s377.py

codex-2026-05-31 S377

Bridge two strands from the prompt:

1. Natural-number operation shadows:
     x+y=z forgets to the complete transitive order x -> z iff x<z.
     x*y=z forgets to the sparse divisibility DAG x -> z iff x|z.

2. Lonely Runner recursion in the threshold denominator n:
     primes have a full unit endpoint skeleton, while composites open
     quotient ladders and nonunit scalar-puncture layers.

The S376 LRC atlas was already computed exactly and stored.  This script
parses that atlas, computes fresh operation-shadow metrics, and runs a small
new exact scan of the largest-proper-divisor ladder channel through the stored
atlas range.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from math import gcd, sqrt
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
NAT = SourceFileLoader(
    "natural_operation_digraphs_s365",
    str(ROOT / "04-computation" / "natural_operation_digraphs_s365.py"),
).load_module()
S356 = SourceFileLoader(
    "lonely_runner_residue_probe_s356",
    str(ROOT / "04-computation" / "lonely_runner_residue_probe_s356.py"),
).load_module()
S360 = SourceFileLoader(
    "lonely_runner_endpoint_protection_s360",
    str(ROOT / "04-computation" / "lonely_runner_endpoint_protection_s360.py"),
).load_module()


N_MIN = 2
N_MAX = 18
OPERATION_MAX = 24
LPD_SCAN_MAX = 18
S376_OUT = ROOT / "05-knowledge" / "results" / "lonely_runner_recursive_metrics_s376.out"


@dataclass(frozen=True)
class OperationRow:
    n: int
    add_edges: int
    mul_edges_unit: int
    mul_edges_nonunit: int
    mul_density: Fraction
    proper_divisors: tuple[int, ...]
    largest_proper_divisor: int | None
    smallest_prime_factor: int | None
    product_sum_product: int | None
    product_sum_seed: tuple[int, ...] | None
    product_sum_ones: int | None
    two_factor_min: int | None
    multifactor_win: bool


@dataclass(frozen=True)
class StoredLRCRow:
    n: int
    factors: str
    patterns: int
    unit_density: Fraction
    puncture_missed: int
    puncture_density: Fraction
    puncture_delta_gcds: str
    ladder_d: int | None
    ladder_gap_ratio: Fraction | None
    ladder_unprotected: int | None


@dataclass(frozen=True)
class BridgeRow:
    operation: OperationRow
    lrc: StoredLRCRow
    ladder_lpd_match: bool | None
    quotient_pressure: Fraction | None


@dataclass(frozen=True)
class LadderScanRow:
    n: int
    lpd: int
    skip: int
    gap_ratio: Fraction
    unprotected: int
    first_unprotected: Fraction | None
    speeds: tuple[int, ...]


def proper_divisors(n: int) -> tuple[int, ...]:
    return tuple(d for d in range(2, n) if n % d == 0)


def largest_proper_divisor(n: int) -> int | None:
    divs = proper_divisors(n)
    return divs[-1] if divs else None


def smallest_prime_factor(n: int) -> int | None:
    lpd = largest_proper_divisor(n)
    if lpd is None:
        return None
    return n // lpd


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


def parse_decimal_fraction(text: str) -> Fraction:
    return Fraction(text)


def parse_ladder_d(label: str | None) -> int | None:
    if not label or label == "-":
        return None
    match = re.search(r"d=(\d+)", label)
    return int(match.group(1)) if match else None


def fmt_float(x: Fraction | None) -> str:
    if x is None:
        return "-"
    return f"{float(x):.6f}"


def fmt_pressure(x: Fraction | None) -> str:
    if x is None:
        return "-"
    return f"{float(x):.1f}"


def fmt_frac(x: Fraction | None) -> str:
    if x is None:
        return "-"
    return S356.fmt_frac(x)


def fmt_optional_int(x: int | None) -> str:
    return "-" if x is None else str(x)


def primitive(speeds: tuple[int, ...]) -> bool:
    g = 0
    for speed in speeds:
        g = gcd(g, speed)
    return g == 1


def pearson(xs: list[float], ys: list[float]) -> float | None:
    if len(xs) != len(ys) or len(xs) < 2:
        return None
    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)
    vx = sum((x - mx) ** 2 for x in xs)
    vy = sum((y - my) ** 2 for y in ys)
    if vx == 0 or vy == 0:
        return None
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    return cov / sqrt(vx * vy)


def operation_rows(max_n: int) -> dict[int, OperationRow]:
    seed_rows = NAT.enumerate_product_sum_seeds(
        max_k=max_n,
        max_product=2 * max_n * max_n,
    )
    rows: dict[int, OperationRow] = {}
    for n in range(N_MIN, max_n + 1):
        add_edges = len(NAT.additive_edges(n))
        mul_edges_unit = len(NAT.product_edges(n, include_unit=True))
        mul_edges_nonunit = len(NAT.product_edges(n, include_unit=False))
        best = seed_rows.get(n, [None])[0]
        two_factor = NAT.two_factor_solutions_for_arity(n)
        two_min = min((p for _, _, p in two_factor), default=None)
        rows[n] = OperationRow(
            n=n,
            add_edges=add_edges,
            mul_edges_unit=mul_edges_unit,
            mul_edges_nonunit=mul_edges_nonunit,
            mul_density=Fraction(mul_edges_unit, add_edges) if add_edges else Fraction(0),
            proper_divisors=proper_divisors(n),
            largest_proper_divisor=largest_proper_divisor(n),
            smallest_prime_factor=smallest_prime_factor(n),
            product_sum_product=best.product if best else None,
            product_sum_seed=best.seed if best else None,
            product_sum_ones=best.ones if best else None,
            two_factor_min=two_min,
            multifactor_win=bool(
                best
                and two_min is not None
                and len(best.seed) > 2
                and best.product < two_min
            ),
        )
    return rows


def load_stored_lrc_rows(path: Path = S376_OUT) -> dict[int, StoredLRCRow]:
    rows: dict[int, StoredLRCRow] = {}
    ladder_d: dict[int, int | None] = {}
    ladder_gap: dict[int, Fraction | None] = {}
    ladder_unprotected: dict[int, int | None] = {}
    current_n: int | None = None
    in_core = False

    core_re = re.compile(
        r"^\s*(\d+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\d+)\s+(\d+)\s+"
        r"([0-9.]+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([0-9.]+)\s+(\S+)\s+(\S+)"
    )
    n_re = re.compile(r"^n=\s*(\d+)\s+gate:")
    ladder_re = re.compile(r"^\s+ladder:\s+(.*)$")

    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("Core recursive LRC metrics"):
            in_core = True
            continue
        if in_core and not line.strip():
            in_core = False
            continue
        if in_core:
            match = core_re.match(line)
            if match:
                n = int(match.group(1))
                rows[n] = StoredLRCRow(
                    n=n,
                    factors=match.group(2),
                    patterns=int(match.group(3)),
                    unit_density=Fraction(int(match.group(6)), n - 1),
                    puncture_missed=int(match.group(10)),
                    puncture_density=parse_decimal_fraction(match.group(11)),
                    puncture_delta_gcds=match.group(13),
                    ladder_d=None,
                    ladder_gap_ratio=None,
                    ladder_unprotected=None,
                )
            continue

        n_match = n_re.match(line)
        if n_match:
            current_n = int(n_match.group(1))
            continue

        ladder_match = ladder_re.match(line)
        if ladder_match and current_n is not None:
            text = ladder_match.group(1)
            if text == "-":
                ladder_d[current_n] = None
                ladder_gap[current_n] = None
                ladder_unprotected[current_n] = None
            else:
                ladder_d[current_n] = parse_ladder_d(text)
                gap_match = re.search(r"gap/th=([0-9.]+)", text)
                unprot_match = re.search(r"unprot=(\d+)", text)
                ladder_gap[current_n] = (
                    parse_decimal_fraction(gap_match.group(1)) if gap_match else None
                )
                ladder_unprotected[current_n] = (
                    int(unprot_match.group(1)) if unprot_match else None
                )

    for n, row in list(rows.items()):
        rows[n] = StoredLRCRow(
            n=row.n,
            factors=row.factors,
            patterns=row.patterns,
            unit_density=row.unit_density,
            puncture_missed=row.puncture_missed,
            puncture_density=row.puncture_density,
            puncture_delta_gcds=row.puncture_delta_gcds,
            ladder_d=ladder_d.get(n),
            ladder_gap_ratio=ladder_gap.get(n),
            ladder_unprotected=ladder_unprotected.get(n),
        )
    return rows


def bridge_rows(operation: dict[int, OperationRow], lrc: dict[int, StoredLRCRow]) -> list[BridgeRow]:
    out: list[BridgeRow] = []
    for n in range(N_MIN, N_MAX + 1):
        op = operation[n]
        row = lrc[n]
        pressure = None
        if row.ladder_gap_ratio and row.ladder_unprotected is not None:
            pressure = Fraction(row.ladder_unprotected, 1) / row.ladder_gap_ratio
        out.append(
            BridgeRow(
                operation=op,
                lrc=row,
                ladder_lpd_match=None if row.ladder_d is None else row.ladder_d == op.largest_proper_divisor,
                quotient_pressure=pressure,
            )
        )
    return out


def best_lpd_ladder(n: int) -> LadderScanRow | None:
    lpd = largest_proper_divisor(n)
    if lpd is None:
        return None
    candidates: list[LadderScanRow] = []
    for skip in range(1, n):
        speeds = tuple(sorted({1} | {lpd * q for q in range(1, n) if q != skip}))
        if len(speeds) != n - 1 or not primitive(speeds):
            continue
        summary = S360.summarize(list(speeds))
        candidates.append(
            LadderScanRow(
                n=n,
                lpd=lpd,
                skip=skip,
                gap_ratio=summary.max_gap / summary.threshold,
                unprotected=summary.unprotected_count,
                first_unprotected=summary.first_unprotected,
                speeds=summary.speeds,
            )
        )
    if not candidates:
        return None
    return min(candidates, key=lambda row: (row.gap_ratio, row.unprotected, row.skip))


def print_operation_shadow_table(operation: dict[int, OperationRow]) -> None:
    print("Natural operation shadow prefixes")
    print(
        "n  addE  mulE+1  mulE>1  mul_density  lpd  spf  "
        "ps_minP  ps_seed  ones  multi_win"
    )
    for n in range(N_MIN, OPERATION_MAX + 1):
        row = operation[n]
        seed = "-" if row.product_sum_seed is None else ".".join(map(str, row.product_sum_seed))
        print(
            f"{n:>2} {row.add_edges:>5} {row.mul_edges_unit:>7} "
            f"{row.mul_edges_nonunit:>7} {fmt_float(row.mul_density):>12} "
            f"{fmt_optional_int(row.largest_proper_divisor):>4} "
            f"{fmt_optional_int(row.smallest_prime_factor):>4} "
            f"{fmt_optional_int(row.product_sum_product):>7} "
            f"{seed:>8} {fmt_optional_int(row.product_sum_ones):>5} "
            f"{str(row.multifactor_win):>9}"
        )
    print()


def print_bridge_table(rows: list[BridgeRow]) -> None:
    print("Operation shadow / stored LRC denominator bridge")
    print(
        "n  factor  mul_density  unit_density  C(n)  dC  moat  "
        "delta_gcd  ladder_d  lpd  match  gap/th   unprot  pressure"
    )
    prev_c = None
    for row in rows:
        lrc = row.lrc
        op = row.operation
        dc = "-" if prev_c is None else str(lrc.patterns - prev_c)
        prev_c = lrc.patterns
        print(
            f"{lrc.n:>2} {lrc.factors:>7} {fmt_float(op.mul_density):>12} "
            f"{fmt_float(lrc.unit_density):>13} {lrc.patterns:>5} {dc:>4} "
            f"{lrc.puncture_missed:>5} {lrc.puncture_delta_gcds:>9} "
            f"{fmt_optional_int(lrc.ladder_d):>9} "
            f"{fmt_optional_int(op.largest_proper_divisor):>4} "
            f"{str(row.ladder_lpd_match):>5} {fmt_float(lrc.ladder_gap_ratio):>8} "
            f"{fmt_optional_int(lrc.ladder_unprotected):>7} "
            f"{fmt_pressure(row.quotient_pressure):>9}"
        )
    print()


def print_lpd_ladder_scan(stored: dict[int, StoredLRCRow]) -> None:
    print(f"Exact largest-proper-divisor ladder scan, n<= {LPD_SCAN_MAX}")
    print("n  factors  lpd  best_skip  gap/th   unprot  first_unprot  stored_best_d")
    for n in range(4, LPD_SCAN_MAX + 1):
        lpd = largest_proper_divisor(n)
        if lpd is None:
            continue
        row = best_lpd_ladder(n)
        stored_d = stored.get(n).ladder_d if n in stored else None
        print(
            f"{n:>2} {factorization(n):>7} {lpd:>4} "
            f"{fmt_optional_int(row.skip if row else None):>9} "
            f"{fmt_float(row.gap_ratio if row else None):>8} "
            f"{fmt_optional_int(row.unprotected if row else None):>7} "
            f"{fmt_frac(row.first_unprotected if row else None):>12} "
            f"{fmt_optional_int(stored_d):>13}"
        )
    print()


def print_transition_table(rows: list[BridgeRow]) -> None:
    print("Transition signatures: operation residue versus LRC layers")
    print("edge  dC  dmul_density  dunit_density  dmoat_density  event")
    for left, right in zip(rows, rows[1:]):
        l0, l1 = left.lrc, right.lrc
        op0, op1 = left.operation, right.operation
        events = []
        if op1.largest_proper_divisor is None:
            events.append("prime")
        else:
            events.append(f"lpd={op1.largest_proper_divisor}")
            if right.lrc.ladder_d == op1.largest_proper_divisor:
                events.append("ladder=lpd")
        if right.lrc.puncture_delta_gcds != "1":
            events.append("nonunit_moat")
        if right.lrc.ladder_gap_ratio and right.lrc.ladder_gap_ratio < Fraction(1, 50):
            events.append("tiny_ladder")
        print(
            f"{l0.n:>2}->{l1.n:<2} {l1.patterns - l0.patterns:>4} "
            f"{fmt_float(op1.mul_density - op0.mul_density):>13} "
            f"{fmt_float(l1.unit_density - l0.unit_density):>14} "
            f"{fmt_float(l1.puncture_density - l0.puncture_density):>14} "
            f"{';'.join(events)}"
        )
    print()


def print_pressure_rankings(rows: list[BridgeRow]) -> None:
    print("Composite pressure rankings")
    pressure_rows = sorted(
        [row for row in rows if row.quotient_pressure is not None],
        key=lambda row: row.quotient_pressure or Fraction(0),
        reverse=True,
    )
    print("By endpoint_defect / gap_ratio:")
    for row in pressure_rows[:8]:
        lrc = row.lrc
        print(
            f"  n={lrc.n:>2} d={lrc.ladder_d} "
            f"gap/th={fmt_float(lrc.ladder_gap_ratio)} "
            f"unprot={lrc.ladder_unprotected:>3} "
            f"pressure={fmt_pressure(row.quotient_pressure)}"
        )
    print()


def print_correlations(rows: list[BridgeRow]) -> None:
    print("Small-n correlation checks (descriptive only)")
    mul_density = [float(row.operation.mul_density) for row in rows]
    complexity = [float(row.lrc.patterns) for row in rows]
    unit_density = [float(row.lrc.unit_density) for row in rows]
    moat_density = [float(row.lrc.puncture_density) for row in rows]
    print(f"corr(mul_density, C(n))={pearson(mul_density, complexity):.6f}")
    print(f"corr(mul_density, unit_density)={pearson(mul_density, unit_density):.6f}")
    print(f"corr(unit_density, moat_density)={pearson(unit_density, moat_density):.6f}")

    comp = [row for row in rows if row.operation.largest_proper_divisor is not None and row.lrc.ladder_gap_ratio]
    lpd_ratio = [
        row.operation.largest_proper_divisor / row.lrc.n
        for row in comp
        if row.operation.largest_proper_divisor is not None
    ]
    ladder_gap = [float(row.lrc.ladder_gap_ratio) for row in comp if row.lrc.ladder_gap_ratio]
    unprotected = [float(row.lrc.ladder_unprotected) for row in comp if row.lrc.ladder_unprotected is not None]
    print(f"corr(lpd/n, ladder_gap_ratio)={pearson(lpd_ratio, ladder_gap):.6f}")
    print(f"corr(lpd/n, ladder_unprotected)={pearson(lpd_ratio, unprotected):.6f}")
    print()


def print_synthesis(rows: list[BridgeRow], operation: dict[int, OperationRow]) -> None:
    lpd_matches = [row for row in rows if row.ladder_lpd_match is not None]
    match_count = sum(1 for row in lpd_matches if row.ladder_lpd_match)
    tiny = [
        row.lrc.n
        for row in rows
        if row.lrc.ladder_gap_ratio and row.lrc.ladder_gap_ratio < Fraction(1, 50)
    ]
    multifactor_wins = [n for n, row in operation.items() if row.multifactor_win]
    print("Synthesis")
    print(
        "  1. Addition supplies the ambient transitive completion; multiplication "
        "is the residue that survives forgetting operation fibers."
    )
    print(
        "  2. LRC composites behave as if the same multiplicative residue opens "
        "quotient channels.  In the stored exact rows n=2..18, best quotient "
        f"ladders match the largest proper divisor in {match_count}/{len(lpd_matches)} "
        "composite cases with a ladder."
    )
    print(
        f"  3. The exact lpd-only ladder scan through n={LPD_SCAN_MAX} extends the same channel "
        "as a concrete family, even where a full all-divisor search has not been "
        "rerun in this script."
    )
    print(
        "  4. Product-sum ones and LRC unit endpoints play the same formal role: "
        "they are additive slack/boundary mass that repairs a multiplicative "
        "defect but also records exactly where the repair can leak."
    )
    print(
        f"  5. Multi-factor product-sum wins begin at n={multifactor_wins[:8]}, "
        "while tiny quotient-ladder gaps below 0.02 occur at "
        f"n={tiny}.  These are the first two places to compare operation-complex "
        "cells with LRC endpoint-transfer matrices."
    )


def main() -> None:
    operation = operation_rows(OPERATION_MAX)
    stored = load_stored_lrc_rows()
    rows = bridge_rows(operation, stored)
    print("LRC operation-shadow bridge (codex-2026-05-31 S377)")
    print("=" * 78)
    print()
    print_operation_shadow_table(operation)
    print_bridge_table(rows)
    print_lpd_ladder_scan(stored)
    print_transition_table(rows)
    print_pressure_rankings(rows)
    print_correlations(rows)
    print_synthesis(rows, operation)


if __name__ == "__main__":
    main()
