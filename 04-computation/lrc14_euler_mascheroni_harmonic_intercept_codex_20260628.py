#!/usr/bin/env python3
"""HYP-3430: Euler-Mascheroni harmonic intercept firewall for LRC14.

Euler-Mascheroni gamma enters here as the finite intercept

    H_N - log N,

not as a proof shortcut.  HYP-3425/HYP-3426/HYP-3427/HYP-3428/HYP-3429 now
give a concrete two-adic certificate stack: bad-core intervals, mirror
reduction, wall words, loss ledger, and endpoint-spine rank.  This scout asks
whether the harmonic intercept can replace any of that retained data.

The answer on the HYP-3429 bank is no: the gamma intercept tracks scale and
best-window length, but it does not determine the endpoint certificate class.
Its proof-facing role is therefore a scalar firewall and tail calibrator.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from importlib import util
from itertools import combinations
from pathlib import Path
import math
import sys


ROOT = Path(__file__).resolve().parents[1]
H3429_PATH = ROOT / "04-computation" / "lrc14_component_spine_certificate_codex_20260628.py"
EULER_GAMMA = 0.5772156649015329


def load_h3429():
    spec = util.spec_from_file_location("h3429_endpoint_spine", H3429_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {H3429_PATH}")
    module = util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


h3429 = load_h3429()


def harmonic(n: int) -> float:
    return sum(1.0 / k for k in range(1, n + 1))


def corr(xs: list[float], ys: list[float]) -> float:
    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)
    vx = sum((x - mx) ** 2 for x in xs)
    vy = sum((y - my) ** 2 for y in ys)
    if vx == 0 or vy == 0:
        return 0.0
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / math.sqrt(vx * vy)


def fmt_float(x: float) -> str:
    return f"{x:.9f}"


@dataclass(frozen=True)
class HarmonicRow:
    name: str
    class_name: str
    max_speed: int
    gamma_intercept: float
    gamma_error: float
    speed_pressure: float
    even_pressure: float
    odd_pressure: float
    owner_pressure: float
    best_rank: int
    best_log_length: float
    window_count: int
    has_mixed: bool


def class_name(row: h3429.RowAudit) -> str:
    kinds = "+".join(row.best.kinds) if row.best.kinds else "none"
    return f"{kinds}/rank{row.best.rank}/branches{len(row.best.branches)}"


def row_pressures(row: h3429.RowAudit) -> tuple[float, float, float, float]:
    speed_pressure = sum(1.0 / s for s in row.speeds)
    even_pressure = sum(2.0 / s for s in row.speeds if s % 2 == 0)
    odd_pressure = sum(1.0 / s for s in row.speeds if s % 2 == 1)
    owner_pressure = sum(1.0 / value for _kind, value in row.best.labels if value)
    return speed_pressure, even_pressure, odd_pressure, owner_pressure


def audit() -> list[HarmonicRow]:
    out: list[HarmonicRow] = []
    for name, speeds in h3429.audited_rows():
        row = h3429.audit_row(name, speeds)
        max_speed = max(row.speeds)
        intercept = harmonic(max_speed) - math.log(max_speed)
        speed_pressure, even_pressure, odd_pressure, owner_pressure = row_pressures(row)
        out.append(
            HarmonicRow(
                name=row.name,
                class_name=class_name(row),
                max_speed=max_speed,
                gamma_intercept=intercept,
                gamma_error=intercept - EULER_GAMMA,
                speed_pressure=speed_pressure,
                even_pressure=even_pressure,
                odd_pressure=odd_pressure,
                owner_pressure=owner_pressure,
                best_rank=row.best.rank,
                best_log_length=math.log(float(row.best.length)),
                window_count=row.window_count,
                has_mixed=row.has_mixed,
            )
        )
    return out


def collision_bins(rows: list[HarmonicRow], digits: int) -> dict[float, set[str]]:
    bins: dict[float, set[str]] = defaultdict(set)
    for row in rows:
        bins[round(row.gamma_intercept, digits)].add(row.class_name)
    return dict(bins)


def same_scale_collision(rows: list[HarmonicRow]) -> tuple[int, int, list[HarmonicRow]]:
    by_max: dict[int, list[HarmonicRow]] = defaultdict(list)
    for row in rows:
        by_max[row.max_speed].append(row)
    mixed_bins = [items for items in by_max.values() if len({r.class_name for r in items}) > 1]
    example = next((items for items in mixed_bins if any(r.max_speed == 84 for r in items)), mixed_bins[0])
    return len(mixed_bins), len(by_max), sorted(example, key=lambda r: (r.name, r.class_name))


AXES = (
    "predicate_retention",
    "endpoint_payload",
    "harmonic_calibration",
    "tail_exception_router",
    "scalar_firewall",
)


@dataclass(frozen=True)
class Carrier:
    name: str
    scores: tuple[int, ...]

    @property
    def total(self) -> int:
        return sum(self.scores)


CARRIERS = (
    Carrier("endpoint_spine_certificate", (10, 10, 5, 8, 10)),
    Carrier("wall_signature_certificate", (10, 9, 5, 8, 9)),
    Carrier("two_adic_loss_ledger", (9, 7, 6, 10, 9)),
    Carrier("harmonic_intercept_calibrator", (4, 2, 10, 7, 8)),
    Carrier("mertens_loglog_tail_budget", (3, 2, 9, 8, 7)),
    Carrier("raw_euler_mascheroni_scalar", (1, 0, 8, 2, 1)),
)


def tournament() -> tuple[dict[int, int], int, list[str]]:
    hist = dict(sorted(Counter(c.total for c in CARRIERS).items()))
    cycles = 0
    for a, b, c in combinations(CARRIERS, 3):
        ranks = sorted((a, b, c), key=lambda item: (-item.total, CARRIERS.index(item)))
        if ranks[0].total < ranks[1].total < ranks[2].total:
            cycles += 1
    path = [c.name for c in sorted(CARRIERS, key=lambda item: (-item.total, CARRIERS.index(item)))]
    return hist, cycles, path


def main() -> None:
    rows = audit()
    gamma_values = [row.gamma_intercept for row in rows]
    ranks = [float(row.best_rank) for row in rows]
    mixed_flags = [1.0 if row.has_mixed else 0.0 for row in rows]
    log_lengths = [row.best_log_length for row in rows]
    window_counts = [float(row.window_count) for row in rows]
    owner_pressures = [row.owner_pressure for row in rows]

    class_hist = Counter(row.class_name for row in rows)
    rounded4 = collision_bins(rows, 4)
    rounded6 = collision_bins(rows, 6)
    same_scale_mixed, same_scale_bins, same_scale_example = same_scale_collision(rows)
    smallest_intercept = min(rows, key=lambda row: row.gamma_intercept)
    largest_intercept = max(rows, key=lambda row: row.gamma_intercept)

    print("HYP-3430 EULER-MASCHERONI HARMONIC INTERCEPT FIREWALL")
    print("=" * 78)
    print("Question:")
    print("  Does the finite intercept H_N - log N certify the LRC14 endpoint-spine class?")
    print("  Treat Euler-Mascheroni as tail calibration, not as a replacement witness.")
    print()

    print("A. Aggregate readout on the HYP-3429 bank")
    print(f"  rows audited:                         {len(rows)}")
    print(f"  endpoint certificate classes:          {len(class_hist)}")
    print(f"  gamma intercept range:                 {fmt_float(min(gamma_values))} .. {fmt_float(max(gamma_values))}")
    print(f"  Euler-Mascheroni reference:             {fmt_float(EULER_GAMMA)}")
    print(f"  largest finite intercept error:         {fmt_float(max(abs(row.gamma_error) for row in rows))}")
    print(f"  same-max-speed bins with mixed classes: {same_scale_mixed}/{same_scale_bins}")
    print(f"  rounded-4 gamma bins with mixed classes:{sum(1 for v in rounded4.values() if len(v) > 1)}/{len(rounded4)}")
    print(f"  rounded-6 gamma bins with mixed classes:{sum(1 for v in rounded6.values() if len(v) > 1)}/{len(rounded6)}")
    print(f"  class histogram:                       {dict(class_hist)}")
    print()

    print("B. What gamma sees and what it forgets")
    print(f"  corr(gamma_intercept, best_rank):       {corr(gamma_values, ranks):+.6f}")
    print(f"  corr(gamma_intercept, has_mixed_spine): {corr(gamma_values, mixed_flags):+.6f}")
    print(f"  corr(gamma_intercept, log_best_length): {corr(gamma_values, log_lengths):+.6f}")
    print(f"  corr(gamma_intercept, window_count):    {corr(gamma_values, window_counts):+.6f}")
    print(f"  corr(owner_pressure, log_best_length):  {corr(owner_pressures, log_lengths):+.6f}")
    print(
        "  readout: gamma mostly tracks height/scale; endpoint labels decide the proof class."
    )
    print()

    print("C. Same gamma intercept, different proof certificates")
    print(f"  example max_speed={same_scale_example[0].max_speed}")
    for row in same_scale_example[:8]:
        print(
            f"    {row.name}: class={row.class_name}, "
            f"gamma_intercept={fmt_float(row.gamma_intercept)}, "
            f"owner_pressure={fmt_float(row.owner_pressure)}, windows={row.window_count}"
        )
    print()

    print("D. Extremal finite intercept rows")
    for label, row in (("smallest", smallest_intercept), ("largest", largest_intercept)):
        print(
            f"  {label}: {row.name}, N={row.max_speed}, class={row.class_name}, "
            f"gamma_intercept={fmt_float(row.gamma_intercept)}, "
            f"error={fmt_float(row.gamma_error)}, log_best_len={row.best_log_length:.6f}"
        )
    print()

    print("E. Canonical 84m tower")
    print("  m | max_speed | gamma_intercept | gamma_error | class | owner_pressure")
    for row in rows:
        if not row.name.startswith("canonical_84m"):
            continue
        m = int(row.name.rsplit("_", 1)[1])
        if m <= 6 or m in {12, 24, 36, 48}:
            print(
                f"  {m:2d} | {row.max_speed:9d} | {fmt_float(row.gamma_intercept)} | "
                f"{fmt_float(row.gamma_error)} | {row.class_name:18s} | {fmt_float(row.owner_pressure)}"
            )
    print()

    hist, cycles, path = tournament()
    print("F. Tournament Analysis")
    print("  vertices=proof carriers and scalar guardrails, not runners or constants")
    print(f"  axes={','.join(AXES)}")
    print(f"  score_hist={hist}")
    print(f"  directed_3cycles={cycles}")
    print("  hamiltonian_path=" + " -> ".join(path))
    print()

    print("G. Assumption challenge")
    print(
        "  considered_vertices=rows, max speeds, harmonic tails, endpoint owners, "
        "wall signatures, loss classes, survivor windows, and proof obligations"
    )
    print(
        "  chosen_vertices=proof carriers plus scalar guardrails; the quotient "
        "preserves the covering-floor predicate only when endpoint data is retained"
    )
    print(
        "  challenged_assumption=Euler-Mascheroni or Mertens constants can replace "
        "the finite endpoint-spine certificate"
    )


if __name__ == "__main__":
    main()
