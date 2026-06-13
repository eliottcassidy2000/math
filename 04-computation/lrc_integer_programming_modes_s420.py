#!/usr/bin/env python3
"""
lrc_integer_programming_modes_s420.py

codex-2026-05-31 S420

Look at the Lonely Runner endpoint program as a small integer program.

The repo's natural-number mode thread split the naturals into horizontal
odd-core chains (x -> x+2) and vertical doubling chains (x -> 2x).  In LRC
language that says a speed column should remember

    speed = 2^h * odd_core.

This script does not solve the full nonlinear LRC feasibility problem.  It
extracts exact finite IP subproblems that any counterexample must satisfy:

* small-denominator invoices from THM-366:
      for every 2 <= m <= n, some chosen speed is divisible by m;
* endpoint-protection invoices from THM-357:
      for every target endpoint t, some chosen speed p has ||p t|| < 1/n.

The point is to expose how horizontal x+2 invoices and vertical x*2 debt
couple in the known n=14 and n=16 frontiers.
"""

from __future__ import annotations

from collections import Counter, defaultdict
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
S360 = SourceFileLoader(
    "lonely_runner_endpoint_protection_s360",
    str(ROOT / "04-computation" / "lonely_runner_endpoint_protection_s360.py"),
).load_module()
S362 = SourceFileLoader(
    "lonely_runner_bohr_descent_s362",
    str(ROOT / "04-computation" / "lonely_runner_bohr_descent_s362.py"),
).load_module()


ONE = Fraction(1, 1)


@dataclass(frozen=True)
class SpeedMode:
    value: int
    odd_core: int
    dyadic_height: int
    odd_chain_index: int


@dataclass(frozen=True)
class SieveIPRow:
    label: str
    n: int
    speeds: tuple[int, ...]
    missing: tuple[int, ...]
    covered: tuple[int, ...]
    forced_moduli: tuple[int, ...]
    dyadic_hist: tuple[tuple[int, int], ...]
    odd_core_hist: tuple[tuple[int, int], ...]


@dataclass(frozen=True)
class CoverColumn:
    speed: int
    mask: int
    size: int


@dataclass(frozen=True)
class SetCoverResult:
    label: str
    n: int
    target_count: int
    candidate_count: int
    active_count: int
    union_count: int
    uncovered_count: int
    forced_columns: tuple[int, ...]
    greedy_size: int | None
    greedy_columns: tuple[int, ...]
    exact_size: int | None
    exact_columns: tuple[int, ...]
    search_calls: int
    max_column_size: int
    dyadic_exposure_hist: tuple[tuple[int, int], ...]
    prime_exposure_hist: tuple[tuple[tuple[int, ...], int], ...]


def fmt(value: Fraction | None) -> str:
    return S356.fmt_frac(value)


def circle(value: Fraction) -> Fraction:
    return value % ONE


def dist_to_integer(value: Fraction) -> Fraction:
    residue = circle(value)
    return min(residue, ONE - residue)


def v2(value: int) -> int:
    value = abs(value)
    if value == 0:
        raise ValueError("v2(0) not used in this finite speed ledger")
    out = 0
    while value % 2 == 0:
        out += 1
        value //= 2
    return out


def prime_factors(value: int) -> tuple[int, ...]:
    value = abs(value)
    out: list[int] = []
    prime = 2
    while prime * prime <= value:
        if value % prime == 0:
            out.append(prime)
            while value % prime == 0:
                value //= prime
        prime += 1 if prime == 2 else 2
    if value > 1:
        out.append(value)
    return tuple(out)


def vp(value: int, prime: int) -> int:
    value = abs(value)
    if value == 0:
        raise ValueError("v_p(0) not used in this finite speed ledger")
    out = 0
    while value % prime == 0:
        out += 1
        value //= prime
    return out


def speed_mode(value: int) -> SpeedMode:
    height = v2(value)
    odd = value >> height
    return SpeedMode(
        value=value,
        odd_core=odd,
        dyadic_height=height,
        odd_chain_index=(odd - 1) // 2,
    )


def mode_histograms(speeds: tuple[int, ...]) -> tuple[Counter[int], Counter[int]]:
    dyadic: Counter[int] = Counter()
    odd_core: Counter[int] = Counter()
    for speed in speeds:
        mode = speed_mode(speed)
        dyadic[mode.dyadic_height] += 1
        odd_core[mode.odd_core] += 1
    return dyadic, odd_core


def compact_hist(counter: Counter[int], limit: int = 8) -> str:
    if not counter:
        return "-"
    items = sorted(counter.items())
    pieces = [f"{key}:{count}" for key, count in items[:limit]]
    if len(items) > limit:
        pieces.append("...")
    return " ".join(pieces)


def primitive_tuple(raw: list[int] | tuple[int, ...]) -> tuple[int, ...]:
    return S356.normalize_speed_set(list(raw))


def initial(n: int) -> tuple[int, ...]:
    return tuple(range(1, n))


def drop_add(n: int, dropped: int, added: int) -> tuple[int, ...]:
    speeds = set(initial(n))
    speeds.remove(dropped)
    speeds.add(added)
    return tuple(sorted(speeds))


def n14_gate_ladder() -> tuple[int, ...]:
    # The S380 high-gate ladder, kept deterministic for cross-session diffs.
    return (1, 14, 28, 42, 56, 70, 98, 112, 126, 140, 154, 168, 182)


def n16_eight_ladder() -> tuple[int, ...]:
    return (1, 8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 120)


def small_denominator_missing(n: int, speeds: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(m for m in range(2, n + 1) if all(speed % m for speed in speeds))


def small_denominator_covered(n: int, speeds: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(m for m in range(2, n + 1) if any(speed % m == 0 for speed in speeds))


def sieve_row(label: str, raw_speeds: tuple[int, ...]) -> SieveIPRow:
    speeds = primitive_tuple(raw_speeds)
    n = len(speeds) + 1
    dyadic, odd_core = mode_histograms(speeds)
    # These are rows that currently have exactly one speed paying them.  In an
    # IP branch they are fragile: deleting that speed immediately opens a
    # THM-366 witness.
    forced = tuple(
        m for m in range(2, n + 1) if sum(1 for speed in speeds if speed % m == 0) == 1
    )
    return SieveIPRow(
        label=label,
        n=n,
        speeds=speeds,
        missing=small_denominator_missing(n, speeds),
        covered=small_denominator_covered(n, speeds),
        forced_moduli=forced,
        dyadic_hist=tuple(sorted(dyadic.items())),
        odd_core_hist=tuple(sorted(odd_core.items())),
    )


def unit_points(n: int) -> tuple[Fraction, ...]:
    return tuple(Fraction(a, n) for a in range(1, n) if gcd(a, n) == 1)


def endpoint_values_for_owner(n: int, owner: int) -> tuple[Fraction, ...]:
    values = {
        circle(Fraction(n * m + sign, n * owner))
        for m in range(owner)
        for sign in (-1, 1)
    }
    return tuple(sorted(values))


def endpoint_values_for_owners(n: int, owners: tuple[int, ...]) -> tuple[Fraction, ...]:
    values: set[Fraction] = set()
    for owner in owners:
        values.update(endpoint_values_for_owner(n, owner))
    return tuple(sorted(values))


def protects_at_threshold(n: int, protector: int, endpoint: Fraction) -> bool:
    return dist_to_integer(protector * endpoint) < Fraction(1, n)


def endpoint_extra_dyadic_depth(n: int, endpoint: Fraction) -> int:
    # Extra 2-adic denominator depth beyond the base threshold denominator.
    base = v2(n) if n % 2 == 0 else 0
    return max(0, v2(endpoint.denominator) - base)


def endpoint_extra_prime_depth(n: int, endpoint: Fraction, primes: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(max(0, vp(endpoint.denominator, prime) - vp(n, prime)) for prime in primes)


def build_cover_columns(
    n: int,
    targets: tuple[Fraction, ...],
    candidates: tuple[int, ...],
) -> tuple[CoverColumn, ...]:
    columns: list[CoverColumn] = []
    for speed in candidates:
        mask = 0
        for index, target in enumerate(targets):
            if protects_at_threshold(n, speed, target):
                mask |= 1 << index
        if mask:
            columns.append(CoverColumn(speed=speed, mask=mask, size=mask.bit_count()))
    return tuple(columns)


def greedy_cover(full_mask: int, columns: tuple[CoverColumn, ...]) -> tuple[int, ...] | None:
    covered = 0
    chosen: list[int] = []
    remaining = list(columns)
    while covered != full_mask:
        best = max(
            remaining,
            key=lambda column: ((column.mask & ~covered).bit_count(), -column.speed),
            default=None,
        )
        if best is None or not (best.mask & ~covered):
            return None
        chosen.append(best.speed)
        covered |= best.mask
        remaining = [column for column in remaining if column.speed != best.speed]
    return tuple(chosen)


def exact_set_cover(
    full_mask: int,
    target_count: int,
    columns: tuple[CoverColumn, ...],
    upper_bound: tuple[int, ...] | None,
) -> tuple[int | None, tuple[int, ...], int]:
    """Branch-and-bound exact set cover for the small rows used here."""

    if not columns:
        return None, tuple(), 0

    speed_to_mask = {column.speed: column.mask for column in columns}
    if upper_bound is not None:
        best_size = len(upper_bound)
        best_cover = upper_bound
    else:
        best_size = len(columns) + 1
        best_cover: tuple[int, ...] = tuple()

    candidates_by_endpoint: list[list[int]] = [[] for _ in range(target_count)]
    for column_index, column in enumerate(columns):
        work = column.mask
        while work:
            bit = work & -work
            endpoint_index = bit.bit_length() - 1
            candidates_by_endpoint[endpoint_index].append(column_index)
            work -= bit

    max_column_size = max(column.size for column in columns)
    seen: dict[int, int] = {}
    calls = 0

    def dfs(covered: int, chosen: tuple[int, ...]) -> None:
        nonlocal best_size, best_cover, calls
        calls += 1
        if len(chosen) >= best_size:
            return
        if covered == full_mask:
            best_size = len(chosen)
            best_cover = chosen
            return
        if seen.get(covered, 10**9) <= len(chosen):
            return
        seen[covered] = len(chosen)

        remaining_count = (full_mask & ~covered).bit_count()
        lower_bound = (remaining_count + max_column_size - 1) // max_column_size
        if len(chosen) + lower_bound >= best_size:
            return

        uncovered = full_mask & ~covered
        best_options: list[int] | None = None
        work = uncovered
        while work:
            bit = work & -work
            endpoint_index = bit.bit_length() - 1
            work -= bit
            options = [
                column_index
                for column_index in candidates_by_endpoint[endpoint_index]
                if columns[column_index].mask & ~covered
            ]
            if best_options is None or len(options) < len(best_options):
                best_options = options
                if len(options) <= 1:
                    break

        if not best_options:
            return

        ordered = sorted(
            best_options,
            key=lambda column_index: (
                (columns[column_index].mask & ~covered).bit_count(),
                -columns[column_index].speed,
            ),
            reverse=True,
        )
        for column_index in ordered:
            column = columns[column_index]
            dfs(covered | column.mask, chosen + (column.speed,))

    # Confirm the supplied upper bound really covers the row.
    if upper_bound is not None:
        mask = 0
        for speed in upper_bound:
            mask |= speed_to_mask.get(speed, 0)
        if mask != full_mask:
            best_size = len(columns) + 1
            best_cover = tuple()

    dfs(0, tuple())
    if best_size == len(columns) + 1:
        return None, tuple(), calls
    return best_size, best_cover, calls


def set_cover_result(
    label: str,
    n: int,
    targets: tuple[Fraction, ...],
    candidates: tuple[int, ...],
    solve_exact: bool = True,
) -> SetCoverResult:
    targets = tuple(sorted(set(targets)))
    candidates = tuple(sorted(set(candidates)))
    full_mask = (1 << len(targets)) - 1
    columns = build_cover_columns(n, targets, candidates)
    union_mask = 0
    endpoint_candidate_count: Counter[int] = Counter()
    for endpoint_index in range(len(targets)):
        bit = 1 << endpoint_index
        for column in columns:
            if column.mask & bit:
                endpoint_candidate_count[endpoint_index] += 1
    for column in columns:
        union_mask |= column.mask

    forced = sorted(
        {
            column.speed
            for endpoint_index, count in endpoint_candidate_count.items()
            if count == 1
            for column in columns
            if column.mask & (1 << endpoint_index)
        }
    )

    uncovered_mask = full_mask & ~union_mask
    dyadic_exposure: Counter[int] = Counter()
    primes = prime_factors(n)
    prime_exposure: Counter[tuple[int, ...]] = Counter()
    work = uncovered_mask
    while work:
        bit = work & -work
        endpoint_index = bit.bit_length() - 1
        dyadic_exposure[endpoint_extra_dyadic_depth(n, targets[endpoint_index])] += 1
        prime_exposure[endpoint_extra_prime_depth(n, targets[endpoint_index], primes)] += 1
        work -= bit

    greedy = None
    exact_size = None
    exact_columns: tuple[int, ...] = tuple()
    calls = 0
    if union_mask == full_mask:
        greedy = greedy_cover(full_mask, columns)
        if solve_exact:
            exact_size, exact_columns, calls = exact_set_cover(
                full_mask,
                len(targets),
                columns,
                greedy,
            )

    return SetCoverResult(
        label=label,
        n=n,
        target_count=len(targets),
        candidate_count=len(candidates),
        active_count=len(columns),
        union_count=union_mask.bit_count(),
        uncovered_count=(full_mask & ~union_mask).bit_count(),
        forced_columns=tuple(forced),
        greedy_size=None if greedy is None else len(greedy),
        greedy_columns=tuple() if greedy is None else greedy,
        exact_size=exact_size,
        exact_columns=exact_columns,
        search_calls=calls,
        max_column_size=max((column.size for column in columns), default=0),
        dyadic_exposure_hist=tuple(sorted(dyadic_exposure.items())),
        prime_exposure_hist=tuple(sorted(prime_exposure.items())),
    )


def print_header(title: str) -> None:
    print("=" * 88)
    print(title)
    print("=" * 88)


def print_mode_table() -> None:
    print_header("SPEED COLUMNS AS ODD-CORE / DOUBLING MODES")
    print(
        "A speed column is written v=2^h*odd_core.  h is the vertical x*2 "
        "height; odd_core moves horizontally by x+2."
    )
    rows = [
        ("n14 initial", initial(14)),
        ("n14 drop13 add14", drop_add(14, 13, 14)),
        ("n14 drop13 add182", drop_add(14, 13, lcm(14, 13))),
        ("n14 S380 gate ladder", n14_gate_ladder()),
        ("n16 initial", initial(16)),
        ("n16 drop15 add16", drop_add(16, 15, 16)),
        ("n16 best 8-ladder", n16_eight_ladder()),
    ]
    print("label                  n  dyadic h-counts        odd-core counts")
    print("-" * 88)
    for label, speeds in rows:
        speeds = primitive_tuple(speeds)
        dyadic, odd_core = mode_histograms(speeds)
        print(
            f"{label:<22} {len(speeds)+1:>2}  "
            f"{compact_hist(dyadic):<22} {compact_hist(odd_core)}"
        )
    print()


def print_sieve_ip_table() -> None:
    print_header("SMALL-DENOMINATOR IP INVOICES")
    print(
        "Rows are m=2..n.  A speed column pays row m exactly when m divides it. "
        "Missing rows give immediate THM-366 witnesses."
    )
    rows = [
        sieve_row("n14 initial", initial(14)),
        sieve_row("n14 drop13 add14", drop_add(14, 13, 14)),
        sieve_row("n14 drop13 add182", drop_add(14, 13, lcm(14, 13))),
        sieve_row("n14 S380 gate ladder", n14_gate_ladder()),
        sieve_row("n16 initial", initial(16)),
        sieve_row("n16 drop15 add16", drop_add(16, 15, 16)),
        sieve_row("n16 best 8-ladder", n16_eight_ladder()),
    ]
    print("label                  missing rows       fragile one-payer rows")
    print("-" * 88)
    for row in rows:
        missing = ",".join(map(str, row.missing)) or "-"
        forced = ",".join(map(str, row.forced_moduli)) or "-"
        print(f"{row.label:<22} {missing:<18} {forced}")
    print()


def print_lrc_audit_table() -> None:
    print_header("EXACT LRC AUDITS OF THE SAME COLUMNS")
    print(
        "These are full interval-union endpoint audits, not just the IP "
        "relaxation.  The key pattern is that sieve repairs produce gap or "
        "endpoint debt, not open covers."
    )
    examples = [
        ("n14 initial", initial(14)),
        ("n14 drop13 add14", drop_add(14, 13, 14)),
        ("n14 drop13 add182", drop_add(14, 13, lcm(14, 13))),
        ("n14 S380 gate ladder", n14_gate_ladder()),
        ("n16 initial", initial(16)),
        ("n16 drop15 add16", drop_add(16, 15, 16)),
        ("n16 best 8-ladder", n16_eight_ladder()),
    ]
    print("label                  class          gap/th     unprotected  peel coreE")
    print("-" * 88)
    for label, speeds in examples:
        report = S360.summarize(list(speeds))
        descent = S362.summarize(list(report.speeds))
        gap_ratio = report.max_gap / report.threshold
        print(
            f"{label:<22} {report.classification:<14} "
            f"{float(gap_ratio):>8.6f} {report.unprotected_count:>11} "
            f"{len(descent.peel_layers):>5} {descent.core_endpoint_count:>5}"
        )
    print()


def print_cover_result(row: SetCoverResult) -> None:
    greedy = "-" if row.greedy_size is None else f"{row.greedy_size}:{row.greedy_columns}"
    exact = "-" if row.exact_size is None else f"{row.exact_size}:{row.exact_columns}"
    forced = ",".join(map(str, row.forced_columns[:16])) or "-"
    if len(row.forced_columns) > 16:
        forced += ",..."
    exposure = " ".join(f"h+{h}:{count}" for h, count in row.dyadic_exposure_hist) or "-"
    primes = prime_factors(row.n)
    prime_exposure = (
        " ".join(
            "{"
            + ",".join(f"{prime}:+{depth}" for prime, depth in zip(primes, vector))
            + f"}}:{count}"
            for vector, count in row.prime_exposure_hist
        )
        or "-"
    )
    print(f"[{row.label}]")
    print(
        f"  n={row.n} targets={row.target_count} candidates={row.candidate_count} "
        f"active={row.active_count} max_col={row.max_column_size}"
    )
    print(
        f"  union={row.union_count}/{row.target_count} "
        f"uncovered={row.uncovered_count} uncovered_dyadic_depths={exposure}"
    )
    print(f"  uncovered_prime_depths={prime_exposure}")
    print(f"  forced_columns={forced}")
    print(f"  greedy={greedy}")
    print(f"  exact={exact} search_calls={row.search_calls}")
    print()


def print_endpoint_cover_ip_table() -> None:
    print_header("ENDPOINT-PROTECTION SET-COVER SUBPROBLEMS")
    print(
        "Each row is an IP relaxation: choose speed columns so every listed "
        "endpoint has at least one strict protector.  The target endpoints "
        "are fixed first, so this is a certificate microscope rather than "
        "the full nonlinear counterexample search."
    )

    rows = [
        set_cover_result(
            "n14 unit layer, candidates 1..28",
            14,
            unit_points(14),
            tuple(range(1, 29)),
        ),
        set_cover_result(
            "n14 owner 14 endpoints, lower 1..13",
            14,
            endpoint_values_for_owner(14, 14),
            tuple(range(1, 14)),
        ),
        set_cover_result(
            "n14 owner 14 endpoints, 1..28",
            14,
            endpoint_values_for_owner(14, 14),
            tuple(range(1, 29)),
        ),
        set_cover_result(
            "n14 S380 ladder endpoints by same ladder",
            14,
            endpoint_values_for_owners(14, n14_gate_ladder()),
            n14_gate_ladder(),
            solve_exact=False,
        ),
        set_cover_result(
            "n16 unit layer, candidates 1..32",
            16,
            unit_points(16),
            tuple(range(1, 33)),
        ),
        set_cover_result(
            "n16 owner 16 endpoints, lower 1..15",
            16,
            endpoint_values_for_owner(16, 16),
            tuple(range(1, 16)),
        ),
        set_cover_result(
            "n16 best 8-ladder endpoints by same ladder",
            16,
            endpoint_values_for_owners(16, n16_eight_ladder()),
            n16_eight_ladder(),
            solve_exact=False,
        ),
    ]
    for row in rows:
        print_cover_result(row)


def print_dual_target() -> None:
    print_header("TWO-MODE IP PROOF TARGET")
    print(
        "Primal variables: x_v in {0,1}, with v=2^h*odd_core.  Horizontal "
        "constraints are small-denominator and odd-core invoices; vertical "
        "constraints are endpoint debt at dyadic heights."
    )
    print()
    print("Counterexample feasibility would require, at minimum:")
    print("  (1)  sum_{m|v} x_v >= 1 for every 2 <= m <= n.")
    print("  (2)  sum_{p: ||p*t|| < 1/n} x_p >= 1 for every active endpoint t.")
    print("  (3)  exactly n-1 primitive speed columns, with gcd(v: x_v=1)=1.")
    print()
    print(
        "The useful dual should weight rows, not speeds.  Unit rows force an "
        "n-gate; the n-gate exports endpoint rows to deeper denominator "
        "layers.  For n=16 the vertical invoice is already sharp: owner 16 "
        "requires nine lower columns (1,3,5,7,8,9,11,13,15).  For n=14 the "
        "horizontal sieve can be paid by 14-gates, but the same-gate ladder "
        "still leaves endpoint rows uncovered at exported 2- and 7-adic "
        "depths."
    )
    print()


def main() -> None:
    print("LRC integer-programming mode lens (codex-2026-05-31 S420)")
    print(
        "The calculation treats natural numbers as odd-core x+2 chains with "
        "vertical x*2 doublings, then builds exact LRC row/column ledgers.\n"
    )
    print_mode_table()
    print_sieve_ip_table()
    print_lrc_audit_table()
    print_endpoint_cover_ip_table()
    print_dual_target()


if __name__ == "__main__":
    main()
