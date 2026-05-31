#!/usr/bin/env python3
"""
lrc_n14_15_16_gate_frontier_s21.py

oracle-2026-05-31 S21

Focused proof-search pass for Lonely Runner at n=14,15,16.

The shared proof pressure is:

1. THM-369 forces every coarse denominator q <= n to be paid by a speed
   divisible by q.  In particular a counterexample must contain an n-gate.
2. The n-gate then owns 2n forbidden endpoints.  These endpoints can be covered
   locally by lower or near-window speeds, but cheap local covers leave coarse
   divisor invoices unpaid or create an imprimitive residue pattern.
3. Paying the unpaid invoices in simple deterministic completions reopens a
   positive Archimedean gap.

This script makes those three statements exact for n=14,15,16.  It is not a
proof; it records a small finite obstruction ledger that future proofs can try
to turn into a Hall/dual certificate.
"""

from __future__ import annotations

from collections import Counter
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

ONE = Fraction(1, 1)


@dataclass(frozen=True)
class CoverAudit:
    n: int
    label: str
    candidate_count: int
    active_count: int
    endpoint_count: int
    union_count: int
    exact_size: int | None
    exact_columns: tuple[int, ...]
    private_count: int
    private_sample: tuple[tuple[Fraction, int], ...]
    gate_cover_gcd: int
    gate_cover_missing_rows: tuple[int, ...]


@dataclass(frozen=True)
class LadderAudit:
    n: int
    scale: int
    skip: int
    gap_ratio: Fraction
    debt: int
    product: Fraction
    layer_hist: tuple[tuple[int, int], ...]


@dataclass(frozen=True)
class CompletionAudit:
    n: int
    label: str
    base_columns: tuple[int, ...]
    completed_speeds: tuple[int, ...]
    missing_rows: tuple[int, ...]
    classification: str
    gap_ratio: Fraction
    debt: int
    product: Fraction


@dataclass(frozen=True)
class FamilyAudit:
    n: int
    label: str
    exact_size: int
    family_count: int
    gcd_hist: tuple[tuple[int, int], ...]
    min_missing_rows: int
    row_free_count: int
    row_free_gcd_hist: tuple[tuple[int, int], ...]
    example: tuple[int, ...]
    example_missing: tuple[int, ...]


def fmt_frac(value: Fraction | None) -> str:
    return S356.fmt_frac(value)


def circle(value: Fraction) -> Fraction:
    return value % ONE


def distance_to_integer(value: Fraction) -> Fraction:
    residue = circle(value)
    return min(residue, ONE - residue)


def primitive_gcd(values: tuple[int, ...] | list[int] | set[int]) -> int:
    out = 0
    for value in values:
        out = gcd(out, value)
    return out


def missing_rows(n: int, speeds: tuple[int, ...] | set[int]) -> tuple[int, ...]:
    return tuple(q for q in range(2, n + 1) if all(speed % q for speed in speeds))


def owner_endpoints(n: int, owner: int) -> tuple[Fraction, ...]:
    values = {
        circle(Fraction(n * m + sign, n * owner))
        for m in range(owner)
        for sign in (-1, 1)
    }
    return tuple(sorted(values))


def protects(n: int, protector: int, endpoint: Fraction) -> bool:
    return distance_to_integer(protector * endpoint) < Fraction(1, n)


def columns_for_endpoints(
    n: int, endpoints: tuple[Fraction, ...], candidates: tuple[int, ...]
) -> tuple[tuple[int, int], ...]:
    cols: list[tuple[int, int]] = []
    for candidate in candidates:
        mask = 0
        for idx, endpoint in enumerate(endpoints):
            if protects(n, candidate, endpoint):
                mask |= 1 << idx
        if mask:
            cols.append((candidate, mask))
    return tuple(cols)


def exact_cover(
    endpoint_count: int, columns: tuple[tuple[int, int], ...]
) -> tuple[int | None, tuple[int, ...]]:
    full = (1 << endpoint_count) - 1
    union = 0
    for _candidate, mask in columns:
        union |= mask
    if union != full:
        return None, ()

    by_endpoint: list[list[int]] = [[] for _ in range(endpoint_count)]
    for idx, (_candidate, mask) in enumerate(columns):
        bitset = mask
        while bitset:
            low = bitset & -bitset
            endpoint_idx = low.bit_length() - 1
            by_endpoint[endpoint_idx].append(idx)
            bitset -= low

    suffix_union = [0] * (len(columns) + 1)
    # Greedy ordering: large columns first.  Rebuild incidence after sorting.
    sorted_cols = tuple(sorted(columns, key=lambda col: (-col[1].bit_count(), col[0])))
    by_endpoint = [[] for _ in range(endpoint_count)]
    for idx, (_candidate, mask) in enumerate(sorted_cols):
        bitset = mask
        while bitset:
            low = bitset & -bitset
            endpoint_idx = low.bit_length() - 1
            by_endpoint[endpoint_idx].append(idx)
            bitset -= low
    for idx in range(len(sorted_cols) - 1, -1, -1):
        suffix_union[idx] = suffix_union[idx + 1] | sorted_cols[idx][1]

    best_size = len(sorted_cols) + 1
    best_columns: tuple[int, ...] = ()

    def dfs(covered: int, chosen: tuple[int, ...]) -> None:
        nonlocal best_size, best_columns
        if len(chosen) >= best_size:
            return
        if covered == full:
            best_size = len(chosen)
            best_columns = tuple(sorted(sorted_cols[idx][0] for idx in chosen))
            return
        if (covered | suffix_union[0]) != full:
            return
        uncovered = full ^ covered
        # Branch on the uncovered endpoint with the fewest available columns.
        best_endpoint = None
        best_options: list[int] | None = None
        bitset = uncovered
        while bitset:
            low = bitset & -bitset
            endpoint_idx = low.bit_length() - 1
            options = [
                idx
                for idx in by_endpoint[endpoint_idx]
                if sorted_cols[idx][1] & uncovered
            ]
            if best_options is None or len(options) < len(best_options):
                best_endpoint = endpoint_idx
                best_options = options
            bitset -= low
        if best_endpoint is None or not best_options:
            return
        best_options.sort(key=lambda idx: (-((sorted_cols[idx][1] & uncovered).bit_count()), sorted_cols[idx][0]))
        for idx in best_options:
            dfs(covered | sorted_cols[idx][1], chosen + (idx,))

    dfs(0, ())
    if best_size == len(sorted_cols) + 1:
        return None, ()
    return best_size, best_columns


def cover_audit(n: int, label: str, candidates: tuple[int, ...]) -> CoverAudit:
    endpoints = owner_endpoints(n, n)
    columns = columns_for_endpoints(n, endpoints, candidates)
    union = 0
    for _candidate, mask in columns:
        union |= mask
    exact_size, exact_columns = exact_cover(len(endpoints), columns)
    private: list[tuple[Fraction, int]] = []
    for idx, endpoint in enumerate(endpoints):
        owners = [candidate for candidate, mask in columns if (mask >> idx) & 1]
        if len(owners) == 1:
            private.append((endpoint, owners[0]))
    gate_cover = tuple(sorted(set(exact_columns) | {n}))
    return CoverAudit(
        n=n,
        label=label,
        candidate_count=len(candidates),
        active_count=len(columns),
        endpoint_count=len(endpoints),
        union_count=union.bit_count(),
        exact_size=exact_size,
        exact_columns=exact_columns,
        private_count=len(private),
        private_sample=tuple(private[:6]),
        gate_cover_gcd=primitive_gcd(gate_cover),
        gate_cover_missing_rows=missing_rows(n, set(gate_cover)),
    )


def summarize_speeds(speeds: tuple[int, ...]) -> tuple[str, Fraction, int, Fraction]:
    summary = S360.summarize(list(speeds))
    gap_ratio = summary.max_gap / summary.threshold
    return summary.classification, gap_ratio, summary.unprotected_count, gap_ratio * summary.unprotected_count


def unprotected_points(speeds: tuple[int, ...]) -> tuple[Fraction, ...]:
    endpoints = {endpoint.value for endpoint in S360.endpoints(speeds)}
    return tuple(
        point
        for point in sorted(endpoints)
        if not any(S360.direct_protects(speeds, speed, point) for speed in speeds)
    )


def debt_layer(n: int, point: Fraction) -> int:
    return point.denominator // gcd(point.denominator, n)


def layer_hist(n: int, speeds: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    return tuple(sorted(Counter(debt_layer(n, point) for point in unprotected_points(speeds)).items()))


def quotient_ladders(n: int) -> tuple[LadderAudit, ...]:
    scales = tuple(d for d in range(2, n + 1) if n % d == 0)
    out: list[LadderAudit] = []
    for scale in scales:
        best: tuple[tuple[Fraction, int, int], tuple[int, ...], tuple[str, Fraction, int, Fraction]] | None = None
        for skip in range(1, n):
            speeds = tuple(sorted({1} | {scale * q for q in range(1, n) if q != skip}))
            if len(speeds) != n - 1 or primitive_gcd(speeds) != 1:
                continue
            summary = summarize_speeds(speeds)
            key = (summary[1], summary[2], skip)
            if best is None or key < best[0]:
                best = (key, speeds, summary)
        if best is None:
            continue
        key, speeds, summary = best
        _classification, gap_ratio, debt, product = summary
        out.append(
            LadderAudit(
                n=n,
                scale=scale,
                skip=key[2],
                gap_ratio=gap_ratio,
                debt=debt,
                product=product,
                layer_hist=layer_hist(n, speeds),
            )
        )
    return tuple(out)


def greedy_complete(n: int, base_columns: tuple[int, ...], limit: int = 240) -> tuple[int, ...]:
    speeds = set(base_columns)
    while len(speeds) < n - 1 and missing_rows(n, speeds):
        misses = missing_rows(n, speeds)
        best: tuple[tuple[int, int], int] | None = None
        for candidate in range(1, limit + 1):
            if candidate in speeds:
                continue
            covered = sum(1 for row in misses if candidate % row == 0)
            if covered == 0:
                continue
            key = (-covered, candidate)
            if best is None or key < best[0]:
                best = (key, candidate)
        if best is None:
            break
        speeds.add(best[1])
    for candidate in range(1, limit + 1):
        if len(speeds) >= n - 1:
            break
        if candidate not in speeds:
            speeds.add(candidate)
    return tuple(sorted(speeds))


def completion_audit(n: int, cover: CoverAudit) -> CompletionAudit:
    base = tuple(sorted(set(cover.exact_columns) | {n}))
    completed = greedy_complete(n, base)
    classification, gap_ratio, debt, product = summarize_speeds(completed)
    return CompletionAudit(
        n=n,
        label=cover.label,
        base_columns=base,
        completed_speeds=completed,
        missing_rows=missing_rows(n, completed),
        classification=classification,
        gap_ratio=gap_ratio,
        debt=debt,
        product=product,
    )


def enumerate_minimum_family(n: int, label: str, candidates: tuple[int, ...], exact_size: int) -> FamilyAudit:
    endpoints = owner_endpoints(n, n)
    columns = columns_for_endpoints(n, endpoints, candidates)
    full = (1 << len(endpoints)) - 1
    from itertools import combinations

    family_count = 0
    gcd_hist: Counter[int] = Counter()
    row_free_gcd_hist: Counter[int] = Counter()
    row_free_count = 0
    min_missing = n
    example: tuple[int, ...] = ()
    example_missing: tuple[int, ...] = ()
    for combo in combinations(columns, exact_size):
        union = 0
        for _candidate, mask in combo:
            union |= mask
        if union != full:
            continue
        speeds = tuple(sorted({n} | {candidate for candidate, _mask in combo}))
        misses = missing_rows(n, set(speeds))
        g = primitive_gcd(speeds)
        family_count += 1
        gcd_hist[g] += 1
        if len(misses) < min_missing:
            min_missing = len(misses)
            example = speeds
            example_missing = misses
        if not misses:
            row_free_count += 1
            row_free_gcd_hist[g] += 1
    return FamilyAudit(
        n=n,
        label=label,
        exact_size=exact_size,
        family_count=family_count,
        gcd_hist=tuple(sorted(gcd_hist.items())),
        min_missing_rows=min_missing,
        row_free_count=row_free_count,
        row_free_gcd_hist=tuple(sorted(row_free_gcd_hist.items())),
        example=example,
        example_missing=example_missing,
    )


def print_cover_audits(audits: tuple[CoverAudit, ...]) -> None:
    print("A. LOCAL n-GATE ENDPOINT INVOICES")
    print("=" * 92)
    print("Exact set cover for the 2n endpoints owned by a single n-gate.")
    print()
    print(
        "  n  candidates       active endpoints union exact columns"
        "                         private gcd missing_rows"
    )
    for audit in audits:
        exact = "-" if audit.exact_size is None else str(audit.exact_size)
        columns = ",".join(str(value) for value in audit.exact_columns) or "-"
        misses = ",".join(str(value) for value in audit.gate_cover_missing_rows) or "-"
        print(
            f"  {audit.n:<2} {audit.label:<14} {audit.active_count:>6}/"
            f"{audit.candidate_count:<3} {audit.endpoint_count:>9} "
            f"{audit.union_count:>5} {exact:>5} {columns:<31} "
            f"{audit.private_count:>7} {audit.gate_cover_gcd:>3} {misses}"
        )
    print()
    print("  Private samples for lower covers:")
    for audit in audits:
        if audit.label != "lower":
            continue
        sample = ", ".join(f"{fmt_frac(point)}->{owner}" for point, owner in audit.private_sample)
        print(f"    n={audit.n}: {sample}")
    print()


def print_family_audits(audits: tuple[FamilyAudit, ...]) -> None:
    print("A2. ALL MINIMUM WINDOW-COVER FAMILIES")
    print("=" * 92)
    print("Enumerate every minimum cover using candidates 1..2n-1, then add the n-gate.")
    print()
    print("  n  size count gcd_hist       min_missing row_free row_free_gcd example")
    for audit in audits:
        gcd_hist = ",".join(f"{g}:{count}" for g, count in audit.gcd_hist) or "-"
        row_gcd = ",".join(f"{g}:{count}" for g, count in audit.row_free_gcd_hist) or "-"
        example = ",".join(str(value) for value in audit.example)
        misses = ",".join(str(value) for value in audit.example_missing) or "-"
        print(
            f"  {audit.n:<2} {audit.exact_size:>4} {audit.family_count:>5} "
            f"{gcd_hist:<14} {audit.min_missing_rows:>11} "
            f"{audit.row_free_count:>8} {row_gcd:<12} {example}  miss={misses}"
        )
    print()


def print_ladders(rows: tuple[LadderAudit, ...]) -> None:
    print("B. BEST QUOTIENT-LADDER EXPORTS")
    print("=" * 92)
    print("For each proper divisor scale, scan skip positions and keep the smallest gap.")
    print()
    print("  n  scale skip gap/th       debt product      debt_layers")
    for row in rows:
        hist = ",".join(f"{layer}:{count}" for layer, count in row.layer_hist[:5])
        if len(row.layer_hist) > 5:
            hist += ",..."
        print(
            f"  {row.n:<2} {row.scale:>5} {row.skip:>4} "
            f"{fmt_frac(row.gap_ratio):>10} {row.debt:>6} "
            f"{fmt_frac(row.product):>12} {hist}"
        )
    print()


def print_completions(rows: tuple[CompletionAudit, ...]) -> None:
    print("C. TWO-STAGE COMPLETION ATTEMPTS")
    print("=" * 92)
    print("Start with an exact local gate cover, then greedily pay all coarse rows q<=n.")
    print()
    print("  n  base    class        gap/th       debt product      missing completed_speeds")
    for row in rows:
        misses = ",".join(str(value) for value in row.missing_rows) or "-"
        speeds = ",".join(str(value) for value in row.completed_speeds)
        print(
            f"  {row.n:<2} {row.label:<7} {row.classification:<12} "
            f"{fmt_frac(row.gap_ratio):>10} {row.debt:>6} "
            f"{fmt_frac(row.product):>12} {misses:<7} {speeds}"
        )
    print()


def main() -> None:
    ns = (14, 15, 16)
    audits: list[CoverAudit] = []
    for n in ns:
        audits.append(cover_audit(n, "lower", tuple(range(1, n))))
        audits.append(cover_audit(n, "window", tuple(range(1, 2 * n))))
    family_rows = [
        enumerate_minimum_family(audit.n, audit.label, tuple(range(1, 2 * audit.n)), audit.exact_size)
        for audit in audits
        if audit.label == "window" and audit.exact_size is not None
    ]

    ladder_rows: list[LadderAudit] = []
    for n in ns:
        ladder_rows.extend(quotient_ladders(n))

    completion_rows = [completion_audit(audit.n, audit) for audit in audits]

    print("LRC n=14,15,16 gate-frontier attempts (oracle-2026-05-31-S21)")
    print()
    print_cover_audits(tuple(audits))
    print_family_audits(tuple(family_rows))
    print_ladders(tuple(ladder_rows))
    print_completions(tuple(completion_rows))

    print("D. SYNTHESIS")
    print("=" * 92)
    print("  1. Local endpoint repair is not the same as a counterexample.")
    print("     The cheapest window covers leave divisor rows unpaid; paying them")
    print("     in deterministic completions gives positive-gap systems in all")
    print("     three cases.")
    print("  2. The minimum-window family audit isolates the obstruction:")
    print("     n=14 can pay all coarse rows at minimum local size only in gcd-2")
    print("     even covers; n=15 minimum local covers always leave at least two")
    print("     coarse rows unpaid; n=16 has a unique minimum local cover, also")
    print("     gcd-2 and still missing rows 7,12,14.")
    print("  3. n=15 behaves like a mixed odd-prime product building, not like the")
    print("     scalar n=14 moat or the pure dyadic n=16 ladder.  Its best quotient")
    print("     ladder product here is the 3-branch value 7/11, while the 5-branch")
    print("     has smaller visible gap but larger endpoint debt.")
    print("  4. A plausible proof target is now a Hall/dual certificate: local")
    print("     gate-endpoint rows plus coarse divisor rows cannot be paid by")
    print("     n-1 primitive columns without leaving either positive gap or")
    print("     unprotected endpoint debt.")


if __name__ == "__main__":
    main()
