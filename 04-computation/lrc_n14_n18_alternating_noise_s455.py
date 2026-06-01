#!/usr/bin/env python3
"""
lrc_n14_n18_alternating_noise_s455.py

codex-2026-06-01 S455

Alternate between the n=14 and n=18 Lonely Runner frontiers.

This pass deliberately uses "noise cards" from nearby LRC literature and repo
threads as hypothesis triggers, but every trigger is forced back into exact
finite computations:

1. local n-gate lower-cover invoices;
2. row-parent and n-gate ladder debt ledgers;
3. one-slot repair scans with and without a multiple of n;
4. Tournament Analysis switchboards at the gate time t=1/n.

The goal is not to prove LRC here.  It is to sharpen proof-shaped invariants
for the first-even seam 7 -> 14 and the square-core seam 9 -> 18.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from itertools import combinations
from math import cos, gcd, pi, prod, sin, sqrt
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
S23 = SourceFileLoader(
    "tournament_analysis_metric_lifts_s23",
    str(ROOT / "04-computation" / "tournament_analysis_metric_lifts_s23.py"),
).load_module()
S356 = SourceFileLoader(
    "lonely_runner_residue_probe_s356",
    str(ROOT / "04-computation" / "lonely_runner_residue_probe_s356.py"),
).load_module()
S360 = SourceFileLoader(
    "lonely_runner_endpoint_protection_s360",
    str(ROOT / "04-computation" / "lonely_runner_endpoint_protection_s360.py"),
).load_module()
S411 = SourceFileLoader(
    "lrc_column_row_modes_s411",
    str(ROOT / "04-computation" / "lrc_column_row_modes_s411.py"),
).load_module()
S420 = SourceFileLoader(
    "lrc_integer_programming_modes_s420",
    str(ROOT / "04-computation" / "lrc_integer_programming_modes_s420.py"),
).load_module()


ONE = Fraction(1, 1)
EPS = 1.0e-12


@dataclass(frozen=True)
class GateAudit:
    n: int
    target_count: int
    exact_size: int | None
    forced: tuple[int, ...]
    covers: tuple[tuple[int, ...], ...]
    free_parts: tuple[tuple[int, ...], ...]
    private_counts: tuple[tuple[int, str, int, int], ...]


@dataclass(frozen=True)
class RowAudit:
    n: int
    label: str
    speeds: tuple[int, ...]
    classification: str
    gap_ratio: Fraction
    unprotected: int
    product: Fraction
    gates_n: int
    parent_multiples: int
    missing: tuple[int, ...]
    fragile: tuple[int, ...]
    depth_hist: tuple[tuple[tuple[int, ...], int], ...]
    frontier_mass: Fraction
    denominator_pressure: int
    top_owners: tuple[tuple[int, int, int], ...]


@dataclass(frozen=True)
class RepairScan:
    n: int
    drop: int
    multiple_added: int
    multiple_gap: Fraction
    multiple_debt: int
    nonmultiple_added: int
    nonmultiple_gap: Fraction
    nonmultiple_debt: int


@dataclass(frozen=True)
class SwitchAudit:
    n: int
    label: str
    t: Fraction
    phase_cycles: int
    safe_cycles: int
    relspeed_cycles: int
    local_gap_cycles: int
    safe_scc: int
    safe_score_hist: tuple[tuple[int, int], ...]


def fmt(value: Fraction | None) -> str:
    return S356.fmt_frac(value)


def fmt_float(value: Fraction) -> str:
    return f"{float(value):.6f}"


def header(title: str) -> None:
    print("=" * 118)
    print(title)
    print("=" * 118)


def circle(value: Fraction) -> Fraction:
    return value % ONE


def prime_factors(value: int) -> tuple[int, ...]:
    out: list[int] = []
    trial = 2
    while trial * trial <= value:
        if value % trial == 0:
            out.append(trial)
            while value % trial == 0:
                value //= trial
        trial += 1 if trial == 2 else 2
    if value > 1:
        out.append(value)
    return tuple(out)


def depth_scale(primes: tuple[int, ...], depth: tuple[int, ...]) -> int:
    return prod(prime**height for prime, height in zip(primes, depth))


def initial(n: int) -> tuple[int, ...]:
    return tuple(range(1, n))


def lpd_ladder(n: int) -> tuple[int, ...]:
    row = S411.best_lpd_ladder(n)
    if row.lpd is None or row.skip is None:
        raise ValueError(f"n={n} has no largest-proper-divisor ladder")
    return S411.ladder(n, row.lpd, row.skip)


def n_gate_ladder(n: int) -> tuple[int, ...]:
    row = S411.best_lpd_ladder(n)
    if row.skip is None:
        raise ValueError(f"n={n} has no row-parent skip")
    speeds = tuple(sorted({1} | {n * q for q in range(1, n) if q != row.skip}))
    if len(speeds) != n - 1 or gcd(*speeds) != 1:
        raise ValueError(f"bad n-gate ladder for n={n}")
    return speeds


def depth_text(primes: tuple[int, ...], depth: tuple[int, ...]) -> str:
    return "{" + ",".join(f"{prime}:+{height}" for prime, height in zip(primes, depth)) + "}"


def depth_hist_text(n: int, hist: tuple[tuple[tuple[int, ...], int], ...]) -> str:
    if not hist:
        return "-"
    primes = prime_factors(n)
    return " ".join(f"{depth_text(primes, depth)}:{count}" for depth, count in hist)


def unprotected_values(speeds: tuple[int, ...]) -> set[Fraction]:
    endpoints = {endpoint.value for endpoint in S360.endpoints(speeds)}
    return {
        value
        for value in endpoints
        if not any(S360.direct_protects(speeds, speed, value) for speed in speeds)
    }


def top_owner_debt(speeds: tuple[int, ...], limit: int = 5) -> tuple[tuple[int, int, int], ...]:
    bad = unprotected_values(speeds)
    by_owner_labels: Counter[int] = Counter()
    by_owner_unique: dict[int, set[Fraction]] = defaultdict(set)
    for endpoint in S360.endpoints(speeds):
        if endpoint.value in bad:
            by_owner_labels[endpoint.speed] += 1
            by_owner_unique[endpoint.speed].add(endpoint.value)
    return tuple(
        (owner, labels, len(by_owner_unique[owner]))
        for owner, labels in by_owner_labels.most_common(limit)
    )


def moduli_ledger(n: int, speeds: tuple[int, ...]) -> tuple[tuple[int, ...], tuple[int, ...]]:
    counts = {
        modulus: sum(1 for speed in speeds if speed % modulus == 0)
        for modulus in range(2, n + 1)
    }
    missing = tuple(modulus for modulus, count in counts.items() if count == 0)
    fragile = tuple(modulus for modulus, count in counts.items() if count == 1)
    return missing, fragile


def row_audit(n: int, label: str, raw_speeds: tuple[int, ...]) -> RowAudit:
    speeds = S356.normalize_speed_set(list(raw_speeds))
    summary = S360.summarize(list(speeds))
    bad = unprotected_values(speeds)
    primes = prime_factors(n)
    depth_hist = Counter(S420.endpoint_extra_prime_depth(n, point, primes) for point in bad)
    frontier_mass = sum(Fraction(1, depth_scale(primes, depth)) for depth in depth_hist.elements())
    pressure = sum(depth_scale(primes, depth) for depth in depth_hist.elements())
    missing, fragile = moduli_ledger(n, speeds)
    return RowAudit(
        n=n,
        label=label,
        speeds=speeds,
        classification=summary.classification,
        gap_ratio=summary.max_gap / summary.threshold,
        unprotected=summary.unprotected_count,
        product=(summary.max_gap / summary.threshold) * summary.unprotected_count,
        gates_n=sum(1 for speed in speeds if speed % n == 0),
        parent_multiples=sum(1 for speed in speeds if speed % (n // 2) == 0),
        missing=missing,
        fragile=fragile,
        depth_hist=tuple(sorted(depth_hist.items())),
        frontier_mass=frontier_mass,
        denominator_pressure=pressure,
        top_owners=top_owner_debt(speeds),
    )


def exact_gate_audit(n: int) -> GateAudit:
    targets = S420.endpoint_values_for_owner(n, n)
    candidates = tuple(range(1, n))
    columns = S420.build_cover_columns(n, targets, candidates)
    result = S420.set_cover_result(f"owner {n} lower cover", n, targets, candidates)
    speed_to_mask = {column.speed: column.mask for column in columns}
    full_mask = (1 << len(targets)) - 1
    covers: list[tuple[int, ...]] = []
    if result.exact_size is not None:
        for combo in combinations(candidates, result.exact_size):
            mask = 0
            for speed in combo:
                mask |= speed_to_mask.get(speed, 0)
            if mask == full_mask:
                covers.append(combo)

    union_without: dict[int, int] = {}
    for column in columns:
        mask = 0
        for other in columns:
            if other.speed != column.speed:
                mask |= other.mask
        union_without[column.speed] = mask

    private_counts: list[tuple[int, str, int, int]] = []
    for column in columns:
        private = (column.mask & ~union_without[column.speed]).bit_count()
        if private or column.speed in result.forced_columns:
            mode = S420.speed_mode(column.speed)
            private_counts.append(
                (
                    column.speed,
                    f"2^{mode.dyadic_height}*{mode.odd_core}",
                    column.size,
                    private,
                )
            )

    forced = result.forced_columns
    free_parts = tuple(tuple(speed for speed in cover if speed not in forced) for cover in covers)
    return GateAudit(
        n=n,
        target_count=len(targets),
        exact_size=result.exact_size,
        forced=forced,
        covers=tuple(covers),
        free_parts=free_parts,
        private_counts=tuple(private_counts),
    )


def one_slot_repair_scan(n: int) -> RepairScan:
    """Drop n-1 from the initial segment and scan one replacement up to n(n-1)."""

    drop = n - 1
    base = tuple(range(1, drop))
    best_multiple: tuple[Fraction, int, int] | None = None
    best_nonmultiple: tuple[Fraction, int, int] | None = None
    for added in range(n, n * (n - 1) + 1):
        speeds = base + (added,)
        if gcd(*speeds) != 1:
            continue
        summary = S360.summarize(list(speeds))
        item = (summary.max_gap / summary.threshold, added, summary.unprotected_count)
        if added % n == 0:
            if best_multiple is None or item < best_multiple:
                best_multiple = item
        else:
            if best_nonmultiple is None or item < best_nonmultiple:
                best_nonmultiple = item
    if best_multiple is None or best_nonmultiple is None:
        raise ValueError(f"repair scan failed for n={n}")
    return RepairScan(
        n=n,
        drop=drop,
        multiple_added=best_multiple[1],
        multiple_gap=best_multiple[0],
        multiple_debt=best_multiple[2],
        nonmultiple_added=best_nonmultiple[1],
        nonmultiple_gap=best_nonmultiple[0],
        nonmultiple_debt=best_nonmultiple[2],
    )


def circ_delta(a: float, b: float) -> float:
    return (b - a) % 1.0


def circ_dist(a: float, b: float) -> float:
    delta = circ_delta(a, b)
    return min(delta, 1.0 - delta)


def chord_dist(a: float, b: float) -> float:
    return sqrt(max(0.0, 2.0 - 2.0 * cos(2.0 * pi * circ_delta(a, b))))


def orient_switch(n: int, switch_value) -> list[list[int]]:
    def beats(i: int, j: int) -> bool:
        value = switch_value(i, j)
        if abs(value) <= EPS:
            return S23.tie_beats(i, j)
        return S23.tie_beats(i, j) if value > 0 else S23.tie_beats(j, i)

    return S23.orient_by_edge_rule(n, beats)


def local_gap_rank(positions: list[float]) -> list[list[int]]:
    n = len(positions)
    order = sorted(range(n), key=lambda index: (positions[index], index))
    scores = [0.0] * n
    for idx, vertex in enumerate(order):
        prev_vertex = order[(idx - 1) % n]
        next_vertex = order[(idx + 1) % n]
        left = circ_delta(positions[prev_vertex], positions[vertex])
        right = circ_delta(positions[vertex], positions[next_vertex])
        scores[vertex] = min(left, right)
    return S23.orient_by_score(scores)


def switch_audit(n: int, label: str, moving_speeds: tuple[int, ...], t: Fraction) -> SwitchAudit:
    speeds = (0,) + tuple(moving_speeds)
    positions = [float((speed * t) % 1) for speed in speeds]
    total = len(speeds)
    phase = S23.stats(S23.phase_halfturn_tournament(positions), compute_h=False)
    safe = S23.stats(
        orient_switch(total, lambda i, j, p=positions: circ_dist(p[i], p[j]) - 1.0 / n),
        compute_h=False,
    )
    relspeed = S23.stats(
        orient_switch(
            total,
            lambda i, j, p=positions, s=speeds: (
                circ_dist(p[i], p[j]) / (1 + abs(s[i] - s[j])) - 1.0 / (2 * n)
            ),
        ),
        compute_h=False,
    )
    local_gap = S23.stats(local_gap_rank(positions), compute_h=False)
    return SwitchAudit(
        n=n,
        label=label,
        t=t,
        phase_cycles=phase.cyclic_triples,
        safe_cycles=safe.cyclic_triples,
        relspeed_cycles=relspeed.cyclic_triples,
        local_gap_cycles=local_gap.cyclic_triples,
        safe_scc=safe.scc_count,
        safe_score_hist=safe.score_hist,
    )


def print_gate_audits(audits: tuple[GateAudit, ...]) -> None:
    header("ROUND 1: LOCAL N-GATE LOWER-COVER INVOICES")
    print("A counterexample branch that spends an n-gate must also pay this local lower-cover invoice.")
    print()
    for audit in audits:
        print(f"n={audit.n}: targets={audit.target_count} exact_size={audit.exact_size} cover_count={len(audit.covers)}")
        print(f"  forced={audit.forced}")
        print(f"  free_parts={audit.free_parts}")
        print("  private endpoint rows by lower column:")
        for speed, mode, covers, private in audit.private_counts:
            print(f"    {speed:>2}  {mode:<7} covers={covers:>2} private={private:>2}")
        print()


def print_row_table(rows: tuple[RowAudit, ...]) -> None:
    header("ROUND 2: ROW-PARENT AND N-GATE DEBT LEDGERS")
    print("product = (max_gap / threshold) * unprotected_endpoint_count.")
    print()
    print(
        "n  label                 class          gap/th      debt product   n-gates parent missing fragile"
    )
    print("-" * 118)
    for row in rows:
        missing = ",".join(map(str, row.missing)) or "-"
        fragile = ",".join(map(str, row.fragile)) or "-"
        print(
            f"{row.n:<2} {row.label:<21} {row.classification:<14} "
            f"{fmt(row.gap_ratio):>8} {row.unprotected:>7} {fmt(row.product):>7} "
            f"{row.gates_n:>8} {row.parent_multiples:>6} {missing:<7} {fragile}"
        )
    print()
    print("Product-depth histograms:")
    for row in rows:
        if row.label in {"initial"}:
            continue
        print(
            f"  n={row.n:<2} {row.label:<21} "
            f"depths={depth_hist_text(row.n, row.depth_hist):<56} "
            f"frontier={fmt(row.frontier_mass):>6} pressure={row.denominator_pressure}"
        )
        owners = ", ".join(f"{owner}:{labels}/{unique}" for owner, labels, unique in row.top_owners)
        print(f"      top owners labels/unique = {owners or '-'}")
    print()


def print_repair_scans(scans: tuple[RepairScan, ...]) -> None:
    header("ROUND 3: ONE-SLOT REPAIR SCAN, MULTIPLE VS NONMULTIPLE")
    print("Scan: drop n-1 from the initial segment and add one speed a in [n, n(n-1)].")
    print("This is a small stress test of 'multiples are required' versus 'multiples create debt'.")
    print()
    print("n  drop  best multiple       gap/th debt   best nonmultiple   gap/th debt")
    print("-" * 86)
    for scan in scans:
        print(
            f"{scan.n:<2} {scan.drop:>4} "
            f"{scan.multiple_added:>14} {fmt(scan.multiple_gap):>8} {scan.multiple_debt:>4} "
            f"{scan.nonmultiple_added:>18} {fmt(scan.nonmultiple_gap):>8} {scan.nonmultiple_debt:>4}"
        )
    print()


def print_switch_table(rows: tuple[SwitchAudit, ...]) -> None:
    header("ROUND 4: TOURNAMENT ANALYSIS SWITCHBOARDS AT T=1/N")
    print("phase is circular halfturn orientation; safe switches edges by pairwise distance >= 1/n.")
    print()
    print("n  label                 t       phaseC safeC relC localC safe_scc safe_score_hist")
    print("-" * 118)
    for row in rows:
        print(
            f"{row.n:<2} {row.label:<21} {fmt(row.t):<7} "
            f"{row.phase_cycles:>6} {row.safe_cycles:>5} {row.relspeed_cycles:>4} "
            f"{row.local_gap_cycles:>6} {row.safe_scc:>8} {row.safe_score_hist}"
        )
    print()


def print_noise_synthesis() -> None:
    header("FORCED-RANDOM NOISE CARDS THAT SURVIVED CONTACT WITH COMPUTATION")
    cards = [
        (
            "covering-radius/zonotope",
            "Treat endpoint protection as a cell-cover invoice.  The n=18 local invoice is smaller than n=14: two bridge choices instead of six.",
        ),
        (
            "mixed-threshold Fourier",
            "Separate unit, half-gate, and bridge thresholds.  At n=18 the bridge is forced into the 6/12 triadic channel.",
        ),
        (
            "covering systems",
            "Small denominator rows are not enough.  The one-slot nonmultiple repairs can be tighter than multiples, but they carry more endpoint debt.",
        ),
        (
            "Tournament Analysis",
            "Initial rows can look ordered under safe-distance switches while ladders create cyclic pairwise crowding.",
        ),
    ]
    for name, finding in cards:
        print(f"  {name:<28} {finding}")
    print()
    print("Working proof guess:")
    print("  n=14 is the wide parity-bridge problem: forced odd fan plus a six-element even bridge fiber.")
    print("  n=18 is the square-core seam: forced unit fan plus half-gate 9 plus only bridges 6 or 12.")
    print("  I would currently attack n=18 first.  Its row-parent product is exactly 1 and its local")
    print("  bridge fiber has size 2, so a dual certificate has fewer bridge cases to break.")
    print("  Then push the resulting triadic debt argument back to n=14 as a 2 x 7 bridge-family")
    print("  inequality over the six even choices.")


def main() -> None:
    print("n=14/n=18 LRC alternating-noise exploration (codex-2026-06-01 S455)")
    print("All LRC endpoint computations are exact over Fraction; tournament switchboards use sampled t=1/n.\n")

    gate_audits = tuple(exact_gate_audit(n) for n in (14, 18))
    print_gate_audits(gate_audits)

    scans = tuple(one_slot_repair_scan(n) for n in (14, 18))
    rows: list[RowAudit] = []
    for n, scan in zip((14, 18), scans):
        rows.extend(
            [
                row_audit(n, "initial", initial(n)),
                row_audit(n, "drop-top add n", tuple(range(1, n - 1)) + (scan.multiple_added,)),
                row_audit(n, "drop-top nonmult", tuple(range(1, n - 1)) + (scan.nonmultiple_added,)),
                row_audit(n, "row-parent ladder", lpd_ladder(n)),
                row_audit(n, "n-gate ladder", n_gate_ladder(n)),
            ]
        )
    print_row_table(tuple(rows))
    print_repair_scans(scans)

    switch_rows: list[SwitchAudit] = []
    for n in (14, 18):
        switch_rows.extend(
            [
                switch_audit(n, "initial", initial(n), Fraction(1, n)),
                switch_audit(n, "row-parent ladder", lpd_ladder(n), Fraction(1, n)),
                switch_audit(n, "n-gate ladder", n_gate_ladder(n), Fraction(1, n)),
            ]
        )
    print_switch_table(tuple(switch_rows))
    print_noise_synthesis()


if __name__ == "__main__":
    main()
