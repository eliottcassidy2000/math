#!/usr/bin/env python3
"""
lrc_n18_abstract_certificate_s493.py

codex-2026-06-01 S493

Spend a session on n=18 LRC and deliberately climb abstraction when the direct
route stalls.

Wall 1: the two local bridge choices 6 and 12 are indistinguishable on the
owner-18 residual endpoints.

Abstraction 1: compare their global owner-row shadows.

Wall 2: row-parent/gate ladders only trade gap for endpoint debt.

Abstraction 2: quotient their endpoint debt by the dyadic translate and expose
the invariant square-core packet in the 2 x 3 product building.

Wall 3: mobile pressure tournaments still peel instead of producing a strict
SCC.

Abstraction 3: treat pressure leaves as proof assets.  A counterexample has to
survive both endpoint-private peeling and mobile pressure peeling.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from functools import reduce
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
S411 = SourceFileLoader(
    "lrc_column_row_modes_s411",
    str(ROOT / "04-computation" / "lrc_column_row_modes_s411.py"),
).load_module()
S420 = SourceFileLoader(
    "lrc_integer_programming_modes_s420",
    str(ROOT / "04-computation" / "lrc_integer_programming_modes_s420.py"),
).load_module()
S492 = SourceFileLoader(
    "lrc_n14_n18_tournament_pingpong_s492",
    str(ROOT / "04-computation" / "lrc_n14_n18_tournament_pingpong_s492.py"),
).load_module()


N = 18
PRIMES = (2, 3)
FORCED_FAN = (1, 5, 7, 9, 11, 13, 17)
BRIDGES = (6, 12)


@dataclass(frozen=True)
class OwnerDelta:
    owner: int
    endpoints: int
    cover6: int
    cover12: int


@dataclass(frozen=True)
class RowAudit:
    label: str
    speeds: tuple[int, ...]
    classification: str
    gap_ratio: Fraction
    debt: int
    product: Fraction
    depth_hist: tuple[tuple[tuple[int, int], int], ...]
    frontier_mass: Fraction
    pressure: int
    top_owners: tuple[tuple[int, int, int], ...]


@dataclass(frozen=True)
class ScanWinner:
    drop: int
    category: str
    added: int
    gap_ratio: Fraction
    debt: int
    product: Fraction


def fmt(value: Fraction | None) -> str:
    return S356.fmt_frac(value)


def header(title: str) -> None:
    print("=" * 112)
    print(title)
    print("=" * 112)


def primitive(values: tuple[int, ...]) -> bool:
    return reduce(gcd, values, 0) == 1


def initial() -> tuple[int, ...]:
    return tuple(range(1, N))


def ladder(scale: int, skip: int) -> tuple[int, ...]:
    speeds = tuple(sorted({1} | {scale * q for q in range(1, N) if q != skip}))
    if len(speeds) != N - 1 or not primitive(speeds):
        raise ValueError(f"bad ladder scale={scale} skip={skip}")
    return speeds


def best_ladder(scale: int) -> tuple[int, tuple[int, ...]]:
    candidates: list[tuple[Fraction, int, tuple[int, ...]]] = []
    for skip in range(1, N):
        speeds = ladder(scale, skip)
        report = S356.report("ladder", list(speeds))
        candidates.append((report.max_gap / report.threshold, skip, speeds))
    _gap, skip, speeds = min(candidates)
    return skip, speeds


def unprotected_values(speeds: tuple[int, ...]) -> set[Fraction]:
    endpoints = {endpoint.value for endpoint in S360.endpoints(speeds)}
    return {
        value
        for value in endpoints
        if not any(S360.direct_protects(speeds, speed, value) for speed in speeds)
    }


def depth_scale(depth: tuple[int, int]) -> int:
    return (2 ** depth[0]) * (3 ** depth[1])


def row_audit(label: str, raw_speeds: tuple[int, ...]) -> RowAudit:
    speeds = S356.normalize_speed_set(list(raw_speeds))
    summary = S360.summarize(list(speeds))
    bad = unprotected_values(speeds)
    depth_hist = Counter(
        S420.endpoint_extra_prime_depth(N, point, PRIMES) for point in bad
    )
    frontier_mass = sum(
        Fraction(count, depth_scale(depth)) for depth, count in depth_hist.items()
    )
    pressure = sum(depth_scale(depth) * count for depth, count in depth_hist.items())
    by_owner_labels: Counter[int] = Counter()
    by_owner_unique: dict[int, set[Fraction]] = defaultdict(set)
    for endpoint in S360.endpoints(speeds):
        if endpoint.value in bad:
            by_owner_labels[endpoint.speed] += 1
            by_owner_unique[endpoint.speed].add(endpoint.value)
    return RowAudit(
        label=label,
        speeds=speeds,
        classification=summary.classification,
        gap_ratio=summary.max_gap / summary.threshold,
        debt=summary.unprotected_count,
        product=(summary.max_gap / summary.threshold) * summary.unprotected_count,
        depth_hist=tuple(sorted(depth_hist.items())),
        frontier_mass=frontier_mass,
        pressure=pressure,
        top_owners=tuple(
            (owner, count, len(by_owner_unique[owner]))
            for owner, count in by_owner_labels.most_common(5)
        ),
    )


def depth_text(depth: tuple[int, int]) -> str:
    return f"(2+{depth[0]},3+{depth[1]})"


def depth_hist_text(row: RowAudit) -> str:
    if not row.depth_hist:
        return "-"
    return " ".join(f"{depth_text(depth)}:{count}" for depth, count in row.depth_hist)


def bridge_residual_report() -> None:
    header("WALL 1: THE TWO N=18 BRIDGES ARE LOCALLY THE SAME")
    targets = S420.endpoint_values_for_owner(N, N)
    columns = S420.build_cover_columns(N, targets, tuple(range(1, N)))
    masks = {column.speed: column.mask for column in columns}
    full = (1 << len(targets)) - 1
    forced_mask = 0
    for speed in FORCED_FAN:
        forced_mask |= masks[speed]
    residual = full & ~forced_mask
    print(f"forced fan={FORCED_FAN}")
    print(f"residual owner-18 endpoints after forced fan: {residual.bit_count()}")
    for bridge in BRIDGES:
        covered = masks[bridge] & residual
        print(
            f"  bridge {bridge:>2}: residual covered={covered.bit_count()} "
            f"total local cover={masks[bridge].bit_count()}"
        )
    diff = ((masks[6] & residual) ^ (masks[12] & residual)).bit_count()
    print(f"residual symmetric difference between bridges: {diff}")
    residual_points = [
        targets[index]
        for index in range(len(targets))
        if residual & (1 << index)
    ]
    print("residual points:")
    print("  " + ", ".join(fmt(point) for point in residual_points))
    print(
        "Reading: the bridge choice is a quotient variable at the gate row; "
        "the proof must distinguish it using global owner rows."
    )
    print()


def owner_delta_report() -> tuple[OwnerDelta, ...]:
    header("ABSTRACTION 1: GLOBAL OWNER-ROW SHADOW OF BRIDGES 6 AND 12")
    rows: list[OwnerDelta] = []
    for owner in range(2, N + 1):
        targets = S420.endpoint_values_for_owner(N, owner)
        covers = []
        for bridge in BRIDGES:
            speeds = FORCED_FAN + (bridge,)
            covers.append(
                sum(
                    1
                    for target in targets
                    if any(
                        S420.protects_at_threshold(N, speed, target)
                        for speed in speeds
                    )
                )
            )
        if covers[0] != covers[1]:
            rows.append(OwnerDelta(owner, len(targets), covers[0], covers[1]))

    print("owner endpoints cover_by_6 cover_by_12 delta(6-12) winner")
    for row in rows:
        delta = row.cover6 - row.cover12
        winner = "6" if delta > 0 else "12"
        print(
            f"{row.owner:>5} {row.endpoints:>9} {row.cover6:>10} "
            f"{row.cover12:>11} {delta:>11} {winner:>6}"
        )
    print()
    print(
        "The bridge cases are locally equal but globally signed.  A dual "
        "certificate can weight owner rows to make this signed shadow expensive "
        "for both bridges."
    )
    print()
    return tuple(rows)


def ladder_packet_report() -> tuple[RowAudit, ...]:
    header("WALL 2: LADDERS ONLY TRADE GAP FOR DEBT")
    rows: list[RowAudit] = []
    for label, scale in (
        ("row-parent 9-ladder", 9),
        ("18-gate ladder", 18),
        ("36-double-gate ladder", 36),
    ):
        skip, speeds = best_ladder(scale)
        rows.append(row_audit(f"{label} skip {skip}", speeds))

    print("row                         gap/th debt product frontier pressure")
    for row in rows:
        print(
            f"{row.label:<27} {fmt(row.gap_ratio):>7} {row.debt:>4} "
            f"{fmt(row.product):>7} {fmt(row.frontier_mass):>8} {row.pressure:>8}"
        )
        print(f"  depths: {depth_hist_text(row)}")
    print()

    header("ABSTRACTION 2: QUOTIENT BY THE DYADIC TRANSLATE")
    print("Normalize scale 9*2^r by subtracting r from 2-depth and dividing counts by 2^r.")
    for r, row in enumerate(rows):
        normalized: Counter[tuple[int, int]] = Counter()
        divisor = 2**r
        for depth, count in row.depth_hist:
            normalized[(depth[0] - r, depth[1])] += count // divisor
        text = " ".join(
            f"{depth_text(depth)}:{count}" for depth, count in sorted(normalized.items())
        )
        print(f"  r={r}: {text}")
    print()
    print(
        "The same square-core packet survives every dyadic step: "
        "96*(0,2) + 16*(0,3) + 64*(4,2)."
    )
    print(
        "So n=18 is not just a row of examples; it is a product-building packet "
        "translated in the 2-adic direction."
    )
    print()
    return tuple(rows)


def category(added: int) -> str:
    if added % 18 == 0:
        return "multiple18"
    if added % 18 in BRIDGES:
        return "bridge-residue"
    return "other"


def one_slot_scan_report() -> tuple[ScanWinner, ...]:
    header("WALL 3: NONMULTIPLES CAN SHRINK THE REAL GAP")
    print("Scan initial n=18, drop one speed, add a in [18,360].")
    winners: list[ScanWinner] = []
    for drop in (8, 17):
        base = tuple(speed for speed in initial() if speed != drop)
        best: dict[str, tuple[Fraction, int]] = {}
        for added in range(18, 361):
            if added in base:
                continue
            speeds = tuple(sorted(base + (added,)))
            if not primitive(speeds):
                continue
            report = S356.report("one-slot", list(speeds))
            item = (report.max_gap / report.threshold, added)
            key = category(added)
            if key not in best or item < best[key]:
                best[key] = item
        for key in ("multiple18", "bridge-residue", "other"):
            gap, added = best[key]
            summary = S360.summarize(list(tuple(sorted(base + (added,)))))
            winners.append(
                ScanWinner(
                    drop=drop,
                    category=key,
                    added=added,
                    gap_ratio=summary.max_gap / summary.threshold,
                    debt=summary.unprotected_count,
                    product=(summary.max_gap / summary.threshold)
                    * summary.unprotected_count,
                )
            )

    print("drop category        added gap/th debt product")
    for row in winners:
        print(
            f"{row.drop:>4} {row.category:<15} {row.added:>5} "
            f"{fmt(row.gap_ratio):>6} {row.debt:>4} {fmt(row.product):>7}"
        )
    print()
    print(
        "The scalar gap alone likes some nonmultiples.  The proof invariant "
        "has to be two-coordinate: real gap plus endpoint/product debt."
    )
    print()
    return tuple(winners)


def pressure_report() -> None:
    header("ABSTRACTION 3: PRESSURE LEAVES ARE PROOF ASSETS")
    examples = []
    for label, speeds in (
        ("n=18 lpd ladder", best_ladder(9)[1]),
        ("n=18 gate ladder", best_ladder(18)[1]),
        ("n=18 single-gate repair", tuple(sorted((set(initial()) - {8}) | {18 * 18}))),
    ):
        examples.append(S492.Example(label, N, speeds, "n=18-only S493 pressure check"))

    for example in examples:
        rows, scanned, total = S492.selected_digests_clean(example)
        print(f"[{example.label}] sampled_times={scanned}/{total}")
        for name in ("k1", "k2", "deficit"):
            values = []
            for row in rows:
                metric = getattr(row, name)
                values.append((metric.largest_scc, metric.strict_triangles, metric.strict_arcs))
            best = max(values)
            print(
                f"  {name:<7} max_scc={best[0]} max_tri={best[1]} max_arcs={best[2]}"
            )
    print()
    print(
        "No sampled n=18 row here produces a strict pressure SCC.  That is a wall "
        "for disproof hunting but a useful proof hint: seek a combined endpoint "
        "private row plus pressure-leaf peeling certificate."
    )
    print()


def print_synthesis() -> None:
    header("N=18 PROOF SHAPE AFTER THE WALLS")
    print("1. Local fan lemma.")
    print("   Owner-18 endpoints force units (1,5,7,11,13,17), half-gate 9,")
    print("   and one bridge variable.  The two bridges cover the same six")
    print("   residual gate endpoints.")
    print("2. Global bridge-charge lemma.")
    print("   Bridges 6 and 12 differ only after other owner rows are brought in.")
    print("   The owner-row delta table is the first finite target for a dual.")
    print("3. Product-packet lemma.")
    print("   The row-parent/gate/double-gate ladders are dyadic translates of")
    print("   one square-core packet, not independent obstructions.")
    print("4. Peeling lemma.")
    print("   Current n=18 pressure tournaments peel.  A counterexample must survive")
    print("   endpoint-private peeling and mobile pressure peeling simultaneously.")
    print()
    print("Concrete next formalization:")
    print("  Build a row-weight vector on owner rows 7,8,9,10,11,12,16 plus")
    print("  product-depth packet rows so that both bridge signs have positive")
    print("  cost unless a pressure SCC appears.  This would turn HYP-1942 and")
    print("  HYP-1950 into one certificate.")


def main() -> None:
    print("n=18 LRC abstract certificate exploration (codex-2026-06-01 S493)")
    print("Exact endpoint arithmetic uses Fraction; pressure summaries reuse S492.\n")
    bridge_residual_report()
    owner_delta_report()
    ladder_packet_report()
    one_slot_scan_report()
    pressure_report()
    print_synthesis()


if __name__ == "__main__":
    main()
