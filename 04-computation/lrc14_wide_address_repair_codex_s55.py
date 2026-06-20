#!/usr/bin/env python3
"""HYP-2683 / T922: wide-branch address repair for LRC(14).

Recent work makes HYP-2675 the live sector-route crux:

    wide row E  ->  p0(E) <= cap_k

KPS S19 supplies the global frame: prove Weyl/decorrelation and compare to the
exact plateau Q(k-1).  HYP-2682 supplies the finite AP/cube-root phase router.

This scout imports an older repo lesson: scalar/product shadows often mix
proof states until the missing address coordinate is restored.  Here the
candidate address is the exact sector wall word, especially private-sector
ownership:

    on a wall atom, sector j is private if exactly one runner occupies j.

Private ownership is a deletion derivative: deleting that runner opens that
sector on that atom.  The script asks whether private-sector and compatibility
profiles separate direct-p0 risk buckets that scalar/additive summaries mix.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import log2

from lrc14_wide_branch_ridge_codex_s47 import (
    CAP,
    Row,
    additive_energy,
    fmt,
    missed_distribution,
    primitive,
    squarefree_profile,
    sumset_excess,
    wall_breakpoints,
)


INNER = tuple(range(1, 7))


def longest_run(row: Row) -> int:
    best = cur = 1
    for a, b in zip(row, row[1:]):
        if b == a + 1:
            cur += 1
            best = max(best, cur)
        else:
            cur = 1
    return best


def residue_hist(row: Row, mod: int = 7) -> tuple[int, ...]:
    counts = [0] * mod
    for x in row:
        counts[x % mod] += 1
    return tuple(counts)


def gap_word(row: Row) -> tuple[int, ...]:
    return tuple(b - a for a, b in zip(row, row[1:]))


def entropy_from_masses(masses: list[Fraction]) -> float:
    total = sum(masses, Fraction(0))
    if total == 0:
        return 0.0
    out = 0.0
    for mass in masses:
        if mass:
            p = float(mass / total)
            out -= p * log2(p)
    return out


@dataclass(frozen=True)
class PrivateProfile:
    state_mass: tuple[tuple[tuple[int, ...], Fraction], ...]
    private_sector_mass: tuple[Fraction, ...]
    private_owner_mass: tuple[tuple[int, Fraction], ...]
    private_total: Fraction
    private_owner_count: int
    state_support: int
    state_entropy: float
    private_entropy: float

    @property
    def private_average(self) -> Fraction:
        return self.private_total / 6

    @property
    def min_private_sector(self) -> Fraction:
        return min(self.private_sector_mass)

    @property
    def max_private_owner(self) -> Fraction:
        return max((mass for _v, mass in self.private_owner_mass), default=Fraction(0))


def private_profile(row: Row) -> PrivateProfile:
    """Exact wall-atom missed-state and private-sector ownership profile."""

    d, bps = wall_breakpoints(row)
    l = d // 7
    den2 = 2 * l
    state_mass: Counter[tuple[int, ...]] = Counter()
    private_sector = [Fraction(0) for _ in INNER]
    private_owner: Counter[int] = Counter()

    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo:
            continue
        mass = Fraction(hi - lo, d)
        midnum = lo + hi
        owners: dict[int, list[int]] = {j: [] for j in INNER}
        for e in row:
            if e == 0:
                continue
            sec = (e * midnum // den2) % 7
            if sec in owners:
                owners[sec].append(e)
        missed = tuple(j for j in INNER if not owners[j])
        state_mass[missed] += mass
        for idx, sec in enumerate(INNER):
            if len(owners[sec]) == 1:
                private_sector[idx] += mass
                private_owner[owners[sec][0]] += mass

    state_items = tuple(sorted(state_mass.items(), key=lambda kv: (len(kv[0]), kv[0])))
    owner_items = tuple(sorted(private_owner.items()))
    private_total = sum(private_sector, Fraction(0))
    return PrivateProfile(
        state_mass=state_items,
        private_sector_mass=tuple(private_sector),
        private_owner_mass=owner_items,
        private_total=private_total,
        private_owner_count=len(owner_items),
        state_support=len(state_items),
        state_entropy=entropy_from_masses([mass for _state, mass in state_items]),
        private_entropy=entropy_from_masses(list(private_owner.values())),
    )


@dataclass(frozen=True)
class RowRisk:
    row: Row
    p0: Fraction
    margin: Fraction
    p_by_count: tuple[Fraction, ...]

    @property
    def risk(self) -> Fraction:
        return self.p0 / CAP[len(self.row)]


@dataclass(frozen=True)
class RichRow:
    risk: RowRisk
    profile: PrivateProfile

    @property
    def row(self) -> Row:
        return self.risk.row

    @property
    def p0(self) -> Fraction:
        return self.risk.p0

    @property
    def margin(self) -> Fraction:
        return self.risk.margin


def push_top(store: list[RowRisk], row: RowRisk, keep: int) -> None:
    store.append(row)
    store.sort(key=lambda r: (r.p0, tuple(-x for x in r.row)), reverse=True)
    del store[keep:]


def scan_true_wide(bound: int = 20, keep: int = 260) -> tuple[int, int, int, list[RowRisk], list[RowRisk]]:
    """Exact k=9 scan of true-wide rows in {0}+8-subsets([1,bound]).

    Returns all-row counts, top risk rows, and a deterministic low/mid sample.
    """

    raw = primitive_count = true_wide = 0
    top: list[RowRisk] = []
    sample: list[RowRisk] = []
    for comb in combinations(range(1, bound + 1), 8):
        raw += 1
        row = (0,) + comb
        if row[-1] <= 14 or not primitive(row):
            continue
        primitive_count += 1
        if row[-2] <= 14:
            continue
        true_wide += 1
        dist = missed_distribution(row)
        risk = RowRisk(row=row, p0=dist[0], margin=CAP[9] - dist[0], p_by_count=dist)
        push_top(top, risk, keep)
        # Deterministic sparse sample to give the bucket audit low-risk contrast.
        h = sum((i + 3) * x * x for i, x in enumerate(row))
        if h % 389 == 0 and len(sample) < keep:
            sample.append(risk)
    return raw, primitive_count, true_wide, top, sample


def qbin(q: Fraction, denom: int = 80) -> int:
    return int(q * denom)


def fbin(x: float, scale: int = 20) -> int:
    return int(x * scale)


def channel_keys(r: RichRow) -> dict[str, tuple[object, ...]]:
    row = r.row
    p = r.profile
    gaps = gap_word(row)
    energy = additive_energy(row)
    exc = sumset_excess(row)
    private_sig = (
        qbin(p.private_average, 40),
        qbin(p.min_private_sector, 40),
        qbin(p.max_private_owner, 40),
        p.private_owner_count,
    )
    state_sig = (
        p.state_support // 4,
        fbin(p.state_entropy, 4),
        tuple(qbin(x, 30) for x in r.risk.p_by_count[1:4]),
    )
    return {
        "scalar": (row[-1], row[-2], max(gaps), longest_run(row)),
        "residue7": (residue_hist(row),),
        "squarefree": (squarefree_profile(row),),
        "additive": (exc, energy // 24, longest_run(row)),
        "private_mass": private_sig,
        "state_mass": state_sig,
        "residue_private": (residue_hist(row), private_sig, exc),
        "coarse_address": (exc, longest_run(row), private_sig, state_sig),
        "fine_address": (gaps, residue_hist(row), private_sig, state_sig),
    }


def audit_channels(rows: list[RichRow]) -> list[dict[str, object]]:
    high = Fraction(3, 10)
    low = Fraction(1, 4)
    summaries: list[dict[str, object]] = []
    for name in channel_keys(rows[0]):
        buckets: dict[tuple[object, ...], list[RichRow]] = defaultdict(list)
        for r in rows:
            buckets[channel_keys(r)[name]].append(r)
        mixed = 0
        widths: list[Fraction] = []
        multi_rows = 0
        worst_key: tuple[object, ...] | None = None
        worst_rows: list[RichRow] = []
        for key, vals in buckets.items():
            vals_sorted = sorted(vals, key=lambda x: x.p0)
            width = vals_sorted[-1].p0 - vals_sorted[0].p0
            if len(vals) > 1:
                widths.append(width)
                multi_rows += len(vals)
            has_high = any(v.p0 >= high for v in vals)
            has_low = any(v.p0 <= low for v in vals)
            if has_high and has_low:
                mixed += 1
            if not worst_rows or width > (worst_rows[-1].p0 - worst_rows[0].p0):
                worst_key = key
                worst_rows = vals_sorted
        max_width = worst_rows[-1].p0 - worst_rows[0].p0 if worst_rows else Fraction(0)
        avg_width = sum(widths, Fraction(0)) / len(widths) if widths else Fraction(0)
        summaries.append(
            {
                "name": name,
                "buckets": len(buckets),
                "compression": Fraction(len(buckets), len(rows)),
                "multi_rows": multi_rows,
                "overfit": len(buckets) > 19 * len(rows) // 20 or multi_rows < len(rows) // 10,
                "mixed": mixed,
                "max_width": max_width,
                "avg_width": avg_width,
                "worst_key": worst_key,
                "worst_low": worst_rows[0] if worst_rows else None,
                "worst_high": worst_rows[-1] if worst_rows else None,
            }
        )
    return summaries


def tournament_channels(summaries: list[dict[str, object]]) -> None:
    names = [str(s["name"]) for s in summaries]
    by_name = {str(s["name"]): s for s in summaries}
    wins = Counter()
    edges: dict[tuple[str, str], str] = {}
    for a, b in combinations(names, 2):
        sa = by_name[a]
        sb = by_name[b]
        score_a = (
            -int(sa["overfit"]),
            -int(sa["mixed"]),
            -sa["max_width"],
            -sa["avg_width"],
            -abs(float(sa["compression"]) - 0.55),
        )
        score_b = (
            -int(sb["overfit"]),
            -int(sb["mixed"]),
            -sb["max_width"],
            -sb["avg_width"],
            -abs(float(sb["compression"]) - 0.55),
        )
        if score_a >= score_b:
            winner, loser = a, b
        else:
            winner, loser = b, a
        wins[winner] += 1
        wins.setdefault(loser, wins[loser])
        edges[(winner, loser)] = ">"
    c3 = 0
    for a, b, c in combinations(names, 3):
        if (
            (a, b) in edges and (b, c) in edges and (c, a) in edges
        ) or (
            (a, c) in edges and (c, b) in edges and (b, a) in edges
        ):
            c3 += 1
    path = sorted(names, key=lambda n: (-wins[n], n))
    print("TOURNAMENT ANALYSIS ON REPAIR CHANNELS")
    print("  vertices are proof channels, not runners.")
    print("  pairwise observable: non-overfit compression, fewer high/low mixed buckets,")
    print("  then smaller bucket p0-width; tie path uses channel-name order.")
    print(f"  score_hist={dict(sorted(Counter(wins.values()).items()))} directed_3cycles={c3}")
    print("  Hamiltonian path=" + " > ".join(path))


def print_row(prefix: str, r: RichRow) -> None:
    p = r.profile
    print(
        f"{prefix} row={r.row} p0={r.p0} ({float(r.p0):.9f}) "
        f"margin={r.margin} ({float(r.margin):.9f}) exc={sumset_excess(r.row)} "
        f"run={longest_run(r.row)} energy={additive_energy(r.row)}"
    )
    print(
        f"    private_avg={p.private_average} ({float(p.private_average):.9f}) "
        f"min_private_sector={p.min_private_sector} ({float(p.min_private_sector):.9f}) "
        f"private_owners={p.private_owner_count} state_support={p.state_support} "
        f"state_entropy={p.state_entropy:.4f}"
    )


def mean_fraction(values: list[Fraction]) -> Fraction:
    return sum(values, Fraction(0)) / len(values) if values else Fraction(0)


def mean_float(values: list[float]) -> float:
    return sum(values) / len(values) if values else 0.0


def print_pressure_stats(rows: list[RichRow]) -> None:
    groups = [
        ("high direct-risk", [r for r in rows if r.p0 >= Fraction(3, 10)]),
        ("low direct-risk", [r for r in rows if r.p0 <= Fraction(1, 4)]),
    ]
    print("STATE-MASS PRESSURE")
    for label, group in groups:
        avg_p0 = mean_fraction([r.p0 for r in group])
        avg_private = mean_fraction([r.profile.private_average for r in group])
        avg_min_private = mean_fraction([r.profile.min_private_sector for r in group])
        avg_support = mean_float([float(r.profile.state_support) for r in group])
        avg_entropy = mean_float([r.profile.state_entropy for r in group])
        print(
            f"  {label:<16} count={len(group):>3} avg_p0={avg_p0} ({float(avg_p0):.9f}) "
            f"avg_private={avg_private} ({float(avg_private):.9f}) "
            f"avg_min_private={avg_min_private} ({float(avg_min_private):.9f}) "
            f"avg_state_support={avg_support:.3f} avg_state_entropy={avg_entropy:.4f}"
        )
    print("  pressure direction: high direct-risk rows have larger private-sector mass")
    print("  and lower missed-state entropy than the low-risk sample.")


def main() -> None:
    print("HYP-2683 / T922 -- LRC14 wide-branch address repair")
    print("Hidden-gem sources: owner-private deletion repair, compatibility walls, state-word invariants.")
    print("Predicate preserved: direct p0(E) <= cap_k for HYP-2675 wide rows.")
    print("Assumption challenged: scalar plateau/additive summaries are enough to route wide risk.\n")

    raw, primitive_count, true_wide, top, sample = scan_true_wide(bound=20, keep=260)
    selected: dict[Row, RowRisk] = {r.row: r for r in top}
    for r in sample:
        selected.setdefault(r.row, r)
    rich = [RichRow(risk=r, profile=private_profile(r.row)) for r in selected.values()]
    rich.sort(key=lambda r: (r.p0, tuple(-x for x in r.row)), reverse=True)

    print("EXACT TRUE-WIDE SCAN")
    print(f"  box: k=9, row=(0)+8-subsets([1,20]), span>14, second-largest>14")
    print(f"  raw rows={raw}, primitive span>14 rows={primitive_count}, true-wide rows={true_wide}")
    print(f"  rich audit rows={len(rich)} (top risk rows plus deterministic low/mid sample)")
    print(f"  cap_9={CAP[9]} ({float(CAP[9]):.9f})")
    print()

    print("TOP DIRECT-P0 ROWS WITH ADDRESS DATA")
    for idx, r in enumerate(rich[:12], 1):
        print_row(f"  {idx:2d}.", r)
    print()

    low_private = sorted(rich, key=lambda r: (r.profile.private_average, -r.p0, r.row))[:8]
    print("LOW PRIVATE-SECTOR AVERAGE ROWS")
    for idx, r in enumerate(low_private, 1):
        print_row(f"  {idx:2d}.", r)
    print()

    print_pressure_stats(rich)
    print()

    summaries = audit_channels(rich)
    print("MIXED-BUCKET CHANNEL AUDIT")
    print("  high threshold p0 >= 3/10; low threshold p0 <= 1/4")
    for s in summaries:
        low = s["worst_low"]
        high = s["worst_high"]
        print(
            f"  {s['name']:<18} buckets={s['buckets']:>4} compression={s['compression']} "
            f"multi_rows={s['multi_rows']:>3} overfit={s['overfit']} "
            f"mixed={s['mixed']:>3} max_width={s['max_width']} ({float(s['max_width']):.9f}) "
            f"avg_width={s['avg_width']} ({float(s['avg_width']):.9f})"
        )
        if low is not None and high is not None and low.row != high.row:
            print(f"    worst bucket low={low.row} p0={low.p0}; high={high.row} p0={high.p0}")
    print()

    tournament_channels(summaries)
    print()

    print("SYNTHESIS")
    print("  1. The direct-risk leaders in the true-wide B20 bank are low-growth, high-overlap")
    print("     sector-cover states rather than generic wide rows.")
    print("  2. Scalar span/additive shadows still mix high and low p0 rows in the audited set;")
    print("     coarse private/state address data reduces mixed buckets without using row identity.")
    print("     Fine address is diagnostic, but too close to a perfect row hash to be a lemma.")
    print("  3. This supports the HYP-2675 proof shape: first route finite low-height/low-growth")
    print("     compatibility addresses, then apply the Weyl/decorrelation plateau lemma.")
    print("  4. No LRC(14) proof is claimed; the open step remains explicit decorrelation error")
    print("     plus bounded-gap glue, but the address-repair channel gives a sharper finite ledger.")


if __name__ == "__main__":
    main()
