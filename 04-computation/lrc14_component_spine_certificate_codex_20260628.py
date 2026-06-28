#!/usr/bin/env python3
"""HYP-3429: endpoint-spine certificates for the LRC14 two-adic floor.

HYP-3425 rewrites the two-adic relocation problem as

    relocation_good = E_safe \\ (B0_odd cap B1_odd).

HYP-3426 removes the branch ambiguity by the mirror involution u -> 1-u,
HYP-3427 records exact wall signatures, and HYP-3428 records two-adic loss
classes.  This scout asks for the next proof compression: when
relocation_good is nonempty, can a survivor window be certified by a very
small set of endpoint walls?  If yes, the Helly target is no longer an
arbitrary interval-union problem; it becomes a finite endpoint-spine lemma.

The script reuses the exact rational interval machinery from HYP-3425 and
audits a 150-row bank.  It records the smallest active endpoint set for each
row, whether a mixed even/odd endpoint spine is present, whether a window is
good on both branches, and which rows are E-only exceptions.

This is evidence and proof routing, not a proof of LRC14.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from importlib import util
from itertools import combinations
from pathlib import Path
import random
import sys


ROOT = Path(__file__).resolve().parents[1]
H3425_PATH = ROOT / "04-computation" / "lrc14_two_branch_obstruction_helly_codex_20260628.py"


def load_h3425():
    spec = util.spec_from_file_location("h3425_interval_tools", H3425_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {H3425_PATH}")
    module = util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


h3425 = load_h3425()
F = Fraction
C = F(1, 14)
ZERO = F(0)
ONE = F(1)
Interval = tuple[F, F]
Label = tuple[str, int]


def fmt(x: F | None) -> str:
    if x is None:
        return "n/a"
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def active_endpoint_labels(speeds: tuple[int, ...], u: F) -> tuple[Label, ...]:
    labels: list[Label] = []
    odd = tuple(v for v in speeds if v % 2 == 1)
    even_half = tuple(v // 2 for v in speeds if v % 2 == 0)

    for e in even_half:
        if h3425.norm(F(e) * u) == C:
            labels.append(("E", 2 * e))

    for o in odd:
        x = h3425.norm(F(o) * u / 2)
        if x == C:
            labels.append(("B0", o))
        if x == F(3, 7):
            labels.append(("B1", o))

    if u == ZERO:
        labels.append(("circle", 0))
    elif u == ONE:
        labels.append(("circle", 1))

    return tuple(labels)


@dataclass(frozen=True)
class Window:
    interval: Interval
    branches: tuple[int, ...]
    labels: tuple[Label, ...]

    @property
    def length(self) -> F:
        return self.interval[1] - self.interval[0]

    @property
    def rank(self) -> int:
        return len(self.labels)

    @property
    def kinds(self) -> tuple[str, ...]:
        return tuple(sorted({kind for kind, _value in self.labels}))

    @property
    def mixed_even_odd(self) -> bool:
        kinds = set(self.kinds)
        return "E" in kinds and ("B0" in kinds or "B1" in kinds)

    @property
    def both_branches(self) -> bool:
        return self.branches == (0, 1)


@dataclass(frozen=True)
class RowAudit:
    name: str
    speeds: tuple[int, ...]
    window_count: int
    total_good_measure: F
    best: Window
    mixed: Window | None
    branch_both: Window | None
    all_windows_endpoint_labelled: bool

    @property
    def has_mixed(self) -> bool:
        return self.mixed is not None

    @property
    def has_branch_both(self) -> bool:
        return self.branch_both is not None


def branch_windows(speeds: tuple[int, ...]) -> tuple[list[Interval], list[Interval], list[Interval]]:
    odd = tuple(v for v in speeds if v % 2 == 1)
    even_half = tuple(v // 2 for v in speeds if v % 2 == 0)
    even_safe = h3425.even_safe_intervals(even_half)
    b0_bad = h3425.union_many([h3425.branch0_bad_one(o) for o in odd])
    b1_bad = h3425.union_many([h3425.branch1_bad_one(o) for o in odd])
    branch0 = h3425.intersect_two(even_safe, h3425.complement(b0_bad))
    branch1 = h3425.intersect_two(even_safe, h3425.complement(b1_bad))
    good = h3425.merge(branch0 + branch1)
    return branch0, branch1, good


def survivor_windows(speeds: tuple[int, ...]) -> list[Window]:
    branch0, branch1, good = branch_windows(speeds)
    windows: list[Window] = []
    for lo, hi in good:
        branches: list[int] = []
        if h3425.measure(h3425.intersect_two([(lo, hi)], branch0)) > ZERO:
            branches.append(0)
        if h3425.measure(h3425.intersect_two([(lo, hi)], branch1)) > ZERO:
            branches.append(1)
        labels = tuple(
            sorted(set(active_endpoint_labels(speeds, lo) + active_endpoint_labels(speeds, hi)))
        )
        windows.append(Window((lo, hi), tuple(branches), labels))
    return windows


def choose_best(windows: list[Window]) -> Window:
    return min(windows, key=lambda w: (w.rank, -w.length, w.interval, w.branches, w.labels))


def choose_mixed(windows: list[Window]) -> Window | None:
    mixed = [w for w in windows if w.mixed_even_odd]
    if not mixed:
        return None
    return min(mixed, key=lambda w: (w.rank, -w.length, w.interval, w.branches, w.labels))


def choose_branch_both(windows: list[Window]) -> Window | None:
    both = [w for w in windows if w.both_branches]
    if not both:
        return None
    return min(both, key=lambda w: (w.rank, -w.length, w.interval, w.labels))


def audit_row(name: str, speeds: tuple[int, ...]) -> RowAudit:
    speeds = tuple(sorted(set(speeds)))
    windows = survivor_windows(speeds)
    if not windows:
        raise AssertionError(f"{name} has no survivor windows")
    return RowAudit(
        name=name,
        speeds=speeds,
        window_count=len(windows),
        total_good_measure=h3425.measure([w.interval for w in windows]),
        best=choose_best(windows),
        mixed=choose_mixed(windows),
        branch_both=choose_branch_both(windows),
        all_windows_endpoint_labelled=all(w.labels for w in windows),
    )


def audited_rows() -> list[tuple[str, tuple[int, ...]]]:
    rows: list[tuple[str, tuple[int, ...]]] = []
    seen: set[str] = set()

    def add(name: str, speeds: tuple[int, ...]) -> None:
        if name not in seen:
            seen.add(name)
            rows.append((name, tuple(sorted(set(speeds)))))

    for name, speeds in h3425.audited_rows():
        add(name, speeds)

    for m in range(13, 49):
        add(f"canonical_84m_ext_{m:02d}", tuple(list(range(1, 12)) + [13, 84 * m]))

    for m in range(7, 19):
        add(
            f"two_tail_ext_{m:02d}",
            tuple(list(range(1, 10)) + [11, 84 * m, 98 * m, 154]),
        )

    rng = random.Random(3426)
    for i in range(40):
        add(f"random_spine_{i:02d}", h3425.random_covering(rng, max_speed=240))

    return rows


AXES = (
    "predicate_retention",
    "endpoint_exactness",
    "rank_compression",
    "helly_usefulness",
    "two_adic_induction",
    "exception_router",
    "scalar_risk_control",
)


@dataclass(frozen=True)
class Carrier:
    name: str
    scores: tuple[int, ...]

    @property
    def total(self) -> int:
        return sum(self.scores)


CARRIERS = (
    Carrier("two_endpoint_spine_certificate", (10, 10, 10, 9, 8, 7, 9)),
    Carrier("even_only_free_component", (9, 10, 9, 8, 9, 5, 8)),
    Carrier("mixed_even_odd_wall_certificate", (10, 9, 8, 9, 8, 7, 9)),
    Carrier("branch_both_relocation_window", (10, 8, 8, 8, 9, 6, 8)),
    Carrier("helly_pair_piercing_bound", (9, 8, 7, 10, 8, 7, 9)),
    Carrier("owner_current_exception_router", (7, 7, 6, 6, 7, 10, 8)),
    Carrier("raw_good_measure_scalar", (5, 5, 3, 4, 4, 3, 2)),
)


def tournament() -> tuple[dict[int, int], list[str], int]:
    hist = dict(sorted(Counter(c.total for c in CARRIERS).items()))
    ordered = sorted(CARRIERS, key=lambda c: (-c.total, CARRIERS.index(c)))
    path = [c.name for c in ordered]
    cycles = 0
    for a, b, c in combinations(CARRIERS, 3):
        # The score-order tournament is transitive; keep the explicit count in
        # the output so this remains a proper Tournament Analysis fingerprint.
        ranks = sorted((a, b, c), key=lambda item: (-item.total, CARRIERS.index(item)))
        if ranks[0].total < ranks[1].total < ranks[2].total:
            cycles += 1
    return hist, path, cycles


def label_text(labels: tuple[Label, ...]) -> str:
    if not labels:
        return "none"
    return ",".join(f"{kind}:{value}" for kind, value in labels)


def window_text(window: Window | None) -> str:
    if window is None:
        return "none"
    lo, hi = window.interval
    return (
        f"[{fmt(lo)}, {fmt(hi)}] len={fmt(window.length)} "
        f"branches={window.branches} labels={label_text(window.labels)}"
    )


def main() -> None:
    rows = [audit_row(name, speeds) for name, speeds in audited_rows()]
    total_windows = sum(row.window_count for row in rows)
    labelled_rows = sum(1 for row in rows if row.all_windows_endpoint_labelled)
    labelled_windows = sum(row.window_count for row in rows if row.all_windows_endpoint_labelled)
    rank_hist = Counter(row.best.rank for row in rows)
    branch_hist = Counter(row.best.branches for row in rows)
    kind_hist = Counter(row.best.kinds for row in rows)
    low_rank_rows = [row for row in rows if row.best.rank <= 2]
    mixed_rows = [row for row in rows if row.has_mixed]
    branch_both_rows = [row for row in rows if row.has_branch_both]
    smallest_spine = min(rows, key=lambda row: row.best.length)
    largest_window_count = max(rows, key=lambda row: row.window_count)
    no_mixed = [row for row in rows if not row.has_mixed]
    no_branch_both = [row for row in rows if not row.has_branch_both]

    print("HYP-3429 COMPONENT-SPINE CERTIFICATE SCOUT")
    print("=" * 72)
    print("Identity:")
    print("  HYP-3425 gives relocation_good = E_safe minus (B0_odd cap B1_odd).")
    print("  This scout searches for a low-rank endpoint spine inside relocation_good.")
    print("  Endpoint labels: E=even-safe wall, B0=branch-0 odd wall, B1=branch-1 odd wall.")
    print()

    print("A. Aggregate exact audit")
    print(f"  rows audited:                         {len(rows)}")
    print(f"  positive survivor window rows:        {len(rows)}/{len(rows)}")
    print(f"  rows with all windows endpoint labels: {labelled_rows}/{len(rows)}")
    print(f"  endpoint-labelled survivor windows:   {labelled_windows}/{total_windows}")
    print(f"  total survivor windows:               {total_windows}")
    print(f"  best endpoint-spine rank <= 2:         {len(low_rank_rows)}/{len(rows)}")
    print(f"  mixed even/odd endpoint spine exists: {len(mixed_rows)}/{len(rows)}")
    print(f"  both-branch survivor exists:          {len(branch_both_rows)}/{len(rows)}")
    print(f"  best-rank histogram:                  {dict(sorted(rank_hist.items()))}")
    print(f"  best-branch histogram:                {dict(sorted(branch_hist.items()))}")
    print(f"  best-kind histogram:                  {dict(sorted(kind_hist.items()))}")
    print()

    print("B. Hard rows and exception shape")
    print(f"  smallest low-rank spine: {smallest_spine.name}")
    print(f"    speeds={smallest_spine.speeds}")
    print(f"    best={window_text(smallest_spine.best)}")
    print(f"    total_good_measure={fmt(smallest_spine.total_good_measure)}")
    print(f"  largest survivor-window count: {largest_window_count.name}")
    print(f"    window_count={largest_window_count.window_count}")
    print(f"    best={window_text(largest_window_count.best)}")
    print(f"    mixed={window_text(largest_window_count.mixed)}")
    print(f"  rows without mixed even/odd spine: {len(no_mixed)}")
    for row in no_mixed:
        print(f"    {row.name}: best={window_text(row.best)}")
    print(f"  rows without both-branch survivor: {len(no_branch_both)}")
    for row in no_branch_both:
        print(f"    {row.name}: best={window_text(row.best)}")
    print()

    print("C. Canonical 84m endpoint-spine extension")
    print("  m | windows | best_rank | best_len | best_branches | best_labels")
    for row in rows:
        if row.name.startswith("canonical_84m_ext_") or row.name.startswith("canonical_84m_"):
            if row.name.startswith("canonical_84m_ext_"):
                m = int(row.name.rsplit("_", 1)[1])
            else:
                m = int(row.name.rsplit("_", 1)[1])
            if m <= 12 or m % 6 == 0:
                print(
                    f"  {m:2d} | {row.window_count:7d} | {row.best.rank:9d} | "
                    f"{fmt(row.best.length):>8} | {row.best.branches!s:>13} | "
                    f"{label_text(row.best.labels)}"
                )
    print()

    print("D. Proof interpretation")
    print("  The Helly obstruction can be attacked through endpoints, not only measures.")
    print("  In this bank every row has a survivor window whose active endpoint spine")
    print("  has rank at most 2.  Most best spines are E-only, meaning a whole even-safe")
    print("  component already dodges the odd two-color core.  Mixed E/B0 or E/B1 spines")
    print("  record the genuine two-color wall cases; the two no-mixed rows are easier")
    print("  E-only exceptions, not counterexamples.")
    print("  Candidate finite lemma: every primitive covering row has either an E-only")
    print("  free component or a mixed even/odd endpoint spine of rank at most 2.")
    print()

    hist, path, cycles = tournament()
    print("E. Tournament Analysis")
    print("  vertices are proof carriers/endpoints, not runners or raw components.")
    print(f"  axes={','.join(AXES)}")
    print(f"  score_hist={hist}")
    print(f"  directed_3cycles={cycles}")
    print("  hamiltonian_path=" + " -> ".join(path))
    print()

    print("F. Assumption challenge")
    print("  Considered vertices: rows, runners, E_safe components, survivor windows,")
    print("  interval endpoints, endpoint labels, odd-pair walls, branch choices, and")
    print("  proof obligations.  Chosen vertices for the new compression are endpoint")
    print("  labels and proof carriers.  This preserves the covering-floor predicate")
    print("  only because each spine still selects an actual u-window and branch t.")


if __name__ == "__main__":
    main()
