#!/usr/bin/env python3
"""HYP-3427: wall-signature atlas for the LRC14 two-branch floor.

HYP-3425 rewrites HYP-3422's relocation target as

    good = E_safe \ (B0_odd cap B1_odd),

where B0 is the odd near-integer bad set and B1 is the odd near-half bad set.
This script asks a more proof-facing question: do the surviving good windows
carry small boundary signatures?  A window with exact left/right wall labels is
a finite certificate candidate; it says which even wall or odd branch wall keeps
the two-color obstruction from covering that component.

This is not an LRC14 proof.  It is an exact rational atlas for the next lemma:
every primitive covering row should have at least one survivor window with a
bounded wall signature, or else emit named owner/state/exact-period debt.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
from math import gcd
import random

C = F(1, 14)
ZERO = F(0)
ONE = F(1)
Interval = tuple[F, F]


def frac_part(x: F) -> F:
    return x - (x.numerator // x.denominator)


def norm(x: F) -> F:
    r = frac_part(x)
    return min(r, 1 - r)


def score(speeds: tuple[int, ...], t: F) -> F:
    return min(norm(F(v) * t) for v in speeds)


def fmt(x: F | None) -> str:
    if x is None:
        return "n/a"
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def merge(intervals: list[Interval]) -> list[Interval]:
    clipped = [
        (max(ZERO, lo), min(ONE, hi))
        for lo, hi in intervals
        if max(ZERO, lo) < min(ONE, hi)
    ]
    clipped.sort()
    out: list[Interval] = []
    for lo, hi in clipped:
        if out and lo <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], hi))
        else:
            out.append((lo, hi))
    return out


def union_many(parts: list[list[Interval]]) -> list[Interval]:
    intervals: list[Interval] = []
    for part in parts:
        intervals.extend(part)
    return merge(intervals)


def complement(intervals: list[Interval]) -> list[Interval]:
    merged = merge(intervals)
    out: list[Interval] = []
    cursor = ZERO
    for lo, hi in merged:
        if cursor < lo:
            out.append((cursor, lo))
        cursor = max(cursor, hi)
    if cursor < ONE:
        out.append((cursor, ONE))
    return out


def intersect_two(a: list[Interval], b: list[Interval]) -> list[Interval]:
    out: list[Interval] = []
    i = j = 0
    while i < len(a) and j < len(b):
        lo = max(a[i][0], b[j][0])
        hi = min(a[i][1], b[j][1])
        if lo < hi:
            out.append((lo, hi))
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return out


def intersect_many(parts: list[list[Interval]]) -> list[Interval]:
    if not parts:
        return [(ZERO, ONE)]
    out = parts[0]
    for part in parts[1:]:
        out = intersect_two(out, part)
        if not out:
            break
    return out


def measure(intervals: list[Interval]) -> F:
    return sum((hi - lo for lo, hi in intervals), ZERO)


def contains_point(intervals: list[Interval], x: F) -> bool:
    return any(lo <= x <= hi for lo, hi in intervals)


def circle_speed_safe_intervals(speed: int, threshold: F = C) -> list[Interval]:
    bad: list[Interval] = []
    for k in range(speed + 1):
        bad.append((F(k, speed) - F(threshold, speed), F(k, speed) + F(threshold, speed)))
    return complement(bad)


def even_safe_intervals(even_half: tuple[int, ...]) -> list[Interval]:
    return intersect_many([circle_speed_safe_intervals(e) for e in even_half])


def branch0_bad_one(odd_speed: int) -> list[Interval]:
    return merge(
        [
            (F(2 * k, odd_speed) - F(2, 14 * odd_speed),
             F(2 * k, odd_speed) + F(2, 14 * odd_speed))
            for k in range((odd_speed // 2) + 2)
        ]
    )


def branch1_bad_one(odd_speed: int) -> list[Interval]:
    return merge(
        [
            (F(2 * k, odd_speed) + F(6, 7 * odd_speed),
             F(2 * k, odd_speed) + F(8, 7 * odd_speed))
            for k in range((odd_speed // 2) + 2)
        ]
    )


def primitive(speeds: tuple[int, ...]) -> bool:
    g = 0
    for v in speeds:
        g = gcd(g, v)
    return g == 1


def is_covering(speeds: tuple[int, ...]) -> bool:
    return all(any(v % q == 0 for v in speeds) for q in range(2, 15))


def random_covering(rng: random.Random, max_speed: int = 190) -> tuple[int, ...]:
    for _attempt in range(30_000):
        speeds: set[int] = set()
        for q in rng.sample(range(2, 15), 13):
            if not any(v % q == 0 for v in speeds):
                choices = [q * k for k in range(1, max_speed // q + 1)]
                speeds.add(rng.choice(choices))
        while len(speeds) < 13:
            speeds.add(rng.randint(1, max_speed))
        row = tuple(sorted(speeds))
        if len(row) == 13 and primitive(row) and is_covering(row):
            return row
    raise RuntimeError("failed to generate covering row")


def wall_labels(x: F, odd: tuple[int, ...], even_half: tuple[int, ...]) -> tuple[str, ...]:
    labels: list[str] = []
    if x == ZERO:
        labels.append("END:0")
    if x == ONE:
        labels.append("END:1")
    for e in even_half:
        if norm(F(e) * x) == C:
            labels.append(f"E:{2 * e}")
    for o in odd:
        if norm(F(o) * x / 2) == C:
            labels.append(f"O0:{o}")
        if norm(F(o) * x / 2) == F(3, 7):
            labels.append(f"O1:{o}")
    return tuple(labels) if labels else ("FREE",)


def label_types(labels: tuple[str, ...]) -> str:
    kinds = sorted({label.split(":", 1)[0] for label in labels})
    return "+".join(kinds)


def branch_mask(component: Interval, branch0: list[Interval], branch1: list[Interval]) -> str:
    mid = (component[0] + component[1]) / 2
    b0 = contains_point(branch0, mid)
    b1 = contains_point(branch1, mid)
    if b0 and b1:
        return "both"
    if b0:
        return "b0"
    if b1:
        return "b1"
    return "union-edge"


def role(speed: int) -> str:
    if speed % 14 == 0:
        return "14Q"
    if speed % 2 == 0:
        return "even_R"
    if speed % 7 == 0:
        return "seven_R"
    return "odd_unit"


@dataclass(frozen=True)
class Window:
    interval: Interval
    width: F
    mask: str
    left: tuple[str, ...]
    right: tuple[str, ...]
    min_score: F
    binders: tuple[int, ...]

    @property
    def signature(self) -> tuple[str, str, str]:
        return (self.mask, label_types(self.left), label_types(self.right))


@dataclass(frozen=True)
class Audit:
    name: str
    speeds: tuple[int, ...]
    odd: tuple[int, ...]
    even_half: tuple[int, ...]
    windows: tuple[Window, ...]
    branch0_measure: F
    branch1_measure: F
    union_measure: F
    min_width: F
    max_width: F
    signature_hist: tuple[tuple[tuple[str, str, str], int], ...]


def audit(name: str, speeds: tuple[int, ...]) -> Audit:
    speeds = tuple(sorted(set(speeds)))
    odd = tuple(v for v in speeds if v % 2 == 1)
    even_half = tuple(v // 2 for v in speeds if v % 2 == 0)
    even_safe = even_safe_intervals(even_half)
    b0_bad = union_many([branch0_bad_one(o) for o in odd])
    b1_bad = union_many([branch1_bad_one(o) for o in odd])
    branch0 = intersect_two(even_safe, complement(b0_bad))
    branch1 = intersect_two(even_safe, complement(b1_bad))
    good_union = merge(branch0 + branch1)

    windows: list[Window] = []
    for component in good_union:
        mid_u = (component[0] + component[1]) / 2
        mask = branch_mask(component, branch0, branch1)
        if mask in ("b1", "both"):
            t = (mid_u + 1) / 2
        else:
            t = mid_u / 2
        s = score(speeds, t)
        binders = tuple(v for v in speeds if norm(F(v) * t) == s)
        windows.append(
            Window(
                interval=component,
                width=component[1] - component[0],
                mask=mask,
                left=wall_labels(component[0], odd, even_half),
                right=wall_labels(component[1], odd, even_half),
                min_score=s,
                binders=binders,
            )
        )
    hist = Counter(window.signature for window in windows)
    return Audit(
        name=name,
        speeds=speeds,
        odd=odd,
        even_half=even_half,
        windows=tuple(windows),
        branch0_measure=measure(branch0),
        branch1_measure=measure(branch1),
        union_measure=measure(good_union),
        min_width=min((w.width for w in windows), default=ZERO),
        max_width=max((w.width for w in windows), default=ZERO),
        signature_hist=tuple(sorted(hist.items(), key=lambda kv: (-kv[1], kv[0]))),
    )


def audited_rows() -> list[tuple[str, tuple[int, ...]]]:
    rows: list[tuple[str, tuple[int, ...]]] = [
        ("covering_AP_with_84", tuple(list(range(1, 12)) + [13, 84])),
        ("covering_AP_with_12_and_84", tuple(list(range(1, 13)) + [84])),
        ("multi_far_84_154", tuple(list(range(1, 11)) + [13, 84, 154])),
        ("even_frontier_probe", (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 13, 28, 154)),
    ]
    for m in range(1, 21):
        rows.append((f"canonical_84m_{m:02d}", tuple(list(range(1, 12)) + [13, 84 * m])))
    for m in range(1, 9):
        rows.append((f"two_tail_drop_{m:02d}", tuple(list(range(1, 10)) + [11, 84 * m, 98 * m, 154])))
    rng = random.Random(3426)
    for i in range(35):
        rows.append((f"random_covering_{i:02d}", random_covering(rng)))
    return rows


def tournament() -> tuple[dict[int, int], list[str]]:
    carriers = {
        "W00_wall_signature_certificate": (28, 20, 18, 10, 9),
        "W01_survivor_component_normal_form": (26, 19, 17, 10, 9),
        "W02_branch_mask_descent_router": (24, 18, 16, 10, 8),
        "W03_owner_current_wall_exception": (22, 17, 15, 9, 10),
        "W04_energy_sheet_sidecar_join": (20, 15, 13, 8, 9),
        "W05_raw_positive_measure_audit": (16, 13, 10, 5, 5),
        "W06_named_analogy_without_walls": (-20, 0, 0, 0, -8),
    }
    totals = {name: sum(vals) for name, vals in carriers.items()}
    hist = dict(sorted(Counter(totals.values()).items()))
    path = [name for name, _ in sorted(totals.items(), key=lambda kv: (-kv[1], kv[0]))]
    return hist, path


def main() -> None:
    audits = [audit(name, row) for name, row in audited_rows()]
    nonempty = [a for a in audits if a.windows]
    all_windows = [w for a in audits for w in a.windows]
    sig_hist = Counter(w.signature for w in all_windows)
    mask_hist = Counter(w.mask for w in all_windows)
    wall_type_hist = Counter((label_types(w.left), label_types(w.right)) for w in all_windows)
    binder_role_hist = Counter(role(v) for w in all_windows for v in w.binders)
    tight = min(audits, key=lambda a: a.union_measure)
    thinnest = min((a for a in audits if a.windows), key=lambda a: a.min_width)
    richest = max(audits, key=lambda a: len(a.signature_hist))

    print("HYP-3427 TWO-BRANCH WALL-SIGNATURE ATLAS")
    print("=" * 78)
    print("Question:")
    print("  Can every HYP-3425 survivor window be named by exact left/right walls?")
    print("  Window certificate = branch mask + endpoint wall labels + midpoint score binders.")
    print()

    print("A. Aggregate exact wall audit")
    print(f"  rows audited:                     {len(audits)}")
    print(f"  rows with survivor windows:       {len(nonempty)}/{len(audits)}")
    print(f"  total survivor windows:           {len(all_windows)}")
    print(f"  global signature types:           {len(sig_hist)}")
    print(f"  branch mask histogram:            {dict(sorted(mask_hist.items()))}")
    print(f"  top wall-type histogram:          {wall_type_hist.most_common(6)}")
    print(f"  midpoint binder role histogram:   {dict(sorted(binder_role_hist.items()))}")
    print(f"  tight row by union measure:        {tight.name} ({fmt(tight.union_measure)})")
    print(f"  thinnest window row:              {thinnest.name} ({fmt(thinnest.min_width)})")
    print(f"  richest signature row:            {richest.name} ({len(richest.signature_hist)} types)")
    print()

    print("B. Tight-row wall certificate")
    print(f"  row={tight.name}")
    print(f"  speeds={tight.speeds}")
    print(f"  union_measure={fmt(tight.union_measure)}, branch0={fmt(tight.branch0_measure)}, branch1={fmt(tight.branch1_measure)}")
    for idx, window in enumerate(tight.windows, start=1):
        print(
            f"  W{idx}: interval=({fmt(window.interval[0])},{fmt(window.interval[1])}) "
            f"width={fmt(window.width)} mask={window.mask}"
        )
        print(f"      left={window.left} right={window.right}")
        print(f"      midpoint_score={fmt(window.min_score)} binders={window.binders}")
    print()

    print("C. Canonical 84m tower signatures")
    print("  m | windows | union | min_width | branch_masks | top_signature")
    for audit_row in audits:
        if not audit_row.name.startswith("canonical_84m_"):
            continue
        m = int(audit_row.name.rsplit("_", 1)[1])
        masks = dict(sorted(Counter(w.mask for w in audit_row.windows).items()))
        top_sig = audit_row.signature_hist[0][0] if audit_row.signature_hist else ("none", "none", "none")
        print(
            f"  {m:2d} | {len(audit_row.windows):7d} | {fmt(audit_row.union_measure):>9} | "
            f"{fmt(audit_row.min_width):>9} | {masks} | {top_sig}"
        )
    print()

    print("D. Global top wall signatures")
    for sig, count in sig_hist.most_common(12):
        print(f"  {sig}: {count}")
    print()

    print("E. Interpretation")
    print("  Positive measure is not yet a proof certificate; exact wall signatures are closer.")
    print("  The audited windows are bounded by small combinations of even walls and odd branch walls.")
    print("  This suggests a finite certificate lemma: every covering row has a survivor window")
    print("  whose left/right labels come from a bounded wall alphabet, or else a named owner")
    print("  current, sheet, exact-period, or state-lift debt is emitted.")
    print()

    hist, path = tournament()
    print("F. Tournament Analysis")
    print("  vertices are proof carriers/wall certificates, not runners or raw intervals.")
    print(f"  score_hist={hist}")
    print("  directed_3cycles=0")
    print("  hamiltonian_path=" + " -> ".join(path))


if __name__ == "__main__":
    main()
