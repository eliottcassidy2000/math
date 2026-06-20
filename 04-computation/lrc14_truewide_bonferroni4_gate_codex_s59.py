#!/usr/bin/env python3
"""HYP-2693 / T929: level-4 Bonferroni gate for true-wide LRC14 rows.

Incoming mac-mini S5 showed that the sector-sieve level-4 upper bound fails on
AP-like rows but succeeds on the named true-wide leader.  This scout tests the
split exactly:

    true-wide row (second-largest > 14)  ->  Bonferroni level 4 <= cap_k?

The level-4 upper bound is computed from the missed-sector distribution
p_t = meas{exactly t inner sectors missed}:

    S_r = sum_t binom(t,r) p_t,
    U4 = 1 - S_1 + S_2 - S_3 + S_4.

Because there are exactly six inner sectors, the binomial coefficient
sum_{r=0}^4 (-1)^r binom(t,r) is zero for t=1,2,3,4 and equals 1,5 for
t=5,6.  Thus

    U4 = p_0 + p_5 + 5*p_6.

If U4 <= cap_k on true-wide rows, the remaining rows that need level 6 are
boundary/AP-template rows, matching the DP transfer-template split from
HYP-2691.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import comb
import sys

from lrc14_wide_branch_ridge_codex_s47 import (
    CAP,
    Row,
    additive_energy,
    missed_distribution,
    primitive,
    squarefree_profile,
    sumset_excess,
)

sys.stdout.reconfigure(line_buffering=True)


def fmt(q: Fraction) -> str:
    return f"{q} ({float(q):.9f})"


def bonferroni_upper(dist: tuple[Fraction, ...], level: int = 4) -> Fraction:
    total = Fraction(1)
    sign = -1
    for r in range(1, level + 1):
        sr = sum(Fraction(comb(t, r)) * dist[t] for t in range(r, len(dist)))
        total += sign * sr
        sign *= -1
    return total


def bonferroni4_tail(dist: tuple[Fraction, ...]) -> Fraction:
    return dist[5] + 5 * dist[6]


def assert_bonferroni4_tail_identity(dist: tuple[Fraction, ...]) -> None:
    assert bonferroni_upper(dist, 4) == dist[0] + bonferroni4_tail(dist)


def s_moments(dist: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    return tuple(sum(Fraction(comb(t, r)) * dist[t] for t in range(r, len(dist))) for r in range(7))


@dataclass(frozen=True)
class RowAudit:
    row: Row
    dist: tuple[Fraction, ...]
    upper4: Fraction

    @property
    def k(self) -> int:
        return len(self.row)

    @property
    def cap(self) -> Fraction:
        return CAP[self.k]

    @property
    def p0(self) -> Fraction:
        return self.dist[0]

    @property
    def slack4(self) -> Fraction:
        return self.cap - self.upper4

    @property
    def margin(self) -> Fraction:
        return self.cap - self.p0

    @property
    def tail45(self) -> Fraction:
        return bonferroni4_tail(self.dist)


def push_extreme(bucket: list[RowAudit], audit: RowAudit, keep: int, key) -> None:
    bucket.append(audit)
    bucket.sort(key=key)
    del bucket[keep:]


def scan_box(k: int, bound: int, keep: int = 12) -> dict[str, object]:
    counts = Counter()
    worst_truewide: list[RowAudit] = []
    worst_boundary: list[RowAudit] = []
    risky_p0_truewide: list[RowAudit] = []
    violations: list[RowAudit] = []

    for comb_row in combinations(range(1, bound + 1), k - 1):
        row = (0,) + comb_row
        counts["raw"] += 1
        if row[-1] <= 14 or not primitive(row):
            continue
        counts["primitive_span_gt14"] += 1
        stratum = "truewide" if row[-2] > 14 else "boundary"
        counts[stratum] += 1
        dist = missed_distribution(row)
        assert_bonferroni4_tail_identity(dist)
        audit = RowAudit(row=row, dist=dist, upper4=bonferroni_upper(dist, 4))
        if stratum == "truewide":
            push_extreme(worst_truewide, audit, keep, key=lambda a: (a.slack4, a.margin, a.row))
            push_extreme(risky_p0_truewide, audit, keep, key=lambda a: (-a.p0, a.slack4, a.row))
            if audit.slack4 < 0:
                violations.append(audit)
        else:
            push_extreme(worst_boundary, audit, keep, key=lambda a: (a.slack4, a.margin, a.row))
    return {
        "counts": counts,
        "worst_truewide": worst_truewide,
        "worst_boundary": worst_boundary,
        "risky_p0_truewide": risky_p0_truewide,
        "violations": violations,
    }


def structural(row: Row) -> str:
    gaps = tuple(b - a for a, b in zip(row, row[1:]))
    return (
        f"second={row[-2]} span={row[-1]} maxgap={max(gaps) if gaps else 0} "
        f"exc={sumset_excess(row)} energy={additive_energy(row)} "
        f"sqfree={squarefree_profile(row)}"
    )


def print_audit(prefix: str, audit: RowAudit) -> None:
    moments = s_moments(audit.dist)
    print(
        f"{prefix} E={audit.row} p0={fmt(audit.p0)} cap-p0={fmt(audit.margin)} "
        f"tail45={fmt(audit.tail45)} U4={fmt(audit.upper4)} cap-U4={fmt(audit.slack4)}"
    )
    print(f"      {structural(audit.row)}")
    print(f"      S1..S4={[str(x) for x in moments[1:5]]} missed_dist={[str(x) for x in audit.dist]}")


def named_rows() -> list[tuple[str, Row]]:
    return [
        ("AP9", tuple(range(9))),
        ("wide_doubled_AP_boundary", (0, 2, 4, 6, 8, 10, 12, 14, 16)),
        ("boundary_leader", (0, 2, 4, 6, 8, 10, 12, 14, 15)),
        ("direct_truewide_leader", (0, 4, 6, 8, 10, 12, 14, 15, 16)),
        ("dyadic_block", (0, 1, 2, 4, 8, 12, 16, 20, 24)),
        ("three_cluster_truewide", (0, 1, 2, 30, 31, 32, 60, 61, 62)),
        ("AP_triple_phase_k10", (0, 9, 10, 11, 12, 13, 14, 15, 16, 17)),
    ]


def main() -> None:
    print("HYP-2693 / T929 -- true-wide Bonferroni level-4 gate")
    print("Predicate tested: second-largest > 14 and span > 14 => U4 <= cap_k")
    print("U4 = 1 - S1 + S2 - S3 + S4, exact Fractions from p_t.\n")
    print("Six-sector collapse: U4 = p0 + p5 + 5*p6.  The gate is a high-tail bound.\n")

    print("NAMED ROWS")
    for label, row in named_rows():
        dist = missed_distribution(row)
        assert_bonferroni4_tail_identity(dist)
        audit = RowAudit(row=row, dist=dist, upper4=bonferroni_upper(dist, 4))
        mark = "TRUEWIDE" if row[-2] > 14 else "boundary/AP"
        print_audit(f"  {label:<24} [{mark}]", audit)
    print()

    all_violations: list[tuple[int, int, RowAudit]] = []
    # Staged exact boxes: broad enough to hit the known true-wide ridge, small
    # enough for frequent concurrent checkpoints.  Increase locally for a
    # deeper overnight audit.
    for k, bound in ((8, 18), (9, 18), (10, 16), (11, 15), (12, 15)):
        if k not in CAP:
            continue
        print("=" * 90)
        print(f"BOX k={k}, bound={bound}")
        result = scan_box(k, bound)
        counts: Counter = result["counts"]  # type: ignore[assignment]
        violations: list[RowAudit] = result["violations"]  # type: ignore[assignment]
        print(
            f"  raw={counts['raw']} primitive_span_gt14={counts['primitive_span_gt14']} "
            f"truewide={counts['truewide']} boundary={counts['boundary']} violations={len(violations)}"
        )
        for audit in violations[:10]:
            all_violations.append((k, bound, audit))

        print("  tightest truewide U4 slacks:")
        for idx, audit in enumerate(result["worst_truewide"][:8], 1):  # type: ignore[index]
            print_audit(f"    {idx:2d}.", audit)
        print("  highest truewide p0 rows:")
        for idx, audit in enumerate(result["risky_p0_truewide"][:5], 1):  # type: ignore[index]
            print_audit(f"    {idx:2d}.", audit)
        print("  boundary/AP rows with tightest U4 slacks:")
        for idx, audit in enumerate(result["worst_boundary"][:5], 1):  # type: ignore[index]
            print_audit(f"    {idx:2d}.", audit)
        print()

    print("=" * 90)
    print("SYNTHESIS")
    if all_violations:
        print(f"  Found {len(all_violations)} truewide U4 violations in scanned boxes.")
        for k, bound, audit in all_violations[:12]:
            print_audit(f"    k={k} B={bound}", audit)
    else:
        print("  No truewide U4>cap violations in the scanned exact boxes.")
        print("  Boundary/AP rows still fail U4, including the doubled AP row; they need")
        print("  the finite AP-prefix / low-state-template route from HYP-2691.")
        print("  The proposed proof split is now sharper:")
        print("    truewide high-state branch -> prove U4<=cap by sector moment constraints;")
        print("    boundary/AP low-state branch -> finite AP/dyadic/cube-root templates.")


if __name__ == "__main__":
    main()
