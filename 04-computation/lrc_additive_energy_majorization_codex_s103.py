#!/usr/bin/env python3
"""
LRC14 additive-energy majorization attack (codex S103).

Goal: sharpen HYP-2885's additive-energy/Fejer extremality route.

Test three candidate theorem shapes on exact bounded banks:

1. Scalar additive energy:
   A(E)=#{a+b=c+d}.  This is known to be maximized by intervals, but a proof
   of p0(E)<=p0(AP) cannot rely on scalar monotonicity if higher A can have
   lower p0 than a lower-A row.

2. Difference-profile majorization:
   Let d_E(h)=#{ unordered pairs {a,b}: |a-b|=h }.  The interval profile
   (k-1,k-2,...,1) should majorize every k-set's sorted difference profile.
   This is the right Fejer object: |Ehat|^2 = k + 2 sum_h d_E(h) cos(2*pi*h*x).

3. Local compressions:
   Replace one element e by e-1 when the slot is empty.  If such moves always
   increased p0 when they increase the Fejer profile, they would give a direct
   compression proof.  Counterexamples route the proof to labelled low-mode
   packets rather than naive compression.

Tournament Analysis discipline:
vertices are proof carriers (scalar energy, AP diff-majorization, pairwise
diff-majorization, local compression, Ly certificate), not runners.  Pairwise
observable is implication reliability on exact banks.  The orientation points
from fewer counterexamples and stronger AP-facing theorem shape to weaker ones.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from math import comb, gcd
from functools import reduce
import sys

sys.stdout.reconfigure(line_buffering=True)


SECTORS = tuple((F(j, 7), F(j + 1, 7)) for j in range(1, 7))


def primitive(E: tuple[int, ...]) -> bool:
    return reduce(gcd, (e for e in E if e), 0) == 1


def fast_profile(E: tuple[int, ...]) -> tuple[F, list[F]]:
    """Return p0 and p[0..6], where p[t]=meas{t inner sectors are empty}.

    E is an offset set containing 0.  Sector 0 is hit by the offset 0, so p0 is
    exactly the seven-sector cover atom for the normalized row.
    """
    E = tuple(sorted(set(E)))
    bp = {F(0), F(1)}
    for e in E:
        if e == 0:
            continue
        for j in range(1, 7):
            for end in (F(j, 7), F(j + 1, 7)):
                m = 0
                while True:
                    xv = (end + m) / e
                    if xv >= 1:
                        break
                    if xv >= 0:
                        bp.add(xv)
                    m += 1
    bps = sorted(b for b in bp if F(0) <= b < F(1))
    nz = [e for e in E if e]
    law = [F(0)] * 7
    for lo, hi in zip(bps, bps[1:] + [F(1)]):
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        hit = set()
        for e in nz:
            fr = (e * mid) % 1
            hit.add((fr.numerator * 7) // fr.denominator)
        free = len([j for j in range(1, 7) if j not in hit])
        law[free] += hi - lo
    return law[0], law


DUAL = {
    8: ([F(1), F(-1), F(1), F(-9, 10), F(3, 5)], 4),
    9: ([F(1), F(-13, 18), F(4, 9), F(-1, 6)], 3),
    10: ([F(1), F(-13, 18), F(4, 9), F(-1, 6)], 3),
}


def ly_from_law(law: list[F], k: int) -> F:
    y, rmax = DUAL[k]
    return sum(y[r] * sum(F(comb(t, r)) * law[t] for t in range(7)) for r in range(rmax + 1))


def additive_energy(E: tuple[int, ...]) -> int:
    sums = Counter(a + b for a in E for b in E)
    return sum(v * v for v in sums.values())


def diff_profile(E: tuple[int, ...]) -> tuple[int, ...]:
    counts = Counter(abs(a - b) for a, b in combinations(E, 2))
    return tuple(sorted(counts.values(), reverse=True))


def majorizes(x: tuple[int, ...], y: tuple[int, ...]) -> bool:
    n = max(len(x), len(y))
    sx = sy = 0
    for i in range(n):
        if i < len(x):
            sx += x[i]
        if i < len(y):
            sy += y[i]
        if sx < sy:
            return False
    return sum(x) == sum(y)


def bank(k: int, b: int) -> list[tuple[int, ...]]:
    rows = []
    for rest in combinations(range(1, b + 1), k - 1):
        E = (0,) + rest
        if primitive(E):
            rows.append(E)
    return rows


def collect(k: int, b: int) -> list[dict]:
    rows = []
    for E in bank(k, b):
        p0, law = fast_profile(E)
        rows.append(
            {
                "E": E,
                "A": additive_energy(E),
                "diff": diff_profile(E),
                "p0": p0,
                "Ly": ly_from_law(law, k),
                "law": law,
            }
        )
    return rows


def fmtq(x: F) -> str:
    return f"{float(x):.6f} ({x})"


def scalar_inversion(rows: list[dict], metric: str) -> tuple[int, tuple | None]:
    """Count envelope inversions A_high>A_low but metric_high<metric_low."""
    grouped: dict[int, list[dict]] = {}
    for r in rows:
        grouped.setdefault(r["A"], []).append(r)
    best_lower = None
    count = 0
    worst = None
    for A in sorted(grouped):
        group = grouped[A]
        if best_lower is not None:
            for r in group:
                if r[metric] < best_lower[metric]:
                    count += 1
                    gap = best_lower[metric] - r[metric]
                    if worst is None or gap > worst[0]:
                        worst = (gap, best_lower, r)
        best_here = max(group, key=lambda r: r[metric])
        if best_lower is None or best_here[metric] > best_lower[metric]:
            best_lower = best_here
    return count, worst


def pairwise_diff_majorization_violations(rows: list[dict], metric: str, limit: int = 8):
    total = 0
    examples = []
    for i, a in enumerate(rows):
        for j, b in enumerate(rows):
            if i == j:
                continue
            if majorizes(a["diff"], b["diff"]) and a[metric] < b[metric]:
                total += 1
                if len(examples) < limit:
                    examples.append((a, b, b[metric] - a[metric]))
    return total, examples


def compression_audit(rows_by_E: dict[tuple[int, ...], dict], metric: str):
    checked = 0
    energy_up = 0
    metric_down = 0
    profile_up = 0
    profile_up_metric_down = 0
    examples = []
    for E, r in rows_by_E.items():
        Eset = set(E)
        for e in E:
            if e == 0 or e - 1 in Eset:
                continue
            B = tuple(sorted((Eset - {e}) | {e - 1}))
            if B not in rows_by_E:
                continue
            checked += 1
            s = rows_by_E[B]
            if s["A"] > r["A"]:
                energy_up += 1
                if s[metric] < r[metric]:
                    metric_down += 1
                    if len(examples) < 6:
                        examples.append((r, s, r[metric] - s[metric]))
            if majorizes(s["diff"], r["diff"]) and s["diff"] != r["diff"]:
                profile_up += 1
                if s[metric] < r[metric]:
                    profile_up_metric_down += 1
                    if len(examples) < 6:
                        examples.append((r, s, r[metric] - s[metric]))
    return {
        "checked": checked,
        "energy_up": energy_up,
        "metric_down": metric_down,
        "profile_up": profile_up,
        "profile_up_metric_down": profile_up_metric_down,
        "examples": examples,
    }


def rank_carriers(stats: dict[str, int]) -> list[tuple[int, str]]:
    # Lower score is better; converted below to a tournament-like order.
    raw = {
        "AP_diff_majorization": stats["ap_diff_fail"],
        "Ly_certificate": stats["ly_over_ap"],
        "pairwise_diff_majorization": stats["pairwise_diff_p0_viol"],
        "local_profile_compression": stats["profile_comp_down"],
        "scalar_additive_energy": stats["scalar_p0_inv"],
    }
    return sorted((v, k) for k, v in raw.items())


def main() -> None:
    print("LRC14 additive-energy majorization attack (codex S103)")
    print("=" * 78)
    configs = [(8, 13), (9, 13), (10, 12)]
    aggregate_stats = {
        "ap_diff_fail": 0,
        "ly_over_ap": 0,
        "pairwise_diff_p0_viol": 0,
        "profile_comp_down": 0,
        "scalar_p0_inv": 0,
    }

    for k, b in configs:
        rows = collect(k, b)
        rows_by_E = {r["E"]: r for r in rows}
        ap = tuple(range(k))
        ap_row = rows_by_E[ap]
        print()
        print(f"--- k={k}, primitive bank E=(0 plus {k-1} from 1..{b}), rows={len(rows)}")
        print(f"AP p0={fmtq(ap_row['p0'])}  Ly={fmtq(ap_row['Ly'])}  A={ap_row['A']}")
        print(f"AP diff profile={ap_row['diff']}")

        ap_fail = [r for r in rows if not majorizes(ap_row["diff"], r["diff"])]
        p0_over_ap = [r for r in rows if r["p0"] > ap_row["p0"]]
        ly_over_ap = [r for r in rows if r["Ly"] > ap_row["Ly"]]
        aggregate_stats["ap_diff_fail"] += len(ap_fail)
        aggregate_stats["ly_over_ap"] += len(ly_over_ap)
        print(f"AP diff-profile majorizes all rows: failures={len(ap_fail)}")
        print(f"AP p0 max in bank: over_AP={len(p0_over_ap)}")
        print(f"AP Ly max in bank: over_AP={len(ly_over_ap)}")

        for metric in ("p0", "Ly"):
            inv_count, worst = scalar_inversion(rows, metric)
            if metric == "p0":
                aggregate_stats["scalar_p0_inv"] += inv_count
            print(f"scalar A monotonicity inversions for {metric}: {inv_count}")
            if worst:
                gap, low, high = worst
                print(
                    f"  worst {metric}: lower-A {low['E']} A={low['A']} {metric}={fmtq(low[metric])}"
                )
                print(
                    f"                 higher-A {high['E']} A={high['A']} {metric}={fmtq(high[metric])}"
                )
                print(f"                 gap={fmtq(gap)}")

        # Pairwise full-profile majorization is stronger than scalar A.  Use the
        # smaller exact banks above so the quadratic audit is still cheap.
        pair_viol, examples = pairwise_diff_majorization_violations(rows, "p0")
        aggregate_stats["pairwise_diff_p0_viol"] += pair_viol
        print(f"pairwise diff-profile majorization => p0 monotone violations: {pair_viol}")
        for a, b2, gap in examples[:3]:
            print(
                f"  diff-majorized row {a['E']} p0={fmtq(a['p0'])} < {b2['E']} p0={fmtq(b2['p0'])}; gap={fmtq(gap)}"
            )

        comp = compression_audit(rows_by_E, "p0")
        aggregate_stats["profile_comp_down"] += comp["profile_up_metric_down"]
        print(
            "one-step compression audit:"
            f" checked={comp['checked']} energy_up={comp['energy_up']} energy_up_p0_down={comp['metric_down']}"
            f" profile_up={comp['profile_up']} profile_up_p0_down={comp['profile_up_metric_down']}"
        )
        for old, new, gap in comp["examples"][:3]:
            print(
                f"  compression {old['E']} -> {new['E']}:"
                f" A {old['A']}->{new['A']}, p0 {fmtq(old['p0'])}->{fmtq(new['p0'])}, down={fmtq(gap)}"
            )

    print()
    print("Tournament Analysis over proof carriers")
    print("=" * 78)
    ranked = rank_carriers(aggregate_stats)
    print("counterexample/obligation scores (lower is better):")
    for score, name in ranked:
        print(f"  {name:30s} {score}")
    print(
        "orientation: AP_diff_majorization -> Ly_certificate -> "
        "pairwise_diff_majorization -> local_profile_compression -> scalar_additive_energy"
    )
    print(
        "readout: scalar A is a correlation heuristic, not a monotone theorem. "
        "The live proof object is AP-facing difference-profile/Fejer majorization "
        "plus labelled sector/Fourier packets to handle p0 signs."
    )


if __name__ == "__main__":
    main()
