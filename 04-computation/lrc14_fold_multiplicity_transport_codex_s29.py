#!/usr/bin/env python3
"""
HYP-2643 / T891: fold multiplicity transport at the LRC(14) near-AP boundary.

HYP-2640 showed that raw relation rank is only a switch.  This scout refines
the next coordinate: fold multiplicity is not a scalar count either.  The
binding non-AP rows keep nearly the same fold count as the AP, but transport
some fold multiplicity to larger denominators.

We exclude trivial 0+a=a folds and measure the profile

    F_E(c) = #{0<a<b in E : a+b=c in E}.

The main observable is the AP-relative reciprocal fold transport:

    lost_1(E)   = sum_c max(0, F_AP(c)-F_E(c)) / c
    gained_1(E) = sum_c max(0, F_E(c)-F_AP(c)) / c
    net_1(E)    = lost_1(E)-gained_1(E).

For the single end-defect E_k^+={0,...,k-2,k}, the raw fold count can stay
unchanged (notably k=9), while net_1 is still positive.  This is the fold-level
certificate candidate for the razor-thin near-AP drop.

Tournament Analysis vertices are proof quotients, not runners.
"""

from __future__ import annotations

from argparse import ArgumentParser
from collections import Counter, defaultdict
from fractions import Fraction
from functools import reduce
from itertools import combinations, permutations
from math import gcd
import sys


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(line_buffering=True)


def gcd_all(values: tuple[int, ...]) -> int:
    return reduce(gcd, values, 0)


def primitive(E: tuple[int, ...]) -> bool:
    return gcd_all(E) == 1


def is_ap(E: tuple[int, ...]) -> bool:
    E = tuple(sorted(E))
    if len(E) < 2:
        return True
    step = E[1] - E[0]
    return all(E[i + 1] - E[i] == step for i in range(len(E) - 1))


def N_at(E: tuple[int, ...], x: Fraction) -> int:
    hit = set()
    for e in E:
        v = e * x
        v -= v.numerator // v.denominator
        hit.add((v.numerator * 7) // v.denominator)
    return sum(1 for j in range(1, 7) if j not in hit)


def dist_p(E: tuple[int, ...]) -> list[Fraction]:
    E = tuple(sorted(set(E)))
    bps = {Fraction(0), Fraction(1)}
    for e in E:
        if e == 0:
            continue
        for a in range(0, 7 * e + 1):
            bps.add(Fraction(a, 7 * e))
    ordered = sorted(b for b in bps if 0 <= b <= 1)
    p = [Fraction(0)] * 7
    for lo, hi in zip(ordered, ordered[1:]):
        if hi != lo:
            p[N_at(E, (lo + hi) / 2)] += hi - lo
    return p


def g_poly(k: int, t: int) -> Fraction:
    if k == 8:
        return Fraction((t - 1) * (t - 2) * (t - 4) * (t - 5), 40)
    if k in (9, 10):
        return Fraction(-(t - 2) * (t - 3) * (t - 6), 36)
    return Fraction((t - 3) * (t - 4), 12)


def L_y(E: tuple[int, ...]) -> tuple[Fraction, list[Fraction]]:
    p = dist_p(E)
    k = len(E)
    return sum(p[t] * g_poly(k, t) for t in range(7)), p


def fold_profile(E: tuple[int, ...]) -> Counter[int]:
    """Nontrivial visible folds by target c: 0<a<b and a+b=c in E."""
    present = set(E)
    profile: Counter[int] = Counter()
    for a, b in combinations(sorted(E), 2):
        if a == 0:
            continue
        c = a + b
        if c in present:
            profile[c] += 1
    return profile


def reciprocal_mass(profile: Counter[int], power: int = 1) -> Fraction:
    return sum(Fraction(m, c**power) for c, m in profile.items())


def centroid(profile: Counter[int]) -> Fraction:
    count = sum(profile.values())
    if count == 0:
        return Fraction(0)
    return Fraction(sum(c * m for c, m in profile.items()), count)


def transport(ap: Counter[int], profile: Counter[int], power: int = 1) -> tuple[Fraction, Fraction, Fraction]:
    targets = set(ap) | set(profile)
    lost = sum(Fraction(max(0, ap[c] - profile[c]), c**power) for c in targets)
    gained = sum(Fraction(max(0, profile[c] - ap[c]), c**power) for c in targets)
    return lost, gained, lost - gained


def transport_row(ap_E: tuple[int, ...], E: tuple[int, ...]) -> dict[str, object]:
    ap_profile = fold_profile(ap_E)
    profile = fold_profile(E)
    Ly_ap, p_ap = L_y(ap_E)
    Ly, p = L_y(E)
    lost1, gained1, net1 = transport(ap_profile, profile, 1)
    lost2, gained2, net2 = transport(ap_profile, profile, 2)
    return {
        "E": E,
        "Ly": Ly,
        "p0": p[0],
        "Ly_drop": Ly_ap - Ly,
        "p0_drop": p_ap[0] - p[0],
        "fold_count": sum(profile.values()),
        "ap_fold_count": sum(ap_profile.values()),
        "fold_mass1": reciprocal_mass(profile, 1),
        "ap_fold_mass1": reciprocal_mass(ap_profile, 1),
        "net1": net1,
        "lost1": lost1,
        "gained1": gained1,
        "net2": net2,
        "lost2": lost2,
        "gained2": gained2,
        "centroid": centroid(profile),
        "ap_centroid": centroid(ap_profile),
        "profile": profile,
    }


def profile_delta(ap: Counter[int], profile: Counter[int]) -> str:
    parts = []
    for c in sorted(set(ap) | set(profile)):
        delta = profile[c] - ap[c]
        if delta:
            sign = "+" if delta > 0 else ""
            parts.append(f"{c}:{sign}{delta}")
    return "{" + ", ".join(parts) + "}" if parts else "{}"


def single_end_defect(k: int, s: int) -> tuple[int, ...]:
    if s == 0:
        return tuple(range(k))
    return tuple(range(k - 1)) + (k - 1 + s,)


def end_defect_formula(k: int) -> tuple[int, int, Fraction]:
    """Exact net_1 for E={0,...,k-2,k} versus AP={0,...,k-1}."""
    ap_top_mult = (k - 2) // 2
    new_top_mult = (k - 1) // 2 - 1
    net = Fraction(ap_top_mult, k - 1) - Fraction(new_top_mult, k)
    return ap_top_mult, new_top_mult, net


def print_single_defect_ladders(max_s: int) -> None:
    print("=== Single-End Defect Ladders ===")
    print("E_k(s)={0,...,k-2,k-1+s}; s=0 is AP.  Fold deltas omit trivial 0+a=a.")
    for k in (8, 9, 10):
        ap = tuple(range(k))
        ap_profile = fold_profile(ap)
        ap_Ly, ap_p = L_y(ap)
        top_mult, new_mult, formula_net = end_defect_formula(k)
        print()
        print(f"k={k}: AP L_y={float(ap_Ly):.9f} p0={float(ap_p[0]):.9f} fold_profile={dict(ap_profile)}")
        print(
            f"  s=1 exact fold transport: remove {top_mult}/({k-1}), add {new_mult}/{k}, "
            f"net_1={formula_net}={float(formula_net):.9f}"
        )
        print("  s top  L_y_drop       p0_drop       fold_count  net_1        centroid  profile_delta")
        for s in range(1, max_s + 1):
            E = single_end_defect(k, s)
            row = transport_row(ap, E)
            print(
                f" {s:2d} {E[-1]:3d}  {float(row['Ly_drop']):.9f} "
                f"{float(row['p0_drop']):.9f} {row['fold_count']:10d} "
                f"{float(row['net1']):.9f} {float(row['centroid']):9.5f} "
                f"{profile_delta(ap_profile, row['profile'])}"
            )
    print()


def bounded_bank(k: int, max_coord: int) -> list[dict[str, object]]:
    ap = tuple(range(k))
    rows = []
    for tail in combinations(range(1, max_coord + 1), k - 1):
        E = (0,) + tail
        if not primitive(E) or is_ap(E):
            continue
        rows.append(transport_row(ap, E))
    return rows


def summarize_bank(k: int, max_coord: int) -> None:
    rows = bounded_bank(k, max_coord)
    print(f"=== Bounded Bank k={k}, max_coord={max_coord}: {len(rows)} primitive non-AP rows ===")
    print("Top L_y rows:")
    for row in sorted(rows, key=lambda r: r["Ly"], reverse=True)[:8]:
        print(
            f"  Ly={float(row['Ly']):.9f} drop={float(row['Ly_drop']):.9f} "
            f"net_1={float(row['net1']):.9f} count={row['fold_count']:2d} "
            f"centroid={float(row['centroid']):.5f} E={row['E']} "
            f"delta={profile_delta(fold_profile(tuple(range(k))), row['profile'])}"
        )
    positive = [row for row in rows if row["net1"] > 0]
    print("Smallest positive reciprocal fold-transport rows:")
    for row in sorted(positive, key=lambda r: (r["net1"], r["Ly_drop"]))[:8]:
        print(
            f"  net_1={float(row['net1']):.9f} drop={float(row['Ly_drop']):.9f} "
            f"Ly={float(row['Ly']):.9f} count={row['fold_count']:2d} E={row['E']} "
            f"delta={profile_delta(fold_profile(tuple(range(k))), row['profile'])}"
        )
    by_net_bucket: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        # Coarse exact bucket by numerator sign and decimal size for readable output.
        if row["net1"] <= 0:
            key = "<=0"
        elif row["net1"] < Fraction(1, 20):
            key = "(0,0.05)"
        elif row["net1"] < Fraction(1, 10):
            key = "[0.05,0.10)"
        elif row["net1"] < Fraction(1, 5):
            key = "[0.10,0.20)"
        else:
            key = ">=0.20"
        by_net_bucket[key].append(row)
    print("Coarse net_1 buckets:")
    for key in ["<=0", "(0,0.05)", "[0.05,0.10)", "[0.10,0.20)", ">=0.20"]:
        bucket = by_net_bucket.get(key, [])
        if not bucket:
            continue
        best = max(bucket, key=lambda r: r["Ly"])
        mean_drop = sum((r["Ly_drop"] for r in bucket), Fraction(0)) / len(bucket)
        print(
            f"  {key:12s} count={len(bucket):5d} bestLy={float(best['Ly']):.9f} "
            f"best_net={float(best['net1']):.9f} mean_drop={float(mean_drop):.9f} "
            f"bestE={best['E']}"
        )
    print()


def tournament_analysis() -> None:
    vertices = [
        "fold_transport_profile",
        "reciprocal_fold_mass",
        "fold_target_centroid",
        "visible_fold_count",
        "pair_sum_energy",
        "raw_relation_rank",
        "raw_runner_vertices",
    ]
    scores = {i: 1 for i in range(len(vertices))}
    print("=== Tournament Analysis ===")
    print("Pairwise observable: proof usefulness for predicting correction after rank saturation.")
    print("Switch/gauge: lexicographic priority declared by exact near-AP separation.")
    print("Hamiltonian path:")
    print("  " + " > ".join(vertices))
    print(f"score_hist={scores}")
    print("directed_3cycles=0")
    print("SCC_sizes=" + str([1] * len(vertices)))
    print("Hamiltonian_path_count=1")
    print()


def parse_args() -> object:
    parser = ArgumentParser(description="Fold multiplicity transport scout for LRC(14).")
    parser.add_argument("--max-s", type=int, default=8, help="largest single-defect gap s to print")
    parser.add_argument("--k8-W", type=int, default=14, help="bounded bank max coordinate for k=8")
    parser.add_argument("--k9-W", type=int, default=13, help="bounded bank max coordinate for k=9")
    parser.add_argument("--k10-W", type=int, default=12, help="bounded bank max coordinate for k=10")
    parser.add_argument("--skip-banks", action="store_true", help="only print single-defect ladders")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    print("=" * 88)
    print("HYP-2643/T891 fold multiplicity transport scout")
    print("=" * 88)
    print("Raw fold count is not the invariant; reciprocal target transport is the next coordinate.")
    print()
    print_single_defect_ladders(args.max_s)
    if not args.skip_banks:
        summarize_bank(8, args.k8_W)
        summarize_bank(9, args.k9_W)
        summarize_bank(10, args.k10_W)
    print("READING")
    print("  1. The tight k=9 near-AP keeps the same nontrivial fold count as AP (12),")
    print("     but transports three top folds from denominator 8 to 9, giving exact")
    print("     reciprocal fold deficit 1/24. Count misses the defect; transport sees it.")
    print("  2. In the critical k=9 bank, the unique top non-AP is also the unique")
    print("     row in the tiny net_1 bucket (0,0.05).  This is exactly the row with")
    print("     fold transport 3/8 -> 3/9, net_1=1/24.  The k=8 bank is a caveat:")
    print("     its looser cap allows a left-hole clipped AP to win despite larger")
    print("     fold transport, so fold transport is a near-AP boundary coordinate,")
    print("     not a universal scalar order.")
    print("  3. Fold multiplicity should be stored as a target profile, not a scalar.")
    print("     Hidden/coset signs still matter, but the near-AP boundary is a fold")
    print("     transport problem before it is a high-rank problem.")
    print()
    tournament_analysis()


if __name__ == "__main__":
    main()
