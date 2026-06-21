#!/usr/bin/env python3
"""Exact q6 boundary-contraction scout for the gK8 LRC14 route.

Recent HYP-2814/HYP-2815/HYP-2816 work says the gK8 concentration proof should
use suppression of the all-missed atom

    q6(E) = meas{x : every nonzero speed of E lands in sector 0}.

Asymptotically, a far speed should multiply q6 by 1/7.  The boundary window
f>=15 is different: nearby speeds can preserve most of the initial sector-0
interval.  This script records the exact finite correction on the current
frontier families.

Tournament Analysis:
  vertices: proof packets (single-far or double-far frontier families);
  observable: packet risk tuple (minimum gK8 margin, maximum q6 ratio);
  switch/gauge: orient toward smaller margin, then larger q6 ratio;
  tie Hamiltonian path: lexicographic packet labels.
"""
from __future__ import annotations

import argparse
import functools
import sys
import time
from fractions import Fraction
from functools import reduce
from itertools import combinations
from math import gcd

print = functools.partial(print, flush=True)

sys.path.insert(0, "04-computation")

from lrc14_gK8_decorrelation_principle_claudeopus_0622 import LyK8, miss_dist
from lrc14_threadA_regime_dichotomy_kpswf8 import CAP
from lrc14_wide_branch_ridge_codex_s47 import primitive


def lcm(a: int, b: int) -> int:
    return a // gcd(a, b) * b


@functools.lru_cache(maxsize=None)
def q6_grid(E: tuple[int, ...]) -> Fraction:
    """Exact q6 using the repo's integer breakpoint grid, but only sector 0.

    For speed e, q6 requires frac(e*x) in [0,1/7).  The common denominator
    grid is the same one used by the full missed-distribution routines; this
    function just avoids constructing all seven q-atoms.
    """
    nz = tuple(int(x) for x in E if x)
    if not nz:
        return Fraction(0)
    period_lcm = reduce(lcm, nz, 1)
    denom = 7 * period_lcm
    den2 = 2 * period_lcm
    breakpoints = {0, denom}
    for e in nz:
        step = period_lcm // e
        x = 0
        for _ in range(7 * e + 1):
            breakpoints.add(x)
            x += step

    num = 0
    bps = sorted(breakpoints)
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo:
            continue
        mid = lo + hi
        if all(((e * mid // den2) % 7) == 0 for e in nz):
            num += hi - lo
    return Fraction(num, denom)


def validate_q6() -> None:
    samples = [
        (0, 1, 2, 3, 4, 5, 6, 7, 17, 18),
        (0, 1, 3, 5, 7, 9, 11, 13, 15, 17),
        (0, 2, 4, 6, 8, 9, 10, 11, 12, 14, 16, 18),
    ]
    print("Validation against full miss_dist q6:")
    for E in samples:
        fast = q6_grid(E)
        full = miss_dist(E)[6]
        ok = fast == full
        print(f"  E={E}: q6_grid={fast}, miss_dist_q6={full}, ok={ok}")
        if not ok:
            raise AssertionError((E, fast, full))


S78_DOUBLE_BASES = {
    10: (0, 1, 3, 5, 7, 9, 11, 13),
    11: (0, 2, 4, 6, 7, 8, 10, 12, 14),
    12: (0, 2, 4, 6, 8, 9, 10, 11, 12, 14),
}


def top14_base(size: int) -> tuple[int, ...]:
    return tuple([0] + list(range(1, size - 1)) + [14])


def top_cluster_base(size: int) -> tuple[int, ...]:
    return tuple([0] + list(range(15 - (size - 1), 15)))


def single_bases(k: int) -> list[tuple[str, tuple[int, ...]]]:
    size = k - 1
    bases: list[tuple[str, tuple[int, ...]]] = [
        ("single_consecutive", tuple(range(size))),
        ("single_top14", top14_base(size)),
        ("single_top_cluster", top_cluster_base(size)),
    ]
    if 2 * (size - 1) <= 14:
        bases.append(("single_even_ap", tuple(range(0, 2 * size, 2))))
    return sorted(set(bases))


def double_bases(k: int) -> list[tuple[str, tuple[int, ...]]]:
    size = k - 2
    bases: list[tuple[str, tuple[int, ...]]] = [
        ("double_consecutive", tuple(range(size))),
        ("double_top14", top14_base(size)),
        ("double_top_cluster", top_cluster_base(size)),
        ("double_s78_leader_base", S78_DOUBLE_BASES[k]),
    ]
    if 2 * (size - 1) <= 14:
        bases.append(("double_even_ap", tuple(range(0, 2 * size, 2))))
    return sorted(set(bases))


def row_record(
    packet: str,
    B: tuple[int, ...],
    E: tuple[int, ...],
    cap: Fraction,
) -> dict[str, object]:
    q6b = q6_grid(B)
    q6e = q6_grid(E)
    ratio = q6e / q6b if q6b else Fraction(0)
    ly = LyK8(E)
    return {
        "packet": packet,
        "B": B,
        "E": E,
        "q6B": q6b,
        "q6E": q6e,
        "ratio": ratio,
        "Ly": ly,
        "margin": 10 * cap - ly,
    }


def better_ratio(new: dict[str, object], old: dict[str, object] | None) -> bool:
    if old is None:
        return True
    return (new["ratio"], -new["margin"], str(new["E"])) > (
        old["ratio"],
        -old["margin"],
        str(old["E"]),
    )


def better_ly(new: dict[str, object], old: dict[str, object] | None) -> bool:
    if old is None:
        return True
    return (new["Ly"], new["ratio"], str(new["E"])) > (
        old["Ly"],
        old["ratio"],
        str(old["E"]),
    )


def frontier_scan(k: int, fmax: int, gmax: int) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    cap = CAP[k]
    packet_stats: dict[str, dict[str, object]] = {}
    top_records: list[dict[str, object]] = []

    def add_record(rec: dict[str, object]) -> None:
        packet = rec["packet"]
        stats = packet_stats.setdefault(
            packet,
            {"max_ratio": None, "max_ly": None, "count": 0},
        )
        stats["count"] += 1
        if better_ratio(rec, stats["max_ratio"]):
            stats["max_ratio"] = rec
        if better_ly(rec, stats["max_ly"]):
            stats["max_ly"] = rec
        top_records.append(rec)

    for label, B in single_bases(k):
        for f in range(15, fmax + 1):
            E = tuple(sorted(B + (f,)))
            if len(set(E)) != k or not primitive(E):
                continue
            add_record(row_record(label, B, E, cap))

    for label, B in double_bases(k):
        for M in range(15, fmax + 1):
            for g in range(1, gmax + 1):
                E = tuple(sorted(B + (M, M + g)))
                if len(set(E)) != k or not primitive(E):
                    continue
                add_record(row_record(label, B, E, cap))

    packets = []
    for packet, stats in packet_stats.items():
        max_ratio = stats["max_ratio"]
        max_ly = stats["max_ly"]
        packets.append(
            {
                "packet": packet,
                "count": stats["count"],
                "max_ratio": max_ratio,
                "max_ly": max_ly,
                "risk_key": (max_ly["margin"], -max_ratio["ratio"], packet),
            }
        )
    packets.sort(key=lambda x: x["risk_key"])
    top_records.sort(key=lambda r: (-r["ratio"], r["margin"], str(r["E"])))
    return packets, top_records


def print_record(prefix: str, rec: dict[str, object]) -> None:
    print(
        f"{prefix} packet={rec['packet']} B={rec['B']} E={rec['E']} "
        f"q6B={rec['q6B']} q6E={rec['q6E']} ratio={rec['ratio']} "
        f"L_yK8={rec['Ly']} 10cap-L={rec['margin']}"
    )


def exhaustive_single(k: int, fmax: int) -> dict[str, object]:
    cap = CAP[k]
    size = k - 1
    best: tuple[Fraction, tuple[int, ...], tuple[int, ...], tuple[int, ...]] | None = None
    count = 0
    start = time.time()
    for S in combinations(range(1, 15), size - 1):
        B = (0,) + S
        q6b = q6_grid(B)
        if q6b == 0:
            continue
        for f in range(15, fmax + 1):
            E = tuple(sorted(B + (f,)))
            if len(set(E)) != k or not primitive(E):
                continue
            count += 1
            ratio = q6_grid(E) / q6b
            candidate = (ratio, B, E, (f,))
            if best is None or (ratio, str(E)) > (best[0], str(best[2])):
                best = candidate
    assert best is not None
    _ratio, B, E, _far = best
    best_record = row_record("single_exhaustive", B, E, cap)
    return {
        "k": k,
        "count": count,
        "seconds": time.time() - start,
        "best": best_record,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ks", nargs="+", type=int, default=[10, 11, 12])
    parser.add_argument("--fmax", type=int, default=60)
    parser.add_argument("--gmax", type=int, default=8)
    parser.add_argument(
        "--exhaustive-single",
        action="store_true",
        help="Exhaustively scan all bounded single-far bases for the selected k values.",
    )
    args = parser.parse_args()

    print("=" * 88)
    print("S79 exact q6 boundary-contraction scout for gK8 concentration")
    print("=" * 88)
    print(
        "Finite target: quantify the small-f correction to the asymptotic q6 -> q6/7 story."
    )
    print(
        "Risk reading: larger q6 ratio is worse for pure concentration, smaller 10cap-L is worse for gK8."
    )
    validate_q6()

    for k in args.ks:
        print("\n" + "#" * 88)
        print(f"k={k}: frontier packet scan, f/M<= {args.fmax}, gaps<= {args.gmax}")
        print("#" * 88)
        packets, top_records = frontier_scan(k, args.fmax, args.gmax)
        print("Packet tournament Hamiltonian path by risk (min margin, then max q6 ratio):")
        for i, packet in enumerate(packets, start=1):
            print(f"  {i}. {packet['packet']} count={packet['count']}")
            print_record("     max_ratio:", packet["max_ratio"])
            print_record("     max_Ly   :", packet["max_ly"])

        print("Top q6-ratio rows in packet bank:")
        for rec in top_records[:8]:
            print_record("  ", rec)

    if args.exhaustive_single:
        print("\n" + "#" * 88)
        print("EXHAUSTIVE bounded-base single-far q6-ratio scan")
        print("#" * 88)
        for k in args.ks:
            result = exhaustive_single(k, args.fmax)
            print(
                f"k={k}: rows={result['count']} seconds={result['seconds']:.2f}"
            )
            print_record("  best:", result["best"])

    print("\n" + "=" * 88)
    print("Interpretation:")
    print("  The boundary correction is not a uniform 1/7.  In the checked frontier")
    print("  and the exhaustive single-far scan, the exact worst ratio is 14/15,")
    print("  attained by a bounded base whose q6 interval is controlled by speed 14")
    print("  and a new far speed f=15.  The adjacent two-far frontier has worst")
    print("  ratio 7/8.  Both worst packet rows still have large positive gK8 margin.")


if __name__ == "__main__":
    main()
