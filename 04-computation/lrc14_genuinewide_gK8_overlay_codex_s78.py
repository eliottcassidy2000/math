#!/usr/bin/env python3
"""Exact gK8 overlay for the S78 genuine-wide bounded-span audit.

The S78 generalized-doublet audit ranks genuine-wide rows by exact `p0=q0`.
Incoming HYP-2809 suggests a sharper proof order: use the Delsarte dual

    L_yK8 = 10*q0 + q3 + 10*q6

on the missed-sector distribution `q_t`.  Since `10*q0 <= L_yK8`, the cap
route follows from `L_yK8 <= 10*cap_k`.

This script reuses the same normalized primitive genuine-wide domain as
`lrc14_genuinewide_generalized_doublet_exact_codex_s78.py`, but ranks rows by
`L_yK8`.  It is still a finite-window audit, not a proof.

Tournament analysis:
  vertices: profile buckets (far_count, gap bucket, parity profile);
  pairwise observable: exact bucket maximum `L_yK8`;
  switch/gauge: orient A -> B iff max_L_yK8(A) >= max_L_yK8(B);
  tie Hamiltonian path: bucket labels sorted by (-max_L_yK8, label).
"""
from __future__ import annotations

import argparse
import functools
import sys
import time
from collections import Counter
from fractions import Fraction
from itertools import combinations

print = functools.partial(print, flush=True)

sys.path.insert(0, "04-computation")

from lrc14_gK8_wide_check_claudeopus_0622 import miss_dist
from lrc14_genuinewide_generalized_doublet_exact_codex_s78 import (
    far_part,
    is_genuine_wide,
    primitive,
    profile,
)
from lrc14_threadA_regime_dichotomy_kpswf8 import CAP, QVAL


def lyk8_from_q(q: list[Fraction]) -> Fraction:
    return 10 * q[0] + q[3] + 10 * q[6]


def row_summary(E: tuple[int, ...], q: list[Fraction], cap: Fraction) -> str:
    ly = lyk8_from_q(q)
    far = far_part(E)
    gap = far[1] - far[0] if len(far) == 2 else None
    gap_txt = f", gap={gap}" if gap is not None else ""
    return (
        f"E={E}, L_yK8={ly}, 10cap-L={10 * cap - ly}, "
        f"p0={q[0]}, cap-p0={cap - q[0]}, "
        f"q3={q[3]}, q6={q[6]}, far={far}{gap_txt}, profile={profile(E)}"
    )


def audit(k: int, span: int, topn: int) -> dict[str, object]:
    start = time.time()
    cap = CAP[k]
    qfloor = QVAL[k]
    total = primitive_count = gw_count = ly_violations = over_q = 0
    by_r: dict[int, tuple[Fraction, tuple[int, ...], list[Fraction]]] = {}
    by_profile: dict[
        tuple[str, str, str], tuple[Fraction, tuple[int, ...], list[Fraction]]
    ] = {}
    top_ly: list[tuple[Fraction, tuple[int, ...], list[Fraction]]] = []
    top_p0: list[tuple[Fraction, tuple[int, ...], list[Fraction]]] = []

    for comb in combinations(range(1, span + 1), k - 1):
        total += 1
        E = (0,) + comb
        if not primitive(E):
            continue
        primitive_count += 1
        if not is_genuine_wide(E):
            continue
        gw_count += 1
        q = miss_dist(E)
        ly = lyk8_from_q(q)
        if ly > 10 * cap:
            ly_violations += 1
        if q[0] > qfloor:
            over_q += 1

        r = len(far_part(E))
        if r not in by_r or ly > by_r[r][0]:
            by_r[r] = (ly, E, q)
        pr = profile(E)
        if pr not in by_profile or ly > by_profile[pr][0]:
            by_profile[pr] = (ly, E, q)

        top_ly.append((ly, E, q))
        top_ly.sort(reverse=True)
        if len(top_ly) > topn:
            top_ly.pop()

        top_p0.append((q[0], E, q))
        top_p0.sort(reverse=True)
        if len(top_p0) > topn:
            top_p0.pop()

    profile_order = sorted(
        by_profile.items(), key=lambda item: (-item[1][0], item[0])
    )
    score_hist = Counter()
    nprof = len(profile_order)
    for i, _item in enumerate(profile_order):
        score_hist[nprof - 1 - i] += 1

    return {
        "k": k,
        "span": span,
        "seconds": time.time() - start,
        "total": total,
        "primitive_count": primitive_count,
        "gw_count": gw_count,
        "over_q": over_q,
        "ly_violations": ly_violations,
        "cap": cap,
        "qfloor": qfloor,
        "by_r": by_r,
        "profile_order": profile_order,
        "score_hist": score_hist,
        "top_ly": top_ly,
        "top_p0": top_p0,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--spans", nargs="+", type=int, default=[18, 20])
    parser.add_argument("--ks", nargs="+", type=int, default=[10, 11, 12])
    parser.add_argument("--topn", type=int, default=8)
    args = parser.parse_args()

    print("=" * 88)
    print("S78 exact gK8 overlay on genuine-wide bounded-span rows")
    print("=" * 88)
    print("Finite-window audit only.  Target inequality: L_yK8 <= 10*cap_k.")
    print(
        "Tournament: vertices=profile buckets; edge A->B iff max_L_yK8(A)>=max_L_yK8(B)."
    )

    for span in args.spans:
        print("\n" + "#" * 88)
        print(f"SPAN <= {span}")
        print("#" * 88)
        for k in args.ks:
            result = audit(k, span, args.topn)
            cap = result["cap"]
            qfloor = result["qfloor"]
            print("\n" + "-" * 88)
            print(
                f"k={k}, span<={span}: total={result['total']}, "
                f"primitive={result['primitive_count']}, genuine_wide={result['gw_count']}, "
                f"over_Q={result['over_q']}, L_yK8_violations={result['ly_violations']}, "
                f"seconds={result['seconds']:.2f}"
            )
            print(f"cap={cap}, 10cap={10 * cap}, Q(k-1)={qfloor}")
            print("Best L_yK8 by far count:")
            for r in sorted(result["by_r"]):
                _ly, E, q = result["by_r"][r]
                print(f"  r={r}: {row_summary(E, q, cap)}")

            best_ly, best_E, best_q = result["top_ly"][0]
            print(f"GLOBAL L_yK8 BEST: {row_summary(best_E, best_q, cap)}")
            print(f"generalized_doublet_Lbest = {len(far_part(best_E)) == 2}")

            print("Top rows by L_yK8:")
            for rank, (_ly, E, q) in enumerate(result["top_ly"], start=1):
                print(f"  {rank}. {row_summary(E, q, cap)}")

            print("Top rows by p0, with L_yK8 overlay:")
            for rank, (_p0, E, q) in enumerate(result["top_p0"], start=1):
                print(f"  {rank}. {row_summary(E, q, cap)}")

            print("Profile tournament Hamiltonian path by max L_yK8:")
            for rank, (pr, (_ly, E, q)) in enumerate(
                result["profile_order"][: args.topn], start=1
            ):
                print(f"  {rank}. {pr}: {row_summary(E, q, cap)}")
            print(
                "Profile tournament score histogram: "
                + ", ".join(
                    f"score{score}:{count}"
                    for score, count in sorted(result["score_hist"].items())
                )
            )


if __name__ == "__main__":
    main()
