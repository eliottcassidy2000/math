#!/usr/bin/env python3
"""Exact bounded-span audit for the generalized-doublet LRC14 route.

codex S78, HYP-2810 scaffold.

Incoming HYP-2807 reframes the mac-mini k=12 odd-bridge obstruction as a
spread doublet: the far part has exactly two elements {M, M+g}.  The published
script for HYP-2807 is useful but samples bases.  This script supplies an exact
bounded-span audit over all normalized rows

    E = {0} union S, |E|=k, max(E)<=W, primitive(E), genuine-wide(E).

It is not a proof of the unbounded theorem.  It is a rigor guardrail: if the
generalized-doublet route is correct, exact span windows should keep their
maxima at far_count r=2, and the dangerous over-Q rows should stay finite and
direct-cap safe.

Tournament analysis:
  vertices: profile buckets (far_count, far_gap bucket, parity profile);
  pairwise observable: exact bucket maximum p0;
  switch/gauge: orient A -> B iff max_p0(A) >= max_p0(B), lexicographic tie;
  tie Hamiltonian path: bucket labels sorted by (-max_p0, label).
"""
from __future__ import annotations

import argparse
import functools
import sys
import time
from collections import Counter, defaultdict
from fractions import Fraction
from functools import reduce
from itertools import combinations
from math import gcd

print = functools.partial(print, flush=True)

sys.path.insert(0, "04-computation")

from lrc14_threadA_regime_dichotomy_kpswf8 import CAP, QVAL, p0_fast


def primitive(E: tuple[int, ...]) -> bool:
    return reduce(gcd, E) == 1


def reprim(E: tuple[int, ...]) -> tuple[int, ...]:
    g = reduce(gcd, E)
    return tuple(x // g for x in E) if g > 1 else E


def is_genuine_wide(E: tuple[int, ...]) -> bool:
    E = tuple(sorted(set(E)))
    if len(E) < 2 or E[0] != 0 or max(E) <= 14 or not primitive(E):
        return False
    for e in E:
        sub = tuple(x for x in E if x != e)
        if len(sub) < 2:
            return False
        sub = reprim(sub)
        if max(sub) - min(sub) <= 14:
            return False
    return True


def far_part(E: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(x for x in E if x > 14)


def far_blocks(far: tuple[int, ...]) -> int:
    if not far:
        return 0
    blocks = 1
    for a, b in zip(far, far[1:]):
        if b != a + 1:
            blocks += 1
    return blocks


def parity_profile(E: tuple[int, ...]) -> str:
    evens = sum(1 for x in E if x % 2 == 0)
    odds = len(E) - evens
    far = far_part(E)
    far_odds = sum(1 for x in far if x % 2 == 1)
    return f"even{evens}_odd{odds}_farodd{far_odds}"


def gap_bucket(far: tuple[int, ...]) -> str:
    if len(far) != 2:
        return "not_r2"
    gap = far[1] - far[0]
    if gap <= 3:
        return f"g{gap}"
    if gap <= 7:
        return "g4_7"
    return "g8plus"


def profile(E: tuple[int, ...]) -> tuple[str, str, str]:
    far = far_part(E)
    return (f"r{len(far)}", gap_bucket(far), parity_profile(E))


def short_row(E: tuple[int, ...], value: Fraction, cap: Fraction, q: Fraction) -> str:
    far = far_part(E)
    gap = far[1] - far[0] if len(far) == 2 else None
    gap_txt = f", gap={gap}" if gap is not None else ""
    return (
        f"E={E}, p0={value}, cap-p0={cap - value}, "
        f"p0-Q={value - q}, far={far}{gap_txt}, blocks={far_blocks(far)}, "
        f"profile={profile(E)}"
    )


def audit(k: int, span: int, topn: int) -> dict[str, object]:
    start = time.time()
    cap = CAP[k]
    q = QVAL[k]
    total = primitive_count = gw_count = over_q = 0
    by_r: dict[int, tuple[Fraction, tuple[int, ...]]] = {}
    by_profile: dict[tuple[str, str, str], tuple[Fraction, tuple[int, ...]]] = {}
    far_gap_counter: Counter[int] = Counter()
    top: list[tuple[Fraction, tuple[int, ...]]] = []

    for comb in combinations(range(1, span + 1), k - 1):
        total += 1
        E = (0,) + comb
        if not primitive(E):
            continue
        primitive_count += 1
        if not is_genuine_wide(E):
            continue
        gw_count += 1
        far = far_part(E)
        if len(far) == 2:
            far_gap_counter[far[1] - far[0]] += 1
        value = p0_fast(E)
        if value > q:
            over_q += 1
        r = len(far)
        if r not in by_r or value > by_r[r][0]:
            by_r[r] = (value, E)
        pr = profile(E)
        if pr not in by_profile or value > by_profile[pr][0]:
            by_profile[pr] = (value, E)
        top.append((value, E))
        top.sort(reverse=True)
        if len(top) > topn:
            top.pop()

    profile_order = sorted(
        by_profile.items(), key=lambda item: (-item[1][0], item[0])
    )
    score_hist = Counter()
    nprof = len(profile_order)
    for i, _item in enumerate(profile_order):
        # Transitive tournament score under the max-p0 orientation.
        score_hist[nprof - 1 - i] += 1

    return {
        "k": k,
        "span": span,
        "seconds": time.time() - start,
        "total": total,
        "primitive_count": primitive_count,
        "gw_count": gw_count,
        "over_q": over_q,
        "cap": cap,
        "q": q,
        "by_r": by_r,
        "by_profile": by_profile,
        "profile_order": profile_order,
        "score_hist": score_hist,
        "far_gap_counter": far_gap_counter,
        "top": top,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--spans", nargs="+", type=int, default=[18, 20, 22])
    parser.add_argument("--ks", nargs="+", type=int, default=[10, 11, 12])
    parser.add_argument("--topn", type=int, default=8)
    args = parser.parse_args()

    print("=" * 88)
    print("HYP-2810 exact bounded-span generalized-doublet audit")
    print("=" * 88)
    print("This is an exact finite-window audit, not an unbounded proof.")
    print(
        "Tournament: vertices=profile buckets; edge A->B iff max_p0(A)>=max_p0(B); "
        "Hamiltonian path is the displayed profile order."
    )

    for span in args.spans:
        print("\n" + "#" * 88)
        print(f"SPAN <= {span}")
        print("#" * 88)
        for k in args.ks:
            result = audit(k, span, args.topn)
            cap = result["cap"]
            q = result["q"]
            print("\n" + "-" * 88)
            print(
                f"k={k}, span<={span}: total={result['total']}, "
                f"primitive={result['primitive_count']}, genuine_wide={result['gw_count']}, "
                f"over_Q={result['over_q']}, seconds={result['seconds']:.2f}"
            )
            print(f"cap={cap}, Q(k-1)={q}, cap-Q={cap - q}")
            print("Best by far count:")
            for r in sorted(result["by_r"]):
                value, E = result["by_r"][r]
                print(f"  r={r}: {short_row(E, value, cap, q)}")
            best_value, best_E = result["top"][0]
            print(f"GLOBAL BEST: {short_row(best_E, best_value, cap, q)}")
            print(f"generalized_doublet_best = {len(far_part(best_E)) == 2}")

            print("Top rows:")
            for rank, (value, E) in enumerate(result["top"], start=1):
                print(f"  {rank}. {short_row(E, value, cap, q)}")

            print("Far-gap counts for r=2:")
            gap_items = sorted(result["far_gap_counter"].items())
            print("  " + ", ".join(f"g={g}:{c}" for g, c in gap_items[:18]))
            if len(gap_items) > 18:
                print(f"  ... {len(gap_items) - 18} more gap values")

            print("Profile tournament Hamiltonian path:")
            for rank, (pr, (value, E)) in enumerate(
                result["profile_order"][: args.topn], start=1
            ):
                print(f"  {rank}. {pr}: max_p0={value}, witness={E}")
            print(
                "Profile tournament score histogram: "
                + ", ".join(
                    f"score{score}:{count}"
                    for score, count in sorted(result["score_hist"].items())
                )
            )


if __name__ == "__main__":
    main()
