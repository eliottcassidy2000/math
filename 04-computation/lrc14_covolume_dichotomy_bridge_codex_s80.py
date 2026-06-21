#!/usr/bin/env python3
"""Exact covolume/peel-depth audit for the LRC14 dichotomy seam.

codex S80.

Incoming mac-mini S23 reframes the remaining structural seam as:

    high coverage above the decorrelated floor should force small relation
    covolume, hence low-complexity remove-one/peel structure.

This script tests that statement on the exact bounded-span genuine-wide rows
used in S78.  It is not a proof.  It is a finite atlas diagnostic that asks
whether the rows with positive `p0-Q(k-1)` are explained by a low relation
depth after removing a small number of elements.

Metrics:
  * full primitive norm squared, the saturated covolume proxy `sum offsets^2`
    after affine normalization;
  * best one-peel and two-peel primitive norm/span;
  * peel depth: least number of removed elements after which the affine
    normalized span is <=14.

Tournament Analysis:
  vertices are proof buckets `(peel_depth, two_peel_span bucket, far_count)`;
  pairwise observable is exact bucket maximum `p0`;
  switch/gauge orients A -> B iff max_p0(A) >= max_p0(B);
  the Hamiltonian path is the displayed bucket order.

Assumption challenge: tournament vertices are not runners.  They are relation
complexity buckets.  This preserves the predicate "can this row beat Q or cap"
and destroys phase data and the exact far placement, so any positive result
must be fed back into endpoint/R-tail ledgers before becoming a proof.
"""
from __future__ import annotations

import argparse
import functools
import math
import statistics
import sys
import time
from collections import Counter
from fractions import Fraction
from itertools import combinations
from typing import Iterable

print = functools.partial(print, flush=True)

sys.path.insert(0, "04-computation")

from lrc14_genuinewide_generalized_doublet_exact_codex_s78 import (  # noqa: E402
    far_part,
    is_genuine_wide,
    primitive,
)
from lrc14_threadA_regime_dichotomy_kpswf8 import CAP, QVAL, p0_fast  # noqa: E402


def affine_normalize(E: Iterable[int]) -> tuple[int, ...]:
    """Translate to min 0 and divide by the gcd of all offsets."""
    S = tuple(sorted(set(int(x) for x in E)))
    if not S:
        return ()
    m = S[0]
    offsets = [x - m for x in S]
    g = 0
    for x in offsets:
        g = math.gcd(g, abs(x))
    if g == 0:
        return (0,)
    return tuple(x // g for x in offsets)


def affine_span(E: Iterable[int]) -> int:
    N = affine_normalize(E)
    return N[-1] - N[0] if N else 0


def primitive_norm2(E: Iterable[int]) -> int:
    N = affine_normalize(E)
    return sum(x * x for x in N)


def best_peel(E: tuple[int, ...], r: int) -> tuple[int, int, tuple[int, ...]]:
    """Return `(best_span, best_norm2, removed)` over all r-element removals."""
    best: tuple[int, int, tuple[int, ...]] | None = None
    for rem in combinations(E, r):
        rem_set = set(rem)
        sub = tuple(x for x in E if x not in rem_set)
        if len(sub) < 2:
            continue
        span = affine_span(sub)
        norm2 = primitive_norm2(sub)
        row = (span, norm2, tuple(rem))
        if best is None or row < best:
            best = row
    if best is None:
        return (10**9, 10**18, ())
    return best


def peel_depth(E: tuple[int, ...], max_depth: int = 4) -> tuple[int, int, int, tuple[int, ...]]:
    """Least removal depth that affine-reduces the row to span <= 14."""
    if affine_span(E) <= 14:
        return (0, affine_span(E), primitive_norm2(E), ())
    for r in range(1, max_depth + 1):
        span, norm2, rem = best_peel(E, r)
        if span <= 14:
            return (r, span, norm2, rem)
    span, norm2, rem = best_peel(E, max_depth)
    return (max_depth + 1, span, norm2, rem)


def bucket(depth: int, span2: int, far_count: int) -> tuple[str, str, str]:
    if depth <= 2:
        d = f"depth{depth}"
    else:
        d = "depth3plus"
    if span2 <= 10:
        s = "two_span_le10"
    elif span2 <= 14:
        s = "two_span_11_14"
    elif span2 <= 18:
        s = "two_span_15_18"
    else:
        s = "two_span_gt18"
    r = f"r{far_count}" if far_count <= 3 else "r4plus"
    return (d, s, r)


def pearson(xs: list[float], ys: list[float]) -> float:
    if len(xs) < 2:
        return 0.0
    mx = statistics.fmean(xs)
    my = statistics.fmean(ys)
    vx = sum((x - mx) ** 2 for x in xs)
    vy = sum((y - my) ** 2 for y in ys)
    if vx == 0 or vy == 0:
        return 0.0
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / math.sqrt(vx * vy)


def row_summary(row: dict[str, object], cap: Fraction, q: Fraction) -> str:
    E = row["E"]
    p0 = row["p0"]
    far = far_part(E)
    return (
        f"E={E}, p0={p0}, p0-Q={p0 - q}, cap-p0={cap - p0}, "
        f"far={far}, depth={row['depth']}, rem={row['rem']}, "
        f"span1={row['span1']}, span2={row['span2']}, "
        f"norm_full={row['norm_full']}, norm1={row['norm1']}, norm2={row['norm2']}, "
        f"bucket={row['bucket']}"
    )


def audit(k: int, span: int, topn: int) -> dict[str, object]:
    start = time.time()
    cap = CAP[k]
    q = QVAL[k]
    rows: list[dict[str, object]] = []
    total = primitive_count = gw_count = over_q = 0
    by_bucket: dict[tuple[str, str, str], tuple[Fraction, dict[str, object]]] = {}

    for comb in combinations(range(1, span + 1), k - 1):
        total += 1
        E = (0,) + comb
        if not primitive(E):
            continue
        primitive_count += 1
        if not is_genuine_wide(E):
            continue
        gw_count += 1
        p0 = p0_fast(E)
        if p0 > q:
            over_q += 1
        span1, norm1, rem1 = best_peel(E, 1)
        span2, norm2, rem2 = best_peel(E, 2)
        depth, dspan, dnorm, drem = peel_depth(E)
        row = {
            "E": E,
            "p0": p0,
            "norm_full": primitive_norm2(E),
            "span1": span1,
            "norm1": norm1,
            "rem1": rem1,
            "span2": span2,
            "norm2": norm2,
            "rem2": rem2,
            "depth": depth,
            "dspan": dspan,
            "dnorm": dnorm,
            "rem": drem,
            "bucket": bucket(depth, span2, len(far_part(E))),
        }
        rows.append(row)
        b = row["bucket"]
        if b not in by_bucket or p0 > by_bucket[b][0]:
            by_bucket[b] = (p0, row)

    profile_order = sorted(by_bucket.items(), key=lambda item: (-item[1][0], item[0]))
    score_hist = Counter()
    nprof = len(profile_order)
    for i, _item in enumerate(profile_order):
        score_hist[nprof - 1 - i] += 1

    over_rows = [r for r in rows if r["p0"] > q]
    xs = [float(r["p0"] - q) for r in rows]
    pos_xs = [float(r["p0"] - q) for r in over_rows]

    def inv_sqrt(key: str, source: list[dict[str, object]]) -> list[float]:
        return [1.0 / math.sqrt(float(r[key])) for r in source if float(r[key]) > 0]

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
        "rows": rows,
        "over_rows": over_rows,
        "top": sorted(rows, key=lambda r: r["p0"], reverse=True)[:topn],
        "profile_order": profile_order,
        "score_hist": score_hist,
        "corr_all_full": pearson(xs, inv_sqrt("norm_full", rows)),
        "corr_all_one": pearson(xs, inv_sqrt("norm1", rows)),
        "corr_all_two": pearson(xs, inv_sqrt("norm2", rows)),
        "corr_pos_full": pearson(pos_xs, inv_sqrt("norm_full", over_rows)),
        "corr_pos_one": pearson(pos_xs, inv_sqrt("norm1", over_rows)),
        "corr_pos_two": pearson(pos_xs, inv_sqrt("norm2", over_rows)),
    }


SPAN20_K12_OVER_Q = [
    (0, 2, 4, 6, 8, 9, 10, 11, 12, 14, 16, 18),
    (0, 2, 4, 6, 7, 8, 9, 10, 12, 14, 16, 18),
    (0, 2, 3, 4, 6, 8, 9, 10, 12, 14, 15, 18),
    (0, 2, 4, 6, 7, 8, 10, 11, 12, 14, 18, 20),
    (0, 2, 3, 4, 6, 8, 9, 10, 12, 14, 16, 18),
    (0, 2, 4, 5, 6, 7, 8, 10, 12, 14, 16, 20),
    (0, 2, 4, 6, 7, 8, 10, 11, 12, 14, 16, 20),
]


def witness_row(E: tuple[int, ...], k: int) -> dict[str, object]:
    p0 = p0_fast(E)
    span1, norm1, rem1 = best_peel(E, 1)
    span2, norm2, rem2 = best_peel(E, 2)
    depth, dspan, dnorm, drem = peel_depth(E)
    row = {
        "E": E,
        "p0": p0,
        "norm_full": primitive_norm2(E),
        "span1": span1,
        "norm1": norm1,
        "rem1": rem1,
        "span2": span2,
        "norm2": norm2,
        "rem2": rem2,
        "depth": depth,
        "dspan": dspan,
        "dnorm": dnorm,
        "rem": drem,
        "bucket": bucket(depth, span2, len(far_part(E))),
    }
    # Keep k in the dict for downstream ad hoc parsing without changing row_summary.
    row["k"] = k
    return row


def print_witness_bank() -> None:
    k = 12
    cap = CAP[k]
    q = QVAL[k]
    print("\n" + "=" * 92)
    print("S78 span<=20 k=12 positive p0-Q witness-bank peel audit")
    print("=" * 92)
    print(
        "These are the seven over-Q rows reported by the S78 span<=20 exact audit; "
        "this targeted bank avoids rerunning the larger span<=20 covolume scan."
    )
    rows = [witness_row(E, k) for E in SPAN20_K12_OVER_Q]
    print(
        "Witness depth histogram: "
        + ", ".join(
            f"d{d}:{c}" for d, c in sorted(Counter(r["depth"] for r in rows).items())
        )
    )
    print(
        "Witness two-peel span: "
        + ", ".join(
            f"{b}:{c}"
            for b, c in sorted(
                Counter("le14" if r["span2"] <= 14 else "gt14" for r in rows).items()
            )
        )
    )
    for rank, row in enumerate(sorted(rows, key=lambda r: r["p0"], reverse=True), start=1):
        print(f"  {rank}. {row_summary(row, cap, q)}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ks", nargs="+", type=int, default=[10, 11, 12])
    parser.add_argument("--span", type=int, default=18)
    parser.add_argument("--topn", type=int, default=8)
    parser.add_argument(
        "--witness-bank-only",
        action="store_true",
        help="only print the S78 span<=20 k=12 positive p0-Q witness-bank peel audit",
    )
    args = parser.parse_args()

    if args.witness_bank_only:
        print_witness_bank()
        return

    print("=" * 92)
    print("S80 covolume/peel-depth bridge for the LRC14 dichotomy seam")
    print("=" * 92)
    print("Exact finite audit only.  Positive p0-Q rows are the structural test cases.")

    for k in args.ks:
        result = audit(k, args.span, args.topn)
        cap = result["cap"]
        q = result["q"]
        print("\n" + "-" * 92)
        print(
            f"k={k}, span<={args.span}: total={result['total']}, "
            f"primitive={result['primitive_count']}, genuine_wide={result['gw_count']}, "
            f"over_Q={result['over_q']}, seconds={result['seconds']:.2f}"
        )
        print(f"cap={cap}, Q(k-1)={q}, cap-Q={cap - q}")
        print(
            "Pearson over all rows: "
            f"p0-Q vs 1/full_norm={result['corr_all_full']:+.4f}, "
            f"1/one_peel_norm={result['corr_all_one']:+.4f}, "
            f"1/two_peel_norm={result['corr_all_two']:+.4f}"
        )
        if result["over_q"]:
            print(
                "Pearson over positive p0-Q rows: "
                f"p0-Q vs 1/full_norm={result['corr_pos_full']:+.4f}, "
                f"1/one_peel_norm={result['corr_pos_one']:+.4f}, "
                f"1/two_peel_norm={result['corr_pos_two']:+.4f}"
            )
        depth_hist = Counter(r["depth"] for r in result["rows"])
        over_depth_hist = Counter(r["depth"] for r in result["over_rows"])
        two_span_over = Counter(
            "le14" if r["span2"] <= 14 else "gt14" for r in result["over_rows"]
        )
        print(
            "Depth histogram all: "
            + ", ".join(f"d{d}:{c}" for d, c in sorted(depth_hist.items()))
        )
        print(
            "Depth histogram over_Q: "
            + (", ".join(f"d{d}:{c}" for d, c in sorted(over_depth_hist.items())) or "none")
        )
        print(
            "Over_Q two-peel span: "
            + (", ".join(f"{b}:{c}" for b, c in sorted(two_span_over.items())) or "none")
        )

        print("Top rows by p0:")
        for rank, row in enumerate(result["top"], start=1):
            print(f"  {rank}. {row_summary(row, cap, q)}")

        print("Bucket tournament Hamiltonian path by max p0:")
        for rank, (b, (_p0, row)) in enumerate(result["profile_order"][: args.topn], start=1):
            print(f"  {rank}. {b}: {row_summary(row, cap, q)}")
        print(
            "Bucket tournament score histogram: "
            + ", ".join(
                f"score{score}:{count}"
                for score, count in sorted(result["score_hist"].items())
            )
        )

    print("\nInterpretation target:")
    print(
        "If every positive p0-Q row has depth<=2 or two-peel span<=14, the failed "
        "one-peel dichotomy can be replaced by a relation-depth dichotomy: "
        "single-far periodicity handles depth 1, while depth 2 belongs to the "
        "generalized-doublet/R-tail atlas."
    )
    print_witness_bank()


if __name__ == "__main__":
    main()
