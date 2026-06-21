#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Corrected two-far diagonal-freeze scout for the LRC14 genuine-wide branch.

Context: ``lrc14_threadA_r2_certificate_macmini_0621s7.py`` contains the right
idea for a two-far frozen tail, but its sanity check fails because the freeze
integral is taken over ``x in [0,1)``.  The exact ``p0_fast`` engine uses the
circle coordinate ``y in [0,7)``:

    sector(e,y) = floor(e*y) mod 7.

For far speeds ``f`` and ``f+g`` with ``f -> infinity``, the correct frozen
pair is

    a uniform in Z/7,  second = a + floor(g*y + theta) mod 7,

where averaging over the carrier offset gives weights ``1-frac(g*y)`` and
``frac(g*y)``.  This script:

1. sanity-checks the corrected exact freeze against large finite far averages;
2. scans frozen tails over bounded bases and gap windows using a fast float
   pass;
3. recomputes the top candidates exactly as Fractions and compares them to
   the Thread-A floor ``Q(k-1)``.

It is a scout/proof-order artifact, not a final theorem: the finite low-gap
window and a proof that large gaps cannot exceed the scanned gap cutoff remain
separate obligations.
"""
from __future__ import annotations

import argparse
import functools
import sys
import time
from fractions import Fraction as F
from functools import reduce
from itertools import combinations
from math import floor, gcd

sys.path.insert(0, "04-computation")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")
print = functools.partial(print, flush=True)

from lrc14_signed_multifar_boundary_hierarchy_codex_s51 import boundary_value_direct
from lrc14_threadA_r2_certificate_macmini_0621s7 import p0_fast
from lrc14_wide_branch_ridge_codex_s47 import CAP

ALL_INNER = 0b1111110
QVAL = {k: boundary_value_direct(tuple(range(k - 1)), 1) for k in CAP}


def base_cells_float(base: tuple[int, ...]) -> list[tuple[float, float, tuple[int, ...]]]:
    """Cells where the bounded base has constant missed inner-sector set."""
    breaks = {0.0, 7.0}
    for e in base:
        if e:
            for j in range(1, 7 * e):
                breaks.add(j / e)
    pts = sorted(breaks)
    cells: list[tuple[float, float, tuple[int, ...]]] = []
    for yl, yr in zip(pts, pts[1:]):
        if yr <= yl:
            continue
        ym = (yl + yr) / 2
        mask = 0
        for e in base:
            if e:
                mask |= 1 << (floor(e * ym + 1e-12) % 7)
        missing = tuple(s for s in range(1, 7) if not ((mask >> s) & 1))
        cells.append((yl, yr, missing))
    return cells


def _cover_count(missing: tuple[int, ...], delta_mod: int) -> float:
    if not missing:
        return 1.0
    cnt = 0
    for a in range(7):
        cov = {a % 7, (a + delta_mod) % 7}
        if all(mm in cov for mm in missing):
            cnt += 1
    return cnt / 7.0


COUNT_FLOAT = {
    (missing, d): _cover_count(missing, d)
    for r in range(3)
    for missing in combinations(range(1, 7), r)
    for d in range(7)
}


def dblock7_float(base: tuple[int, ...], g: int) -> float:
    """Fast floating corrected frozen two-far p0 limit."""
    total = 0.0
    for yl, yr, missing in base_cells_float(base):
        t = len(missing)
        width = yr - yl
        if t == 0:
            total += width
            continue
        if t > 2:
            continue
        cur = yl
        while cur < yr - 1e-15:
            d0 = floor(g * cur + 1e-12)
            nxt = min(yr, (d0 + 1) / g)
            if nxt <= cur + 1e-15:
                nxt = yr
            cnt_a = COUNT_FLOAT[(missing, d0 % 7)]
            cnt_b = COUNT_FLOAT[(missing, (d0 + 1) % 7)]
            subw = nxt - cur
            int_phi = g * (nxt * nxt - cur * cur) / 2 - d0 * subw
            total += cnt_a * (subw - int_phi) + cnt_b * int_phi
            cur = nxt
    return total / 7.0


def base_cells_exact(base: tuple[int, ...]) -> list[tuple[F, F, tuple[int, ...]]]:
    breaks = {F(0), F(7)}
    for e in base:
        if e:
            for j in range(1, 7 * e):
                breaks.add(F(j, e))
    pts = sorted(breaks)
    cells: list[tuple[F, F, tuple[int, ...]]] = []
    for yl, yr in zip(pts, pts[1:]):
        if yr <= yl:
            continue
        ym = (yl + yr) / 2
        mask = 0
        for e in base:
            if e:
                mask |= 1 << (int((e * ym).__floor__()) % 7)
        missing = tuple(s for s in range(1, 7) if not ((mask >> s) & 1))
        cells.append((yl, yr, missing))
    return cells


def _cover_count_exact(missing: tuple[int, ...], delta_mod: int) -> F:
    if not missing:
        return F(1)
    cnt = 0
    for a in range(7):
        cov = {a % 7, (a + delta_mod) % 7}
        if all(mm in cov for mm in missing):
            cnt += 1
    return F(cnt, 7)


COUNT_EXACT = {
    (missing, d): _cover_count_exact(missing, d)
    for r in range(3)
    for missing in combinations(range(1, 7), r)
    for d in range(7)
}


def dblock7_exact(base: tuple[int, ...], g: int) -> F:
    """Exact corrected frozen two-far p0 limit."""
    total = F(0)
    for yl, yr, missing in base_cells_exact(base):
        t = len(missing)
        width = yr - yl
        if t == 0:
            total += width
            continue
        if t > 2:
            continue
        cur = yl
        while cur < yr:
            d0 = int((g * cur).__floor__())
            nxt = min(yr, F(d0 + 1, g))
            if nxt <= cur:
                nxt = yr
            cnt_a = COUNT_EXACT[(missing, d0 % 7)]
            cnt_b = COUNT_EXACT[(missing, (d0 + 1) % 7)]
            subw = nxt - cur
            int_phi = g * (nxt * nxt - cur * cur) / 2 - d0 * subw
            total += cnt_a * (subw - int_phi) + cnt_b * int_phi
            cur = nxt
    return total / 7


def finite_far_average(base: tuple[int, ...], g: int, start: int = 1200, span: int = 49) -> F:
    s = F(0)
    n = 0
    for f in range(start, start + span):
        s += p0_fast(tuple(sorted(set(base) | {f, f + g})))
        n += 1
    return s / n


def dense_ok(base: tuple[int, ...], max_gap: int | None) -> bool:
    if max_gap is None:
        return True
    return all((b - a) <= max_gap for a, b in zip(base, base[1:]))


def scan(k: int, gap_max: int, topn: int, dense_diff: int | None) -> list[tuple[float, tuple[int, ...], int]]:
    b = k - 2
    heap: list[tuple[float, tuple[int, ...], int]] = []
    best = -1.0
    checked = 0
    t0 = time.time()
    for rest in combinations(range(1, 15), b - 1):
        base = (0,) + rest
        if not dense_ok(base, dense_diff):
            continue
        cells = base_cells_float(base)
        # Inline variant using precomputed cells for speed.
        for g in range(1, gap_max + 1):
            total = 0.0
            for yl, yr, missing in cells:
                mt = len(missing)
                width = yr - yl
                if mt == 0:
                    total += width
                    continue
                if mt > 2:
                    continue
                cur = yl
                while cur < yr - 1e-15:
                    d0 = floor(g * cur + 1e-12)
                    nxt = min(yr, (d0 + 1) / g)
                    if nxt <= cur + 1e-15:
                        nxt = yr
                    cnt_a = COUNT_FLOAT[(missing, d0 % 7)]
                    cnt_b = COUNT_FLOAT[(missing, (d0 + 1) % 7)]
                    subw = nxt - cur
                    int_phi = g * (nxt * nxt - cur * cur) / 2 - d0 * subw
                    total += cnt_a * (subw - int_phi) + cnt_b * int_phi
                    cur = nxt
            val = total / 7.0
            checked += 1
            if val > best:
                best = val
                print(
                    f"    new float leader k={k} val={val:.9f} base={base} g={g} "
                    f"room_vs_Q={float(QVAL[k]) - val:+.9f}"
                )
            heap.append((val, base, g))
    heap.sort(reverse=True, key=lambda x: x[0])
    print(
        f"  scanned k={k}: candidates={checked}, kept_top={topn}, "
        f"time={time.time()-t0:.1f}s"
    )
    return heap[:topn]


def ceil_fraction(x: F) -> int:
    return (x.numerator + x.denominator - 1) // x.denominator


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gap-max", type=int, default=24)
    ap.add_argument("--topn", type=int, default=12)
    ap.add_argument("--dense-diff", type=int, default=None)
    ap.add_argument("--ks", default="10,11,12", help="comma-separated k values to scan")
    ap.add_argument("--no-sanity", action="store_true")
    args = ap.parse_args()
    ks = tuple(int(x) for x in args.ks.split(",") if x.strip())

    print("=" * 96)
    print("LRC14 genuine-wide corrected two-far diagonal-freeze tail scout (codex-S77)")
    print("=" * 96)
    print("Correct coordinate: y in [0,7), sector(e,y)=floor(e*y) mod 7.")
    print("Old x in [0,1) Dblock undercounts and fails large-f sanity checks.")
    print(f"gap_max={args.gap_max}, topn={args.topn}, dense_diff={args.dense_diff}")
    print()

    if not args.no_sanity:
        print("-" * 96)
        print("Corrected exact freeze vs finite large-f average")
        sanity = [
            ((0, 1, 2, 3, 4, 5, 6, 7), 1),
            ((0, 1, 2, 3, 4, 5, 6, 7), 2),
            ((0, 1, 2, 3, 4, 5, 6, 7, 8), 1),
            ((0, 2, 4, 6, 8, 10, 12, 14), 1),
        ]
        for base, g in sanity:
            exact = dblock7_exact(base, g)
            avg = finite_far_average(base, g)
            print(
                f"  base={base} g={g}: D7={exact} ({float(exact):.9f}) "
                f"large-f-avg={float(avg):.9f} diff={float(avg-exact):+.9e}"
            )
        print()

    for k in ks:
        print("-" * 96)
        print(f"k={k}: Q(k-1)={QVAL[k]} ({float(QVAL[k]):.9f})")
        top = scan(k, args.gap_max, args.topn, args.dense_diff)
        print("  exact top candidates:")
        for i, (val, base, g) in enumerate(top, 1):
            exact = dblock7_exact(base, g)
            room = QVAL[k] - exact
            if room > 0:
                tail_f = ceil_fraction(F(7, 1) / room)
                tail_note = f"D7+7/f<Q once f>{tail_f}"
            else:
                tail_note = "D7 >= Q: frozen tail itself obstructs Q-floor"
            print(
                f"    {i:02d}. D7={exact} ({float(exact):.9f}) "
                f"float={val:.9f} room={room} ({float(room):+.9f}) "
                f"base={base} g={g} :: {tail_note}"
            )
        print()


if __name__ == "__main__":
    main()
