#!/usr/bin/env python3
"""
lrc_torsion_leakage_census_monad_s3.py

monad-explorer-2026-06-06-S3

TORSION-LEAKAGE CORRESPONDENCE for the composite-n Lonely Runner full-cell
system (deepens HYP-1832 / S377 / S367).

MODEL (exact reproduction of S367 lonely_runner_k13_scalar_gauge):
  - n runners, k = n-1 speeds v_1..v_k in Z/n.
  - "cells" = the distinct floor-bin patterns bins_i(alpha)=floor(n*frac(i*alpha))
    as alpha ranges over (0,1); breakpoints a/(n*i).
  - candidate = (shift s in Z/n, cell pattern p).  candidate_count = n*|patterns|.
  - coordinate i with residue r BLOCKS candidate (s,p) iff
        (s*r + bins_i(p)) mod n in {0, n-1}.
    (the speed-i runner is within 1/n of the origin -> not a lonely time)
  - leak(v) = # candidates blocked by NO coordinate.
  - scalar ramp v_i = m*i (the AP / gauge) blocks ALL cells: leak = 0.
    Gauge: adding m*i reparametrizes; we normalize v_1 = 0.

CLAIM under test (the seed):
  the minimum-leak NON-scalar defects are supported in the subgroup
    (n/p) Z/n = ker( Z/n -> Z/(n/p) )   for a prime p | n,
  i.e. the p-torsion that PROJECTS TO ZERO in the (n/p)-runner base.
  Prediction: the global minimum support-1 leak is at residue r = n/p* where
  p* = smallest prime factor of n (the largest torsion generator), at the
  product-sum / resonance coordinates.

Outputs the full support-1 leak census per composite n and tabulates the
minimizers against the torsion prediction.  VALIDATES against S377:
  n=14 coord6 res7 -> missed 56 ;  n=15 -> 120 at coords 6,14 res 5/10.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction
from math import gcd
from sympy import factorint

ONE = Fraction(1, 1)


def cell_bins(n, k, alpha):
    return tuple(int((n * ((i * alpha) % ONE)) // ONE) for i in range(1, k + 1))


def cell_patterns(n, k):
    breaks = {Fraction(0), ONE}
    for i in range(1, k + 1):
        for a in range(n * i + 1):
            breaks.add(Fraction(a, n * i))
    ordered = sorted(breaks)
    out = []
    seen = set()
    for lo, hi in zip(ordered, ordered[1:]):
        if lo == hi:
            continue
        bins = cell_bins(n, k, (lo + hi) / 2)
        if bins in seen:
            continue
        seen.add(bins)
        out.append(bins)
    return out


def build_masks(n):
    """Return (patterns, masks[i][r] as Python int bitmask over candidates, cand_count)."""
    k = n - 1
    patterns = cell_patterns(n, k)
    P = len(patterns)
    cand_count = n * P
    # masks[i][r]: bit (s*P + p_idx) set iff coord i, residue r blocks (s, p)
    masks = [[0] * n for _ in range(k)]
    # exact S367 brute form: coord i, residue r blocks (s,p) iff (s*r+bins_i)%n in {0,n-1}
    for s in range(n):
        for p_idx, bins in enumerate(patterns):
            cidx = s * P + p_idx
            bit = 1 << cidx
            for i in range(k):
                bv = bins[i]
                for r in range(n):
                    if (s * r + bv) % n in (0, n - 1):
                        masks[i][r] |= bit
    return patterns, masks, cand_count


def leak_of(masks, v, cand_count):
    """v: tuple length k of residues. leak = candidates blocked by no coord."""
    blocked = 0
    for i, r in enumerate(v):
        blocked |= masks[i][r]
    return cand_count - bin(blocked).count("1")


def torsion_predict(n):
    """Return dict: residue r -> (prime p, base m=n/p) for each nonzero element of
    (n/p)Z/n that projects to 0 mod (n/p). These are the predicted leak residues."""
    pred = {}
    for p in factorint(n):
        m = n // p  # base; ker(mod m) = mZ/n = {0, m, 2m, ..., (p-1)m}
        for j in range(1, p):
            r = (j * m) % n
            pred.setdefault(r, []).append((p, m))
    return pred


def analyze(n):
    primes = sorted(factorint(n))
    pstar = primes[0]
    patterns, masks, cand_count = build_masks(n)
    k = n - 1
    base_leak = leak_of(masks, tuple([0] * k), cand_count)  # all-zero scalar ramp
    # support-1 census: v = r at coordinate i (i index 0..k-1 => coordinate i+1), others 0
    results = []  # (leak, coord(1-based), residue)
    for i in range(k):
        for r in range(1, n):
            v = [0] * k
            v[i] = r
            lk = leak_of(masks, tuple(v), cand_count)
            results.append((lk, i + 1, r))
    results.sort()
    min_leak = results[0][0]
    minimizers = [(c, r) for (lk, c, r) in results if lk == min_leak]
    # residues achieving the global min
    min_res = sorted(set(r for (c, r) in minimizers))
    pred = torsion_predict(n)
    # per-residue best leak (min over coordinates)
    best_by_res = {}
    for (lk, c, r) in results:
        if r not in best_by_res or lk < best_by_res[r][0]:
            best_by_res[r] = (lk, c)
    return dict(
        n=n, primes=primes, pstar=pstar, cand_count=cand_count, npat=len(patterns),
        base_leak=base_leak, min_leak=min_leak, minimizers=minimizers,
        min_res=min_res, pred=pred, best_by_res=best_by_res,
    )


def main():
    composites = [6, 10, 12, 14, 15, 18, 20, 21, 22]
    print("=" * 78)
    print("TORSION-LEAKAGE CENSUS (S367 full-cell model), support-1 defects")
    print("=" * 78)
    summary = []
    for n in composites:
        A = analyze(n)
        print(f"\n### n={n}  factor={dict(factorint(n))}  p*={A['pstar']}  "
              f"cells={A['npat']} cands={A['cand_count']}  base(all-0) leak={A['base_leak']}")
        print(f"  predicted torsion leak-residues (r: [(p,base)]):")
        for r in sorted(A['pred']):
            print(f"      r={r:3d}  {A['pred'][r]}   (r mod base = "
                  f"{[ (r % m) for (_,m) in A['pred'][r]]})")
        print(f"  GLOBAL min support-1 leak = {A['min_leak']}")
        print(f"     achieved at (coord,res): {A['minimizers'][:12]}"
              + (" ..." if len(A['minimizers']) > 12 else ""))
        print(f"     minimizing residues = {A['min_res']}")
        # Is the global minimizer residue == n/p* (largest torsion gen)?
        npstar = n // A['pstar']
        pred_res = sorted(A['pred'].keys())
        in_torsion = all(r in A['pred'] for r in A['min_res'])
        hits_npstar = npstar in A['min_res']
        print(f"     n/p* = {npstar}.  min-res all torsion? {in_torsion}.  "
              f"n/p* among minimizers? {hits_npstar}")
        # best leak per torsion residue, sorted
        print(f"  best leak per residue (top 10 smallest):")
        ranked = sorted(A['best_by_res'].items(), key=lambda kv: kv[1][0])[:10]
        for r, (lk, c) in ranked:
            tag = "TORSION" + str(A['pred'][r]) if r in A['pred'] else "non-tors"
            print(f"      res={r:3d}  leak={lk:5d}  @coord {c:2d}  [{tag}]")
        summary.append((n, A['min_leak'], A['min_res'], npstar, in_torsion, hits_npstar))

    print("\n" + "=" * 78)
    print("SUMMARY")
    print("=" * 78)
    print(f"{'n':>3} {'minleak':>8} {'min_res':>14} {'n/p*':>5} {'allTors':>8} {'n/p*?':>6}")
    for (n, ml, mr, nps, it, hn) in summary:
        print(f"{n:>3} {ml:>8} {str(mr):>14} {nps:>5} {str(it):>8} {str(hn):>6}")


if __name__ == "__main__":
    main()
