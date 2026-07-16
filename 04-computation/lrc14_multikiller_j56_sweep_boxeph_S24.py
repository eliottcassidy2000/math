#!/usr/bin/env python3
"""
THE j = 5, 6 MULTI-KILLER SWEEP (the decidable compact-core sweep, leg 1)
boxeph-2026-07-16-S24 -- independent implementation of the THM-883 (mac-mini
fragmentation lemma) box sweep, completing the j = 5, 6 residual their S114
launched but did not report, and cross-validating j = 2, 3, 4.

SETTING (all derived independently from the Fragmentation Lemma, radius 1/13):
  S = P u W, P = S cap {1..12}, W = {w_1 < ... < w_j} outliers >= 13.
  If M(S) < 1/13 then the open arc systems B_w (arcs of half-width 1/(13w) at
  the w-grid) of the REMAINING killers must cover the closed good set G_i of
  P u {w_1..w_i} at every stage i, giving the NECESSARY bounds
     (stage bound)  w_{i+1} <= 2m / ((13-2m) * lmax(G_i)),   m = j-i remaining
     (last stage)   every component of G_{j-1} sits INSIDE one open arc of
                    w_j  =>  w_j < 2/(13 * lmax(G_{j-1})).
  LRC(<=13) guarantees every G_i is nonempty; moreover |P u W_i| = 13-m speeds
  => M >= 1/(14-m) > 1/13, so lmax(G_i) > 0 and every bound is finite:
  the configuration space is a FINITE box. Sweeping it exactly and finding no
  M < 1/13 config proves: every 13-set (covering or not) with j in {2..6}
  outliers has M >= 1/13.

SOUNDNESS PROTOCOL (mirrors mac-mini's fast sweep):
  - in-tree good sets are computed in floats with INWARD-conservative rounding
    (arcs expanded by EPS): the tracked good set is a SUBSET of the true one,
    so lmax is an UNDER-estimate and every w-range an OVER-enumeration;
  - every leaf configuration is re-verified from scratch in EXACT Fractions
    (emptiness of the closed good set at 1/13);
  - float endpoints are rationals with denominator <= 13*wmax ~ 1e5, whose
    pairwise gaps are >> EPS = 1e-9: the conservative rounding can never flip
    interval order, only directions at exact ties.

Pure Python, no numpy.
"""

import sys
from fractions import Fraction as Fr
from itertools import combinations
from math import floor, gcd

EPS = 1e-9
RADIUS_NUM, RADIUS_DEN = 1, 13          # radius 1/13

# ------------------------------------------------------------ float side

def subtract_arcs_f(gs, w):
    """gs: sorted list of (lo,hi) float good intervals in [0,1].
    Subtract open arcs of speed w at radius 1/13, arcs EXPANDED by EPS
    (=> resulting good set is inward-conservative)."""
    r = 1.0 / (13.0 * w) + EPS
    out = []
    for lo, hi in gs:
        # arcs with centers k/w intersecting [lo,hi]
        k0 = int(floor(lo * w))          # center candidates k0-1 .. k1+1
        k1 = int(floor(hi * w)) + 1
        cur = lo
        for k in range(k0 - 1, k1 + 2):
            c = k / w
            alo, ahi = c - r, c + r
            if ahi <= cur:
                continue
            if alo >= hi:
                break
            if alo > cur:
                out.append((cur, alo))
            cur = max(cur, ahi)
            if cur >= hi:
                break
        if cur < hi:
            out.append((cur, hi))
    return [(a, b) for a, b in out if b - a > 0.0]

def good_set_f(speeds):
    gs = [(0.0, 1.0)]
    for v in sorted(speeds):
        gs = subtract_arcs_f(gs, v)
        if not gs:
            return gs
    return gs

def lmax_f(gs):
    return max((b - a for a, b in gs), default=0.0)

# ------------------------------------------------------------ exact side

def subtract_arcs_q(gs, w):
    """exact: subtract OPEN arcs of speed w radius 1/13 from closed intervals.
    Closed good set: keep endpoints (closed complement of open arcs)."""
    r = Fr(RADIUS_NUM, RADIUS_DEN * w)
    out = []
    for lo, hi in gs:
        k0 = int(lo * w)
        k1 = int(hi * w) + 1
        cur = lo
        for k in range(k0 - 1, k1 + 2):
            c = Fr(k, w)
            alo, ahi = c - r, c + r
            if ahi <= cur:
                continue
            if alo >= hi:
                break
            if alo >= cur:
                out.append((cur, alo))     # closed remnant [cur, alo]
            cur = max(cur, ahi)
            if cur > hi:
                break
        if cur <= hi:
            out.append((cur, hi))
    # degenerate points [x,x] are legitimate members of the CLOSED good set
    return out

def good_set_q(speeds):
    gs = [(Fr(0), Fr(1))]
    for v in sorted(speeds):
        gs = subtract_arcs_q(gs, v)
        if not gs:
            return gs
    return gs

def exact_nonempty(speeds):
    """closed good set at 1/13 nonempty? (M >= 1/13)"""
    return bool(good_set_q(speeds))

def is_covering(S):
    return all(any(v % q == 0 for v in S) for q in range(2, 15))

# ------------------------------------------------------------ the tree

class Sweep:
    def __init__(self, j, wcap_hard=10**7, shard=0, nshards=1):
        self.j = j
        self.branches = 0
        self.leaves = 0
        self.violations = []
        self.noncov_below = []
        self.wcap_hard = wcap_hard
        self.maxdepth_w = 0
        self.shard = shard
        self.nshards = nshards

    def run(self):
        j = self.j
        shapes = [P for idx, P in enumerate(combinations(range(1, 13), 13 - j))
                  if idx % self.nshards == self.shard]
        for si, P in enumerate(shapes):
            gs = good_set_f(P)
            self.recurse(list(P), [], gs, 12, j)
            if (si + 1) % 20 == 0:
                print(f"    [j={j} shard {self.shard}/{self.nshards}] "
                      f"{si+1}/{len(shapes)} shapes, {self.branches} branches, "
                      f"{self.leaves} leaves", flush=True)
        return self

    def recurse(self, P, W, gs, w_prev, m):
        """m killers remaining; gs = inward-conservative good set of P u W."""
        if m == 0:
            self.leaf(P, W)
            return
        ell = lmax_f(gs)
        if ell <= 4 * EPS:
            # float set degenerate -- fall back to exact lmax (rare)
            eq = good_set_q(P + W)
            ell_exact = max((float(b - a) for a, b in eq), default=0.0)
            if ell_exact == 0.0:
                # good set is isolated points; components have zero length:
                # the last-stage swallow bound is vacuous ONLY if m>1 arcs
                # must still cover points -- but stage bound needs ell>0.
                # LRC(<=13) guarantees this cannot happen (M >= 1/(14-m) > 1/13
                # gives a positive-length component); flag loudly if reached.
                print(f"!! ZERO-LENGTH lmax at P={P} W={W} -- theory violated?",
                      flush=True)
                return
            ell = ell_exact * (1 - 1e-9)
            gs = [(float(a), float(b)) for a, b in eq]   # refresh components
        if m == 1:
            hi = 2.0 / (13.0 * ell * (1 - 1e-9))
        else:
            hi = 2.0 * m / ((13 - 2 * m) * ell * (1 - 1e-9))
        hi = int(hi + 1e-9)
        if hi > self.wcap_hard:
            print(f"!! w-cap exceeded at P={P} W={W}: hi={hi}", flush=True)
            hi = self.wcap_hard
        lo = max(13, w_prev + 1)
        if m == 1:
            # FAST last stage: w swallows the good set iff every component
            # [a,b] fits inside one open arc: exists integer k in
            # (w*b - 1/13, w*a + 1/13). Conservative widening by DELTA, then
            # exact confirmation in the leaf. O(#components) per w.
            DELTA = 1e-6
            comps = gs
            for w in range(lo, hi + 1):
                self.branches += 1
                ok = True
                for a, b in comps:
                    kmin = floor(w * b - 1.0 / 13.0 - DELTA) + 1
                    if kmin > w * a + 1.0 / 13.0 + DELTA:
                        ok = False
                        break
                if ok:
                    self.leaf(P, W + [w])
                if w > self.maxdepth_w:
                    self.maxdepth_w = w
            return
        for w in range(lo, hi + 1):
            self.branches += 1
            gs2 = subtract_arcs_f(gs, w)
            self.recurse(P, W + [w], gs2, w, m - 1)
            if w > self.maxdepth_w:
                self.maxdepth_w = w

    def leaf(self, P, W):
        """candidate config: float says killers might cover G at 1/13.
        Exact confirmation from scratch."""
        self.leaves += 1
        S = P + W
        if not exact_nonempty(S):
            entry = (tuple(P), tuple(W), is_covering(S))
            if entry[2]:
                self.violations.append(entry)
                print(f"  ** COVERING VIOLATION (M < 1/13): P={P} W={W}",
                      flush=True)
            else:
                self.noncov_below.append(entry)
                print(f"  (non-covering config with M < 1/13: P={P} W={W})",
                      flush=True)

def main():
    args = sys.argv[1:]
    shard, nshards = 0, 1
    if "--shard" in args:
        i = args.index("--shard")
        shard, nshards = int(args[i + 1]), int(args[i + 2])
        args = args[:i] + args[i + 3:]
    js = [int(x) for x in args] or [2, 3, 4, 5, 6]
    print(f"SHARD {shard}/{nshards}", flush=True)
    print("j | shapes | branches | float-leaves | covering M<1/13 | noncov M<1/13 | max w",
          flush=True)
    for j in js:
        sw = Sweep(j, shard=shard, nshards=nshards).run()
        from math import comb
        print(f"RESULT j={j} shard={shard}/{nshards} | branches {sw.branches} | "
              f"leaves {sw.leaves} | covering-viol {len(sw.violations)} | "
              f"noncov-below {len(sw.noncov_below)} | maxw {sw.maxdepth_w}",
              flush=True)
        if sw.violations:
            print("  !! THM-726/883 COVERING VIOLATIONS FOUND -- report to fleet",
                  flush=True)
    print(f"done shard {shard}/{nshards}", flush=True)

if __name__ == "__main__":
    main()
