#!/usr/bin/env python3
"""
lrc14_final17_k8-direct-net-bound_kps-S6-wf.py   (kind-pasteur-2026-06-18-S6)

ANGLE: k8-direct-net-bound.

TARGET (the last gap of LRC(14)):  for every integer co-offset set E with 0 in E,
|E|=k, 8<=k<=12, prove   mu_{1/7}(E) >= thr_k,
where thr_k (EXACT, decreasing) = 1 - min_{|P|=13-k} meas(G_P).

Equivalently (since 1 - mu_{1/7}(E) = meas{x : ALL k gaps of {frac(e_i x)} are <= 1/7}):
   prove   meas{ x in [0,1) : maxgap_E(x) <= 1/7 }  <=  1 - thr_k =: cap_k.

cap_k  =  min_{|P|=13-k} meas(G_P).

This module:
  (1) reproduces the EXACT mu_theta engine (from prompt) and the safe_set / G_P engine,
  (2) recomputes thr_k and cap_k EXACTLY,
  (3) verifies the consecutive values mu_{1/7}(consec_k) and the binding k=8 margins,
  (4) develops the DIRECT net bound: bound meas{all gaps<=1/7} from the
      "every (k-1)-subset's big gap must be split" constraint, WITHOUT assuming
      'consecutive minimizes'.

EXACT Fractions throughout.  stdlib only.
"""
import sys, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)

ONE7 = F(1,7)
HALF14 = F(1,14)

# ----------------------------------------------------------------------------
# EXACT mu_theta engine (verbatim port from prompt; verified consec_8 -> 691/735)
# ----------------------------------------------------------------------------
def mu_theta(E, theta):
    E = sorted(set(E)); n = len(E); bp = set([F(0), F(1)])
    for i in range(n):
        for j in range(i+1, n):
            d = E[j]-E[i]
            for m in range(0, d+1): bp.add(F(m, d))
    bp = sorted(b for b in bp if 0 <= b <= 1); total = F(0)
    for a, b in zip(bp, bp[1:]):
        if b <= a: continue
        mid = (a+b)/2; order = sorted(range(n), key=lambda i: (E[i]*mid) % 1)
        ks = [(E[order[t]]*mid).__floor__() for t in range(n)]; subs = []
        for t in range(n):
            o1 = order[t]; o2 = order[(t+1) % n]; k1 = ks[t]; k2 = ks[(t+1) % n]; wrap = 1 if t == n-1 else 0
            s = E[o2]-E[o1]; c = F(k1-k2+wrap)
            if s == 0:
                if c > theta: subs.append((a, b))
            elif s > 0:
                lo = max(a, (theta-c)/s);  subs.append((lo, b)) if lo < b else None
            else:
                hi = min(b, (theta-c)/s);  subs.append((a, hi)) if a < hi else None
        subs.sort(); cur = cb = None
        for lo, hi in subs:
            if cur is None: cur, cb = lo, hi
            elif lo <= cb: cb = max(cb, hi)
            else: total += cb-cur; cur, cb = lo, hi
        if cur is not None: total += cb-cur
    return total

def mu17(E): return mu_theta(E, ONE7)

# ----------------------------------------------------------------------------
# G_P engine (the small-part safe set; meas(G_P) and thr_k / cap_k)
# ----------------------------------------------------------------------------
def merge(iv):
    iv = sorted(iv); out = []
    for a, b in iv:
        if out and a <= out[-1][1]: out[-1] = (out[-1][0], max(out[-1][1], b))
        else: out.append((a, b))
    return out
def meas(arcs): return sum((b-a for a, b in arcs), F(0))
def complement(arcs):
    arcs = merge(arcs); out = []; prev = F(0)
    for a, b in arcs:
        if a > prev: out.append((prev, a))
        prev = max(prev, b)
    if prev < 1: out.append((prev, F(1)))
    return out
def danger_arcs(u, h=HALF14):
    # {x : ||u x|| < 1/14}: union over residues of (j/u - h/u, j/u + h/u) mod 1
    iv = []
    for j in range(u):
        c = F(j, u); a = (c - h/u) % 1; b = (c + h/u) % 1
        if a < b: iv.append((a, b))
        else: iv.append((a, F(1))); iv.append((F(0), b))
    return iv
def safe_set(P, h=HALF14):
    # G_P = {x : ||p x|| >= 1/14 for all p in P}
    if not P: return [(F(0), F(1))]
    return complement(merge([iv for u in P for iv in danger_arcs(u, h)]))
def measGP(P): return meas(safe_set(P))

def compute_thr_cap():
    thr = {}; cap = {}
    for k in range(8, 14):
        psz = 13 - k
        m = min(measGP(list(P)) for P in itertools.combinations(range(1, 14), psz))
        cap[k] = m
        thr[k] = 1 - m
    return thr, cap

# ----------------------------------------------------------------------------
# MAIN
# ----------------------------------------------------------------------------
if __name__ == "__main__":
    print("="*78)
    print("LRC(14)  k8-direct-net-bound   (kind-pasteur-S6-wf)")
    print("="*78)

    # ---- sanity: engine reproduces known values ----
    print("\n[0] Engine sanity checks (EXACT):")
    for k in (7, 8):
        E = list(range(k))
        v = mu17(E)
        print(f"    mu_1/7(consec_{k}) = {v} = {float(v):.6f}")
    assert mu17(list(range(7))) == F(1), "consec_7 must be 1"
    assert mu17(list(range(8))) == F(691, 735), "consec_8 must be 691/735"
    print("    OK: consec_7 -> 1, consec_8 -> 691/735")

    # ---- thr_k / cap_k EXACT ----
    print("\n[1] thr_k and cap_k (EXACT):")
    thr, cap = compute_thr_cap()
    consec_mu = {}
    for k in range(8, 14):
        E = list(range(k)); cv = mu17(E); consec_mu[k] = cv
        slack = cv - thr[k]
        print(f"    k={k:2d}: cap_k={str(cap[k]):>14s}={float(cap[k]):.5f}  "
              f"thr_k={str(thr[k]):>14s}={float(thr[k]):.5f}  "
              f"mu_consec={str(cv):>14s}={float(cv):.5f}  "
              f"slack={float(slack):+.5f}  {'OK' if slack>=0 else 'FAIL'}")
    # Cross-check against prompt's stated thr/consec values:
    expect_thr8 = F(3637, 5880); expect_consec8 = F(691, 735)
    assert thr[8] == expect_thr8, f"thr_8 mismatch: {thr[8]} vs {expect_thr8}"
    assert consec_mu[8] == expect_consec8
    print(f"    cross-check thr_8={thr[8]} (==3637/5880), 1-thr_8=cap_8={cap[8]}={float(cap[8]):.5f}")
    print(f"    BINDING k=8: need meas(all gaps<=1/7) <= cap_8 = {cap[8]} = {float(cap[8]):.6f}")

    # ------------------------------------------------------------------
    # [2] Direct net-set N(E) = {x : all gaps <= 1/7}, meas = 1 - mu17(E).
    #     We need a DIRECT upper bound on meas(N(E)).  First confirm the
    #     pair-condition: N(E) subset {x : ||(e_i - e_j) x|| in nice ranges}.
    # ------------------------------------------------------------------
    print("\n[2] Net-set net_k(E) := meas(N(E)) = 1 - mu17(E)  (sanity vs cap):")
    for k in range(8, 13):
        E = list(range(k)); netv = 1 - consec_mu[k]
        print(f"    k={k:2d}: meas(N(consec))={str(netv):>14s}={float(netv):.5f}  "
              f"cap_k={float(cap[k]):.5f}  {'<=cap OK' if netv<=cap[k] else 'OVER'}")

    # ------------------------------------------------------------------
    # [3] SINGLE-DIFFERENCE NET BOUND.
    #
    #  Claim (geometric):  if all k gaps of {frac(e_i x)} are <= 1/7, then for
    #  any pair i,j the two points frac(e_i x), frac(e_j x) are separated by at
    #  most  ceil( ||..|| )... more precisely, the ARC between two of the points
    #  is covered by some of the <=1/7 gaps.  For ADJACENT-in-E differences this
    #  gives ||(e_i-e_j) x|| constrained.
    #
    #  Cleanest usable fact:  Let d = e_{(t+1)} - e_{(t)} be the difference of two
    #  E-values.  As x varies, frac(e_i x) - frac(e_j x) = frac(d x) or frac(dx)-1.
    #  But the CIRCLE distance between the two image points is exactly ||d x||.
    #  If x in N(E), the two image points frac(e_i x), frac(e_j x) are both among
    #  the k near-equispaced points, so their circle distance ||d x|| is a SUM of
    #  consecutive gaps -> it is a multiple-ish of ~1/k but each step <= 1/7.
    #
    #  We DON'T yet have a single clean inequality; instead measure, per difference
    #  d, meas{ x in N(E) }  vs  meas{ x : ||d x|| <= 1/7 } and see containment.
    #  Goal: find a difference d (or a small family) s.t. N(E) is contained in a
    #  SMALL union of ||dx|| in [explicit] intervals, giving meas(N(E)) <= small.
    # ------------------------------------------------------------------
    def net_set_intervals(E):
        """Return N(E) = {x: all gaps<=1/7} as a merged interval list (EXACT)."""
        E = sorted(set(E)); n = len(E); bp = set([F(0), F(1)])
        for i in range(n):
            for j in range(i+1, n):
                d = E[j]-E[i]
                for m in range(0, d+1): bp.add(F(m, d))
        bp = sorted(b for b in bp if 0 <= b <= 1); arcs = []
        for a, b in zip(bp, bp[1:]):
            if b <= a: continue
            mid = (a+b)/2
            pts = sorted((E[i]*mid) % 1 for i in range(n))
            gaps = [pts[t+1]-pts[t] for t in range(n-1)] + [pts[0]+1-pts[-1]]
            # On this whole cell the COMBINATORIAL order is fixed; maxgap as a
            # function of x is piecewise linear; but "all gaps<=1/7" is a union of
            # linear conditions.  Reuse mu_theta complement at finer breakpoints:
            # Easiest: net = [0,1] minus the mu17-good region; recompute via theta-cells.
            arcs.append((a, b, gaps))
        return arcs

    # meas(N(E)) directly = 1 - mu17; just expose net intervals for k=8 consec:
    E8 = list(range(8))
    netmeas = 1 - mu17(E8)
    print(f"\n[3] Single-difference containment probe, E=consec_8, meas(N)={netmeas}={float(netmeas):.5f}")
