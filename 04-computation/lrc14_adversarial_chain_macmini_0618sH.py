#!/usr/bin/env python3
"""
ANGLE H — ADVERSARIAL re-derivation of the ENTIRE LRC(14) reduction chain.
Independent engine. Exact rational arithmetic (fractions.Fraction).

We re-derive from scratch and cross-check every load-bearing claim:
  (M)   M(S) = max_tau min_{v in S} ||v tau||  (exact, via candidate breakpoints)
  (G)   meas(G_P), G_P = {tau: ||p tau|| >= 1/14  for all p in P}
  (mu)  mu_theta(E) = meas{x: circular maxgap of {frac(e x): e in E} > theta}
  (S7)  meas(S7(E)) = meas{x: each sector [j/7,(j+1)/7) hit by some frac(e x)}
  (M7)  M7(k) closed form
  (link2) GLOBAL WITNESS: gap>1/7 at x  ==>  there is a point safe for all of S? (re-derive)
  (link3) N(E) subset S7(E):  maxgap<=1/7 ==> every sector hit.   exact / boundary.
  (link4) k<=7 pigeonhole: mu_{1/7}=1.

This file = the independent control. Outputs to results/.
"""
from fractions import Fraction as F
from itertools import combinations
import sys

# ----------------------------------------------------------------------
# EXACT M(S): M(S) = max_tau min_{v in S} ||v tau||.
# ||v tau|| is piecewise linear in tau with breakpoints at multiples of
# 1/(2v) (where ||v tau|| has a corner: at k/v it's 0, at (k+1/2)/v it's 1/2).
# The function f(tau)=min_v ||v tau|| is piecewise linear; its max over
# [0,1) is attained at a point where two runners' sawtooths cross, OR at a
# local max of a single sawtooth.  Standard: candidate tau are rationals
# a/q with q = v+w or |v-w| crossings, but the robust exact method:
# M(S) = max over candidate tau in the finite set of "crossing" rationals.
#
# Cleanest exact approach (used in canon): M(S) is rational a/Q where Q
# divides 2*lcm or appears as v_i + v_j etc.  We instead sample the EXACT
# breakpoints of min_v||v tau||.  The maximum of a min of sawtooths occurs
# at tau where ||v tau|| = ||w tau|| for some pair (the binding pair, THM-524),
# i.e.  v tau - a = +-(w tau - b)  =>  tau = (a +- b)/(v +- w).
# Collect all such tau in (0,1), evaluate f, take max.  Also include the
# single-runner peaks tau=(2k+1)/(2v).
# ----------------------------------------------------------------------

def frac_dist(x):
    """||x|| = distance to nearest integer, exact for Fraction x."""
    r = x - int(x)            # in (-1,1)
    if r < 0:
        r += 1                # in [0,1)
    return r if r <= F(1,2) else 1 - r

def M_of_S(S):
    S = sorted(set(S))
    cands = set()
    # single-runner peaks
    for v in S:
        for k in range(0, v):
            cands.add(F(2*k+1, 2*v))
    # pairwise crossings  tau=(a+-b)/(v+-w)
    for v, w in combinations(S, 2):
        for den in (v+w, abs(v-w)):
            if den == 0:
                continue
            for num in range(0, den+1):
                cands.add(F(num, den))
    best = F(0)
    arg = None
    for t in cands:
        if t < 0 or t >= 1:
            continue
        val = min(frac_dist(v*t) for v in S)
        if val > best:
            best = val
            arg = t
    return best, arg

# ----------------------------------------------------------------------
# EXACT meas(G_P):  G_P = { tau in [0,1): ||p tau|| >= 1/14 for all p in P }.
# For a single p, the unsafe set is union over k of (k/p - 1/(14p), k/p + 1/(14p)),
# i.e. ||p tau|| < 1/14.  The safe set for p is a union of p closed arcs each of
# length (1 - 2/14)/p = (6/7)/p... wait: per period 1/p, safe fraction is
# 1 - 2*(1/14) = 6/7 of length 1/p.  We compute the intersection exactly by
# tracking breakpoints.  All breakpoints are rationals of form (k +- 1/14)/p
# = (14k +- 1)/(14 p).  Put everything over common denom and sweep.
# ----------------------------------------------------------------------

def measure_safe(P, half=F(1,14)):
    """meas{tau in [0,1): ||p tau|| >= half for all p in P}.  half=1/14 default."""
    P = sorted(set(P))
    if not P:
        return F(1)
    # breakpoints: tau where some ||p tau|| = half, i.e. p tau = k +- half
    bps = set([F(0), F(1)])
    for p in P:
        # p tau = k +- half  => tau = (k +- half)/p, k=0..p-1 (and edge)
        kmax = p
        for k in range(0, kmax+1):
            for s in (half, -half):
                t = (F(k) + s) / p
                if F(0) <= t <= F(1):
                    bps.add(t)
    bps = sorted(bps)
    total = F(0)
    for i in range(len(bps)-1):
        a, b = bps[i], bps[i+1]
        if b == a:
            continue
        mid = (a+b)/2
        if all(frac_dist(p*mid) >= half for p in P):
            total += (b-a)
    return total

# ----------------------------------------------------------------------
# EXACT mu_theta(E): measure of x where circular maxgap of points
# {frac(e x): e in E} exceeds theta.
# Breakpoints in x: the cyclic ORDER of the points {frac(e x)} changes only
# when two points coincide: frac(e_i x)=frac(e_j x) => (e_i-e_j) x in Z =>
# x = m/(e_i-e_j).  Also a gap can equal theta at x where frac(e_i x) and the
# next point differ by exactly theta or 1-... ; handle by sampling midpoints
# of order-cells AND checking the maxgap there, but maxgap as a function of x
# is continuous & piecewise *linear* within an order cell, so its sup over a
# cell can occur at endpoints.  The set {maxgap>theta} within a cell is an
# interval (maxgap linear in x). We must find where maxgap crosses theta:
#    a particular gap g = frac(e_b x) - frac(e_a x) (mod 1, the gap between
#    cyclically adjacent points) is an affine function (e_b-e_a) x + const...
#    but the "const" jumps by integers; within an order cell each gap is
#    ( (e_b - e_a) x ) mod 1 with a FIXED integer offset => affine in x.
# So crossing maxgap=theta is at x where some gap = theta:  (e_b-e_a)x + c = theta.
# Robust exact method: collect ALL candidate breakpoints
#   x = m/(e_i - e_j)          (order changes)
#   x = (theta + n)/(e_i - e_j)  and  ((1-theta)+n)/(e_i-e_j) ...
# Simpler & fully rigorous: the maxgap>theta indicator is piecewise constant
# between consecutive elements of the breakpoint set
#   B = { (n + s)/d : d = e_i - e_j (i!=j), n in Z, s in {0, theta, -theta} } cap (0,1)
# Evaluate indicator at midpoints. This is the order-cell + gap=theta refinement
# used by the canon engines.  We implement it exactly.
# ----------------------------------------------------------------------

def circular_maxgap(pts):
    """pts: list of Fraction in [0,1). Return max circular gap (Fraction)."""
    xs = sorted(set(pts))
    if len(xs) == 1:
        return F(1)
    g = F(0)
    for i in range(len(xs)):
        nxt = xs[(i+1) % len(xs)]
        gap = nxt - xs[i] if i+1 < len(xs) else (xs[0] + 1 - xs[i])
        if gap > g:
            g = gap
    return g

def mu_theta(E, theta):
    """meas{x in (0,1): circular maxgap{frac(e x): e in E} > theta}, EXACT."""
    E = sorted(set(E))
    diffs = set()
    for a, b in combinations(E, 2):
        d = abs(a-b)
        if d:
            diffs.add(d)
    # also differences with 0-offset point if 0 in E covered above
    bps = set([F(0), F(1)])
    for d in diffs:
        # order-change points  n/d
        for n in range(0, d+1):
            bps.add(F(n, d))
        # gap = theta crossings:  a gap is affine with slope d' = some diff;
        # to be safe add (n + theta)/d and (n - theta)/d for each diff d
        # (covers gap value hitting theta). Use theta as Fraction.
        nlo = -2
        nhi = d + 2
        for n in range(nlo, nhi+1):
            for s in (theta, -theta, 1-theta, theta-1):
                t = F(n) + s
                t = t / d
                if F(0) < t < F(1):
                    bps.add(t)
    bps = sorted(b for b in bps if F(0) <= b <= F(1))
    total = F(0)
    for i in range(len(bps)-1):
        a, b = bps[i], bps[i+1]
        if b == a:
            continue
        mid = (a+b)/2
        pts = [frac_dist_signed(e*mid) for e in E]
        if circular_maxgap(pts) > theta:
            total += (b-a)
    return total

def frac_dist_signed(x):
    """frac(x) in [0,1) for Fraction x (NOT distance-to-int; the actual fractional part)."""
    r = x - int(x)
    if r < 0:
        r += 1
    return r

# ----------------------------------------------------------------------
# EXACT meas(S7(E)): each of 7 sectors [j/7,(j+1)/7) hit by some frac(e x).
# Breakpoints in x: frac(e x) crosses a sector boundary j/7, i.e.
#   e x = m + j/7  => x = (m + j/7)/e = (7m + j)/(7 e).
# Indicator piecewise constant between these. Evaluate at midpoints.
# ----------------------------------------------------------------------

def meas_S7(E):
    E = sorted(set(E))
    bps = set([F(0), F(1)])
    for e in E:
        if e == 0:
            continue
        # e x = m + j/7
        mmax = e
        for m in range(0, mmax+1):
            for j in range(0, 7):
                t = F(7*m + j, 7*e)
                if F(0) <= t <= F(1):
                    bps.add(t)
    bps = sorted(b for b in bps if F(0) <= b <= F(1))
    total = F(0)
    for i in range(len(bps)-1):
        a, b = bps[i], bps[i+1]
        if b == a:
            continue
        mid = (a+b)/2
        hit = set()
        for e in E:
            fx = frac_dist_signed(e*mid)
            sec = int(fx * 7)  # which sector [j/7,(j+1)/7)
            if sec == 7:
                sec = 6
            hit.add(sec)
        if len(hit) == 7:
            total += (b-a)
    return total

def M7(k):
    """closed-form independent main term."""
    from math import comb
    s = F(0)
    for t in range(0, 7):
        s += F((-1)**t) * comb(6, t) * F(7-t, 7)**(k-1)
    return s

# ----------------------------------------------------------------------
if __name__ == "__main__":
    print("="*70)
    print("ANGLE H — independent re-derivation. Anchor cross-checks.")
    print("="*70)

    # --- ANCHOR 1: mu_{2/7} consecutive small k (canon: mu_4=19/21, mu_6=4/7? 9/14?)
    print("\n[mu_{2/7} consecutive]  canon: mu_4=19/21, mu_5=9/14, mu_6=4/7")
    for k in [3,4,5,6,7]:
        E = list(range(k))
        v = mu_theta(E, F(2,7))
        print(f"  k={k}: mu_2/7 = {v} = {float(v):.5f}")

    # --- ANCHOR 2: mu_{1/7} consecutive (canon: mu_7=1, mu_8=691/735, mu_13=477/1078)
    print("\n[mu_{1/7} consecutive]  canon: mu_7=1, mu_8=691/735, mu_9=247/294,")
    print("                              mu_10=38/49, mu_11=1381/2205, mu_13=477/1078")
    for k in [7,8,9,10,11,12,13]:
        E = list(range(k))
        v = mu_theta(E, F(1,7))
        print(f"  k={k}: mu_1/7 = {v} = {float(v):.5f}")

    # --- ANCHOR 3: cap_k = min_{|P|=13-k} meas(G_P)
    # canon cap_8=2243/5880, cap_13=1, m_P=14249/252252
    print("\n[meas(G_P) min over |P|=psz]  canon row psz=1..10:")
    print("   6/7,66/91,55/91,1979/4004,2243/5880,3029/10780,45107/229320,2479/17640,10601/114660,14249/252252")
    for psz in range(1, 11):
        best = None; argP=None
        for P in combinations(range(1,14), psz):
            m = measure_safe(list(P))
            if best is None or m < best:
                best = m; argP = P
        print(f"  |P|={psz}: min meas(G_P) = {best} = {float(best):.6f}   at {argP}")
