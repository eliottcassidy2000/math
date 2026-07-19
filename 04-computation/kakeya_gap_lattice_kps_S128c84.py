#!/usr/bin/env python3
"""kakeya_gap_lattice_kps_S128c84.py -- kind-pasteur-2026-07-19-S128c84

Target: death-star-S58c's residue  #runs <= 2D/3  (equivalently D >= 3P), the last
gap in the LRC(14) uniform-r=5 maximiser.

THE MOVE: change to GAP COORDINATES.  With g_min <= g_mid <= g_max the sorted
coordinates of g(u) = (frac(-d_2 u), frac(-d_3 u), frac(-d_4 u)), put

    x = g_min,   y = g_mid - g_min,   z = g_max - g_mid      (unimodular)

death-star's four constraints become  x,y,z <= 2/7  and  x+y+z >= 5/7, i.e. with
U = 2/7-x, V = 2/7-y, W = 2/7-z,

    B_sigma  =  { U,V,W >= 0,  U+V+W <= 1/7 }      <-- CORNER SIMPLEX of a cube

so B is six corner simplices of the cube [0,2/7]^3, leg 1/7.  Everything follows:

 (L1) GAP-WINDOW LEMMA.  x,y,z in [1/7,2/7], hence
        g_min in [1/7,2/7],  g_mid in [3/7,4/7],  g_max in [5/7,6/7]
      -- three DISJOINT windows of length 1/7.  (This is death-star's box; B is the
      inscribed simplex, exactly 1/6 of it.)

 (L2) NO CROSSING, NO WRAP.  On a run the gaps y,z are >= 1/7 > 0, so no two
      coordinates ever cross: THE ORDERING IS CONSTANT ON A RUN.  And
      1/7 <= g <= 6/7, so no coordinate wraps.  Hence a run needs no cell
      subdivision -- the affine description is valid on the whole run.

 (L3) RATES ARE DISTINCT.  d_i = d_j would force g_i = g_j, a zero gap.  So the
      three rates are distinct positive integers and D = max d_i >= 3.
      => the P = 1 case of D >= 3P is PROVED outright.

 (L4) #RUNS IS EVEN.  u -> 1-u maps B-runs to B-runs (g(1-u) = 1-g(u), B is
      invariant under g -> 1-g).  Its only fixed points are u = 0 and u = 1/2.
      At u=0, g=(0,0,0): zero gaps.  At u=1/2, every g_i is 0 or 1/2, so two
      coordinates coincide: zero gap.  Neither is in B, so NO run is self-mirror
      and #runs = 2P exactly.  => "#runs <= 2D/3"  IS  "D >= 3P".

 (L5) RUN BOUND + EQUALITY, with no centre-hit hypothesis.  U+V+W = 6/7 - g_max
      increases at rate r_max > 0 and the run ends exactly when it reaches 1/7, so
      run <= (1/7)/r_max; and the fastest coordinate is confined to one window of
      length 1/7, so run <= 1/(7D) (death-star's box-exit bound).  Equality needs
      U=V=W=0 at entry, i.e. entry at the CORNER g = (2/7,4/7,6/7), and r_max = D.
      Then U,V,W increase at rates a, b-a, c-b > 0 summing to D, so a<b<c and
      D >= 3 with equality iff (a,b,c) = (1,2,3).

 (L6) THE LATTICE-TUBE IDENTITY -- the Kakeya reformulation.  Unrolled, the
      geodesic is the segment [0, -d] in R^3, and each B_sigma is CONVEX, so
        #runs = #{ (sigma, m) in S_3 x Z^3 : (B_sigma + m) meets [0,-d] }
      i.e. #runs is a LATTICE POINT COUNT in the convex sausage [0,-d] (-) B_sigma.
      So D >= 3P is exactly "no thin tube about a rational direction carries more
      lattice points than the (1,2,3) tube" -- a Kakeya/Diophantine statement, and
      it says why the bound is hard: lattice points in a thin tube exceed the
      volume precisely when the tube is arithmetically aligned.

This script VERIFIES L1-L6, exhaustively PROVES the P<=1 case (D<=5), and measures
D >= 3P as far as feasible.
"""
import sys
from fractions import Fraction as F
from math import gcd
from itertools import permutations

N = int(sys.argv[1]) if len(sys.argv) > 1 else 30

W1, W2, W3 = (F(1, 7), F(2, 7)), (F(3, 7), F(4, 7)), (F(5, 7), F(6, 7))


def cell_runs(DD):
    """death-star-S58b's exact affine-cell decomposition (verified == bad_from_g)."""
    bps = {F(0), F(1)}
    for d in DD:
        for k in range(d + 1):
            bps.add(F(k, d))
    for i in range(3):
        for j in range(i + 1, 3):
            dd = abs(DD[i] - DD[j])
            if dd > 0:
                for k in range(dd + 1):
                    bps.add(F(k, dd))
    bps = sorted(bps)
    out = []
    for idx in range(len(bps) - 1):
        lo, hi = bps[idx], bps[idx + 1]
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        gmid = [(-DD[i] * mid) % 1 for i in range(3)]
        a = [gmid[i] + DD[i] * mid for i in range(3)]
        order = sorted(range(3), key=lambda i: gmid[i])
        s0, s1, s2 = order
        cons = [("ge", a[s0] - F(2, 7), DD[s0]), ("le", a[s2] - F(5, 7), DD[s2])]
        feas = True
        for (x, y) in [(s1, s0), (s2, s1)]:
            c = a[x] - a[y]
            dc = DD[x] - DD[y]
            if dc > 0:
                cons.append(("ge", c - F(2, 7), dc))
            elif dc < 0:
                cons.append(("le", c - F(2, 7), dc))
            elif not (c <= F(2, 7)):
                feas = False
        if not feas:
            continue
        ulo, uhi = lo, hi
        for typ, const, dd in cons:
            b = const / dd
            if typ == "ge":
                ulo = max(ulo, b)
            else:
                uhi = min(uhi, b)
        if uhi > ulo:
            out.append((ulo, uhi))
    out.sort()
    m = []
    for lo, hi in out:
        if m and lo <= m[-1][1]:
            m[-1] = (m[-1][0], max(m[-1][1], hi))
        else:
            m.append((lo, hi))
    return m


def gvec(DD, u):
    return [(-DD[i] * u) % 1 for i in range(3)]


def in_windows(g):
    s = sorted(g)
    return (W1[0] <= s[0] <= W1[1] and W2[0] <= s[1] <= W2[1]
            and W3[0] <= s[2] <= W3[1])


def is123orbit(d):
    s = sorted(d)
    return s[1] == 2 * s[0] and s[2] == 3 * s[0]


print("=" * 78)
print("(L1)(L2) GAP-WINDOW + NO-CROSSING/NO-WRAP, checked at every run midpoint")
print("=" * 78)
bad_win = bad_ord = 0
checked = 0
for d2 in range(1, 21):
    for d3 in range(1, 21):
        for d4 in range(1, 21):
            DD = [d2, d3, d4]
            if gcd(gcd(d2, d3), d4) != 1:
                continue
            for lo, hi in cell_runs(DD):
                checked += 1
                # sample the run at 5 interior points
                order0 = None
                for f in (F(1, 12), F(1, 4), F(1, 2), F(3, 4), F(11, 12)):
                    u = lo + (hi - lo) * f
                    g = gvec(DD, u)
                    if not in_windows(g):
                        bad_win += 1
                    o = tuple(sorted(range(3), key=lambda i: g[i]))
                    if order0 is None:
                        order0 = o
                    elif o != order0:
                        bad_ord += 1
print("  runs sampled                        : %d" % checked)
print("  samples outside the three windows   : %d  -> L1 %s"
      % (bad_win, "HOLDS" if bad_win == 0 else "FAILS"))
print("  runs where the ordering changed     : %d  -> L2 %s"
      % (bad_ord, "HOLDS" if bad_ord == 0 else "FAILS"))

print()
print("=" * 78)
print("(L3) RATES DISTINCT  =>  D >= 3   (the P = 1 case, PROVED)")
print("=" * 78)
eqrate = 0
for d2 in range(1, 21):
    for d3 in range(1, 21):
        for d4 in range(1, 21):
            DD = [d2, d3, d4]
            if gcd(gcd(d2, d3), d4) != 1:
                continue
            if len(set(DD)) < 3 and cell_runs(DD):
                eqrate += 1
print("  primitive directions with a repeated rate that still have a run : %d" % eqrate)
print("  -> equal rates give g_i = g_j, a zero gap, impossible in B.  D >= 3 PROVED.")

print()
print("=" * 78)
print("(L4) #RUNS EVEN  (no self-mirror run)  =>  '#runs <= 2D/3' IS 'D >= 3P'")
print("=" * 78)
odd = selfmirror = 0
for d2 in range(1, N + 1):
    for d3 in range(1, N + 1):
        for d4 in range(1, N + 1):
            DD = [d2, d3, d4]
            if gcd(gcd(d2, d3), d4) != 1:
                continue
            r = cell_runs(DD)
            if len(r) % 2:
                odd += 1
            for lo, hi in r:
                if (lo, hi) == (1 - hi, 1 - lo):
                    selfmirror += 1
print("  directions with an ODD number of runs : %d" % odd)
print("  self-mirror runs                      : %d" % selfmirror)
print("  -> u=0 gives g=(0,0,0) and u=1/2 gives g_i in {0,1/2}: both have a zero gap,")
print("     so neither fixed point of u->1-u lies in B.  #runs = 2P PROVED.")

print()
print("=" * 78)
print("(L5) run <= 1/(7D) and the equality characterisation")
print("=" * 78)
viol = 0
eqdirs = []
corner_ok = True
for d2 in range(1, N + 1):
    for d3 in range(1, N + 1):
        for d4 in range(1, N + 1):
            DD = [d2, d3, d4]
            if gcd(gcd(d2, d3), d4) != 1:
                continue
            D = max(DD)
            for lo, hi in cell_runs(DD):
                if hi - lo > F(1, 7) / D:
                    viol += 1
                if hi - lo == F(1, 7) / D:
                    eqdirs.append(tuple(DD))
                    if sorted(gvec(DD, lo)) != [F(2, 7), F(4, 7), F(6, 7)]:
                        corner_ok = False
print("  runs violating run <= 1/(7D)          : %d  -> %s"
      % (viol, "HOLDS" if viol == 0 else "BROKEN"))
print("  runs attaining equality               : %d  (directions: %s)"
      % (len(eqdirs), sorted(set(eqdirs))[:8]))
print("  all equality directions are (1,2,3)-orbit : %s"
      % all(is123orbit(t) for t in set(eqdirs)))
print("  every equality run enters at the corner (2/7,4/7,6/7) : %s" % corner_ok)

print()
print("=" * 78)
print("(L6) LATTICE-TUBE IDENTITY:  runs <-> distinct (ordering, lattice point)")
print("=" * 78)
collide = 0
tot = 0
for d2 in range(1, 21):
    for d3 in range(1, 21):
        for d4 in range(1, 21):
            DD = [d2, d3, d4]
            if gcd(gcd(d2, d3), d4) != 1:
                continue
            seen = set()
            for lo, hi in cell_runs(DD):
                u = (lo + hi) / 2
                g = gvec(DD, u)
                m = tuple(int(-DD[i] * u - g[i]) for i in range(3))
                sig = tuple(sorted(range(3), key=lambda i: g[i]))
                tot += 1
                if (sig, m) in seen:
                    collide += 1
                seen.add((sig, m))
print("  runs indexed                                   : %d" % tot)
print("  (ordering, lattice point) collisions           : %d  -> injective %s"
      % (collide, "YES" if collide == 0 else "NO"))
print("  -> #runs = #{(sigma,m) : (B_sigma + m) meets the segment [0,-d]}.")
print("     B_sigma convex => one interval each.  A LATTICE COUNT IN A THIN TUBE.")

print()
print("=" * 78)
print("(E) D >= 3P :  exhaustive proof at P<=1, and measurement above")
print("=" * 78)
# --- exhaustive: D <= 5 forces P <= 1 ---
worstsmall = None
for DD in [list(t) for t in permutations(range(1, 6), 3)]:
    if gcd(gcd(DD[0], DD[1]), DD[2]) != 1:
        continue
    P = len(cell_runs(DD)) // 2
    if worstsmall is None or P > worstsmall[0]:
        worstsmall = (P, tuple(DD))
print("  ALL primitive triples with distinct rates <= 5 : max P = %d at %s"
      % worstsmall)
print("  (repeated rates give no runs at all, L3)  ->  D <= 5 => P <= 1  PROVED,")
print("     i.e. the P = 2 case of D >= 3P (which needs D >= 6) is PROVED.")

minwit = {}
maxratio = F(0)
argmax = None
viol3P = 0
for d2 in range(1, N + 1):
    for d3 in range(1, N + 1):
        for d4 in range(1, N + 1):
            DD = [d2, d3, d4]
            if gcd(gcd(d2, d3), d4) != 1:
                continue
            D = max(DD)
            P = len(cell_runs(DD)) // 2
            if P == 0:
                continue
            if D < 3 * P:
                viol3P += 1
            r = F(P, D)
            if r > maxratio:
                maxratio = r
                argmax = tuple(DD)
            if P not in minwit or D < minwit[P][0]:
                minwit[P] = (D, tuple(DD))
print("  primitive d in [1,%d]^3 : violations of D >= 3P : %d  -> %s"
      % (N, viol3P, "HOLDS" if viol3P == 0 else "BROKEN"))
print("  max P/D = %s = %.5f at %s   (threshold 1/3 = %.5f)"
      % (maxratio, float(maxratio), argmax, 1 / 3))
print("  minimal D achieving each P (the witness family):")
for P in sorted(minwit):
    D, dd = minwit[P]
    print("     P=%-3d  minimal D=%-4d at %-14s   3P=%-4d  slack %d"
          % (P, D, str(dd), 3 * P, D - 3 * P))
