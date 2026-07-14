#!/usr/bin/env python3
"""
lrc14_F_telescoping_xi_thm744_opus_S277.py
==========================================
opus-2026-07-13-S277.  THM-744: the F-TELESCOPING refinement of Xi (THM-742's partial-wrap caps),
using the EXACT segment-end heights.

THE STRUCTURE.  Per maximal exposed segment of a slope -j boundary line, the first-order
sawtooth contribution telescopes: with F(x) = x(1-x)/2 (an antiderivative of psi({x}) = 1/2-{x}
along the height march, CONTINUOUS across wraps since F(0)=F(1)=0),

    (j/W) * Sum_crossings psi(s0)  =  orient * [F(x_start) - F(x_end)]  +  O(j/W)-corrections,

where x_start, x_end are the segment's END HEIGHTS -- exact rationals from the arrangement
(the meeting heights of the pair events that terminate the segment).  Hence three valid
per-segment first-order costs, and Xi may take the least:

    cost_seg = min( j*du/2,  1/4,  |F(x_s)-F(x_e)| (* + 2j/W -> C2 inflation *) ).

DEEPER (the loop structure, logged as the bridge to cross-line cancellation): along each closed
boundary loop of R the segments CHAIN -- at a same-orientation vertex the F(x*) terms of the
two joined segments CANCEL; only run-birth/death vertices (the r = +-1/7 pair events) survive,
weighted +-2F(x*), and r=0 swap vertices are FREE.  The SIGNED total Xi_signed =
|Sum_seg orient*(F(x_s)-F(x_e))| measures how much cross-segment cancellation is available
arrangement-side; computed below as a diagnostic (it is the quantity a joint estimate would bound).

This script: extends the segment walk with exact end heights; computes Xi_old, Xi_new (with the
exact C2 inflation from the DF-selected segments), Xi_signed; re-verifies the full THM-742/743
bound with the new constants over the W-spread (exact, zero-violation check); reports new W0.
"""
from fractions import Fraction as F

LAM = F(1, 14)

def normalize(ivs):
    ivs = sorted((a, b) for a, b in ivs if b > a)
    out = []
    for a, b in ivs:
        if out and a <= out[-1][1]:
            if b > out[-1][1]: out[-1] = (out[-1][0], b)
        else: out.append((a, b))
    return out

def intersect(A, B):
    out = []; i = j = 0
    while i < len(A) and j < len(B):
        a1, b1 = A[i]; a2, b2 = B[j]
        lo, hi = max(a1, a2), min(b1, b2)
        if hi > lo: out.append((lo, hi))
        if b1 <= b2: i += 1
        else: j += 1
    return out

def measure(A): return sum(b - a for a, b in A)

def safe_set(w):
    return [((k + LAM) / w, (k + 1 - LAM) / w) for k in range(w)]

def good_intervals(speeds):
    cur = [(F(0), F(1))]
    for w in sorted(speeds): cur = intersect(cur, safe_set(w))
    return cur

def A_J(u, J):
    ivs = []
    for j in J:
        c = (-j * u) % 1
        a, b = c - LAM, c + LAM
        if a < 0: ivs += [(a + 1, F(1)), (F(0), b)]
        elif b > 1: ivs += [(a, F(1)), (F(0), b - 1)]
        else: ivs.append((a, b))
    return 1 - measure(normalize(ivs))

def breakpoints(J, extra):
    pts = set(extra)
    diffs = sorted({abs(a - b) for a in J for b in J if a != b})
    for d in diffs:
        for k in range(d + 1):
            for off in (F(0), F(1, 7), -F(1, 7)):
                u = (F(k) + off) / d
                if 0 <= u <= 1: pts.add(u)
    pts.add(F(0)); pts.add(F(1))
    return sorted(pts)

def integrate_AJ(J, dom):
    pts = breakpoints(J, [e for iv in dom for e in iv])
    tot = F(0)
    def integ(a, b, fa, fb, d=0):
        m = (a + b) / 2; fm = A_J(m, J)
        assert fm * 2 == fa + fb or d <= 12
        if fm * 2 == fa + fb:
            return (fa + fb) * (b - a) / 2
        return integ(a, m, fa, fm, d + 1) + integ(m, b, fm, fb, d + 1)
    for a, b in zip(pts, pts[1:]):
        if a == b: continue
        for da, db in dom:
            lo, hi = max(a, da), min(b, db)
            if hi > lo: tot += integ(lo, hi, A_J(lo, J), A_J(hi, J))
    return tot, pts

def Ffun(x): return x * (1 - x) / 2

def segments_with_heights(J, GB, pts):
    """maximal exposed segments per sloped line, with exact end u's and heights."""
    pieces = []
    for a, b in zip(pts, pts[1:]):
        if a == b: continue
        for da, db in GB:
            lo, hi = max(a, da), min(b, db)
            if hi > lo: pieces.append((lo, hi))
    pieces.sort()
    segs = []
    for j in J:
        if j == 0: continue
        for sign in (1, -1):
            run = None; prev_hi = None
            for lo, hi in pieces:
                mid = (lo + hi) / 2
                x = (-j * mid + F(sign, 14)) % 1
                exposed = True
                for j2 in J:
                    if j2 == j: continue
                    d = (x - (-j2 * mid)) % 1
                    if (d < LAM or d > 1 - LAM) and d != LAM and d != 1 - LAM:
                        exposed = False; break
                if exposed:
                    if run is not None and prev_hi == lo:
                        run = (run[0], hi)
                    else:
                        if run is not None: segs.append((j, sign, run[0], run[1]))
                        run = (lo, hi)
                    prev_hi = hi
                else:
                    if run is not None:
                        segs.append((j, sign, run[0], run[1])); run = None
            if run is not None: segs.append((j, sign, run[0], run[1]))
    return segs

def xi_variants(segs):
    xi_old = F(0); xi_new = F(0); infl = 0; signed = F(0); n_df = 0
    for j, sign, u1, u2 in segs:
        du = u2 - u1
        x1 = (-j * u1 + F(sign, 14)) % 1
        x2 = (-j * u2 + F(sign, 14)) % 1
        dF = abs(Ffun(x1) - Ffun(x2))
        c_old = min(F(j) * du / 2, F(1, 4))
        c_new = min(c_old, dF)
        xi_old += c_old
        xi_new += c_new
        if c_new == dF and dF < c_old:
            infl += 2 * j; n_df += 1
        orient = -1 if sign == 1 else 1
        signed += orient * (Ffun(x1) - Ffun(x2))
    return xi_old, xi_new, infl, n_df, abs(signed)

CFG = [
    ("shape 1 {1}u{W..W+11}", [1], list(range(12)), F(3690), 339),
    ("shape 2 {1,2,3}u{W..W+9}", [1, 2, 3], list(range(10)), F(1971), 513),
]

print("=" * 104)
print("THM-744: F-telescoping refinement of Xi (exact segment-end heights).  F(x) = x(1-x)/2")
print("=" * 104)

for name, Bset, J, C2_743, oldW0 in CFG:
    GB = good_intervals(Bset)
    area, pts = integrate_AJ(J, GB)
    segs = segments_with_heights(J, GB, pts)
    xi_old, xi_new, infl, n_df, xi_signed = xi_variants(segs)
    C1_old = 2 * (len(GB) + xi_old)
    C1_new = 2 * (len(GB) + xi_new)
    C2_new = C2_743 + infl
    print("\n%s:  Area = %.6f ;  %d segments" % (name, float(area), len(segs)))
    print("   Xi: old = %.4f  ->  new = %.4f   (%d segments use the dF form; C2 inflation +%d)"
          % (float(xi_old), float(xi_new), n_df, infl))
    print("   Xi_SIGNED (loop-telescoped diagnostic) = %.4f   <- the cross-line-cancellation target"
          % float(xi_signed))
    print("   C1: %.2f -> %.2f ;  C2: %d -> %d" % (float(C1_old), float(C1_new), C2_743, C2_new))
    print("   %6s %12s %14s %8s" % ("W", "|L-Area|", "new bound", "ok"))
    viol = 0
    for W in [10, 20, 30, 50, 90, 150, 250, 400, 800]:
        body = list(Bset) + [W + j for j in J]
        if len(set(body)) < 13: continue
        Lx = measure(good_intervals(body))
        err = abs(Lx - area)
        bnd = C1_new / W + C2_new / (W * W)
        ok = err <= bnd
        viol += 0 if ok else 1
        print("   %6d %12.6f %14.6f %8s" % (W, float(err), float(bnd), ok))
    lo, hi = 1, 10 ** 7
    while lo < hi:
        mid = (lo + hi) // 2
        if C1_new / mid + C2_new / F(mid * mid) < area: hi = mid
        else: lo = mid + 1
    print("   NEW W0 = %d  (THM-743: %d)   violations: %d" % (lo, oldW0, viol))

print("\ndone.")
