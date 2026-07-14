#!/usr/bin/env python3
"""
lrc14_rho_residual_identity_thm745_opus_S278.py
===============================================
opus-2026-07-13-S278.  THM-745: the EXACT RESIDUAL IDENTITY and the EUCLIDEAN TOWER.

THE IDENTITY (proved by direct expansion; verified EXACTLY per segment below).  For a height
march h_{m+1} = h_m - alpha (mod 1), alpha = j/W, psi(h) = 1/2 - h, F(h) = h(1-h)/2:

  per NON-WRAP step:  alpha psi(h_m) = -[F(h_{m+1}) - F(h_m)] - alpha^2/2
  per WRAP step (h_m < alpha):  alpha psi(h_m) = -[F(h_{m+1}) - F(h_m)] - alpha^2/2 + (alpha - h_m)

so over any consecutive run of K+1 crossings:

  (j/W) Sum psi(h_m)  =  -[F(h_last+1) - F(h_first)]  -  (K+1) alpha^2/2  +  Sum_wraps (alpha - h_m).

CONSEQUENCES.
 (1) DETERMINISTIC CANCELLATION (microscopic Raabe): #wraps ~ (K+1)alpha, each wrap term ~ alpha/2,
     so Sum_wraps ~ (K+1)alpha^2/2 -- cancelling the -(K+1)alpha^2/2 stream to O(alpha) per segment.
 (2) THE TOWER: the wrap heights h (in [0, alpha)) follow the RETURN MAP of the rotation = rotation
     by (W mod j)/j on the rescaled interval -- the EUCLIDEAN/CF descent of the pair (j, W).  The
     residual of the level-1 sawtooth is a level-2 sawtooth on the CF convergent's clock:
     perspectives all the way down (Mode A / Stern-Brocot / Ostrowski -- the repo's oldest threads).
 (3) THE TERMINAL PER-LINE BOUND: |L - Area| <= [2(#comp(G_B) + Xi_signed)]/W + C2_745/W^2 with
     C2_745 = C2(THM-743) + Sum_seg (2 j_seg + j_seg^2 du_seg / 2 + j_seg^2)  -- C1 collapses to
     ~2.1 (shape 1); the irreducible per-line residual pots are measured and named.

This script: (a) verifies the identity EXACTLY (Fractions) on every segment at W = 90;
(b) computes the signed residual totals Sum orient*rho at several W vs the absolute bound
(measuring the joint-cancellation prize); (c) assembles + verifies the terminal bound over the
W-spread (zero violations) and solves the new W0.
"""
from fractions import Fraction as F

LAM = F(1, 14)
HALF = F(1, 2)

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
        if fm * 2 == fa + fb: return (fa + fb) * (b - a) / 2
        return integ(a, m, fa, fm, d + 1) + integ(m, b, fm, fb, d + 1)
    for a, b in zip(pts, pts[1:]):
        if a == b: continue
        for da, db in dom:
            lo, hi = max(a, da), min(b, db)
            if hi > lo: tot += integ(lo, hi, A_J(lo, J), A_J(hi, J))
    return tot, pts

def Ffun(x): return x * (1 - x) / 2

def segments_with_heights(J, GB, pts):
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
                    if run is not None and prev_hi == lo: run = (run[0], hi)
                    else:
                        if run is not None: segs.append((j, sign, run[0], run[1]))
                        run = (lo, hi)
                    prev_hi = hi
                else:
                    if run is not None: segs.append((j, sign, run[0], run[1])); run = None
            if run is not None: segs.append((j, sign, run[0], run[1]))
    return segs

def rho_decompose(j, sign, u1, u2, W):
    """exact identity pieces over the full-inside strips of a segment; returns None if empty."""
    import math
    m1 = (u1 * W).__ceil__()
    m2 = (u2 * W).__floor__() - 1
    if m2 < m1: return None
    alpha = F(j, W)
    h = (F(sign, 14) - F(j * m1, W)) % 1
    lhs = F(0); wrapsum = F(0); K1 = m2 - m1 + 1
    hfirst = h
    for m in range(m1, m2 + 1):
        lhs += alpha * (HALF - h)
        if h < alpha:
            wrapsum += alpha - h
            h = h - alpha + 1
        else:
            h = h - alpha
    hlastnext = h
    dF = Ffun(hlastnext) - Ffun(hfirst)
    rho = -K1 * alpha * alpha / 2 + wrapsum
    rhs = -dF + rho
    return lhs, rhs, dF, rho, -K1 * alpha * alpha / 2, wrapsum

CFG = [
    ("shape 1 {1}u{W..W+11}", [1], list(range(12)), F(3690), 336),
    ("shape 2 {1,2,3}u{W..W+9}", [1, 2, 3], list(range(10)), F(1971), 462),
]

print("=" * 104)
print("THM-745: exact residual identity + Euclidean tower.  All Fractions.")
print("=" * 104)

for name, Bset, J, C2_743, oldW0 in CFG:
    GB = good_intervals(Bset)
    area, pts = integrate_AJ(J, GB)
    segs = segments_with_heights(J, GB, pts)
    # (a) identity verification at W=90
    W = 90
    nver = 0; nfail = 0
    for (j, sign, u1, u2) in segs:
        r = rho_decompose(j, sign, u1, u2, W)
        if r is None: continue
        lhs, rhs = r[0], r[1]
        if lhs == rhs: nver += 1
        else: nfail += 1
    print("\n%s: identity check at W=90: %d segments verified EXACTLY, %d failures"
          % (name, nver, nfail))
    # (b) signed residual totals
    print("   %6s %16s %14s %14s %16s" % ("W", "Sum orient*rho", "Sum |rho|", "bound", "cancel factor"))
    for W in [90, 150, 250, 400, 800]:
        tot_signed = F(0); tot_abs = F(0)
        for (j, sign, u1, u2) in segs:
            r = rho_decompose(j, sign, u1, u2, W)
            if r is None: continue
            rho = r[3]
            orient = -1 if sign == 1 else 1
            tot_signed += orient * rho
            tot_abs += abs(rho)
        naive = sum(2 * F(j, W) for (j, s, a, b) in segs)
        cf = float(tot_abs / abs(tot_signed)) if tot_signed != 0 else float('inf')
        print("   %6d %16.8f %14.8f %14.8f %16.1f" % (W, float(tot_signed), float(tot_abs), float(naive), cf))
    # (c) terminal bound
    Xi_signed = F(0)
    infl = F(0)
    for (j, sign, u1, u2) in segs:
        du = u2 - u1
        x1 = (-j * u1 + F(sign, 14)) % 1
        x2 = (-j * u2 + F(sign, 14)) % 1
        orient = -1 if sign == 1 else 1
        Xi_signed += orient * (Ffun(x1) - Ffun(x2))
        infl += 2 * j + F(j * j, 1) * du / 2 + j * j  # ends + fluctuation + quadratic margin
    Xi_signed = abs(Xi_signed)
    C1 = 2 * (len(GB) + Xi_signed)
    C2 = C2_743 + infl
    print("   Xi_signed = %.4f ;  C1 = %.3f ;  C2 = %.0f (743: %.0f + inflation %.0f)"
          % (float(Xi_signed), float(C1), float(C2), float(C2_743), float(infl)))
    viol = 0
    for W in [10, 20, 30, 50, 90, 150, 250, 400, 800]:
        body = list(Bset) + [W + j for j in J]
        if len(set(body)) < 13: continue
        Lx = measure(good_intervals(body))
        err = abs(Lx - area)
        bnd = C1 / W + C2 / (W * W)
        if err > bnd: viol += 1
    lo, hi = 1, 10 ** 7
    while lo < hi:
        mid = (lo + hi) // 2
        if C1 / mid + C2 / F(mid * mid) < area: hi = mid
        else: lo = mid + 1
    print("   TERMINAL BOUND: violations over W-spread = %d ;  W0 = %d (was %d)" % (viol, lo, oldW0))

print("\ndone.")
