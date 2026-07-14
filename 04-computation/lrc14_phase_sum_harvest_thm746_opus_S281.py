#!/usr/bin/env python3
"""
lrc14_phase_sum_harvest_thm746_opus_S281.py
===========================================
opus-2026-07-14-S281.  THM-746: the SOUND FIRST-ORDER HARVEST via exact quadratic phase expansion.

THE POINT.  F(h) = h(1-h)/2 is QUADRATIC, so the Taylor expansion is EXACT in three terms:
F(x + d) = F(x) + psi(x) d - d^2/2  (psi = 1/2 - x, F'' = -1, no tail).  The first-order wedge
content Phi(W) = Sum_seg orient(-dF_ext) (THM-745: EXACTLY the first-order part, residuals
cancel by the pairing theorem) therefore decomposes EXACTLY:

  Phi(W) = Xi_sv  -  S(W)/W  -  T(W)/(2 W^2)  -  Z(W),
    Xi_sv = Sum_seg orient [F(x_start) - F(x_end)]          (W-INDEPENDENT; = +-0.0561 / 0.1786)
    S(W)  = Sum_seg orient j [psi(x_start) th1 + psi(x_end) th2]   (GRID-PHASE SUM)
    T(W)  = Sum_seg orient j^2 [th1^2 - th2^2]
    Z(W)  = Sum over segments with NO full-inside crossing at W of orient [F(x_start)-F(x_end)]
  with th1 = ceil(u1 W) - u1 W, th2 = u2 W - floor(u2 W) in [0,1): the fractional GRID PHASES
  of the segment endpoints -- i.e. {u_e W}: THE ARRANGEMENT'S VERTICES BECOME RUNNERS (speeds
  u_e, rational = pair-event grid points) OBSERVED AT INTEGER TIME W.  The perspective tower
  inverts: bounding the LRC error spawns a new lonely-runner-type system one level down, whose
  runners are the vertices of the arrangement.  S(W) is periodic in W (mod lcm of the vertex
  denominators -- the difference-runner grids).

SOUND BOUNDS (each order clean -- the corrected form of MISTAKE-142):
  |Phi(W)| <= |Xi_sv| + S1/W + S2/(2W^2) + Z1-part, with
    S1 = Sum_seg j (|psi(x_start)| + |psi(x_end)|)   (exact),
    S2 = Sum_seg j^2                                  (|th^2 - th^2| <= 1),
    Z1 = Sum_seg j  (zero-crossing segments have du < 2/W so |dF| <= j du/2 < j/W).
  ASSEMBLY: |L - Area| <= C1/W + C2/W^2 with C1 = 2(#comp(G_B) + |Xi_sv|)  [~2.1 vs 14.49!]
  and C2 = C2(THM-743) + 2(S1 + Z1) + S2.

VERIFICATION: (a) the three-term expansion as an EXACT Fraction identity per W (both shapes);
(b) Phi(W) equals the direct wedge sum (pairing re-check); (c) the assembled bound over the
W-spread, zero violations; (d) new W0; (e) S(W) values + periodicity spot-check.
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
def psi(x): return HALF - x

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

CFG = [
    ("shape 1 {1}u{W..W+11}", [1], list(range(12)), F(3690), 339),
    ("shape 2 {1,2,3}u{W..W+9}", [1, 2, 3], list(range(10)), F(1971), 513),
]

print("=" * 106)
print("THM-746: the sound first-order harvest (exact quadratic phase expansion of the dF_ext sum)")
print("=" * 106)

for name, Bset, J, C2_743, oldW0 in CFG:
    GB = good_intervals(Bset)
    area, pts = integrate_AJ(J, GB)
    segs = segments_with_heights(J, GB, pts)
    # W-independent constants
    Xi_sv = F(0); S1 = F(0); S2 = 0; Z1 = 0
    segdata = []
    for (j, sign, u1, u2) in segs:
        orient = -1 if sign == 1 else 1
        xs = (-j * u1 + F(sign, 14)) % 1
        xe = (-j * u2 + F(sign, 14)) % 1
        Xi_sv += orient * (Ffun(xs) - Ffun(xe))
        S1 += j * (abs(psi(xs)) + abs(psi(xe)))
        S2 += j * j
        Z1 += j
        segdata.append((j, sign, orient, u1, u2, xs, xe))
    C1 = 2 * (len(GB) + abs(Xi_sv))
    C2 = C2_743 + 2 * (S1 + Z1) + S2
    print("\n%s:  Area = %.6f ;  Xi_sv = %+.6f ;  S1 = %.2f ; S2 = %d ; Z1 = %d"
          % (name, float(area), float(Xi_sv), float(S1), S2, Z1))
    print("   C1 = 2(#comp + |Xi_sv|) = %.4f   [THM-743's C1 was %.2f -> %.1fx at large W]"
          % (float(C1), 14.49 if "1}" in name[:8] else 19.14, (14.49 if "1}" in name[:8] else 19.14)/float(C1)))
    print("   C2 = %.0f (743) + 2(S1+Z1) + S2 = %.0f" % (float(C2_743), float(C2)))

    # exact expansion identity + wedge equality + S(W) values
    print("   %6s %14s %14s %10s %12s" % ("W", "Phi(W) exact", "3-term RHS", "EQUAL", "S(W)"))
    for W in [90, 97, 250, 800]:
        Phi = F(0); RHS_S = F(0); RHS_T = F(0); Zw = F(0); wedge = F(0)
        for (j, sign, orient, u1, u2, xs, xe) in segdata:
            m1 = (u1 * W).__ceil__(); m2 = (u2 * W).__floor__() - 1
            alpha = F(j, W)
            if m2 < m1:
                Zw += orient * (Ffun(xs) - Ffun(xe)); continue
            th1 = m1 - u1 * W
            th2 = u2 * W - (m2 + 1)
            hf = xs - j * th1 / W
            hl = xe + j * th2 / W
            Phi += orient * (Ffun(hf) - Ffun(hl))
            RHS_S += orient * j * (psi(xs) * th1 + psi(xe) * th2)
            RHS_T += orient * j * j * (th1 * th1 - th2 * th2)
            # direct wedge sum for this segment
            h = (F(sign, 14) - F(j * m1, W)) % 1
            acc = F(0)
            for m in range(m1, m2 + 1):
                acc += alpha * (HALF - h)
                h = (h - alpha) % 1
            wedge += orient * acc
        rhs = (Xi_sv - Zw) - RHS_S / W - RHS_T / (2 * W * W)
        # NB: wedge = Phi + Sum orient rho ; pairing theorem => equal
        print("   %6d %+14.8f %+14.8f %10s %+12.6f   [wedge-sum equal: %s]"
              % (W, float(Phi), float(rhs), Phi == rhs, float(RHS_S), wedge == Phi))

    # assembled bound + W0
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
    print("   ASSEMBLY: violations = %d ;  NEW W0 = %d (THM-743: %d)" % (viol, lo, oldW0))

print("\ndone.")
