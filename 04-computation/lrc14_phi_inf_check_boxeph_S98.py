#!/usr/bin/env python3
"""
COMPUTE Phi_inf(frame) AND CHECK sqrt(C) < Phi_inf  (boxeph-2026-07-18-S98)

S97 pinned the density route's resonant peel: Q_s(d) = C*d^2 with C a FIXED
frame constant, so the Cauchy-Schwarz proxy Error ~ sqrt(C).  The row closes iff
the two-scale Error stays below the main term Phi_inf(frame).  Here we compute
the GROUND TRUTH directly from good-set measures (no Fourier reconstruction):

  good set  G(E) = { t in [0,1) : ||v t|| >= theta  for all v in E },  theta = 1/14
  Phi(E)          = |G(E)|                                  (exact, rational)
  Phi_inf(E')     = (1 - 2 theta) * Phi(E') = (6/7) Phi(E') (far elt equidistributes)
  Error(d)        = Phi(E' u {d}) - Phi_inf(E')             (the TRUE two-scale error)

Check: max_d |Error(d)| < Phi_inf   AND   sqrt(C) vs Phi_inf.
sqrt(C) is the LOSSY CS proxy (per-section, no 1/2pi, no g-hat weights); the
direct Error(d) is the honest quantity.

Also (TASK 2 seed): test the fixed-frame scaling law on the REAL deep-well frame
{1..12} with far element 182 -- does Q_s = C12 * d^2 hold, and what is C12?

Pure Python, exact Fractions for the good-set measures.
"""

from fractions import Fraction as Fr
from math import lcm, gcd, pi, sqrt
import cmath


# ---------- exact good-set measure ---------------------------------------
def bad_arcs(v, theta):
    """arcs (lo,hi) in [0,1) where ||v t|| < theta:  t in (k/v - theta/v, k/v + theta/v)."""
    out = []
    for k in range(v):
        c = Fr(k, v)
        out.append((c - theta / v, c + theta / v))
    return out


def union_measure(arcs):
    """measure of a union of (lo,hi) arcs on the circle [0,1) (wrap-aware), exact."""
    # normalise into [0,1), splitting wrap-around arcs
    segs = []
    for lo, hi in arcs:
        lo %= 1
        hi = lo + (hi - lo)          # keep width; hi may exceed 1
        if hi - lo >= 1:
            return Fr(1)
        hi_m = hi % 1
        if hi_m < lo or hi > 1:
            segs.append((lo, Fr(1)))
            segs.append((Fr(0), hi_m))
        else:
            segs.append((lo, hi))
    segs.sort()
    tot = Fr(0)
    cur_lo, cur_hi = None, None
    for lo, hi in segs:
        if cur_hi is None:
            cur_lo, cur_hi = lo, hi
        elif lo <= cur_hi:
            if hi > cur_hi:
                cur_hi = hi
        else:
            tot += cur_hi - cur_lo
            cur_lo, cur_hi = lo, hi
    if cur_hi is not None:
        tot += cur_hi - cur_lo
    return tot


def good_measure(E, theta=Fr(1, 14)):
    arcs = []
    for v in E:
        arcs += bad_arcs(v, theta)
    return 1 - union_measure(arcs)


# ---------- TASK 1: Phi_inf(frame {1..6}) and the Error scan --------------
theta = Fr(1, 14)
frame = [1, 2, 3, 4, 5, 6]
Phi_frame = good_measure(frame, theta)
Phi_inf = (1 - 2 * theta) * Phi_frame          # (6/7) Phi_frame
print("=" * 74)
print("TASK 1 -- Phi_inf(frame {1..6}) and the true two-scale Error")
print("=" * 74)
print(f"theta = 1/14 ;  Phi(frame {{1..6}}) = {Phi_frame} = {float(Phi_frame):.5f}")
print(f"Phi_inf = (6/7)*Phi_frame = {Phi_inf} = {float(Phi_inf):.5f}")
print()
# scan far elements d; find the worst (min Phi, max |Error|)
worst_lo = (None, Fr(1))     # min Phi(d)
worst_err = (None, Fr(0))    # max |Error|
rows = []
for d in range(7, 421):
    Phi_d = good_measure(frame + [d], theta)
    err = Phi_d - Phi_inf
    if Phi_d < worst_lo[1]:
        worst_lo = (d, Phi_d)
    if abs(err) > worst_err[1]:
        worst_err = (d, abs(err))
    rows.append((d, Phi_d, err))
print(f"scanned d = 7..420 :  min_d Phi(d) = {float(worst_lo[1]):.5f} at d={worst_lo[0]}"
      f"  (still > 0 => good set never empty)")
print(f"max_d |Error(d)| = {float(worst_err[1]):.5f} at d={worst_err[0]}")
print()
# a few representative rows incl. resonant (mult of 60=lcm) and generic
for d in [7, 42, 60, 120, 180, 210, 420]:
    Phi_d = good_measure(frame + [d], theta)
    err = Phi_d - Phi_inf
    res = "RESONANT(60|d)" if d % 60 == 0 else ("res(gcd>1)" if gcd(d, 60) > 1 else "generic")
    print(f"  d={d:>4} {res:>15}:  Phi(d)={float(Phi_d):.5f}  Error={float(err):+.5f}"
          f"  |Error|/Phi_inf={float(abs(err)/Phi_inf):.3f}")
print()

# the CS proxy sqrt(C) from S97 (worst section), for the comparison the owner asked
sqrtC = 0.3694
print(f"THE CHECK the owner asked:  sqrt(C) (CS proxy, S97) = {sqrtC:.4f}   vs   "
      f"Phi_inf = {float(Phi_inf):.4f}")
print(f"  sqrt(C) < Phi_inf ?  {sqrtC < float(Phi_inf)}   "
      f"(proxy verdict) ;  BUT sqrt(C) omits 1/2pi and the g-hat weights -->")
print(f"  honest CS bound ~ (1/2pi)*sqrt(C)*||ghat|| = {sqrtC/(2*pi):.4f}*||ghat|| ; "
      f"and the TRUE max|Error| above is the ground truth.")
print(f"  TRUE max|Error| = {float(worst_err[1]):.5f}  <  Phi_inf = {float(Phi_inf):.5f} ?  "
      f"{worst_err[1] < Phi_inf}")


# ---------- TASK 2 seed: scaling law on the REAL deep-well frame {1..12} ---
print()
print("=" * 74)
print("TASK 2 -- fixed-frame scaling law on the deep-well frame {1..12}, far elt d")
print("=" * 74)


def endpoints(E, s):
    bps = sorted(set(Fr(k, 7 * e) for e in E for k in range(7 * e)) | {Fr(0), Fr(1)})
    pts, prev_in, first_in = [], None, None
    for i in range(len(bps) - 1):
        mid = (bps[i] + bps[i + 1]) / 2
        occ = set(int((e * mid % 1) * 7) for e in E)
        cur = (len(occ) == 6) and (s not in occ)
        if prev_in is None:
            first_in = cur
        else:
            if cur and not prev_in:
                pts.append([bps[i], +1])
            elif prev_in and not cur:
                pts.append([bps[i], -1])
        prev_in = cur
    if prev_in != first_in:
        pts.append([Fr(0), +1 if (first_in and not prev_in) else -1])
    return [(p, sg) for p, sg in pts]


def Qs_over_d2(fr, d):
    """Q_s(d)/d^2 at the far element, worst section, for frame `fr` + far d."""
    E = sorted(fr + [d])
    P = 7 * lcm(*E)
    best = 0.0
    for s in range(7):
        pts = endpoints(E, s)
        if not pts:
            continue
        pos = [int(p * P) for p, sg in pts]
        sgn = [sg for p, sg in pts]
        tot = 0.0
        for a in range(len(pts)):
            for b in range(len(pts)):
                x = float((d * (Fr(pos[a] - pos[b], P))) % 1)
                tot += sgn[a] * sgn[b] * x * (1 - x)
        Q = -2 * pi * pi * tot
        best = max(best, Q / d ** 2)
    return best


print("frame {1..12}:  Q_s(d)/d^2 -> C12 (fixed) ?   [note: the deep-well far elt is 182]")
print(f"{'d':>6} {'Qs(d)/d^2':>11}")
for d in [182, 300, 500, 800, 1400]:
    c = Qs_over_d2(list(range(1, 13)), d)
    tag = "  <- deep-well far element" if d == 182 else ""
    print(f"{d:>6} {c:>11.5f}{tag}")
print("If C12 stabilises, the deep-well resonance obeys the SAME fixed-frame scaling")
print("law -- Q_s = C12*d^2 with C12 an invariant of the {1..12} frame's section comb.")
