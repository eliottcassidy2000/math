#!/usr/bin/env python3
"""
klein-2026-07-02-S118 (HYP-4023) - pin the SPREAD PAIR-FLOOR statement + constants.

Target (kps pairCredits first term, pristine window):
  pairmass(w1,w2,a,b) = sum over w2-teeth d of |window cap d cap teeth(w1)|
Claim structure:
  (1) PER-TOOTH IDENTITY: each w2-tooth overlaps AT MOST ONE w1-tooth (near-equal:
      2h/w2 < (1-2h)/w1), and the overlap is the trapezoid
      T(r) = max(0, min(2h/w1, 2h/w2, (h(w1+w2) - |r|)/(w1*w2))),  r = n*w2 - m*w1.
  (2) THE WALK: r_{m+1} = r_m - D (mod w2), D = w2 - w1.
  (3) FLOOR: DL >= 1 => pairmass >= 4h^2 L - errors; pin the exact error constants.

This script verifies (1),(2) exactly (Fractions) and measures the floor quality over
adversarial phases to choose the Lean statement's constants.
"""
from fractions import Fraction as F
import random, math

h = F(1, 14)

def teeth_in(w, a, b):
    """List of (lo,hi) teeth of runner w meeting [a,b] (exact)."""
    import math
    mlo = math.ceil(w * a) - 1
    mhi = math.floor(w * b) + 1
    return [((F(m) - h) / w, (F(m) + h) / w, m) for m in range(mlo, mhi + 1)]

def interval_cap(p, q):
    lo, hi = max(p[0], q[0]), min(p[1], q[1])
    return max(F(0), hi - lo)

def pairmass(w1, w2, a, b):
    """Exact sum over w2-teeth of |[a,b] cap d cap teeth(w1)| ."""
    T1 = teeth_in(w1, a, b)
    tot = F(0)
    per_tooth = []
    for lo2, hi2, m in teeth_in(w2, a, b):
        dlo, dhi = max(lo2, a), min(hi2, b)
        if dlo >= dhi:
            per_tooth.append((m, F(0), 0)); continue
        s = F(0); npart = 0
        for lo1, hi1, n in T1:
            c = interval_cap((max(lo1, dlo), min(hi1, dhi)), (dlo, dhi))
            ov = max(F(0), min(hi1, dhi) - max(lo1, dlo))
            if ov > 0: npart += 1
            s += ov
        per_tooth.append((m, s, npart))
        tot += s
    return tot, per_tooth

def trap(r, w1, w2):
    return max(F(0), min(2*h/w1, 2*h/w2, (h*(w1+w2) - abs(F(r)))/(w1*w2)))

# (1)+(2) verify per-tooth identity + walk on random cases (interior teeth)
random.seed(7)
bad = 0
for trial in range(200):
    w1 = random.randint(50, 400)
    D = random.randint(1, max(1, w1 // 30))
    w2 = w1 + D
    a = F(random.randint(0, 10**6), 10**6)
    L = F(random.randint(50, 500), 100 * w2)  # window a few teeth to many
    b = a + L
    tot, per = pairmass(w1, w2, a, b)
    for m, s, npart in per:
        lo2 = (F(m) - h)/w2
        hi2 = (F(m) + h)/w2
        if lo2 < a or hi2 > b:   # only interior teeth for the identity
            continue
        if npart > 1: bad += 1; print("TWO-PARTNER at", w1, w2, m); break
        # identity: s == T(r) for the nearest n
        n = round(F(m) * w1 / w2)
        rbest = min((abs(n2*w2 - m*w1), n2*w2 - m*w1) for n2 in (n-1, n, n+1))[1]
        if s != trap(rbest, w1, w2):
            bad += 1; print("IDENT FAIL", w1, w2, m, s, trap(rbest, w1, w2)); break
print(f"per-tooth identity + single-partner: {'PASS' if bad==0 else 'FAIL'} over 200 random cases")

# (3) floor quality: adversarial phase sweep. For each (w1,D,L) find the MIN pairmass
# over window positions (adversary chooses a) and compare to 4h^2 L and candidate floors.
print("\nfloor sweep: min over phases of pairmass vs L/49, as DL varies")
print(f"{'w1':>6} {'D':>4} {'DL':>6} {'min/L':>10} {'1/49':>8} {'ratio':>7}")
for (w1, D, Lnum) in [(97, 2, 1), (97, 2, 2), (97, 2, 3), (199, 3, 2), (199, 5, 2),
                      (313, 4, 1), (313, 4, 2), (500, 7, 3), (1000, 11, 2), (150, 1, 1), (150, 1, 2)]:
    w2 = w1 + D
    L = F(Lnum, D)          # DL = Lnum exactly
    best = None
    for k in range(160):    # phase sweep (adversary): shift window by k/(160 D)
        a = F(k, 160 * D) + F(1, 7 * w1 * (k + 3))  # irregular offsets
        tot, _ = pairmass(w1, w2, a, a + L)
        if best is None or tot < best: best = tot
    print(f"{w1:>6} {D:>4} {float(D*L):>6.2f} {float(best/L):>10.6f} {float(F(1,49)):>8.6f} {float(best/L*49):>7.3f}")
