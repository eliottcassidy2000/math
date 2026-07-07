#!/usr/bin/env python3
"""
klein-2026-07-06-S152 (HYP-4711).  THE r=13 COARSE CRUX: v = a + L*AP is lonely.

The coarse descent M(v) >= M(K) - A/L (kps-S52/S53) fails ONLY when the coarse part K is the
AP {1..13} (M(K)=1/14, no slack). That crux family is  v_i = a_i + L*i  (coarse=AP, fine=a).
[mac-mini's 'escape families' v_i=i+L*k_i with varying k have r<=12 clusters => already closed
 by the coarse descent + LRC(<=13); they were 'escape' only for the finite covering (wrong object).]

CANDIDATE MECHANISM (klein-S152): the dilated AP L*{1..13} has phi(14)=6 optimal witnesses
    t_c = c/(14L),  c in units(14) = {1,3,5,9,11,13},
each binding exactly ONE antipodal pair {i_+, i_-} (i_+ c = +1, i_- c = -1 mod 14).
A small shift t = c/(14L) + delta clears that pair (keeps both >= 1/14) iff
    a_{i_+}/i_+ >= a_{i_-}/i_-.
The conjugate witness c' = 14 - c binds the SAME pair with i_+ <-> i_- swapped, giving the
REVERSE inequality.  So for every pair, one of {c, 14-c} works  =>  M(v) >= 1/14 ALWAYS.

This script VERIFIES the mechanism with EXACT arithmetic (Fraction):
  * construct the analytic delta-interval per c, grid-search it, EXACTLY evaluate min_i ||v_i t||,
  * confirm some c gives min >= 1/14 (a valid witness => M(v) >= 1/14, a rigorous LOWER bound),
  * confirm the analytic slope condition a_{i+}/i+ >= a_{i-}/i- predicts the winning c,
  * stress it: adversarial a, dilated AP (K=d*{1..13}), permuted AP, boundary (all slopes equal).
No Qcap issue: we only certify a LOWER bound via an explicit witness.
"""
from fractions import Fraction as F
import itertools, random

def norm(x):
    """distance to nearest integer, exact."""
    fx = x - (x.numerator // x.denominator)  # in [0,1)
    return min(fx, 1 - fx)

def minreach(v, t):
    return min(norm(v_i * t) for v_i in v)

UNITS14 = [1, 3, 5, 9, 11, 13]
TH = F(1, 14)

def binding_pair(c):
    """i_+ with i_+*c = 1 mod 14, i_- with i_-*c = 13 mod 14 (i in 1..13)."""
    ip = im = None
    for i in range(1, 14):
        if (i * c) % 14 == 1: ip = i
        if (i * c) % 14 == 13: im = i
    return ip, im

def certify_ap_escape(a, L, dil=1, perm=None):
    """
    v_i = a[i] + L*dil*kk[i], kk = perm of {1..13} (default identity = AP).
    Try each c-witness with analytic delta-interval + grid; return (ok, best_min, winning_c, details).
    """
    kk = list(range(1, 14)) if perm is None else perm
    v = [F(a[i] + L * dil * kk[i]) for i in range(13)]
    best = F(0); best_c = None; det = []
    for c in UNITS14:
        # AP witness in s = (dil-scaled) coordinate: for K = dil*{1..13}, min_i ||dil*kk_i * s|| = 1/14
        # at s = c/(14*dil).  Then t = s/L.  Base t0 = c/(14*dil*L).
        t0 = F(c, 14 * dil * L)
        # binding pair is on the AP index: which kk_i hits +-1/14 at this witness
        # ||dil*kk_i*c/(14*dil)|| = ||kk_i*c/14||, binding when kk_i*c = +-1 mod 14
        ip = im = None
        for idx in range(13):
            r = (kk[idx] * c) % 14
            if r == 1: ip = idx
            if r == 13: im = idx
        if ip is None or im is None:
            continue
        # analytic delta interval (linearized): need
        #   a_ip*c/(14 L) + L*dil*kk_ip*delta >= 0   and   a_im*c/(14 L)+L*dil*kk_im*delta <= 0
        denom_p = F(L * dil * kk[ip]); denom_m = F(L * dil * kk[im])
        dmin = -F(a[ip] * c, 14 * L) / denom_p
        dmax = -F(a[im] * c, 14 * L) / denom_m
        slope_ok = (F(a[ip], kk[ip]) >= F(a[im], kk[im]))  # predicted winner
        if dmin <= dmax:
            # grid search inside [dmin, dmax]
            cmax = F(0)
            N = 200
            for j in range(N + 1):
                delta = dmin + (dmax - dmin) * F(j, N)
                mr = minreach(v, t0 + delta)
                if mr > cmax: cmax = mr
            det.append((c, ip, im, slope_ok, True, cmax))
            if cmax > best: best, best_c = cmax, c
        else:
            det.append((c, ip, im, slope_ok, False, None))
    return best >= TH, best, best_c, det

def rough_M(v, Qs):
    """rough upper-ish check: best min over t=p/Q for Q in Qs (to see if family is ~1/14 = sharp)."""
    b = F(0)
    for Q in Qs:
        for p in range(1, Q):
            mr = minreach(v, F(p, Q))
            if mr > b: b = mr
    return b

random.seed(12345)
print("=" * 78)
print("r=13 COARSE CRUX  v = a + L*AP  is LONELY (conjugate-witness delta-adjustment)")
print("=" * 78)

# ---- 1. base case K = {1..13}, random small a, several L ----
print("\n(1) BASE: K={1..13}, v_i = a_i + L*i, random a in [-A,A].")
for L in [183, 500, 3600]:
    for A in [1, 3, 6]:
        fails = 0; sharp = 0; trials = 200
        for _ in range(trials):
            a = [random.randint(-A, A) for _ in range(13)]
            ok, best, bc, det = certify_ap_escape(a, L)
            if not ok: fails += 1
        print(f"   L={L:5d} A={A}:  certified M>=1/14 : {trials-fails}/{trials}"
              + ("   <-- FAIL" if fails else ""))

# ---- 2. adversarial a: try to break it (all a_i equal sign/magnitude, monotone, etc.) ----
print("\n(2) ADVERSARIAL a (structured):")
L = 3600
adv = {
    "all +A":      [6]*13,
    "all -A":      [-6]*13,
    "monotone up": list(range(-6, 7)),
    "monotone dn": list(range(6, -7, -1)),
    "alt +-":      [6 if i%2 else -6 for i in range(13)],
    "one spike":   [0]*6 + [6] + [0]*6,
    "slopes equal (a_i = i)": list(range(1,14)),   # a_i/i all = 1 -> boundary
    "slopes eq neg (a_i=-i)": [-i for i in range(1,14)],
}
for name, a in adv.items():
    ok, best, bc, det = certify_ap_escape(a, L)
    print(f"   {name:24s}: certified={ok}  best_min={float(best):.5f} (1/14={float(TH):.5f}) win_c={bc}")

# ---- 3. verify the SLOPE CONDITION predicts the winning c (a_i+/i+ >= a_i-/i-) ----
print("\n(3) SLOPE-CONDITION check: does a_{i+}/i+ >= a_{i-}/i- predict a valid c? (100 random a)")
L = 3600; mism = 0; nofeasible = 0
for _ in range(100):
    a = [random.randint(-6, 6) for _ in range(13)]
    ok, best, bc, det = certify_ap_escape(a, L)
    # among c with slope_ok True, at least one should be feasible (dmin<=dmax) & clear
    slope_true_feasible = any(d[3] and d[4] and d[5] is not None and d[5] >= TH for d in det)
    if not ok: nofeasible += 1
    if ok and not slope_true_feasible: mism += 1
print(f"   families not certified: {nofeasible}/100 ; certified-but-slope-didn't-predict: {mism}/100")

# ---- 4. DILATED AP  K = d*{1..13} ----
print("\n(4) DILATED AP: K = d*{1..13}, v_i = a_i + L*d*i.")
L = 3600
for d in [1, 2, 3, 5, 7]:
    fails = 0; trials = 100
    for _ in range(trials):
        a = [random.randint(-6, 6) for _ in range(13)]
        ok, best, bc, det = certify_ap_escape(a, L, dil=d)
        if not ok: fails += 1
    print(f"   d={d}:  certified {trials-fails}/{trials}" + ("  <-- FAIL" if fails else ""))

# ---- 5. PERMUTED AP  K = perm of {1..13} ----
print("\n(5) PERMUTED AP: K = random permutation of {1..13}, v_i = a_i + L*K_i.")
L = 3600; fails = 0; trials = 200
for _ in range(trials):
    perm = list(range(1,14)); random.shuffle(perm)
    a = [random.randint(-6, 6) for _ in range(13)]
    ok, best, bc, det = certify_ap_escape(a, L, perm=perm)
    if not ok: fails += 1
print(f"   permuted:  certified {trials-fails}/{trials}" + ("  <-- FAIL" if fails else ""))

# ---- 6. is the family genuinely SHARP (M ~ 1/14)?  sample rough_M for a few ----
print("\n(6) SHARPNESS: rough max-min over small-denominator t (is M ~ 1/14, i.e. genuinely tight?)")
L = 500
for _ in range(5):
    a = [random.randint(-6,6) for _ in range(13)]
    v = [F(a[i] + L*(i+1)) for i in range(13)]
    # denominators near 14L capture AP witnesses; add a few multiples for refinement
    Qs = [14*L, 28*L, 7*L]
    rm = rough_M(v, Qs)
    ok, best, bc, det = certify_ap_escape(a, L)
    print(f"   a={a}  cert_lb={float(best):.5f}  rough_M(smallQ)~{float(rm):.5f}")

print("\nDONE.")
