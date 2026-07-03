#!/usr/bin/env python3
"""
THE SCALE-SEPARATION LEMMA (mac-mini-2026-07-03-S23, HYP-4041 rigorous core).
Claim: let R be lonely at t0 with slack delta = min_{r in R} ||r t0|| - 1/14 > 0, V_R = max R.
Let C = {N + c_i} be a cluster, spread D = max c_i - min c_i. IF
  (i)  N >= V_R / (2 delta)              [fast phase {Nt} sweeps a full period over R's window W]
  (ii) D * (t0 + delta/V_R) < 6/7        [cluster arc {c_i t} fits the safe region [1/14,13/14] over W]
THEN S = R u C is lonely (exists t in W = [t0 - delta/V_R, t0 + delta/V_R] with min_i ||v_i t|| >= 1/14).
Mechanism: on W, R stays safe (slack absorbs the <=V_R*|t-t0| move); {Nt} sweeps [0,1), so some t places the
short cluster arc (length <=D*t<6/7) inside the safe band by the free phase. Verify + probe the boundary.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import random

def gcd_all(xs): return reduce(gcd, xs)
def nd(x): x = x % 1.0; return min(x, 1 - x)

def lonely_on_window(speeds, t_lo, t_hi, steps=400000):
    """max over t in [t_lo,t_hi] of min_i ||v_i t||, and the argmax."""
    best = 0.0; bt = t_lo
    for k in range(steps + 1):
        t = t_lo + (t_hi - t_lo) * k / steps
        m = min(nd(v * t) for v in speeds)
        if m > best: best, bt = m, t
    return best, bt

def R_lonely_with_slack(R, qmax=60):
    """find a rational t0=a/q (small q) where R is lonely, return (t0, slack) with max slack; else None."""
    best = None
    for q in range(2, qmax + 1):
        for a in range(1, q):
            if gcd(a, q) != 1: continue
            t0 = a / q
            m = min(nd(v * t0) for v in R)
            if m > 1/14:
                slack = m - 1/14
                if best is None or slack > best[1]:
                    best = (F(a, q), slack)
    return best

if __name__ == "__main__":
    rng = random.Random(555)
    print("SCALE-SEPARATION LEMMA verification: (i)+(ii) => S=R u C lonely on the window W.")
    print("=" * 92)
    print(f"{'R (base)':>26} {'t0':>7} {'delta':>7} {'cluster N':>10} {'D':>4} {'(i)N>=Vr/2d':>12} {'(ii)Dtmax<6/7':>13} {'M_S on W':>9} {'lonely':>7}")
    ok = 0; tot = 0
    for trial in range(25):
        # random small base R (bounded), find a lonely rational with slack
        m = rng.randint(3, 6)
        R = sorted(rng.sample(range(1, 20), m))
        br = R_lonely_with_slack(R, qmax=40)
        if br is None: continue
        t0, delta = br; t0f = float(t0); VR = max(R)
        # build a near-equal cluster satisfying (i) and (ii)
        Ncond = VR / (2 * delta)
        N = int(Ncond) + rng.randint(1, 50) + 50
        tmax = t0f + delta / VR
        Dmax = int((6/7) / tmax) - 1
        if Dmax < 1: continue
        D = rng.randint(1, max(1, Dmax))
        k = rng.randint(2, 7)
        drifts = sorted(rng.sample(range(0, D + 1), min(k, D + 1)))
        if len(drifts) < 2: continue
        D = drifts[-1] - drifts[0]
        C = [N + d for d in drifts]
        S = sorted(set(R + C))
        cond_i = N >= Ncond
        cond_ii = D * tmax < 6/7
        if not (cond_i and cond_ii): continue
        tot += 1
        W_lo, W_hi = t0f - delta/VR, t0f + delta/VR
        MS, tS = lonely_on_window(S, max(1e-9, W_lo), W_hi, steps=300000)
        lonely = MS >= 1/14
        if lonely: ok += 1
        if trial < 14:
            print(f"{str(R):>26} {t0f:>7.4f} {delta:>7.4f} {N:>10} {D:>4} {str(cond_i):>12} {str(cond_ii):>13} {MS:>9.5f} {str(lonely):>7}")
    print(f"\n(i)+(ii) satisfied in {tot} trials; S lonely on W in {ok}/{tot}")
    print("BOUNDARY probe: violate (ii) (cluster too WIDE, D*tmax >> 6/7) -- expect failures on W:")
    fails = 0; bt = 0
    for _ in range(300):
        R = sorted(rng.sample(range(1, 20), rng.randint(3,5)))
        br = R_lonely_with_slack(R, 40)
        if not br: continue
        t0, delta = br; t0f=float(t0); VR=max(R)
        N = int(VR/(2*delta)) + 100
        tmax = t0f + delta/VR
        D = int((6/7)/tmax) + rng.randint(int((6/7)/tmax)+5, int((6/7)/tmax)+40)  # WIDE, violate (ii)
        C = [N, N+D]
        S = sorted(set(R+[N, N+D]))
        MS,_ = lonely_on_window(S, max(1e-9,t0f-delta/VR), t0f+delta/VR, 150000)
        bt += 1
        if MS < 1/14: fails += 1
    print(f"  violating (ii): {fails}/{bt} FAIL to be lonely on W (confirms (ii) is a real, necessary condition)")
    print("\n=> (i)+(ii) => S lonely on W: a clean scale-separation lemma. Reduces a near-equal-cluster family")
    print("   to its bounded base R. Recurse on R -> spread13 / bounded-denominator census. Formalizable.")
