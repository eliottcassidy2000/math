#!/usr/bin/env python3
"""Cross-lane test: death-star's THM-3002 Bernstein capacity criterion applied to
the two readings of the eq(27) weight alpha = 2457/6592.

THM-3002 (death-star): with d_m = floor(gamma*m) + D0, the epoch of length R has
d_i = d_{R+i} for i = 0..R-1, and the necessary criterion is

    sum_{i<=t} C(d_i, t-i) * 2^(t-i)  >=  C(R-1, t)      for every t.

THE TWO READINGS of alpha (klein-S428, 07-reflections/eq27-is-a-logit-gate-...):
  R1  gamma = alpha                    = 2457/6592 = 0.3727245...
  R2  alpha = (C-1)/C = gamma/(1+gamma)  =>  gamma = alpha/(1-alpha)
                                       = 2457/4135 = 0.5941959...,  C = 6592/4135

RESULT
  * POSITIVE CONTROL: this implementation reproduces death-star's stated
    trichotomy at gamma = 1/2 -- ample at R = 8, 16 and deficient by R = 64 --
    which is what licenses the readings below.
  * R1 IS REFUTED.  gamma = 2457/6592 is DEFICIENT at every R >= 16 and every
    D0 in {0,1,2}, exponentially so (min ratio 0.008 at R=32, 0.000 at R=64).
  * R2 SURVIVES.  gamma = 2457/4135 is UNIFORMLY AMPLE through R = 1024, with
    the binding index t = 1 and min ratio ~1.19, stable in R -- exactly the
    signature death-star reports for gamma >= 3/5.
  * The criterion does NOT separate 2457/4135 from 3/5 (both bind at t=1 with
    ratio -> 2*gamma).  So it eliminates one reading; it cannot pin the constant.
  * CORRECTION (2026-07-31, after death-star pushed to R = 2048).  An earlier
    version of this file warned "do not extrapolate the finite-R threshold
    sequence", on the grounds that at R = 1024 the binding index for 2457/4135 is
    already t = 1.  THAT WARNING WAS WRONG, and the extrapolation it dismissed was
    right.  The finite-R thresholds 0.5313, 0.5606, 0.5758, 0.5849 (R = 32..256)
    continue 0.59065 (R=512), 0.59393 (R=1024) and converge to ~0.5982, matching
    the asymptotic two-ray entropy value 0.59799.  Since 2457/4135 = 0.594196 lies
    BELOW that limit, it is eventually deficient, and death-star reports it dying
    at R = 2048.  Ampleness at R <= 1024 with binding index t = 1 is therefore NOT
    the ample signature I took it for -- t = 1 binds first and the large-t
    constraint overtakes it only later.
  * NET EFFECT ON THE TWO READINGS: R1 (gamma = 2457/6592) is refuted here at
    R >= 16; R2 (gamma = 2457/4135) is refuted by death-star at R = 2048.  BOTH
    readings of the eq(27) weight are dead as deadline slopes.  gamma = 3/5 is the
    first round rate above the threshold and death-star has now closed every
    dyadic epoch through R = 64 there, giving C <= 8/5 for all n <= 127.

At t = 1 the criterion reads 2*(floor(gamma*R)+D0) >= R-1, i.e. ratio -> 2*gamma,
so t=1 is satisfied for every gamma > 1/2; ampleness is decided by large t.

Reproduce: python3 04-computation/amm12592_capacity_criterion_eliminates_the_raw_weight_reading_klein.py
"""

from fractions import Fraction as Fr
from math import comb, floor


def min_ratio(gamma, R, D0=0):
    """min over t of  [sum_{i<=t} C(d_i,t-i) 2^(t-i)] / C(R-1,t), and the argmin."""
    d = [floor(gamma * (R + i)) + D0 for i in range(R)]
    best, at = None, None
    for t in range(1, R):
        S = 0
        for i in range(0, min(t, R - 1) + 1):
            if t - i <= d[i]:
                S += comb(d[i], t - i) * (1 << (t - i))
        r = Fr(S, comb(R - 1, t))
        if best is None or r < best:
            best, at = r, t
    return best, at


def rule(s):
    print("=" * 74)
    print(s)
    print("=" * 74)


CANDS = [("1/2 (death-star control)", Fr(1, 2)),
         ("2457/6592  = READING R1", Fr(2457, 6592)),
         ("2457/4135  = READING R2", Fr(2457, 4135)),
         ("3/5", Fr(3, 5)),
         ("2/3", Fr(2, 3))]


def main():
    rule("A. CRITERION LEDGER (min ratio; >=1 is AMPLE)")
    for D0 in (0, 1, 2):
        print(f"  --- D0 = {D0} ---")
        print(f"    {'gamma':28s}" + "".join(f"R={R:<13d}" for R in (8, 16, 32, 64)))
        for name, g in CANDS:
            row = ""
            for R in (8, 16, 32, 64):
                r, _ = min_ratio(float(g), R, D0)
                row += f"{float(r):6.3f} {'AMPLE' if r >= 1 else 'DEFIC'}  "
            print(f"    {name:28s}{row}")
        print()

    control_ok = True
    for R, want in ((8, True), (16, True), (64, False)):
        r, _ = min_ratio(0.5, R, 0)
        control_ok &= ((r >= 1) == want)
    print(f"  POSITIVE CONTROL (gamma=1/2 ample at R=8,16 and dead by R=64): {control_ok}")

    r1_dead = all(min_ratio(2457 / 6592, R, D0)[0] < 1
                  for R in (16, 32, 64) for D0 in (0, 1, 2))
    print(f"  READING R1 (gamma=2457/6592) deficient at every R>=16, D0<=2: {r1_dead}")

    print()
    rule("B. DEEP CHECK: does READING R2 survive large R?")
    r2_ok = True
    for R in (256, 512, 1024):
        for name, g in (("2457/4135 (R2)", 2457 / 4135), ("3/5", 0.6)):
            r, t = min_ratio(g, R)
            r2_ok &= (r >= 1)
            print(f"    R={R:5d}  {name:16s} min ratio={float(r):8.4f} at t={t:4d}"
                  f"  {'AMPLE' if r >= 1 else 'DEFICIENT'}   (2*gamma = {2 * g:.4f})")
    print(f"  READING R2 uniformly ample through R=1024, binding at t=1: {r2_ok}")

    print()
    rule("C. FINITE-R THRESHOLDS (do NOT extrapolate these)")
    for R in (32, 64, 128, 256):
        lo, hi = 0.50, 0.65
        for _ in range(30):
            mid = (lo + hi) / 2
            if min_ratio(mid, R)[0] >= 1:
                hi = mid
            else:
                lo = mid
        print(f"    R={R:4d}: finite-R threshold gamma ~ {hi:.6f}")
    print("    These converge to ~0.5982 (death-star: 0.59065 at R=512, 0.59393 at")
    print("    R=1024), matching the asymptotic entropy value 0.59799.  Since")
    print("    2457/4135 = 0.594196 is BELOW that, it is eventually deficient and")
    print("    dies at R=2048.  An earlier version of this file wrongly warned")
    print("    against this extrapolation; the extrapolation was right.")

    print()
    rule(f"SUMMARY  control={control_ok}  R1_refuted={r1_dead}  R2_survives={r2_ok}")


if __name__ == "__main__":
    main()
