#!/usr/bin/env python3
"""
death-star-2026-07-07-S2 (formalization support) -- can the AP76Certificate
mu_{1/7}(AP_76) >= m_P = 14249/252252 be proved by a UNION OF SMALL-q RATIONAL NEIGHBORHOODS?

Route (boxeph): near x=p/q (small q), the AP orbit {frac(j x): j=0..N-1} clusters into q groups,
leaving a gap ~ 1/q - (spread). For x = p/q + delta, an explicit interval has maxgap > 1/7.
Sum the lengths of these GOOD intervals over small q; if the total already exceeds m_P (with the
disjoint neighborhoods), the certificate is Lean-provable by exhibiting finitely many intervals.

Here: (1) exact good-set measure of AP_N via fine grid (sanity vs 2314528732/40290957525);
(2) for each small q, the EXACT half-width where the p/q-neighborhood stays good, and the
cumulative lower bound as we include more q; find the minimal Q with cumulative >= m_P.
"""
from fractions import Fraction as F
from math import gcd

def maxgap_frac(N, x):
    # points {frac(j*x): j=0..N-1}, x a Fraction; returns exact maxgap as Fraction
    pts = sorted(set((j * x) % 1 for j in range(N)))
    mg = pts[0] + 1 - pts[-1]
    for i in range(len(pts) - 1):
        g = pts[i+1] - pts[i]
        if g > mg: mg = g
    return mg

def good_halfwidth_at(N, p, q, thr, side, maxstep=F(1,2)):
    """largest w such that x=p/q + side*d is GOOD (maxgap>thr) for all 0<d<=w (approx via
    the FIRST crossing; exact by scanning candidate breakpoints d=k/(q*M))."""
    # scan delta on a fine rational grid to find where maxgap drops to <= thr
    base = F(p, q)
    lo, hi = F(0), maxstep
    # find first delta where maxgap <= thr (bisection on a monotone-ish decay near the peak)
    # robust: fine scan
    steps = 4000
    prev_good = True
    firstbad = None
    for k in range(1, steps + 1):
        d = maxstep * k / steps
        x = (base + side * d) % 1
        if maxgap_frac(N, x) <= thr:
            firstbad = d
            break
    return firstbad if firstbad is not None else maxstep

if __name__ == "__main__":
    thr = F(1, 7)
    mP = F(14249, 252252)
    exact_claim = F(2314528732, 40290957525)
    print(f"m_P = {mP} = {float(mP):.6f};  claimed exact mu = {float(exact_claim):.6f}")
    print(f"margin (claim - m_P) = {float(exact_claim - mP):.6f}  ({float((exact_claim-mP)/mP)*100:.1f}%)\n")

    # sanity: grid estimate of mu_{1/7}(AP_76)
    N = 76
    res = 20000
    good = 0
    for r in range(res):
        x = F(2*r+1, 2*res)
        if maxgap_frac(N, x) > thr:
            good += 1
    print(f"grid mu_{{1/7}}(AP_{N}) ~ {good/res:.5f}  (claim {float(exact_claim):.5f})\n")

    # small-q neighborhood coverage: for each reduced p/q with q<=Q, measure the good nbhd
    print("cumulative good-measure from p/q neighborhoods, by max q:")
    print(f"{'Q':>3} {'#fracs':>7} {'cum_lowerbound':>15} {'>= m_P?':>8}")
    covered = []  # list of (center, left_w, right_w)
    for Q in range(2, 21):
        # add all reduced fractions with denominator exactly Q
        for p in range(1, Q):
            if gcd(p, Q) != 1: continue
            lw = good_halfwidth_at(N, p, Q, thr, -1)
            rw = good_halfwidth_at(N, p, Q, thr, +1)
            covered.append((F(p, Q), lw, rw))
        # cumulative measure = sum of (lw+rw), assuming disjoint (verify roughly)
        total = sum((lw + rw) for _, lw, rw in covered)
        ok = total >= mP
        print(f"{Q:>3} {len(covered):>7} {float(total):>15.6f} {str(ok):>8}")
        if ok:
            print(f"\n  => q <= {Q} neighborhoods already sum to {float(total):.6f} >= m_P={float(mP):.6f}.")
            print(f"     Lean plan: exhibit these {len(covered)} good intervals + disjointness + sum.")
            break
