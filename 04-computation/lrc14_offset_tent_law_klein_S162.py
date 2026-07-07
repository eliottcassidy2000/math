#!/usr/bin/env python3
"""
klein-2026-07-07-S162 -- THM-645: THE OFFSET TENT LAW, its verification, and the 2-anchor
applications (PA_2 exact minimizer-candidates; half-shift quasi-independence census).

THE LAW (derivation in THM-645; generalizes THM-638): for coprime q1,q2 >= 1, theta = c/s
in lowest terms, windows A_i = {x : frac(q_i x) in (alpha_i, alpha_i + theta]} at rational
offsets alpha_i, with r_i = c*q_i mod s and the OFFSET PHASE psi = frac(alpha_2 q1 - alpha_1 q2):
    meas(A_1 ∩ A_2) = theta^2 + (s*Lambda(psi) - r1*r2) / (s^2 q1 q2),
    Lambda(psi) = s * |arc(psi, psi + r1/s] ∩ arc(0, r2/s]|   (mod-1 arc overlap; a tent).
Peak Lambda = min(r1,r2) at psi = 0 (THM-638 same-sign); valley (r1+r2-s)+ (mixed-sign);
piecewise-linear in psi with slopes ±s between.
HALF-SHIFT COROLLARY: alpha_1 = 0, alpha_2 = 1/2 => psi = frac(q1/2): 0 if q1 EVEN
(peak, positively correlated) and 1/2 if q1 ODD. For odd-q directions with r1+r2 <= s and
r's <= s/2 the half-shifted mass is theta^2 - r1r2/(s^2 q1q2) < theta^2: NEGATIVE correlation.
The even/odd duality appears as the tent phase.

THIS SCRIPT:
 1. exhaustive verification of the law (theta = 1/7 and 2/7; offsets in {0, 1/2, 1/3, 1/4, 2/7};
    q <= 40 coprime pairs) against the exact interval engine.
 2. exact-flavored PA_2(AP_k) table, k = 8..13 (grid, two resolutions): the conjectured
    minimizer values of boxeph's discharge object vs the T_k bars; plus P(g0>1/7),
    P(g1/2>1/7), and the overlap (T ∩ sT) for the AP (all-odd? {1..k} is mixed-parity!)
    and for the all-odd AP {1,3,..,2k-1} (where PA_2 = meas(T u T-1/2) exactly).
 3. half-shift quasi-independence census: R_s(E) = meas(T ∩ sT)/(meas T)^2 across shapes
    (T = {gap@0 > 1/7}); skeleton-vs-bulk reading.
"""
import numpy as np
from fractions import Fraction as F
from math import gcd

def intervals_off(q, theta: F, alpha: F):
    out = []
    for j in range(q):
        lo = (j + alpha)/q; hi = lo + theta/q
        out.append((lo, hi))
    return out

def m_exact_offsets(q1, q2, theta: F, a1: F, a2: F):
    A = intervals_off(q1, theta, a1); B = intervals_off(q2, theta, a2)
    tot = F(0)
    for (x0, x1) in A:
        for (y0, y1) in B:
            # wrap both into [0,1) pieces then intersect (windows may cross 1)
            for (u0, u1) in [(x0 % 1, x0 % 1 + (x1-x0))]:
                for (v0, v1) in [(y0 % 1, y0 % 1 + (y1-y0))]:
                    # split at 1
                    pa = [(u0, min(u1, F(1)))] + ([(F(0), u1-1)] if u1 > 1 else [])
                    pb = [(v0, min(v1, F(1)))] + ([(F(0), v1-1)] if v1 > 1 else [])
                    for (p0, p1) in pa:
                        for (s0, s1) in pb:
                            lo, hi = max(p0, s0), min(p1, s1)
                            if hi > lo: tot += hi - lo
    return tot

def law_offsets(q1, q2, theta: F, a1: F, a2: F):
    c, s = theta.numerator, theta.denominator
    r1, r2 = (c*q1) % s, (c*q2) % s
    psi = (a2*q1 - a1*q2) % 1   # convention pinned by exhaustive test: 0 violations
    # Lambda = s * overlap of arc(psi, psi + r1/s] with arc(0, r2/s] on the circle
    lo1, len1 = psi, F(r1, s)
    ov = F(0)
    # circle overlap of (lo1, lo1+len1] with (0, r2/s]
    for shift in (F(0), F(-1)):
        a_lo, a_hi = lo1 + shift, lo1 + shift + len1
        lo, hi = max(a_lo, F(0)), min(a_hi, F(r2, s))
        if hi > lo: ov += hi - lo
    Lam = s * ov
    return theta*theta + (s*Lam - r1*r2) / (s*s*q1*q2)

def grid(NG): return (np.arange(NG)+0.5)/NG

def tail_sets(E, x, t=1/7):
    A = np.asarray(E, float)
    P = (A[None, :]*x[:, None]) % 1.0
    g0 = P.min(axis=1) + (1.0 - P).min(axis=1)
    Q = (P - 0.5) % 1.0
    g1 = Q.min(axis=1) + (1.0 - Q).min(axis=1)
    return g0 > t, g1 > t

if __name__ == "__main__":
    rng = np.random.default_rng(16200)

    print("=== 1. THE OFFSET TENT LAW: exhaustive verification ===")
    for theta in (F(1,7), F(2,7)):
        bad = tested = 0
        for a1, a2 in [(F(0), F(1,2)), (F(0), F(1,3)), (F(1,4), F(1,2)), (F(0), F(2,7)), (F(1,3), F(1,4))]:
            for q1 in range(1, 41):
                for q2 in range(1, 41):
                    if q1 == q2 and q1 > 1: continue
                    if gcd(q1, q2) != 1: continue
                    tested += 1
                    if m_exact_offsets(q1, q2, theta, a1, a2) != law_offsets(q1, q2, theta, a1, a2):
                        bad += 1
        print(f"  theta={theta}: {tested} (pair,offset) cases: {bad} violations")
    print("  half-shift corollary check (psi by parity of q1):")
    for (q1, q2) in [(3, 5), (4, 7), (5, 8), (6, 11)]:
        m = m_exact_offsets(q1, q2, F(1,7), F(0), F(1,2))
        print(f"    (q1,q2)=({q1},{q2}) q1 {'even' if q1%2==0 else 'odd'}: m = {m} = {float(m):.6f}  vs theta^2 = {float(F(1,49)):.6f}"
              f"  ({'>= (peak side)' if m >= F(1,49) else '< (valley side)'})")

    print("\n=== 2. PA_2 minimizer-candidate table (AP_k and all-odd AP), vs T_k ===")
    TK = {8: 0.6185, 9: 0.5057, 10: 0.3956, 11: 0.2747, 12: 0.1429, 13: 0.0565}
    x = grid(120017)
    for k in range(8, 14):
        E = list(range(1, k+1))
        T0, T1 = tail_sets(E, x)
        pa2 = float((T0 | T1).mean())
        Eodd = list(range(1, 2*k, 2))
        S0, S1 = tail_sets(Eodd, x)
        pa2o = float((S0 | S1).mean())
        # all-odd identity check: T1 for odd family = T0 shifted by 1/2
        print(f"  k={k:>2}: PA2(AP_k) = {pa2:.4f}  [P(g0)={T0.mean():.4f}, P(g1/2)={T1.mean():.4f}, overlap={float((T0&T1).mean()):.4f}]"
              f"   PA2(odd-AP) = {pa2o:.4f}   T_k = {TK[k]}   margins {pa2-TK[k]:+.3f}/{pa2o-TK[k]:+.3f}")

    print("\n=== 3. half-shift quasi-independence census: R_s = P(T ∩ sT)/P(T)^2 ===")
    for nm, E in [("AP13", list(range(1,14))), ("odd-AP13", list(range(1,27,2))),
                  ("record", [2,4,6,8,10,11,12,13,14,16,18,20,22]),
                  ("spread", [6,10,14,18,22,26,30,34,38,42,46,50,58]),
                  ("primes", [2,3,5,7,11,13,17,19,23,29,31,37,41])]:
        T0, _ = tail_sets(E, x)
        sT0 = np.roll(T0, len(x)//2)   # T shifted by 1/2 (grid is uniform, len even-ish)
        p = T0.mean(); joint = (T0 & sT0).mean()
        print(f"  {nm:>10}: P(T) = {p:.4f}   P(T ∩ sT) = {joint:.4f}   R_s = {joint/(p*p):.3f}"
              f"   [union = {float((T0|sT0).mean()):.4f} vs 2P-P^2 = {2*p-p*p:.4f}]")
