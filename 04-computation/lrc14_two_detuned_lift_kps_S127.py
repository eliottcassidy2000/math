# -*- coding: utf-8 -*-
# kind-pasteur-2026-07-10-S127: the mod-2g lift for the (2,2) half-harmonic residual.
#
# The last detuned residual: v = g*H u {d1,d2} with q1=q2=2 (both d ≡ g/2 mod g, congruent half-harmonics).
# The branch count at g SATURATES (two q=2 coords each cover g/2 branches => fill [0,g)); no branch clears both.
#
# monad's THM-678 mod-2g lift is a 2-adic DESCENT (canon THM-678-multi-detuned-counting-dispatch.md):
# at 2g the two detuned become q=4 (good), BUT every ODD harmonic multiplier m becomes a NEW half-harmonic
# of 2g (g*m ≡ g mod 2g => q=2). So the obstruction descends one level; full closure is at the open core.
#
# This script VERIFIES the split: the lift CLOSES iff there are no odd harmonic multipliers (k=0), and
# DESCENDS (fails at 2g) when k>=1. It is the evidence behind LRCTwoDetunedLift.lonely14_of_two_detuned_lift2g.
from math import gcd
from fractions import Fraction as F

def q_at(s, Q):
    return Q // gcd(s, Q)

def Nq(q):
    return q // 7 + 1                       # bad-count numerator; per-coord density = Nq/q

def generic_fires_at(speeds, Q):
    """Does opus/kps generic counting fire at modulus Q? (sum of N_i/q_i over Q-non-multiples < 1)"""
    det = [s for s in speeds if s % Q != 0]
    S = sum(F(Nq(q_at(s, Q)), q_at(s, Q)) for s in det)
    return S < 1, len(det), S

def is_lonely(speeds):
    """Exact: exists t in [0,1) with all ||s t|| >= 1/14 ?"""
    bps = {F(0), F(1)}
    for s in speeds:
        s = abs(s)
        for k in range(s):
            bps.add(F(14 * k + 1, 14 * s)); bps.add(F(14 * k + 13, 14 * s))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    for a, b in zip(bps, bps[1:]):
        mid = (a + b) / 2
        if all(F(1, 14) <= (s * mid) % 1 <= F(13, 14) for s in speeds):
            return True
    return False

def main():
    print("(2,2) half-harmonic family, mod-2g lift, vary #odd multipliers k:")
    g = 6                                    # g even, g/2 = 3
    for k in range(0, 4):
        odd_mults  = [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21][:k]
        even_mults = [2, 4, 8, 10, 14, 16, 20, 22, 26, 28, 32][:11 - k]
        v = [g * m for m in odd_mults + even_mults] + [(g // 2) * 1, (g // 2) * 5]
        assert len(v) == 13
        qg = [q_at(s, g) for s in [(g // 2) * 1, (g // 2) * 5]]
        fires, ndet, S = generic_fires_at(v, 2 * g)
        print(f"  k={k}: detuned q@g={qg}; @2g detuned-count={ndet}, sum N/q={S}={float(S):.3f}, LIFT FIRES={fires}")
    print("=> k=0 (no odd multiplier): lift fires at 2g (2 detuned at q=4). k>=1: fails (odd mults -> q=2 at 2g).")
    print("=> mod-2g lift CLOSES iff no odd harmonic multiplier; else DESCENDS (monad THM-678, open core).")
    print()
    # k=0 base case is genuinely lonely (LRC(14) true here):
    g = 2
    v = [g * m for m in [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22]] + [1, 3]   # harmonics 4k' (even mult), detuned 1,3 odd
    print(f"k=0 base-case family {sorted(v)}: lonely={is_lonely(v)}  (closed in Lean by lonely14_of_two_detuned_lift2g)")

if __name__ == "__main__":
    main()
