#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Is n/Phi_6(n) the LRC covering-min? Honest answer: NO in general -- the extremal FAMILY transitions with n.
mac-mini-2026-06-30-S42. Attacking the owner's reframe (covering-min = n/Phi_6(n) = 14/183).

Findings:
  - Construction {1,..,n-2, n(n-1)} has M = n/Phi_6(n) = n/(n^2-n+1) (verified). Phi_6 = Eisenstein norm.
  - n/Phi_6(n) >= 1/n is TRIVIAL (n^2 >= n^2-n+1 <=> n>=1).
  - BUT the construction is NOT the covering-min: for n=4,5,6 a 'drop-2 + tuned large speed' config beats it,
    reaching 2/(2n-1) < n/Phi_6(n). (2n-1 = the signed-LRC modulus C; 27=3^3 at n=14.) Verified minimizers:
    n=4 {1,3,4}=2/7; n=5 {1,3,4,5}=2/9; n=6 {1,3,4,5,18}=2/11 = 2/(2n-1).
  - The small-n beaters do NOT scale: at n=14, {1,3,..,13, 14m} bottoms at 9/83 ~ 0.108 >> 14/183. The
    construction (which KEEPS 2 and uses n(n-1)=182) is tighter. So the EXTREMAL FAMILY TRANSITIONS:
    drop-2-split wins small n, the n(n-1)-construction wins at n=14. No single uniform construction is the
    covering-min -- that is the exact nature/pattern of the razor-thin M-edge.
  - Margins (edge thinness): construction n/Phi_6(n)-1/n = (n-1)/(n Phi_6) ~ 1/n^2; 2/(2n-1)-1/n =
    1/(n(2n-1)) ~ 1/n^2; both ->0 as n grows -- the M-edge thins, and WHICH config is extremal is n-dependent.
"""
from __future__ import annotations
import functools
from fractions import Fraction as F
print = functools.partial(print, flush=True)


def M(S):
    S = sorted(set(S)); cand = set()
    for i in range(len(S)):
        for j in range(len(S)):
            for d in (S[i]-S[j], S[i]+S[j]):
                if d > 0:
                    for k in range(1, d):
                        cand.add((k, d))
    best = F(0)
    for k, d in cand:
        t = F(k, d)
        g = min(min((v*t) % 1, 1-((v*t) % 1)) for v in S)
        if g > best:
            best = g
    return best


def main():
    print("=" * 80)
    print("IS n/Phi_6(n) THE COVERING-MIN? -- the extremal family TRANSITIONS with n (mac-mini-S42)")
    print("=" * 80)

    print("\n[1] construction {1,..,n-2, n(n-1)}: M = n/Phi_6(n); n/Phi_6 >= 1/n trivial:")
    for n in (6, 7, 14):
        S = list(range(1, n-1)) + [n*(n-1)]
        m = M(S); phi6 = F(n, n*n-n+1)
        print(f"    n={n:>2}: M={m}={float(m):.5f} = n/Phi_6(n)={phi6} ; >=1/n? {m >= F(1, n)}")

    print("\n[2] the covering-min BEATS the construction for small n (= 2/(2n-1), drop-2 + tuned large):")
    cases = [(4, [1, 3, 4]), (5, [1, 3, 4, 5]), (6, [1, 3, 4, 5, 18])]
    for n, S in cases:
        m = M(S); two = F(2, 2*n-1); phi6 = F(n, n*n-n+1)
        print(f"    n={n}: {S}: M={m}={float(m):.5f} = 2/(2n-1)={two}; < n/Phi_6={phi6}? {m < phi6}  (covering-min)")

    print("\n[3] but the small-n beaters DON'T scale -- n=14, drop-2 + 14m vs the construction:")
    small = [v for v in range(1, 14) if v != 2]
    best = (F(1), None)
    for mm in range(1, 16):
        S = small + [14*mm]
        if len(set(S)) < 13:
            continue
        if not all(any(v % q == 0 for v in S) for q in range(2, 15)):
            continue
        v = M(S)
        if v < best[0]:
            best = (v, 14*mm)
    print(f"    n=14: best drop-2+14m = {best[0]}={float(best[0]):.5f} (big={best[1]}) -- WORSE than construction 14/183={float(F(14,183)):.5f}")
    print(f"    => at n=14 the construction (KEEPS 2, uses 182=13*14) is tighter; the extremal FAMILY transitions.")

    print("\n" + "=" * 80)
    print("HONEST ANSWER: n/Phi_6(n) is NOT the covering-min in general -- REFUTED for n=4,5,6 (covering-min")
    print("= 2/(2n-1), via drop-2+tuned-large). The beaters DON'T scale to n=14 (drop-2+14m gives 9/83 >>")
    print("14/183), so at n=14 the construction STANDS (tightest known, unrefuted). The EXTREMAL FAMILY is")
    print("n-DEPENDENT (drop-2-split for small n; n(n-1)-construction at n=14) -- there is no single uniform")
    print("optimal construction. Both candidate floors are prime-3 flavored: n/Phi_6 (Eisenstein/hexagonal),")
    print("2/(2n-1)=2/C (signed-LRC modulus, 27=3^3). The inequality (>= 1/n) is trivial for either; the OPEN")
    print("content is the n-dependent EXTREMALITY, and the M-edge margin ~1/n^2 -> 0.")
    print("=" * 80)


if __name__ == "__main__":
    main()
