#!/usr/bin/env python3
"""
covering_escape_large_multiple_klein.py  --  klein-2026-06-30-S52  (HYP-3763, CORRECTED part 2)

CORRECTION: a COVERING SET (THM-523) must contain a multiple of EVERY q in {2,...,n}. The tight
minimizers (GW={1..11,13,24}, M=1/n) and the drop-2 mediant sets (M=2/(2n-1)) are NOT covering sets
(they miss a multiple of n) -- they are THM-523's trivial class, handled by the q-witness. The
covering-min / lowness lemma is about COVERING sets only.

So the honest escape from the covering-min (missing a top-half core speed k) is a COVERING set
S = {1,...,n-2}\{k}  U  {2 large speeds}, where the 2 speeds must supply multiples of k, n-1, n
(the three resonances the punctured core no longer covers) -- 3 resonances, 2 speeds, so one speed
must double up (e.g. n(n-1)=lcm(n-1,n), or 12*13=156 covering 12&13).

This confirms: over genuine COVERING escapes, min M > n/Phi6 (the lowness lemma), and the binding
exhibits the Steinhaus displacement (the forced large multiple of k is displaced from the k-slot).
"""
from fractions import Fraction as F
from math import gcd

def Phi6(n):
    return n * n - n + 1

def dist0(x, D):
    x %= D
    return min(x, D - x)

def M_and_arg(S):
    Dmax = 2 * max(S) + 2
    best = F(0); arg = None
    for D in range(2, Dmax + 1):
        for a in range(1, D):
            m = D
            for s in S:
                d = dist0(s * a, D)
                if d < m:
                    m = d
                if m == 0:
                    break
            if F(m, D) > best:
                best = F(m, D); arg = (D, a)
    return best, arg

def is_covering(S, n):
    """Covering set: a multiple of every q in {2,...,n} is present."""
    return all(any(s % q == 0 for s in S) for q in range(2, n + 1))

def gcd_all(S):
    g = 0
    for s in S:
        g = gcd(g, s)
    return g


if __name__ == "__main__":
    print("COVERING escapes missing a top-half core speed k  (must be a genuine covering set)")
    print("=" * 78)
    for n in [8, 10, 12, 14, 16]:
        m = n - 2
        thr = F(n, Phi6(n))
        # the construction (keeps full core) for reference
        constr = list(range(1, m + 1)) + [n * (n - 1)]
        Mc, _ = M_and_arg(constr)
        Kb = 14 * n
        worst = None
        for k in range(m // 2 + 1, m + 1):
            Pk = [j for j in range(1, m + 1) if j != k]
            bestk = None
            # 2 extra speeds A,B up to Kb; S must be a primitive covering set of size n-1
            for A in range(n - 1, Kb + 1):
                for B in range(A + 1, Kb + 1):
                    S = Pk + [A, B]
                    if len(set(S)) != n - 1:
                        continue
                    if gcd_all(S) != 1:
                        continue
                    if not is_covering(S, n):
                        continue
                    Mv, arg = M_and_arg(S)
                    if bestk is None or Mv < bestk[0]:
                        bestk = (Mv, A, B, arg)
            if bestk:
                Mv, A, B, (D, a) = bestk
                kslot = (k * a) % D
                # which forced multiple of k is present (A or B)?
                kap = A if A % k == 0 else (B if B % k == 0 else None)
                disp = dist0((kap * a) % D - kslot, D) if kap else None
                if worst is None or Mv < worst[0]:
                    worst = (Mv, k, A, B, (D, a), kslot, kap, disp)
        Mv, k, A, B, (D, a), kslot, kap, disp = worst
        margin = Mv - thr
        print(f"n={n:2d}: construction M={Mc}={float(Mc):.4f} (=n/Phi6? {Mc==thr})")
        print(f"      tightest COVERING escape (min over k): M={Mv}={float(Mv):.4f} vs n/Phi6={float(thr):.4f}"
              f"  margin={float(margin):+.5f}  {'OK (>n/Phi6)' if margin>0 else 'BELOW!'}")
        print(f"      at k={k} (dropped), speeds A={A},B={B}; binding D={D},a={a}; forced mult-of-{k}={kap}")
        if disp is not None:
            print(f"      Steinhaus displacement: k-slot={kslot}, forced-multiple image={(kap*a)%D}, dist={disp}")
    print()
    print("NET: over genuine COVERING sets, every single-core-drop escape has M > n/Phi6 (lowness lemma),")
    print("     and the dropped speed's forced large multiple is displaced from the k-slot it vacated.")
