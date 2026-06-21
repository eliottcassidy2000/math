#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_q108_threadB_corner_monotone_kpswf4.py   (kind-pasteur 2026-06-21, THREAD B / part c)

The L_y survival-basis decomposition (opus routeB, mac-mini U4=1-G_1+G_5+4G_6):
   L_y functional = sum_t phi(t) p_t  with  phi(t) = 1 - t + sum_{a=1..5} w_a (t-a)_+ ,
where t = N = #missed inner sectors and p_t = P(N=t).  Equivalently
   L_y(E) = 1 - E[N] + sum_a w_a * C_a(E),     C_a(E) := E[(N-a)_+].
The corner functionals C_a are the "monotone cuts".  EMPIRICAL (routeB):
   - C_a is CONSEC-MAXIMIZED for a >= 3 (k=8) and a >= 2 (k=9,10): 0 exceeders.
   - C_1 (= E[(N-1)_+]) is NOT consec-maximized (458/75/15 exceeders).
   - E[N] is consec-MINIMIZED.
So the win is NOT a single monotone functional: the LINEAR term (-E[N]) wants consec
(min mean miss), the LOW corner C_1 fights consec, the HIGH corners C_a (a>=3) want
consec.  The L_y WEIGHTS make the signed sum come out consec-max.

THIS SCRIPT rigorously characterizes each piece:
 (1) E[N] = S_1 = sum over inner sectors of P(s missed).  Is consec the UNIQUE minimizer
     of E[N] (S_1)?  Test exhaustively.  [If yes -> the linear part is a clean lemma.]
 (2) High corner C_a (a>=3): is consec the UNIQUE maximizer?  Test exhaustively k=8,9,10.
     C_a(E) = E[(N-a)_+] = sum_{t>a} (t-a) p_t = "deep-miss mass": measure-weighted excess
     of missed sectors beyond a.  Consec wins => consec produces the HEAVIEST deep-miss
     tail.  This is the OPPOSITE intuition to "consec equidistributes": consec ALSO has
     the largest probability of MANY simultaneous misses.  Confirm + locate the mechanism.
 (3) Is there a SINGLE monotone scalar (a convex combination of the corners with the
     CORRECT signs) that is BOTH consec-extremal AND closes the cap?  i.e. is the L_y
     vector phi the minimal such certificate, or is there a cleaner all-nonneg one?
 (4) HONEST QUESTION: can the consec-max of L_y be REDUCED to:
        [E[N] consec-min]  AND  [each high-corner C_a (a in active set) consec-max] ?
     i.e. does dominance hold COORDINATEWISE on (-E[N], C_3, C_4, C_5) so that any
     nonneg-weighted combo is consec-max -- EXCEPT the C_1 term?  Quantify how much the
     C_1 anomaly costs and whether the L_y weight on C_1 (w_1) is small enough that the
     other terms dominate.
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)

def dist_p(E):
    """p[t] = measure{x : exactly t of the 6 inner sectors {1..6} are missed by E}.
       (sector 0 always hit since 0 in E.)  Returns p[0..6]."""
    E = sorted(set(E))
    bps = {F(0), F(1)}
    for e in E:
        if e == 0:
            continue
        for a in range(0, 7 * e + 1):
            bps.add(F(a, 7 * e))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    p = [F(0)] * 7
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        hit = set()
        for e in E:
            v = (e * mid) % 1
            hit.add((v.numerator * 7) // v.denominator)
        t = sum(1 for j in range(1, 7) if j not in hit)
        p[t] += (hi - lo)
    return p

def EN(p):
    return sum(t * p[t] for t in range(7))

def corner(p, a):
    """C_a = E[(N-a)_+] = sum_{t>a} (t-a) p_t."""
    return sum((t - a) * p[t] for t in range(a + 1, 7))

def measS7_from_p(p):
    return p[0]

def consec(k):
    return list(range(k))

def primitive(E):
    return reduce(gcd, [e for e in E if e != 0], 0) == 1

def main():
    print("#" * 84)
    print("# THREAD B part(c): the survival-basis CORNER decomposition of L_y")
    print("#" * 84)
    print("  L_y(E) = 1 - E[N] + sum_a w_a C_a(E),  C_a = E[(N-a)_+].")
    print("  Question: is consec-max of L_y reducible to clean per-corner monotone lemmas?")

    for k in (8, 9, 10):
        print("\n" + "=" * 84)
        print(f"k = {k}")
        print("=" * 84)
        C = consec(k)
        pc = dist_p(C)
        ENc = EN(pc)
        cornersC = {a: corner(pc, a) for a in range(1, 6)}
        print(f"  consec[0..{k-1}]: measS7={float(pc[0]):.5f}  E[N]={float(ENc):.5f}")
        print(f"    corners C_a: " + ", ".join(f"C_{a}={float(cornersC[a]):.5f}" for a in range(1, 6)))

        # bank: all primitive k-sets with 0, in {0..span}
        span = k + 5
        bank = [[0] + list(r) for r in itertools.combinations(range(1, span + 1), k - 1)]
        bank = [E for E in bank if primitive(E)]
        # (1) E[N] minimizer
        en_lower = 0
        en_min = ENc
        en_min_set = C
        # (2) corner maximizers
        corner_higher = {a: 0 for a in range(1, 6)}
        corner_max = {a: cornersC[a] for a in range(1, 6)}
        corner_max_set = {a: C for a in range(1, 6)}
        # (4) coordinatewise dominance on (-E[N], C_3,C_4,C_5): does ANY set beat consec
        #     on ALL of {lower E[N]} and {>= each high corner}? (i.e. Pareto-dominate)
        pareto_beats = 0
        for E in bank:
            p = dist_p(E)
            en = EN(p)
            if en < ENc - F(1, 10**15):
                en_lower += 1
                if en < en_min:
                    en_min = en; en_min_set = E
            cor = {a: corner(p, a) for a in range(1, 6)}
            for a in range(1, 6):
                if cor[a] > cornersC[a] + F(1, 10**15):
                    corner_higher[a] += 1
                    if cor[a] > corner_max[a]:
                        corner_max[a] = cor[a]; corner_max_set[a] = E
            # Pareto on (-E[N], C_3, C_4, C_5): all >= consec and one strictly >
            ge = (en <= ENc) and all(cor[a] >= cornersC[a] for a in (3, 4, 5))
            strict = (en < ENc) or any(cor[a] > cornersC[a] for a in (3, 4, 5))
            if ge and strict:
                pareto_beats += 1

        print(f"  (1) E[N]: #shapes with LOWER E[N] than consec = {en_lower}  "
              f"({'consec is UNIQUE min' if en_lower == 0 else 'NOT min'})"
              f"  global-min E[N]={float(en_min)}")
        if en_lower:
            print(f"      lowest-E[N] set = {en_min_set}")
        for a in range(1, 6):
            ch = corner_higher[a]
            print(f"  (2) C_{a}=E[(N-{a})_+]: #shapes EXCEEDING consec = {ch}  "
                  f"({'consec UNIQUE max' if ch == 0 else 'NOT max'})"
                  + (f"  top set {corner_max_set[a]} C={float(corner_max[a]):.5f}" if ch else ""))
        print(f"  (4) Pareto on (-E[N], C_3,C_4,C_5): #shapes weakly-dominating consec & "
              f"strict somewhere = {pareto_beats}  "
              f"({'consec PARETO-OPTIMAL on these 4 coords' if pareto_beats == 0 else 'dominated!'})")

    print("\n" + "=" * 84)
    print("INTERPRETATION")
    print("=" * 84)
    print("  If for every k: E[N] consec-MIN (unique) AND C_a (a>=3) consec-MAX (unique),")
    print("  AND consec is Pareto-optimal on (-E[N],C_3,C_4,C_5), then EVERY nonneg-weighted")
    print("  combination  -alpha*E[N] + sum_{a>=3} w_a C_a  (alpha,w_a >= 0) is consec-max.")
    print("  The ONLY obstruction to a clean monotone proof is the LOW corner C_1 (and the")
    print("  +C_2 term at k=8): L_y puts weight w_1>0 on C_1 where consec is NOT max. So the")
    print("  consec-max of L_y is NOT coordinatewise -- it needs the C_1 deficit to be")
    print("  OUTWEIGHED by the (-E[N]) + high-corner surplus.  Quantify that net below.")

    # Quantify the C_1 anomaly cost vs the rest, on the WORST C_1-exceeding shapes.
    for k in (8,):
        print(f"\n  --- net-balance audit at k={k} (the binding even-band case) ---")
        C = consec(k); pc = dist_p(C)
        # L_y weights (from routeB): k=8 -> phi=(1,0,0,1/10,0,0,1); decompose
        phi = [F((t-1)*(t-2)*(t-4)*(t-5), 40) for t in range(7)]
        w = {a: phi[a-1] - 2*phi[a] + phi[a+1] for a in range(1, 6)}
        beta = phi[1] - phi[0]
        Lc = sum(pc[t]*phi[t] for t in range(7))
        span = k + 5
        bank = [[0]+list(r) for r in itertools.combinations(range(1, span+1), k-1)]
        bank = [E for E in bank if primitive(E)]
        # find shape with max C_1 surplus; show its full term ledger vs consec
        worst = None
        for E in bank:
            p = dist_p(E)
            dc1 = corner(p,1) - corner(pc,1)
            if worst is None or dc1 > worst[0]:
                worst = (dc1, E, p)
        dc1, E, p = worst
        print(f"      max-C_1 shape {E}: C_1 surplus over consec = {float(dc1):.5f}")
        print(f"      term ledger (L_y term = weight * functional), shape minus consec:")
        print(f"        d(-E[N]) contribution: {float(beta*(EN(p)-EN(pc))):+.6f}  (beta={beta})")
        for a in range(1, 6):
            dca = corner(p,a) - corner(pc,a)
            print(f"        w_{a}={str(w[a]):>5} * dC_{a}={float(dca):+.6f}  -> {float(w[a]*dca):+.6f}")
        L_E = sum(p[t]*phi[t] for t in range(7))
        print(f"      total L_y(shape)-L_y(consec) = {float(L_E - Lc):+.6f}  "
              f"(<0 confirms consec still wins despite C_1 deficit)")

    print("\nDONE.")

if __name__ == "__main__":
    main()
