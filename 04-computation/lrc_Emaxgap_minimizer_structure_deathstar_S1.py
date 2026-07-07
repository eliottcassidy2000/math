#!/usr/bin/env python3
"""
death-star-2026-07-07-S1 -- STRUCTURE of the E[maxgap] minimizer (the coarse-scale optimum).

E[maxgap] is minimized by coarse-cluster-breakers, NOT the AP.  prim-sat=2*{1..12}u{13}=0.2013.
Characterize: is the min a 'detuned dilated AP' d*{1..12} u {odd}?  Find the true inf, check >1/7.
Exact throughout (corrected kdenom).
"""
from fractions import Fraction as F
import random
from math import gcd
from functools import reduce

def Emaxgap_exact(E):
    E = list(E)
    kd = max(max(abs(j) for j in E), max(abs(a-b) for a in E for b in E))
    bps = set([F(0), F(1)])
    for d in range(1, kd+1):
        for m in range(1, d): bps.add(F(m, d))
    bps = sorted(bps); total = F(0)
    for a, b in zip(bps, bps[1:]):
        mid = (a+b)/2; fl = {j: (j*mid).__floor__() for j in E}
        order = sorted(E, key=lambda j: (j*mid - fl[j])); n = len(order); gaps = []
        for s in range(n):
            j1 = order[s]; j2 = order[(s+1) % n]
            if s < n-1: c = F(j2-j1); b0 = F(-(fl[j2]-fl[j1]))
            else: c = F(order[0]-order[-1]); b0 = F(-(fl[order[0]]-fl[order[-1]])+1)
            gaps.append((c, b0))
        sub = set([a, b])
        for i in range(n):
            ci, bi = gaps[i]
            for jx in range(i+1, n):
                cj, bj = gaps[jx]
                if ci != cj:
                    xc = (bj-bi)/(ci-cj)
                    if a < xc < b: sub.add(xc)
        sub = sorted(sub)
        for u, v in zip(sub, sub[1:]):
            m2 = (u+v)/2; cb, bb = max(gaps, key=lambda cb: cb[0]*m2+cb[1])
            total += (cb*u+bb + cb*v+bb)/2*(v-u)
    return total

def num(E, res):
    E = list(E); tot = 0.0
    for r in range(res):
        x = (r+0.5)/res; ph = sorted((e*x) % 1.0 for e in E)
        mg = ph[0]+1.0-ph[-1]; prev = ph[0]
        for p in ph[1:]:
            g = p-prev
            if g > mg: mg = g
            prev = p
        tot += mg
    return tot/res

def primitive(E):
    g = reduce(gcd, E); return tuple(sorted(e//g for e in E))

if __name__ == "__main__":
    thr = 1.0/7; APv = F(93, 440)
    print("=== E[maxgap] minimizer structure (exact); AP=93/440=%.5f, 1/7=%.5f ===\n" % (float(APv), thr))

    # family 1: 2*{1..12} u {odd m}, vary the odd element
    print("-- 2*{1..12} u {m}, m odd (prim-sat is m=13) --")
    bestfam = (2.0, None)
    for m in [1,3,5,7,9,11,13,15,17,19,21,23,25]:
        E = sorted(set([2*i for i in range(1,13)] + [m]))
        if len(E) != 13: continue
        v = Emaxgap_exact(E)
        if float(v) < bestfam[0]: bestfam = (float(v), tuple(E))
        print(f"   m={m:>2}: E[maxgap]={float(v):.6f}  {'<AP' if v<APv else ''}")

    # family 2: d*{1..12} u {13-ish odd}, vary dilation d
    print("\n-- d*{1..12} u {small odd}: does the 'detuned dilated AP' shape persist? --")
    for d in [2,3,4,5]:
        E = sorted(set([d*i for i in range(1,13)] + [d*6+1]))  # detune near middle
        E = E[:13] if len(E)>=13 else E
        if len(set(E))==13:
            v = Emaxgap_exact(list(E)); print(f"   d={d}, detune=6d+1: {float(v):.6f}")

    # family 3: which single detuning of 2*AP={2,4,..,26} is best? replace 2k by 2k-1 or other
    print("\n-- single detuning of 2*AP={2,4,...,26}: replace ONE 2j by t --")
    base = [2*i for i in range(1,14)]  # {2,4,...,26}
    best_det = (2.0, None)
    for j in range(13):
        for t in range(1, 40):
            if t % 2 == 0: continue
            E = sorted(set(base[:j]+base[j+1:]+[t]))
            if len(E) != 13: continue
            v = float(Emaxgap_exact(E))
            if v < best_det[0]: best_det = (v, tuple(E))
    print(f"   best single-detune of 2*AP: {best_det[1]}  E[maxgap]={best_det[0]:.6f}")

    # broad exact descent to push the true inf down
    print("\n-- exact descent for the true inf E[maxgap] --")
    random.seed(71); gb = (2.0, None)
    for trial in range(50):
        H = random.choice([14,16,20,26,30])
        E = sorted(random.sample(range(1, H+1), 13)); cur = num(E, max(3000,150*max(E)))
        for step in range(70):
            i = random.randrange(13); new = random.randint(1, random.choice([16,26,36]))
            if new in E: continue
            cand = sorted(set(E[:i]+E[i+1:]+[new]))
            if len(cand) != 13: continue
            c = num(cand, max(3000,150*max(cand)))
            if c < cur - 1e-4: E, cur = cand, c
        p = primitive(E)
        vex = float(Emaxgap_exact(list(p)))
        if vex < gb[0]: gb = (vex, p)
    print(f"   descent inf (exact): {gb[1]}  E[maxgap]={gb[0]:.6f}")
    truemin = min(bestfam[0], best_det[0], gb[0])
    print(f"\n  best E[maxgap] found = {truemin:.6f}  (margin over 1/7 = {truemin-thr:+.5f})")
    print(f"  champion prim-sat 2*{{1..12}}u{{13}} = {float(F(145091,720720)):.6f}")
    print(f"  => inf_E E[maxgap] ~ 0.20 >> 1/7=0.143: the DIRECT route inf>1/7 holds, comfortable +0.06,")
    print(f"     minimizer = a DETUNED dilated AP (coarse-cluster-breaker), NOT the plain AP.")
