#!/usr/bin/env python3
"""
Multi-straddle moment counts: how do several coprime return levels combine? (mac-mini-S151)
================================================================================
THM-1740: for a SINGLE straddle (one charge mult r, opposite mult 1), M* = r*m0,
m0 = (p+|n|)/gcd(p,|n|) the coprime charge-pair return level.  HYP-8560: what happens with
SEVERAL straddles, each with its own return level?

A straddle = a (pos charge p, neg charge n) pair; its level is L = r*(p+|n|)/gcd(p,|n|).
A pattern with charges {+p, -a, -b} has two straddles (p,a),(p,b) sharing +p; {+p,+q,-a} has
(p,a),(q,a) sharing -a.  We measure M* and compare it to candidate combination laws of the
per-straddle levels L_i:
    max L_i      (straddles independent, biggest wins)
    lcm L_i      (simultaneous return -- primes combine MULTIPLICATIVELY)
    sum L_i      (levels stack)
The answer tells us the arithmetic of how the moment ideal saturates across straddles.
"""
import sympy as sp
from math import factorial, gcd
from functools import reduce
def lcm(a, b): return a*b//gcd(a, b)

def moments(monos, coeffs, M):
    P = {}
    for i, (a, b) in enumerate(monos): P[(a, b)] = P.get((a, b), 0) + coeffs[i]
    def mul(X, Y):
        o = {}
        for (a1, b1), c1 in X.items():
            for (a2, b2), c2 in Y.items():
                k = (a1+a2, b1+b2); o[k] = o.get(k, 0) + c1*c2
        return o
    Pm = {(0, 0): sp.Integer(1)}; g = []
    for m in range(1, M+1):
        Pm = mul(Pm, P)
        g.append(sp.expand(sum(c*factorial(a) for (a, b), c in Pm.items() if a == b)))
    return g

def minM(monos, coeffs, Mmax):
    pos = [coeffs[i] for i, (a, b) in enumerate(monos) if a > b]
    neg = [coeffs[i] for i, (a, b) in enumerate(monos) if a < b]
    if not pos or not neg: return 0
    w = sp.Symbol('w'); allg = moments(monos, coeffs, Mmax)
    for M in range(1, Mmax+1):
        ok = True
        for cp in pos:
            for cn in neg:
                G = sp.groebner(list(allg[:M]) + [1 - w*cp*cn], *(list(coeffs)+[w]), order='grevlex')
                if not any(x.is_number and x != 0 for x in G.exprs): ok = False; break
            if not ok: break
        if ok: return M
    return None

def straddle_levels(monos):
    """per-straddle level L = mult(neg charge involved... )*(p+|n|)/gcd; here mult=1 patterns."""
    chg = [a-b for a, b in monos]
    from collections import Counter
    mult = Counter(chg)
    P = [c for c in set(chg) if c > 0]; N = [-c for c in set(chg) if c < 0]
    levels = {}
    for p in P:
        for n in N:
            g = gcd(p, n); base = (p+n)//g
            # multiplicity of the straddle: how many terms at the busier of the two charges
            r = max(mult[p], mult[-n])
            levels[(p, -n)] = r*base
    return levels

print("=" * 78)
print("MULTI-STRADDLE PATTERNS: M* vs max / lcm / sum of per-straddle levels")
print("=" * 78)
print(f"{'charges':>16} {'straddle levels':>26} {'M*':>4} {'max':>4} {'lcm':>5} {'sum':>4} {'fits':>8}")
pats = [
    [(2, 0), (0, 1), (0, 3)],       # +2,-1,-3 : (2,1)=3, (2,3)=5
    [(3, 0), (0, 1), (0, 2)],       # +3,-1,-2 : (3,1)=4, (3,2)=5
    [(2, 0), (0, 1), (0, 5)],       # +2,-1,-5 : (2,1)=3, (2,5)=7
    [(3, 0), (0, 1), (0, 4)],       # +3,-1,-4 : (3,1)=4, (3,4)=7
    [(2, 0), (3, 0), (0, 1)],       # +2,+3,-1 : (2,1)=3, (3,1)=4
    [(2, 0), (5, 0), (0, 1)],       # +2,+5,-1 : (2,1)=3, (5,1)=6
    [(3, 0), (5, 0), (0, 1)],       # +3,+5,-1 : (3,1)=4, (5,1)=6
    [(2, 0), (0, 1), (0, 2)],       # +2,-1,-2 : (2,1)=3, (2,2)=2
]
for monos in pats:
    coeffs = sp.symbols('c0:%d' % len(monos))
    Ms = minM(list(monos), list(coeffs), Mmax=20)
    ch = sorted({a-b for a, b in monos})
    lv = straddle_levels(monos); vals = sorted(lv.values())
    mx = max(vals); lc = reduce(lcm, vals); sm = sum(vals)
    fit = []
    if Ms == mx: fit.append("max")
    if Ms == lc: fit.append("lcm")
    if Ms == sm: fit.append("sum")
    print(f"{str(ch):>16} {str(dict(lv)):>26} {str(Ms):>4} {mx:>4} {lc:>5} {sm:>4} "
          f"{','.join(fit) if fit else 'OTHER':>8}")

print()
print("=" * 78)
print("MULTIPLICITY x MULTI-STRADDLE: does each straddle carry its own r?")
print("=" * 78)
print(f"{'pattern':>34} {'charges':>14} {'levels':>18} {'M*':>4} {'max':>4} {'lcm':>5}")
for monos, nm in (
    ([(2, 0), (0, 1), (1, 2), (0, 3)], "+2,-1x2(rad1,3),-3"),   # straddle(2,1) r=2 L=6, (2,3) L=5
    ([(3, 0), (0, 1), (0, 2), (1, 3)], "+3,-1,-2x2"),           # (3,1)=4,(3,2) r=2 L=10
):
    coeffs = sp.symbols('c0:%d' % len(monos))
    Ms = minM(list(monos), list(coeffs), Mmax=24)
    ch = sorted({a-b for a, b in monos})
    lv = straddle_levels(monos); vals = sorted(lv.values())
    mx = max(vals); lc = reduce(lcm, vals)
    print(f"{nm:>34} {str(ch):>14} {str(vals):>18} {str(Ms):>4} {mx:>4} {lc:>5}")

print()
print("SUMMARY")
print("  Read the 'fits' column: if M* = LCM of the per-straddle levels, the straddles combine")
print("  MULTIPLICATIVELY over primes -- a simultaneous-return / covering structure. If M* = max,")
print("  they are independent.  Either way this is the moment-ideal analogue of a return-time")
print("  lattice, and the combination law is the point of contact with the repo's LRC / covering")
print("  / three-gap threads (coprime speeds, first-return times, lcm of periods).")
