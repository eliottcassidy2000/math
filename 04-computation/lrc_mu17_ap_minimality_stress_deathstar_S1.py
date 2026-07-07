#!/usr/bin/env python3
"""
death-star-2026-07-07-S1 -- STRESS-TEST the load-bearing claim: AP minimizes mu_1/7.

mu_1/7(E) = meas{x in [0,1] : circular maxgap of {frac(e_i x)} > 1/7}.  Route-1's density
floor is mu_1/7 >= mu_1/7(AP_k) = 477/1078 (k=13), CONDITIONAL on AP-minimality (opus-S130).
Since E[maxgap] AP-minimality was just REFUTED (HYP-4777), I re-verify mu_1/7 AP-minimality
EXACTLY (not just numeric descent that missed prim-sat for E[maxgap]).

Exact mu_1/7 via the same three-gap piecewise-linear breakpoints as E[maxgap]:
within each linear region maxgap(x)=c*x+b0, {maxgap>1/7} is a sub-interval (exact length).
"""
from fractions import Fraction as F
import random
from math import gcd
from functools import reduce

THR = F(1, 7)

def regions(E):
    """yield (u, v, c, b0) : on (u,v), maxgap(x) = c*x + b0 (linear, exact)."""
    E = list(E)
    kdenom = max(max(abs(j) for j in E), max(abs(a-b) for a in E for b in E))
    bps = set([F(0), F(1)])
    for d in range(1, kdenom+1):
        for m in range(1, d):
            bps.add(F(m, d))
    bps = sorted(bps)
    for a, b in zip(bps, bps[1:]):
        mid = (a+b)/2
        fl = {j: (j*mid).__floor__() for j in E}
        order = sorted(E, key=lambda j: (j*mid - fl[j]))
        n = len(order); gaps = []
        for s in range(n):
            j1 = order[s]; j2 = order[(s+1) % n]
            if s < n-1:
                c = F(j2-j1); b0 = F(-(fl[j2]-fl[j1]))
            else:
                c = F(order[0]-order[-1]); b0 = F(-(fl[order[0]]-fl[order[-1]])+1)
            gaps.append((c, b0))
        subbp = set([a, b])
        for i in range(n):
            ci, bi = gaps[i]
            for jx in range(i+1, n):
                cj, bj = gaps[jx]
                if ci != cj:
                    xc = (bj-bi)/(ci-cj)
                    if a < xc < b: subbp.add(xc)
        subbp = sorted(subbp)
        for u, v in zip(subbp, subbp[1:]):
            m2 = (u+v)/2
            cbest, bbest = max(gaps, key=lambda cb: cb[0]*m2+cb[1])
            yield u, v, cbest, bbest

def Emaxgap_exact(E):
    return sum((c*u+b0 + c*v+b0)/2*(v-u) for u, v, c, b0 in regions(E))

def mu17_exact(E):
    tot = F(0)
    for u, v, c, b0 in regions(E):
        # length of {x in (u,v): c*x+b0 > 1/7}
        if c == 0:
            if b0 > THR: tot += (v-u)
        elif c > 0:
            xc = (THR-b0)/c
            lo = xc if xc > u else u
            if lo < v: tot += (v-lo)
        else:
            xc = (THR-b0)/c
            hi = xc if xc < v else v
            if hi > u: tot += (hi-u)
    return tot

def stats_num(E, res):
    E = list(E); thr = 1.0/7; cnt = 0
    for r in range(res):
        x = (r + 0.5) / res
        ph = sorted((e * x) % 1.0 for e in E)
        mg = ph[0] + 1.0 - ph[-1]; prev = ph[0]
        for p in ph[1:]:
            g = p - prev
            if g > mg: mg = g
            prev = p
        if mg > thr: cnt += 1
    return cnt/res

def primitive(E):
    g = reduce(gcd, E); return tuple(sorted(e // g for e in E))

if __name__ == "__main__":
    K = 13
    AP = list(range(1, 14))
    mu_ap = mu17_exact(AP)
    print("=== mu_1/7 AP-minimality stress test, K=13 ===\n")
    print(f"mu_1/7(AP_13) EXACT = {mu_ap} = {float(mu_ap):.6f}   (opus-S130 claims 477/1078={float(F(477,1078)):.6f})")
    print(f"cross-check: match = {mu_ap == F(477,1078)}\n")

    # structured families: does any dip BELOW the AP for mu_1/7?
    cands = {
        "prim-sat 2*{1..12}+13": [2,4,6,8,10,12,13,14,16,18,20,22,24],
        "S57-adversarial-.42": [2,6,8,10,11,12,14,16,18,20,22,26,42],
        "min-E[mg] {..29}": [1,3,5,6,7,8,9,10,11,13,15,20,29],
        "{1..12,14}": list(range(1,13))+[14],
        "odd {1,3..25}": list(range(1,26,2)),
        "3*{1..12}+13": [3*i for i in range(1,13)]+[13],
        "{1,2,4,8,16,..}+fill": [1,2,4,8,16,32,64,3,5,7,9,11,13],
    }
    print(f"{'family':>28} {'mu_1/7 exact':>16} {'~':>8} {'vs AP':>10}")
    below = []
    for nm, E in cands.items():
        if len(set(E)) != K:
            print(f"{nm:>28}  (not 13 distinct, skip)"); continue
        mu = mu17_exact(E)
        tag = "BELOW AP!" if mu < mu_ap else "above/eq"
        if mu < mu_ap: below.append((nm, E, mu))
        print(f"{nm:>28} {str(mu):>16} {float(mu):>8.5f} {tag:>10}")

    # hard adversarial descent to MINIMIZE mu_1/7 (numeric screen, exact-confirm)
    print("\n-- adversarial descent minimizing mu_1/7 (numeric screen) --")
    random.seed(17); gb = (2.0, None)
    for trial in range(80):
        H = random.choice([14,16,20,26,34,48])
        E = sorted(random.sample(range(1, H+1), K)); cur = stats_num(E, max(3000, 200*max(E)))
        for step in range(70):
            i = random.randrange(K); new = random.randint(1, random.choice([16,26,40,60]))
            if new in E: continue
            cand = sorted(set(E[:i]+E[i+1:]+[new]))
            if len(cand) != K: continue
            c = stats_num(cand, max(3000, 200*max(cand)))
            if c < cur - 1.5e-3: E, cur = cand, c
        if cur < gb[0]: gb = (cur, primitive(E))
    winner = list(gb[1]); mu_w = mu17_exact(winner)
    print(f"  descent winner: {gb[1]}")
    print(f"    numeric mu={gb[0]:.4f}   EXACT mu_1/7={mu_w}={float(mu_w):.6f}")
    print(f"    vs AP {float(mu_ap):.6f}: {'BELOW AP -- AP-min REFUTED' if mu_w < mu_ap else 'AP-min holds (winner >= AP)'}")
    if below:
        print("\n  STRUCTURED families below AP:", [(nm, float(m)) for nm,_,m in below])
    else:
        print("\n  No structured family dipped below the AP for mu_1/7.")
