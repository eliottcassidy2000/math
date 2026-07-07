#!/usr/bin/env python3
"""
death-star-2026-07-07-S1 -- TRUE minimizer of E[maxgap] AND of mu_1/7, K=13.

Findings to establish:
 (1) AP-minimality of E[maxgap] is FALSE (prim-sat 0.2013 < AP 0.2114). Find true inf, check >1/7.
 (2) Is mu_1/7 = P(maxgap>1/7) STILL AP-minimized?  (Different functional, maybe different min.)
     => reverse-Markov reduction may LOSE the AP-min structure that mu_1/7 has.
E[maxgap] and mu_1/7 are dilation-invariant. Numeric screen (pure python), exact-confirm winner.
"""
from fractions import Fraction as F
import random
from math import gcd
from functools import reduce

def Emaxgap_exact(E):
    E = list(E)
    kdenom = max(max(abs(j) for j in E), max(abs(a-b) for a in E for b in E))
    bps = set([F(0), F(1)])
    for d in range(1, kdenom+1):
        for m in range(1, d):
            bps.add(F(m, d))
    bps = sorted(bps)
    total = F(0)
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
            total += (cbest*u+bbest + cbest*v+bbest)/2*(v-u)
    return total

def stats_num(E, res):
    """returns (E[maxgap], mu_1/7) numerically."""
    E = list(E); thr = 1.0/7; sm = 0.0; cnt = 0
    for r in range(res):
        x = (r + 0.5) / res
        ph = sorted((e * x) % 1.0 for e in E)
        mg = ph[0] + 1.0 - ph[-1]; prev = ph[0]
        for p in ph[1:]:
            g = p - prev
            if g > mg: mg = g
            prev = p
        sm += mg
        if mg > thr: cnt += 1
    return sm/res, cnt/res

def primitive(E):
    g = reduce(gcd, E)
    return tuple(sorted(e // g for e in E))

def screen_res(E):
    return max(3000, min(24000, 200 * max(E)))

if __name__ == "__main__":
    K = 13; thr = 1.0/7; APmg = F(93,440)
    print("=== E[maxgap] & mu_1/7 minimizers, K=13 (1/7=%.5f, AP E[mg]=93/440=%.5f) ===\n" % (thr, float(APmg)))

    cands = {
        "AP {1..13}": list(range(1,14)),
        "prim-sat {even..24,13}": [2,4,6,8,10,12,13,14,16,18,20,22,24],
        "odd {1,3..25}": list(range(1,26,2)),
        "{1..12,14}": list(range(1,13))+[14],
        "{1..6,8..14}": list(range(1,7))+list(range(8,15)),
        "primes<=41": [2,3,5,7,11,13,17,19,23,29,31,37,41],
        "{2,4,6,8,10,12,14,16,18,20,22,24,13}": [2,4,6,8,10,12,14,16,18,20,22,24,13],
    }
    print(f"{'family':>36} {'E[mg]~':>8} {'mu17~':>7}  {'E[mg] flag':>10}")
    for nm, E in cands.items():
        if len(set(E)) != K: continue
        mg, mu = stats_num(E, screen_res(E))
        fl = "BELOW AP" if mg < float(APmg)-1e-4 else ("=AP" if abs(mg-float(APmg))<1e-4 else "above")
        print(f"{nm:>36} {mg:>8.5f} {mu:>7.4f}  {fl:>10}")

    # descent to MIN E[maxgap]
    print("\n-- descent: minimize E[maxgap] --")
    random.seed(3); gb = (2.0, None)
    for trial in range(60):
        H = random.choice([14,16,20,26,34])
        E = sorted(random.sample(range(1, H+1), K)); cur = stats_num(E, screen_res(E))[0]
        for step in range(60):
            i = random.randrange(K); new = random.randint(1, random.choice([16,26,40]))
            if new in E: continue
            cand = sorted(set(E[:i]+E[i+1:]+[new]))
            if len(cand) != K: continue
            c = stats_num(cand, screen_res(cand))[0]
            if c < cur - 1e-4: E, cur = cand, c
        if cur < gb[0]: gb = (cur, primitive(E))
    ex = Emaxgap_exact(list(gb[1]))
    print(f"  min-E[maxgap] winner: {gb[1]}")
    print(f"    numeric={gb[0]:.5f}  EXACT={ex}={float(ex):.6f}  margin_over_1/7={float(ex)-thr:+.5f}")

    # descent to MIN mu_1/7
    print("\n-- descent: minimize mu_1/7 --")
    random.seed(5); gm = (2.0, None)
    for trial in range(60):
        H = random.choice([14,16,20,26,34])
        E = sorted(random.sample(range(1, H+1), K)); cur = stats_num(E, screen_res(E))[1]
        for step in range(60):
            i = random.randrange(K); new = random.randint(1, random.choice([16,26,40]))
            if new in E: continue
            cand = sorted(set(E[:i]+E[i+1:]+[new]))
            if len(cand) != K: continue
            c = stats_num(cand, screen_res(cand))[1]
            if c < cur - 1e-4: E, cur = cand, c
        if cur < gm[0]: gm = (cur, primitive(E))
    mu_ap, _ = stats_num(list(range(1,14)), 24000)[::-1][0], None
    apmu = stats_num(list(range(1,14)), 24000)[1]
    print(f"  min-mu_1/7 winner: {gm[1]}  mu={gm[0]:.4f}   (AP mu_1/7~{apmu:.4f})")
    print(f"    => AP {'IS' if gm[0]>=apmu-2e-3 else 'is NOT'} the mu_1/7 minimizer (numeric).")
