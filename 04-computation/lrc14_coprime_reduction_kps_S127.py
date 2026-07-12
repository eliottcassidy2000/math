# -*- coding: utf-8 -*-
# kps-2026-07-11-S127cont45: the COPRIME-REDUCTION lemma for Route B (spread-bulk anti-concentration).
#
# CONTEXT (klein-S258 finish map, opus-S239): LRC(14) reduces to Route B -- every spread divisor-complete (DC)
# family clears at some non-14 q in [15,29] (=> M>1/14). opus-S238 analyzed PRIME moduli via fold-classes and
# found NO small-prime shortcut (the coprime family is FULL at a prime). THIS exploits COMPOSITE moduli:
#
# STRUCTURAL LEMMA (coprime reduction, CORRECTED -- the (A) check below caught the q|v edge case). At
# modulus q with danger arc {0,+-1}, at any UNIT multiplier p (gcd(p,q)=1), a runner v is SAFE UNLESS
# it is a multiple of q (then v*p===0 always). Precisely, with g=gcd(v,q):
#   - g=1 (coprime): danger iff p === +-v^{-1} (v*p===+-1);
#   - 1<g<q: SAFE (v*p===0 mod g so cannot be +-1; and q'|p fails for q'=q/g>1, a unit p, so not 0);
#   - g=q (q|v): v*p===0 always => DANGER at every p.
#   So at a unit p, bandCount = #{coprime runners at +-1} + #{runners divisible by q}. Hence:
#   >> CLEARING via a unit multiplier  <=>  (q nmid every v_i)  AND  the COPRIME-to-q runners miss a fold-class.
#   >> SUFFICIENT: if [q nmid v_i for all i]  AND  [#{i: gcd(v_i,q)=1} <= phi(q)/2 - 1]  then CLEARS.
# For DC families at composite q (2^a*3^b*...), MOST runners share a small factor with q => FEW coprime
# runners => the anti-concentration SHRINKS from 13 runners to a handful. This is the composite handle
# opus's prime analysis (full coprime family) could not see.
import random
from math import gcd
from functools import reduce
from sympy import totient

def is_DC(v):  return all(any(x % d == 0 for x in v) for d in range(2, 15))
def prim(v):   return reduce(gcd, v) == 1
def longest_run(v):
    v = sorted(set(v)); b = m = 1
    for i in range(1, len(v)):
        if v[i] == v[i-1] + 1: m += 1; b = max(b, m)
        else: m = 1
    return b
def clears_at(v, q):
    lo = -(-q // 14)
    return any(all(lo <= (vi * p) % q <= q - lo for vi in v) for p in range(1, q))
def phi(q): return int(totient(q))

def rand_spread_DC(rng, vmax):
    for _ in range(20000):
        v = sorted(rng.sample(range(1, vmax + 1), 13))
        if prim(v) and is_DC(v) and longest_run(v) <= 7:
            return v
    return None

def main():
    rng = random.Random(20260711)
    window = [q for q in range(15, 30) if q % 14]
    fams = []
    for vmax in (40, 60, 80, 120, 200):
        for _ in range(700):
            v = rand_spread_DC(rng, vmax)
            if v: fams.append(v)
    # dedup
    fams = [list(t) for t in {tuple(v) for v in fams}]
    print(f"random spread primitive DC families (varied Vmax<=200, longest-run<=7): {len(fams)}")

    # (A) VERIFY the CORRECTED lemma: at a unit p, a runner with 1 < gcd(v,q) < q (shares a factor but is
    #     NOT a multiple of q) is never in danger. (multiples of q ARE in danger -- excluded here.)
    viol = 0
    for v in fams[:400]:
        for q in window:
            for p in range(1, q):
                if gcd(p, q) != 1: continue
                for vi in v:
                    g = gcd(vi, q)
                    if 1 < g < q and (vi * p) % q in (0, 1, q - 1):
                        viol += 1
    print(f"(A) CORRECTED LEMMA CHECK: runner with 1<gcd(v,q)<q in danger at a UNIT multiplier: {viol} (expect 0)")

    # (B) the CORRECTED sufficient guarantee: some q with [q nmid every v_i] AND [#coprime <= phi(q)/2 - 1]
    thresh = {q: phi(q) // 2 - 1 for q in window}
    print(f"(B) per-q guarantee: q nmid every v_i  AND  #coprime <= phi(q)/2 - 1 = {thresh}")
    def guar_q(v, q):
        return all(vi % q != 0 for vi in v) and sum(1 for vi in v if gcd(vi, q) == 1) <= thresh[q]
    guaranteed = 0; witness_q = {}
    for v in fams:
        qs = [q for q in window if guar_q(v, q)]
        if qs:
            guaranteed += 1
            witness_q[qs[0]] = witness_q.get(qs[0], 0) + 1
    print(f"    provably cleared by the CORRECTED coprime-guarantee: {guaranteed}/{len(fams)} ({100*guaranteed/len(fams):.1f}%)")
    print(f"    first-witness-q distribution: {dict(sorted(witness_q.items()))}")

    # (C) min #coprime-to-q over the window per family (how small does the effective family get?)
    mincop = [min(sum(1 for vi in v if gcd(vi, q) == 1) for q in window) for v in fams]
    from collections import Counter
    print(f"(C) min-over-window #coprime-to-q distribution: {dict(sorted(Counter(mincop).items()))}")
    print(f"    => effective anti-concentration family shrinks from 13 to a median of {sorted(mincop)[len(mincop)//2]} coprime runners")

    # (D) residual: families NOT caught by the guarantee -- do they still clear (via fold-class miss / non-unit p)?
    notguar = [v for v in fams if not any(guar_q(v, q) for q in window)]
    stillclear = sum(1 for v in notguar if any(clears_at(v, q) for q in window))
    print(f"(D) residual (guarantee misses): {len(notguar)}; of these, clear anyway in window: {stillclear}/{len(notguar)}")
    if notguar:
        ex = notguar[0]
        cq = [q for q in window if clears_at(ex, q)]
        print(f"    example residual family: {ex}")
        print(f"       clears at q={cq}; #coprime-to-q per window q: {[(q,sum(1 for vi in ex if gcd(vi,q)==1)) for q in cq]}")

if __name__ == "__main__":
    main()
