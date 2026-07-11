# -*- coding: utf-8 -*-
# kps-2026-07-11-S127cont27: ADVERSARIAL hunt for the hB5 floor.
# hB5 (the single Lean obligation): every residual covering family has some q with B5(v,q) > 0.
# If max_q B5(v,q) <= 0 for some residual family, hB5 is FALSE at depth-5 (would need higher depth) -- a
# potential counterexample. Hill-climb to MINIMIZE max_q B5 over covering families (distinct, some |v|>=23),
# to find the true floor and whether depth-5 always suffices.
from math import comb
import random

def bandCount(v, q, p):
    c = 0
    for vi in v:
        r = (vi * p) % q
        if not (q <= 14 * r <= 13 * q): c += 1
    return c
def B5(v, q):
    S = [0]*6
    for p in range(1, q):
        b = bandCount(v, q, p)
        for d in range(6): S[d] += comb(b, d)
    return sum((-1)**d * S[d] for d in range(6))
def maxB5(v, qmax=320):
    best = -10**9; bq = 0
    for q in range(14, qmax+1):
        b = B5(v, q)
        if b > best: best = b; bq = q
    return best, bq
def covering_ok(v):
    # residual-ish: distinct |.|, some |v|>=23, no easy small-q witness at q<=13 (ratio>13 not peelable trivially)
    if len(set(abs(x) for x in v)) != 13: return False
    if max(abs(x) for x in v) < 23: return False
    return True

def main():
    random.seed(2026)
    # seed from the known tight near-AP shapes and hill-climb to minimize maxB5
    seeds = [list(range(1,13))+[26], list(range(1,12))+[13,26], list(range(1,12))+[13,36],
             [1,2,3,4,5,6,7,8,9,10,11,12,23]]
    globalmin = (10**9, None, 0)
    for seed in seeds:
        v = sorted(seed)
        cur, cq = maxB5(v)
        for it in range(140):
            i = random.randrange(13)
            cand = v[:]
            cand[i] = v[i] + random.choice([-2,-1,1,2,3,-3])
            if cand[i] < 1: continue
            cand = sorted(cand)
            if not covering_ok(cand): continue
            m, mq = maxB5(cand)
            if m < cur:                       # descend toward smaller max-B5 (harder family)
                v, cur, cq = cand, m, mq
        if cur < globalmin[0]:
            globalmin = (cur, v, cq)
        print(f"  seed {seed[:4]}...: min maxB5 reached = {cur} (v={v}, @q={cq})")
    print()
    m,v,q = globalmin
    print(f"  GLOBAL adversarial floor found: max_q B5 = {m}  at v={v} (@q={q})")
    print(f"  hB5 (some q with B5>0) holds on this family? {'YES (B5>0)' if m>0 else 'VIOLATED at depth-5 -- '+str(v)}")
    # also random adversaries
    worst = (10**9, None, 0)
    tested = 0
    for _ in range(250):
        base = sorted(random.sample(range(1,40), 12)) + [random.randint(23,60)]
        v = sorted(set(base))
        if len(v) != 13 or not covering_ok(v): continue
        tested += 1
        mm, mq = maxB5(v)
        if mm < worst[0]: worst = (mm, v, mq)
    print(f"\n  {tested} random residual families: worst (min) max_q B5 = {worst[0]} at v={worst[1]} (@q={worst[2]})")
    print(f"  => depth-5 B5>0 achieved on ALL tested ? {'YES' if min(m,worst[0])>0 else 'NO'}")
    print(f"  FLOOR = {min(m,worst[0])} (> 0 => hB5 empirically robust; near-AP shapes are the binding case)")

if __name__ == '__main__':
    main()
