#!/usr/bin/env python3
"""
Second fully-independent cross-check, using a DIFFERENT algorithm:
 - Monte Carlo estimate of prof_t(B) (random x in [0,1)) to confirm the exact
   interval-scan profiles (no shared breakpoint logic).
 - Direct Monte Carlo of P_r(B) itself: sample x AND r far-runner phases
   theta_1..theta_r uniform in [0,1); a "sector" for far runner is uniform 0..6;
   estimate P( B-sectors(x) UNION runner-sectors == all 7 ).  This is EXACTLY the
   r->infinity decorrelated model (far runners' sectors are iid uniform), so it
   must converge to P_r(B).  Confirms the formula P_r=sum prof_t c_t independently.
"""
import random
from math import comb

DEN = 7
def sector(e, x):
    return int(((e * x) % 1.0) * DEN) % DEN

def prof_mc(B, N=4_000_000):
    cnt = [0]*7
    for _ in range(N):
        x = random.random()
        hit = len({sector(e, x) for e in B})
        cnt[7-hit] += 1
    return [c/N for c in cnt]

def c_t(t, r):
    return sum((-1)**i * comb(t,i) * (1 - i/7)**r for i in range(t+1))

def P_r_mc(B, r, N=4_000_000):
    """Direct decorrelated model: r far runners contribute iid uniform sectors."""
    good = 0
    for _ in range(N):
        x = random.random()
        sectors = {sector(e, x) for e in B}
        for _ in range(r):
            sectors.add(random.randrange(7))
        if len(sectors) == 7:
            good += 1
    return good / N

random.seed(12345)

tests = {
    (9,2): (0,1,2,3,4,5,6),
    (10,2): (0,1,2,3,4,5,6,7),
}
exact = {
    (9,2): 1268/5145,    # = 0.246453
    (10,2): 14368/36015, # = 0.398945
}

print("Monte-Carlo cross-check of P_r(B) at the claimed argmax cores")
print("="*70)
for (k,r), B in tests.items():
    # build P_r from MC profile (route 1)
    pr = prof_mc(B)
    P_from_prof = sum(pr[t]*c_t(t,r) for t in range(7))
    # direct decorrelated MC (route 2)
    P_direct = P_r_mc(B, r)
    print(f"\n(k={k},r={r})  B={B}")
    print(f"   exact P_r(B)            = {exact[(k,r)]:.6f}")
    print(f"   MC via prof*c_t         = {P_from_prof:.6f}")
    print(f"   MC direct decorrelated  = {P_direct:.6f}")
    print(f"   MC profile t=0..6       = {[round(p,4) for p in pr]}")
