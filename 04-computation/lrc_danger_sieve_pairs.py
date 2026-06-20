from fractions import Fraction as F
from itertools import combinations

delta = F(1,14)

def danger_measure(T):
    if not T: return F(1)
    bps = set([F(0), F(1)])
    for v in T:
        for k in range(0, v+1):
            for off in (delta, 1-delta, F(0)):
                t = (k+off)/v
                if 0 <= t <= 1: bps.add(t)
    bps = sorted(bps)
    tot = F(0)
    for i in range(len(bps)-1):
        x0,x1 = bps[i],bps[i+1]
        if x1<=x0: continue
        mid = (x0+x1)/2
        ok = all(((v*mid)%1 < delta or (v*mid)%1 > 1-delta) for v in T)
        if ok: tot += x1-x0
    return tot

# pairwise: depends on gcd and relation structure
print("=== PAIRS: meas{both in danger} ===")
print("indep upper bound for pair = (1/7)^2 =", F(1,49), "=", float(F(1,49)))
pairs = [(1,2),(1,3),(2,3),(1,7),(2,5),(3,5),(1,13),(6,7),(7,8),(1,14),(7,14),(1,8),(5,9)]
for (a,b) in pairs:
    m = danger_measure([a,b])
    print(f"pair ({a},{b}): {m} = {float(m):.5f}   ratio to 1/49 = {float(m*49):.3f}")
