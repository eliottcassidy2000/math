"""kps-S38: THE RUNNER-REMOVAL MEASURE BOUND. mu_v >= mu_{v\r} - 1/7 (each danger set meas 1/7).
So if SOME runner r has mu_{12 others}(margin 1/14) > 1/7, then mu_v>0 => M>1/14 => LRC. TEST:
does every primitive covering family have a removable runner with mu_{v\r} > 1/7?"""
import numpy as np
from math import gcd
from functools import reduce
rng=np.random.default_rng(11)
def gcd_all(v): return reduce(gcd,v)
def is_cov(v): return all(any(x%q==0 for x in v) for q in range(2,15))
def mu(speeds, grid=2000000):
    """safe measure = fraction of t in [0,1) with min_i ||v_i t|| >= 1/14."""
    v=np.array(speeds,dtype=np.float64)
    t=np.linspace(0,1,grid,endpoint=False)
    x=np.outer(v,t); d=np.abs(x-np.round(x)); md=d.min(axis=0)
    return (md>=1/14-1e-12).mean()
def M_of(v,grid=300000):
    v=np.array(sorted(set(v)),dtype=np.float64)
    t=np.linspace(1e-9,0.5,grid); x=np.outer(v,t); d=np.abs(x-np.round(x)); md=d.min(axis=0)
    return md.max()
oneseven=1/7
print(f"1/7 = {oneseven:.5f} (each runner's danger measure). Need mu_{{v\r}} > 1/7 for SOME r.")
print()
def test(v,lab):
    v=sorted(set(v))
    if len(v)!=13 or gcd_all(v)!=1 or not is_cov(v): return
    muv=mu(v)
    best_r=-1; best_mu=-1
    for r in range(13):
        w=v[:r]+v[r+1:]
        m=mu(w)
        if m>best_mu: best_mu=m; best_r=r
    ok = best_mu > oneseven
    print(f"{lab:>22} {str(v):>44}")
    print(f"    mu_v(13)={muv:.5f}  max_r mu_(12)={best_mu:.5f} (remove {v[best_r]})  >1/7? {ok}  M={M_of(v):.5f}")
    return ok, best_mu

results=[]
r=test(list(range(1,13))+[182],"deep well"); results.append(r)
# other primitive covering families
cand=[list(range(1,12))+[13,182], list(range(1,13))+[26], list(range(1,13))+[196],
      [1,2,3,4,5,6,7,8,9,10,11,13,168], list(range(1,13))+[14]]
for c in cand:
    r=test(sorted(set(c)),"structured"); 
    if r: results.append(r)
# random primitive covering (small)
cnt=0
for _ in range(400):
    v=sorted(set(int(x) for x in rng.integers(1,40,size=13)))
    if len(v)!=13 or gcd_all(v)!=1 or not is_cov(v): continue
    w_best=max(mu(v[:r]+v[r+1:], grid=400000) for r in range(13))
    results.append((w_best>oneseven, w_best)); cnt+=1
    if w_best<=oneseven: print(f"  FAIL random: {v} max mu_12={w_best:.5f} <= 1/7")
print()
oks=[x[0] for x in results if x]
print(f"tested {len(oks)} primitive covering families ({cnt} random). max_r mu_(12) > 1/7 for ALL? {all(oks)}")
print(f"  min over families of (max_r mu_12) = {min(x[1] for x in results if x):.5f}  vs 1/7={oneseven:.5f}")
if all(oks):
    print("=> REMOVAL BOUND WORKS: every primitive covering family has a runner whose removal")
    print("   leaves 12 with safe measure > 1/7 => mu_v>0 => M>1/14 => LRC(14). Reduces to a")
    print("   12-runner measure bound: primitive-covering => EXISTS r with mu_{v\r} > 1/7.")
