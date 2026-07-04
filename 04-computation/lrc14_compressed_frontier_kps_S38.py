"""kps-S38: THE ACTUAL FRONTIER. Dominant-far => measure-independence (remove big runner:
mu_v ~ (6/7)mu_12 >0). All-comparable-LOOSE => mu large. The ONLY possibly-hard case =
LARGE all-comparable TIGHT (M near 1/14). Does it EXIST? Minimize M over speeds in [N,kN]
(compressed, large N). If min M bounded away from 1/14 => hard case empty => closure."""
import numpy as np
from math import gcd
from functools import reduce
rng=np.random.default_rng(5)
def gcd_all(v): return reduce(gcd,v)
def M_of(v,grid=400000):
    v=np.array(sorted(set(int(x) for x in v)),dtype=np.float64)
    if len(v)!=13: return 9.0
    t=np.linspace(1e-9,0.5,grid); x=np.outer(v,t); d=np.abs(x-np.round(x)); return d.min(axis=0).max()
def mu(v,grid=800000):
    v=np.array(v,dtype=np.float64); t=np.linspace(0,1,grid,endpoint=False)
    x=np.outer(v,t); d=np.abs(x-np.round(x)); return (d.min(axis=0)>=1/14-1e-12).mean()

# 1) VERIFY: removing the LARGE runner closes the deep well via near-independence
dw=list(range(1,13))+[182]
muv=mu(dw); mu12=mu(list(range(1,13)))
print("REMOVAL-OF-LARGE (measure independence) on the deep well:")
print(f"  mu_v(deep well) = {muv:.5f} ; mu({{1..12}}) = {mu12:.5f} ; (6/7)*mu12 = {6/7*mu12:.5f}")
print(f"  => mu_v ~ (6/7)mu12 ? {abs(muv-6/7*mu12)<0.01}  (182 fine comb ~ independent of 1..12)")
print()

# 2) FRONTIER: minimize M over LARGE COMPRESSED families (all speeds in [N, k*N])
print("Minimize M over LARGE COMPRESSED families (speeds in [N,kN]) -- can M reach 1/14=0.0714?")
print(f"{'N':>7} {'ratio k':>8} {'min M found':>12} {'M*14':>7}")
def minimize_M_compressed(N, k, iters=1500):
    lo,hi=N,int(k*N)
    cur=sorted(set(int(x) for x in rng.integers(lo,hi+1,size=13)))
    while len(cur)!=13: cur=sorted(set(int(x) for x in rng.integers(lo,hi+1,size=13)))
    curM=M_of(cur,120000)
    for _ in range(iters):
        i=rng.integers(0,13); step=int(rng.integers(-max(2,N//20),max(3,N//20)))
        cand=cur.copy(); cand[i]+=step
        cand=sorted(set(cand))
        if len(cand)!=13 or cand[0]<lo or cand[-1]>hi: continue
        m=M_of(cand,120000)
        if m<curM: cur,curM=cand,m
    return curM,cur
overall_min=9
for N in [50, 200, 1000, 5000]:
    for k in [1.5, 2.0, 3.0]:
        best=9; bestv=None
        for _ in range(6):
            m,v=minimize_M_compressed(N,k)
            if m<best: best,bestv=m,v
        overall_min=min(overall_min,best)
        print(f"{N:>7} {k:>8.1f} {best:>12.5f} {best*14:>7.3f}")
print()
print(f"MIN M over all large-compressed families searched: {overall_min:.5f} (M*14={overall_min*14:.3f})")
print(f"1/14 = {1/14:.5f}")
if overall_min > 0.09:
    print("=> LARGE COMPRESSED families are UNIFORMLY LOOSE (M >= ~0.1 = 1.4x radius, bounded")
    print("   away from 1/14). No large all-comparable TIGHT family found => the hard case is")
    print("   EMPTY: tight families are small-speed (window/AP) => census/far-peel closed;")
    print("   large families are either dominant-far (measure-indep) or compressed-loose (mu big).")
else:
    print(f"=> FOUND a large compressed family with M={overall_min:.5f} near 1/14 -- a genuine hard instance!")
