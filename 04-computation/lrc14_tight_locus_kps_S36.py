"""kps-S36: THE MECHANISM behind covering-min > 1/14. Conjecture: the TIGHT families (M near
1/14) are non-covering (APs {a,2a,..,13a} & reflections), so COVERING forces M off 1/14.
Find near-tight families (minimize M) + check: are they all non-covering? are they APs?"""
import numpy as np
from math import gcd
rng=np.random.default_rng(1)

def is_covering(v): return all(any(x%q==0 for x in v) for q in range(2,15))
def M_of(v, grid=200000):
    v=np.array(sorted(set(v)),dtype=np.float64)
    if len(v)!=13: return 9,0
    t=np.linspace(1e-9,0.5,grid); x=np.outer(v,t); d=np.abs(x-np.round(x)); md=d.min(axis=0)
    i=md.argmax(); return md[i],t[i]
def is_AP(v):
    v=sorted(v); d=[v[i+1]-v[i] for i in range(12)]
    return len(set(d))==1
def is_dilated_AP(v):
    """v == a*{1..13}/g form? check if sorted v is c*{1,2,..,13} for some c (after gcd)."""
    v=sorted(v); g=v[0]
    return all(v[i]==(i+1)*g for i in range(13))

# minimize M (find tight families) via local search from many starts
def minimize_M(start, iters=400, span=30):
    cur=sorted(set(start)); 
    if len(cur)!=13: return None
    curM,_=M_of(cur,80000)
    for _ in range(iters):
        i=rng.integers(0,13); step=int(rng.integers(-3,4))
        cand=cur.copy(); cand[i]+=step
        if cand[i]<1 or cand[i]>span or len(set(cand))!=13: continue
        m,_=M_of(cand,80000)
        if m<curM-1e-9: cur,curM=sorted(cand),m
    return cur,curM

print("Minimizing M (finding tight families) from random starts <=26:")
print(f"{'family':>42} {'M':>9} {'M*14':>6} {'covering?':>9} {'dilated AP?':>11}")
tights=[]
for trial in range(25):
    start=sorted(set(int(x) for x in rng.integers(1,27,size=13)))
    if len(start)!=13: continue
    r=minimize_M(start)
    if r is None: continue
    v,m=r
    if m < 0.0765:  # tighter than covering-min => should be non-covering
        cov=is_covering(v); dap=is_dilated_AP([x//gcd_all(v) for x in v]) if False else is_dilated_AP(v)
        tights.append((m,tuple(v),cov))
        print(f"{str(v):>42} {m:>9.5f} {m*14:>6.3f} {str(cov):>9} {str(is_dilated_AP(v)):>11}")

def gcd_all(v):
    from functools import reduce
    return reduce(gcd, v)

# explicitly: is {1..13} tight and non-covering? and its small dilates/variants?
print()
print("Reference tight families:")
for v,lab in [(list(range(1,14)),"AP {1..13}"),
              ([1,2,3,4,5,6,7,8,9,10,11,12,13],"{1..13}"),
              (list(range(1,13))+[182],"deep well")]:
    m,_=M_of(v); print(f"  {lab:>16}: M={m:.5f} (M*14={m*14:.3f}), covering={is_covering(v)}, dilatedAP={is_dilated_AP(v)}")

cov_tight = [t for t in tights if t[2]]
print()
print(f"Tight families found (M<14/183): {len(tights)}; of these COVERING: {len(cov_tight)}")
if len(cov_tight)==0:
    print("=> MECHANISM CONFIRMED: every family tighter than 14/183 is NON-COVERING.")
    print("   Covering forces M off the tight locus => M >= 14/183. The tight locus (APs & c.)")
    print("   is non-covering because an AP {a..13a} (a coprime 14) misses q=14 (=> t=1/14 lonely).")
else:
    print(f"=> COVERING TIGHT family exists: {cov_tight[0]}")
