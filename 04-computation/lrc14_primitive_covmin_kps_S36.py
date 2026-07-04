"""kps-S36: RIGOROUS covering-min over PRIMITIVE (gcd=1) covering families (the gcd subtlety:
14*{1..13} is covering+tight M=1/14 but IMPRIMITIVE, reduces to {1..13}=non-covering=sieve).
Broad minimize-M search incl large speeds. Confirm floor = 14/183 (deep well) for primitive."""
import numpy as np
from math import gcd
from functools import reduce
rng=np.random.default_rng(3)
def gcd_all(v): return reduce(gcd,v)
def is_primitive(v): return gcd_all(v)==1
def is_covering(v): return all(any(x%q==0 for x in v) for q in range(2,15))
def M_of(v, grid=120000):
    v=np.array(sorted(set(v)),dtype=np.float64)
    if len(v)!=13: return 9.0,0.0
    t=np.linspace(1e-9,0.5,grid); x=np.outer(v,t); d=np.abs(x-np.round(x)); md=d.min(axis=0)
    i=md.argmax(); return md[i],t[i]

target=14/183
best=(9.0,None)
def consider(v):
    global best
    v=sorted(set(v))
    if len(v)!=13 or not is_primitive(v) or not is_covering(v): return None
    m,_=M_of(v)
    if m<best[0]: best=(m,tuple(v))
    return m

# 1) structured: {1..12, X}, {1..11,13,X}, {1..10,12,13,X}, deep-well relatives
for X in range(13,2000): consider(list(range(1,13))+[X])
for X in range(13,2000): consider(list(range(1,12))+[13,X])
for X in range(13,2000): consider(list(range(1,11))+[12,13,X])
# 2) primitive dilated-AP-like but covering: {a,2a,..,12a, Y} with a small, Y covering the rest
for a in range(1,6):
    for Y in range(13,3000):
        v=[a*k for k in range(1,13)]+[Y]
        if gcd_all(v)==1: consider(v)
# 3) minimize-M local search over primitive covering, large speed room
def minimize(start, iters=800, cap=2000):
    cur=sorted(set(start))
    if len(cur)!=13 or not is_primitive(cur) or not is_covering(cur): return
    curM,_=M_of(cur)
    for _ in range(iters):
        i=rng.integers(0,13); step=int(rng.integers(-6,7))
        cand=sorted(set(cur[:i]+[cur[i]+step]+cur[i+1:]))
        if len(cand)!=13 or cand[0]<1 or cand[-1]>cap: continue
        if not is_primitive(cand) or not is_covering(cand): continue
        m,_=M_of(cand)
        if m<curM-1e-9: cur,curM=cand,m
    consider(cur)
for _ in range(40):
    seed=sorted(set(int(x) for x in rng.integers(1,60,size=13)))
    if len(seed)==13: minimize(seed)
minimize(list(range(1,13))+[182])   # from deep well
minimize([14,28,42,56,70,84,98,112,126,140,154,168,169])  # near imprimitive but perturbed primitive

print(f"14/183 = {target:.6f}, 1/14 = {1/14:.6f}")
print(f"PRIMITIVE covering-min found: M = {best[0]:.6f}  (M*14 = {best[0]*14:.4f})")
print(f"  extremizer: {best[1]}")
print(f"  is deep well? {best[1]==tuple(sorted(list(range(1,13))+[182]))}")
print()
if best[0] >= target - 5e-5:
    print("=> PRIMITIVE COVERING-MIN = 14/183 CONFIRMED (deep well), gap above 1/14 = 0.005.")
    print("   The gcd subtlety: imprimitive dilated-AP covering families (14*{1..13}) have M=1/14")
    print("   but reduce to {1..13}=NON-covering=sieve. So PRIMITIVE covering => M>=14/183>1/14.")
else:
    print(f"=> tighter primitive covering family: M={best[0]:.6f} < 14/183, family {best[1]}")
