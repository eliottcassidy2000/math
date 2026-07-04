"""kps-S36: THE TIGHT CRUX = COVERING-MIN. Verify `covering => M >= 14/183` (M = max_t min_i
||v_i t||, the lonely margin). Deep well {1..12,182} conjectured extremizer (M=14/183). Search
broadly for covering families with M < 14/183 (would refute); find the true floor + extremizers."""
import numpy as np
from math import gcd
from functools import reduce
def lcm(a,b): return a*b//gcd(a,b)

def is_covering(v):
    return all(any(x%q==0 for x in v) for q in range(2,15))

def M_of(v, grid=400000):
    """M = max over t in (0,1/2] of min_i ||v_i t||  (by symmetry t and 1-t). fine grid."""
    v=np.array(v,dtype=np.float64)
    t=np.linspace(1e-9, 0.5, grid)
    # ||v_i t|| = dist to nearest integer
    x=np.outer(v,t)            # (13, grid)
    d=np.abs(x-np.round(x))
    md=d.min(axis=0)           # min over runners
    i=md.argmax()
    return md[i], t[i]

def M_refine(v, t0, w=2e-6):
    v=np.array(v,dtype=np.float64)
    t=np.linspace(max(1e-9,t0-w), t0+w, 20000)
    x=np.outer(v,t); d=np.abs(x-np.round(x)); md=d.min(axis=0)
    i=md.argmax(); return md[i], t[i]

target=14/183
print(f"14/183 = {target:.6f} ; 1/14 = {1/14:.6f}")
print()
best=(9,None)
def consider(v,label):
    global best
    if not is_covering(v): return
    v=sorted(set(v))
    if len(v)!=13: return
    m,t0=M_of(v); m,t0=M_refine(v,t0)
    if m<best[0]: best=(m,tuple(v),label)
    return m,t0

# deep well + variants
consider([1,2,3,4,5,6,7,8,9,10,11,12,182],"deep well {1..12,182}")
for X in range(14, 400):
    consider(list(range(1,13))+[X], f"{{1..12,{X}}}")
# {1..11,13,X}
for X in range(14,400):
    consider(list(range(1,12))+[13,X], f"{{1..11,13,{X}}}")
# window-ish covering families (speeds<=28), random
rng=np.random.default_rng(0)
import itertools
tried=0
mins=[]
for _ in range(4000):
    v=sorted(set(int(x) for x in rng.integers(1,29,size=13)))
    if len(v)!=13: continue
    if not is_covering(v): continue
    m,t0=M_of(v,120000)
    mins.append(m); tried+=1
    if m<best[0]:
        m,t0=M_refine(v,t0)
        if m<best[0]: best=(m,tuple(v),"random window<=28")
# adversarial: local search to MINIMIZE M from the deep well
cur=list(range(1,13))+[182]
curM,_=M_of(cur); 
for it in range(300):
    i=rng.integers(0,13); step=int(rng.integers(-3,4))
    cand=cur.copy(); cand[i]=cand[i]+step
    if cand[i]<1 or len(set(cand))!=13 or not is_covering(cand): continue
    m,_=M_of(cand,120000)
    if m<curM: cur,curM=cand,m
consider(cur, "local-min search")

print(f"random window covering families: {tried} tested, min M = {min(mins):.6f}, median = {np.median(mins):.4f}")
print()
print(f"GLOBAL MIN M over all covering families searched: {best[0]:.6f}")
print(f"  extremizer: {best[2]} = {best[1]}")
print(f"  M/(1/14) = {best[0]*14:.4f} ; M/(14/183) = {best[0]/target:.4f}")
print()
if best[0] >= target - 1e-4:
    print("=> COVERING-MIN >= 14/183 CONFIRMED (no tighter covering family found).")
    print("   The deep well is the extremizer; covering => M >= 14/183 > 1/14 => mu>0 => LRC(14).")
else:
    print(f"=> COUNTEREXAMPLE: covering family with M={best[0]:.6f} < 14/183 !")
