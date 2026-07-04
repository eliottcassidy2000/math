"""kps-S37 (fast): tight locus + KEY test: do tight families optimize on the 14-grid t=k/14?
That is the hinge of the rigidity proof: primitive-tight => optimal on 14-grid; covering =>
14-grid unsafe (runner div 14 at observer) => covering can't be primitive-tight => M>1/14."""
import numpy as np
from math import gcd
from functools import reduce
def gcd_all(v): return reduce(gcd,v)
def is_cov(v): return all(any(x%q==0 for x in v) for q in range(2,15))
def M_and_t(v, grid=60000):
    v=np.array(sorted(set(int(x) for x in v)),dtype=np.float64)
    if len(v)!=13: return 9.0,0.0
    t=np.linspace(1e-9,0.5,grid); x=np.outer(v,t); d=np.abs(x-np.round(x)); md=d.min(axis=0)
    i=md.argmax(); return md[i],t[i]
def on_14grid(t0):
    # is t0 approx k/14?
    return min(abs(t0-k/14) for k in range(1,8)) < 2e-3
inv14=1/14
tight=[]
def consider(v,lab):
    v=sorted(set(int(x) for x in v))
    if len(v)!=13: return
    m,t0=M_and_t(v)
    if m<=inv14+1.5e-4: tight.append((round(m,6),tuple(v),is_cov(v),gcd_all(v),round(t0,5),on_14grid(t0),lab))

# structured: APs and dilates
for a in range(1,13):
    for dd in range(1,13):
        consider([a+k*dd for k in range(13)], f"AP({a},{dd})")
for c in range(1,7): consider([c*k for k in range(1,14)], f"{c}*{{1..13}}")
# random small
rng=np.random.default_rng(7)
for _ in range(2500):
    v=tuple(sorted(rng.choice(range(1,30),size=13,replace=False)))
    consider(v,"rand")

uniq={}
for m,v,cov,g,t0,grid14,lab in tight: uniq[v]=(m,cov,g,t0,grid14,lab)
print(f"1/14={inv14:.6f}. tight (M<=1/14+1.5e-4): {len(uniq)}")
print(f"{'M':>9} {'cov':>4} {'gcd':>4} {'t*':>8} {'on14grid':>9}  family")
for v,(m,cov,g,t0,grid14,lab) in sorted(uniq.items(),key=lambda kv:(kv[1][0],kv[0])):
    print(f"{m:>9.6f} {str(cov):>4} {g:>4} {t0:>8.4f} {str(grid14):>9}  {list(v)}  [{lab}]")
prim=[(v,d) for v,d in uniq.items() if d[2]==1]
print()
print(f"PRIMITIVE tight: {len(prim)}")
for v,d in prim: print(f"   {list(v)}  M={d[0]:.6f} t*={d[3]:.4f} on14grid={d[4]} cov={d[1]}")
# KEY: primitive tight all {1..13}? all optimize on 14-grid? all non-covering?
allprim_1to13 = all(list(v)==list(range(1,14)) for v,_ in prim)
allprim_grid = all(d[4] for _,d in prim)
anyprim_cov = any(d[1] for _,d in prim)
print()
print(f"every primitive tight == {{1..13}}?  {allprim_1to13}")
print(f"every primitive tight optimizes on 14-grid?  {allprim_grid}")
print(f"any primitive tight covering?  {anyprim_cov}")
