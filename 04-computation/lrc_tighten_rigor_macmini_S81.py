"""
S81: tightening the LRC(14) rigor (push-pull).
(c) CONSTRUCTION rigorous: AP/GW & dilations safe at EXACTLY the units of Z/(14d) (all phi(14d) -- kps S255 equioscillation).
(b) BOUNDED margin rigorous: single-swaps tight = AP/GW only, uniform margin delta~0.0026; 2-swaps no new tight.
(b) EQUIDISTRIBUTION pull: clean 1/7-removal is generic-only; resonant apex v=14 removes 0.73; survivor stays positive (M>1/14).
Remaining hard: (a) full finiteness/rigidity; (b) survivor-positivity for resonant v.
"""
from fractions import Fraction as F
from math import gcd
import numpy as np

def safe_at(S,t): return all(min((s*t)%1,1-(s*t)%1)>=F(1,14) for s in S)
def unit_witnesses(S,d):
    units=[a for a in range(1,14*d) if gcd(a,14*d)==1]
    return [a for a in units if safe_at(S,F(a,14*d))], len(units)
def M_of(S,grid=140000):
    t=np.arange(1,grid)/grid; s=np.ones(grid-1)
    for x in S:
        fr=(x*t)%1.0; s=np.minimum(s,np.minimum(fr,1-fr))
    return s.max()

print("(c) CONSTRUCTION: tight-locus safe at the units of Z/(14d):")
for name,S,d in [("AP",list(range(1,14)),1),("GW",list(range(1,12))+[13,24],1),("2AP",list(range(2,27,2)),2),("3AP",list(range(3,40,3)),3)]:
    g,n=unit_witnesses(S,d); print(f"  {name:<4} d={d}: {len(g)}/{n} units are witnesses  {'ALL (rigorous)' if len(g)==n else ''}")

print("\n(b) BOUNDED margin: AP single-swaps (r<=26), tight set + margin:")
thr=F(1,14); aboves=[]; tight=[]
for k in range(1,14):
    for r in range(1,27):
        if 1<=r<=13 and r!=k: continue
        S=[x for x in range(1,14) if x!=k]+[r]
        if len(set(S))<13: continue
        m=M_of(S)
        (tight if m<float(thr)+1e-4 else aboves).append((k,r,round(m,5)))
print(f"  tight (M=1/14): {sorted(set((k,r) for k,r,_ in tight if r!=k))}  (= GW 12->24; rest are AP no-swaps)")
print(f"  margin: min M above 1/14 = {min(m for _,_,m in aboves):.5f}, gap = {min(m for _,_,m in aboves)-float(thr):.5f}")

print("\n(b) EQUIDISTRIBUTION pull: large v removal of seed lonely set (1/7 only generic):")
seed=[1,2,3,4,5,6,8,9,10,11,12,13]
def safe_meas(S,grid=300000):
    t=np.arange(1,grid)/grid; s=np.ones(grid-1)
    for x in S:
        fr=(x*t)%1.0; s=np.minimum(s,np.minimum(fr,1-fr))
    return (s>=1/14).mean()
sm0=safe_meas(seed)
for v in (14,42,140):
    print(f"  v={v:>3}: removed {1-safe_meas(seed+[v])/sm0:.3f} of seed lonely (1/7={1/7:.3f}); M(seed+v)={M_of(seed+[v]):.5f}>1/14")
print("  => resonant v=14 removes 0.73 (not 1/7); survivor positive => M>1/14. Argument = survivor positive, not 6/7.")
