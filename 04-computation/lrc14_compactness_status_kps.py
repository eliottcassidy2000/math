#!/usr/bin/env python3
"""
lrc14_compactness_status_kps  (kind-pasteur, PROVE side)

PRECISE status of the COMPACTNESS LEMMA: "small L => bounded-lcm perturbation of
a tight config".

We found small-L configs with UNBOUNDED lcm (e.g. {1..11,13,150}, {1,3,4,5,7,8,9,
10,11,12,13,17,25} lcm 30M). So "bounded lcm" is FALSE as stated. BUT the relevant
question for inf L>0 is NOT bounded lcm of every small-L config; it is whether the
INFIMUM is attained / approached only by bounded configs. THM-518 decoupling says:
if one entry -> infinity, L -> (6/7)*meas(bounded core), which is BOUNDED BELOW by a
positive constant determined by the 12-element core. So large entries do NOT drive
L to 0.

This script makes the decoupling lower bound EXPLICIT and quantitative:
  For a config S = core(12 fixed elements) + {w}, as w->inf,
     L(S) -> (6/7) * meas(complement of danger(core)) = (6/7)*meas(G_core).
  We verify this for several cores and show the limit is >= 1/1260 always, i.e.
  pushing the perturbing element to infinity NEVER beats the bounded optimum.

Then we test 2-fold growth: S = (11 fixed) + {w, w'} with both large; does L stay
bounded below? We sample large coordinated pairs.

CONCLUSION target: state exactly which configs evade "bounded lcm" and why they
still cannot drive inf L below 1/1260 (decoupling floor).
"""
from fractions import Fraction as F
import math

def danger(v):
    out=[]; w=F(1,14*v)
    for k in range(v+1):
        lo=F(k,v)-w; hi=F(k,v)+w
        if lo<0: out += [(F(0),hi),(1+lo,F(1))]
        elif hi>1: out += [(lo,F(1)),(F(0),hi-1)]
        else: out.append((lo,hi))
    return [(x,y) for x,y in out if y>x]
def union(arcs):
    arcs=sorted((x,y) for x,y in arcs if y>x)
    if not arcs: return []
    res=[]; cl,ch=arcs[0]
    for lo,hi in arcs[1:]:
        if lo<=ch: ch=max(ch,hi)
        else: res.append((cl,ch)); cl,ch=lo,hi
    res.append((cl,ch)); return res
def total(a): return sum(y-x for x,y in union(a))
def complement(arcs):
    u=union(arcs); res=[]; prev=F(0)
    for lo,hi in u:
        if lo>prev: res.append((prev,lo))
        prev=max(prev,hi)
    if prev<1: res.append((prev,F(1)))
    return res
def L_exact(S):
    arcs=[]
    for v in set(S): arcs+=danger(v)
    return F(1)-total(arcs)
def L_float(S):
    arcs=[]
    for v in set(S):
        inv=1.0/(14*v)
        for k in range(v+1):
            lo=max((14*k-1)*inv,0.0); hi=min((14*k+1)*inv,1.0)
            if lo<hi: arcs.append((lo,hi))
    arcs.sort(); tot=0.0; cl,ch=arcs[0]
    for lo,hi in arcs[1:]:
        if lo<=ch: ch=hi if hi>ch else ch
        else: tot+=ch-cl; cl,ch=lo,hi
    tot+=ch-cl
    return 1.0-tot

# decoupling limit for cores = {1..13}\{e}
print("DECOUPLING FLOOR: core = {1..13}\\{e}, add w->inf.  L -> (6/7)*meas(G_e)")
print(f"  {'e':>3s} {'meas(G_e)':>12s} {'(6/7)meas(G_e)=floor':>22s} {'>=1/1260?':>10s}")
base=list(range(1,14))
floors=[]
for e in base:
    core=[x for x in base if x!=e]
    mg=L_exact(core)  # = meas(G_e)
    floor=F(6,7)*mg
    floors.append(floor)
    print(f"  {e:>3d} {str(mg):>12s} {str(floor):>22s} {'yes' if floor>=F(1,1260) else 'NO':>10s}")
print(f"  MIN decoupling floor over e = {min(floors)} = {float(min(floors)):.6e}")
print(f"  (this is the smallest L achievable by sending ONE element to infinity)")
print()

# verify the limit numerically for e=12 (smallest core gap)
print("Numerical check e=12: L({1..11,13,w}) for large w should approach (6/7)*426/35035:")
core=[x for x in base if x!=12]
lim=float(F(6,7)*F(426,35035))
for w in [500,1000,5000,20000,100000]:
    print(f"   w={w:7d}: L={L_float(core+[w]):.8e}   target {lim:.8e}")
print()

# two-element growth: 11 fixed + {w,w'} large, can L dip below 1/1260?
print("TWO-ELEMENT growth: drop {12,13}, add {w, w2} large; smallest L over scan:")
import random
random.seed(7)
best=(2.0,None)
fixed=[x for x in base if x not in (12,13)]
for _ in range(200000):
    w=random.randint(14,2000); w2=random.randint(14,2000)
    if w==w2 or w in fixed or w2 in fixed: continue
    lf=L_float(fixed+[w,w2])
    if lf<best[0]: best=(lf,(w,w2))
print(f"   best L={best[0]:.6e} at {best[1]}  (vs 1/1260={1/1260:.6e})")
if best[0]<1/1260*1.2:
    S=fixed+list(best[1]); print(f"   exact L = {L_exact(S)}")
print()
print("Drop {6,12}, add {w,w2} (6 and 12 are the two smallest-gap drops):")
best2=(2.0,None); fixed2=[x for x in base if x not in (6,12)]
for _ in range(200000):
    w=random.randint(14,2000); w2=random.randint(14,2000)
    if w==w2 or w in fixed2 or w2 in fixed2: continue
    lf=L_float(fixed2+[w,w2])
    if lf<best2[0]: best2=(lf,(w,w2))
print(f"   best L={best2[0]:.6e} at {best2[1]}")
if best2[0]<1/1260*1.2:
    print(f"   exact L = {L_exact(fixed2+list(best2[1]))}")
