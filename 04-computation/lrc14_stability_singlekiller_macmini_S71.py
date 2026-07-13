#!/usr/bin/env python3
"""mac-mini-S71: toward HYP-2566 (closed-form covering-min). Attack the SINGLE-KILLER case.
The one gap in THM-724 is near-tight NON-dilated cores with large s. Route: quantitative
LRC(13) STABILITY -- if M_core is near 1/13, the core is near a dilated AP c*{1..12}, where the
shallow witness (base 13c) applies (STABLY). Test both legs:
(A) LRC13 STABILITY: does M_core <= 1/13 + eta force the 12-core near a dilated AP?
(B) SHALLOW-WITNESS STABILITY: for a PERTURBED dilated core + covering killer, does the
    base-13c witness still give M(S) >= 14/183?
If both hold, single-killer closes in closed form modulo the stability constants."""
from fractions import Fraction as F
from math import gcd
target=F(14,183); onethird=F(1,13)
def resd(x,q): r=x%q; return min(r,q-r)
def M_exact(S,Qmax):
    best=F(0); arg=None
    for q in range(2,Qmax+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            mind=q
            for v in S:
                d=resd(a*v,q)
                if d<mind: mind=d
                if d==0: break
            if mind>0 and F(mind,q)>best: best=F(mind,q); arg=(q,a)
    return best,arg
def M_core(C,Qmax=250):
    best=F(0); arg=None
    for q in range(2,Qmax+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            mind=min(resd(a*v,q) for v in C)
            if mind>0 and F(mind,q)>best: best=F(mind,q); arg=(q,a)
    return best,arg
def dist_to_dilatedAP(C):
    """min over c of |C - c*{1..12}| (as multiset L1 after best alignment); crude: try c=round(min)"""
    C=sorted(C); best=(10**9,None)
    for c in range(1, C[0]+3):
        ap=[c*k for k in range(1,13)]
        d=sum(abs(x-y) for x,y in zip(C,ap))
        if d<best[0]: best=(d,c)
    return best

print(f"14/183={float(target):.6f}  1/13={float(onethird):.6f}\n")
print("(A) LRC13 STABILITY: near-tight 12-cores (M_core close to 1/13) vs distance to dilated AP")
print("    scan 12-subsets of small integers, bucket by M_core, report min L1-distance to c*{1..12}:")
import itertools, random
random.seed(5)
buckets={}
for _ in range(6000):
    C=sorted(random.sample(range(1,40),12))
    if gcd(*C[:2])==0: pass
    Mc,_=M_core(C,120)
    key = "tight[1/13,1/13+.001]" if onethird<=Mc<=onethird+F(1,1000) else \
          "near(1/13+.001,1/12]" if onethird<Mc<=F(1,12) else \
          ">1/12" if Mc>F(1,12) else "<1/13"
    d,c=dist_to_dilatedAP(C)
    buckets.setdefault(key,[]).append(d)
for k in ["tight[1/13,1/13+.001]","near(1/13+.001,1/12]",">1/12","<1/13"]:
    v=buckets.get(k,[])
    if v: print(f"    {k:26s}: {len(v)} cores, L1-dist to dilated-AP  min={min(v)} median={sorted(v)[len(v)//2]}")

print("\n(B) SHALLOW-WITNESS STABILITY: perturbed dilated core c*{1..12}+delta + covering killer.")
print("    Does M(S) stay >= 14/183 (via base-13c-type witness)?")
belowB=[]; testedB=0
for c in range(2,7):
    base=[c*k for k in range(1,13)]
    for j in range(12):
        for delta in [-2,-1,1,2]:
            core=sorted(base[:j]+[base[j]+delta]+base[j+1:])
            if len(set(core))!=12 or core[0]<1: continue
            for killer in [182,364,c*13*14, 2*c*13*7]:
                S=sorted(set(core+[killer]))
                if len(S)!=13: continue
                g=0
                for x in S: g=gcd(g,x)
                if g!=1: continue
                if not all(any(v%q==0 for v in S) for q in range(2,15)): continue
                testedB+=1
                M,arg=M_exact(S,min(2*max(S),400))
                if M<target: belowB.append((M,S,arg))
print(f"    tested {testedB} perturbed-dilated-core primitive covering single-killer configs.")
print(f"    below 14/183: {'NONE' if not belowB else belowB[:6]}")
if not belowB:
    print("    => shallow-witness is STABLE: perturbed dilated cores stay >= 14/183.")
    print("    Combined with kps-S127 (primitivity forces s=1 at exact tight) + balance-slack")
    print("    (M_core>1/13+eps), the SINGLE-KILLER case closes modulo the stability window.")
