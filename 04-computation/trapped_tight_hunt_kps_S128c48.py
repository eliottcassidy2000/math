#!/usr/bin/env python3
"""trapped_tight_hunt_kps_S128c48.py -- ADVERSARIAL: try to REFUTE THM-995's conjecture by
finding a TRAPPED family that is TIGHT (M = 1/14).  Two prongs:
 (A) the dilation mechanism: c*V dilates of any tight V fail non-clusterable (common factor c);
     confirm, and confirm primitive tight reps fail some other hypothesis.
 (B) direct hunt: perturb/search near the tight configs for a trapped M<=1/14+eps family.
"""
import sys, random
from fractions import Fraction as F
from math import gcd
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)
random.seed(999)
def nd(x):
    fx=x-int(x)
    if fx<0: fx+=1
    return min(fx,1-fx)
def M_exact(V):
    cand=set()
    n=len(V)
    for i in range(n):
        vi=abs(V[i])
        for a in range(1,2*vi): cand.add(F(a,2*vi))
    for i in range(n):
        for j in range(i+1,n):
            for s in (abs(V[i]+V[j]),abs(V[i]-V[j])):
                if s==0: continue
                for a in range(1,s): cand.add(F(a,s))
    best=F(0)
    for t in cand:
        if not (0<t<1): continue
        m=min(nd(v*t) for v in V)
        if m>best: best=m
    return best
def gap(V): A=[abs(v) for v in V]; return any(a>13*b for a in A for b in A)
def compressed(V): A=[abs(v) for v in V]; return all(any(j!=i and A[i]<=13*A[j] for j in range(len(A))) for i in range(len(A)))
def distinct(V): A=[abs(v) for v in V]; return len(set(A))==len(A)
def maxge23(V): return max(abs(v) for v in V)>=23
def noncluster(V):
    n=len(V)
    for i0 in range(n):
        rest=[abs(V[j]) for j in range(n) if j!=i0]
        g=reduce(gcd,rest,0)
        if g>=2 and V[i0]%g!=0: return False
    return True
def covering(V):
    for q in range(2,15):
        if min(nd(v*F(1,q)) for v in V)>=F(1,14): return False
    return True
def comm_res(V):
    for d in range(2,max(abs(v) for v in V)+1):
        for a in range(d):
            if all((v-a)%d==0 for v in V): return False
    return True
HYPS=[("cov",covering),("gap",gap),("comp",compressed),("dist",distinct),("m23",maxge23),("ncl",noncluster),("ncr",comm_res)]
def fails(V): return [nm for nm,fn in HYPS if not fn(V)]

print("== (A) dilation mechanism: dilates of tight families fail non-clusterable ==")
tights=[("AP{1..13}",list(range(1,14))),("sporadic{1..11,13,24}",list(range(1,12))+[13,24])]
for name,V in tights:
    for c in (1,2,3,5):
        cV=[c*v for v in V]
        f=fails(cV)
        m=M_exact(cV)
        tag = "TIGHT" if m==F(1,14) else "M=%s"%m
        print("  %-22s x%d: %s  fails=%s"%(name,c,tag,f))
print()
print("== (B) direct adversarial hunt: any TRAPPED family with M <= 1/14 + 1e-9 ? ==")
def is_trapped(V): return len(fails(V))==0
found=None; scanned=0
# hunt near tight structure: base = sporadic-like, perturb one element to seek covering+trapped
bases=[list(range(1,12))+[13,24], list(range(1,13))+[24], list(range(1,12))+[14,27]]
for base in bases:
    for trial in range(3000):
        V=base[:]
        # perturb 1-3 coords
        for _ in range(random.randint(1,3)):
            i=random.randrange(13); V[i]=V[i]+random.choice([-2,-1,1,2,3])
        if len(set(abs(x) for x in V))<13 or any(x<=0 for x in V): continue
        scanned+=1
        if is_trapped(V):
            m=M_exact(V)
            if m <= F(1,14)+F(1,10**6):
                found=(V,m); break
    if found: break
if found:
    print("  *** REFUTED: trapped tight family found: %s M=%s ***"%(list(found[0]),found[1]))
else:
    print("  no trapped family with M near 1/14 in %d perturbations of the tight bases"%scanned)
    print("  -> THM-995 conjecture SURVIVES the adversarial hunt (near-tight perturbations all escape the cut)")
print("DONE")
