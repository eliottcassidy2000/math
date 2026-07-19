"""Does a covering non-AP family exist with SMALL margin (M just above 1/13) AND a FAR-from-AP core
(Hamming distance >6 from every dilated AP)? If NO (small margin => near-AP), the favorable shape
holds: far cores have comfortable margin and a crude bound suffices; near cores are done by
THM-1004/5/6. death-star-S58h."""
from fractions import Fraction as F
from math import gcd
import random, itertools
def Mexact(V):
    mx=max(V);Qs=set()
    for i in range(len(V)):
        for j in range(i,len(V)):
            s=V[i]+V[j]
            for d in range(2,2*mx+1):
                if s%d==0:Qs.add(d)
    best=F(0)
    for q in Qs:
        for a in range(1,q):
            if gcd(a,q)!=1:continue
            m=min(min((v*a)%q,q-((v*a)%q)) for v in V)
            if F(m,q)>best:best=F(m,q)
    return best
def covers(V,hi=14): return all(any(v%k==0 for v in V) for k in range(2,hi+1))
def ham_core(C):
    best=12
    for d in range(1,max(C)//12+2):
        ap=set(d*i for i in range(1,13)); best=min(best,12-len(set(C)&ap))
    return best
random.seed(11)
THR=F(1,13); tested=0
smallfar=[]  # covering non-AP, margin<0.01, Hamming>6
minmargin_far=(None,99)  # smallest margin among far cores
for trial in range(6000):
    # build a 13-family: core mostly-random, plus a far element to help covering
    base=sorted(random.sample(range(1,40),12))
    w=random.choice([84,156,182,168,364,120,110,132,140,273,productf] if False else [84,156,182,168,120,110,132,140,364,78,132])
    V=sorted(set(base+[w]))
    if len(V)!=13: continue
    if not covers(V,13): continue   # M<1/13 requires covers 2..13
    tested+=1
    M=Mexact(V)
    if M<=THR: 
        print("!! covering M<=1/13 non-deepwell:",V,"M=",M,flush=True); continue
    marg=float(M-THR)
    C=V[:]; C.remove(max(V)); h=ham_core(C)
    if h>6 and marg<minmargin_far[1]: minmargin_far=(V,marg)
    if marg<0.01 and h>6:
        smallfar.append((V,float(M),marg,h))
        if len(smallfar)<=8: print("SMALL-MARGIN FAR CORE:",V,"M=%.5f margin=%+.5f H=%d"%(float(M),marg,h),flush=True)
print("tested=%d covering families; small-margin(<0.01) FAR(H>6) cores found=%d"%(tested,len(smallfar)),flush=True)
print("min margin among FAR cores: margin=%.5f at %s"%(minmargin_far[1] if minmargin_far[0] else -1, minmargin_far[0]),flush=True)
