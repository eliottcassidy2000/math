# death-star-2026-07-17-S56 (THM-996 Part III): AP vs GW census-twin / spectrum-split.
# liveCount twins at threshold; full depth histogram splits below it. + general resonance-confinement.
from math import gcd
from collections import Counter
from fractions import Fraction as F
def minnum(v,q,p): r=(v*p)%q; return min(r,q-r)
def histogram(fam,q):
    h=Counter()
    for p in range(1,q): h[min(minnum(v,q,p) for v in fam)]+=1
    return h
def live_count(fam,n,q): return sum(1 for p in range(1,q) if all(q<=n*((v*p)%q)<=(n-1)*q for v in fam))
def M_exact(fam,Qmax):
    best=F(0); arg=None
    for q in range(2,Qmax+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            m=min(minnum(v,q,a) for v in fam); val=F(m,q)
            if val>best: best=val; arg=(a,q)
    return best,arg
if __name__=="__main__":
    n=14; AP=list(range(1,14)); GW=[1,2,3,4,5,6,7,8,9,10,11,13,24]; DW=list(range(1,13))+[182]
    for nm,fam,Qm in [("AP",AP,60),("GW",GW,120)]:
        print(nm,"M=",M_exact(fam,Qm))
    print("Part I: nonresonant live for AP/GW to q=560:",
          [q for q in range(2,561) if q%14 and (live_count(AP,n,q) or live_count(GW,n,q))])
    print("Deep-well off-resonance live@183:",live_count(DW,n,183))
    for q in [14,28,42,70,98,112,196]:
        print(f"q={q}: liveAP={live_count(AP,n,q)} liveGW={live_count(GW,n,q)} hist_same={histogram(AP,q)==histogram(GW,q)}")
