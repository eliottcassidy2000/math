from fractions import Fraction as F
from math import gcd
def Mexact(V):
    mx=max(V); Qs=set()
    for i in range(len(V)):
        for j in range(i,len(V)):
            s=V[i]+V[j]
            for d in range(2,2*mx+1):
                if s%d==0: Qs.add(d)
    best=F(0); bq=None
    for q in sorted(Qs):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            m=q
            for v in V:
                r=(v*a)%q; d=r if r<q-r else q-r
                if d<m: m=d
            if F(m,q)>best or (F(m,q)==best and bq and m<bq[0]):
                best=F(m,q); bq=(m,q,a)
    return best,bq
def show(name,V):
    M,(val,q,a)=Mexact(V)
    R=sorted((v*a)%q for v in V); gaps=[R[i+1]-R[i] for i in range(len(R)-1)]
    cov14 = any(v%14==0 for v in V)
    print("%-30s M=%s(%.5f) val=%d q=%d a=%d cover14?%s"%(name,M,float(M),val,q,a,cov14))
    print("   residues:",R); print("   gaps:",gaps," #gaps<val=",sum(1 for g in gaps if g<val))
show("{1..11,13,24}",[1,2,3,4,5,6,7,8,9,10,11,13,24])
show("deepwell {1..12,182}",list(range(1,13))+[182])
# hunt: val>=2, non-covering, non-AP core, M<1/13 ?
print("--- hunt val>=2 non-cover14 non-AP-core M<1/13 (small) ---")
found=0
import itertools
base=list(range(1,13))
for w in range(13,120):
    if w%14==0: continue  # skip covering-of-14
    for drop in range(1,13):
        C=[x for x in base if x!=drop]  # non-AP core (missing one of 1..12)
        V=sorted(C+[w])
        if len(set(V))<13: continue
        M,(val,q,a)=Mexact(V)
        if M<F(1,13) and val>=2:
            print("  FOUND",V,"M=%s val=%d q=%d"%(M,val,q)); found+=1
            if found>=5: break
    if found>=5: break
if not found: print("  none found (val>=2, non-cover14, one-hole core, w<120)")
