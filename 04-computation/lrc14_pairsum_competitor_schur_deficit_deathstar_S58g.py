"""The pair-sum competitor q'=v_i+v_j (THM-724 edge pair) and its margin over 1/13 for covering
non-AP cores (death-star-S58g). Finding: margin M-1/13 tracks the core's Schur deficit 66-T(core)
(THM-730 quantitative inverse, localized to the competitor). Smallest margin at deficit 1 (near-AP),
which is PROVEN M>=2/25>1/13 by THM-1004/5/6. Uniform bound (margin>0 for all deficit>=1) = the
Freiman/E3 crux, OPEN (THM-1028 residual = near-tight fragmented far cores)."""
from fractions import Fraction as F
from math import gcd
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
def schur(C): S=set(C); return sum(1 for a in C for b in C if a+b in S)
for name,C,far in [("{1..11,13}",list(range(1,12))+[13],156),("{1..11,14}",list(range(1,12))+[14],156),
                   ("{1..10,12,13}",list(range(1,11))+[12,13],110),("{1..10,13,110}",list(range(1,11))+[13,110],156),
                   ("AP{1..12}",list(range(1,13)),182)]:
    V=sorted(C+[far]); M=Mexact(V)
    print("%-16s T=%d deficit=%d M=%s margin=%+.5f"%(name,schur(C),66-schur(C),M,float(M-F(1,13))))
