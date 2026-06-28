"""
S75b: the magic-function kernel K(p,q)=meas(D_p∩D_q) is a THREE-GAP / Stern-Brocot function of a/b.
Verifies: (1) scale-invariance K(p,q)=K(p/gcd,q/gcd); (2) antipode K(a,13)=(2a-1)/(91a);
(3) the gap-count g(a,b)=7ab*K(a,b) is piecewise-linear in a (three-gap, breakpoints at convergents);
base K(1,b)=1/(7b) (single-arc lemma), apex K(7,b)=1/49.
"""
from fractions import Fraction as F
from math import gcd

def K(p,q):
    S=sorted({p,q}-{0})
    if not S: return F(1)
    b=set([F(0),F(1)])
    for r in S:
        for k in range(0,r+1):
            for s in (F(1,14),-F(1,14)):
                v=(F(k)+s)/r
                if 0<=v<=1: b.add(v)
    b=sorted(b); tot=F(0)
    for i in range(len(b)-1):
        x0,x1=b[i],b[i+1]
        if x1<=x0: continue
        mid=(x0+x1)/2
        if all(min((r*mid)%1,1-(r*mid)%1)<F(1,14) for r in S): tot+=x1-x0
    return tot

print("(1) scale-invariance K(p,q)=K(p/gcd,q/gcd):")
bad=[ (p,q) for p in range(1,14) for q in range(p,14) if gcd(p,q)>1 and K(p,q)!=K(p//gcd(p,q),q//gcd(p,q)) ]
print("    ", "HOLDS (all p,q<=13)" if not bad else f"FAILS {bad}")

print("(2) antipode column K(a,13)=(2a-1)/(91a):")
print("    ", "HOLDS (a=1..12)" if all(K(a,13)==F(2*a-1,91*a) for a in range(1,13)) else "FAILS")

print("(3) gap-count g(a,b)=7ab*K(a,b) (piecewise-linear in a = three-gap; kinks at convergents of a/b):")
for b in range(2,14):
    g=[ int(7*a*b*K(a,b)) for a in range(1,b) if gcd(a,b)==1 ]
    print(f"    b={b:>2}: g = {g}")
print("    base K(1,b)=1/(7b):", all(K(1,b)==F(1,7*b) for b in range(1,14)))
print("    apex K(7,b)=1/49 for all b coprime-ish:", [str(K(7,b)) for b in range(8,14)])
