from fractions import Fraction as F
from math import gcd, lcm
import numpy as np
def nrm(x): f=x-(x.numerator//x.denominator); return min(f,1-f)
def M_exact(S):
    C=set()
    for i in range(len(S)):
        for j in range(len(S)):
            for d in (S[i]+S[j],abs(S[i]-S[j])):
                if d:
                    for m in range(d+1): C.add(F(m,d))
        for m in range(2*S[i]+1): C.add(F(m,2*S[i]))
    best=F(0)
    for t in C:
        if 0<t<1:
            v=min(nrm(s*t) for s in S)
            if v>best: best=v
    return best
def M_grid(S,Q):
    t=np.arange(1,Q)/Q; f=np.ones(Q-1)
    for v in S: f=np.minimum(f,np.abs(((v*t+0.5)%1)-0.5))
    return f.max()
def is_cov(S): return all(any(s%q==0 for s in S) for q in range(2,15))
for nm,S in [("{1..11,13,84}",[1,2,3,4,5,6,7,8,9,10,11,13,84]),("{1..12,182}",list(range(1,13))+[182]),
             ("{1..10,12,13,?} drop11+force",None)]:
    if S is None: continue
    Me=M_exact(S); Mg=M_grid(S,600000)
    print(f"{nm:>16}: M_exact={Me}={float(Me):.6f}  grid={Mg:.6f}  covering={is_cov(S)}")
print(f"7/89={float(F(7,89)):.6f}, 14/183={float(F(14,183)):.6f}; 14/183<7/89: {F(14,183)<F(7,89)}")
# focused: drop k (a binding unit), add lcm(k,14) -- the tightest displacements
print("\ndrop a UNIT k (binding), add killer lcm(k,14): M and displacement (units bind tightest?):")
ap=list(range(1,14)); units=[1,3,5,9,11,13]
res=[]
for k in units:
    S=[x for x in ap if x!=k]+[lcm(k,14)]
    if len(set(S))==13 and is_cov(S):
        M=M_exact(S); res.append((M,k,lcm(k,14)))
for M,k,w in sorted(res):
    print(f"   drop unit {k:>2}, add {w:>3}: M={M}={float(M):.6f}  delta={float(M-F(1,14)):.6f}")
