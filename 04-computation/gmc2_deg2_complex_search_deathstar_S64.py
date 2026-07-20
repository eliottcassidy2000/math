from fractions import Fraction as Fr
from math import factorial
from itertools import product
# complex arithmetic via Gaussian rationals: represent as (re,im) Fraction pairs
class G:
    __slots__=('r','i')
    def __init__(s,r=0,i=0): s.r=Fr(r); s.i=Fr(i)
    def __add__(a,b): return G(a.r+b.r,a.i+b.i)
    def __sub__(a,b): return G(a.r-b.r,a.i-b.i)
    def __mul__(a,b): return G(a.r*b.r-a.i*b.i, a.r*b.i+a.i*b.r)
    def __eq__(a,b): return a.r==b.r and a.i==b.i
    def __repr__(s): return f"{s.r}{'+' if s.i>=0 else ''}{s.i}i"
    def iszero(s): return s.r==0 and s.i==0
def pmul(p,q):
    d={}
    for i,u in enumerate(p):
        for j,v in enumerate(q):
            k=i+j; d[k]=d.get(k,G())+u*v
    n=max(d)+1 if d else 1; o=[G()]*n
    for k,v in d.items(): o[k]=v
    return o
def padd(*ps):
    n=max(len(p) for p in ps); o=[G() for _ in range(n)]
    for p in ps:
        for i,u in enumerate(p): o[i]=o[i]+u
    return o
def Er(p): 
    s=G()
    for j,c in enumerate(p): s=s+G(c.r*factorial(j),c.i*factorial(j))
    return s
def EPm(a,b,c,M):
    ac=pmul(a,c); alpha=[G()]+ac; res=[]
    for m in range(1,M+1):
        tot=[G()]
        for k in range(m//2+1):
            coef=Fr(factorial(m),factorial(k)**2*factorial(m-2*k))
            term=[G(1)]
            for _ in range(k): term=pmul(term,alpha)
            for _ in range(m-2*k): term=pmul(term,b)
            tot=padd(tot,[G(coef*u.r,coef*u.i) for u in term])
        res.append(Er(tot))
    return res
def charges_nonzero(a,b,c):  # two-sided if a!=0 and c!=0, OR b!=0 with... one-sided = b==0 & ac==0
    az=all(x.iszero() for x in a); cz=all(x.iszero() for x in c); bz=all(x.iszero() for x in b)
    ac=pmul(a,c); acz=all(x.iszero() for x in ac)
    onesided = bz and acz
    return not onesided

# SEARCH complex deg-2 b, a=c=1: b0,b1,b2 in Gaussian integers {0,±1,±i,±1±i,±2,±2i}
vals=[G(0),G(1),G(-1),G(0,1),G(0,-1),G(1,1),G(-1,-1),G(2),G(0,2)]
M0=8; found=0; checked=0
for co in product(vals,repeat=3):
    b=list(co)
    if all(x.iszero() for x in b): continue
    checked+=1
    E=EPm([G(1)],b,[G(1)],M0)
    if all(x.iszero() for x in E):
        if charges_nonzero([G(1)],b,[G(1)]):
            found+=1
            if found<=5: print("  COMPLEX deg-2 nullcone (a=c=1):", [str(x) for x in b])
print(f"deg-2 complex (a=c=1, b Gaussian-int box): checked {checked}, NON-one-sided nullcone members: {found}")

# SEARCH deg-3 b, a=c=1, smaller box
vals3=[G(0),G(1),G(-1),G(0,1),G(0,-1),G(1,1)]
found3=0; checked3=0
for co in product(vals3,repeat=4):
    b=list(co)
    if all(x.iszero() for x in b): continue
    checked3+=1
    E=EPm([G(1)],b,[G(1)],M0)
    if all(x.iszero() for x in E):
        if charges_nonzero([G(1)],b,[G(1)]): found3+=1
print(f"deg-3 complex (a=c=1, b box): checked {checked3}, NON-one-sided nullcone members: {found3}")
print("=> NC2 holds (no counterexample) for complex deg-2 and deg-3." if found+found3==0 else "!!! COUNTEREXAMPLE")
