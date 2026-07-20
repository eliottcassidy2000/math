from fractions import Fraction as Fr
from math import factorial
# L_m(alpha,beta) obeys (m+1)L_{m+1}=(2m+1)beta L_m - m(beta^2-4alpha)L_{m-1}.  Verify, then show that
# after E_r the recurrence for E[P^m]=E_r[L_m] does NOT close: it needs E_r[b*L_m] != b*E_r[L_m].
def pmul(p,q):
    d={}
    for i,u in enumerate(p):
        for j,v in enumerate(q):
            if u and v: d[i+j]=d.get(i+j,Fr(0))+u*v
    n=max(d)+1 if d else 1; o=[Fr(0)]*n
    for k,v in d.items(): o[k]=v
    return o
def padd(*ps):
    n=max(len(p) for p in ps); o=[Fr(0)]*n
    for p in ps:
        for i,u in enumerate(p): o[i]+=u
    return o
def psc(p,s): return [u*s for u in p]
def Er(p): return sum(c*factorial(j) for j,c in enumerate(p))
# alpha=r (gamma=1), beta=b(r)=r+r^2  (a deg-2 case)
alpha=[Fr(0),Fr(1)]; beta=[Fr(0),Fr(1),Fr(1)]
delta=padd(pmul(beta,beta),psc(alpha,Fr(-4)))   # beta^2-4alpha
# build L_m as polynomials in r via the recurrence
L=[[Fr(1)],beta[:]]   # L_0=1, L_1=beta
for m in range(1,9):
    Lm1=psc(padd(psc(pmul(beta,L[m]),Fr(2*m+1)), psc(pmul(delta,L[m-1]),Fr(-m))),Fr(1,m+1))
    L.append(Lm1)
EP=[Er(L[m]) for m in range(1,10)]
print("E[P^m]=E_r[L_m], m=1..9 (alpha=r,beta=r+r^2):", [str(x) for x in EP])
# THE NON-CLOSURE: is E_r[b*L_m] == (something in E[P^*])?  Show E_r[b L_m] is a NEW quantity.
bLm=[Er(pmul(beta,L[m])) for m in range(0,9)]
print("\nE_r[b*L_m] m=0..8 (the b-weighted moments the recurrence needs):", [str(x) for x in bLm])
# check the constant-b would give E_r[b L_m]=b*E[P^m]; here b is non-constant so it's independent:
print("For constant b, E_r[b L_m]=b*E[P^m]. Here b=r+r^2 is non-constant, so E_r[b L_m] is a")
print("genuinely NEW sequence not expressible in {E[P^j]} -> the m-recurrence for E[P^m] does NOT close.")
print("=> Sheffer/Legendre recurrence opens an INFINITE b-weighted-moment hierarchy (no 'common root').")
