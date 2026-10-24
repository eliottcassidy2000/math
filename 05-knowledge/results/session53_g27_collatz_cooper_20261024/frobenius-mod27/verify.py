import json
from math import comb
M=80

def require(ok,msg):
    if not ok:raise RuntimeError(msg)
def add(a,b):return [x+y for x,y in zip(a,b)]
def scale(a,k):return [k*x for x in a]
def mul(a,b):return [sum(a[j]*b[n-j] for j in range(n+1)) for n in range(len(a))]
def inv(a):
    require(a[0]==1,'nonunit constant');r=[1]
    for n in range(1,len(a)):r.append(-sum(a[j]*r[n-j] for j in range(1,n+1)))
    return r
def divint(a,k):
    require(all(x%k==0 for x in a),'inexact integer division');return [x//k for x in a]
def chi(d):return 0 if d%7==0 else (1 if d%7 in (1,2,4) else -1)
raw=[1]+[2*sum(chi(d) for d in range(1,n+1) if n%d==0) for n in range(1,3*M+1)]
raw2=mul(raw,raw);raw3=mul(raw2,raw)
th=raw[:M+1];ith=inv(th);one=[1]+[0]*M
thp=[th[n//3] if n%3==0 else 0 for n in range(M+1)]
U=lambda a:[a[3*n] for n in range(M+1)]
P1=scale(mul(U(raw),ith),3)
P2=scale(mul(U(raw2),mul(ith,ith)),3)
P3=scale(mul(U(raw3),mul(mul(ith,ith),ith)),3)
e1=P1;e2=divint(add(mul(P1,P1),scale(P2,-1)),2)
e3=divint(add(add(mul(mul(P1,P1),P1),scale(mul(P1,P2),-3)),scale(P3,2)),6)
z=mul(thp,ith)
actual=[scale(one,27),add(scale(e1,9),scale(z,-27)),add(scale(e2,3),scale(mul(z,e1),-9)),add(e3,scale(mul(z,e2),-3)),scale(mul(z,e3),-1)]
# Independent eta expressions for t and s via h.
eta=[1]+[0]*M
for n in range(1,M+1):
    for j in range(M,n-1,-1):eta[j]-=eta[j-n]
eta7=[eta[n//7] if n%7==0 else 0 for n in range(M+1)]
a=mul(eta7,inv(eta));a2=mul(a,a);a4=mul(a2,a2);h=[0]+a4[:-1];h2=mul(h,h)
den=add(add(one,scale(h,13)),scale(h2,49));iden=inv(den)
t=mul(h,iden);s=mul(add(one,scale(h2,-49)),iden)
expected=[scale(one,27),scale(one,0),scale(one,-18),scale(s,-8),add(scale(t,24),scale(one,-1))]
require(actual==expected,'T3 exact q-series identity')
# Excellent-lift congruence in q coordinates, avoiding formal inversion.
F=mul(th,th);Fp=mul(thp,thp)
trunc=add(add(one,scale(t,4)),scale(mul(t,t),48))
require(all(x%27==0 for x in add(F,scale(mul(trunc,Fp),-1))),'F excellent congruence mod27')
# Symbolic b,t polynomial identity, with exponent pairs (b,t).
def pa(a,b):
    r=dict(a)
    for m,c in b.items():r[m]=r.get(m,0)+c
    return {m:c for m,c in r.items() if c}
def ps(a,k):return {m:k*c for m,c in a.items()}
def pm(a,b):
    r={}
    for (i,j),c in a.items():
        for (k,l),d in b.items():r[i+k,j+l]=r.get((i+k,j+l),0)+c*d
    return r
Q={(0,0):1,(0,1):-26,(0,2):-27};X={(0,0):1,(0,1):6,(1,0):9}
expr=pa(pm(Q,pa(pm({(0,0):1,(0,1):-24},pm(X,X)),ps(X,8))),{(0,0):18})
rem={m:c%27 for m,c in expr.items() if c%27}
require(rem=={(1,1):9,(1,0):9,(0,1):18,(0,3):18},'symbolic b residual')
s7=lambda n:sum(comb(n,k)**2*comb(n+k,k)*comb(2*k,n) for k in range(n+1) if 2*k>=n)
sharp=s7(3)-s7(1)
require(sharp==756 and sharp%27==0 and sharp%81!=0,'sharpness')
print(json.dumps({'sharpness_t3_coefficient':sharp,'status':'PASS','q_bound':M,'exact_T3_coefficient_checks':5,'symbolic_residual':'9*((1+t)*b-t-t^3) modulo27','F_over_Fsigma_mod27':'1+4*t+21*t^2','proof_dependency':'inherited modularity and Fricke degree bounds, plus eta/theta parameter identities'},indent=2))
