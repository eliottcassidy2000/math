from functools import reduce
from fractions import Fraction
from itertools import product
from math import gcd
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent /
                      "session328_fermat_gluing_and_phase_repair_20261114"))
from fermat_profiles import strong, canon, extend, hp_all, ends, rooted, glue, profile, submasks

def require(x):
    if not x: raise RuntimeError('S332 replay failed')
def conv(A):
    m=(1<<len(A))-1
    return tuple(m^x^(1<<i) for i,x in enumerate(A))
def induced(A,verts):
    return tuple(sum(((A[i]>>j)&1)<<k for k,j in enumerate(verts)) for i in verts)
def state(A,v):
    F=hp_all(A);m=len(F)-1;d=m^(1<<v);S=ends(conv(A))
    L=F[d]+sum(F[d^B]*sum(S[B][i] for i in range(len(A)) if (B&A[v])>>i&1) for B in submasks(d) if B)
    return F[m],F[d],L,ends(A)[m][v]
def seam(A,v,key,j=1):
    n=len(A); B=glue(A,v,j)
    r={'O':v,'I':n,'T':n+j}.get(key,key)
    return B,r
def word(A,v,seq):
    for a in seq:A,v=seam(A,v,a)
    return A,v
def step(Z,j,k):
    H,h,L,E=Z;d=1<<j;F=H+d*L
    if k==0:return F,h,2*d*L-h,h
    if k==j+1:return F,H,F,H
    t=j-k
    return F,H+(d//2)*L if j>=2 else H+L-E,(2<<t)*H+(1<<k)*3**t*L-(E if k==1 else 0),(1<<(k-1))*L
def decode(A):
    degree=[r.bit_count() for r in A]
    high=[i for i,d in enumerate(degree) if d>=6]
    if 2 in degree:
        mode="OIT"
        cand=[i for i,d in enumerate(degree) if d==5 and any(A[i]>>j&1 for j in high)]
        require(len(cand)==1)
        marked=cand[0]; verts=sorted(high+[marked])
    else:
        mode="IOT";verts=high
        cand=[i for i in high if degree[i]==6]
        require(len(cand)==1);marked=cand[0]
    return induced(A,verts),verts.index(marked),mode

B=(2,4,1);v=2;L=[('OIT',109),('IOT',109)]
for m,H,e in [(1,109,31),(2,4003,1137)]:
    q=[]
    for mword in product(('OIT','IOT'),repeat=m):
        A,r=B,v
        for w in mword:A,r=word(A,r,w)
        require(strong(A) and state(A,r)==(H,e,H,e) and len(A)==3+6*m)
        q.append(profile(A)[1])
    require(q==([4580,4580] if m==1 else [1671956,1698644,1631444,1671956]))
for mword in product(('OIT','IOT'),repeat=5):
    A,r=B,v
    for w in mword:A,r=word(A,r,w)
    for w in reversed(mword):
        A,r,decoded=decode(A)
        require(decoded==w)
    require((A,r)==(B,v))
print('literal H,e through depth two and all 32 inverse degree codes at depth five')

def ctx(A):
    C=induced(A,list(range(1,len(A))))
    _,_,q,p,m,_=list(rooted(conv(A)))[0]
    return profile(C)[1],q,p,m
require(ctx(word((0,),0,'OIT')[0])==(128,596,395,87))
require(ctx(word((0,),0,'IOT')[0])==(256,644,395,183))
def bu(A,u):
    F=hp_all(A); m=len(F)-1
    return sum(F[m^S]*p[u] for S,p in enumerate(ends(A)))
for Z in ('OIT','IOT'):
    U,r=word(B,v,Z)
    A,t=seam(U,r,'T',1)
    p=profile(word(A,t,'OIT')[0])[1]
    q=profile(word(A,t,'IOT')[0])[1]
    require(p-q==48*(4*bu(U,r)-2*profile(U)[1])<0)
print('universal observer-side context constants verified and strict signs checked')

T={(0,)};cases=0
for n in range(1,6):
    if n>1:T={canon(extend(A,c)) for A in T for c in range(1<<(n-1))}
    for A in T:
        if not strong(A):continue
        for v in range(n):
            X=state(A,v); H,h,L,E=X
            for j in (1,2,4):
                U=glue(A,v,j)
                for k in range(j+2):
                    Y=state(U,n+k-1 if k else v)
                    require(Y==step(X,j,k));a,b,c,d=Y
                    g=gcd(H,L,h) if k==0 else gcd(H,L,E) if k==1 else gcd(H,L)
                    require(gcd(*Y)==gcd(a,c,b)==gcd(a,c,d)==g);cases+=1
print(f'exact balanced valuation law checked on {cases} literal root transitions')

def twoport(X,j,k):
    return step(X,j,0 if k=='O' else 1)
anchor=(Fraction(1),Fraction(2,3),Fraction(1),Fraction(1,2))
W=Fraction(0);by={}
for j1,j2,j3 in product((1,2,4),repeat=3):
    for p1,p2 in product('OI',repeat=2):
        f=twoport(twoport(twoport(anchor,j1,p1),j2,p2),j3,'O')[0]
        term=2/f;W+=term;by[p1+p2]=by.get(p1+p2,0)+term
require(W < Fraction(959,1000))
require(Fraction(959,1000)+3*Fraction(79,10000)*Fraction(6,5)**2 <1)
print('short three-step inverse scale subtotal =',format(float(W),'.12f'),'with exact <959/1000 check')

D=((39,-8),(11,-2)); identity=((1,0),(0,1))
def mm(A,B,p):
    return tuple(tuple(sum(A[i][k]*B[k][j] for k in range(2))%p for j in range(2)) for i in range(2))
def mx(A,x,p):
    return tuple(sum(v*w for v,w in zip(row,x))%p for row in A)
def matpow(A,n,p):
    R=identity
    while n:
        if n%2:R=mm(R,A,p)
        A=mm(A,A,p);n//=2
    return R
x=(3,1)
require(len({mx(matpow(D,n,7),x,7) for n in range(48)})==48 and matpow(D,48,7)==identity)
require(matpow(D,48,49)==((22,35),(7,36)))
require([n for n in range(8) if mx(matpow(D,n,7),x,7)[0]==0]==[3])
require(matpow(D,8,49)==((17,42),(28,24)) and mx(matpow(D,3,7),x,7)==(0,4))
require(mx(matpow(D,35,7**5),x,7**5)[0]%7**4==0 and mx(matpow(D,35,7**5),x,7**5)[0]!=0)
print('scalar 7-adic complete-cycle and simple zero-phase derivatives checked')

def primitive_dead(Z):
    H,h,L,E=Z
    return gcd(*Z)==1 and gcd(H,L,h)>1 and gcd(H,L,E)>1
T={(0,)}
for n in range(1,7):
    if n>1:T={canon(extend(A,c)) for A in T for c in range(1<<(n-1))}
    for A in T:
        if strong(A):
            require(not any(primitive_dead(state(A,v)) for v in range(n)))
S=(6,24,42,81,69,27,39)
require(strong(S) and state(S,1)==(135,33,165,25) and primitive_dead(state(S,1)))
def inflation(A,w,m,v):
    p=[]
    for i in range(len(A)):p.extend([i]*(m if i==w else 1))
    U=tuple(sum(1<<j for j,b in enumerate(p) if (a==b and i<j) or (a!=b and A[a]>>b&1)) for i,a in enumerate(p))
    return U,p.index(v)
for m in [1,2,5,9]:
    A,v=inflation(S,0,m,1)
    Z=state(A,v);H,h,L,E=Z
    require(strong(A) and h==16*(1<<m)+1 and E==16*(1<<m)-7)
    if m%4==1:
        require(tuple(x%15 for x in Z)==(0,3,0,10) and primitive_dead(Z))
print('minimum primitive dead root through order seven and three strong inflation controls verified')
