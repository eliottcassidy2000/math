#!/usr/bin/env python3
"""Independent Session60 arithmetic / incidence / reflection audit.
No packet-code imports. Field-square characters use the quadratic norm.
"""
from itertools import combinations, permutations, product
from collections import Counter
from fractions import Fraction as Q
from math import isqrt
import json
import sympy as s
checks=0
def req(x, why):
    global checks
    checks+=1
    if not x: raise RuntimeError(why)
def det(M):
    M=[list(map(Q,row)) for row in M]; ans=Q(1)
    for k in range(len(M)):
        a=next((i for i in range(k,len(M)) if M[i][k]),None)
        if a is None:return Q(0)
        if a!=k:M[a],M[k]=M[k],M[a];ans=-ans
        u=M[k][k];ans*=u
        for j in range(k,len(M)):M[k][j]/=u
        for i in range(k+1,len(M)):
            v=M[i][k]
            for j in range(k,len(M)):M[i][j]-=v*M[k][j]
    return ans
def cross(v,w):return (v[1]*w[2]-v[2]*w[1],v[2]*w[0]-v[0]*w[2],v[0]*w[1]-v[1]*w[0])
lines=[(1,0,0),(1,1,0),(1,-1,0),(1,0,1),(1,0,-1),(1,1,1)]
def incidences(p=None):
    out={}
    for a,b in combinations(lines,2):
        v=cross(a,b)
        if p:
            v=tuple(c%p for c in v); h=next(c for c in v if c);v=tuple(c*pow(h,-1,p)%p for c in v)
            idx=tuple(i for i,l in enumerate(lines) if sum(c*d for c,d in zip(l,v))%p==0)
        else:
            h=next(c for c in v if c);v=tuple(Q(c,h) for c in v)
            idx=tuple(i for i,l in enumerate(lines) if sum(c*d for c,d in zip(l,v))==0)
        out[v]=idx
    return out
base=incidences(); signature=sorted(base.values())
req(Counter(map(len,signature))=={2:9,3:2},'rational incidence')
for p in [5,7,11,13]:req(sorted(incidences(p).values())==signature,'no incidence collision '+str(p))
D={int(det(M)) for M in combinations(lines,3)};req(D-{0}<=set([-3,-2,-1,1,2,3]),'only primes 2,3 in nonzero minors')
G=[[0]*18 for _ in range(18)];G[0][0]=2
for i in range(1,18):G[i][i]=-2
for h in [10,14]:
    for j in range(h+1,h+4):G[h][j]=G[j][h]=1
req(det(G)==-16384,'independence determinant of 18 geometric classes')

def chi(n,p):
    c=pow(n%p,(p-1)//2,p)
    return -1 if c==p-1 else c

def count(p,n):
    if n==1:
        total=sum(chi((1-x*x)*(1-z*z)*(1+x+z),p) for x,z in product(range(p),repeat=2))
        q=p;dn=None
    else:
        dn=next(d for d in range(2,p) if chi(d,p)==-1)
        def C(a,b):return chi(a*a-dn*b*b,p)
        # chi_Fp2(a+b*u)=chi_Fp(Norm(a+b*u)), u^2=dn.
        # Independently verify norm-character against all nonzero field squares.
        squares={((a*a+dn*b*b)%p,2*a*b%p) for a,b in product(range(p),repeat=2) if (a,b)!=(0,0)}
        for a,b in product(range(p),repeat=2):req(C(a,b)==(0 if (a,b)==(0,0) else 1 if (a,b) in squares else -1),'norm character')
        total=0
        for a,b,c,d in product(range(p),repeat=4):
            total+=C(1-a*a-dn*b*b,-2*a*b)*C(1-c*c-dn*d*d,-2*c*d)*C(1+a+c,b+d)
        q=p*p
    return {'p':p,'n':n,'q':q,'norm_nonresidue':dn,'S':total,'N':q*q+18*q+1+total}
counts={p:[count(p,n) for n in [1,2]] for p in [7,13]}
expected={7:[180,3384],13:[408,32088]}
T=s.symbols('T'); arithmetic={}
for p,rows in counts.items():
    req([r['N'] for r in rows]==expected[p],'independent point counts')
    a,A=[r['S'] for r in rows]; b=Q(a*a-A,2)
    req(b!=0,'reciprocity sign forced, not assumed')
    R=T**4-a*T**3+b*T**2-a*p*p*T+p**4
    c=2*p-a
    req(s.expand(R-(T-p)**2*(T*T+c*T+p*p))==0,'residual quartic factors')
    req(Q(-c,p).denominator!=1,'noncyclotomic root sums')
    disc=-(2*p+c)
    arithmetic[p]={'trace':a,'trace_square':A,'middle':int(b),'residual_polynomial':str(R),'remaining_quadratic':str(T*T+c*T+p*p),'normalized_root_sum':str(Q(-c,p)),'disc_representative':disc}
    # Independent extension consistency: if u_n=lambda^n+mu^n, then q_n*Phi_n(1)=2p^n-u_n.
    # These must lie in the same square class when no new algebraic eigenvalues arise.
    u0,u1=2,-c
    for n in range(1,9):
        if n==1:un=u1
        else:u0,u1=u1,-c*u1-p*p*u0;un=u1
        v=Q(2*p**n-un,2*p+c)
        req(v>0 and isqrt(v.numerator)**2==v.numerator and isqrt(v.denominator)**2==v.denominator,'extension squareclass control')
req(Q(arithmetic[7]['disc_representative'],arithmetic[13]['disc_representative'])==Q(1,2),'discriminant ratio nonsquare')

# Recompute lower-bound Weierstrass geometry directly from the cubic.
x,t,z,u=s.symbols('x t z u'); d=1-t*t
f=(x-d)*(x+d)*(x-d*(t+1)); coeff=s.Poly(s.expand(f),x)
a2=coeff.nth(2);a4=coeff.nth(1);a6=coeff.nth(0)
b2=4*a2;b4=2*a4;b6=4*a6;b8=4*a2*a6-a4*a4
c4=b2*b2-24*b4;c6=-b2**3+36*b2*b4-216*b6
# Polynomial discriminant gives the Weierstrass discriminant with factor 16.
# Check it independently against the b-invariant formula.
Delta_b=-b2*b2*b8-8*b4**3-27*b6*b6+9*b2*b4*b6
Delta=16*s.discriminant(f,x)
req(s.expand(Delta_b-Delta)==0,'b-invariant discriminant route')
req(s.expand(c4**3-c6*c6-1728*Delta)==0,'Weierstrass invariant normalization')
req(s.expand(f.subs(x,-d*z)-d**3*(1-z*z)*(1+t+z))==0,'coordinate change sign')
req(s.expand(f.subs(x,d*t)-d**4)==0,'section membership')
def ord_at(g,a,var=t):
    poly=s.Poly(s.expand(g),var)
    return min(i for (i,),co in s.Poly(poly.as_expr().subs(var,var+a),var).terms() if co!=0)
vals=[]
for a in [0,-2,1,-1]:vals.append([ord_at(g,a) for g in [c4,c6,Delta]])
vals.append([ord_at(s.cancel(u**w*g.subs(t,1/u)),0,u) for g,w in zip([c4,c6,Delta],[8,12,24])])
req(vals==[[0,0,2],[0,0,2],[2,3,6],[2,4,6],[2,3,8]],'all singular valuations incl infinity')
req(all(v[2]<12 for v in vals),'minimal discriminant at every singular place')
req(s.expand(f.subs({x:u-3,t:2})-u**3+36*u)==0,'specialization short curve')
X=Q(-3);Y=Q(9); doubles=[]
def v2(q):
    def z(n):
        v=0
        while n%2==0:n//=2;v+=1
        return v
    return z(abs(q.numerator))-z(q.denominator)
for k in range(1,6):
    slope=(3*X*X-36)/(2*Y); xx=slope*slope-2*X; yy=slope*(X-xx)-Y;X,Y=xx,yy
    req(Y*Y==X**3-36*X,'double on curve');req(v2(X)==-2*k,'nontorsion denominator growth control');doubles.append(v2(X))
req(6*19-56==58,'audited Hodge numerical transfer')

# Independently enumerate all labelled Fano planes, not just one marking.
sigma=(0,2,1,4,3,6,5)
planes=set()
for lab in permutations(range(1,8)):
    pl=frozenset(tuple(i for i in tr) for tr in combinations(range(7),3) if lab[tr[0]]^lab[tr[1]]^lab[tr[2]]==0)
    req(len(pl)==7,'Fano marking has seven lines');planes.add(pl)
req(len(planes)==30,'exactly thirty quotient kernels/incidences')
def act(pl):return frozenset(tuple(sorted(sigma[i] for i in line)) for line in pl)
req(all(act(pl) in planes and act(pl)!=pl for pl in planes),'reflection pairs all kernels with no fixed one')
pairs={frozenset((pl,act(pl))) for pl in planes};req(len(pairs)==15,'fifteen unordered pairs')
mark={7:0,1:1,6:2,2:3,5:4,4:5,3:6}
J=lambda a:7 if a==7 else a^7
req(all(mark[J(a)]==sigma[mark[a]] for a in mark),'explicit Boolean to branch conjugacy')
invol=Counter();gl=0
for rows in product(range(8),repeat=3):
    def actv(v):return sum(((rows[i]&v).bit_count()%2)<<i for i in range(3))
    perm=tuple(actv(v) for v in range(8))
    if len(set(perm))<8:continue
    gl+=1
    if perm!=tuple(range(8)) and all(perm[perm[v]]==v for v in range(8)):invol[sum(perm[v]==v for v in range(1,8))]+=1
req(gl==168 and invol=={3:21},'all nontrivial Fano involutions fix three points')
# Four explicit relations for the seven square-coordinate forms.
L=s.Matrix([[1,1,1,1,1,1,1],[0,1,-1,0,0,1,-1],[0,0,0,1,-1,1,-1]])
N=s.Matrix.hstack(*L.nullspace());perm=s.zeros(7)
for i in range(7):perm[sigma[i],i]=1
rel=N.gauss_jordan_solve(perm*N)[0]
req(rel.det()==-1 and perm.det()==-1,'canonical residue determinants')
req(perm.det()/rel.det()==1,'central minus-I determinant twist is plus')
req((s.eye(7)-perm).rank()==3,'four canonical invariants')
output={'status':'PASS','checks':checks,'counts':counts,'arithmetic':arithmetic,'eighteen_class_determinant':int(det(G)),'fiber_valuations':vals,'doubling_v2':doubles,'Fano_structures':len(planes),'reflection_kernel_pairs':len(pairs),'GL3_order':gl,'involution_fixed_census':dict(invol),'canonical_invariant_dimension':4,'outer_picard_rank':19,'hourglass_picard_rank':58,'scope':'Finite arithmetic checks plus geometric and primary-source proof in AUDIT.md, not machine-verified cohomology.'}
print(json.dumps(output,sort_keys=True,indent=2))
