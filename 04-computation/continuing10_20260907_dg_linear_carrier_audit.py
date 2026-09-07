#!/usr/bin/env python3
"""Independent sparse-polynomial, rational-rank and two-chart DG audit."""
from fractions import Fraction as Q
from itertools import combinations_with_replacement
from pathlib import Path
from hashlib import sha256
import sys
import sympy as s
sys.stdout.reconfigure(encoding='utf-8',newline='\n')
G=0
def need(value,label):
    global G
    G+=1
    if not bool(value):raise ArithmeticError(label)
def zero(value,label):need(s.cancel(value)==0,label)

# Native Q[x,t] dictionaries, separate from the producer's symbolic Jacobians.
def clean(P):return {k:Q(v) for k,v in P.items() if v}
def add(P,R):
    Z=dict(P)
    for k,v in R.items():Z[k]=Z.get(k,0)+v
    return clean(Z)
def scale(P,c):return clean({k:v*c for k,v in P.items()})
def mul(P,R):
    Z={}
    for (i,j),v in P.items():
        for (k,l),w in R.items():Z[i+k,j+l]=Z.get((i+k,j+l),0)+v*w
    return clean(Z)
def power(P,n):
    Z={(0,0):Q(1)}
    for _ in range(n):Z=mul(Z,P)
    return Z
def derivative(P,axis):
    Z={}
    for k,v in P.items():
        if not k[axis]:continue
        j=list(k);j[axis]-=1;Z[tuple(j)]=k[axis]*v
    return clean(Z)
def bracket(P,R):return add(mul(derivative(P,0),derivative(R,1)),scale(mul(derivative(P,1),derivative(R,0)),-1))
def make_F(A,B):return add({(i,1):v for i,v in enumerate(A) if v},{(i,0):v for i,v in enumerate(B) if v})

def rank(matrix):
    M=[[Q(v) for v in row] for row in matrix]
    if not M:return 0
    i=0
    for j in range(len(M[0])):
        pivot=next((k for k in range(i,len(M)) if M[k][j]),None)
        if pivot is None:continue
        M[i],M[pivot]=M[pivot],M[i];d=M[i][j]
        M[i]=[v/d for v in M[i]]
        for k in range(i+1,len(M)):
            if M[k][j]:
                d=M[k][j];M[k]=[a-d*b for a,b in zip(M[k],M[i])]
        i+=1
        if i==len(M):break
    return i

# Complete negative-Laurent coefficient matrices, built monomial by monomial.
layers=[]
for N in range(0,17):
    columns=[]
    for i in range(N+1):
        columns.append({(2-i,0):-1,(4-i,1):-1})  # x^i t
    columns.extend({(-i,0):1} for i in range(N+1))
    bad=sorted({k for C in columns for k in C if k[0]<0})
    matrix=[[C.get(k,0) for C in columns] for k in bad]
    dimension=len(columns)-rank(matrix)
    need(dimension==min(N+2,6),'native full Laurent nullity, with arbitrary x degree in the finite control')
    layers.append([N,dimension])

# Unpruned polynomial-mate systems and exact image controls on fresh cases.
cases=[([3],[1,-2,0,1],True),([0,1],[2,0,1],False),([1,0,1],[2,-3,1],False),
       ([-8,12,-6,1],[0,1],False),([81,-108,54,-12,1],[0,-12,1],False),
       ([4,4,-3,-2,1],[0,-2,1],False)]
systems=[]
for ci,(A,B,success) in enumerate(cases):
    F=make_F(A,B)
    for nt in [0,1,2,4]:
        monomials=[(i,j) for j in range(nt+1) for i in range(9)]
        columns=[bracket(F,{m:Q(1)}) for m in monomials]
        rows=sorted({(0,0)}|{m for P in columns for m in P})
        matrix=[[P.get(m,0) for P in columns] for m in rows]
        target=[int(m==(0,0)) for m in rows]
        attained=rank(matrix)==rank([row+[v] for row,v in zip(matrix,target)])
        need(attained==success,'fresh full polynomial coefficient image system')
        systems.append([ci,nt,attained])
    for n in [1,2,3,5,7]:
        q={(i,0):Q((-1)**i,i+1) for i in range(1,6)}
        trial=add(power(F,n),q)
        want=scale(mul({(i,0):v for i,v in enumerate(A) if v},derivative(q,0)),-1)
        need(bracket(F,trial)==want,'unbounded descent identity on native full powers and independent residual')

# Complete small root-multiplicity bank checks all degree cases of criticality.
def poly_roots(roots):
    P=[Q(1)]
    for r in roots:
        Z=[Q(0)]*(len(P)+1)
        for i,c in enumerate(P):Z[i]-=r*c;Z[i+1]+=c
        P=Z
    return P
def value(A,x):return sum(a*x**i for i,a in enumerate(A))
critical=[]
for degree in range(5):
    for roots in combinations_with_replacement([-2,0,3],degree):
        A=poly_roots(roots);coeff=A+[Q(0)]*(5-len(A));a0,a1,a2,a3,a4=coeff
        B=[Q(7),a3,a4];Ap=[i*A[i] for i in range(1,len(A))];Bp=[a3,2*a4]
        affine=any(value(Ap,r)!=0 or value(Bp,r)==0 for r in set(roots))
        boundary=(a4==0 and (a3!=0 or a1==0))
        expected=degree==4 and all(roots.count(r)>=2 and r!=0 for r in set(roots))
        need((not affine and not boundary)==expected,'entire small multiplicity/zero-root critical-point bank')
        critical.append([list(roots),affine,boundary])

# Whole-parameter identities, using explicit special-fibre factors and inverses.
x,t,r,b,z,h,k,a,c,Y=s.symbols('x t r b z h k a c Y')
A=a*(x-h)**2*(x-k)**2;B=a*(x*x-2*(h+k)*x)+c
F=A*t+B
Fi=s.cancel(F.subs({x:1/r,t:-r*r-r**4*b},simultaneous=True))
zero(s.diff(F,x).subs(x,h)+2*a*k,'all-parameter first affine root gradient')
zero(s.diff(F,x).subs(x,k)+2*a*h,'all-parameter second affine root gradient')
zero(s.diff(Fi,b).subs(r,0)+a,'boundary separating gradient is a nonzero constant')
# Compute the simple-pole coefficient by differentiation of the regular factor.
res=s.diff(-1/(a*(x-k)**2),x).subs(x,h)
zero(res-2/(a*(h-k)**3),'independent local residue extraction, without a residue API')
four=F.subs(k,h);four_i=Fi.subs(k,h);f0=c-3*a*h*h;w=1-h*r
E2=h*h*(3-h*r)+b*w**3
zero(four_i-f0+a*w*E2,'complete second-chart special fibre factorization')
zero(E2.subs(r,1/h)-2*h*h,'the two component factors are disjoint in the boundary chart')
zero(w*(-h*h-b*w*w)/(2*h*h)-(1-E2/(2*h*h)),'explicit inverse of 1-hr in the E2 coordinate ring')
zero(four-f0-a*(x-h)*((x-h)**3*t+x-3*h),'complete first-chart factorization')
zero(((x-h)**3*t+x-3*h).subs(x,h)+2*h,'first-chart component disjointness and reducedness')

# Elimination recovers E2 rather than merely testing its proposed parametrization.
BZ=-2*h**5*z**3-7*h**4*z*z-8*h**3*z-3*h*h
el=s.resultant(r*(1+h*z)-z,b-BZ,z)
zero(el-E2,'native two-equation elimination gives the exact E2 equation')
zero(E2.subs({r:z/(1+h*z),b:BZ},simultaneous=True),'parameter lies on the entire boundary-chart component')
zero((r/(1-h*r)).subs(r,z/(1+h*z))-z,'parameter and inverse agree')
zero((1+h*(r/(1-h*r)))-1/(1-h*r),'two-chart localization transition')
zero((1+h*z)-h*z-1,'the two parameter domains cover A1')
zero(BZ.subs(z,0)+3*h*h,'exact E2 boundary point')
zero((h+1/z).subs(z,-1/h),'the omitted boundary-chart parameter belongs to U0 with x0')

# The full generic fibre also has the stated A1 completion.
X=h+1/z;TZ=(Y-f0)*z**4/a+2*h*z**3-z*z
zero(four.subs({x:X,t:TZ},simultaneous=True)-Y,'generic fibre source parametrization')
generic_b=s.cancel(-X*X-X**4*TZ)
need(s.denom(generic_b)==a,'generic boundary coordinate has no z pole')
zero(generic_b.subs(z,0)-(c-6*a*h*h-Y)/a,'generic fibre includes the actual D point')

# Rational mate and divisor orders: a rational correction cannot repair both.
g0=1/(3*a*(x-h)**3)
zero(s.diff(four,x)*s.diff(g0,t)-s.diff(four,t)*s.diff(g0,x)-1,'actual rational conjugate')
zero(g0.subs(x,X)-z**3/(3*a),'rational conjugate regular on E2')
zero(s.limit((four-f0)/(x-h),x,h)+2*a*h,'special value has order1 on E1')
zero(s.limit(g0*(x-h)**3,x,h)-1/(3*a),'rational mate has order minus3 on E1')
zero(s.diff(four_i,b).subs({r:0,b:-3*h*h})+a,'special value has order1 on E2 at its D point')
# Clearing E1's leading pole forces a nonzero degree-three principal coefficient.
Hminus3=8*a*a*h**3/3
zero(1/(3*a)+Hminus3/(-2*a*h)**3,'necessary third-order principal coefficient at E1')
need(Hminus3!=0,'symbolically nonzero correction principal coefficient for a*h nonzero')

# Scope controls: affine-plane positive pair, vanishing volume and outer powers.
U0=t+x**3;mate=-x
zero(s.diff(U0,x)*s.diff(mate,t)-s.diff(U0,t)*s.diff(mate,x)-1,'positive affine-plane polynomial pair')
zero(s.limit(r*mate.subs(x,1/r),r,0)+1,'that mate has a genuine boundary pole')
jac_transition=s.det(s.Matrix([[s.diff(1/r,r),s.diff(1/r,b)],[s.diff(-r*r-r**4*b,r),s.diff(-r*r-r**4*b,b)]]))
zero(jac_transition-r*r,'global source form vanishes to order2 on D')
for n in [2,3,5]:
    zero(s.diff(four**n,x)*s.diff(g0,t)-s.diff(four**n,t)*s.diff(g0,x)-n*four**(n-1),'outer chain rule retains the nonconstant factor')

print('GLOBAL_LAYER native Laurent matrices through degree16:',layers)
print('ALL_DEGREE_IMAGE independent sparse identities;24 unpruned rational coefficient systems:',systems)
print('CRITICALITY complete finite multiplicity bank:',len(critical),'root multisets; quartic double/double or quadruple nonzero roots exactly')
print('SPECIAL_FIBRE both charts; E2 resultant, inverse of1-hr, A1 localization cover, generic A1 boundary and both divisor orders PASS')
print('SCOPE constant source Jacobian on the specified DG surface; omega has divisor2D; no arbitrary-coordinate or JC conclusion')
print('PASS',G,'always-active exact gates; raw LF; no producer imports')
