#!/usr/bin/env python3
"""Exact regression/hostile controls; no claims based only on finite search."""
import json
import math
import random
from collections import Counter
import sympy as s

x, y, t = s.symbols('x y t')
counts = Counter()

def gate(ok, label):
    counts[label] += 1
    if not ok:
        raise RuntimeError('Failed gate: ' + label)

def equal(a, b, label):
    gate(s.expand(a-b) == 0, label)

def bracket(P, Q):
    return s.expand(s.diff(P,x)*s.diff(Q,y)-s.diff(P,y)*s.diff(Q,x))

def component(P, p, q, weight, modulus=None):
    return s.Add(*(coef*x**i*y**j for (i,j),coef in s.Poly(P,x,y).terms()
                   if ((p*i-q*j-weight) % modulus == 0 if modulus
                       else p*i-q*j == weight)))

def mono(D):
    return [(i,j) for i in range(D+1) for j in range(D+1-i)]

rng = random.Random(114)
A = 2+3*t-2*t**2
B = 5-2*t+7*t**2
for p in range(1,5):
    for q in range(1,5):
        if math.gcd(p,q) != 1:
            continue
        T=x**q*y**p
        exactQ=y*B.subs(t,T)
        exactP=x*A.subs(t,T)
        expected=(A*B+q*t*s.diff(A,t)*B+p*t*A*s.diff(B,t)).subs(t,T)
        equal(bracket(exactP,exactQ),expected,'general_weighted_bracket')
        arbitraryP=sum(rng.randrange(-4,5)*x**i*y**j for i,j in mono(7))
        pp=component(arbitraryP,p,q,p)
        equal(component(bracket(arbitraryP,exactQ),p,q,0),
              bracket(pp,exactQ),'one_coordinate_channel_isolation')
        arbitraryQ=sum(rng.randrange(-4,5)*x**i*y**j for i,j in mono(7))
        qq=component(arbitraryQ,p,q,-q)
        equal(component(bracket(exactP,arbitraryQ),p,q,0),
              bracket(exactP,qq),'symmetric_channel_isolation')

for p in range(1,7):
    for q in range(1,7):
        if math.gcd(p,q) != 1:
            continue
        for d in range(1,13):
            M=max(q*(d-1),p*d+q)
            N=max(q*d+p,p*(d-1))
            for extra in (1,2,5):
                m=M+extra
                supported=[(i,j) for i,j in mono(d) if (p*i-q*j+q)%m==0]
                gate(all(p*i-q*j==-q for i,j in supported),'finite_Q_no_alias')
                m=N+extra
                supported=[(i,j) for i,j in mono(d) if (p*i-q*j-p)%m==0]
                gate(all(p*i-q*j==p for i,j in supported),'finite_P_no_alias')
        for a in range(6):
            for b in range(6):
                gate(1+q*a+p*b>0,'positive_top_coefficient')

# Exact support at the old full-map boundary, including the asymmetric case.
for w in range(1,11):
    for D in range(2,31):
        m=w*D+1
        ps={(i,j) for i,j in mono(D) if (i-w*j-1)%m==0}
        qs={(i,j) for i,j in mono(D) if (i-w*j+w)%m==0}
        expectedP={(i,j) for i,j in mono(D) if i-w*j==1}|{(0,D)}
        expectedQ={(i,j) for i,j in mono(D) if i-w*j==-w}
        if w==1:
            expectedQ.add((D,0))
        gate(ps==expectedP,'weighted_boundary_P_support')
        gate(qs==expectedQ,'weighted_boundary_Q_support')
        if w>1:
            gate(m>max(w*(D-1),D+w),'weighted_boundary_uses_one_coordinate')
        equal(bracket(2*x+3*y**D,5*y),10,'boundary_shear_is_Keller')

alpha,beta=s.symbols('alpha beta')
for D in range(2,10):
    T=x*y
    P=x*A.subs(t,T)+alpha*y**D
    Q=y*B.subs(t,T)+beta*x**D
    H=(A*B+t*(s.diff(A,t)*B+A*s.diff(B,t))-alpha*beta*D**2*t**(D-1)).subs(t,T)
    low=-alpha*D*y**(D+1)*s.diff(B,t).subs(t,T)
    high=-beta*D*x**(D+1)*s.diff(A,t).subs(t,T)
    equal(bracket(P,Q),H+low+high,'balanced_boundary_full_formula')
    equal(component(bracket(P,Q),1,1,0),H,'balanced_zero_channel')
    equal(component(bracket(P,Q),1,1,-(D+1)),low,'balanced_negative_channel')
    equal(component(bracket(P,Q),1,1,D+1),high,'balanced_positive_channel')
    equal(bracket(2*x+alpha*y**D,3*y+beta*x**D),
          6-alpha*beta*D**2*(x*y)**(D-1),'two_alias_obstruction')

for m in range(3,61):
    ps={(i,j) for i,j in mono(m) if (i-j-1)%m==0}
    qs={(i,j) for i,j in mono(m) if (i-j+1)%m==0}
    gate(ps=={(i,j) for i,j in mono(m) if i-j==1}|{(0,m-1)},'balanced_extra_shell_P')
    gate(qs=={(i,j) for i,j in mono(m) if i-j==-1}|{(m-1,0)},'balanced_extra_shell_Q')

# The companion is unrestricted: it may have no finite equivariance and much
# larger degree than the single symmetric coordinate's degree bound.
for p,q in ((1,1),(1,3),(2,3),(3,2)):
    for degree in (2,7,19,40):
        P=2*x+y**degree+3*y+7
        Q=3*y
        equal(bracket(P,Q),6,'unrestricted_companion_Keller')
        equal(component(Q,p,q,-q),Q,'single_coordinate_exact_weight')
        equal(P.subs({x:(x-(y/3)**degree-y-7)/2,y:y/3},simultaneous=True),
              x,'explicit_inverse_first_coordinate')
        equal(Q.subs(y,y/3),y,'explicit_inverse_second_coordinate')

for D in range(2,16):
    m=D+1
    Q=y+x**D
    equal(component(Q,1,1,-1,m),Q,'strict_bound_cyclic_hostile')
    gate(component(Q,1,1,-1)!=Q,'strict_bound_not_exact_weight')
    equal(bracket(x,Q),1,'strict_bound_hostile_still_Keller')

# D=1 is a genuine exception to the claimed triangular form, not an omission.
for w in range(1,15):
    m=w+1
    P=x+y; Q=x+2*y
    equal(component(P,1,w,1,m),P,'degree_one_scalar_P')
    equal(component(Q,1,w,-w,m),Q,'degree_one_scalar_Q')
    equal(bracket(P,Q),1,'degree_one_nondiagonal_Keller')

# Sharpness of the balanced degree-m shell at m=3, and a broader nontriangular
# family without claiming its degree is minimal for every modulus.
for m in range(3,8):
    P=x+y**(m-1)
    Q=y+P**(m-1)
    equal(bracket(P,Q),1,'two_shears_Keller')
    equal(component(P,1,1,1,m),s.expand(P),'two_shears_equivariant_P')
    equal(component(Q,1,1,-1,m),s.expand(Q),'two_shears_equivariant_Q')
    gate(s.Poly(Q,x,y).total_degree()==(m-1)**2,'two_shears_degree')
    gate(s.diff(P,y)!=0 and s.diff(Q,x)!=0,'two_shears_not_coordinate_triangular')

print(json.dumps({'status':'PASS','total_gates':sum(counts.values()),
                  'counts':dict(sorted(counts.items()))},indent=2,sort_keys=True))
