#!/usr/bin/env python3
"""Exact arithmetic audits for the independent tensor proof. No assert checks."""
from fractions import Fraction as F
from itertools import product
from random import Random
import json
rng=Random(139526)
checks={}
def check(cond,label):
    if not cond:raise RuntimeError(label)
    checks[label]=True

def mul(u,v):
    a,b=u;c,d=v
    return a*c-b*d,a*d+b*c+b*d

def add(u,v):return u[0]+v[0],u[1]+v[1]
def scale(c,u):return c*u[0],c*u[1]
def conj(u):return u[0]+u[1],-u[1]
def norm(u):return u[0]**2+u[0]*u[1]+u[1]**2
z=(0,1);zc=(1,-1)
check(mul(z,z)==(-1,1),'zeta_quadratic')
check(conj(z)==zc,'zeta_conjugation')
check(mul(z,zc)==(1,0),'zeta_unit')
check(mul(mul(z,z),z)==(-1,0),'zeta_cube')
for u in product(range(-3,4),repeat=2):
    check(mul(u,conj(u))==(norm(u),0),'norm_identity_'+str(u))
    if u!=(0,0):check(norm(u)>=1,'nonzero_norm_'+str(u))
for u in [(1,0),(0,1),(1,1)]:
    check(norm(u)%2==1,'F4_unit_norm_'+str(u))
    check(any(tuple(t%2 for t in mul(u,v))==(1,0) for v in [(1,0),(0,1),(1,1)]),'F4_inverse_'+str(u))
check(all((a*a+a+1)%2==1 for a in [0,1]),'F4_irreducible_quadratic')
# Exact matrices and Q2=M^T M.
Q=[[[F(1),0],[0,F(1)]],[[F(1),F(1,2)],[F(1,2),F(1)]],[[F(1),F(1)],[F(1),F(2)]]]
M=[[1,1],[0,1]]
check([[sum(M[k][i]*M[k][j] for k in range(2)) for j in range(2)] for i in range(2)]==Q[2],'Q2_integral_Gram')
check([M[i][1]%2 for i in range(2)]==[1,1],'Q2_parity_two_odd_slices')
for x in [F(0),F(1,7),F(1,2),F(6,7),F(1)]:
    weights=[(1-x)**2,2*x*(1-x),x*x]
    reconstructed=[[sum(weights[k]*Q[k][i][j] for k in range(3)) for j in range(2)] for i in range(2)]
    check(reconstructed==[[1,x],[x,1+x*x]],'Bernstein_matrix_'+str(x))
    check(sum(weights[k]*[1,1,2][k] for k in range(3))==1+x*x,'Bernstein_target_'+str(x))

def evaluate(c,zs):
    out=(0,0)
    for mask,coef in enumerate(c):
        term=(1,0)
        for i,zv in enumerate(zs):
            if mask>>i&1:term=mul(term,zv)
        out=add(out,scale(coef,term))
    return out

def quadratic(c,alpha):
    value=F(0)
    for i,a in enumerate(c):
        if not a:continue
        for j,b in enumerate(c):
            if not b:continue
            v=F(a*b)
            for k,typ in enumerate(alpha):v*=Q[typ][(i>>k)&1][(j>>k)&1]
            value+=v
    return value

phase_cases=0
for r in range(6):
    for trial in range(3):
        c=[2*rng.randrange(-2,3) for _ in range(1<<r)]
        c[-1]=2*rng.randrange(-2,3)+1
        vals=[evaluate(c,zs) for zs in product([z,zc],repeat=r)]
        check(all(norm(v)>=1 and tuple(t%2 for t in v)!=(0,0) for v in vals),'phase_nonvanishing_r%d_t%d'%(r,trial))
        mean=sum(F(norm(v)) for v in vals)/(1<<r)
        check(mean==quadratic(c,[1]*r),'phase_Gram_identity_r%d_t%d'%(r,trial))
        phase_cases+=1

coefficient_cases=0
for r in range(5):
    for trial in range(4):
        c=[2*rng.randrange(-2,3) for _ in range(1<<r)]
        c[-1]=2*rng.randrange(-2,3)+1
        for alpha in product(range(3),repeat=r):
            check(quadratic(c,alpha)>=2**alpha.count(2),'tensor_coefficient_r%d_t%d_a%s'%(r,trial,alpha))
            coefficient_cases+=1
# Dropping parity is genuinely invalid: P=x-1, at x=1.
check(quadratic([-1,1],[2])==1<2,'negative_control_non_even_constant')
# Sharp unit-top polynomials attain every target coefficient.
for r in range(5):
    c=[0]*(1<<r);c[-1]=1
    for alpha in product(range(3),repeat=r):
        check(quadratic(c,alpha)==2**alpha.count(2),'sharp_top_monomial_r%d_a%s'%(r,alpha))
print(json.dumps({'status':'PASS','checks_passed':len(checks),'phase_cases':phase_cases,'tensor_coefficient_cases':coefficient_cases,'negative_control':'P=x-1 fails when even-proper-coefficient hypothesis is removed','all_checks_active_under_optimization':True},sort_keys=True,indent=2))
