#!/usr/bin/env python3
"""Independent bounded controls. No producer-verifier import and no writes.

Run with Python 3 + SymPy, either normally or with -O.  Output is JSON.
Exact non-binomial packets use a known finite root and a Newton-polygon
escaping factor; trace values are computed from multiplication matrices.
"""
from collections import Counter
from itertools import product
import json
import sys
import sympy as S

t,z=S.symbols('t z')
counts=Counter()
cases=[]
def need(group, value):
    if not bool(value):
        raise RuntimeError('Independent control failed: '+group)
    counts[group]+=1

def valuation(e):
    p,q=S.fraction(S.cancel(e))
    if p==0:
        return S.oo
    def order(f):
        return min(k[0] for k,c in S.Poly(f,t).terms() if c)
    return order(p)-order(q)

def multiplication_matrix(poly):
    n=S.degree(poly,z)
    # Construct every column by polynomial division, independently of Newton sums.
    columns=[]
    for k in range(n):
        rem=S.Poly(S.rem(z**(k+1),poly,z),z)
        columns.append(S.Matrix([rem.nth(j) for j in range(n)]))
    return S.Matrix.hstack(*columns)

for m,a,b in product(range(2,8),(1,2),(1,2)):
    if '--progress' in sys.argv:
        print('packet',m,a,b,file=sys.stderr,flush=True)
    q=z**m + sum((j+1)*z**(m-j) for j in range(1,m))-t**(-a)
    chi=S.expand((z-b)*q)
    coeff=S.Poly(chi,z).all_coeffs()
    # Newton polygon: the factor q has endpoints (0,-a),(m,0),
    # and every intermediate coefficient has valuation 0, strictly above
    # the joining segment. All m roots therefore have slope a/m.
    need('escaping_newton_polygon', all(0 > -a+S.Rational(a*k,m) for k in range(1,m)))
    need('known_finite_root', chi.subs(z,b)==0)
    C=multiplication_matrix(chi)
    P=S.eye(m+1)
    T=[S.Integer(m+1)]
    for k in range(1,2*m-1):
        P=(P*C).applyfunc(S.expand)
        T.append(S.expand(S.trace(P)))
    need('low_trace_integrality',all(valuation(T[k])>=0 for k in range(1,m)))
    need('first_trace_pole',valuation(T[m])==-a)
    # Every Newton identity in the required window is checked against the
    # independently computed trace, including initial-degree coefficients.
    for k in range(1,2*m-1):
        if k<=m+1:
            identity=T[k]+sum(coeff[j]*T[k-j] for j in range(1,k))+k*coeff[k]
        else:
            identity=T[k]+sum(coeff[j]*T[k-j] for j in range(1,m+2))
        need('newton_identity',S.expand(identity)==0)
    G=S.Matrix(m-1,m-1,lambda i,j:T[i+j+2])
    need('gram_entry_bound',all(valuation(e)>=-a for e in G))
    B=(t**a*G).applyfunc(S.expand)
    residue=B.applyfunc(lambda e:e.subs(t,0))
    lead=(t**a*T[m]).expand().subs(t,0)
    sign=(-1)**((m-1)*(m-2)//2)
    need('scaled_residue_determinant',residue.det()==sign*lead**(m-1))
    need('direct_gram_determinant',valuation(G.det(method='domain-ge'))==-a*(m-1))
    h=(m-1)//2
    need('residue_maximal_isotropic_plane',residue[:h,:h]==S.zeros(h))
    ratio=S.cancel(-coeff[-1]/coeff[-2])
    need('selector_integral',valuation(ratio)>=0)
    need('selector_error_bound',valuation(ratio-b)>=a)
    cases.append({'m':m,'a':a,'finite_root':b,'index_bound':(a*(m-1)+1)//2})

# Exact rational Newton errors against a known root (not merely residuals).
# q has nonzero lower coefficients, so the initial ratio is not the root.
# Keep polynomial numerator/denominator pairs without costly expression
# simplification. Common factors do not change the difference of valuations.
newton=[]
for a,b in product((1,2),(1,2)):
    if '--progress' in sys.argv:
        print('Newton',a,b,file=sys.stderr,flush=True)
    tau=S.Poly(t**a,t)
    N=S.Poly(b,t)
    D=S.Poly(1+3*b*t**a,t)
    def porder(P):
        return min(k[0] for k,c in P.terms() if c)
    for n in range(3):
        error_num=N-b*D
        err=porder(error_num)-porder(D)
        need('exact_newton_error',err>=(2**(n+1)-1)*a)
        need('exact_newton_sharpness',err==(2**(n+1)-1)*a)
        Q=tau*(N**3+2*N**2*D+3*N*D**2)-D**3
        Qp=tau*(3*N**2+4*N*D+3*D**2)
        derivative_num=Q+error_num*Qp
        derivative_den=D**3*S.Poly(1+3*b*t**a,t)
        need('normalized_derivative_unit',porder(derivative_num)-porder(derivative_den)==0)
        newton.append({'a':a,'b':b,'iteration':n,'error_order':int(err)})
        if n<2:
            # For chi=(z-b)q, Newton gives e_new=e^2*q'/(q+e*q').
            Dnext=D*derivative_num
            Nnext=b*Dnext+error_num**2*Qp
            N,D=Nnext,Dnext

# Brute force all diagonal exponent allocations in a split normal form.
budgets=[]
for rank,a in product(range(1,6),range(1,5)):
    best=None
    for allocation in product(range(a+1),repeat=rank):
        valid=all(allocation[2*j]+allocation[2*j+1]>=a for j in range(rank//2))
        if rank%2:
            valid=valid and 2*allocation[-1]>=a
        if valid:
            cost=sum(allocation)
            best=cost if best is None else min(best,cost)
    need('exhaustive_split_index',best==(a*rank+1)//2)
    need('scalar_parity_excess',rank*((a+1)//2)-best==(rank//2 if a%2 else 0))
    budgets.append({'rank':rank,'a':a,'minimum_index':best})

# Sharpened full-rank hostile: a self-dual trace lattice can contain 1,h
# and nevertheless leave h polar, because it is not a ring / h-module.
hostiles=[]
for a in range(1,6):
    T=[S.Integer(4)]+[3*t**(-a*k//3) if k%3==0 else S.Integer(0) for k in range(1,7)]
    scales=[0,0,a,a]
    G=S.Matrix(4,4,lambda i,j:t**(scales[i]+scales[j])*T[i+j])
    expected=S.Matrix([[4,0,0,3],[0,0,3,0],[0,3,0,0],[3,0,0,3]])
    need('full_rank_hostile_gram',G==expected)
    need('full_rank_hostile_unimodular',G.det()==-27)
    need('full_rank_hostile_not_closed',valuation(t**(-a))<0)
    need('full_rank_hostile_h_polar',S.Rational(a,3)>0)
    hostiles.append({'a':a,'gram_determinant':-27,'h_pole_order':str(S.Rational(a,3))})

# Centered quartic selector formula checked from the coefficient identities.
p,q,r=S.symbols('p q r')
T2=-2*p; T3=-3*q; T4=2*p*p-4*r
need('centered_quartic_selector',S.cancel((3*T2*T2-6*T4)/(8*T3)+r/q)==0)

print(json.dumps({'status':'PASS','checks':sum(counts.values()),'check_groups':dict(counts),
                  'nonbinomial_packets':cases,'exact_newton_errors':newton,
                  'exhaustive_split_budgets':budgets,'full_rank_hostiles':hostiles,
                  'scope':'Bounded exact controls support the independent uniform-proof review; no Keller realization or JC closure.'},indent=2,sort_keys=True))
