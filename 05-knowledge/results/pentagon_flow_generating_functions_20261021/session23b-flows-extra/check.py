#!/usr/bin/env python3
"""Exact certificates for R5 flow series; no fitted numerators or floating point."""
from collections import Counter
from functools import lru_cache
from itertools import combinations
import json
import sympy as s


def require(b,label):
    if not b:raise RuntimeError(label)

q=s.symbols('q')
# Cone in the saturated integer lattice with coordinates t,y0,...,y4.
rays=[(1,0,0,0,0,0)]
for i in range(5):rays.append((1,)+tuple(int(j==i)for j in range(5)))
for i in range(5):rays.append((1,)+tuple(int(j in(i,(i+2)%5))for j in range(5)))
rays.append((2,1,1,1,1,1))
facets=[]
for i in range(5):facets.append((0,)+tuple(int(j==i)for j in range(5)))
for i in range(5):facets.append((1,)+tuple(-int(j in(i,(i+1)%5))for j in range(5)))
require(all(sum(a*b for a,b in zip(n,v))>=0 for n in facets for v in rays),'rays feasible')
# Exhaust all active-constraint sets on t=1. The polytope is bounded since 0<=y_i<=1.
vertices=set(); nonsingular=0
for active in combinations(range(10),5):
    M=s.Matrix([(1,0,0,0,0,0)]+[facets[i]for i in active])
    if not M.det():continue
    nonsingular+=1
    v=M.inv()*s.Matrix([1,0,0,0,0,0])
    if all(sum(a*b for a,b in zip(n,v))>=0 for n in facets):vertices.add(tuple(v))
expected={tuple(s.Rational(x,r[0])for x in r)for r in rays}
require(vertices==expected and len(vertices)==12,'complete exact vertex enumeration')
zeros=[frozenset(i for i,v in enumerate(rays) if sum(a*b for a,b in zip(n,v))==0)for n in facets]
require(all(s.Matrix.hstack(*(s.Matrix(rays[i])for i in z)).rank()==5 for z in zeros),'all ten inequalities are facets')

@lru_cache(None)
def triang(face):
    face=frozenset(face)
    if len(face)==1:return (tuple(face),)
    apex=min(face)
    proper={face&z for z in zeros if face&z and face&z!=face}
    facet_faces=[f for f in proper if not any(f<g for g in proper)]
    rank=s.Matrix.hstack(*(s.Matrix(rays[i])for i in face)).rank()
    require(all(s.Matrix.hstack(*(s.Matrix(rays[i])for i in f)).rank()==rank-1 for f in facet_faces),'recursive face codimension')
    out=[]
    for f in sorted(facet_faces,key=lambda x:tuple(sorted(x))):
        if apex not in f:
            out.extend(tuple(sorted((apex,)+cell))for cell in triang(tuple(sorted(f))))
    return tuple(out)

cells=triang(tuple(range(12)))
require(len(cells)==25 and len(set(cells))==25,'pulling simplex count')
# Positive combination gives an interior vector, checked not on any simplex hyperplane.
w=s.Matrix([sum(2**i*r[j]for i,r in enumerate(rays))for j in range(6)])
require(all(sum(a*b for a,b in zip(n,w))>0 for n in facets),'strictly interior generic vector')
terms=[]; records=[]; groups={}
for cell in cells:
    M=s.Matrix.hstack(*(s.Matrix(rays[i])for i in cell))
    det=M.det();require(abs(det)==1,'unimodularity in circulation lattice')
    alpha=M.inv()*w;require(all(a!=0 for a in alpha),'generic offset')
    weights=tuple(5*rays[i][0]-sum(rays[i][1:])for i in cell)
    open_indices=[i for i,a in zip(cell,alpha)if a<0]
    shift=sum(wt for wt,a in zip(weights,alpha)if a<0)
    terms.append(q**shift/s.prod(1-q**wt for wt in weights))
    key=tuple(sorted(weights));groups.setdefault(key,Counter())[shift]+=1
    records.append({'rays':cell,'determinant':int(det),'offset_coordinates':[int(a)for a in alpha],'open_ray_indices':open_indices,'weight_shift':shift,'ray_weights':weights})
P=1+2*q**3+3*q**4+q**5+3*q**6+q**7+3*q**8+2*q**9+q**12
R=P/((1-q**3)**3*(1-q**4)**2*(1-q**5))
require(s.cancel(sum(terms)-R)==0,'raw R5 rational identity from all simplices')
A=1+q**3+3*q**4+2*q**6+4*q**7+q**8+3*q**10
B=2*q**3+2*q**6+2*q**7+q**9+2*q**10+q**13
require(groups[(3,3,4,4,5,5)]==Counter({0:1,3:1,4:3,6:2,7:4,8:1,10:3}),'15-cone group')
require(groups[(3,3,3,4,5,5)]==Counter({3:2,6:2,7:2,9:1,10:2,13:1}),'10-cone group')
require(s.expand(A*(1-q**3)+B*(1-q**4)-P*(1-q**5))==0,'short numerator cancellation certificate')
N=1+q**4+q**5+q**6+q**7+q**8+q**12
K=N/((1-q**3)*(1-q**5)*(1-q**6)*(1-q**8))
geo=1/(1-1/q)*(1/(1-q**2)*(1/((1-q**6)*(1-q**3))-q**2/((1-q**8)*(1-q**5)))-1/q/(1-q)*(1/((1-q**5)*(1-q**3))-q/((1-q**6)*(1-q**4))))
require(s.cancel(geo-K)==0,'elementary anti-fixed geometric sum')
F=(R+4/(1-q**5)**2)/5
# Known interior-translation reciprocity is just a cross-check, not used as a proof of R or K.
require(all(s.cancel(f.subs(q,1/q)-q**10*f)==0 for f in(R,F,K,(F+K)/2)),'R5 reciprocal degree-gap control')

def gf_coeffs(num,weights,bound):
    coeff=[0]*(bound+1)
    for power,c in s.Poly(num,q).terms():
        if power[0]<=bound:coeff[power[0]]=int(c)
    for wt in weights:
        for k in range(wt,bound+1):coeff[k]+=coeff[k-wt]
    return coeff

# Independent integer-point counting directly in the pentagon inequalities.
M=40
raw=[0]*(M+1)
for t in range(2*M//5+1):
    for a in range(t+1):
      for b in range(t-a+1):
       for c in range(t-b+1):
        for d in range(t-c+1):
         for e in range(min(t-d,t-a)+1):
            weight=5*t-a-b-c-d-e
            if weight<=M:raw[weight]+=1
require(raw==gf_coeffs(P,(3,3,3,4,4,5),M),'independent raw pentagon count through 40')
MK=80
fixed=[0]*(MK+1)
for t in range(2*MK//5+1):
    for b in range(t//2+1):
      for a in range(t-b+1):
       for c in range(t-a+1):
        weight=5*t-2*a-2*b-c
        if weight<=MK:fixed[weight]+=1
require(fixed==gf_coeffs(N,(3,5,6,8),MK),'independent anti-fixed cone sum through 80')
orb=[(raw[m]+4*(m//5+1 if m%5==0 else 0))//5 for m in range(M+1)]
require(all((raw[m]+4*(m//5+1 if m%5==0 else 0))%5==0 for m in range(M+1)),'Burnside integrality')
require(all(orb[m]>=fixed[m]and (orb[m]-fixed[m])%2==0 for m in range(M+1)),'converse orbit integrality')
# Reduced cyclotomic denominators and leading constants, exact arithmetic.
require(s.gcd(P,s.prod(1-q**w for w in(3,3,3,4,4,5)))==1,'raw denominator reduced')
require(s.gcd(N,s.prod(1-q**w for w in(3,5,6,8)))==1,'fixed denominator reduced')
raw_leading=s.Rational(P.subs(q,1),s.factorial(5)*3**3*4**2*5)
fixed_leading=s.Rational(N.subs(q,1),s.factorial(3)*3*5*6*8)
require(raw_leading==s.Rational(17,259200),'raw leading coefficient')
require(fixed_leading==s.Rational(7,4320),'fixed leading coefficient')

out={
 'status':'PASS',
 'vertex_certificate':{'active_sets_tested':252,'nonsingular_active_sets':nonsingular,'vertices':[[str(x)for x in v]for v in sorted(vertices)],'facet_zero_ray_sets':[sorted(z)for z in zeros]},
 'triangulation':{'rays':rays,'generic_interior_vector':[int(x)for x in w],'simplex_count':25,'all_unimodular':True,'cells':records,'grouped_shifts':[{'weights':key,'shift_counts':dict(sorted(val.items()))}for key,val in sorted(groups.items())]},
 'rational_formulas':{'raw_numerator':str(P),'raw_denominator_weights':[3,3,3,4,4,5],'anti_fixed_numerator':str(N),'anti_fixed_denominator_weights':[3,5,6,8],'ordinary_orbit_formula':'(R(q)+4/(1-q^5)^2)/5'},
 'independent_counts':{'raw_total_use_bound':M,'raw_coefficients':raw,'ordinary_orbit_coefficients':orb,'fixed_total_use_bound':MK,'fixed_coefficients':fixed},
 'asymptotic_constants':{'raw_m5':str(raw_leading),'ordinary_orbit_m5':str(raw_leading/5),'anti_fixed_m3':str(fixed_leading),'fixed_fraction_m_minus2':str(fixed_leading/(raw_leading/5))},
 'scope':'Exact R5 certificate only. No fitted numerator, floating point triangulation, or general R_n classification.'
}
print(json.dumps(out,sort_keys=True,indent=2))
