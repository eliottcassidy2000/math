#!/usr/bin/env python3
"""Exact controls for the simplex/sphere/cube completion packet.
Standard-library only. No repository writes or imports. Tests supplement proofs.
"""
from fractions import Fraction as F
from itertools import product, combinations
from math import comb
import json
from pathlib import Path

COUNTS = {}
def check(name, cond):
    if not cond:
        raise RuntimeError(name)
    COUNTS[name] = COUNTS.get(name, 0) + 1

def chi(a, x):
    return -1 if (a & x).bit_count() % 2 else 1

def even_masks(m):
    return [a for a in range(1 << m) if a.bit_count() % 2 == 0]

def moments(m, weights):
    # Representatives have last sign +1; global-sign orbits are exactly these.
    return {a: sum(w * chi(a, x) for x, w in enumerate(weights)) for a in even_masks(m)}

def invert(m, mus):
    n = 1 << (m-1)
    return [sum(mus[a] * chi(a, x) for a in mus)/n for x in range(n)]

def matmul(a,b):
    return [[sum(a[i][k]*b[k][j] for k in range(len(b))) for j in range(len(b[0]))] for i in range(len(a))]
def transpose(a): return list(map(list,zip(*a)))
def hadamard(n): return [[chi(i,j) for j in range(n)] for i in range(n)]
def householder(v):
    s=sum(x*x for x in v)
    return [[F(int(i==j))-F(2*v[i]*v[j],s) for j in range(len(v))] for i in range(len(v))]

# Every even-character record except one can miss maximally disjoint laws.
for m in range(2,7):
    n=1<<(m-1); chars=even_masks(m)
    for a in chars[1:]:
        plus=[F(1+chi(a,x),n) for x in range(n)]
        minus=[F(1-chi(a,x),n) for x in range(n)]
        mp,mm=moments(m,plus),moments(m,minus)
        check('parity_normalization',sum(plus)==sum(minus)==1)
        check('parity_disjoint',all(p*q==0 for p,q in zip(plus,minus)))
        for b in chars:
            check('parity_moments',mp[b]-mm[b]==(2 if a==b else 0))
    # Independent matrix/Walsh and atom-inversion paths, five rational laws.
    for seed in range(1,6):
        raw=[(seed*(x+3)**2+7*x+1)%23+1 for x in range(n)]
        w=[F(v,sum(raw)) for v in raw]
        mus=moments(m,w)
        check('walsh_inverse',invert(m,mus)==w)
        h=[[chi(a,x) for x in range(n)] for a in chars]
        moment=[[mus[a^b] for b in chars] for a in chars]
        gram=[[sum(w[x]*h[i][x]*h[j][x] for x in range(n)) for j in range(n)] for i in range(n)]
        check('moment_gram',moment==gram)
        for x in range(n):
            hv=[h[i][x] for i in range(n)]
            check('moment_eigen',all(sum(moment[i][j]*hv[j] for j in range(n))==n*w[x]*hv[i] for i in range(n)))

# Triangle four-character spectrum, in a complete 9^3 rational grid.
triangle_psd_only=0; triangle_boolean=0
for a,b,c in product([F(i,4) for i in range(-4,5)],repeat=3):
    lam=[1+a+b+c,1+a-b-c,1-a+b-c,1-a-b+c]
    mus={0:F(1),3:a,5:b,6:c}
    w=invert(3,mus)
    check('triangle_spectrum',sorted(lam)==sorted(4*z for z in w))
    boolean=min(lam)>=0
    ordinary=(abs(a)<=1 and abs(b)<=1 and abs(c)<=1 and 1+2*a*b*c-a*a-b*b-c*c>=0)
    if boolean: triangle_boolean+=1
    if ordinary and not boolean: triangle_psd_only+=1
    check('triangle_boolean_implies_psd',not boolean or ordinary)
check('triangle_hostile',1+3*F(-1,2)<0 and 1+2*F(-1,2)**3-3*F(-1,2)**2==0)

# Full characters are essential: any one-row deletion can miss negative mass.
missing_row=[]
for m in range(3,7):
    n=1<<(m-1); c=F(-1,n-2)
    mus={a:(F(1) if a==0 else c) for a in even_masks(m)}
    w=invert(m,mus)
    check('missing_row_bad_mass',w[0]==F(-1,n*(n-2)))
    check('missing_row_other_masses',all(v==F(n-1,n*(n-2)) for v in w[1:]))
    check('missing_row_submatrix_psd',1+(n-2)*c==0 and 1-c>0)
    check('missing_row_full_negative',1+(n-1)*c<0)
    missing_row.append({'roots':m,'matrix_size':n,'nonconstant_moment':str(c),'negative_orbit_mass':str(w[0]),'other_orbit_mass':str(w[1])})

# The four-root pair projection remains realizable when its proposed quartic is not.
m=4;n=8
uniform=[F(1,n)]*n
balanced=[F(1,3) if x.bit_count()==2 else F(0) for x in range(n)]
w=[(p+q)/2 for p,q in zip(uniform,balanced)]
mu=moments(4,w)
check('four_root_pair_control',all(mu[a]==F(-1,6) for a in mu if a.bit_count()==2))
check('four_root_quartic_control',mu[15]==F(1,2))

# Inherited even-star coefficient and new one-character consumer discrepancy.
stars=[]
for m in range(2,11,2):
    full=(1<<m)-1
    values=[min(x.bit_count(),m-x.bit_count()) for x in range(1<<m)]
    coeff=F(sum(values[x]*chi(full,x) for x in range(1<<m)),1<<m)
    r=m//2; formula=F((-1)**r*comb(2*r-2,r-1),2**(2*r-1))
    check('star_top_coefficient',coeff==formula)
    ep=F(sum(values[x] for x in range(1<<m) if chi(full,x)==1),1<<(m-1))
    em=F(sum(values[x] for x in range(1<<m) if chi(full,x)==-1),1<<(m-1))
    check('star_parity_gap',ep-em==2*coeff)
    stars.append({'roots':m,'coefficient':str(coeff),'plus_cost':str(ep),'minus_cost':str(em)})

# Metric speed of a parity fiber; no floating point/arcsin numeric test needed.
for m in range(2,7):
    n=1<<(m-1); a=even_masks(m)[1]
    for t in [F(-3,4),F(-1,2),F(0),F(1,3),F(4,5)]:
        fisher=sum((F(chi(a,x),n)**2)/F(1+t*chi(a,x),n) for x in range(n))
        check('fisher_parity_speed',fisher==1/(1-t*t))

# Squared-entry orthogonal couplings generally fail to compose.
for d in range(2,7):
    u=householder(list(range(1,d+1)));v=householder([2*i+1 for i in range(d)])
    vu=matmul(v,u)
    bu=[[x*x for x in row] for row in u];bv=[[x*x for x in row] for row in v]
    classical=matmul(bv,bu)
    defect=[[vu[i][j]**2-classical[i][j] for j in range(d)] for i in range(d)]
    for i,j in product(range(d),repeat=2):
        cross=2*sum(v[i][k]*v[i][l]*u[k][j]*u[l][j] for k,l in combinations(range(d),2))
        check('coupling_composition_identity',defect[i][j]==cross)
    check('coupling_defect_row_zero',all(sum(row)==0 for row in defect))
    check('coupling_defect_column_zero',all(sum(row)==0 for row in transpose(defect)))
walsh=[]
for k in range(1,6):
    n=1<<k; h=hadamard(n); hh=matmul(h,h)
    check('walsh_involution',all(hh[i][j]==n*int(i==j) for i,j in product(range(n),repeat=2)))
    # U=H/sqrt(n), entrywise square J/n; B(UU)=I, B(U)^2=J/n.
    tv=sum(abs(F(int(j==0))-F(1,n)) for j in range(n))/2
    check('walsh_composition_tv',tv==1-F(1,n))
    walsh.append({'size':n,'total_variation_gap':str(tv)})

# Regular-simplex Gram feasibility versus odd-cardinality Boolean parity.
equicorrelation=[]
for m in range(2,13):
    allowed=[x for x in range(1<<m) if x.bit_count() in {m//2,(m+1)//2}]
    c=F(-1,m if m%2 else m-1)
    for i,j in combinations(range(m),2):
        actual=F(sum(chi((1<<i)|(1<<j),x) for x in allowed),len(allowed))
        check('equicorrelation_attainment',actual==c)
    minimum_square=F(sum((m-2*x.bit_count())**2 for x in allowed),len(allowed))
    check('equicorrelation_parity',minimum_square==(m%2))
    check('equicorrelation_psd_gap',c-F(-1,m-1)==(F(1,m*(m-1)) if m%2 else 0))
    equicorrelation.append({'roots':m,'boolean_minimum':str(c),'psd_minimum':str(F(-1,m-1))})

# q17 valid cut gap -> certified correlation separation, inherited finite bound.
edges=2601; maxcut_upper=1730; witness=F(-1,3)
valid_sum=edges-2*maxcut_upper; witness_sum=edges*witness
check('q17_gap',valid_sum-witness_sum==8)
check('q17_l2_denominator',51*51==edges)

out={'status':'PASS','exact_check_total':sum(COUNTS.values()),'checks':COUNTS,
     'triangle_grid':{'triples':729,'boolean_feasible':triangle_boolean,'ordinary_psd_nonboolean':triangle_psd_only},
     'missing_row_hostiles':missing_row,'star_controls':stars,'walsh_composition_controls':walsh,
     'equicorrelation_boundary':equicorrelation,
     'q17_correlation_distance_lower_bounds':{'l1':'8','l2':'8/51','linf':'8/2601'},
     'scope':'Finite exact controls supplement the elementary proofs. The q17 cut bound is inherited, not re-proved by this script. No major conjecture is claimed solved.'}
print(json.dumps(out,indent=2,sort_keys=True))
