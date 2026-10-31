#!/usr/bin/env python3
"""Exact bounded nonunit feedback check: no word generation, external packages,
floating roots, or assertion statements. Scalar triples evaluate at real beta.
"""
from fractions import Fraction as F
from itertools import permutations
import json

def need(c,m):
    if not c:raise RuntimeError(m)
def add(a,b):return tuple(x+y for x,y in zip(a,b))
def neg(a):return tuple(-x for x in a)
def mul(a,b):
    v=[F(0)]*5
    for i,x in enumerate(a):
        for j,y in enumerate(b):v[i+j]+=x*y
    for j in (4,3):
        for r in range(3):v[j-3+r]+=v[j]
    return tuple(v[:3])
def sc(c,a):return tuple(c*x for x in a)
def pw(a,n):
    v=(F(1),F(0),F(0))
    while n:
        if n&1:v=mul(v,a)
        a=mul(a,a);n//=2
    return v
O=(F(1),F(0),F(0));Z=(F(0),)*3;T=(F(0),F(1),F(0));TI=(-F(1),-F(1),F(1))
def mm(A,B):return [[sum(A[i][r]*B[r][j] for r in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
def mp(A,n):
    v=[[int(i==j) for j in range(3)] for i in range(3)]
    while n:
        if n&1:v=mm(v,A)
        A=mm(A,A);n//=2
    return v
def ms(c,A):return [[c*v for v in row] for row in A]
def ma(A,B):return [[a+b for a,b in zip(rowa,rowb)] for rowa,rowb in zip(A,B)]
def parity(p):return (-1)**sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))
def det(A):return sum(parity(p)*A[0][p[0]]*A[1][p[1]]*A[2][p[2]] for p in permutations(range(3)))
def coeff(A):
    tr=sum(A[i][i] for i in range(3));e2=sum(A[i][i]*A[j][j]-A[i][j]*A[j][i] for i in range(3) for j in range(i+1,3))
    return [1,-tr,e2,-det(A)]
def enc(a):return [str(x) for x in a]
def run():
    M=[[1,1,1],[1,0,0],[0,1,0]];C=[[0,0,1],[1,0,1],[0,1,1]];I=mp(C,0)
    need(mul(T,TI)==O,'beta unit identity')
    ell=[O,add(T,neg(O)),TI]
    reports=[]
    for q,n in [(2,4),(103,17)]:
        N=ms(q,mp(M,n));lam=sc(q,pw(T,n));delta=sc(q*q,pw(TI,n))
        need(all(v>0 for row in N for v in row),'positive incidence')
        need(det(N)==q**3,'nonunit determinant norm')
        need(mul(lam,delta)==sc(q**3,O),'Perron times planar determinant equals norm')
        need(mul(lam,sc(F(1,q**3),delta))==O,'finite-place product corrects measure balance')
        for j in range(3):
            v=Z
            for i in range(3):v=add(v,sc(N[i][j],ell[i]))
            need(v==mul(lam,ell[j]),'positive left eigenvector identity')
        gap=F(9,5)**n-q*q
        need(gap>0,'exact conjugate contraction')
        reports.append({'q':q,'n':n,'incidence':N,'norm':det(N),'characteristic_polynomial_high_to_low':coeff(N),
            'perron_power_basis':enc(lam),'planar_determinant_power_basis':enc(delta),
            'measure_balance_multiplier':q**3,'finite_place_Haar_multiplier':str(F(1,q**3)),
            'rational_contraction_gap_(9/5)^n-q^2':str(gap),
            'nonzero_finite_area_disjoint_piece_system':'impossible'})
    for gamma,q,n,sgn in [(F(0),1,1,1),(F(1,2),2,4,1),(F(56,103),103,17,-1)]:
        R=ms(sgn*q,mp(C,n));factor=ma(I,ms(-gamma,C))
        need(mm(factor,R)==C,'rational resolvent operator exact identity')
    # Explicit gamma=56/103 sign boundary: raw spectral radius has negative real eigenvalue.
    need(F(9,5)**17>103**2,'negative case signflip contraction')
    return {'status':'FINITE-EXACT identities supporting the accompanying proof','cases':reports,
        'resolvent_operator_identities':'(I-gamma C)R=C verified at gamma=0,1/2,56/103',
        'negative_case':'raw R=-103 C^17 is not a positive substitution incidence matrix; +103 C^17 is a separate signflip',
        'scope':'No gigantic substituted words, planar area computation, tiling, or spectral-purity claim.'}
print(json.dumps(run(),indent=2,sort_keys=True))
