#!/usr/bin/env python3
"""Exact algebra controls; no physical-state inference is tested."""
import json, itertools, math
from collections import Counter
import sympy as s
from sympy.matrices.normalforms import smith_normal_form

def require(x,msg):
    if not x: raise RuntimeError(msg)
def snf(M):
    D=smith_normal_form(M,domain=s.ZZ)
    return [abs(int(D[i,i])) for i in range(M.rows)]
def skew(n,mask):
    S=s.zeros(n)
    for j,(a,b) in enumerate(itertools.combinations(range(n),2)):
        S[a,b]=1 if (mask>>j)&1 else -1
        S[b,a]=-S[a,b]
    return S
hist=Counter(); witnesses={}
for mask in range(64):
    S=skew(4,mask); K=s.zeros(8); K[:4,4:]=S; K[4:,:4]=S.T
    d=int(S.det()); hist[d]+=1
    require(d%2==1 and math.isqrt(d)**2==d,'odd square determinant')
    require(K==K.T and all(K[i,i]%2==0 for i in range(8)),'even symmetric')
    require(abs(K.det())==d*d,'BF determinant')
    R=s.diag(-1,-1,-1,-1,1,1,1,1)
    require(R.T*K*R==-K,'converse equivalence')
    witnesses.setdefault(str(d),{'mask':mask,'S':list(map(list,S.tolist())),'S_snf':snf(S),'BF_snf':snf(K)})
K=s.Matrix([[3,4],[4,3]]); t=s.Matrix([1,1]); KB=s.Matrix([[0,1],[1,0]])
require(K==KB+3*t*t.T,'three flux attachment')
require(K==-s.eye(2)+4*t*t.T,'reverse Jain alternative')
require((t.T*K.inv()*t)[0]==s.Rational(2,7),'filling')
require(K.det()==-7,'determinant')
# Gradient holonomy on every simple oriented cyclic listing for n <= 6.
grad=0
for n in range(3,7):
    for c in itertools.product(range(3),repeat=n):
        require(sum(c[(i+1)%n]-c[i] for i in range(n))%3==0,'exact gradient circulation')
        grad+=1
# q_+ and q_- on Z/256 are well-defined but are not isometric.
for x in range(256):
    require(((x+256)**2-x*x)%512==0,'q periodicity')
require(not any((u*u+1)%512==0 for u in range(1,512,2)),'opposite q isometry obstruction')
# Odd K=[3] naive spin fails on representatives 0 and 3.
require((s.Rational(3*3,2*3)%1)==s.Rational(1,2),'odd K representative ambiguity')
# A triangle grounded Laplacian gives even K and Z/3, with q(e1)=1/3.
L=s.Matrix([[2,-1],[-1,2]])
require(snf(L)==[1,3] and L.inv()[0,0]/2==s.Rational(1,3),'triangle discriminant')
receipt={'status':'PASS','scope':'exact finite algebra only','n4_tournament_count':64,'skew_det_histogram':dict(hist),'witnesses':witnesses,'two_over_seven':{'K':K.tolist(),'inverse':[[str(x) for x in row] for row in K.inv().tolist()],'nu':'2/7','SNF':snf(K),'charge_e1':'1/7','representative_spin_e1':'exp(-3 pi i /7)'},'gradient_checks':grad,'z256_opposite_refinements_nonisometric':True,'odd_K_naive_spin_ambiguity_confirmed':True}
print(json.dumps(receipt,indent=2,default=int))
