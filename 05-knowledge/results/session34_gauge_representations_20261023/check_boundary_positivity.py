#!/usr/bin/env python3
import itertools,json
from fractions import Fraction as F
import sympy as S
G=list(itertools.permutations(range(3)));idx={g:i for i,g in enumerate(G)};e=(0,1,2)
def mul(a,b):return tuple(a[b[i]] for i in range(3))
def inv(a):return tuple(a.index(i) for i in range(3))
def sign(g):return (-1)**sum(g[i]>g[j] for i in range(3) for j in range(i+1,3))
def chi(g):return sum(g[i]==i for i in range(3))-1
def w(g,eps):return F(1)+eps*sign(g)+F(1,4)*chi(g)
checks=0
def ck(x,m):
 global checks
 checks+=1
 if not x:raise ValueError(m)
rows=[]
for eps in [F(1,4),F(-1,4)]:
 vals=[w(g,eps) for g in G]
 ck(min(vals)>0 and sum(vals)==6,'positive normalized density')
 K=S.Matrix([[S.Rational(w(mul(inv(g),h),eps))/6 for h in G] for g in G])
 sig=S.Matrix([sign(g) for g in G]);std=S.Matrix([chi(g) for g in G])
 ck(K*sig==S.Rational(eps)*sig,'sign eigenvector')
 ck(K*std==S.Rational(1,8)*std,'standard eigenvector')
 expected={S.Integer(1):1,S.Rational(eps):1,S.Rational(1,8):4}
 ck(K.eigenvals()==expected,'full eigenvalues')
 Z=F(0)
 for a,b,c in itertools.product(G,repeat=3):
  # Tree-gauge fields g01=g02=g03=e, g12=a,g13=b,g23=c.
  Z+=w(a,eps)*w(b,eps)*w(c,eps)*w(mul(mul(a,c),inv(b)),eps)
 Z/=216
 ck(Z==F(1029,1024),'tetrahedral partition')
 rows.append({'epsilon':str(eps),'class_values_identity_transposition_threecycle':[str(w(g,eps)) for g in [e,(1,0,2),(1,2,0)]],'tetrahedron_Z':str(Z),'reflection_form_sign':str(eps),'spectrum':{str(k):v for k,v in expected.items()}})
# Generator divided by log2: each transposition rate1/3, each nonidentity3cycle2/3.
rate={g:(F(1,3) if sign(g)==-1 else F(2,3)) for g in G if g!=e}
L=S.zeros(6)
for i,g in enumerate(G):
 for a,r in rate.items():L[i,idx[mul(a,g)]]+=S.Rational(r);L[i,i]-=S.Rational(r)
ck(L.eigenvals()=={S.Integer(0):1,S.Integer(-2):1,S.Integer(-3):4},'Markov generator')
ck(all(L[i,j]>=0 for i in range(6) for j in range(6) if i!=j),'generator positivity')
ck(all(sum(L.row(i))==0 for i in range(6)),'generator conservation')
# Exp(log2 L) checked by exact spectral projectors, avoiding transcendental evaluation.
P0=(L+2*S.eye(6))*(L+3*S.eye(6))/6
P2=-L*(L+3*S.eye(6))/2
P3=L*(L+2*S.eye(6))/3
Kplus=S.Matrix([[S.Rational(w(mul(inv(g),h),F(1,4)))/6 for h in G] for g in G])
ck(P0+P2/4+P3/8==Kplus,'positive continuous semigroup time1')
print(json.dumps({'status':'PASS','checks':checks,'weights':rows,'generator_rates_divided_by_log2':{'each_transposition':'1/3','each_threecycle':'2/3'},'scope':'finite one-step reflected-boundary kernel; equality of all closed triangulated orientable surface amplitudes follows from character formula and even face count'},indent=2,sort_keys=True))
