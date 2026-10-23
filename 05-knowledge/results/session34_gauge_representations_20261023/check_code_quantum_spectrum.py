#!/usr/bin/env python3
"""Exact finite representation reduction; not a continuum field-theory proof."""
import itertools,json
from collections import Counter
from fractions import Fraction
import sympy as S
checks=0
def ck(p,msg):
 global checks
 checks+=1
 if not p:raise ValueError(msg)
def dot(a,b):return (a&b).bit_count()%2
cols=list(range(1,8))
def syndrome(m):
 s=0
 for i,a in enumerate(cols):
  if m>>i&1:s^=a
 return s
C={m for m in range(128) if syndrome(m)==0}
D={sum(dot(b,a)<<i for i,a in enumerate(cols)) for b in range(8)}
ck(len(C)==16 and len(D)==8,'dimensions')
ck(all(dot(c,d)==0 for c in C for d in D),'orthogonal codes')
remaining=set(range(128));cosets=[]
while remaining:
 q=min(remaining);O={q^d for d in D};remaining-=O;cosets.append(O)
counts=Counter(tuple(sorted(Counter(q.bit_count() for q in O).items())) for O in cosets)
expected={((0,1),(4,7)):1,((1,1),(3,4),(5,3)):7,((2,3),(4,4),(6,1)):7,((3,7),(7,1)):1}
ck(dict(counts)==expected,'sector enumerators')
# Verify full finite Fourier actions entrywise using integers.
# X_i character q = (-1)^q_i character q; zero-syndrome projector Fourier kernel.
for p in range(128):
 for q in range(128):
  coeff=sum((-1)**dot(p^q,m) for m in C)
  ck(coeff==(16 if (p^q) in D else 0),'potential Fourier kernel')
for q in range(128):
 ck(sum(1-(-1)**((q>>i)&1) for i in range(7))==2*q.bit_count(),'kinetic Fourier')
t,z,c=S.symbols('t z c');rows=[];coeffs=[]
for weights,mult in sorted(counts.items()):
 vals=[w for w,n in weights]; ns=[n for w,n in weights]
 B=S.diag(*[1+2*t*w for w in vals])-S.ones(len(vals),1)*S.Matrix([ns])/8
 pol=S.factor(B.charpoly(z).as_expr())
 ck(S.diff(pol,z).subs({z:0,t:0})!=0,'simple low branch')
 trial=7*t-7*t*t+c*t**3
 ser=S.Poly(S.series(pol.subs(z,trial),t,0,4).removeO(),t)
 ck(all(ser.coeff_monomial(t**k)==0 for k in range(3)),'common first two orders')
 sol=S.solve(ser.coeff_monomial(t**3),c)
 ck(len(sol)==1,'unique cubic coefficient')
 # Direct 8x8 Fourier-sector characteristic polynomial.
 seq=[w for w,n in weights for _ in range(n)]
 M=S.diag(*[1+2*t*w for w in seq])-S.ones(8)/8
 fact=pol*S.prod((z-1-2*t*w)**(n-1) for w,n in weights)
 ck(S.expand(M.charpoly(z).as_expr()-fact)==0,'sector factorization')
 coeffs.append(sol[0]);rows.append({'weights':weights,'sector_multiplicity':mult,'reduced_charpoly':str(pol),'low_branch':'7t-7t^2+('+str(sol[0])+')t^3+O(t^4)','extra_eigenvalues':[[str(1+2*t*w),n-1] for w,n in weights if n>1]})
ck(coeffs==[-42,-6,6,42],'cubic coefficients')
# Invariant sector, exact full8x8 quotient directly in syndrome coordinates.
H=t*(8*S.eye(8)-S.ones(8))+S.diag(0,*([1]*7))
char=(z-1-8*t)**6*(z*z-(1+8*t)*z+7*t)
ck(S.expand(H.charpoly(z).as_expr()-char)==0,'invariant spectrum')
# Seven weight3 codewords have cube-corner graph eigenvalues7,1,-1,-7.
lines=[c for c in C if c.bit_count()==3]
cayley=Counter(sum((-1)**dot(next(iter(O)),c) for c in lines) for O in cosets)
ck(cayley==Counter({7:1,1:7,-1:7,-7:1}),'third order effective graph')
print(json.dumps({'status':'PASS','checks':checks,'sectors':rows,'effective_code_graph_spectrum':dict(cayley),'full_gap_asymptotic':'36 h^3/J^2 + O(h^4/J^3), fixed h>0, J/h -> infinity','invariant_gap':'sqrt((J-6h)^2+28h^2)','scope':'finite128-state quantum Hamiltonian; code symmetry is not identified with physical Yang-Mills gauge'},indent=2,sort_keys=True))
