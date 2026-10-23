from itertools import product
from math import comb
import json
checks=0

def need(c,msg):
 global checks
 checks+=1
 if not c: raise RuntimeError(msg)

def points(r):return list(product(range(3),repeat=r))
def coeff(f,r,m):
 out={}
 for a in points(r):
  out[a]=sum((-1)**(sum(a)-sum(b))*__import__('functools').reduce(lambda x,y:x*y,(comb(x,y) for x,y in zip(a,b)),1)*f[b] for b in product(*(range(x+1) for x in a)))%m
 return out

def evaluate(c,a,m):
 return sum(v*__import__('functools').reduce(lambda x,y:x*y,(comb(x,y) if y<=x else 0 for x,y in zip(a,b)),1) for b,v in c.items())%m

def shift(a,v):return tuple((x+y)%3 for x,y in zip(a,v))
def kernel(f,r):return [v for v in points(r) if all(f[shift(a,v)]==f[a] for a in points(r))]
def audit(f,r,m):
 c=coeff(f,r,m);need(all(evaluate(c,a,m)==f[a] for a in points(r)),'integer Newton inversion')
 for i in range(r):
  v=tuple(int(j==i) for j in range(r));literal=any(f[a]!=f[shift(a,v)] for a in points(r))
  need(literal==any(x and a[i] for a,x in c.items()),'essential axis')
 K=kernel(f,r);need(all(shift(u,v) in K for u,v in product(K,repeat=2)),'translation subgroup')
 classes={frozenset(shift(a,v) for v in K) for a in points(r)}
 need(len(classes)*len(K)==3**r,'retain all quotient lifts')
 need(all(len({f[a] for a in cl})==1 for cl in classes),'quotient sufficient')

for values in product(range(9),repeat=3):audit(dict(zip(points(1),values)),1,9)
for values in product(range(2),repeat=9):audit(dict(zip(points(2),values)),2,9)

def T(v):return v[1:]+v[:1]
def D(v,m):return tuple((a-b)%m for a,b in zip(T(v),v))
for k in range(1,7):
 m=3**k
 for v in ((1,0,0),(0,1,0),(0,0,1)):
  w=v
  for j in range(1,2*k+2):
   w=D(w,m)
   if j==2*k:need(any(w),'sharp lower index')
  need(not any(w),'nilpotency upper index')
for a,b in product(range(8),repeat=2):
 v=(a,b,(-a-b)%8)
 for _ in range(12):v=D(v,8)
 need(any(v)==bool(a or b),'zero-sum action injective at modulus8')
f={a:sum(a)%3 for a in points(2)};c=coeff(f,2,9);K=kernel(f,2)
need(set(K)=={(0,0),(1,2),(2,1)},'rank-one period observer')
need(c[1,0]!=0 and c[0,1]!=0,'Newton support indices spuriously span rank two')
hostile_coefficients={str(a):v for a,v in c.items()}
for vals in product(range(9),repeat=3):
 f=dict(zip(points(1),vals));c=coeff(f,1,9);d=coeff(dict(zip(points(1),T(vals))),1,9)
 want=((c[0,]+c[1,])%9,(c[1,]+c[2,])%9,(-3*c[1,]-2*c[2,])%9)
 need(tuple(d[a] for a in points(1))==want,'exact cyclic wrap matrix')
print(json.dumps({'status':'PASS','checks':checks,'one_axis_mod9_tables':729,'two_axis_boolean_tables':512,'sharp_single_axis_nilpotency_indices':{str(3**k):2*k+1 for k in range(1,7)},'newton_support_hostile':{'kernel':K,'coefficients':hostile_coefficients }},indent=2))
