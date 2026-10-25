#!/usr/bin/env python3
"""Independent controls; imports no character-packet implementation."""
from itertools import combinations, product
from math import comb
from collections import Counter
import json
import sympy as S
N=0

def check(ok, name):
 global N
 N+=1
 if not ok: raise RuntimeError(name)

rows=[]
for n in range(4,17):
 degree=n-6
 hilbert=0 if degree<0 else sum((-1)**j*comb(n-3,j)*comb(degree-2*j+n-1,n-1) for j in range(min(n-3,degree//2)+1))
 chars=sum(comb(n,w)*comb(w//2-1,2) for w in range(6,n+1,2))
 # Squarefree root monomials over three degree-two base variables.
 parity=sum(comb(n,a)*comb((degree-a)//2+2,2) for a in range(max(0,degree+1)) if a<=n and (degree-a)%2==0)
 check(hilbert==chars==parity, 'three dimensions')
 rows.append((n,hilbert))
# Complement characters and complete canonical linear representation at n=7.
all7=127
canonical=[all7^(1<<i) for i in range(7)]
check(all(x.bit_count()==6 for x in canonical),'six-line complementary characters')
for j in range(7):
 values=[(-1)**((c>>j)&1) for c in canonical]
 check(values.count(1)==1 and values[j]==1,'one K3 form fixed under deletion')
check(sum(1 for I in combinations(canonical,3) if I[0]^I[1]^I[2]==0)==0,'no full-deck alternating invariant')
# Hamming kernel, represented by even words: Fano labels 1,...,7.
H=[]
for mask in range(128):
 syndrome=0
 for i in range(7):
  if mask>>i&1: syndrome^=i+1
 if syndrome==0 and mask.bit_count()%2==0: H.append(mask)
check(len(H)==8,'kernel representatives')
characters=[tuple((-1)**((c&m).bit_count()%2) for m in H) for c in canonical]
check(len(set(characters))==7,'distinct restricted characters')
for k,m in enumerate(H):
 check(sum(c[k] for c in characters)==(7 if m==0 else -1),'reduced regular trace')
inv=0
for I in combinations(range(7),3):
 if all(characters[I[0]][j]*characters[I[1]][j]*characters[I[2]][j]==1 for j in range(8)): inv+=1
check(inv==7,'seven restricted alternating invariants')
# D4 local blowup directly: w^2=xy(x+y), x=s,y=st,w=sv.
x,y,w,s,t,v=S.symbols('x y w s t v')
f=w*w-x*y*(x+y)
check(S.expand(f.subs({x:s,y:s*t,w:s*v},simultaneous=True)-s*s*(v*v-s*t*(1+t)))==0,'D4 first blowup')
# Exceptional plus strict branch has three distinct marked directions 0,-1,infinity.
check(S.factor(s*t*(1+t))==s*t*(t+1),'D4 exceptional is branch')
# Scalar i on the quadratic Veronese coordinates of C²/{±1} gives minus all roots.
a,b=S.symbols('a b')
for q in (a*a,a*b,b*b):
 check(S.expand(q.subs({a:S.I*a,b:S.I*b},simultaneous=True)+q)==0,'dependent triple fixed exceptional')
D=S.Matrix([[-2,1,1,1],[1,-2,0,0],[1,0,-2,0],[1,0,0,-2]])
check(D.det()==4 and all(e<0 for e in D.eigenvals()),'D4 resolution configuration')
print(json.dumps({'status':'PASS','checks':N,'n_pg':rows,'exterior_cube_dimensions':[0,7],'scope':'Finite controls plus separate proof audit; no producer imports'},sort_keys=True,indent=2))
