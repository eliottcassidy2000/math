#!/usr/bin/env python3
from itertools import product
from collections import deque
from pathlib import Path
import hashlib,json
# An independent positive-coordinate implementation; no author code imported.
M={
 'O':((1,0,1,0),(1,2,1,2),(0,0,0,0),(0,0,1,1)),
 'Z':((2,1,1,1),(0,1,0,0),(0,0,1,1),(0,0,0,0)),
 'J':((0,1,0,1),(0,0,0,0),(1,0,1,1),(0,0,1,0))}
def require(v,msg):
 if not v:raise RuntimeError(msg)
def mv(A,y):return tuple(sum(a*b for a,b in zip(r,y)) for r in A)
def move(y,w):
 for c in w:y=mv(M[c],y)
 return y
def raw(y):
 a,f,e,p=y;return (a+2*e+p,e,a+f+2*e+2*p,e+p)
def F(y):return sum(a*b for a,b in zip(y,(2,1,3,3)))
def U(y):return raw(y)[2]
def ell(y):return y[0]+y[3]
def h(y):return y[1]-y[2]-y[3]
def c(y):return y[1]+y[2]+y[3]
seed=(0,0,0,1)
def inventory(cutoff):
 out={seed};q=deque([seed])
 while q:
  y=q.popleft()
  for ch in M:
   z=move(y,ch)
   if U(z)<=cutoff and z not in out:out.add(z);q.append(z)
 return out
def limit(N):
 L=0
 while 3**(L+2)<=N*2**(L+1):L+=1
 return L
def suffixes(L):
 out=set()
 for i in range(L+1):
  out.update(('O'*i+'J','Z'*i+'J','O'*i+'JJ'))
  for j in range(L+1):
   out.update(('O'*i+'Z'*j,'Z'*i+'JO'+'Z'*j,'O'*i+'JJO'+'Z'*j))
 return {w for w in out if len(w)<=L}
# Symbolic identities are tested on the four coordinate basis vectors.
for j in range(4):
 y=tuple(int(i==j) for i in range(4));D,E,u,d=raw(y)
 require(move(y,'JZ')==move(y,'OJ'),'rewrite on basis')
 require(F(move(y,'O'))==F(move(y,'J'))==2*F(y)-ell(y),'binary O/J')
 require(F(move(y,'Z'))==2*F(y)+h(y),'binary Z')
 require((h(move(y,'Z')),c(move(y,'Z')))==(h(y),c(y)),'Z sidecars')
 require((h(move(y,'J')),c(move(y,'J')))==(-D,D),'J sidecars')
 require((h(move(y,'O')),c(move(y,'O')))==(y[0]+2*y[1]+y[3],2*u-D),'O sidecars')
 for ch in M:require(2*F(move(y,ch))>=3*F(y),'coefficient growth')
N=1024;states=inventory(N); L=limit(N);stats=[]
for C in (1,2,3,4,5,8,12,20):
 cores={y for y in states if U(y)<=4*C};forms=suffixes(L)
 represented={move(y,w) for y in cores for w in forms}
 qualifying={y for y in states if F(y)<=N and (ell(y)<=C or abs(h(y))<=C)}
 require(qualifying<=represented,'independent normal-form union fails')
 parents={F(y) for y in qualifying if F(y)>4*C}
 bound=6*len(cores)*(L+1)**2
 require(len(qualifying)<=bound,'profile bound')
 # Scalar-first union retains every actual profile in each fiber.
 actual={F(y) for y in states if 4*C<F(y)<=N and any(abs(F(move(y,ch))-2*F(y))<=C for ch in M)}
 require(parents==actual,'arbitrary rebase parent union')
 stats.append({'C':C,'core_profiles':len(cores),'qualifying_profiles':len(qualifying),'rebased_parent_counts':len(parents),'proven_bound':bound})
# Test all recipes, including different recipes with the same numerical profile.
words=0
for n in range(9):
 for tup in product('OZJ',repeat=n):
  w=''.join(tup);y=move(seed,w)
  require(3**(n+1)<=F(y)*2**n,'all-recipe length bound')
  require(ell(y)>0 and ell(y)%2 and h(y)%2 and c(y)>0,'physical parity/positivity')
  words+=1
# Sharp boundary n>4C is intentional; at n=4C equality cannot trigger the strict proof.
# Independently recover the familiar hostile parent after taking all source rebases.
fiber=[y for y in states if F(y)==25]
minoffset=min(abs(F(move(y,ch))-50) for y in fiber for ch in M)
require(minoffset==5,'hostile fiber25')
# Independent automaton count for the added JZ-avoidance core bound.
def fib(n):
 a,b=0,1
 for _ in range(n):a,b=b,a+b
 return a
x,y=1,0; cumulative=0; bigram=[]
for n in range(61):
 total=x+y;cumulative+=total
 require(total==fib(2*n+2),'JZ avoidance Fibonacci identity')
 require(cumulative==fib(2*n+3)-1,'cumulative Fibonacci identity')
 bigram.append(total)
 # x = empty or ending in O/Z; y = ending in J.
 x,y=2*x+y,x+y
for C in range(1,257):
 # Since 4C<=1024 this is a complete bounded core census.
 count=sum(U(y)<=4*C for y in states)
 require(count<=fib(2*limit(8*C)+3)-1,'polynomial core bound')
proof=Path(__file__).with_name('PROOF-frozen.md').read_bytes()
print(json.dumps({'status':'PASS','proof_sha256':hashlib.sha256(proof).hexdigest(),'inventory_cover_bound':N,'profiles':len(states),'length_bound_at_1024':L,'all_recipe_checks':words,'JZ_avoidance_lengths_checked':61,'complete_core_bounds_checked':256,'normal_form_union_checks':stats,'parent25_profiles':len(fiber),'parent25_minimum_absolute_offset':minoffset,'scope':'all-order proof audited independently; numerical replay verifies finite controls only'},sort_keys=True,indent=2))
