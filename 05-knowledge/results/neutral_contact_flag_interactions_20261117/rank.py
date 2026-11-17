from collections import Counter
from itertools import combinations_with_replacement,product
from fractions import Fraction as F
from transfer import R
S=list(combinations_with_replacement(range(5),3));idx={a:i for i,a in enumerate(S)};d=len(S)
T=[]
for b in (0,1):
 q=[]
 for p in S:
  v=[0]*d
  for i,c in enumerate(p):
   w=(1+3*b)if not c else ((c-1)&2)+1+b
   t=tuple(sorted(p[:i]+(w,)+p[i+1:]))
   v[idx[t]]+=1
  q.append(v)
 T.append(q)
def move(v,A):return tuple(sum(v[i]*A[i][j]for i in range(d))for j in range(d))
def functional(l,A):return tuple(sum(a*b for a,b in zip(r,l))for r in A)
def basis(seed,func):
 B={};stack=list(seed)
 while stack:
  v=list(map(F,stack.pop()))
  for i,r in sorted(B.items()):
   if v[i]: q=v[i];v=[x-q*y for x,y in zip(v,r)]
  if not any(v):continue
  i=next(i for i,t in enumerate(v)if t);a=v[i]
  v=tuple(t/a for t in v); B[i]=v
  stack.extend(func(v,A) for A in T)
 return list(B.values())
v0=tuple(int(i==0)for i in range(d))
reach=basis([v0],move)
print('reachable',len(reach),'of',d)
for who in (0,1):
 w=[]
 for p in S:
  v=set(__import__('itertools').permutations(p))
  w.append(sum((F(R[(who,0),s]+1,len(v))for s in v),F(0)))
 obs=basis([w],functional)
 G=[]
 for u in reach:
  G.append(tuple(sum(a*b for a,b in zip(u,v))for v in obs))
 B={}
 for row in G:
  w=list(row)
  for i,r in sorted(B.items()):
   if w[i]:x=w[i];w=[a-x*b for a,b in zip(w,r)]
  if any(w):
   i=next(j for j,v in enumerate(w)if v);B[i]=tuple(x/w[i]for x in w)
 print('order',who,'observable',len(obs),'Hankel',len(B))
def center_observable(who):
 return [sum((F(R[(who,0),s]+1,len(set(__import__('itertools').permutations(p))))for s in set(__import__('itertools').permutations(p))),F(0))for p in S]
def project(op):return tuple(sum(a*b for a,b in zip(u,op))for u in reach)
from collections import deque
for order in (0,1,2):
 start=[(center_observable(w),str(w)+':')for w in range(2)] if order==2 else [(center_observable(order),'')]
 B={};queues=deque(start);saved=[]
 def reduce(v,key=True):
  v=list(v);wk=[F(0)]*len(saved)
  for i,(r,c) in sorted(B.items()):
   x=v[i]
   if x:v=[a-x*b for a,b in zip(v,r)];wk=[a+x*b for a,b in zip(wk,c)]
  return v,wk
 while queues:
  v,w=queues.popleft();p,co=reduce(project(v))
  if not any(p):continue
  i=next(j for j,u in enumerate(p)if u);z=p[i]
  for j,(r,t) in list(B.items()):B[j]=(r,t+[F(0)])
  B[i]=(tuple(u/z for u in p),[-u/z for u in co]+[1/z])
  saved.append((v,w))
  queues.extend((functional(v,A),(str(b)+w) if order!=2 else w[0:2]+str(b)+w[2:])for b,A in enumerate(T))
 print('read order',order, 'suffix words', [w for v,w in saved])
 for b,A in enumerate(T):
  M=[]
  for v,w in saved:
   p,c=reduce(project(functional(v,A)))
   assert not any(p);M.append(c)
  print('M',b, [' '.join(map(str,x))for x in M])
def comm(v): return tuple(x-y for x,y in zip(move(move(v,T[0]),T[1]),move(move(v,T[1]),T[0])))
for wh in (0,1):
  obs=basis([center_observable(wh)],functional)
  X=[tuple(sum(a*b for a,b in zip(comm(u),q))for q in obs)for u in reach]
  # plain row rank
  U={}
  for r in X:
   v=list(r)
   for j,q in sorted(U.items()):
    if v[j]:w=v[j];v=[a-w*b for a,b in zip(v,q)]
   if any(v):
    j=next(j for j,x in enumerate(v)if x);U[j]=[x/v[j]for x in v]
  print('commutator section rank',wh,len(U))
