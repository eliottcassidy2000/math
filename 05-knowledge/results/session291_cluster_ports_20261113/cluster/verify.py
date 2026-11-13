#!/usr/bin/env python3
"""Exact independent graph-path and squarefree pivot tests; no external packages."""
from itertools import permutations, combinations
from collections import Counter
import json, random

ONE={0:1}
def add(a,b):
 d=a.copy()
 for m,c in b.items(): d[m]=d.get(m,0)+c
 return {m:c for m,c in d.items() if c}
def neg(a): return {m:-c for m,c in a.items()}
def sub(a,b): return add(a,neg(b))
def mul(a,b):
 d={}
 for m,c in a.items():
  for n,e in b.items():
   if not m&n: d[m|n]=d.get(m|n,0)+c*e
 return {m:c for m,c in d.items() if c}
def inv(a,n):
 if a.get(0)!=1: raise ValueError('unit must have constant one')
 z=sub(ONE,a); term=ONE; result=ONE
 for _ in range(n): term=mul(term,z); result=add(result,term)
 if mul(a,result)!=ONE: raise ValueError('bad inverse')
 return result

def adj(n,bits):
 a=[[0]*n for _ in range(n)]
 for t,(i,j) in enumerate(combinations(range(n),2)):
  if bits>>t&1:a[i][j]=1
  else:a[j][i]=1
 return a

def mutate(a,k):
 n=len(a); b=[r[:] for r in a]
 for i in range(n):
  for j in range(n):
   if i!=j and i!=k and j!=k:b[i][j]+=a[i][k]*a[k][j]
 for i in range(n):b[i][k],b[k][i]=b[k][i],b[i][k]
 for i,j in combinations(range(n),2):
  m=min(b[i][j],b[j][i]);b[i][j]-=m;b[j][i]-=m
 return b

def tournament(a):
 n=len(a)
 return all(a[i][i]==0 for i in range(n)) and all(a[i][j]+a[j][i]==1 and a[i][j] in (0,1) and a[j][i] in (0,1) for i,j in combinations(range(n),2))

def paths(a):
 n=len(a); p=[[{} for _ in range(n)] for _ in range(n)]
 def walk(s,t,m):
  p[s][t][m]=p[s][t].get(m,0)+1
  for v in range(n):
   if not m>>v&1 and a[t][v]:walk(s,v,m|1<<v)
 for s in range(n):walk(s,s,1<<s)
 return p

def hp(a):
 p=paths(a);m=(1<<len(a))-1
 return sum(x.get(m,0) for r in p for x in r)

def eliminate(a,order,check_exchange=False):
 n=len(a);N=3*n;q=[[{} for _ in range(N)] for _ in range(N)]
 for i in range(n):
  q[n+i][i]={1<<i:1};q[i][2*n+i]=ONE
  for j in range(n):
   if a[i][j]:q[i][j]={1<<j:1}
 active=list(range(N));exchanges=0
 for k in order:
  rest=[i for i in active if i!=k]
  if check_exchange:
   others=[i for i in rest if i<n]
   if others:
    u=others[0]
    aa=sub(ONE,q[k][k]);dd=sub(ONE,q[u][u]);bc=mul(q[k][u],q[u][k]);D=sub(mul(aa,dd),bc)
    if mul(aa,dd)!=add(D,bc):raise ValueError('exchange failed')
    di=inv(D,n)
    # Compare the displayed two-by-two inverse to two sequential pivots,
    # at every surviving pair, not just its defining determinant identity.
    sa=inv(aa,n)
    first=lambda i,j:add(q[i][j],mul(mul(q[i][k],sa),q[k][j]))
    su=inv(sub(ONE,first(u,u)),n)
    for i in rest:
     if i==u:continue
     for j in rest:
      if j==u:continue
      num=add(add(mul(mul(q[i][k],dd),q[k][j]),mul(mul(q[i][k],q[k][u]),q[u][j])),add(mul(mul(q[i][u],q[u][k]),q[k][j]),mul(mul(q[i][u],aa),q[u][j])))
      block=add(q[i][j],mul(num,di))
      sequential=add(first(i,j),mul(mul(first(i,u),su),first(u,j)))
      if block!=sequential:raise ValueError('two-pivot numerator mismatch')
    exchanges+=1
  star=inv(sub(ONE,q[k][k]),n)
  old=[[e.copy() for e in row] for row in q]
  for i in rest:
   for j in rest:q[i][j]=add(old[i][j],mul(mul(old[i][k],star),old[k][j]))
  active=rest
 p=[[q[n+i][2*n+j] for j in range(n)] for i in range(n)]
 if any(c<0 for r in p for e in r for c in e.values()):raise ValueError('positivity')
 return p,exchanges

def scc(a):
 n=len(a);r=[row[:] for row in a]
 for i in range(n):r[i][i]=1
 for k in range(n):
  for i in range(n):
   for j in range(n):r[i][j]|=r[i][k]&r[k][j]
 groups=[];todo=set(range(n))
 while todo:
  i=min(todo);g=tuple(j for j in sorted(todo) if r[i][j] and r[j][i]);groups.append(g);todo-=set(g)
 groups.sort(key=lambda g:-sum(a[g[0]][h[0]] for h in groups if h!=g))
 return groups

def key(a):return tuple(tuple(r) for r in a)
def orbit(a):
 seen={key(a)};todo=[a]
 while todo:
  b=todo.pop();n=len(b)
  for k in range(n):
   if sum(b[k]) in (0,n-1):
    c=mutate(b,k);ck=key(c)
    if ck not in seen:seen.add(ck);todo.append(c)
 return seen

def skew_mutate(b,k):
 n=len(b)
 return tuple(tuple(-b[i][j] if i==k or j==k else b[i][j]+max(b[i][k],0)*max(b[k][j],0)-max(-b[i][k],0)*max(-b[k][j],0) for j in range(n)) for i in range(n))

def pentagon(a,u,v):
 n=len(a);b=tuple(tuple(a[i][j]-a[j][i] for j in range(n)) for i in range(n));old=b
 for k in (u,v,u,v,u): b=skew_mutate(b,k)
 p=list(range(n));p[u],p[v]=p[v],p[u]
 if b!=tuple(tuple(old[p[i]][p[j]] for j in range(n)) for i in range(n)):raise ValueError('pentagon')

def main():
 counts={};checks=0;components=0;pentagons=0
 for n in range(1,6):
  valid=0;orbit_seen=set()
  for bits in range(1<<(n*(n-1)//2)):
   a=adj(n,bits);h=hp(a)
   for u,v in combinations(range(n),2):pentagon(a,u,v);pentagons+=1
   for k in range(n):
    b=mutate(a,k);v=tournament(b)
    if v!=(sum(a[k]) in (0,n-1)):raise ValueError('FZ gate')
    checks+=1
    if v:
     valid+=1
     if hp(b)!=h:raise ValueError('H changed')
   if key(a) not in orbit_seen:
    oo=orbit(a);orbit_seen |=oo;cc=scc(a)
    nons=[i for i,g in enumerate(cc) if len(g)>1]
    expected=n if not nons else 1+nons[0]+len(cc)-1-nons[-1]
    if len(oo)!=expected:raise ValueError(('component',n,bits,len(oo),expected))
    components+=1
  counts[str(n)]=valid
 rng=random.Random(291);graphs=0;orders=0;exchanges=0
 for n in range(1,7):
  vals=range(1<<(n*(n-1)//2)) if n<=4 else sorted(set([0,(1<<(n*(n-1)//2))-1]+[rng.randrange(1<<(n*(n-1)//2)) for _ in range(20)]))
  for bits in vals:
   a=adj(n,bits);direct=paths(a)
   oo=list(permutations(range(n))) if n<=3 else [tuple(range(n)),tuple(reversed(range(n)))]+[tuple(rng.sample(range(n),n)) for _ in range(3)]
   for order in oo:
    p,e=eliminate(a,order,True)
    if p!=direct:raise ValueError(('path atlas',n,bits,order))
    orders+=1;exchanges+=e
   graphs+=1
 # Honest nonpositive principal determinant: cyclic triangle determinant=1-x0x1x2.
 cyclic=[[0,1,0],[0,0,1],[1,0,0]]
 if hp(cyclic)!=3:raise ValueError('cycle control')
 print(json.dumps({'status':'PASS','pentagon_tests':pentagons,'fz_vertex_tests':checks,'fz_valid_counts':counts,'labelled_components':components,'atlas_graphs':graphs,'elimination_orders':orders,'two_pivot_exchange_checks':exchanges,'cycle3_determinant':{'constant':1,'x0x1x2':-1},'scope':'FZ all labelled n<=5; atlas all labelled n<=4 and seeded samples n5,n6'},sort_keys=True,indent=2))
if __name__=='__main__':main()
