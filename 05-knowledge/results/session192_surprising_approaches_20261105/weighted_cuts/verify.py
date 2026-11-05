#!/usr/bin/env python3
"""Exact ear capacities, flows, zonotope strata and convex hostile controls."""
from itertools import combinations, permutations
from collections import deque
from fractions import Fraction
import json, hashlib
from pathlib import Path

def need(test,msg):
 if not test: raise RuntimeError(msg)
def bits(mask):
 while mask:
  b=mask&-mask; yield b.bit_length()-1; mask-=b

def rows_from_code(n,code):
 r=[0]*n
 for k,(a,b) in enumerate(combinations(range(n),2)):
  if code>>k&1:r[a]|=1<<b
  else:r[b]|=1<<a
 return r

def extend(r,S):
 n=len(r);return [v|((0 if S>>i&1 else 1)<<n) for i,v in enumerate(r)]+[S]

def tables(r):
 n=len(r);N=1<<n;full=N-1
 e=[[0]*n for _ in range(N)];s=[[0]*n for _ in range(N)]
 incoming=[sum(1<<j for j in range(n)if r[j]>>i&1)for i in range(n)]
 for i in range(n):e[1<<i][i]=s[1<<i][i]=1
 for S in range(1,N):
  for i in bits(S):
   A=S^(1<<i)
   if A:
    e[S][i]=sum(e[A][j]for j in bits(A&incoming[i]))
    s[S][i]=sum(s[A][j]for j in bits(A&r[i]))
 return s,e

def hp(r):
 if not r:return 1
 return sum(tables(r)[1][-1])

def boundary(r):
 n=len(r);s,e=tables(r);full=(1<<n)-1
 Q=[[0]*n for _ in r]
 for S in range(1,full):
  for i in bits(S):
   if e[S][i]:
    for j in bits(full^S):Q[i][j]+=e[S][i]*s[full^S][j]
 C=[[Q[i][j]+Q[j][i]for j in range(n)]for i in range(n)]
 H=sum(e[full]);field=[]
 for i in range(n):
  a=s[full][i]-e[full][i]+sum(Q[j][i] for j in range(n))
  need(a==sum(C[i][j] for j in bits(r[i])),'outgoing capacity identity')
  need(all(C[i][j]%2==0 for j in range(n)),'even capacities')
  field.append(a-sum(C[i])//2)
 return H,C,field

def reachability(r):
 n=len(r);reach=[]
 for i in range(n):
  S=1<<i
  while True:
   A=S
   for j in bits(S):A|=r[j]
   if A==S:break
   S=A
  reach.append(S)
 return reach

def sccs(r):
 n=len(r);R=reachability(r);remaining=set(range(n));out=[]
 while remaining:
  i=min(remaining);C={j for j in remaining if R[i]>>j&1 and R[j]>>i&1};out.append(C);remaining-=C
 out.sort(key=lambda C:-max(R[i].bit_count() for i in C))
 return out

def strong(r):return bool(r) and all(x==(1<<len(r))-1 for x in reachability(r))
def cut(r,C,S):return sum(C[i][j] for i in bits(S)for j in bits(r[i]&~S))

def flow(r,C,a,b):
 """Edmonds-Karp with integral arc certificates. r fixes directed capacities."""
 n=len(r);cap=[[C[i][j] if r[i]>>j&1 else 0 for j in range(n)] for i in range(n)];res=[v[:]for v in cap];value=0
 while True:
  p=[-1]*n;p[a]=a;q=deque([a])
  while q and p[b]<0:
   u=q.popleft()
   for v in range(n):
    if res[u][v]>0 and p[v]<0:p[v]=u;q.append(v)
  if p[b]<0:break
  bottleneck=10**100;v=b
  while v!=a:bottleneck=min(bottleneck,res[p[v]][v]);v=p[v]
  v=b
  while v!=a:u=p[v];res[u][v]-=bottleneck;res[v][u]+=bottleneck;v=u
  value+=bottleneck
 S=sum(1<<i for i in range(n)if p[i]>=0)
 need(S>>a&1 and not S>>b&1,'flow pins')
 need(cut(r,C,S)==value,'maxflow equals cut')
 F=[[0]*n for _ in r]
 for i in range(n):
  for j in bits(r[i]):F[i][j]=cap[i][j]-res[i][j];need(0<=F[i][j]<=cap[i][j],'capacity')
 for i in range(n):
  net=sum(F[i])-sum(F[j][i] for j in range(n))
  need(net==(value if i==a else -value if i==b else 0),'flow conservation')
 return value,S,F

def minstrong(r,C):
 cc=sccs(r);best=None;records=[]
 for a in sorted(cc[0]):
  for b in sorted(cc[-1]):
   if a==b:continue
   val,S,F=flow(r,C,a,b);records.append({'source':a,'sink':b,'flow':val,'mask':S,'positive_flow_arcs':[[i,j,F[i][j]] for i in range(len(r)) for j in range(len(r)) if F[i][j]]})
   if best is None or val<best[0]:best=(val,S)
 return best,records

def literal_boundary(r):
 n=len(r);Q=[[0]*n for _ in r];H=0
 for p in permutations(range(n)):
  bad=[k for k in range(n-1)if not r[p[k]]>>p[k+1]&1]
  H+=not bad
  for k in range(n-1):
   if not bad or bad==[k]:Q[p[k]][p[k+1]]+=1
 return H,[[Q[i][j]+Q[j][i] for j in range(n)]for i in range(n)]

def main():
 counts={'parents':0,'cuts':0,'flows':0};minimum_capacity=None;zero_capacity=None;triangle_counter=None;minimum_counter=None
 for n in range(2,6):
  for code in range(1<<(n*(n-1)//2)):
   r=rows_from_code(n,code);H,C,h=boundary(r);cc=sccs(r);full=(1<<n)-1
   need(literal_boundary(r)==(H,C),'independent literal exposed pairs')
   zeros=[];valid=[];tri=[]
   for S in range(1<<n):
    d=cut(r,C,S);need(H+d==hp(extend(r,S)),'direct child DP')
    need(2*d==sum(C[i][j] for i in bits(S)for j in bits(full^S))+2*sum(h[i] for i in bits(S)),'zonotope slack')
    if not d:zeros.append(S)
    if strong(extend(r,S)):
     valid.append((d,S));tri.append((sum((r[i]&~S).bit_count()for i in bits(S)),d,S))
    counts['cuts']+=1
   suffixes=[0];suffix=0
   for component in reversed(cc):suffix|=sum(1<<i for i in component);suffixes.append(suffix)
   need(sorted(zeros)==sorted(suffixes),'zero cuts exactly SCC suffixes')
   best,records=minstrong(r,C);need(best[0]==min(valid)[0],'mincut strong minimum');counts['flows']+=len(records)
   if triangle_counter is None:
    for x in tri:
     for y in tri:
      if x[0]<y[0] and x[1]>y[1]:triangle_counter={'rows':r,'fewer_triangles':x,'more_triangles':y}
   mt=min(v[0]for v in tri);md=min(v[1]for v in tri)
   if minimum_counter is None and all(d>md for t,d,S in tri if t==mt):minimum_counter={'rows':r,'triangle_minimum':mt,'increment_minimum':md,'states':tri}
   for i,j in combinations(range(n),2):
    keep=[v for v in range(n) if v not in (i,j)]
    smaller=[sum(1<<b for b,v in enumerate(keep) if r[u]>>v&1) for u in keep]
    need(C[i][j]>=2*hp(smaller),'capacity lower bound')
    minimum_capacity=C[i][j]if minimum_capacity is None else min(minimum_capacity,C[i][j])
    if not C[i][j] and zero_capacity is None:zero_capacity={'rows':r,'pair':[i,j]}
   counts['parents']+=1
 core=[6,8,26,17,3];H,C,h=boundary(core);states=[{'mask':S,'H':H+cut(core,C,S)}for S in range(1,31)]
 lower=next(x for x in states if x['H']==37);upper=next(x for x in states if x['H']==43)
 # Mixture has 39 but cannot be a physical one-birth extension.
 need(Fraction(2,3)*lower['H']+Fraction(1,3)*upper['H']==39,'fractional hostile')
 need(39 not in [x['H']for x in states],'physical 39 gap')
 best,records=minstrong(core,C)
 out={'status':'PASS','counts':counts,'minimum_capacity':minimum_capacity,'zero_capacity':zero_capacity,'triangle_order_hostile':triangle_counter,'different_argmin_hostile':minimum_counter,'core13':{'rows':core,'H':H,'capacities':C,'field':h,'minimum_increment':best,'flow_records':records,'strong_spectrum':sorted(set(x['H']for x in states)),'convex_39':{'lower':lower,'upper':upper,'lower_weight':'2/3','upper_weight':'1/3'}}}
 raw=json.dumps(out,indent=2,sort_keys=True)+'\n';Path(__file__).with_name('receipt.json').write_text(raw)
 print(json.dumps({'status':'PASS',**counts,'minimum_capacity':minimum_capacity,'triangle_order_hostile':triangle_counter,'different_argmin_hostile':minimum_counter,'receipt_sha256':hashlib.sha256(raw.encode()).hexdigest()},sort_keys=True))
if __name__=='__main__':main()
