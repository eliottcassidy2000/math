#!/usr/bin/env python3
import itertools,json
from collections import Counter

def check(condition, message):
 if not condition:
  raise RuntimeError(message)

def tournament(n,mask):
 a=[[False]*n for _ in range(n)]
 for k,(i,j) in enumerate(itertools.combinations(range(n),2)):
  if mask>>k&1:a[i][j]=True
  else:a[j][i]=True
 return a

def cycles(a):
 n=len(a);out=[]
 for k in range(3,n+1,2):
  for s in itertools.combinations(range(n),k):
   for tail in itertools.permutations(s[1:]):
    c=(s[0],)+tail
    if all(a[c[i]][c[(i+1)%k]] for i in range(k)):out.append(c)
 return out

def scc(a):
 n=len(a);r=[row[:] for row in a]
 for i in range(n):r[i][i]=True
 for k in range(n):
  for i in range(n):
   for j in range(n):r[i][j]|=r[i][k] and r[k][j]
 return {frozenset(j for j in range(n) if r[i][j] and r[j][i]) for i in range(n)}

def graph(c):
 return [(i,j) for i,j in itertools.combinations(range(len(c)),2) if set(c[i])&set(c[j])]

def components(v,e):
 p=list(range(v))
 def root(i):
  while p[i]!=i:i=p[i]
  return i
 for i,j in e:p[root(i)]=root(j)
 return {frozenset(j for j in range(v) if root(j)==root(i)) for i in range(v)}

def hamiltonian_count(a):
 n=len(a); dp=[[0]*n for _ in range(1<<n)]
 for j in range(n):dp[1<<j][j]=1
 for mask in range(1,1<<n):
  for j in range(n):
   if mask>>j&1:
    prev=mask^(1<<j)
    if prev:dp[mask][j]=sum(dp[prev][i] for i in range(n) if prev>>i&1 and a[i][j])
 return sum(dp[-1])

def long_edge(n):
 a=[[i<j for j in range(n)] for i in range(n)];a[0][-1]=False;a[-1][0]=True
 return a

def ordered_sum(blocks):
 n=sum(map(len,blocks));a=[[i<j for j in range(n)] for i in range(n)];offset=0
 for b in blocks:
  for i,row in enumerate(b):
   for j,x in enumerate(row):a[offset+i][offset+j]=x
  offset+=len(b)
 return a

summary={}
for n in range(1,6):
 hist=Counter();strong_hist=Counter()
 for mask in range(1<<(n*(n-1)//2)):
  a=tournament(n,mask);cs=cycles(a);es=graph(cs);cc=components(len(cs),es);ss=scc(a)
  supports={frozenset(x for i in c for x in cs[i]) for c in cc}
  check(supports=={s for s in ss if len(s)>1}, f"SCC correspondence failed: n={n}, mask={mask}")
  forest=len(es)==len(cs)-len(cc)
  check(forest==all(len(s)<=4 for s in ss), f"Forest classification failed: n={n}, mask={mask}")
  hist[(len(cs),len(es))]+=1
  if len(ss)==1:
   strong_hist[len(cs)]+=1
   if n==4:check(len(cs)==2 and len(es)==1, f"Strong4 shape failed: mask={mask}")
   if n==5:check(len(cs)>=4 and len(es)==len(cs)*(len(cs)-1)//2, f"Strong5 shape failed: mask={mask}")
 summary[str(n)]={'tournaments':sum(hist.values()),'strong_cycle_counts':dict(strong_hist)}
# Explicit reversed-long-edge tournament and its closed form.
for n in range(3,9):
 a=[[i<j for j in range(n)] for i in range(n)];a[0][-1]=False;a[-1][0]=True
 cs=cycles(a);q=2**(n-3)
 check(len(cs)==q and len(graph(cs))==q*(q-1)//2, f"Long-edge clique formula failed: n={n}")
 check(all(0 in c and n-1 in c for c in cs), f"Long-edge common support failed: n={n}")
 check(hamiltonian_count(a)==1+2*q, f"Long-edge Hamiltonian formula failed: n={n}")
summary['long_edge_family_verified_through_n']=8
check(hamiltonian_count(ordered_sum([long_edge(5)]+[long_edge(3)]*3))==243, "Order14 Hamiltonian count failed")
summary['first_nonforest_sparse']={'tournament_order':14,'conflict_vertices':7,'conflict_edges':6,'H':243,'SCC_sizes':[5,3,3,3]}
print(json.dumps(summary,indent=2,sort_keys=True))
