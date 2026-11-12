"""Independent permutation/run audit; imports no producer code."""
from itertools import permutations
from math import comb,factorial
import json
from pathlib import Path

def graph(n,incoming):
 G=[[False]*(n+1) for _ in range(n+1)]
 for i in range(n):
  G[0][i+1]=i not in incoming;G[i+1][0]=i in incoming
  for j in range(n):G[i+1][j+1]=0<(j-i)%n<=(n-1)//2
 return G

def valid(G,p):return all(G[a][b] for a,b in zip(p,p[1:]))
def H(G):return sum(valid(G,p) for p in permutations(range(len(G))))
def cover(G):
 n=len(G);c=[0]*(n+1)
 for p in permutations(range(n)):
  good=sum(G[a][b] for a,b in zip(p,p[1:]));bad=n-1-good
  for extra in range(good+1):c[bad+extra+1]+=comb(good,extra)
 return [x//factorial(i) for i,x in enumerate(c)]
def endpoint(G,ending=False):
 n=len(G);p=[0]*(n+1)
 for tail in permutations(range(1,n)):
  for k in range(1,n+1):
   a=(tail[:k-1]+(0,)) if ending else ((0,)+tail[:k-1])
   if valid(G,a) and valid(G,tail[k-1:]):p[k]+=1
 return p

def triangle_sub(G):
 n=len(G); R=[r+[True,False] for r in G]
 R += [[False]*n+[False,True],[True]*n+[False,False]]
 return R

def subset_count(G):
 n=len(G);dp={((1<<i),i):1 for i in range(n)}
 for mask in range(1,1<<n):
  for last in range(n):
   value=dp.get((mask,last),0)
   if value:
    for nxt in range(n):
     if not(mask>>nxt&1) and G[last][nxt]:
      key=(mask|1<<nxt,nxt);dp[key]=dp.get(key,0)+value
 return sum(dp.get(((1<<n)-1,i),0) for i in range(n))
out={}
for n,As in [(5,({0,1},{0,2})),(7,({0,2},{0,1}))]:
 graphs=[graph(n,A) for A in As];cs=[cover(G) for G in graphs]
 starts=[endpoint(G) for G in graphs];ends=[endpoint(G,True) for G in graphs]
 assert starts[0]==starts[1] and ends[0]==ends[1]
 vals=[subset_count(triangle_sub(G)) for G in graphs]
 assert vals==[3*c[1]+2*c[2] for c in cs]
 out[str(n)]={'covers':cs,'start_by_length':starts[0],'end_by_length':ends[0],'triangle_counts':vals}
assert [b-a for a,b in zip(*out['5']['covers'])]==[0,2,2,0,0,0,0]
assert [b-a for a,b in zip(*out['7']['covers'])]==[0,2,-10,-12,-2,0,0,0,0]
assert out['7']['triangle_counts'][1]-out['7']['triangle_counts'][0]==-14
out['status']='PASS'; print(json.dumps(out,indent=2))
