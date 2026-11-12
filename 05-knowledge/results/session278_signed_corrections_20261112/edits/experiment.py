from itertools import permutations,combinations
from math import comb
import random,json
R=random.Random(278)
def path_count(A):
 n=len(A);dp={(1<<i,i):1 for i in range(n)}
 for mask in range(1,1<<n):
  for i in range(n):
   v=dp.get((mask,i),0)
   if v:
    for j in range(n):
     if not mask>>j&1 and A[i][j]:dp[(mask|1<<j,j)]=dp.get((mask|1<<j,j),0)+v
 return sum(dp.get(((1<<n)-1,i),0) for i in range(n))
def data(A,E):
 n=len(A);U=[r[:] for r in A];old=set(E);new={(v,u) for u,v in E}
 for u,v in E:U[v][u]=1
 hist={}
 for p in permutations(range(n)):
  arcs=list(zip(p,p[1:]))
  if all(U[u][v] for u,v in arcs):
   a=sum(e in old for e in arcs);b=sum(e in new for e in arcs)
   hist[a,b]=hist.get((a,b),0)+1
 B=[r[:] for r in A]
 for u,v in E:B[u][v]=0;B[v][u]=1
 h0=sum(v for (a,b),v in hist.items() if b==0);h1=sum(v for (a,b),v in hist.items() if a==0)
 assert (h0,h1)==(path_count(A),path_count(B))
 rows=[]
 for r in range(n):
  def approx(axis):return sum(v*sum((-1)**j*comb((a,b)[axis],j) for j in range(min(r,(a,b)[axis])+1)) for (a,b),v in hist.items())
  S0=approx(1);S1=approx(0)
  error=sum(v*((comb(a-1,r) if a>r else 0)-(comb(b-1,r) if b>r else 0)) for (a,b),v in hist.items())*(-1)**r
  assert S1-S0-(h1-h0)==error
  rows.append([r,S1-S0,error])
 return {'n':n,'rows':[''.join(map(str,row)) for row in A],'edits':E,'counts':[h0,h1],'joint_histogram':[[a,b,v] for (a,b),v in sorted(hist.items())],'truncations_r_delta_error':rows}
checks=0;examples=[]
for n in range(3,8):
 for t in range(60 if n<7 else 15):
  A=[[0]*n for _ in range(n)]
  for u,v in combinations(range(n),2):
   if R.randrange(2):A[u][v]=1
   else:A[v][u]=1
  arcs=[(u,v) for u in range(n) for v in range(n) if A[u][v]]
  E=R.sample(arcs,min(4,len(arcs)));d=data(A,E);checks+=1
  if len(examples)<8 and abs(d['truncations_r_delta_error'][1][2])>=4:examples.append(d)
print(json.dumps({'seed':278,'independent_DP_controls':checks,'nontrivial_correction_examples':examples},indent=2))
