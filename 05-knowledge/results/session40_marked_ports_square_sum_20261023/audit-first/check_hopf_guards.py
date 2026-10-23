#!/usr/bin/env python3
import itertools,math,json
from pathlib import Path

def require(b,msg):
 if not b:raise RuntimeError(msg)
def rank3(rows):
 a=[list(r) for r in rows];r=0
 if not a:return 0
 for j in range(len(a[0])):
  p=next((k for k in range(r,len(a)) if a[k][j]%3),None)
  if p is None:continue
  a[r],a[p]=a[p],a[r];u=pow(a[r][j]%3,-1,3);a[r]=[(x*u)%3 for x in a[r]]
  for k in range(len(a)):
   if k!=r:
    q=a[k][j];a[k]=[(x-q*y)%3 for x,y in zip(a[k],a[r])]
  r+=1
 return r
checks=0
for r in range(2,13):
 N=2**r;support={0:1,N//4:4,N//2:6,3*N//4:4,N:1}
 for j in range(N+1):require(math.comb(N,j)%8==support.get(j,0),'five tap');checks+=1
for r in range(1,9):
 N=2**r
 for j in range(1,N):
  b=math.comb(N,j);v=(b&-b).bit_length()-1
  require(v==r-((j&-j).bit_length()-1),'binomial valuation');checks+=1
allE=list(itertools.combinations(range(4),2));assignments=0;graphrows=[]
for mask in range(64):
 E=[e for i,e in enumerate(allE) if mask>>i&1];index={e:i for i,e in enumerate(E)};m=len(E);T=[]
 for a,b,c in itertools.combinations(range(4),3):
  if all(e in index for e in ((a,b),(b,c),(a,c))):
   row=[0]*m;row[index[(a,b)]]=row[index[(b,c)]]=1;row[index[(a,c)]]=2;T.append(row)
 parent=list(range(4))
 def root(i):
  while parent[i]!=i:i=parent[i]
  return i
 for a,b in E:parent[root(a)]=root(b)
 comp=len({root(i) for i in range(4)});beta=m-4+comp;rank=rank3(T);passing=0
 for v in itertools.product(range(3),repeat=m):
  assignments+=1
  if all(sum(x*y for x,y in zip(t,v))%3==0 for t in T):passing+=1
 require(passing==3**(m-rank),'triangle solution count');checks+=1
 require(0<=rank<=beta,'triangle subspace bound');checks+=1
 graphrows.append({'mask':mask,'beta':beta,'triangle_rank':rank,'hidden_phase_dimension':beta-rank})
# Exact Pauli hostile encoded as 4B=(real,imaginary); orthogonal triples.
u=(0,0,1);v=(1,0,0)
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def cross(a,b):return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
for sign in(-1,1):
 w=(0,sign,0);require([dot(a,b) for a,b in ((u,v),(v,w),(w,u))]==[0,0,0],'pair mags');require((1+dot(u,v)+dot(v,w)+dot(w,u),dot(u,cross(v,w)))==(1,sign),'Bargmann sign');checks+=2
out={'status':'PASS','checks':checks,'phase_assignments':assignments,'mod8_max_depth':4096,'graph_rows':graphrows}
Path(__file__).with_name('hopf-guards-audit-receipt.json').write_text(json.dumps(out,indent=2)+'\n');print({k:v for k,v in out.items() if k!='graph_rows'})
