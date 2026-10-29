from itertools import combinations,product
import json

def require(x,m):
 if not x:raise RuntimeError(m)
def rank(M,p=101):
 A=[list(row) for row in M];r=0
 for j in range(len(A[0]) if A else 0):
  k=next((k for k in range(r,len(A)) if A[k][j]%p),None)
  if k is None:continue
  A[r],A[k]=A[k],A[r];q=pow(A[r][j]%p,-1,p);A[r]=[(v*q)%p for v in A[r]]
  for k in range(r+1,len(A)):
   z=A[k][j]%p
   if z:A[k]=[(v-z*w)%p for v,w in zip(A[k],A[r])]
  r+=1
 return r
cases=0
for m in range(1,5):
 es=list(combinations(range(m),2)); states=list(product((-1,0,1),repeat=m));idx={s:i for i,s in enumerate(states)}
 for mask in range(1<<len(es)):
  edges={e for k,e in enumerate(es) if mask>>k&1};L=[]
  for i in range(m):
   for others in product((-1,0,1),repeat=m-1):
    s=list(others);s.insert(i,0);t=s.copy();t[i]=1
    parity=sum(s[j] for j in range(i) if s[j]>=0 and (j,i) in edges)%2
    row=[0]*len(states);row[idx[tuple(s)]]=1;row[idx[tuple(t)]]=(-1)**parity;L.append(row)
  indep=[]
  for sm in range(1<<m):
   if not any(sm>>i&1 and sm>>j&1 for i,j in edges):indep.append(sm)
  # Explicit independent kernel vectors, one disjoint support per independent set.
  for sm in indep:
   f=[(-1)**sum(v==1 for v in s) if sum((v>=0)<<i for i,v in enumerate(s))==sm else 0 for s in states]
   require(all(sum(a*b for a,b in zip(row,f))==0 for row in L),'integer kernel')
  r=rank(L);require(r==3**m-len(indep),'modular lower bound sharp')
  cases+=1
print(json.dumps({'all_graphs_orders_1_to_4':cases,'rank_method':'F101 boundary elimination plus independent integer kernel witnesses','status':'pass'}))
