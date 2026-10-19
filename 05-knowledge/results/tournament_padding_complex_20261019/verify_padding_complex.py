from itertools import permutations
from functools import lru_cache
import json
checks=0
def ck(x,m):
 global checks
 checks+=1
 if not x:raise RuntimeError(m)
def pairs(n):return [(i,j)for i in range(n)for j in range(i+1,n)]
@lru_cache(None)
def edge(n,mask,i,j):
 if i>j:return 1-edge(n,mask,j,i)
 return (mask>>pairs(n).index((i,j)))&1
@lru_cache(None)
def relabel(n,mask,order):
 return sum(edge(n,mask,order[i],order[j])<<k for k,(i,j)in enumerate(pairs(len(order))))
@lru_cache(None)
def canon(n,mask):
 return min(relabel(n,mask,p) for p in permutations(range(n)))
classes=[]
for n in range(6):classes.append(sorted({canon(n,m) for m in range(1<<(n*(n-1)//2))}))
index=[{m:i for i,m in enumerate(cs)}for cs in classes]
d=[];h=[[]]
for n in range(5):
 rows=[]
 for T in classes[n+1]:
  row=0
  for v in range(n+1):
   sub=relabel(n+1,T,tuple(i for i in range(n+1)if i!=v))
   row^=1<<index[n][canon(n,sub)]
  rows.append(row)
 d.append(rows)
for n in range(1,6):
 rows=[]
 for T in classes[n-1]:
  # New source is last vertex; old->source is false.
  aug=sum(edge(n-1,T,i,j)<<k for k,(i,j)in enumerate(pairs(n))if j<n-1)
  rows.append(1<<index[n][canon(n,aug)])
 h.append(rows)
def multiply(A,B):
 out=[]
 for row in A:
  z=0
  for j in range(len(B)):
   if row>>j&1:z^=B[j]
  out.append(z)
 return out
def rank(rows):
 basis={}
 for x in rows:
  while x:
   p=x.bit_length()-1
   if p not in basis:basis[p]=x;break
   x^=basis[p]
 return len(basis)
for n in range(4):ck(all(x==0 for x in multiply(d[n+1],d[n])),('d2',n))
for n in range(5):
 left=multiply(d[n-1],h[n]) if n else [0]
 right=multiply(h[n+1],d[n])
 ck([x^y for x,y in zip(left,right)]==[1<<i for i in range(len(classes[n]))],('contraction',n))
ranks=[rank(x)for x in d]
for n,r in enumerate(ranks):
 ck(r==sum((-1)**(n-j)*len(classes[j])for j in range(n+1)),('rank',n))
for n,cs in enumerate(classes):
 for T in cs:
  aut=sum(relabel(n,T,p)==T for p in permutations(range(n)))
  ck(aut%2==1,('odd automorphism',n,T))
print(json.dumps({'status':'PASS','checks':checks,'class_counts':[len(x)for x in classes],'direct_ranks':ranks,'predicted_rank_d5':9,'scope':'Exact class-space matrices through five-vertex hosts'},indent=2))
