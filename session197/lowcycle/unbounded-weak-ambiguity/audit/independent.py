from functools import lru_cache
from itertools import permutations
from math import prod,factorial
from pathlib import Path
import json
P=Path(__file__).parent

def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def make(js):
 rows=[0];root=0
 for j in js:
  n=len(rows);new=list(range(n,n+j+1));rows += [0]*(j+1)
  for a in range(n):
   for b in new:
    if b==n and a==root:rows[b]|=1<<a
    else:rows[a]|=1<<b
  # The gadget root is last in its transitive order.
  order=list(range(n+1,n+j+1))+[n]
  for i,a in enumerate(order):
   for b in order[i+1:]:rows[a]|=1<<b
  root=n
 return rows,root

def inventory(rows,root):
 n=len(rows);full=(1<<n)-1
 rev=[sum(1<<j for j in range(n)if rows[j]>>i&1)for i in range(n)]
 def counter(arcs):
  @lru_cache(None)
  def f(M,v):
   if M==1<<v:return 1
   N=M^(1<<v)
   return sum(f(N,u)for u in range(n)if N>>u&1 and arcs[u]>>v&1)
  def h(M):return sum(f(M,v)for v in range(n)if M>>v&1)if M else 1
  return f,h
 end,H=counter(rows);start,_=counter(rev)
 R=sum(H(full^M)*start(M,root)for M in range(1<<n)if M>>root&1)
 return H(full),end(full,root),R,start(full,root),H(full^(1<<root))
source=json.loads((P/'source'/'weak_ambiguity_certificate.json').read_text())
checks=literal=0
for record in source['physical']['witnesses']:
 m=record['m'];last=record['last'];js=[j for j in range(1,m+1)if j!=last]+[last]
 A,r=make(js);n=len(A)
 for i in range(n):A[i]|=1<<n
 A.append(0)
 need(A==record['rows']and r==record['root'],'independent adjacency')
 h=prod(2**j+1 for j in range(1,m+1));d=h//(2**last+1)
 need(inventory(A,r)==(h,0,2*h,1,d),'independent predecessor DP')
 if m<=2:
  total=start=0;deletion=0
  for p in permutations(range(len(A))):
   if all(A[u]>>v&1 for u,v in zip(p,p[1:])):total+=1;start+=p[0]==r
  for p in permutations([v for v in range(len(A))if v!=r]):
   deletion+=all(A[u]>>v&1 for u,v in zip(p,p[1:]))
  need((total,start,deletion)==(h,1,d),'literal graph controls');literal+=1
 checks+=1
orders=0
for m in range(1,8):
 Pm=prod(2**j+1 for j in range(1,m+1));seen=set()
 for js in permutations(range(1,m+1)):
  h,e,R,d,n=1,1,1,1,1
  for j in js:h,e,R,d,n=h+2**j*R,h,h+2**j*R,h,n+j+1
  need(h==R==Pm and e==d==Pm//(2**js[-1]+1)and n+1==2+m*(m+3)//2,'all reorderings')
  seen.add(d);orders+=1
 need(len(seen)==m,'all-ordering deletion image')
print(json.dumps(dict(status='PASS',independent_graphs=checks,literal_controls=literal,all_gadget_orderings=orders,max_reordering_m=7),sort_keys=True,indent=2))
