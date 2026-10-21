"""Independent finite controls for the coprime-factor primeness extension."""
from itertools import product
import json

def skewsets(n):
 for signs in product((1,-1),repeat=(n-1)//2):
  yield {a*s%n for a,s in zip(range(1,(n+1)//2),signs)}
def rows(m,n,A,B):
 pts=list(product(range(m),range(n)));out=[]
 for a,b in pts:
  r=0
  for i,(c,d) in enumerate(pts):
   x,y=(c-a)%m,(d-b)%n
   sign=(y in B) if x==0 else ((x not in A) if y==0 else (x in A))
   if (x or y) and sign:r|=1<<i
  out.append(r)
 return out

def closure(rows,d):
 full=(1<<len(rows))-1;M=1|(1<<d);todo=[d]
 while todo:
  v=todo.pop();D=(rows[0]^rows[v])&~M&full;M|=D
  while D:
   b=D&-D;D-=b;todo.append(b.bit_length()-1)
 return M

def control(m,n):
 results=[]
 for A in skewsets(m):
  for B in skewsets(n):
   r=rows(m,n,A,B)
   proper=sorted({closure(r,d).bit_count() for d in range(1,m*n)}-{m*n})
   results.append(proper)
 return {'m':m,'n':n,'factor_pairs':len(results),'prime_count':sum(not x for x in results),'proper_closure_sizes':sorted({k for x in results for k in x})}
result=[control(m,n) for m,n in [(9,5),(3,5),(3,3)]]
if result[0]['prime_count']!=64 or result[1]['prime_count']!=8 or result[2]['prime_count']!=0:raise RuntimeError('unexpected extension control')
r=rows(3,3,{1},{1});diag={0,4,8}
if any(len({bool(r[x]>>y&1) for y in diag})!=1 for x in range(9) if x not in diag):raise RuntimeError('diagonal is not a module')
print(json.dumps({'finite_controls':result,'noncoprime_hostile_diagonal':sorted(diag)},indent=2))
