"""Independent permutation-and-cut boundary checks; no imports from producers."""
from itertools import permutations, combinations, product
from collections import Counter
from pathlib import Path
import json

def need(ok, message):
 if not ok:raise ArithmeticError(message)

def topo(A):
 n=len(A); D=[set() for _ in A]
 for u in range(n):
  for v in range(n):
   if A[u][v] and v!=(u+1)%n:D[(u+1)%n].add(v)
 pending=set(range(n)); order=[]
 while pending:
  free=[v for v in pending if not any(v in D[u] for u in pending)]
  if not free:return None
  v=min(free);pending.remove(v);order.append(v)
 return {v:i for i,v in enumerate(order)}

def covers(A,k):
 n=len(A)
 for seq in permutations(range(n)):
  for cuts in combinations(range(1,n),k-1):
   boundaries=(0,)+cuts+(n,)
   if all(A[seq[j]][seq[j+1]] for j in range(n-1) if j+1 not in cuts):
    yield tuple(seq[boundaries[i]] for i in range(k)),tuple(seq[boundaries[i+1]-1] for i in range(k))

rec={'dag_tournaments':{},'all_k_observers':0,'bridge_cases':0}
for n in range(3,6):
 cycle={(i,(i+1)%n) for i in range(n)}
 free=[(i,j) for i in range(n) for j in range(i+1,n) if (i,j) not in cycle and (j,i) not in cycle]
 count=0
 for bits in product((0,1),repeat=len(free)):
  A=[[0]*n for _ in range(n)]
  for u,v in cycle:A[u][v]=1
  for (u,v),b in zip(free,bits):A[u if b else v][v if b else u]=1
  f=topo(A)
  if f is None:continue
  count+=1
  for k in range(1,n+1):
   diag={}
   for starts in permutations(range(n),k):
    ends=tuple(((min(starts,key=lambda t:(t-s)%n or n))-1)%n for s in starts)
    diag[starts]=ends
   need(len(set(diag.values()))==len(diag), (n,k,"matching bijection"))
   M=Counter(covers(A,k))
   need(all(M[s,e]==1 for s,e in diag.items()), (n,k,"unit matching"))
   for (s,e),v in M.items():
    w=sum(f[(x+1)%n] for x in e)-sum(f[x] for x in s)
    need(w>=0 and (w>0 or e==diag[s]),(n,k,s,e,w))
   rec['all_k_observers']+=1
 rec['dag_tournaments'][n]=count

# R={a=0,s=1,t=2}; vary s-t; Q={r=3,...,c=n-1}.
for m in range(2,5):
 n=3+m;r=3;c=n-1
 for st,q in product((0,1),(0,1)):
  forced={(1,0),(0,2),(1,2) if st else (2,1),(c,r)}
  forced|={(u,c) for u in range(4,c)}
  forced|={(u,v) if (u,v)==(q,r) else (v,u) for u in range(3) for v in range(3,n)}
  free=[(u,v) for u in range(3,n) for v in range(u+1,n) if (u,v) not in forced and (v,u) not in forced]
  for bits in product((0,1),repeat=len(free)):
   A=[[0]*n for _ in range(n)]
   for u,v in forced:A[u][v]=1
   for (u,v),b in zip(free,bits):A[u if b else v][v if b else u]=1
   coefficients={(0,r):1,(r,0):-1,(1,r):-1,(r,1):1}
   residual=Counter()
   for starts,ends in covers(A,2):residual[starts]+=coefficients.get(ends,0)
   need(all(v==0 for v in residual.values()),(m,st,q,bits,residual))
   rec['bridge_cases']+=1
rec['status']='PASS'
Path(__file__).with_name('receipt.json').write_text(json.dumps(rec,indent=2)+'\n')
print(json.dumps(rec,indent=2))
