from itertools import combinations, permutations
from check import tiling, chord, sandwich, subset_tables

def fixed_cycles(T):
 n=len(T); count={}
 for k in range(3,n+1,2):
  for U in combinations(range(n),k):
   root=U[0]
   for rest in permutations(U[1:]):
    p=(root,)+rest
    if not all(T[p[i]][p[(i+1)%k]] for i in range(k)): continue
    transformed=[n-1-v for v in p[::-1]]
    j=transformed.index(root) if root in transformed else -1
    if j>=0 and tuple(transformed[j:]+transformed[:j])==p: count[k]=count.get(k,0)+1
 return count
for k in range(1,5):
 m=(k-1)*(k-2)//2 if k>=2 else 0
 data=[]
 for x in range(1<<m):
  A=tiling(k,x); X=chord(A); _,end=subset_tables(X);a=k-1
  Y=chord(sandwich(A)); by_length=fixed_cycles(Y)
  for t in range(1,k+1):
   actual=by_length.get(2*t+1,0)
   target=sum(end[u][a] for u in range(1,1<<k) if u>>a&1 and u.bit_count()==t)
   if actual!=target: raise RuntimeError(('fixed cycle count',k,x,t,actual,target))
  data.append((x,by_length))
 print('k=',k,'mask, fixed-cycles-by-length=',data)
