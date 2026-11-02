#!/usr/bin/env python3
"""Independent unordered-cover / ordered-run coefficient audit; stdlib only."""
from functools import lru_cache
from itertools import combinations
from math import comb,factorial
import json

def need(c,m):
 if not c:raise RuntimeError(m)
def decode(n,code):
 a=[[0]*n for _ in range(n)]
 for k,(i,j) in enumerate(combinations(range(n),2)):
  a[i][j]=(code>>k)&1;a[j][i]=1-a[i][j]
 return tuple(map(tuple,a))
def H(a):
 n=len(a)
 if not n:return 1
 @lru_cache(None)
 def walk(last,mask):
  if not mask:return 1
  return sum(walk(v,mask^(1<<v)) for v in range(n) if mask>>v&1 and a[last][v])
 return sum(walk(v,((1<<n)-1)^(1<<v)) for v in range(n))
def induced(a,mask):
 vs=[v for v in range(len(a)) if mask>>v&1]
 return tuple(tuple(a[i][j] for j in vs) for i in vs)
def covers(a):
 n=len(a);full=(1<<n)-1
 hs=[H(induced(a,m)) for m in range(1<<n)]
 @lru_cache(None)
 def p(mask,k):
  if not mask:return int(k==0)
  if k==0:return 0
  least=mask&-mask;ans=0;part=mask
  while part:
   if part&least:ans+=hs[part]*p(mask^part,k-1)
   part=(part-1)&mask
  return ans
 return [p(full,k) for k in range(n+1)]
def substitute(q,v,t):
 external=[i for i in range(len(q)) if i!=v];n=len(t)
 a=[[0]*(n+len(external)) for _ in range(n+len(external))]
 for i in range(len(a)):
  for j in range(len(a)):
   a[i][j]=t[i][j] if i<n and j<n else q[v if i<n else external[i-n]][v if j<n else external[j-n]]
 return tuple(map(tuple,a))
def cyclic(q,v):return any(q[v][a] and q[a][b] and q[b][v] for a in range(len(q)) for b in range(len(q)))
def main():
 source_count=0;cover_inequalities=0;triangle_checks=0;equalities=[]
 triangle=decode(3,2)
 small=[]
 for n in range(1,6):
  for code in range(1<<(n*(n-1)//2)):
   t=decode(n,code);h=H(t);ps=covers(t);source_count+=1
   if n<=3:small.append((t,h,ps))
   for k in range(1,n+1):
    need(factorial(k)*ps[k]>=h*comb(n-1,k-1),'unordered path-cut lower bound')
    cover_inequalities+=1
   value=H(substitute(triangle,0,t));p2=ps[2] if n>=2 else 0
   need(value==3*h+2*p2,'triangle coefficient convention')
   bound=(n+2)*h+int(n%2==0)
   need(value>=bound,'cyclic-port amplification')
   if value==bound:equalities.append([n,code,h,value])
   triangle_checks+=1
 contexts=cyclic_contexts=kernel_checks=0
 for m in range(1,5):
  for code in range(1<<(m*(m-1)//2)):
   q=decode(m,code)
   for v in range(m):
    contexts+=1;ws=[]
    for k in range(1,m+1):
     independent=tuple((0,)*k for _ in range(k))
     labelled=H(substitute(q,v,independent))
     need(labelled%factorial(k)==0,'twin-label factorial')
     ws.append(labelled//factorial(k))
    need(ws[0]==H(q),'W1 convention')
    if cyclic(q,v):
     cyclic_contexts+=1
     need(ws[0]>=3 and ws[1]>=1,'cyclic-port coefficients')
    for t,h,ps in small:
     need(H(substitute(q,v,t))==sum(ws[k-1]*factorial(k)*ps[k] for k in range(1,min(len(t),m)+1)),'exact run-cover kernel')
     kernel_checks+=1
 print(json.dumps({'status':'PASS','source_tournaments_n_le5':source_count,'cover_inequalities':cover_inequalities,
  'triangle_substitutions':triangle_checks,'marked_quotient_contexts_n_le4':contexts,'cyclic_contexts':cyclic_contexts,
  'exact_kernel_controls':kernel_checks,'sharp_lower_bound_examples':equalities[:8],
  'convention':'p_k counts unordered nonempty directed path covers; W_k counts words with k indistinguishable marked runs; coefficient is k!*W_k.',
  'scope':'Finite controls only; path-cut injection and triangle-twin insertion give the all-order proof.'},indent=2,sort_keys=True))
if __name__=='__main__':main()
