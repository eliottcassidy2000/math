#!/usr/bin/env python3
"""Exact insertion-circuit checks. Standard library only; assertions not used."""
import itertools,json
from collections import Counter
from explore import tables,kernel

def require(p,msg):
 if not p:raise RuntimeError(msg)

def tournament(n,bits):
 a=[0]*n
 for k,(u,v) in enumerate(itertools.combinations(range(n),2)):
  if bits>>k&1:a[u]|=1<<v
  else:a[v]|=1<<u
 return a

def extend(a,b):
 n=len(a); out=a[:]+[0]
 for u in range(n):
  if b>>u&1:out[n]|=1<<u
  else:out[u]|=1<<n
 return out

def literal(a):
 return sum(all(a[p[i]]>>p[i+1]&1 for i in range(len(a)-1)) for p in itertools.permutations(range(len(a))))

def scc_count(a):
 n=len(a);reach=[a[u]|1<<u for u in range(n)]
 for k in range(n):
  for u in range(n):
   if reach[u]>>k&1:reach[u]|=reach[k]
 return len({tuple(v for v in range(n) if reach[u]>>v&1 and reach[v]>>u&1) for u in range(n)})

def main():
 receipts=[]; total_states=0; brute_states=0
 for n in range(1,6):
  signatures=Counter(); minima=[];checked=0
  for bits in range(1<<(n*(n-1)//2)):
   a=tournament(n,bits);H,K=kernel(a);pairs=list(itertools.combinations(range(n),2))
   w={(u,v):(K[u][v]+K[v][u])//2 for u,v in pairs}
   require(all(K[u][v]+K[v][u]>0 and (K[u][v]+K[v][u])%2==0 for u,v in pairs),'weight positivity/parity')
   vals=[]
   for b in range(1<<n):
    e=extend(a,b);E,_=tables(e);actual=sum(E[-1])
    cut=sum(z for (u,v),z in w.items() if (a[u]>>v&1 and b>>u&1 and not b>>v&1) or (a[v]>>u&1 and b>>v&1 and not b>>u&1))
    require(actual==H+2*cut,'cut identity')
    if n<=4:
     require(actual==literal(e),'literal permutation mismatch');brute_states+=1
    vals.append(actual);checked+=1
   require(sum(v==H for v in vals)==scc_count(a)+1,'SCC zero-cut count')
   for u,v in pairs:
    require(vals[1<<u]+vals[1<<v]-vals[(1<<u)|(1<<v)]-H==2*w[u,v],'tomography')
   for b,val in enumerate(vals):
    pred=H+sum((vals[1<<u]-H) for u in range(n) if b>>u&1)-2*sum(z for (u,v),z in w.items() if b>>u&1 and b>>v&1)
    require(pred==val,'quadratic reconstruction')
   require(2*sum(vals)==(1<<n)*(2*H+sum(w.values())),'uniform mean')
   signatures[tuple(vals)]+=1
  total_states+=checked
  receipts.append({'order':n,'tournaments':1<<(n*(n-1)//2),'extension_states':checked,'distinct_labeled_response_functions':len(signatures),'max_response_fiber':max(signatures.values())})
 transitive=[]
 for n in range(2,11):
  a=[sum(1<<v for v in range(u+1,n)) for u in range(n)];H,K=kernel(a)
  for u,v in itertools.combinations(range(n),2):
   d=v-u;expected=1 if d==1 else 1<<(d-2)
   require(K[u][v]+K[v][u]==2*expected,'transitive weight formula')
  transitive.append(n)
 c=[2,4,1];cr=[4,1,2]
 f=lambda a:[sum(tables(extend(a,b))[0][-1]) for b in range(8)]
 require(f(c)==f(cr),'cycle information-loss witness')
 print(json.dumps({'status':'PASS','exhaustive':receipts,'extension_states':total_states,'independent_permutation_states':brute_states,'transitive_orders_checked':transitive,'opposite_cycles_response':f(c)},indent=2,sort_keys=True))
if __name__=='__main__':main()
