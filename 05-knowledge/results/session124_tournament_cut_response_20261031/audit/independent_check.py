"""Independent literal two-path-cover and all-cut audit; explicit exceptions."""
from pathlib import Path
from itertools import permutations,product
import json,hashlib,sys
from math import gcd
SOURCE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(SOURCE))
from cut_response import compile_response

def graph(n,code):
 edges=[[False]*n for _ in range(n)];j=0
 for u in range(n):
  for v in range(u+1,n):
   edges[u][v]=bool(code&(1<<j));edges[v][u]=not edges[u][v];j+=1
 return edges

def hamiltonian(A):
 n=len(A)
 if not n:return 1
 d={(1<<u,u):1 for u in range(n)}
 for mask in range(1,1<<n):
  for last in range(n):
   x=d.get((mask,last),0)
   if not x:continue
   for nxt in range(n):
    if not mask>>nxt&1 and A[last][nxt]:
     key=(mask|(1<<nxt),nxt);d[key]=d.get(key,0)+x
 return sum(d.get(((1<<n)-1,v),0) for v in range(n))

def literal(A):
 n=len(A);H=0;first=[0]*n;last=[0]*n;W=[[0]*n for _ in range(n)]
 for p in permutations(range(n)):
  bad=[j for j in range(n-1) if not A[p[j]][p[j+1]]]
  if not bad:
   H+=1
   if p:first[p[0]]+=1;last[p[-1]]+=1
  for split in range(1,n):
   if all(j==split-1 for j in bad):W[p[split-1]][p[split]]+=1
 return H,first,last,W

def child(A,mask):
 n=len(A)
 return [row+[not bool(mask>>u&1)] for u,row in enumerate(A)]+[[bool(mask>>u&1) for u in range(n)]+[False]]

def suffixes(A):
 n=len(A);R=[[i==j or A[i][j] for j in range(n)] for i in range(n)]
 for k in range(n):
  for i in range(n):
   for j in range(n):R[i][j]=R[i][j] or (R[i][k] and R[k][j])
 comps=[];unused=set(range(n))
 while unused:
  u=min(unused);C={v for v in unused if R[u][v] and R[v][u]};comps.append(C);unused-=C
 # Sort components by number of other components dominated; source first.
 original=list(comps)
 comps.sort(key=lambda C:-sum(A[next(iter(C))][next(iter(D))] for D in original if D!=C))
 result={0};mask=0
 for C in reversed(comps):mask|=sum(1<<u for u in C);result.add(mask)
 return result

checks=dict(parents=0,cut_dp=0,literal_cover_tables=0,pair_deletion_bounds=0,endpoint_balances=0,forced_minima=0,mixed_differences=0,scc_suffixes=0,positive_minima=0,gain_content=0)
digest=hashlib.sha256()
for n in range(6):
 for code in range(1<<(n*(n-1)//2)):
  A=graph(n,code);H,starts,ends,W=literal(A)
  out=tuple(sum(1<<v for v in range(n) if A[u][v]) for u in range(n));r=compile_response(out)
  if (r.h,r.starts,r.ends,r.covers)!=(H,starts,ends,W):raise RuntimeError('literal table')
  checks['parents']+=1;checks['literal_cover_tables']+=1
  b=[[0]*n for _ in range(n)]
  for u in range(n):
   balance=sum(W[u][v] for v in range(n) if A[u][v])-sum(W[v][u] for v in range(n) if A[v][u])
   if balance!=starts[u]-ends[u]:raise RuntimeError('balance')
   checks['endpoint_balances']+=1
   for v in range(u+1,n):
    twice=W[u][v]+W[v][u]
    remain=[i for i in range(n) if i not in (u,v)]
    old=hamiltonian([[A[i][j] for j in remain] for i in remain])
    if twice%2 or twice<2*old:raise RuntimeError('positive integral weight')
    b[u][v]=b[v][u]=twice//2;checks['pair_deletion_bounds']+=1
  if r.weights!=b:raise RuntimeError('weight table')
  values=[];mincuts=set()
  for mask in range(1<<n):
   capacity=sum(b[u][v] for u in range(n) for v in range(n) if mask>>u&1 and not mask>>v&1 and A[u][v])
   actual=hamiltonian(child(A,mask))
   if actual!=H+2*capacity or r.value(mask)!=actual:raise RuntimeError(('cut',n,code,mask))
   values.append(actual);checks['cut_dp']+=1
   if actual==H:mincuts.add(mask)
  if r.all_values()!=values:raise RuntimeError('all-values recurrence')
  positive=[x for x in values if x>H];answer=r.minimum_positive()
  if positive:
   if answer is None or answer[0]!=min(positive) or values[answer[1]]!=answer[0]:raise RuntimeError('minimum positive')
  elif answer is not None:raise RuntimeError('spurious positive gain')
  checks['positive_minima']+=1
  divisor=0
  for value in values:divisor=gcd(divisor,value-H)
  if r.gain_content()!=divisor:raise RuntimeError('gain gcd')
  checks['gain_content']+=1
  if mincuts!=suffixes(A):raise RuntimeError('SCC suffix minimizers')
  checks['scc_suffixes']+=1
  for u in range(n):
   for v in range(u+1,n):
    for mask in range(1<<n):
     if mask&((1<<u)|(1<<v)):continue
     if values[mask|1<<u]+values[mask|1<<v]-values[mask]-values[mask|1<<u|1<<v]!=2*b[u][v]:raise RuntimeError('mixed difference')
     checks['mixed_differences']+=1
  constraints=product(range(3),repeat=n) if n<=4 else [tuple(u%3 for u in range(n))]
  for marks in constraints:
   inside=sum(1<<u for u in range(n) if marks[u]==1);outside=sum(1<<u for u in range(n) if marks[u]==2)
   true=min(values[x] for x in range(1<<n) if x&inside==inside and not x&outside)
   value,x=r.minimum(inside,outside)
   if value!=true or values[x]!=true or x&inside!=inside or x&outside:raise RuntimeError('mincut')
   checks['forced_minima']+=1
  digest.update(json.dumps([n,code,H,W,values],separators=(',',':')).encode()+b'\n')
# Independent verification of two free-order block insertions via all comparison words.
blockchecks=0
for n in range(7):
 for word in product(range(4),repeat=n):
  count=0
  for position in range(n+1):
   # Low bit: old vertex beats u; high bit: old vertex beats v.
   left=(1,1) if position==0 else (word[position-1]&1,word[position-1]>>1)
   right=(0,0) if position==n else (word[position]&1,word[position]>>1)
   count+=left[0] and not right[1]
   count+=left[1] and not right[0]
  if count<2:raise RuntimeError('block lower bound')
  blockchecks+=1
print(json.dumps(dict(status='PASS',checks=checks,block_comparison_words=blockchecks,data_sha256=digest.hexdigest(),source_sha256=hashlib.sha256((SOURCE/'cut_response.py').read_bytes()).hexdigest()),sort_keys=True))
