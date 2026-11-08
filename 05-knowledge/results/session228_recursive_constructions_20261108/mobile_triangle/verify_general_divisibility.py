#!/usr/bin/env python3
import json,hashlib,math
from itertools import combinations,product
from pathlib import Path
from inherited import profile,tournaments,require

def induced(a,vs):return [sum(1<<j for j,w in enumerate(vs) if a[u]>>w&1) for u in vs]
def hq(a):
 n=len(a);dp=[[0]*n for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for mask in range(1,1<<n):
  for i in range(n):
   todo=a[i]&~mask
   while todo:
    bit=todo&-todo;todo-=bit;dp[mask|bit][bit.bit_length()-1]+=dp[mask][i]
 h=list(map(sum,dp));h[0]=1;full=(1<<n)-1
 return h[full],sum(h[s]*h[full^s] for s in range(full+1))
def components(a):
 n=len(a);reach=[row|(1<<i) for i,row in enumerate(a)]
 for k in range(n):
  for i in range(n):
   if reach[i]>>k&1:reach[i]|=reach[k]
 left=set(range(n));comps=[]
 while left:
  u=min(left);vs=[v for v in left if reach[u]>>v&1 and reach[v]>>u&1]
  comps.append(vs);left.difference_update(vs)
 # The count test needs only the components, not their order.
 return comps

def bound(gs):
 G=math.prod(gs);L=1;compgcd=0
 for size in range(min(2,len(gs))+1):
  for S in combinations(range(len(gs)),size):
   removed=math.prod(gs[i] for i in S);L=math.lcm(L,removed);compgcd=math.gcd(compgcd,G//removed)
 require(G%L==0 and G//L==compgcd,'lcm and gcd complement equality')
 return G//L
controls=0;nontrivial=0
for n in range(1,7):
 for a in tournaments(n):
  dgraph=induced(a,list(range(1,n)));gs=[]
  for vs in components(dgraph):
   H,Q=hq(induced(dgraph,vs));gs.append(math.gcd(H,Q))
  K=bound(gs);D=profile(a,0)[0]
  require(D%K==0,'universal SCC bound');controls+=1;nontrivial+=K>1
# Exact arithmetic controls include heterogeneous valuations and the m=0,1,2 boundaries.
synthetic=0
for m in range(6):
 for gs in product((1,3,9,13,39),repeat=m):
  K=bound(gs)
  for p in (3,13):
   vals=[]
   for g in gs:
    e=0
    while g%p==0:g//=p;e+=1
    vals.append(e)
   vals.sort(reverse=True);want=sum(vals[2:]);z=K;got=0
   while z%p==0:z//=p;got+=1
   require(got==want,'valuation formula')
  synthetic+=1
# Larger ordered-triangle examples make the divisor nontrivial beyond the small census.
sharp=[]
for t in range(2,15):
 require(bound([3]*t)==3**(t-2),'triangle specialization')
 sharp.append({'t':t,'K':bound([3]*t)})
require(bound([3,9,13,39])==3,'heterogeneous bound')
files=['GENERAL_DIVISIBILITY.md','verify_general_divisibility.py','inherited.py']
print(json.dumps({'status':'PASS','root0_labelled_orders1_to6':controls,'nontrivial_small_bounds':nontrivial,'synthetic_gcd_lists':synthetic,'triangle_specializations':sharp,'sha256':{f:hashlib.sha256(Path(__file__).with_name(f).read_bytes()).hexdigest() for f in files}},sort_keys=True,indent=2))
