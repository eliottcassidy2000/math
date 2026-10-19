#!/usr/bin/env python3
"""Exact controls for monotone bins, arbitrary chain-surplus gadgets, and endpoint obstruction."""
from collections import defaultdict,Counter
from itertools import combinations,product
from functools import lru_cache
import hashlib,json

def require(ok,msg):
 if not ok:raise RuntimeError(msg)
def close(p):
 p=list(p)
 for k in range(len(p)):
  for v in range(len(p)):
   if p[v]>>k&1:p[v]|=p[k]
 return tuple(p)
def posets(n):
 edges=list(combinations(range(n),2));seen=set()
 for mask in range(1<<len(edges)):
  p=[0]*n
  for k,(u,v) in enumerate(edges):
   if mask>>k&1:p[v]|=1<<u
  q=close(p)
  if q not in seen:seen.add(q);yield q
def extensions(p):
 out=[];n=len(p)
 def rec(mask,w):
  if len(w)==n:out.append(tuple(w));return
  for v in range(n):
   if not mask>>v&1 and not p[v]&~mask:rec(mask|1<<v,w+[v])
 rec(0,[]);return out

def strategy(q):
 E=sum(q);require(E>0 and q[0]==q[1] and all(q[i]>=q[i+1]>=0 for i in range(len(q)-1)),'Invalid monotone-bin law')
 if 2*q[0]==E:return (0,1)
 prefix=0
 for ell,x in enumerate(q):
  if 2*(prefix+x)>=E:break
  prefix+=x
 require(1<=ell<len(q)-1 and x<=prefix,'Adjacent-cut boundary failed')
 return (ell-1,ell) if 2*prefix+x>=E else (ell,ell-1)
def check(q):
 cuts=strategy(q);C=[[0,0],[0,0]]
 for i,n in enumerate(q):C[i<=cuts[0]][i<=cuts[1]]+=n
 E=sum(q)
 for row in C:
  for cell in row:
   if cell:require(E>=sum(row)+cell,'Gold deficit failed')
 return cuts,C
stats=Counter();byk=Counter();digest=hashlib.sha256()
# Independent algebraic controls, including zero-bin boundary cases.
for length in range(3,8):
 for A in range(1,6):
  for tail in product(range(A+1),repeat=length-2):
   if any(tail[i]<tail[i+1] for i in range(len(tail)-1)):continue
   q=(A,A)+tail;check(q);stats['integer_monotone_laws']+=1
for n in range(3,7):
 for p in posets(n):
  stats['posets']+=1;U=tuple(sum(1<<v for v in range(n) if p[v]>>u&1) for u in range(n));words=extensions(p)
  for dual in [False,True]:
   D,V=(U,p) if dual else (p,U)
   for u in range(n):
    for v in range(n):
     if u==v or D[u]!=D[v] or V[v]&~V[u]:continue
     extras=[i for i in range(n) if (V[u]&~V[v])>>i&1]
     if not extras or any(not(D[a]>>b&1 or D[b]>>a&1) for a,b in combinations(extras,2)):continue
     extras.sort(key=lambda a:D[a].bit_count());chain=[u]+extras;marked=set(chain+[v]);fibres=defaultdict(lambda:[0]*(len(chain)+1))
     for word in words:
      oriented=tuple(reversed(word)) if dual else word;r={a:i for i,a in enumerate(oriented)};binidx=sum(r[a]<r[v] for a in chain)
      skeleton=tuple(-1 if a in marked else a for a in word);fibres[skeleton][binidx]+=1
     ambient=[0]*(len(chain)+1)
     for skeleton,q in fibres.items():
      cuts,C=check(q);stats['skeleton_laws']+=1;digest.update(repr((p,dual,u,v,extras,skeleton,q,cuts)).encode())
      for i,a in enumerate(q):ambient[i]+=a
     check(ambient);stats['marked_gadgets']+=1;byk[len(extras)]+=1
# Six-point actual-poset endpoint obstruction for two incomparable extras.
p=(0,0,1,7,3,1);words=extensions(p);E=len(words);require(E==37,'Obstruction total mismatch')
@lru_cache(None)
def dp(mask):
 if mask==63:return 1
 return sum(dp(mask|1<<v) for v in range(6) if not mask>>v&1 and not p[v]&~mask)
require(dp(0)==E,'Independent obstruction total mismatch')
endpoints=(0,1,2,5);pairs=list(combinations(endpoints,2));ranks=[{a:i for i,a in enumerate(w)} for w in words];rows=[]
for first in pairs:
 branch_rows=[]
 for answer in [False,True]:
  branch=[r for r in ranks if (r[first[0]]<r[first[1]])==answer]
  if not branch:continue
  options=[]
  for second in pairs:
   a=sum(r[second[0]]<r[second[1]] for r in branch);options.append(max(a,len(branch)-a))
  best=min(options);branch_rows.append({'answer':answer,'survivors':len(branch),'best_larger_second_cell':best,'deficit':E-len(branch)-best})
 margin=min(b['deficit'] for b in branch_rows);require(margin<0,'Restricted strategy unexpectedly succeeds');rows.append({'first_pair':first,'branches':branch_rows,'best_worst_deficit':margin})
require(max(r['best_worst_deficit'] for r in rows)==-3,'Optimal restricted deficit mismatch')
C=[[0,0],[0,0]]
for r in ranks:C[r[4]<r[5]][r[2]<r[5]]+=1
margin=min(E-sum(C[i])-C[i][j] for i in range(2) for j in range(2) if C[i][j]);require(margin==4,'Unrestricted Gold control failed')
print(json.dumps({'status':'EXACT_CONTROLS_PASSED','scope':'All-size theorem proved symbolically; finite controls through six vertices and monotone-bin lengths through seven.','stats':dict(stats),'gadget_counts_by_extra_chain_length':dict(sorted(byk.items())),'digest':digest.hexdigest(),'incomparable_extra_endpoint_obstruction':{'pred':p,'u':0,'v':1,'extras':[2,5],'extensions':E,'optimized_first_pair_rows':rows,'best_restricted_deficit':-3,'unrestricted_strategy':{'comparisons':[[4,5],[2,5]],'table_false_true':C,'minimum_deficit':margin}}},indent=2))
