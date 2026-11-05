#!/usr/bin/env python3
"""Physical controls for permanent-module compression and capped continuations."""
import itertools,json,hashlib,random
from functools import lru_cache

def need(p,s):
 if not p:raise RuntimeError(s)
def bounds(N,k):
 a=b=0
 while 3**(a+1)<=N:a+=1
 while 5**(b+1)<=N**3:b+=1
 return a,b,b+1,a+b+k+(a+k+1)*(b+1)
def valid(r):
 for i in range(len(r)):
  need(not(r[i]>>i&1),'loop')
  for j in range(i+1,len(r)):need((r[i]>>j&1)+(r[j]>>i&1)==1,'orientation')
def induced(r,V):return tuple(sum((r[u]>>v&1)<<j for j,v in enumerate(V)) for u in V)
def scc(r):
 n=len(r);reach=[r[i]|1<<i for i in range(n)]
 for v in range(n):
  for u in range(n):
   if reach[u]>>v&1:reach[u]|=reach[v]
 seen=set();parts=[]
 for u in range(n):
  if u in seen:continue
  C=[v for v in range(n)if reach[u]>>v&1 and reach[v]>>u&1];seen.update(C);parts.append(C)
 parts.sort(key=lambda C:-reach[C[0]].bit_count())
 for i,C in enumerate(parts):
  for D in parts[i+1:]:need(all(r[u]>>v&1 for u in C for v in D),'SCC total order')
 return parts
@lru_cache(None)
def raw(r):
 n=len(r)
 if not n:return 1
 d=[[0]*n for _ in range(1<<n)]
 for v in range(n):d[1<<v][v]=1
 for S in range(1,1<<n):
  for v,x in enumerate(d[S]):
   if not x:continue
   opts=r[v]&~S
   while opts:
    bit=opts&-opts;opts-=bit;d[S|bit][bit.bit_length()-1]+=x
 return sum(d[-1])
@lru_cache(None)
def hp(r):
 h=1
 for C in scc(r):
  if len(C)>1:h*=raw(induced(r,C))
 return h

def birth(r,p,op):
 side,mask,drop=op;n=len(r);rr=list(r)+[0]
 for v in range(n):
  selected=any(v==u and mask>>j&1 for j,u in enumerate(p))
  win=selected if side==0 else not selected
  if win:rr[n]|=1<<v
  else:rr[v]|=1<<n
 pp=list(p)+[n];del pp[drop]
 return tuple(rr),tuple(pp)
def compress(r,p,N):
 M=bounds(N,len(p))[2];keep=[];run=[]
 def flush():
  keep.extend(run[:M]);run.clear()
 for C in scc(r):
  if len(C)==1 and C[0]not in p:run.extend(C)
  else:flush();keep.extend(C)
 flush();keep.sort();ix={v:i for i,v in enumerate(keep)}
 out=induced(r,keep);ports=tuple(ix[v]for v in p)
 need(hp(out)==hp(r),'current exact H preserved')
 if hp(r)<=N:need(len(out)<=bounds(N,len(p))[3],'bounded state order')
 return out,ports,keep

def alphabet(k):return tuple(itertools.product(range(2),range(1<<k),range(k+1)))
def step(state,op,N):
 if state is None:return None
 r,p=birth(*state,op)
 if hp(r)>N:return None
 rr,pp,_=compress(r,p,N);return rr,pp

def trans(n):return tuple(sum(1<<j for j in range(i+1,n))for i in range(n))
def ordinal(parts):
 sizes=[len(r)for r in parts];n=sum(sizes);out=[];offset=0
 for r in parts:
  later=((1<<n)-1)^((1<<(offset+len(r)))-1)
  out.extend((v<<offset)|later for v in r);offset+=len(r)
 return tuple(out)

def main():
 digest=hashlib.sha256();module_controls=0;word_controls=0;compressions=0
 # Each long block is genuinely unexposed. Exterior pieces include nontrivial SCCs.
 cases=[(1,ordinal([trans(5),trans(1)]),(5,),list(range(5))),
        (1,ordinal([trans(1),trans(5),trans(1)]),(0,6),list(range(1,6))),
        (3,ordinal([(2,4,1),trans(5),trans(1)]),(0,8),list(range(3,8))),
        (3,ordinal([trans(1),trans(5),(2,4,1)]),(0,6),list(range(1,6))),
        (5,ordinal([trans(6),(2,4,1)]),(6,),list(range(6))),
        (3,ordinal([trans(5),(2,4,1)]),(),list(range(5)))]
 for N,r,p,B in cases:
  valid(r);rr,pp,keep=compress(r,p,N);need(len(rr)<len(r),'actual compression triggered');compressions+=1
  for word in itertools.product(alphabet(len(p)),repeat=2):
   x,q=r,p;y,qq=rr,pp;retained=keep[:]
   for op in word:
    z=len(x);x,q=birth(x,q,op);y,qq=birth(y,qq,op);retained.append(z)
    need(induced(x,retained)==y,'same physical retained graph under continuation')
    need(tuple(retained[v]for v in qq)==q,'corresponding ordered ports')
    outside=set(range(len(x)))-set(B)
    need(all(len({(x[v]>>b)&1 for b in B})==1 for v in outside),'module survives every birth')
    parts=scc(x);hits=[C for C in parts if set(C)&set(B)]
    need((all(len(C)==1 for C in hits))or(len(hits)==1 and set(B)<=set(hits[0])),'no partial SCC merge')
    need(min(hp(x),N+1)==min(hp(y),N+1),'all tested capped continuation outputs agree')
    digest.update(str((N,x,q,y,qq)).encode());module_controls+=1
 # Exhaustive words test repeated compression, not only one initial replacement.
 for N,seed,ports,depth in [(1,trans(2),(0,),5),(3,(2,4,1),(0,),4),(5,(2,4,1),(0,1),3),(1,trans(1),(),8)]:
  for word in itertools.product(alphabet(len(ports)),repeat=depth):
   r,p=seed,ports;cr,cp,_=compress(r,p,N);state=(cr,cp)
   for op in word:
    r,p=birth(r,p,op);state=step(state,op,N)
    got=N+1 if state is None else hp(state[0])
    need(got==min(hp(r),N+1),'deterministic machine exact output')
    word_controls+=1
 # Complete small machines, no arbitrary search cutoff.
 machines=[]
 for k in [0,1,2]:
  N=1;seed=trans(max(1,k));ports=tuple(range(k));r,p,_=compress(seed,ports,N);initial=(r,p);seen={initial};todo=[initial]
  for state in todo:
   for op in alphabet(k):
    out=step(state,op,N)
    if out not in seen:seen.add(out);todo.append(out)
  machines.append({'cap':N,'ports':k,'states_including_reachable_sink':len(seen),'max_graph_order':max(len(s[0])for s in seen if s is not None)})
 rng=random.Random(197555);raw_checks=0
 for n in range(1,8):
  for _ in range(12):
   r=[0]*n
   for i in range(n):
    for j in range(i+1,n):
     if rng.randrange(2):r[i]|=1<<j
     else:r[j]|=1<<i
   r=tuple(r);need(raw(r)==hp(r),'independent whole-graph/SCC count');raw_checks+=1
 for N in range(1,1001):
  a,b,M,V=bounds(N,2)
  need(3**a<=N<3**(a+1)and 5**b<=N**3<5**M,'exact integer thresholds')
 result={'status':'PASS','initial_long_block_compressions':compressions,'module_prefix_controls':module_controls,'iterated_machine_prefix_controls':word_controls,'complete_small_machines':machines,'whole_graph_SCC_count_controls':raw_checks,'threshold_caps_checked':1000,'physical_controls_sha256':digest.hexdigest()}
 print(json.dumps(result,sort_keys=True,indent=2))
if __name__=='__main__':main()
