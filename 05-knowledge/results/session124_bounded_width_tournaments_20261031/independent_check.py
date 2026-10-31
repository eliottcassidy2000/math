"""Exact independent adjacency DP and exhaustive small edge-subset path covers."""
from pathlib import Path
from itertools import product,permutations
from collections import Counter
import json,hashlib
import transfer
ROOT=Path(__file__).resolve().parent

def graph(word,k):
 out=[0]*len(word)
 for j,mask in enumerate(word):
  r=min(j,k)
  if mask<0 or mask>=1<<r:raise ValueError('invalid physical cut')
  for i in range(j):
   backward=j-i<=k and bool(mask&(1<<(i-max(0,j-k))))
   a,b=(j,i) if backward else (i,j)
   out[a]|=1<<b
 return out

def paths(out):
 n=len(out);dp=[[0]*n for _ in range(1<<n)]
 for v in range(n):dp[1<<v][v]=1
 for mask in range(1,1<<n):
  for end,value in enumerate(dp[mask]):
   allowed=out[end]&~mask
   if value:
    while allowed:
     bit=allowed&-allowed;allowed-=bit
     dp[mask|bit][bit.bit_length()-1]+=value
 return sum(dp[-1])

def covers(word,k):
 # Enumerate actual edge subsets, reconstruct full components, then project endpoints.
 out=graph(word,k);n=len(out);edges=[(a,b) for a in range(n) for b in range(n) if out[a]>>b&1]
 result=Counter();offset=max(0,n-k)
 for chosen in range(1<<len(edges)):
  succ=[-1]*n;pred=[-1]*n;valid=True
  for j,(a,b) in enumerate(edges):
   if chosen>>j&1:
    if succ[a]>=0 or pred[b]>=0:valid=False;break
    succ[a]=b;pred[b]=a
  if not valid:continue
  pairs=[];seen=set()
  for start in range(n):
   if pred[start]>=0:continue
   end=start;seen.add(end)
   while succ[end]>=0:end=succ[end];seen.add(end)
   pairs.append((start-offset if start>=offset else -1,end-offset if end>=offset else -1))
  if len(seen)!=n:continue # A cyclic component has no start.
  if sum(a==-1 for a,b in pairs)>1:continue
  result[tuple(sorted(pairs))]+=1
 return result

def endpoint_step(state,r,k,mask):
 # Materialize each endpoint-paired component as an actual short directed path.
 succ0={};nodes=set();starts=[];ends=[];anon=-2
 for a,b in state:
  if a<0:a=anon;anon-=1
  if b<0:b=anon;anon-=1
  starts.append(a);ends.append(b);nodes.update((a,b))
  if a!=b:succ0[a]=b
 incoming=[None]+[e for e in ends if e<0 or not mask>>e&1]
 outgoing=[None]+[a for a in starts if a>=0 and mask>>a&1]
 result=Counter()
 for end in incoming:
  for start in outgoing:
   succ=dict(succ0);allnodes=nodes|{r}
   if end is not None:succ[end]=r
   if start is not None:succ[r]=start
   pred=set(succ.values());seen=set();parts=[]
   def project(a):
    if a<0:return -1
    return a-1 if r==k else a
   for a in allnodes-pred:
    b=a;seen.add(b)
    while b in succ:b=succ[b];seen.add(b)
    parts.append((project(a),project(b)))
   if seen!=allnodes:continue
   if sum(a==-1 for a,b in parts)>1:continue
   result[tuple(sorted(parts))]+=1
 return result
independent_sizes=[];transition_checks=0
for k in range(5):
 states={()}
 for r in range(k):
  states={v for u in states for a in range(1<<r) for v in endpoint_step(u,r,k,a)}
 todo=list(states)
 while todo:
  u=todo.pop()
  for a in range(1<<k):
   step=endpoint_step(u,k,k,a)
   if dict(step)!=dict(transfer.step(u,k,k,a)):raise RuntimeError('endpoint transition')
   transition_checks+=1
   for v in step:
    if v not in states:states.add(v);todo.append(v)
 if states!=set(transfer.carrier(k)):raise RuntimeError('carrier state set')
 independent_sizes.append(len(states))
small=0
for k,maxn in ((0,4),(1,5),(2,5),(3,4)):
 for n in range(1,maxn+1):
  for w in product(*(range(1<<min(i,k)) for i in range(n))):
   v={():1}
   for i,a in enumerate(w):v=transfer.advance(v,min(i,k),k,a)
   if dict(v)!=dict(covers(w,k)):raise RuntimeError(('full cover state',k,w))
   small+=1
checks=0;by_order={};seen=set()
for n in range(1,9):
 total=0
 for w in product(*(range(1<<min(i,2)) for i in range(n))):
  h=paths(graph(w,2))
  if h!=transfer.count(w,2) or h==47:raise RuntimeError(('width2',w,h))
  seen.add(h);checks+=1;total+=1
 by_order[n]=total
if checks!=10923:raise RuntimeError('universe size')
ws=[];got=set();max_n=0
for line in (ROOT/'width3_witnesses.tsv').read_text().splitlines():
 h,*w=map(int,line.split());ws.append((h,w));got.add(h);max_n=max(max_n,len(w))
 if paths(graph(w,3))!=h or transfer.count(w,3)!=h:raise RuntimeError(('witness',h,w))
if len(ws)!=825 or len(got)!=825 or got!={h for h in range(1,1654,2) if h not in (7,21)} or max_n>9:raise RuntimeError('coverage')
w=(0,0,0,1,3,2,1);out=graph(w,3)
permutation47=sum(all(out[a]>>b&1 for a,b in zip(p,p[1:])) for p in permutations(range(7)))
if permutation47!=47 or transfer.count(w,3)!=47:raise RuntimeError('47 control')
rejected=0
if paths(out)==49:raise RuntimeError('mutated count accepted')
rejected+=1
try:graph((0,2),3)
except ValueError:rejected+=1
else:raise RuntimeError('invalid mask accepted')
mut=list(w);mut[-1]^=1
if paths(graph(mut,3))==47:raise RuntimeError('chosen edge mutation not hostile')
rejected+=1
states=[len(transfer.carrier(k)) for k in range(5)]
if states!=[2,6,25,116,567]:raise RuntimeError('carrier size')
if not 47**3<5**8:raise RuntimeError('Busch cutoff')
print(json.dumps(dict(full_edge_subset_cover_comparisons=small,width2_checks=checks,width2_order_counts=by_order,width2_47_absent=True,width3_witnesses=len(ws),width3_allowed_prefix=1653,maximum_witness_order=max_n,witness47_word=w,witness47_literal_permutations=permutation47,carrier_states=states,independently_reconstructed_carriers=independent_sizes,endpoint_transition_checks=transition_checks,rejected_mutations=rejected,witness_sha256=hashlib.sha256((ROOT/'width3_witnesses.tsv').read_bytes()).hexdigest()),sort_keys=True))
