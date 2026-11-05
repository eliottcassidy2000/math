"""Independent SCC/endpoint-DP and literal graph-map audit. No author imports."""
from pathlib import Path
from functools import lru_cache
from collections import Counter
import hashlib,json
P=Path(__file__).resolve().parent

def induced(T,vs):return tuple(sum(1<<j for j,v in enumerate(vs)if T[u]>>v&1)for u in vs)
def converse(T):return tuple(sum(1<<j for j in range(len(T))if T[j]>>i&1)for i in range(len(T)))
def path_counts(T):
 n=len(T);full=(1<<n)-1
 @lru_cache(None)
 def start(S,v):
  U=S^(1<<v)
  if not U:return 1
  ans=0;bits=U&T[v]
  while bits:
   b=bits&-bits;bits-=b;ans+=start(U,b.bit_length()-1)
  return ans
 return [start(full,v)for v in range(n)]
def cycle_count(T):
 n=len(T)
 if n<3:return 0
 @lru_cache(None)
 def finish(S,v):
  if not S:return int(bool(T[v]&1))
  ans=0;bits=S&T[v]
  while bits:
   b=bits&-bits;bits-=b;ans+=finish(S^b,b.bit_length()-1)
  return ans
 return finish(((1<<n)-1)^1,0)
def scc(T,vs):
 # Kosaraju, using literal adjacency and a separate reverse traversal.
 seen=set();finish=[];vs=set(vs)
 def visit(v):
  seen.add(v)
  for u in sorted(vs):
   if u not in seen and T[v]>>u&1:visit(u)
  finish.append(v)
 for v in sorted(vs):
  if v not in seen:visit(v)
 seen=set();out=[]
 def reverse(v,part):
  seen.add(v);part.append(v)
  for u in sorted(vs):
   if u not in seen and T[u]>>v&1:reverse(u,part)
 for v in reversed(finish):
  if v not in seen:
   part=[];reverse(v,part);out.append(sorted(part))
 for a,A in enumerate(out):
  for B in out[a+1:]:
   if not all(T[x]>>y&1 for x in A for y in B):raise RuntimeError('SCC order')
 return out

def add(T,outneighbors):
 n=len(T);A=list(T)+[sum(1<<v for v in outneighbors)]
 for v in range(n):
  if v not in outneighbors:A[v]|=1<<n
 return tuple(A),n

def literal(kind,V,b,v):
 A,a=add(V,{b});A,y=add(A,{b})
 if kind=='N':root=a
 elif kind=='M':A,z=add(A,{a});root=y
 elif kind=='Q':A,root=add(A,{a})
 else:raise RuntimeError('kind')
 # All tail vertices are appended in increasing dominance order. A bit1
 # means the tail vertex beats the retained root; others are sinks.
 for bit in bin(v-1)[2:][::-1]if v>1 else'':
  A,_=add(A,{root}if bit=='1'else set())
 return A,root

def recover(T,r):
 blocks=scc(T,[v for v in range(len(T))if v!=r]);U=blocks[0]
 if cycle_count(induced(T,U))!=1:return None
 if any(len(x)!=1 for x in blocks[1:]):raise RuntimeError('non-singleton tail')
 tail=[x[0]for x in blocks[1:]];param=1+sum(1<<i for i,x in enumerate(tail)if T[x]>>r&1)
 entries=[x for x in U if T[r]>>x&1]
 if len(entries)==2:
  uniqueout={x:[y for y in U if T[x]>>y&1]for x in U}
  possibilities=[(x,y)for x in entries for y in entries if x!=y and uniqueout[x]==[y]]
  if possibilities:
   if len(possibilities)!=1:raise RuntimeError('adjacent uniqueness')
   a,b=possibilities[0];V=[x for x in U if x!=a];kind='N';mapping=V+[r,a]+tail
   if param<2:raise RuntimeError('adjacent empty tail')
  else:
   possibilities=[(u,b,c)for u in entries for c in entries if u!=c for b in U if b not in entries and uniqueout[u]==[b]and uniqueout[b]==[c]]
   if len(possibilities)!=1:raise RuntimeError('outer pair uniqueness')
   u,b,c=possibilities[0];V=[x for x in U if x not in(u,b)];kind='M';mapping=V+[b,r,u]+tail;b=c
 elif len(entries)==1:
  a=entries[0];possibilities=[]
  for y in U:
   if y==a:continue
   ys=[v for v in U if T[y]>>v&1]
   if len(ys)==1 and set(v for v in U if T[a]>>v&1)=={y,ys[0]}:possibilities.append((y,ys[0]))
  if len(possibilities)!=1:raise RuntimeError('clone uniqueness')
  y,b=possibilities[0];V=[x for x in U if x not in(a,y)];kind='Q';mapping=V+[a,y,r]+tail
 else:raise RuntimeError('entry number')
 source=induced(T,V);source_root=V.index(b)
 if len(V)<2 or path_counts(source)[source_root]!=1:raise RuntimeError('source membership')
 A,root=literal(kind,source,source_root,param)
 if mapping[root]!=r or induced(T,mapping)!=A:raise RuntimeError(('literal map',kind,T,r,param))
 return kind,len(source),param

def main():
 stats={n:Counter()for n in range(5,11)};records=[];roots=maps=0
 for line in(P/'c2-inputs.csv').read_text().splitlines():
  z=list(map(int,line.split(',')));n,C=z[:2];T=tuple(z[2:])
  if len(T)!=n or C!=2 or cycle_count(T)!=2:raise RuntimeError('bank cycle count')
  ss=path_counts(T);found=Counter()
  for r,s in enumerate(ss):
   if s!=2:continue
   roots+=1;rec=recover(T,r)
   if rec:
    maps+=1;found[rec[0]]+=1;stats[n]['root_'+rec[0]]+=1
  stats[n]['C2']+=1;stats[n]['Gamma']+=bool(found);stats[n]['qualifying']+=sum(found.values())
  records.append((n,tuple(T),sum(ss),dict(found)))
 expectedC2=[2,5,18,60,191,589];expectedGamma=[1,5,18,57,169,482];expectedRoots=[3,10,31,91,258,715]
 for i,n in enumerate(range(5,11)):
  if [stats[n][k]for k in('C2','Gamma','qualifying')]!=[expectedC2[i],expectedGamma[i],expectedRoots[i]]:raise RuntimeError(('inventory',n,stats[n]))
 hostiles=[]
 for H,T in[(77,(2,4,9,19,103,79,143,63)),(289,(2,4,9,83,39,79,151,319,639,255))]:
  result=[]
  for A in(T,converse(T)):
   ss=path_counts(A)
   if sum(ss)!=H or cycle_count(A)!=2:raise RuntimeError('hostile count')
   rs=[r for r,s in enumerate(ss)if s==2];data=[]
   for r in rs:
    blocks=scc(A,[x for x in range(len(A))if x!=r]);U=blocks[0]
    data.append(dict(root=r,initial=U,cycles=cycle_count(induced(A,U)),compiled=recover(A,r)))
   result.append(data)
  if any(x['compiled']for x in result[0])or H==289 and any(x['compiled']for x in result[1]):raise RuntimeError('hostile qualification')
  hostiles.append(dict(H=H,orientations=result))
 # Exact formal coefficient expansion of canonical rooted GF.
 coeff=[0]*61;den=[1,-5,7,-2];num={5:3,6:-5,7:2}
 for n in range(61):coeff[n]=num.get(n,0)-sum(den[j]*coeff[n-j]for j in range(1,4)if n>=j)
 if coeff[5:11]!=expectedRoots:raise RuntimeError('rooted GF')
 text=json.dumps(dict(status='PASS',complete_C2_graphs=len(records),start_two_roots=roots,physical_inverse_maps=maps,inventory={n:dict(c)for n,c in stats.items()},hostiles=hostiles,rooted_GF_coefficients_through60=coeff,input_sha256=hashlib.sha256((P/'c2-inputs.csv').read_bytes()).hexdigest()),sort_keys=True,indent=2)+'\n'
 print(text,end='')
if __name__=='__main__':main()
