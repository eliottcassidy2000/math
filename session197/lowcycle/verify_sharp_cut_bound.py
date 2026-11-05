#!/usr/bin/env python3
from pathlib import Path
import itertools as it,json
P=Path(__file__).resolve().parent

def need(x,s):
 if not x:raise RuntimeError(s)
def graphs(n):
 pairs=list(it.combinations(range(n),2))
 for bits in range(1<<len(pairs)):
  r=[0]*n
  for j,(u,v)in enumerate(pairs):
   if bits>>j&1:r[u]|=1<<v
   else:r[v]|=1<<u
  yield r
def reach(r):
 a=[x|1<<i for i,x in enumerate(r)]
 for k in range(len(r)):
  for i in range(len(r)):
   if a[i]>>k&1:a[i]|=a[k]
 return a
def strong(r):return all(x==(1<<len(r))-1 for x in reach(r))
def induce(r,vs):return [sum(1<<j for j,v in enumerate(vs)if r[u]>>v&1)for u in vs]
def scc(r):
 a=reach(r);left=set(range(len(r)));parts=[]
 while left:
  u=min(left);S={v for v in left if a[u]>>v&1 and a[v]>>u&1};parts.append(S);left-=S
 return sorted(parts,key=lambda S:-sum(bool(r[next(iter(S))]>>next(iter(V))&1)for V in parts if V!=S))
def cycles(r):
 n=len(r)
 if n<3:return []
 return [(0,)+w for w in it.permutations(range(1,n))if all(r[u]>>v&1 for u,v in zip((0,)+w,w+(0,)))]
def endpoints(r,vs):
 starts={u:0 for u in vs};ends=starts.copy()
 for w in it.permutations(vs):
  if all(r[u]>>v&1 for u,v in zip(w,w[1:])):starts[w[0]]+=1;ends[w[-1]]+=1
 return starts,ends
def bridges(r):
 ans=[]
 for u in range(len(r)):
  for v in range(len(r)):
   if r[u]>>v&1:
    q=r.copy();q[u]^=1<<v
    if not strong(q):ans.append((u,v))
 return ans

def main():
 total=cuts=products=0
 for n in range(1,7):
  for r in graphs(n):
   total+=1
   if n<3 or not strong(r):continue
   cs=cycles(r);C=len(cs);br=None
   for v in range(n):
    old=[u for u in range(n)if u!=v];small=induce(r,old);parts=scc(small)
    if len(parts)<2:continue
    cuts+=1;parts=[{old[u]for u in S}for S in parts]
    starts,_=endpoints(r,parts[0]);_,ends=endpoints(r,parts[-1])
    A=sum(h for u,h in starts.items()if r[v]>>u&1);B=sum(h for u,h in ends.items()if r[u]>>v&1);prod=A*B
    for S in parts[1:-1]:prod*=sum(endpoints(r,S)[0].values())
    need(prod==C,'cut product');products+=1
    if br is None:br=bridges(r)
    need(br,'smaller bridge-free separating tournament')
    if A<=2 or B<=2:need(br,'endpoint low implies bridge')
 W=[14,116,120,114,33,65,17];cs=cycles(W)
 need(len(cs)==9 and strong(W)and not bridges(W),'sharp witness')
 small=induce(W,list(range(1,7)));parts=scc(small);need(parts==[{0,1,2},{3,4,5}],'witness SCCs')
 separators=[v for v in range(7)if not strong(induce(W,[u for u in range(7)if u!=v]))];need(separators==[0],'witness separator')
 result={'status':'PASS','complete_labelled_tournaments':total,'separating_vertex_cases':cuts,'literal_cut_products':products,'witness_rows':W,'witness_cycles':cs,'witness_strong_arc_deletions':21,'witness_separating_vertices':separators,'all_order_minimum_cycles':9,'minimum_order':7}
 data=json.dumps(result,sort_keys=True,indent=2)+'\n';(P/'sharp_cut_bound_certificate.json').write_text(data);print(data,end='')
if __name__=='__main__':main()
