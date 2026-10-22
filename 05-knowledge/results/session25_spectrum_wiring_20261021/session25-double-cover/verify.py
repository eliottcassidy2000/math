from itertools import combinations, permutations
from collections import Counter, defaultdict
import json

def parity(p):
 return sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))%2

def run(n):
 pairs=list(combinations(range(n),2)); ix={p:i for i,p in enumerate(pairs)}; N=len(pairs); full=(1<<N)-1; m=(n-1)*(n-2)//2
 groups=[]
 for p in permutations(range(n)):
  perm=[]; flip=0
  for i,(a,b) in enumerate(pairs):
   x,y=p[a],p[b]; j=ix[tuple(sorted((x,y)))]; perm.append(j)
   if x>y: flip|=1<<j
  vals=[0]*(1<<N)
  for t in range(1,1<<N):
   bit=t&-t; vals[t]=vals[t^bit]^(1<<perm[bit.bit_length()-1])
  for eps in range(2):
   ff=flip^(full if eps else 0); tr=[x^ff for x in vals]; ch=parity(p)^(eps*(N%2))
   assert all(((t.bit_count()+tr[t].bit_count())%2)==ch for t in range(1<<N))
   groups.append((ch,tr))
 edges=set(); pathof={}
 for p in permutations(range(n)):
  pm=0; fixed=0
  for a,b in zip(p,p[1:]):
   j=ix[tuple(sorted((a,b)))]; pm|=1<<j
   if a<b: fixed|=1<<j
  free=full^pm; sub=free
  while True:
   t=fixed|sub; u=t^free; e=tuple(sorted((t,u))); edges.add(e);pathof[e]=p
   if sub==0: break
   sub=(sub-1)&free
 assert len(edges)==__import__('math').factorial(n)*2**(m-1)
 vc={}; vk={}
 for t in range(1<<N):
  vc[t]=min(tr[t] for ch,tr in groups); vk[t]=min(tr[t] for ch,tr in groups if not ch)
 fibres=defaultdict(set)
 for t in vc: fibres[vc[t]].add(vk[t])
 assert all(len(f)==2 for f in fibres.values())
 ec={}; ek={}; dk={}; dc={}; colors={}; inversions=0
 for e in edges:
  a,b=e; orbit=[]; korbit=[]; blue=False
  for ch,tr in groups:
   aa,bb=tr[a],tr[b]; ee=tuple(sorted((aa,bb))); orbit.append(ee)
   if not ch: korbit.append(ee)
   if aa==b and bb==a: inversions+=1
  ec[e]=min(orbit); ek[e]=min(korbit)
  p=pathof[e]; r={p[i]:p[n-1-i] for i in range(n)}
  def arc(t,a,b): return bool(t&(1<<ix[tuple(sorted((a,b)))]))==(a<b)
  blue=all(arc(a,x,y)!=arc(a,r[x],r[y]) for x,y in pairs)
  colors[e]=blue
 assert inversions==0
 ef=defaultdict(set)
 for e in edges:
  ef[ec[e]].add(ek[e]); assert colors[e]==colors[ec[e]]
 assert all(len(f)==2 for f in ef.values())
 # Darts preserve the endpoint distinction even over base loops.
 for e in edges:
  for a,b in (e,e[::-1]):
   d=(a,b); dc[d]=min((tr[a],tr[b]) for ch,tr in groups); dk[d]=min((tr[a],tr[b]) for ch,tr in groups if not ch)
 # Each base dart and choice of its K vertex has one K dart.
 df=defaultdict(set)
 for d in dc: df[(dc[d],vk[d[0]])].add(dk[d])
 assert all(len(s)==1 for s in df.values())
 assert len(df)==2*len(set(dc.values()))
 assert len(set(dc.values()))==2*len(ef)
 # All lifted endpoint signs obey prescribed voltage, including loops.
 for a,b in set(ek.values()): assert (a.bit_count()+b.bit_count())%2==m%2
 base=set(ec.values()); loops=sum(vc[a]==vc[b] for a,b in base)
 # Check actual base edge colors agree with cube <complement,reflection> count.
 blue=sum(colors[e] for e in base); black=len(base)-blue; dim=(m+(n-1)//2)//2
 assert blue==2**(dim-1)
 assert black==(2**m-2**dim)//4
 return dict(n=n,labelled_vertices=1<<N,labelled_edges=len(edges),base_vertices=len(fibres),cover_vertices=len(set(vk.values())),base_edges=len(base),cover_edges=len(set(ek.values())),base_loops=loops,blue_edges=blue,fundamental_black_edges=black,base_darts=len(set(dc.values())),cover_darts=len(set(dk.values())),edge_inversions=inversions,voltage=m%2)

def comp(vertices,edges,start):
 seen={start};todo=[start]
 while todo:
  v=todo.pop()
  for a,b in edges:
   if a==v and b not in seen: seen.add(b);todo.append(b)
   if b==v and a not in seen: seen.add(a);todo.append(a)
 return seen

def bip(edges,vs):
 color={}
 for root in vs:
  if root in color: continue
  color[root]=0;todo=[root]
  while todo:
   v=todo.pop()
   for a,b in edges:
    if v not in (a,b): continue
    w=b if a==v else a
    if w in color and color[w]==color[v]: return False
    if w not in color: color[w]=1-color[v];todo.append(w)
 return True

def bridge_test():
 checked=0;cases=0
 for n in range(1,5):
  pairs=list(combinations(range(n),2))+[(i,i) for i in range(n)]
  for bits in range(1<<len(pairs)):
   edges=[pairs[i] for i in range(len(pairs)) if bits>>i&1]
   for duplicate in (False,True):
    gg=edges*(2 if duplicate else 1)
    h=[((a,s),(b,1-s)) for a,b in gg for s in (0,1)]
    for i,(a,b) in enumerate(gg):
     cut=gg[:i]+gg[i+1:];ca=comp(range(n),cut,a);cb=comp(range(n),cut,b)
     expect=b not in ca and (bip(cut,ca) or bip(cut,cb))
     for s in (0,1):
      j=2*i+s;hc=h[:j]+h[j+1:];x,y=h[j]
      assert (y not in comp([],hc,x))==expect
     checked+=1
    cases+=1
 return dict(graph_cases=cases,base_edge_checks=checked)
if __name__=='__main__':
 for n in range(3,6): print(json.dumps(run(n)),flush=True)
 print(json.dumps(bridge_test()),flush=True)
