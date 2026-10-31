from graph_tools import *
from functools import lru_cache
import time,json

def paths(n,cycle=False,count=True):
 a=graph(n)
 if len(comps(a))>1:return dict(count=0,witness=None,states=0)
 degs=[len(a[v]) for v in a]
 if n>1 and (min(degs)==0 or sum(d==1 for d in degs)>2):return dict(count=0,witness=None,states=0)
 if cycle and min(degs)<2:return dict(count=0,witness=None,states=0)
 bits=[sum(1<<(w-1) for w in a[v]) for v in a]
 full=(1<<n)-1;total=0;witness=None;nodes=0
 starts=[0] if cycle else list(range(n))
 for start in starts:
  if not cycle and n>1 and 1 in degs and degs[start]!=1:continue
  def walk(v,mask,trail):
   nonlocal witness,nodes
   nodes+=1
   if mask==full:
    if cycle and not(bits[v]>>start&1):return 0
    if witness is None:witness=[x+1 for x in trail]
    return 1
   avail=full^mask;nexts=bits[v]&avail
   if not nexts:return 0
   # Every unseen vertex must have future incidence; at most one can
   # have only one future incidence in the remaining path problem.
   future=avail|(1<<v)
   ones=0
   for j in range(n):
    if avail>>j&1:
     deg=(bits[j]&future).bit_count()
     if cycle and (bits[j]>>start&1):deg+=1
     if deg==0:return 0
     if deg==1:
      ones+=1
      if ones>(0 if cycle else 1):return 0
   # Connectivity of remaining induced graph plus endpoint.
   seen=1<<v;todo=seen
   while todo:
    q=todo&-todo;todo-=q
    new=bits[q.bit_length()-1]&future&~seen
    seen|=new;todo|=new
   if seen!=future:return 0
   ans=0
   choices=[]
   while nexts:
    q=nexts&-nexts;nexts-=q;j=q.bit_length()-1
    choices.append(((bits[j]&avail).bit_count(),j,q))
   for _,j,q in sorted(choices):
    ans+=walk(j,mask|q,trail+[j])
    if ans and not count:return ans
   return ans
  total+=walk(start,1<<start,[start])
  if total and not count:break
 # When leaves exist, every undirected path is searched twice if two
 # leaves, and once if one leaf. Without leaves, twice (except singleton).
 if count:
  if cycle:total//=2
  elif n>1 and sum(d==1 for d in degs)!=1:total//=2
 return dict(count=total if count else None,exists=bool(witness),witness=witness,states=nodes)
