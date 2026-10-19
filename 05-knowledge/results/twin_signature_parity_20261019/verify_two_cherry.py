#!/usr/bin/env python3
"""Exact source-side deletion/reversal certificates, not random host evidence."""
from verify_twin_signature import *
from math import factorial
import json

def forest_aut(edges,n,state):
 a=[[] for _ in range(n)];out=[0]*n
 for (u,v),s in zip(edges,state):
  if not s:continue
  a[u].append(v);a[v].append(u)
  if s==1:out[u]|=1<<v
  else:out[v]|=1<<u
 seen=set();codes=[];group=1
 for v in range(n):
  if v in seen:continue
  vv=[v];seen.add(v)
  for u in vv:
   for w in a[u]:
    if w not in seen:seen.add(w);vv.append(w)
  ix={u:i for i,u in enumerate(vv)}
  aa=[[ix[w] for w in a[u]] for u in vv]
  oo=[sum(1<<ix[w] for w in a[u] if out[u]>>w&1) for u in vv]
  code,g=ocanon_aut(aa,oo);codes.append(code);group*=g
 for m in Counter(codes).values():group*=factorial(m)
 return group,sum(not x for x in a)

def main():
 e=((0,1),(0,2),(0,3),(0,4),(1,5),(1,6),(2,7),(7,8));n=9;a=adj_of(e,n)
 require(aut(a)==4 and afree(a)==16 and not hereditary(a),'tree identities')
 direct=0;reductions=0;asymmetric=0
 for mask in range(256):
  state=[2 if mask>>i&1 else 1 for i in range(8)]
  g,_=forest_aut(e,n,state);asymmetric+=g==1
  if state[2]==state[3] or state[4]==state[5]:
   require(g%2==0,'direct swap fails');direct+=1;continue
  x=state[:];x[2]=3-x[2]
  require(forest_aut(e,n,x)[0]%2==0,'first reversed term not even')
  x=state[:];x[2]=0;y=x[:];y[4]=3-y[4]
  require(forest_aut(e,n,y)[0]%2==0,'second reversed term not even')
  x[4]=0;g,isolates=forest_aut(e,n,x)
  require(g%2==0 and isolates==2,'terminal isolate swap fails');reductions+=1
 require((direct,reductions,asymmetric)==(192,64,64),'orientation route totals')
 a,out=orient_data(e,n,40)
 require(ocanon_aut(a,out)[1]==1,'asymmetric concrete orientation')
 require(embedding_count(e,n,40,transitive(n))==260,'transitive count')
 require(brute_emb(e,n,40,out)==1,'self-host count')
 census=[];expected={7:(2,0),8:(6,0),9:(14,3),10:(38,12),11:(94,33)}
 for n,es in trees(11).items():
  if n<7:continue
  total=outside=0
  for e in es:
   a=adj_of(e,n)
   two=sum(sum(len(a[w])==1 for w in a[v])>=2 for v in range(n))>=2
   total+=two;outside+=two and not hereditary(a) and afree(a)>0
  require((total,outside)==expected[n],f'census at {n}')
  census.append({'n':n,'two_cherry_trees':total,'outside_hereditary_and_auto_zero':outside})
 print(json.dumps({'status':'PASS','direct_automorphism_routes':direct,'two_deletion_routes':reductions,'asymmetric_labelled_orientations':asymmetric,'asymmetric_orientation_classes':16,'census':census},indent=2))
if __name__=='__main__':main()
