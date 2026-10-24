#!/usr/bin/env python3
"""Direct complete-state Kruskal checker, independent of helper-sheet BFS."""
from pathlib import Path
from itertools import combinations
from collections import Counter,defaultdict
import json
HERE=Path(__file__).resolve().parent
D=json.loads((HERE/'inputs.json').read_text());CORE=D['cores'];D={int(k):v for k,v in D['records'].items()};FREE=sorted(set(D)-set(CORE)-{368,776,3139})
ADD=lambda a,b:((((a>>4)+(b>>4))%128)<<4)|((a^b)&15)
NEG=lambda q:((-(q>>4)%128)<<4)|(q&15)
def charge(ss):
 q=0
 for s in ss:q=ADD(q,D[s]['label'])
 return q
def check(ok,why):
 if not ok:raise RuntimeError(why)
def main():
 out=[];hostiles=[]
 for k in (0,1):
  for sign in (1,-1):
   fibres=defaultdict(list)
   for C in combinations(CORE,k):
    for H in combinations(FREE,2):
     q=ADD(charge(C) if sign==1 else NEG(charge(C)),charge(H));fibres[q].append((C,H))
   hist=Counter()
   for q,states in sorted(fibres.items()):
    edges=[]
    for i,A in enumerate(states):
     a=set(A[0]+A[1])
     for j in range(i):
      b=set(states[j][0]+states[j][1]);w=sum(D[s]['cost'] for s in a^b);edges.append((w,j,i))
    edges.sort();p=list(range(len(states)))
    def find(x):
     while p[x]!=x:x=p[x]
     return x
    bottleneck=0;selected=[]
    for w,i,j in edges:
     a,b=find(i),find(j)
     if a!=b:p[a]=b;bottleneck=w;selected.append((w,i,j))
    check(len(selected)==len(states)-1,'complete nonempty graph spans')
    hist[bottleneck]+=1
    if k==1 and bottleneck==33:
     p=list(range(len(states)))
     for w,i,j in edges:
      if w<=32:
       a,b=find(i),find(j);p[a]=b
     comps=defaultdict(list)
     for i in range(len(states)):comps[find(i)].append(i)
     cross=[(w,i,j) for w,i,j in edges if find(i)!=find(j)]
     check(len(comps)>1 and min(w for w,i,j in cross)==33,'sharp32 obstruction actual complete fibre')
     hostiles.append(dict(sign=sign,q=q,states=[dict(core=list(c),helpers=list(h)) for c,h in states],components32=list(comps.values()),minimum_crossing_clock=min(w for w,i,j in cross),crossing_witness=min(cross)))
   out.append(dict(grade=k,sign=sign,states=sum(map(len,fibres.values())),fibres=len(fibres),minimum_fibre_size=min(map(len,fibres.values())),bottleneck_histogram=dict(sorted(hist.items()))))
 check(max(int(w) for x in out for w in x['bottleneck_histogram'])==33,'uniform low0/1 sharp33')
 receipt=dict(summary=out,sharp33_hostiles=hostiles)
 print(json.dumps(receipt,sort_keys=True,indent=2))
if __name__=='__main__':main()
