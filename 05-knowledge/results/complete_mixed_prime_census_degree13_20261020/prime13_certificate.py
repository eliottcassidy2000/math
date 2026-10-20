#!/usr/bin/env python3
"""Exact prime-degree classification control using all skew sets, no group database."""
from itertools import permutations
from collections import Counter
import json
from pathlib import Path

def require(v,msg):
 if not v:raise RuntimeError(msg)

def main():
 p=13;rows=[]
 for bits in range(1<<6):
  S=tuple(sorted(i if bits>>(i-1)&1 else p-i for i in range(1,7)))
  A=[[int((v-u)%p in S) for v in range(p)] for u in range(p)]
  R=tuple(v for v in range(1,p) if v not in S)
  sig={v:tuple(A[v][s] for s in S) for v in R}
  require(len(set(sig.values()))==len(R),'in-neighbor signatures must be distinct')
  inv={s:v for v,s in sig.items()};auts=[]
  for image in permutations(S):
   if any(A[u][v]!=A[x][y] for u,x in zip(S,image) for v,y in zip(S,image)):continue
   q=[None]*p;q[0]=0
   for u,v in zip(S,image):q[u]=v
   for u in R:
    # Find v with A[v][q(s)] = A[u][s] for all s in S.
    target=tuple(A[u][S[image.index(s)]] for s in S)
    q[u]=inv.get(target)
   if None in q or len(set(q))!=p:continue
   if all(A[u][v]==A[q[u]][q[v]] for u in range(p) for v in range(p)):auts.append(tuple(q))
  aff=[tuple(a*i%p for i in range(p)) for a in range(1,p) if {a*s%p for s in S}==set(S)]
  require(set(auts)==set(aff),'all root stabilizers must be affine')
  rows.append({'skew_set':S,'root_automorphisms':auts,'automorphism_order':p*len(auts)})
 require(Counter(x['automorphism_order'] for x in rows)=={13:60,39:4},'prime13 histogram')
 data={'status':'PASS','scope':'all 64 cyclic tournaments at degree13; each root stabilizer checked over6! out-neighbor bijections with uniquely forced in-neighbor images','automorphism_order_histogram':dict(Counter(x['automorphism_order'] for x in rows)),'rows':rows}
 Path(__file__).with_suffix('.json').write_text(json.dumps(data,indent=2)+'\n');print(data['automorphism_order_histogram'])
if __name__=='__main__':main()
