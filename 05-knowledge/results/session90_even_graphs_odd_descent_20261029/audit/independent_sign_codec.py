#!/usr/bin/env python3
"""Verify exact Royle sign under generic sparse-graph encoding through n=5."""
from independent_sparse_audit import p
from itertools import combinations,permutations
import json
C=0
def ck(x,msg):
 global C
 C+=1
 if not x:raise RuntimeError(msg)
def main():
 objects=maps=0;odd=even=0
 for n in range(6):
  es=list(combinations(range(n),2));idx={e:i for i,e in enumerate(es)}
  for mask in range(1<<len(es)):
   g=p.encode_graph(n,mask);objects+=1
   for pi in permutations(range(n)):
    moved=sum(((mask>>j)&1)<<idx[tuple(sorted((pi[a],pi[b])))] for j,(a,b) in enumerate(es))
    if moved!=mask:continue
    source=sum(pi[a]>pi[b] for j,(a,b) in enumerate(es) if mask>>j&1)%2
    q=list(range(len(g)))
    if n:
     for i in range(n):
      for j in range(11):q[12+11*i+j]=12+11*pi[i]+j
    ck(all({q[w] for w in g[v]}==g[q[v]] for v in range(len(g))),'lifted graph automorphism')
    target=sum(q[v]>q[w] for v in range(len(g)) for w in g[v] if v<w)%2
    ck(target==source,'Royle endpoint sign preserved')
    maps+=1;odd+=source;even+=1-source
 print(json.dumps(dict(status='PASS',checks=C,labelled_graphs=objects,lifted_automorphisms=maps,positive_maps=even,negative_maps=odd),indent=2))
if __name__=='__main__':main()
