#!/usr/bin/env python3
"""Pruefer-universe audit of tree nonedges <-> unordered pairs of tree edges."""
from itertools import product,permutations,combinations
from collections import deque
import random,json
C=0
def ck(x,msg):
 global C
 C+=1
 if not x:raise RuntimeError(msg)
def tree(n,word):
 if n==1:return [set()]
 d=[1]*n
 for v in word:d[v]+=1
 g=[set() for _ in range(n)]
 for v in word:
  u=next(i for i in range(n) if d[i]==1)
  g[u].add(v);g[v].add(u);d[u]-=1;d[v]-=1
 a,b=[i for i in range(n) if d[i]==1];g[a].add(b);g[b].add(a)
 return g
def path(g,u,v):
 prev={u:None};todo=deque([u])
 while v not in prev:
  x=todo.popleft()
  for y in g[x]:
   if y not in prev:prev[y]=x;todo.append(y)
 out=[v]
 while out[-1]!=u:out.append(prev[out[-1]])
 return tuple(reversed(out))
def edge(a,b):return tuple(sorted((a,b)))
def boundary(g,u,v):
 p=path(g,u,v)
 if len(p)<3:raise ValueError('not a nonedge')
 return edge(*p[:2]),edge(*p[-2:])
def outer(g,e,f):
 # Remove e and find which endpoint component does not meet f.
 todo=[e[0]];seen={e[0]}
 while todo:
  x=todo.pop()
  for y in g[x]:
   if edge(x,y)==e:continue
   if y not in seen:seen.add(y);todo.append(y)
 return e[1] if f[0] in seen else e[0]
def main():
 rng=random.Random(901);rows=[];roundtrips=0
 for n in range(1,8):
  trees=0
  for word in product(range(n),repeat=max(0,n-2)):
   g=tree(n,word);es=[(u,v) for u in range(n) for v in g[u] if u<v];non=[(u,v) for u in range(n) for v in range(n) if u!=v and v not in g[u]]
   images={boundary(g,u,v) for u,v in non}
   ck(images=={(e,f) for e in es for f in es if e!=f},'ordered boundary bijection')
   for u,v in non:
    e,f=boundary(g,u,v)
    ck((outer(g,e,f),outer(g,f,e))==(u,v),'outer-endpoint inverse')
    ck(boundary(g,v,u)==(f,e),'reversal')
   # Natural under arbitrary tree isomorphism, independent per tree.
   q=list(range(n));rng.shuffle(q);h=[set() for _ in g]
   for u in range(n):h[q[u]]={q[v] for v in g[u]}
   for u,v in non:
    e,f=boundary(g,u,v)
    ck(boundary(h,q[u],q[v])==(edge(q[e[0]],q[e[1]]),edge(q[f[0]],q[f[1]])),'equivariance')
   if n<=4:
    pairs=list(combinations(range(n),2));idx={e:i for i,e in enumerate(pairs)};eidx={e:i for i,e in enumerate(es)};ep=list(combinations(range(len(es)),2));epidx={e:i for i,e in enumerate(ep)}
    for bits in range(1<<len(pairs)):
     orient=sum(((bits>>idx[e])&1)<<i for i,e in enumerate(es));small=0
     for u,v in combinations(range(n),2):
      if v in g[u]:continue
      e,f=boundary(g,u,v);a,b=eidx[e],eidx[f];bit=(bits>>idx[u,v])&1
      if a>b:a,b=b,a;bit^=1
      small|=bit<<epidx[a,b]
     recovered=sum(((orient>>i)&1)<<idx[e] for i,e in enumerate(es))
     for i,(a,b) in enumerate(ep):
      u,v=outer(g,es[a],es[b]),outer(g,es[b],es[a]);bit=(small>>i)&1
      if u>v:u,v=v,u;bit^=1
      recovered|=bit<<idx[u,v]
     ck(recovered==bits,'full tournament codec');roundtrips+=1
   trees+=1
  ck(trees==(1 if n==1 else n**(n-2)),'Cayley universe size');rows.append([n,trees])
 # Hostile to accidental oddness assertion for the oriented tree itself:
 # all arcs from star center 0 to three leaves, S3 acts before completion.
 ck(sum(q[0]==0 for q in permutations(range(4)))==6,'oriented outstar Aut may be even')
 print(json.dumps(dict(status='PASS',checks=C,all_labelled_trees=rows,full_tournament_roundtrips=roundtrips),indent=2))
if __name__=='__main__':main()
