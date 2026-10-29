#!/usr/bin/env python3
"""Independent structural decode and AHU symmetry audit of producer graph codecs.
Usage: python3 independent_sparse_audit.py [producer.py]
"""
import importlib.util,itertools,json,math,random,sys
from pathlib import Path
producer=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parents[1]/'trijection'/'check_sparse_carrier.py'
spec=importlib.util.spec_from_file_location('producer',producer);p=importlib.util.module_from_spec(spec);spec.loader.exec_module(p)
N=0
def ck(x,msg):
 global N
 N+=1
 if not x:raise RuntimeError(msg)
def comps(g):
 unseen=set(range(len(g)));out=[]
 while unseen:
  block={unseen.pop()};todo=list(block)
  while todo:
   v=todo.pop()
   for w in g[v]&unseen:unseen.remove(w);block.add(w);todo.append(w)
  out.append(block)
 return out

def tree_signature(g):
 n=len(g);alive=set(range(n));deg=list(map(len,g))
 while len(alive)>2:
  leaves=[v for v in alive if deg[v]<=1]
  ck(bool(leaves),'tree has leaves')
  alive.difference_update(leaves)
  for v in leaves:
   for w in g[v]:deg[w]-=1
 def rooted(v,parent):
  child=[rooted(w,v) for w in g[v] if w!=parent]
  sig=tuple(sorted(a for a,b in child)); aut=math.prod(b for a,b in child)
  for a in set(sig):aut*=math.factorial(sig.count(a))
  return sig,aut
 if len(alive)==1:
  return ('vertex',rooted(next(iter(alive)),-1)[0]),rooted(next(iter(alive)),-1)[1]
 u,v=alive;a,aa=rooted(u,v);b,bb=rooted(v,u)
 return ('edge',tuple(sorted((a,b)))),aa*bb*(2 if a==b else 1)

def decode(g,kind):
 if len(g)==1:return (0,()),1
 core={v for v in range(len(g)) if len(g[v])>=4};roles={}
 for v in core:
  lengths=[]
  for w in g[v]-core:
   prev=v;length=1
   while len(g[w])==2:
    nxt=next(iter(g[w]-{prev}));prev,w=w,nxt;length+=1
   ck(len(g[w])==1,'tag terminates at leaf')
   lengths.append(length)
  roles[v]=tuple(sorted(lengths))
 if kind=='tournament':
  c=sorted(v for v in core if roles[v]==(1,2,3));n=len(c)
  tail=[next(w for w in g[v]&core if roles[w]==(1,2,4)) for v in c]
  head=[next(w for w in g[v]&core if roles[w]==(1,2,5)) for v in c]
  adj=tuple(tuple(head[j] in g[tail[i]] for j in range(n)) for i in range(n))
 else:
  c=sorted(v for v in core if roles[v]==(1,2,3,4));n=len(c)
  adj=tuple(tuple(c[j] in g[c[i]] for j in range(n)) for i in range(n))
 aut=sum(all(adj[i][j]==adj[q[i]][q[j]] for i in range(n) for j in range(n)) for q in itertools.permutations(range(n)))
 canon=min(tuple(adj[q[i]][q[j]] for i in range(n) for j in range(n)) for q in itertools.permutations(range(n)))
 return (n,canon),aut

def main():
 rng=random.Random(5029);rows={}
 for kind,enc in [('tournament',p.encode_tournament),('graph',p.encode_graph)]:
  counts=[]
  for n in range(5):
   classes={};seenaut={}
   for mask in range(1<<(n*(n-1)//2)):
    g=enc(n,mask);key,aut=decode(g,kind);classes[key]=1;seenaut[key]=aut
    ck(sum(map(len,g))==2*(len(g)-1),'Euler edge count')
    q=list(range(len(g)));rng.shuffle(q);h=[set() for _ in g]
    for v in range(len(g)):h[q[v]]={q[w] for w in g[v]}
    ck(decode(h,kind)==(key,aut),'arbitrary labels preserve decode')
   counts.append(len(classes))
  rows[kind]=counts
 ck(rows['tournament']==[1,1,1,2,4],'tournament isoclasses')
 ck(rows['graph']==[1,1,2,4,11],'graph isoclasses')
 for k in range(1,129):
  codes=set()
  for j in range(k):
   g=p.padding_tree(k,j);code,aut=tree_signature(g)
   ck(aut==1,'padding rigid by tree automorphism recursion')
   ck(code not in codes,'padding types distinct');codes.add(code)
   ck(max(map(len,g))==3,'max degree three')
 # Same-size sparse graphs at n3: all3 labelled two-edge masks form one orbit.
 n=3;ss=set()
 for edges in itertools.combinations(range(3),2):
  g=[set() for _ in range(n)]
  for e in edges:
   a,b=list(itertools.combinations(range(n),2))[e];g[a].add(b);g[b].add(a)
  code,aut=tree_signature(g);ss.add(code)
 ck(len(ss)==1,'one sparse3 type versus two tournament types')
 print(json.dumps(dict(status='PASS',checks=N,decoded_class_counts=rows,padding_families_through=128,padding_trees_checked=sum(range(1,129))),indent=2))
if __name__=='__main__':main()
