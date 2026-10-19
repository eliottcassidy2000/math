#!/usr/bin/env python3
"""Independent exact move replay for the whole hostile fibers k2,3,4."""
import json
from collections import deque
checks=0;moves_checked=0

def check(x,msg):
 global checks
 checks+=1
 if not x:raise RuntimeError(msg)
def edge(x,y):return tuple(sorted((x,y)))
def frozen(k,m):return {j:edge(0,j) for j in range(k+1,m+1)if j!=m-1}|{m-1:edge(1,m)}
def components(edges,vertices):
 g={v:[]for v in vertices}
 for a,b in edges:
  if a in g and b in g:g[a].append(b);g[b].append(a)
 unseen=set(g);out=[]
 while unseen:
  s={unseen.pop()};q=list(s)
  while q:
   for v in g[q.pop()]:
    if v in unseen:unseen.remove(v);s.add(v);q.append(v)
  out.append(s)
 return out
def validate(t,k,m,F):
 check(set(t)==set(range(1,m+1)),'one edge per color')
 check(all(0<=a<b<=m and b-a==d for d,(a,b)in t.items()),'labels and differences')
 check(all(t[d]==e for d,e in F.items()),'frozen tail')
 check(len(components(t.values(),range(m+1)))==1,'tree connectivity')
def enumerate_fiber(k,m):
 F=frozen(k,m);out=[]
 def rec(d,t):
  if d>k:out.append(t.copy());return
  comps=components(t.values(),range(m+1));part={v:i for i,c in enumerate(comps)for v in c}
  for a in range(m-d+1):
   if part[a]!=part[a+d]:t[d]=(a,a+d);rec(d+1,t);del t[d]
 rec(1,F.copy());return out

def transform(source,target,k,m):
 global moves_checked
 F=frozen(k,m);t=source.copy();b=m-1;H=next(c for c in components(F.values(),range(m+1))if 0 in c);protected=set(F.values());steps=0
 def move(new):
  nonlocal t,steps
  global moves_checked
  changed=[d for d in t if t[d]!=new[d]]
  check(1<=len(changed)<=2 and all(d<=k for d in changed),'move width and frozen colors')
  check(protected<=set(new.values()),'protected targets')
  validate(new,k,m,F);t=new;steps+=1;moves_checked+=1
 def attachment():
  a=[d for d,e in t.items()if b in e];check(len(a)==1,'buffer leaf');return a[0]
 for e,g in sorted(target.items()):
  if e>k or b in g:continue
  if t[e]==g:protected.add(g);continue
  d=attachment()
  if d!=e:
   rem=[x for color,x in t.items()if color not in(d,e)];cs=components(rem,[v for v in range(m+1)if v!=b]);check(len(cs)==2,'attachment cut')
   S=next(c for c in cs if not(c&H));check(S and S<=set(range(2,k+1)),'small cut side');z=max(S)
   new=t.copy();new[d]=edge(z,z+d);new[e]=edge(b-e,b);move(new)
  check(attachment()==e,'desired buffer color')
  # Recover the old path between the target endpoints without the buffer.
  graph={v:[]for v in range(m+1)if v!=b}
  for color,(a,c)in t.items():
   if b not in(a,c):graph[a].append((c,color));graph[c].append((a,color))
  start,end=g;prev={start:None};q=deque([start])
  while q and end not in prev:
   v=q.popleft()
   for w,color in graph[v]:
    if w not in prev:prev[w]=(v,color);q.append(w)
  check(end in prev,'target fundamental path');path=[];v=end
  while v!=start:v0,color=prev[v];path.append(color);v=v0
  choices=[c for c in path if t[c]not in protected and t[c]not in set(target.values())]
  # A target-tree cycle is impossible, so an edge outside the entire target exists.
  check(bool(choices),'non-target fundamental-cycle edge');f=choices[0]
  new=t.copy();new[e]=g;new[f]=edge(b-f,b);move(new);protected.add(g)
 d=attachment()
 if t[d]!=target[d]:new=t.copy();new[d]=target[d];move(new)
 check(t==target,'target reached');check(steps<=2*k-1,'diameter bound');return steps
reports=[]
for k in(2,3,4):
 m=2*k+2;F=frozen(k,m);ts=enumerate_fiber(k,m)
 for t in ts:validate(t,k,m,F)
 # Verify each abstract color contact graph, using F contractions.
 b=m-1;vertices=[v for v in range(m+1)if v!=b]
 for d in range(1,k+1):
  es=list(F.values())+[(a,a+d)for a in range(m-d+1)if b not in(a,a+d)]
  check(len(components(es,vertices))==1,'connected color contact graph')
  check(b-d in next(c for c in components(F.values(),range(m+1))if 0 in c),'buffer contact')
 maxsteps=0
 for source in ts:
  for target in ts:maxsteps=max(maxsteps,transform(source,target,k,m))
 reports.append({'k':k,'m':m,'completions':len(ts),'ordered_pairs':len(ts)**2,'max_constructed_moves':maxsteps,'proved_bound':2*k-1})
print(json.dumps({'status':'passed','checks':checks,'moves_replayed':moves_checked,'full_fiber_controls':reports,'scope':'Exact exhaustive controls for k2,3,4; universal theorem follows from the cut/buffer proof. Abstract tree type may change.'},indent=2,sort_keys=True))
