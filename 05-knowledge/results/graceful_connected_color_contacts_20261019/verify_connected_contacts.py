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
def enumerate_fiber(k,m,F=None):
 F=frozen(k,m) if F is None else F;out=[]
 def rec(d,t):
  if d>k:out.append(t.copy());return
  comps=components(t.values(),range(m+1));part={v:i for i,c in enumerate(comps)for v in c}
  for a in range(m-d+1):
   if part[a]!=part[a+d]:t[d]=(a,a+d);rec(d+1,t);del t[d]
 rec(1,F.copy());return out

def contact_connected(d,F,m):
 return len(components(list(F.values())+[(a,a+d)for a in range(m-d+1)],range(m+1)))==1

def transform(source,target,k,m,F):
 global moves_checked
 bad=[d for d in range(1,k+1)if not contact_connected(d,F,m)];check(len(bad)<=1,'contact hypothesis')
 order=bad+[d for d in range(1,k+1)if d not in bad];t=source.copy();steps=0
 for e in order:
  g=target[e]
  if t[e]==g:continue
  graph={v:[]for v in range(m+1)}
  for color,(a,b)in t.items():graph[a].append((b,color));graph[b].append((a,color))
  start,end=g;prev={start:None};q=deque([start])
  while q and end not in prev:
   v=q.popleft()
   for w,c in graph[v]:
    if w not in prev:prev[w]=(v,c);q.append(w)
  path=[];v=end
  while v!=start:u,c=prev[v];path.append(c);v=u
  new=t.copy()
  if e in path:new[e]=g
  else:
   choices=[c for c in path if t[c]not in set(target.values())];check(bool(choices),'non-target cycle edge');f=choices[0]
   check(f<=k and contact_connected(f,F,m),'connected repair color')
   rem=[x for c,x in t.items()if c not in(e,f)]+[g];cs=components(rem,range(m+1));check(len(cs)==2,'repair cut');S=cs[0]
   candidates=[(a,a+f)for a in range(m-f+1)if(a in S)!=((a+f)in S)];check(bool(candidates),'repair contact')
   new[e]=g;new[f]=candidates[0]
  changed=[c for c in t if t[c]!=new[c]];check(1<=len(changed)<=2 and all(c<=k for c in changed),'width')
  check({x for x in t.values()if x in set(target.values())}<=set(new.values()),'monotone target overlap')
  validate(new,k,m,F);t=new;steps+=1;moves_checked+=1
 check(t==target and steps<=k,'target and bound');return steps

from itertools import product
reports=[];boundary=[]
for m in(3,4):
 F={4:(0,4)}if m==4 else{};ts=enumerate_fiber(3,m,F);diam=0
 for t in ts:validate(t,3,m,F)
 for source in ts:
  dist={tuple(source[d]for d in range(1,m+1)):0};queue=[source]
  for u in queue:
   key=tuple(u[d]for d in range(1,m+1))
   for v in ts:
    vk=tuple(v[d]for d in range(1,m+1))
    if vk not in dist and sum(u[d]!=v[d]for d in range(1,m+1))<=2:dist[vk]=dist[key]+1;queue.append(v)
  check(len(dist)==len(ts),'boundary connectivity');diam=max(diam,max(dist.values()))
 boundary.append({'m':m,'states':len(ts),'diameter':diam})
for m in range(5,8):
 fibers=states=pairs=maximum=0
 for starts in product(*(range(m-d+1)for d in range(4,m+1))):
  F={d:(a,a+d)for d,a in zip(range(4,m+1),starts)}
  if len(components(F.values(),range(m+1)))!=4:continue
  check(contact_connected(1,F,m)and contact_connected(2,F,m),'one and two connected')
  ts=enumerate_fiber(3,m,F);fibers+=1;states+=len(ts);pairs+=len(ts)**2
  for t in ts:validate(t,3,m,F)
  for source in ts:
   for target in ts:maximum=max(maximum,transform(source,target,3,m,F))
 reports.append({'m':m,'acyclic_tails':fibers,'states':states,'ordered_pairs':pairs,'max_moves':maximum})
F={4:(0,4),5:(1,6),6:(0,6),7:(0,7)}
check([contact_connected(d,F,7)for d in(1,2,3)]==[True,True,False],'specific one-bad-color example')
hostile=[]
for k in(2,3,4):
 m=2*k+2;F=frozen(k,m);ts=enumerate_fiber(k,m,F);mx=0
 for source in ts:
  for target in ts:mx=max(mx,transform(source,target,k,m,F))
 hostile.append({'k':k,'states':len(ts),'ordered_pairs':len(ts)**2,'max_moves':mx,'bound':k})
print(json.dumps({'status':'passed','checks':checks,'moves_replayed':moves_checked,'boundaries':boundary,'three_color_fibers':reports,'improved_hostile_controls':hostile},indent=2,sort_keys=True))
