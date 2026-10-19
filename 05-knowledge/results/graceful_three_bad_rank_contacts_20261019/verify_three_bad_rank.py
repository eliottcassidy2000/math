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

def transform_two_bad(source,target,k,m,F):
 global moves_checked
 bad=[d for d in range(1,k+1)if not contact_connected(d,F,m)]
 check(len(bad)<=2,'at most two bad');good=set(range(1,k+1))-set(bad);t=source.copy();steps=0
 def commit(new,protected):
  nonlocal t,steps
  global moves_checked
  changed=[d for d in t if t[d]!=new[d]]
  if not changed:return
  check(len(changed)<=2 and set(changed)<=set(range(1,k+1)),'width and colors');check(protected<=set(new.values()),'protected set')
  validate(new,k,m,F);t=new;steps+=1;moves_checked+=1
 def install(e,g,protected):
  if t[e]==g:return
  graph={v:[]for v in range(m+1)}
  for color,(a,b)in t.items():graph[a].append((b,color));graph[b].append((a,color))
  start,end=g;prev={start:None};todo=[start]
  for v in todo:
   for w,c in graph[v]:
    if w not in prev:prev[w]=(v,c);todo.append(w)
  path=[];v=end
  while v!=start:u,c=prev[v];path.append(c);v=u
  new=t.copy()
  if e in path:new[e]=g
  else:
   opts=[c for c in path if c in good and t[c]not in protected];check(bool(opts),'connected cycle repair edge');f=opts[0]
   rem=[x for c,x in t.items()if c not in(e,f)]+[g];cs=components(rem,range(m+1));check(len(cs)==2,'two-component repair');S=cs[0]
   opts=[(a,a+f)for a in range(m-f+1)if(a in S)!=((a+f)in S)];check(bool(opts),'connected crossing contact');new[e]=g;new[f]=opts[0]
  commit(new,protected)
 if len(bad)==2:
  a,b=bad;cs=components(F.values(),range(m+1));part={v:i for i,c in enumerate(cs)for v in c}
  def typ(e):return tuple(sorted((part[e[0]],part[e[1]])))
  A={};B={}
  for color,bank in((a,A),(b,B)):
   for v in range(m-color+1):
    g=(v,v+color);z=typ(g)
    if z[0]!=z[1]:bank.setdefault(z,g)
  x,y,u,v=typ(t[a]),typ(t[b]),typ(target[a]),typ(target[b]);plan=[]
  check(x!=y and u!=v,'acyclic partial choices')
  if u!=y:plan=[(a,target[a]),(b,target[b])]
  elif v!=x:plan=[(b,target[b]),(a,target[a])]
  elif set(A)-{x,y}:
   w=min(set(A)-{x,y});plan=[(a,A[w]),(b,target[b]),(a,target[a])]
  elif set(B)-{x,y}:
   w=min(set(B)-{x,y});plan=[(b,B[w]),(a,target[a]),(b,target[b])]
  else:
   new=t.copy();new[a]=target[a];new[b]=target[b];commit(new,set(F.values()))
  for e,g in plan:
   other=b if e==a else a
   check(typ(g)!=typ(t[other]),'proposed bad forest acyclic');install(e,g,set(F.values())|{t[other]})
  check(t[a]==target[a] and t[b]==target[b] and steps<=3,'bad phase')
 else:
  for e in bad:install(e,target[e],set(F.values()))
 for e in sorted(good):install(e,target[e],set(t.values())&set(target.values()))
 check(t==target and steps<=k+1,'target and total bound');return steps

branches={'direct_single':0,'common_A':0,'parallel_swap':0,'BC_type_swap':0}
def transform_three_bad(source,target,k,m,F):
 global moves_checked
 bad=[d for d in range(1,k+1)if not contact_connected(d,F,m)]
 if len(bad)<=2:return transform_two_bad(source,target,k,m,F)
 check(len(bad)==3,'three bad');rank=lambda d: len(components(F.values(),range(m+1)))-len(components(list(F.values())+[(a,a+d)for a in range(m-d+1)],range(m+1)))
 A=next((d for d in bad if rank(d)>=3),None);check(A is not None,'rank-three color');B,C=[d for d in bad if d!=A];good=set(range(1,k+1))-set(bad);t=source.copy();steps=0
 def commit(new,protected):
  nonlocal t,steps
  global moves_checked
  changed=[d for d in t if t[d]!=new[d]]
  if not changed:return
  check(len(changed)<=2 and set(changed)<=set(range(1,k+1)),'width/colors');check(protected<=set(new.values()),'protected')
  validate(new,k,m,F);t=new;steps+=1;moves_checked+=1
 def install(e,g,protected):
  if t[e]==g:return
  graph={v:[]for v in range(m+1)}
  for c,(u,v)in t.items():graph[u].append((v,c));graph[v].append((u,c))
  start,end=g;prev={start:None};todo=[start]
  for v in todo:
   for w,c in graph[v]:
    if w not in prev:prev[w]=(v,c);todo.append(w)
  path=[];v=end
  while v!=start:u,c=prev[v];path.append(c);v=u
  new=t.copy()
  if e in path:new[e]=g
  else:
   opts=[c for c in path if c in good and t[c]not in protected];check(bool(opts),'good cycle edge');f=opts[0]
   rem=[x for c,x in t.items()if c not in(e,f)]+[g];cs=components(rem,range(m+1));check(len(cs)==2,'repair cut');S=cs[0]
   opts=[(a,a+f)for a in range(m-f+1)if(a in S)!=((a+f)in S)];check(bool(opts),'repair contact');new[e]=g;new[f]=opts[0]
  commit(new,protected)
 def single_bad(e,g):
  candidate={d:(g if d==e else t[d])for d in bad};check(len(components(list(F.values())+list(candidate.values()),range(m+1)))==k-2,'proposed bad forest')
  install(e,g,set(F.values())|{t[d]for d in bad if d!=e})
 def type_map(es):
  cs=components(es,range(m+1));p={v:i for i,c in enumerate(cs)for v in c}
  return lambda e:tuple(sorted((p[e[0]],p[e[1]])))
 def change(e,g):
  if t[e]==g:return
  other=C if e==B else B;typ=type_map(list(F.values())+[t[other]]);x,y=typ(t[e]),typ(g);a=typ(t[A]);check(x[0]!=x[1] and y[0]!=y[1],'nonloops')
  if x==y or a!=y:branches['direct_single']+=1;single_bad(e,g);return
  bank=[(v,v+A)for v in range(m-A+1)];common=[z for z in bank if typ(z)[0]!=typ(z)[1] and typ(z)not in(x,y)]
  if common:branches['common_A']+=1;single_bad(A,common[0]);single_bad(e,g);return
  candidates=[z for z in bank if typ(z)==x];check(bool(candidates),'parallel A alternative');branches['parallel_swap']+=1
  before={frozenset(c)for c in components(list(F.values())+[t[d]for d in bad],range(m+1))};new=t.copy();new[A]=candidates[0];new[e]=g
  after={frozenset(c)for c in components(list(F.values())+[new[d]for d in bad],range(m+1))};check(before==after,'swap closure');commit(new,set(F.values())|{t[other]})
 typ=type_map(F.values());banks={}
 for d in(B,C):
  bank={}
  for v in range(m-d+1):
   e=(v,v+d);z=typ(e)
   if z[0]!=z[1]:bank.setdefault(z,e)
  banks[d]=bank
 x,y,u,v=typ(t[B]),typ(t[C]),typ(target[B]),typ(target[C]);check(x!=y and u!=v,'partial forest');plan=[]
 if u!=y:plan=[(B,target[B]),(C,target[C])]
 elif v!=x:plan=[(C,target[C]),(B,target[B])]
 elif set(banks[B])-{x,y}:
  w=min(set(banks[B])-{x,y});plan=[(B,banks[B][w]),(C,target[C]),(B,target[B])]
 elif set(banks[C])-{x,y}:
  w=min(set(banks[C])-{x,y});plan=[(C,banks[C][w]),(B,target[B]),(C,target[C])]
 else:
  branches['BC_type_swap']+=1;new=t.copy();new[B]=target[B];new[C]=target[C];commit(new,set(F.values())|{t[A]})
 for e,g in plan:change(e,g)
 check(t[B]==target[B]and t[C]==target[C]and steps<=6,'BC phase');single_bad(A,target[A]);check(steps<=7,'bad phase')
 for e in sorted(good):install(e,target[e],set(t.values())&set(target.values()))
 check(t==target and steps<=k+4,'target and bound');return steps

from itertools import product
boundaries=[]
for m in(5,6):
 F={6:(0,6)}if m==6 else{};ts=enumerate_fiber(5,m,F);keys=[tuple(t[d]for d in range(1,m+1))for t in ts];adj=[[j for j,b in enumerate(keys)if sum(x!=y for x,y in zip(a,b))<=2]for a in keys];diam=0
 for i in range(len(ts)):
  dist={i:0};todo=[i]
  for u in todo:
   for v in adj[u]:
    if v not in dist:dist[v]=dist[u]+1;todo.append(v)
  check(len(dist)==len(ts),'boundary connected');diam=max(diam,max(dist.values()))
 boundaries.append({'m':m,'states':len(ts),'diameter':diam})
reports=[]
for m in(7,8):
 tails=states=paths=mx=0
 for starts in product(*(range(m-d+1)for d in range(6,m+1))):
  F={d:(a,a+d)for d,a in zip(range(6,m+1),starts)}
  if len(components(F.values(),range(m+1)))!=6:continue
  check(contact_connected(1,F,m)and contact_connected(2,F,m),'good colors')
  bad=[d for d in range(1,6)if not contact_connected(d,F,m)]
  if len(bad)!=3:continue
  check(6-len(components(list(F.values())+[(a,a+3)for a in range(m-2)],range(m+1)))>=3,'residue rank')
  ts=enumerate_fiber(5,m,F);tails+=1;states+=len(ts);targets=[ts[0],ts[len(ts)//2],ts[-1]]
  for source in ts:
   for target in targets:mx=max(mx,transform_three_bad(source,target,5,m,F));paths+=1
 reports.append({'m':m,'three_bad_tails':tails,'all_source_states':states,'source_target_checks':paths,'max_moves':mx,'proved_bound':9})
# Abstract rank-two hostile: two bad forests, requiring all three colors to change.
star=[[(0,1),(0,2)],[(0,2),(0,3)],[(0,3),(0,1)]];valid=[]
for choice in product(*star):
 if len(components(choice,range(4)))==1 and len(set(choice))==3:valid.append(choice)
check(len(valid)==2 and sum(x!=y for x,y in zip(*valid))==3,'generic cyclic-star obstruction')
# Abstract parallel-contact branch: colors A/B exchange two parallel classes
# after contracting C. Colored contacts may be parallel in an F quotient.
A_edges=[(0,1),(0,2),(0,3)];C_edge=(0,1)
old=[(0,3),(0,2),C_edge];new=[(0,2),(0,3),C_edge]
check(6-len(components(A_edges,range(6)))==3,'abstract A rank three')
check({frozenset(c)for c in components(old,range(6))}=={frozenset(c)for c in components(new,range(6))},'abstract swap closure')
check(len(components(old+[(0,4),(0,5)],range(6)))==1 and len(components(new+[(0,4),(0,5)],range(6)))==1,'abstract full-tree swap')
check(sum(x!=y for x,y in zip(old,new))==2,'abstract two-color width')
print(json.dumps({'status':'passed','checks':checks,'moves_replayed':moves_checked,'boundaries':boundaries,'three_bad_controls':reports,'repair_branches':branches,'generic_rank_two_bad_selections':len(valid),'abstract_parallel_swap_geometry':'passed'},indent=2,sort_keys=True))
