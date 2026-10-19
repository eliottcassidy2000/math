#!/usr/bin/env python3
"""Standard-library deterministic construction and independent per-move replay."""
import json
from itertools import combinations
checks=0
stats={'single_partial':0,'double_partial':0,'two_single_repair':0,'parallel_swap':0,'full_moves':0,'paths':0,'states':0}
def check(ok,msg):
 global checks
 checks+=1
 if not ok:raise RuntimeError(msg)
def partition(edges,m):
 p=list(range(m+1))
 def root(a):
  while p[a]!=a:a=p[a]
  return a
 for a,b in edges:
  a,b=root(a),root(b)
  if a!=b:p[b]=a
 groups={}
 for a in range(m+1):groups.setdefault(root(a),[]).append(a)
 return tuple(sorted(tuple(v) for v in groups.values()))
def independent(edges,m):return len(partition(edges,m))==m+1-len(edges)
def partial_path(source,target,colors,classes,F,m):
 if not colors:return []
 akey=colors[0];rest=colors[1:];t=source.copy();out=[]
 def commit(new):
  changed=[d for d in colors if t[d]!=new[d]]
  if not changed:return
  check(1<=len(changed)<=2,'partial width')
  check(independent(list(F.values())+list(new.values()),m),'partial independent')
  if len(changed)==2:
   check(partition(list(F.values())+list(t.values()),m)==partition(list(F.values())+list(new.values()),m),'double closure')
   stats['double_partial']+=1
  else:stats['single_partial']+=1
  t.clear();t.update(new);out.append(t.copy())
 tail=partial_path({d:source[d] for d in rest},{d:target[d] for d in rest},rest,classes,F,m)
 for step in tail:
  changed=[d for d in rest if step[d]!=t[d]]
  proposal=t|step
  if len(changed)==2 or independent(list(F.values())+list(proposal.values()),m):commit(proposal);continue
  check(len(changed)==1,'recursive width')
  d=changed[0]
  for anew in classes[akey]:
   trial=proposal|{akey:anew}
   if independent(list(F.values())+list(trial.values()),m):break
  else:raise RuntimeError('rank repair unavailable')
  first=t|{akey:anew}
  if independent(list(F.values())+list(first.values()),m):
   stats['two_single_repair']+=1;commit(first);commit(trial)
  else:stats['parallel_swap']+=1;commit(trial)
 commit(t|{akey:target[akey]})
 check(t==target,'partial exact target')
 check(len(out)<=2**len(colors)-1,'partial bound')
 return out

def validate(t,r,m,F):
 check(set(t)==set(range(1,m+1)),'all colors')
 check(all(0<=a<b<=m and b-a==d for d,(a,b) in t.items()),'difference labels')
 check(all(t[d]==e for d,e in F.items()),'frozen tail')
 check(len(partition(t.values(),m))==1,'full spanning tree')
def transform(source,target,r,m,F):
 classes={d:[(a,a+d) for a in range(m-d+1)] for d in range(1,r+1)}
 validate(source,r,m,F);validate(target,r,m,F)
 t=source.copy();moves=0
 def commit(new):
  nonlocal moves,t
  changed=[d for d in t if t[d]!=new[d]]
  if not changed:return
  check(len(changed)<=2,'full width two');validate(new,r,m,F)
  moves+=1;stats['full_moves']+=1;t=new.copy()
 if r:
  path=partial_path({d:t[d] for d in range(2,r+1)},{d:target[d] for d in range(2,r+1)},list(range(2,r+1)),classes,F,m)
  for part in path:
   new=F|part;comps=partition(new.values(),m)
   check(len(comps)==2,'partial has two components')
   side=set(comps[0])
   def crosses(e):return (e[0] in side)!=(e[1] in side)
   chosen=t[1] if crosses(t[1]) else next(e for e in classes[1] if crosses(e))
   commit(new|{1:chosen})
  commit(t|{1:target[1]})
 check(t==target,'full exact target');check(moves<=(2**(r-1) if r else 0),'full bound');stats['paths']+=1
 return moves

def enumerate_fiber(r,m,F):
 out=[]
 def visit(d,t):
  if d>r:out.append(t.copy());return
  for a in range(m-d+1):
   e=(a,a+d)
   if independent(list(t.values())+[e],m):t[d]=e;visit(d+1,t);del t[d]
 visit(1,F.copy());return out

def main():
 cases=[(m,m,{}) for m in range(8)]
 cases += [(k,2*k+2,{d:(0,d) for d in range(k+1,2*k+3) if d!=2*k+1}|{2*k+1:(1,2*k+2)}) for k in (2,3,4)]
 cases += [(r,7,{d:(0,d) for d in range(r+1,8)}) for r in range(1,7)]
 receipts=[]
 for r,m,F in cases:
  states=enumerate_fiber(r,m,F);stats['states']+=len(states)
  targets=states if m<=5 else [states[i] for i in sorted({0,len(states)//2,len(states)-1})]
  worst=0
  for source in states:
   for target in targets:worst=max(worst,transform(source,target,r,m,F))
  receipts.append({'r':r,'m':m,'fixed':F,'states':len(states),'targets':len(targets),'max_moves':worst})
 print(json.dumps({'checks':checks,'statistics':stats,'cases':receipts},sort_keys=True,indent=2))
if __name__=='__main__':main()
