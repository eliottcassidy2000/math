"""Exact finite-state target compiler. No completeness assertion."""
from probe import make
from collections import deque
import json,sys,time

def compile_target(h,budget=None):
 init,step,val=make(2,[[0,1],[0,0]])
 cap=h+1
 def trunc(v):return tuple(tuple(min(cap,x) for x in row) for row in v)
 queue=deque([init]);parent={init:None};expanded=0
 while queue:
  v=queue.popleft();expanded+=1
  if val(v)==h:
   w=[];u=v
   while parent[u] is not None:
    u,m=parent[u];w.append(m)
   return {'status':'FOUND','target':h,'word':w[::-1],'expanded':expanded,'states':len(parent)}
  if budget is not None and expanded>=budget:return {'status':'INCOMPLETE','target':h,'expanded':expanded,'states':len(parent)}
  if val(v)>h:continue # Hamilton-path count is induced-subtournament monotone.
  for m in range(4):
   u=trunc(step(v,m))
   if u not in parent:parent[u]=(v,m);queue.append(u)
 return {'status':'ABSENT','target':h,'expanded':expanded,'states':len(parent)}
if __name__=='__main__':
 for h in map(int,sys.argv[1:] or [1,3,5,7,9,21,23]):
  t=time.time();r=compile_target(h,200000);r['seconds']=time.time()-t;print(json.dumps(r),flush=True)
