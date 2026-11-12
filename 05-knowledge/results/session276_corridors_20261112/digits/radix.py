#!/usr/bin/env python3
"""Exact positive radix inverse for a physical central sandwich family.
Not a Chan--Pak solver. Reports all solutions of the specified polynomial.
"""
from dataclasses import dataclass
import json

def g(a):
 if not a:return 0
 q,r=divmod(a,2)
 return g(q)+(q+(q&1) if r==0 else 0)
@dataclass(frozen=True)
class State:
 c:int;h:int;q:int;u:int;v:int;e:int;f:int
 def value(self,A,B):return self.c+self.h*(g(A)+g(B))+self.q*A*B+self.u*A+self.v*B+self.e*(A&1)+self.f*(B&1)
 def split(self,r,s):
  return State(self.c+self.u*r+self.v*s+self.q*r*s+self.e*r+self.f*s,self.h,4*self.q,2*self.u+2*self.q*s+self.h*(1-r),2*self.v+2*self.q*r+self.h*(1-s),self.h*(1-r),self.h*(1-s))

def inverse(N,state):
 if min(state.h,state.q,state.u,state.v)<1 or min(state.c,state.e,state.f)<0:raise ValueError('positive family required')
 stats={'nodes':0,'max_depth':0,'constant_rejections':0,'forced_zero_leaves':0}
 out=[]
 def dfs(S,lowA,lowB,k):
  stats['nodes']+=1;stats['max_depth']=max(stats['max_depth'],k)
  budget=N-S.c
  if budget<0:stats['constant_rejections']+=1;return
  if budget==0:out.append((lowA,lowB));return
  if S.u>budget and S.v>budget:stats['forced_zero_leaves']+=1;return
  rr=[0] if S.u>budget else [0,1];ss=[0] if S.v>budget else [0,1]
  for r in rr:
   for s in ss:dfs(S.split(r,s),lowA+(r<<k),lowB+(s<<k),k+1)
 dfs(state,0,0,0)
 return sorted(out),stats

C5=State(41,15,120,42,63,6,9)
C9=State(15381,3357,50004,25915,20732,1865,1492)

def run():
 checks=0
 for state in [C5,C9]:
  for a in range(16):
   for b in range(16):
    r=a&1;s=b&1
    if state.value(a,b)!=state.split(r,s).value(a//2,b//2):raise RuntimeError('split mismatch')
    checks+=1
  for N in list(range(state.c,state.c+250))+[state.value(a,b) for a,b in [(0,9),(7,3),(12,13),(33,27)]]:
   actual,stats=inverse(N,state)
   budget=N-state.c
   expected=[]
   for A in range(max(-1,budget//state.u)+1):
    for B in range(max(-1,budget//(state.v+state.q*A))+1):
     if state.value(A,B)==N:expected.append((A,B))
   if actual!=expected:raise RuntimeError(('inverse mismatch',N,actual,expected))
   if budget>=0 and stats['max_depth']>(budget//min(state.u,state.v)+1).bit_length():raise RuntimeError('depth bound')
   checks+=1
 examples=[]
 for S in [C5,C9]:
  for A,B in [(0,0),(7,3),(33,27),(63,62)]:
   N=S.value(A,B);sol,stats=inverse(N,S)
   examples.append(dict(core=S.c,N=N,solutions=sol,stats=stats))
 print(json.dumps(dict(status='PASS',checks=checks,examples=examples,scope='exact finite decision tree for fixed physical central sandwich, not universal spectrum coverage'),indent=2))
if __name__=='__main__':run()
