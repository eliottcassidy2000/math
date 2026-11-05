"""Uniform finite control for strong-C<=2 outputs of the physical two-port language."""
from twoport_core import shapes,transitions,brute,H,physical,step
from collections import deque
from itertools import permutations,combinations
import json


def need(x,s):
 if not x:raise RuntimeError(s)
def cycle_readout(v,S):
 return sum(x for x,shape in zip(v,shapes(2))
            if len(shape)==1 and shape[0][0]<2 and S>>shape[0][0]&1
            and (shape[0][1]==2 or not S>>shape[0][1]&1))
def capped_step(v,S,d,cap):
 out=[0]*27
 for i,x in enumerate(v):
  if x:
   for j,m in transitions(2,S,d)[i]:out[j]=min(cap,out[j]+m*x)
 return tuple(out)

def machine():
 seed=brute([2,0],[0,1]);states=[seed];seen={seed:0};edges=[]
 for v in states:
  row=[]
  for S in range(4):
   C=cycle_readout(v,S)
   for d in range(3):
    if C>2:row.append([-1,3]);continue
    w=capped_step(v,S,d,3)
    if w not in seen:seen[w]=len(states);states.append(w)
    row.append([seen[w],C])
  edges.append(row)
 # Add an absorbing rejecting state. Edge output1 means new graph is strong C1/2.
 reject=len(states)
 E=[[(q if q>=0 else reject,int(C in (1,2))) for q,C in row] for row in edges]+[[(reject,0)]*12]
 cls=[0]*len(E);rounds=0
 while True:
  keys={};new=[]
  for row in E:
   key=tuple((a,cls[q]) for q,a in row)
   if key not in keys:keys[key]=len(keys)
   new.append(keys[key])
  if new==cls:break
  cls=new;rounds+=1
 quotient=[None]*len(set(cls))
 for p,row in enumerate(E):
  result=[[cls[q],a] for q,a in row]
  if quotient[cls[p]] is None:quotient[cls[p]]=result
  else:need(quotient[cls[p]]==result,'unstable quotient')
 need(len(states)==1369 and len(quotient)==39,'frozen finite-control dimensions')
 return dict(cap=3,raw_states=[list(v) for v in states],raw_edges=edges,
             class_of_state=cls,quotient=quotient,initial_class=cls[0],reject_class=cls[reject],refinement_rounds=rounds)

def bounded_python(X):
 seed=brute([2,0],[0,1]);states=[seed];seen={seed};values=set();edges=reject=0
 for v in states:
  for S in range(4):
   C=cycle_readout(v,S)
   for d in range(3):
    if C>2:reject+=1;continue
    w=capped_step(v,S,d,X+1)
    if H(w,2)>X:continue
    edges+=1
    if C:values.add(H(w,2))
    if w not in seen:seen.add(w);states.append(w)
 return dict(bound=X,states=len(states),admissible_edges=edges,cycle_rejected_edges=reject,support=sorted(values))

if __name__=='__main__':print(json.dumps(machine(),separators=(',',':'),sort_keys=True))
