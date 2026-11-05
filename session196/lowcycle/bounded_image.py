"""Exact arbitrary-length bounded image for the two-apex model.
Nonnegative path/cycle transfers permit finite coordinate saturation.
No order cutoff is used. Rejecting H>B prefixes is safe by path insertion.
"""
from apex_transfer import *
from collections import deque
import json,time,sys
B=int(sys.argv[1]) if len(sys.argv)>1 else 67
raw,t,o=build_cycle_dfa();partition,t,o=minimize_dfa(t,o)
live={i for i,c in enumerate(o) if c in (1,2)}
while True:
 new=live|{i for i,row in enumerate(t) if any(j in live for j in row)}
 if new==live:break
 live=new
initial=(partition[0],(1,)+(0,)*26,(1,)+(0,)*26)
states=[initial];index={initial:0};parents=[None];witnesses={};ntrans=0;start=time.time();maxdepth=0;depths=[0]
for idx,(cs,ab,ba) in enumerate(states):
 h=finish(ab,'path_ab')+finish(ba,'path_ba')
 assert h<=B
 if o[cs] in (1,2):witnesses.setdefault(h,idx)
 for a,ns in enumerate(t[cs]):
  if ns not in live:continue
  nab=step(ab,'path_ab',a,B+1);nba=step(ba,'path_ba',a,B+1)
  if finish(nab,'path_ab')+finish(nba,'path_ba')>B:continue
  nxt=(ns,nab,nba);ntrans+=1
  if nxt not in index:
   index[nxt]=len(states);states.append(nxt);parents.append((idx,a));depths.append(depths[idx]+1)
 if idx and idx%100000==0:print('processed',idx,'states',len(states),'elapsed',time.time()-start,flush=True)
def word(i):
 out=[]
 while parents[i] is not None:i,a=parents[i];out.append(a)
 return out[::-1]
missing=[n for n in range(3,B+1,2) if n not in witnesses]
print('B',B,'exact states',len(states),'transitions',ntrans,'max discovery depth',max(depths),'elapsed',time.time()-start)
print('missing strong counts:',missing)
packet={'bound':B,'exact_arbitrary_length':True,'states':len(states),'transitions':ntrans,'max_discovery_depth':max(depths),'missing_strong_odd':missing,'witnesses':{str(h):word(i) for h,i in witnesses.items()}}
json.dump(packet,open(__file__.replace('.py',f'_{B}.json'),'w'),indent=2)
