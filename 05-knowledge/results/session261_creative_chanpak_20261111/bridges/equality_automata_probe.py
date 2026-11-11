"""Validate explicit finite-core Hamilton and domination state machines, not full basis synthesis."""
from collections import defaultdict
from itertools import permutations,product
import random,json
from pathlib import Path
from path_cover_probe import H,sub,graph

def build(word):
 # core indices z=0,a=1,b=2; bit i means spine vertex points to core i.
 n=3+len(word);r=[0]*n;r[0]|=2;r[1]|=4;r[2]|=1
 for k,c in enumerate(word,3):
  for i in range(3):
   if c>>i&1:r[k]|=1<<i
   else:r[i]|=1<<k
  for j in range(3,k):r[j]|=1<<k
 return tuple(r)

def h_automaton(word,core):
 total=0
 def core_arc(i,j):return (i,j) in [(0,1),(1,2),(2,0)]
 for order in permutations(core):
  r=len(core);d={(0,)*(r+1):1}
  for c in word:
   dd=defaultdict(int)
   for state,mass in d.items():
    for i in range(r+1):
     if state[i]==0 and i>0 and c>>order[i-1]&1:continue
     ns=list(state);ns[i]=1 if i==r or c>>order[i]&1 else 2
     dd[tuple(ns)]+=mass
   d=dd
  for state,mass in d.items():
   if 2 in state:continue
   if any(state[i]==0 and not core_arc(order[i-1],order[i]) for i in range(1,r)):continue
   total+=mass
 return total

def d_automaton(word):
 total=0
 for C in range(4):
  sa=bool(C&1);sb=bool(C&2);d={(False,False,bool(C),bool(C)):1}
  for c in word:
   isA=not bool(c&1) # z->v
   adj_b=isA and bool(c&4) # v->b
   adj_a=not isA and not bool(c&2) # a->v
   dd=defaultdict(int)
   for (seen,pending,da,db),mass in d.items():
    # Select current spine vertex.
    if isA:dd[(True,pending,da,db or adj_b)]+=mass
    else:dd[(seen,False,da or adj_a,db)]+=mass
    # Leave current vertex unselected.
    if isA:dd[(seen,pending or not(sb and adj_b),da,db)]+=mass
    elif seen or (sa and adj_a):dd[(seen,pending,da,db)]+=mass
   d=dd
  total+=sum(m for (seen,pending,da,db),m in d.items() if not pending and da and db)
 return total

def main():
 rng=random.Random(2611712);words=[w for n in range(4) for w in product(range(8),repeat=n)]
 words +=[tuple(rng.randrange(8) for _ in range(rng.randrange(4,8))) for _ in range(100)]
 for w in words:
  r=build(w);h=H(r);h0=H(sub(r,list(range(1,len(r)))));dg=graph(r,0)[-1]
  if (h_automaton(w,(0,1,2)),h_automaton(w,(1,2)),d_automaton(w))!=(h,h0,dg):raise RuntimeError(('automaton',w))
  if h-h0-dg+1<0:raise RuntimeError('slack')
 # The necessary equality skeleton's counterexample.
 # U=(a,x,b,c) with x->a, a->{b,c}, b->{c,x}, c->x; root z->a,x and b,c->z.
 r=(12,1,26,18,3);s=H(r)-H(sub(r,[0,1,2,3]))-graph(r,4)[-1]+1
 if s!=2:raise RuntimeError(('skeleton counterexample',s))
 return {'scope':'Every eight-letter word of lengths0..3, plus100 seeded words of lengths4..7; all three automata compared to literal subset Hamilton DP and direct domination enumeration','seed':2611712,'word_cases':len(words),'state_bounds':{'H_triangle':486,'H_deleted_root':54,'domination':64,'constant':1,'signed_slack_total':605},'skeleton_counterexample_slack':s,'basis_synthesis':'NOT RUN','result':'PASS'}
if __name__=='__main__':
 out=main();Path(__file__).with_name('equality_automata_receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
