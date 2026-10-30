"""Exact B2 DAG audit, no SAT solver and no assertion statements required."""
from itertools import combinations, product
from collections import Counter
import json

def require(ok,msg):
 if not ok: raise RuntimeError(msg)
def norm(t,mask): return t^mask if t&1 else t
def inputs(n): return tuple(sum(((x>>i)&1)<<x for x in range(1<<n)) for i in range(n))
def defect(n,k,v=None,c=0):
 if v is None:v=(1<<n)-1
 return sum((c^((x&((1<<k)-1)).bit_count()&1)^(x==v))<<x for x in range(1<<n))
def normalized_reach(n,depth):
 mask=(1<<(1<<n))-1; inp=inputs(n)
 states={tuple(sorted(inp))}; mins={t:0 for t in inp}; mins[0]=0; sizes=[]
 for d in range(1,depth+1):
  nxt=set()
  for state in states:
   for a,b in combinations(state,2):
    for raw in [a&b,a&(mask^b),(mask^a)&b,a|b,a^b]:
     c=norm(raw,mask)
     if c==0 or c in state:continue
     mins.setdefault(c,d);nxt.add(tuple(sorted(state+(c,))))
  states=nxt;sizes.append(len(states))
 return mins,sizes

def raw_gate(a,b,t,mask):
 return ((mask^a)&(mask^b) if t&1 else 0)|((mask^a)&b if t&2 else 0)|(a&(mask^b) if t&4 else 0)|(a&b if t&8 else 0)
def raw_reach_three():
 n=3; mask=255; states={inputs(n)}; reached=set(inputs(n));sizes=[]
 for depth in range(1,4):
  nxt=set()
  for state in states:
   outs={raw_gate(a,b,t,mask) for a,b in product(state,repeat=2) for t in range(16)}
   reached.update(outs)
   if depth<3:nxt.update(state+(z,) for z in outs)
  sizes.append(len(nxt) if depth<3 else None);states=nxt
 return reached,sizes

def upper_value(n,k,x):
 # Translation to defect at all ones and final affine constants are handled
 # by truth-table relabeling; here verify the canonical representative.
 if k==0:return int(x==(1<<n)-1)
 a=(x&1)&(1-int((x>>1)==(1<<(n-1))-1))
 return a^(((x>>1)&((1<<(k-1))-1)).bit_count()&1)

def main():
 mins,sizes=normalized_reach(3,3)
 require(sizes==[15,180,2168], 'normalized state counts')
 require(len(mins)==116,'normalized reach count')
 raw,rawsizes=raw_reach_three()
 require({norm(t,255) for t in raw}==set(mins),'independent raw DAG mismatch')
 costs=[mins.get(norm(defect(3,k),255),4) for k in range(4)]
 require(costs==[2,2,3,4],'base complexity')
 # Validate the upper-bound formula on every input and every k through n=10.
 evaluations=0
 for n in range(3,11):
  for k in range(n+1):
   t=defect(n,k)
   for x in range(1<<n):
    require(upper_value(n,k,x)==((t>>x)&1),'construction mismatch')
    evaluations+=1
 # Every translated/complemented three-variable member is covered, not just
 # one representative; input translations and permutations preserve size
 # for functions depending on at least two inputs.
 for k in range(4):
  for v in range(8):
   for c in range(2):
    got=mins.get(norm(defect(3,k,v,c),255),4)
    require(got==costs[k],'translation orbit mismatch')
 # Boundary: two-variable parity with a point defect is one B2 gate.
 require(defect(2,2)==6^8,'two-variable boundary')
 print(json.dumps({'status':'FINITE-EXACT','basis':'all 16 binary Boolean gates','normalized_state_counts':sizes,'normalized_functions_through_three_gates':len(mins),'raw_ordered_state_counts':rawsizes,'raw_functions_through_three_gates':len(raw),'three_variable_costs_by_affine_support':costs,'construction_evaluations_n3_through_n10':evaluations,'independent_enumerators_agree':True},indent=2,sort_keys=True))
if __name__=='__main__':main()
