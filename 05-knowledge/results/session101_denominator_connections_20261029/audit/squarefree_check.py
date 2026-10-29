#!/usr/bin/env python3
"""Independent modular audit; no imports from producer scripts."""
import json
N=10201
checks=0
def check(b):
 global checks
 checks+=1
 if not b: raise RuntimeError('audit failed')

def orbit(C,m,start=0):
 out=[];seen=set();x=start
 while x not in seen:
  seen.add(x);out.append(x);x=(27*x+C)*pow(128,-1,m)%m
 return out

out={}
for C,expected,limit in [(19,[0,3795,5214],4985),(31,[0,8601,8780],8600)]:
 cyc=orbit(C,N);check(len(cyc)==N)
 inv={x:i for i,x in enumerate(cyc)}
 # Prefix numerator coefficients independent of halving denominators, which are units.
 anchors=sorted([inv[0],inv[-pow(3,-1,N)%N],inv[-5*pow(9,-1,N)%N]])
 check(anchors==expected)
 alpha=27*pow(128,-1,N)%N;beta=C*pow(128,-1,N)%N;c=(alpha-1)//101
 for i,z in enumerate(cyc):
  check(z==beta*(i+101*c*i*(i-1)//2)%N)
  y=z*pow(beta,-1,N)%N
  check(i==(y-101*c*y*(y-1)//2)%N)
 def count(r):
  forbidden={(-t)%N for t in range(r+1)}
  for a in anchors[1:]:forbidden.update((a-t)%N for t in range(r))
  return sorted(set(range(N))-forbidden)
 gaps=[anchors[1],anchors[2]-anchors[1],N-anchors[2]-1]
 for r in list(range(21))+[1000,limit-1,limit,limit+1,N]:
  check(len(count(r))==sum(max(g-r,0) for g in gaps))
 survivors=count(limit);check(len(survivors)==1);check(not count(limit+1))
 out[str(C)]={'anchors':anchors,'last_possible_blocks':limit,'unique_last_source':cyc[survivors[0]],'gap_formula_terms':gaps}

# All iterates modulo49 from fixed point plus7; direct phase arithmetic.
f=31*pow(101,-1,49)%49
cyc=orbit(31,49,(f+7)%49)
check(len(cyc)==6)
for x in cyc:
 phases=[x,(3*x+1)*pow(2,-1,49)%49,(9*x+5)*pow(16,-1,49)%49]
 check(all(y!=0 for y in phases))
 check(phases[0]%7!=0 and phases[1]%7!=0 and phases[2]%7==0)
out['seven_repair_cycle']=cyc
out['checks']=checks
print(json.dumps(out,indent=2,sort_keys=True))
