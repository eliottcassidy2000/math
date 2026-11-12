#!/usr/bin/env python3
"""Exact declared cyclic-family experiment; not an iso-class census."""
import json,math
from itertools import combinations
from explore import coefficients

def profile(n,S):
 rows=[sum(1<<((i+j)%n) for j in S) for i in range(n)]
 h,Q,u,v,d,e,b,s,t=coefficients(rows)
 W=[[b[i][j]+b[j][i] for j in range(n)] for i in range(n)];r=(Q-2*h)//n
 records=[];unit_pairs=0
 for k in range(2,(n+1)//2):
  vals=sorted({h+r*k-sum(W[i][j] for i,j in combinations(C,2)) for C in combinations(range(n),k)})
  g=math.gcd(*(a-vals[0] for a in vals)) if len(vals)>1 else 0
  longest=run=1;start=best=vals[0]
  for a,z in zip(vals,vals[1:]):
   if z-a==2:run+=1;unit_pairs+=1
   else:run=1;start=z
   if run>longest:longest=run;best=start
  records.append(dict(k=k,counts=len(vals),gcd=g,run=longest,start=best,lo=vals[0],hi=vals[-1]))
 return dict(n=n,S=S,h=h,weights=sorted(set(W[0][1:])),fibers=records),unit_pairs

out={'carousel':[],'all_cyclic':[]}
for n in range(5,18,2):out['carousel'].append(profile(n,list(range(1,(n+1)//2)))[0])
for n in [9,11,13]:
 hits=0;best=1;examples=[]
 for code in range(1<<((n-1)//2)):
  S=[j if code>>(j-1)&1 else n-j for j in range(1,(n+1)//2)]
  record,count=profile(n,S);hits+=count
  run=max(x['run'] for x in record['fibers']);best=max(best,run)
  if count:examples.append(record)
 out['all_cyclic'].append(dict(n=n,orientations=1<<((n-1)//2),maximum_run=best,unit_pairs_with_presentation_multiplicity=hits,positive_presentations=examples))
print(json.dumps(out,sort_keys=True,indent=2))
