#!/usr/bin/env python3
from itertools import product
from pathlib import Path
import json

checks=0

def gate(q):
 global checks
 if not q:raise RuntimeError('Exact count gate failed')
 checks+=1

def formula(p,m,D):
 N=p**(m-1);q,t=divmod(D+1,N)
 return 1+t*(p**(q+1)-1)+(N-t)*(p**q-1)

for p in (2,3,5):
 for m in range(1,6):
  for D in range(7):
   N=p**(m-1)
   sectors=sum(len({i%N for i,a in enumerate(v) if a})<=1
               for v in product(range(p),repeat=D+1))
   gate(sectors==formula(p,m,D))

# Receipts from the frozen complete lift universes; the exact totals are pinned there.
for p,m,D,count in ((2,2,3,7),(2,3,3,5),(2,3,4,7),(2,4,2,4),
                    (3,2,3,13),(3,3,2,7),(5,2,2,13),(2,4,3,5)):
 gate(formula(p,m,D)==count)

# The fixed compatible lifts 1 and X^2 add to the wrong lift modulo4.
a=[1,0,1];sq=[0]*5
for i,x in enumerate(a):
 for j,y in enumerate(a):sq[i+j]=(sq[i+j]+x*y)%4
phi=[1,0,0,0,1]
gate(sq!=phi)
b=[1,2,1];sq=[0]*5
for i,x in enumerate(b):
 for j,y in enumerate(b):sq[i+j]=(sq[i+j]+x*y)%4
phi=[1,0,2,0,1]
gate(sq==phi)
# A fresh dense-convolution path counts compatible lifts with multiplicities.
def conv(a,b,q):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]=(c[i+j]+x*y)%q
 return c

universes=[]
for p,m,D in ((2,3,4),(3,2,3),(5,2,2),(2,4,3)):
 q=p**m;fibres={};total=0
 for a in product(range(q),repeat=D+1):
  b=[1]
  for _ in range(p):b=conv(b,a,q)
  c=[0]*(p*D+1)
  for j,x in enumerate(a):c[p*j]=x
  ok=b==c
  f=tuple(x%p for x in a)
  gate(not ok or len({i%(p**(m-1)) for i,x in enumerate(f) if x})<=1)
  if ok:fibres[f]=fibres.get(f,0)+1
  total+=1
 expected={f for f in product(range(p),repeat=D+1)
           if len({i%(p**(m-1)) for i,x in enumerate(f) if x})<=1}
 gate(set(fibres)==expected)
 for multiplicity in fibres.values():gate(multiplicity==1)
 gate(sum(fibres.values())==formula(p,m,D))
 universes.append({'p':p,'m':m,'D':D,'complete_coefficient_lifts':total,
                   'compatible_lifts':sum(fibres.values()),
                   'admissible_reductions':len(fibres),'maximum_fibre_size':max(fibres.values())})
print(json.dumps({'status':'PASS','exact_gates':checks,
                  'sector_universe':{'p':[2,3,5],'m':'1..5','D':'0..6'},
                  'stored_complete_lift_boxes_checked':8,
                  'fresh_multiplicity_universes':universes,
                  'fresh_complete_lifts':sum(u['complete_coefficient_lifts'] for u in universes)},
                 sort_keys=True,indent=2))
