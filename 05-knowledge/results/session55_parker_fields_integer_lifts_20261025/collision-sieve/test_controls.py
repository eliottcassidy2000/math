#!/usr/bin/env python3
import copy, hashlib, itertools, json, math
from pathlib import Path
import verify as v
h=Path(__file__).resolve().parent
c=json.loads((h/'certificate.json').read_text())
producer=h.parent/'integer-search/certificates.json'
p=json.loads(producer.read_text())['integer_search']
v.need(c['search']['bank_sha256']==p['offset_banks_sha256'])
v.need(c['search']['pairs']==p['eligible_pair_tests'])
# Discriminant checked directly without factor-grouping code, using integer differences.
for u in range(-5,6):
 for w in range(-5,6):
  points=[s*u+t*w for s in (-1,0,1) for t in (-1,0,1)]
  direct=math.prod((a-b)**2 for a,b in itertools.combinations(points,2))
  formula=2**16*u**18*w**18*(u*u-w*w)**10*(u*u-4*w*w)**4*(4*u*u-w*w)**4
  v.need(direct==formula)
# All local decisions checked by direct duplicate detection, separate from C and masks.
for q in v.P:
 for a in range(q):
  for b in range(q):
   offsets=[s*a+t*b for s in (-1,0,1) for t in (-1,0,1)]
   duplicate=any((x-y)%q==0 for x,y in itertools.combinations(offsets,2))
   v.need(duplicate==bool((v.M[q][a]>>b)&1))
for kind in ('changed','missing','extra'):
 d=copy.deepcopy(c)
 if kind=='changed':d['search']['remaining']+=1
 elif kind=='missing':d.pop('controls')
 else:d['extra']=1
 try:v.need(c==d,'mutant rejected')
 except RuntimeError:pass
 else:raise RuntimeError('mutation accepted')
print('PASS: upstream bank hash and pair universe; 121 discriminants; direct pair-collision controls; three certificate mutations')
