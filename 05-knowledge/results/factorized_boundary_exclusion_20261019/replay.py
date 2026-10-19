#!/usr/bin/env python3
"""Independent integer-only replay of one factorized boundary exclusion."""
from itertools import permutations
from collections import Counter
from pathlib import Path
import json,hashlib
ROOT=Path(__file__).resolve().parent
def require(condition,message):
 if not condition:raise RuntimeError(message)
p=(0,28,0,0,12,4,13,93,125)
order=(2,3,0,4,5,6,1,7,8)
pairs=[(a,b) if order.index(a)<order.index(b) else (b,a) for a in range(9) for b in range(a+1,9) if not(p[a]>>b&1 or p[b]>>a&1)]
counts={}; ext=0
for w in permutations(range(9)):
 m=0;valid=True
 for v in w:
  if p[v]&~m:valid=False;break
  m|=1<<v
 if not valid:continue
 ext+=1;rank={v:i for i,v in enumerate(w)}
 h=min(rank[0],rank[1]);k=max(rank[0],rank[1])+1
 values=[3*int(rank[a]<rank[b])-2 for a,b in pairs]
 for i in range(h+1):
  for j in range(k,10):
   key=(w[:i],w[j:]);col=counts.setdefault(key,[0]*len(pairs))
   for t,value in enumerate(values):col[t]+=value
us=sorted({u for u,v in counts},key=lambda w:(len(w),w));vs=sorted({v for u,v in counts},key=lambda w:(len(w),w))
require(ext==540 and len(counts)==175 and len(us)==len(vs)==16,'Core or generator cardinality mismatch')
require(us[1]==(2,) and us[11]==(2,3,4,5),'Named multiplier variables mismatch')
require(pairs==[(0,1),(2,0),(3,0),(0,4),(0,5),(5,1),(6,1),(1,7),(1,8),(2,3),(3,5),(4,5),(4,6),(5,6),(5,7),(7,8)],'Oriented pair mapping mismatch')
file=ROOT/'small_certificate.json'
certificate=json.loads(file.read_text()); polynomial=Counter();positive_constant=False
for term in certificate:
 row=term['row'];m=tuple(term['multiplier']);coef=term['coefficient']
 require(type(coef) is int and coef>0 and type(row) is int and 0<=row<len(pairs),'Invalid term coefficient or row')
 require(m in ((),(0,),(10,)),'Invalid multiplier word index')
 positive_constant |= not m
 for (u,v),col in counts.items():
  mon=tuple(sorted(m+(() if not u else (us.index(u)-1,))+(() if not v else (14+vs.index(v),))))
  polynomial[mon]+=coef*col[row]
require(positive_constant and max(polynomial.values())<=0,'Certificate polynomial has a positive coefficient or no constant multiplier')
require(len(certificate)==20,'Unexpected certificate term count')
require(sum(v<0 for v in polynomial.values())==486 and sum(v==0 for v in polynomial.values())==9,'Expanded coefficient counts mismatch')
report={'status':'PASS','scope':'All nonnegative factorized boundary weights for the specified core and majority order','core':p,'order':order,'extensions':ext,'generators':len(counts),'prefix_variables':len(us)-1,'suffix_variables':len(vs)-1,'certificate':file.name,'terms':len(certificate),'max_polynomial_coefficient':max(polynomial.values()),'strictly_negative_coefficients':sum(v<0 for v in polynomial.values()),'zero_coefficients':sum(v==0 for v in polynomial.values()),'positive_constant_multiplier':bool(positive_constant),'named_multiplier_variables':{'0':list(us[1]),'10':list(us[11])},'oriented_pairs':pairs,'sha256':hashlib.sha256(file.read_bytes()).hexdigest()}
(ROOT/'expanded_polynomial.json').write_text(json.dumps([{'monomial':m,'coefficient':v} for m,v in sorted(polynomial.items())],indent=2)+'\n');print(json.dumps(report,indent=2))
