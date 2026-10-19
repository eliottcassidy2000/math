#!/usr/bin/env python3
"""Exact positive and hostile controls for the simple-resonance lemma."""
import json
rows=[]
for p,r in [(3,2),(7,3),(31,5)]:
 def it(x,k,m):
  for _ in range(k):x=(x*x+2*x)%m
  return x
 q=p*p;mod=p**3
 base=[x for x in range(p,q,p)]
 if not all(it(x,r,q)==x for x in base):raise ValueError('positive control failed')
 lifts=[x+k*q for x in base for k in range(p)]
 if not all(it(x,r,mod)!=x for x in lifts):raise ValueError('nonlift control failed')
 if ((1+p)*(-p)+(-p)**2)!=-p:raise ValueError('unit-hypothesis hostile failed')
 rows.append({'p':p,'r':r,'map':'t^2+2t','fixed_point':0,'multiplier':2,'noncentral_mod_p2_points':len(base),'lifts_mod_p3_without_period_dividing_r':len(lifts),'lambda_minus_one_hypothesis_hostile':'(1+p)t+t^2 has the second fixed point -p','status':'PASS'})
print(json.dumps(rows,indent=2))
