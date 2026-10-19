#!/usr/bin/env python3
"""Finite hypotheses for two infinite good-prime parameter families."""
import json
rows=[]
for p,c,expected in [(3,2,2),(5,1,3)]:
 cycles=set()
 for x in range(p):
  seen=[]
  while x not in seen:
   seen.append(x);x=(x*x+c)%p
  cyc=seen[seen.index(x):]
  cycles.add(min(tuple(cyc[k:]+cyc[:k]) for k in range(len(cyc))))
 if len(cycles)!=1:raise ValueError('not a unique cycle')
 cyc=next(iter(cycles));mult=1
 for x in cyc:mult=mult*2*x%p
 if len(cyc)!=expected or mult!=0:raise ValueError('critical-cycle certificate failed')
 rows.append({'prime':p,'parameter_residue':c,'unique_cycle':list(cyc),'multiplier_mod_p':mult,'consequence':'All affine Q_p-periodic points lie on the unique lifted cycle; rational periodic points, if any, have its exact period.'})
print(json.dumps(rows,indent=2))
