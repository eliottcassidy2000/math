"""Actual boundary: the ordinal profile dichotomy fails for strong Q4.

This is a mixed-profile witness, NOT a zero-rigidity counterexample.
"""
import sys,json,hashlib
from pathcover_probe import profile,induce
from reference_paths import tournament,c3_blowup,extend
q=tournament(4,4);T=c3_blowup(q);B=induce(T,list(range(1,12)));f=profile(T);rows=[]
for S in [867,995,996]:
 U=extend(B,S);ff=profile(U);d=[b-a for a,b in zip(f,ff)]
 if not min(d)<0<max(d):raise RuntimeError('mixed profile')
 rows.append({'cut':S,'parent_masks':B,'child_masks':U,'original_H':f[0],'new_H':ff[0],'difference':d,'C3_singleton_context_gap':3*d[0]+d[1],'C3_triangle_context_gap':sum(a*b for a,b in zip(d,[99,261,216,36]))})
for row in rows:print(json.dumps(row,sort_keys=True))
print('PASS sha256='+hashlib.sha256(json.dumps(rows,separators=(',',':')).encode()).hexdigest())
