#!/usr/bin/env python3
"""Exact finite controls for the proved local-budget theorem; not JC tests."""
import json
from fractions import Fraction
from pathlib import Path

def require(c,msg):
    if not c: raise RuntimeError(msg)

def series(edges,b,steps):
    n=len(b); out=[dict(v) for v in b]; layer=[dict(v) for v in b]
    for _ in range(steps):
        nxt=[{} for i in range(n)]
        for s,v,r,j,c in edges:
            for (a,k),q in layer[s].items():
                key=(r+a,j+k)
                nxt[v][key]=nxt[v].get(key,0)+q*c
        for v in range(n):
            for key,c in nxt[v].items(): out[v][key]=out[v].get(key,0)+c
        layer=nxt
    return out

good=series([(0,1,1,1,1),(1,0,1,0,1)],[{(0,0):1},{}],60)
require(all(2*j-r<=0 for (r,j),a in good[0].items() if a),'good output budget')
require(all(2*j-r<=1 for (r,j),a in good[1].items() if a),'good state offset')
require(len(good[0])==31,'return paths')
bad=series([(0,1,1,2,1),(1,0,1,0,1)],[{(0,0):1},{}],60)
violations=[(r,j,2*j-r) for (r,j),a in bad[0].items() if a and 2*j-r>5]
require(violations[0]==(6,6,6),'first bad row at B5')
signed=series([(0,0,1,1,1),(1,1,1,1,1)],[{(0,0):1},{(0,0):1}],60)
require(signed[0]==signed[1],'exact signed cancellation')
res={'status':'PASS','universe':'specified two-state recurrences through 60 transitions',
     'positive_control_max_weight':max(2*j-r for r,j in good[0]),
     'offset_state_max_weight':max(2*j-r for r,j in good[1]),
     'hostile_first_B5_violation':violations[0],
     'hostile_max_weight':max(2*j-r for r,j in bad[0]),
     'signed_cancellation_checked_terms':len(signed[0]),
     'scope':'controls only; all-row statements are proved in FILTERED_TRANSFER.md'}
print(json.dumps(res,indent=2))
