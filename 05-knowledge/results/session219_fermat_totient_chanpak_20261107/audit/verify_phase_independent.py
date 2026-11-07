#!/usr/bin/env python3
"""Independent zero-phase certificate by solving for y, not a pair census."""
import json
from pathlib import Path
P,Q=641,6700417
zeros={}
for p in [P,Q]:
    powers={pow(2,j,p):j for j in range(64)}
    assert len(powers)==64
    rows=[]
    for x,r in powers.items():
        if x==1:continue
        y=(1-2*pow(x-1,-1,p))%p
        if y in powers:rows.append((r,powers[y]))
    zeros[p]=sorted(rows)
assert zeros[Q]==[(1,32),(32,1)]
assert set(zeros[Q])<=set(zeros[P])
# Linear matrices act on column state D,E,U,d.
matrices={'O':[[1,0,0,0],[0,0,0,0],[0,0,2,0],[0,0,0,1]],'Z':[[1,-1,1,0],[0,0,0,1],[0,0,2,-1],[0,0,0,1]],'J':[[1,-1,1,0],[1,-1,0,0],[1,0,1,0],[1,0,0,0]]}
def apply(M,v):return [sum(a*b for a,b in zip(row,v)) for row in M]
v=[1,0,2,1]
for c in 'OOZJOZJJOJJOZ':v=apply(matrices[c],v)
assert v==[4909,859,6013,859]
v20=v[:]
for _ in range(20):v20=apply(matrices['O'],v20)
h=v20[0]-v20[1]+v20[2]
assert h==941*Q and h%P==380
# Verify generating identity algebra by matching coefficients.
# (x+1)(2x^2-4x+3)=2x^3-2x^2-x+3 = 2+(x-1)(2x^2-1).
receipt={'status':'PASS','method':'solve y=1−2/(x−1) in each prime field, then inverse lookup in order64 subgroup','zero_phases':zeros,'escape_state':v,'initial_h':h,'initial_mod641':h%P,'scope':'physical realization imports Session200 O/Z/J model; no independent literal graph enumeration'}
Path(__file__).with_name('phase_independent_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
