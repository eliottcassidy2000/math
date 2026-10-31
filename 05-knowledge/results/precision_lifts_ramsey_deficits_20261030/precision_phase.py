#!/usr/bin/env python3
"""Paid all-one prefixes aligning one actual suffix with retained precision."""
from fractions import Fraction as F
from pathlib import Path
from math import gcd, comb
import json

E=9
q=5**E
w=(2,2,3,1,2,2,1,1,2,1,4)
lam=F(5,8)*(1-F(1,5**11))


def h(x):
    if x==0:return E
    j=0
    while j<E and x%5==0:j+=1;x//=5
    return j


rows=[]
for R in range(1,E+1):
    mod=5**R
    end=pow(3,11,mod)*pow(pow(2,21,mod),-1,mod)%mod
    step=3*pow(2,-1,mod)%mod
    x,k=end,0
    while x!=1:
        x=x*step%mod;k+=1
        if k>=2*5**(R-1):raise ValueError('phase not in the cyclic subgroup')
    m=pow(3*pow(2,-1,q)%q,k,q);z=0
    initial=sum(k//(2*5**(j-1)) for j in range(1,E+1))
    charges=[]
    for a in w:
        precision=max(R,h(z))
        m=m*3*pow(2,-a,q)%q
        z=(3*z+2**a-2)*pow(2,-a,q)%q
        charges.append(min(h(z),h((m-1)%q)))
    if z or h((m-1)%q)<precision:raise ValueError('incompatible terminal lift')
    total=initial+sum(charges[:-1])+E
    rows.append(dict(R=R,ones_prefix=k,clock=k+sum(w),initial_charge=initial,
        suffix_actual_charges=charges,abstract_total=total,
        endpoint_price=str(lam*(k+sum(w))-total),actual_terminal_content=charges[-1]))

# Independent literal integer recurrence for the two minimal-precision hostiles.
for R, k in ((1,0),(2,4)):
    P,B,Q,S=1,0,1,0
    for a in (1,)*k+w:
        oldz=B-P+Q
        P,B,Q=3*P,3*B+Q,Q*2**a
        c=h(gcd(P-Q,B));S+=c
    abstract=S-c+E
    if rows[R-1]['abstract_total']!=abstract or h(P-Q)<max(R,h(oldz)):
        raise ValueError('literal hostile mismatch')
    suffix = '' if R==1 else '_r2'
    audit=json.loads(Path(f'adaptive9{suffix}_audit.json').read_text())
    if F(-audit['terminal'],audit['D'])!=lam*(k+21)-abstract:
        raise ValueError('hostile fails to attain certified terminal value')

print(json.dumps(dict(status='FINITE-EXACT witnesses; optimality only for R1 and R2 with audited potentials',
    suffix=w,rows=rows,R2_fixed_length_universe=comb(24,14),
    R2_unique_last_lift_witness={'ones_prefix':4,'actual_total':5,'abstract_total':12},pass_checks=True),indent=2))
