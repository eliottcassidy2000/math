#!/usr/bin/env python3
"""Finite-width digital reflection and recovery of the existing Mersenne fan."""
import json
from pathlib import Path
import sys
from source_inverse import f
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'square_sum_spectrum_overnight_20261031c'))
from audit_port_words import adjacency,count


def reverse(x,k):return sum(((x>>j)&1)<<(k-1-j) for j in range(k))


def mersenne_address(k):
    if k<4:raise ValueError('fan family begins at k=4')
    return (1<<(k+1))-(3*(1<<(k//2)) if k%2==0 else 5*(1<<((k-1)//2)))-2


checks=0
for k in range(1,17):
    for x in range(1<<k):
        rhs=f(x)+reverse(x,k)-x+x%2-(x>>(k-1))
        if f((1<<k)-1-x)!=rhs:raise RuntimeError(('reflection',k,x))
        checks+=1
fan_checks=0
for k in range(4,129):
    address=mersenne_address(k)
    if f(address)!=(1<<k)-1:raise RuntimeError(('fan address',k))
    word=[('S',int(b)) for b in bin(address)[2:]]
    rows=adjacency(word)
    if len(rows)!=k+3 or rows[0]!=(1<<len(rows))-2:raise RuntimeError('one explicit source padding vertex')
    losing={j+1 for j in range(k+1) if not(rows[1]>>(j+2)&1)}
    expected={k//2,k+1-k//2,k+1}
    if losing!=expected:raise RuntimeError('not the prior THM-4748 fan')
    # Removing the named dominating source leaves the exact earlier fan.
    reduced=[row>>1 for row in rows[1:]]
    if count(reduced,[0])!=(1<<k)-1:raise RuntimeError('actual fan count')
    fan_checks+=1
for k in range(129,4097):
    if f(mersenne_address(k))!=(1<<k)-1:raise RuntimeError('large exact identity')
print(json.dumps({'status':'passed','reflection_checks':checks,'actual_prior_fan_checks':fan_checks,
                  'fan_identity_through_exponent':4096,'scope':'The Mersenne existence family is inherited THM-4748; this is its exact digital-address and source-padding identification.'},indent=2,sort_keys=True))
