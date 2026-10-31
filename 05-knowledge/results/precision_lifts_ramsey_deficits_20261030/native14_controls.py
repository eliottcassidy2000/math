#!/usr/bin/env python3
"""Exact arithmetic and port checks for the necessary five-form grammar."""
from fractions import Fraction as F
from math import ceil
import json

eta=F(183105426,48828125)
d=F(3,5**13)
q=5**14
def need(x):
    if not x:raise ValueError('native14 arithmetic or port failure')
def compose(first,second):
    return second[0]*first[0]%q,(second[0]*first[1]+second[1])%q
def inverse(a):
    u=pow(a[0],-1,q)
    return u,-u*a[1]%q
def ones(e):
    a=pow(3*pow(2,-1,q)%q,2*5**(e-1),q)
    return a,(a-1)%q

need(2*(eta-2)-1>1)
need(2*(eta-1)-2>1)
patterns=[]
for target in (F(1),F(4,5)):
    permitted=[]
    for a11 in range(4):
        for a12 in range(4):
            lower=eta-3+F(31,125)*a11+F(6,25)*a12
            if lower>=target:continue
            cap=(target-lower)/d
            max_clock=ceil(cap)-1
            need(lower+d*max_clock<target<=lower+d*(max_clock+1))
            permitted.append((a11,a12))
            if a11+a12==0:
                placements=['bare']
            else:
                placements=['prefix','suffix']
            end=(1,0) if target==1 else (1+5**13,5**13)
            for place in placements:
                block=ones(11 if a11 else 12) if place!='bare' else (1,0)
                core=(compose(inverse(block),end) if place=='prefix' else
                      compose(end,inverse(block)) if place=='suffix' else end)
                actual=(compose(block,core) if place=='prefix' else
                        compose(core,block) if place=='suffix' else core)
                need(actual==end)
                need(core[0]%5**11==1 and core[1]%5**11==0)
                patterns.append(dict(target=str(target),a11=a11,a12=a12,placement=place,
                    strict_core_clock_cap=str(cap),max_integer_core_clock=max_clock,
                    required_core_map_mod_5_14=core))
    need(permitted==([(0,0),(0,1),(1,0)] if target==1 else [(0,0)]))

print(json.dumps(dict(status='FINITE-EXACT necessary grammar arithmetic; no large native-core enumeration',
    eta10=str(eta),slope_increment=str(d),modulus=q,patterns=patterns,
    all_word_claim='proved necessary form, not existence or exclusion of the cores',
    pass_checks=True),indent=2))
