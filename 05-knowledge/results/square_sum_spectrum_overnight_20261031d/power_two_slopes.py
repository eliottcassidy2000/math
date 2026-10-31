#!/usr/bin/env python3
"""A complete bounded-state audit of the actual power-of-two slope boundary."""
from collections import deque
import json
from pathlib import Path
HERE=Path(__file__).resolve().parent
profiles=[tuple(map(int,s.split())) for s in (HERE.parent/'square_sum_spectrum_overnight_20261031c'/'D_profiles.tsv').read_text().splitlines()]


def children(a):
    F,q,g,h=a
    yield F+q+g,2*q,0,h
    yield F,2*q+h,h,h
    for H,E,p,Q,s,t,B,mask in profiles:
        yield H*F+p*q+(H-s)*g+(B-H)*h,E*q+Q*h,(H-t)*h,H*h


def classified(a):
    F,q,g,h=a
    if not q or q&(q-1):return True
    exponent=q.bit_length()-1;r=0;z=h
    while z>1:
        if z%5:return False
        z//=5;r+=1
    expected=[]
    if exponent>=5*r:
        G=1 if r==0 and exponent==0 else 0
        expected.append((q+h-G,q,G,h))
    if r>=1 and exponent>=5*r-1:
        G=2 if r==1 and exponent==4 else 0
        expected.append((q+h//5-G,q,G,h))
    return a in expected


seen={(1,0,0,1)};pending=deque(seen);power_states=[]
while pending:
    a=pending.popleft()
    if not classified(a):raise RuntimeError(('power slope',a))
    if a[1] and not(a[1]&(a[1]-1)):power_states.append(a)
    for b in children(a):
        if b[1]>256 or b[3]>25:continue
        if b[0]>b[1]+b[3]-b[2]:raise RuntimeError('positive interior invariant')
        if b not in seen:seen.add(b);pending.append(b)
print(json.dumps({'status':'passed','q_bound':256,'h_bound':25,'all_word_states':len(seen),
                  'power_states':sorted(power_states),'scope':'Complete actual-word closure under monotone q,h bounds; the all-parameter statement has a separate predecessor proof.'},indent=2,sort_keys=True))
