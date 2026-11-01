#!/usr/bin/env python3
"""Literal-integer controls for the explicit phase exit and exact reset word."""
from itertools import product
from math import comb
import json,random
from primegame_basin import PRIMES,valuations,integer,measures,require

RACK=((17,91),(78,85),(19,51),(23,38),(29,33),(77,29),(95,23),
      (77,19),(1,17),(11,13),(13,11),(15,2),(1,7),(55,1))
def raw(n):
    for i,(u,v) in enumerate(RACK):
        if n%v==0: return n//v*u,'ABCDEFGHIJKLMN'[i]
    raise RuntimeError('total')

counts={'phase_cycle2':0,'phase_cycle6':0,'exit_standard':0,'exit_cycle2':0}
max_phase=max_total=checks=0

def verify(x):
    global max_phase,max_total,checks
    T,Q,k=measures(x)
    if not k: return
    B=(T+1)*comb(Q+2,2)*comb(k+6,6)
    n=integer(x);seen={};elapsed=0
    while True:
        y,s=valuations(n);t,q,h=measures(y)
        require(s==1 and t==T and q<=Q and h<=k,'phase invariant');checks+=1
        require(elapsed<=B,'phase bound')
        if h==0: break
        if n in seen:
            period=elapsed-seen[n]
            require(period in (2,6),'phase period');checks+=1
            counts['phase_cycle'+str(period)]+=1
            max_phase=max(max_phase,elapsed);max_total=max(max_total,elapsed)
            return
        seen[n]=elapsed
        n,letter=raw(n);elapsed+=1
        require(letter in 'ABCDEFGHIJK','positive-phase priority');checks+=1
    require(y[0]==T and y[1]==y[2]==0,'I exit shape');checks+=1
    d=y[3];word='L'*T+'M'*d+'N'+'EF'*T+'K'
    phase_steps=elapsed
    for intended in word:
        n,letter=raw(n);elapsed+=1
        require(letter==intended,'exact reset word');checks+=1
    require(n==5**(T+1)*7**T*13,'exact candidate reset');checks+=1
    require(elapsed<=B+3*T+Q+2,'raw-entry bound');checks+=1
    counts['exit_standard' if T else 'exit_cycle2']+=1
    max_phase=max(max_phase,phase_steps);max_total=max(max_total,elapsed)

cases=0
for x in product(range(2),repeat=10):
    verify(x);cases+=1
for q in PRIMES[4:]:
    for a,b,c,d in product(range(4),repeat=4):
        x=[a,b,c,d,0,0,0,0,0,0];x[PRIMES.index(q)]=1
        verify(tuple(x));cases+=1
rng=random.Random(13509)
for _ in range(200):
    verify(tuple(rng.randrange(7) for _ in range(10)));cases+=1
print(json.dumps(dict(status='FINITE-EXACT PASS',inputs=cases,positive_control_inputs=sum(counts.values()),
                     cases=counts,checks=checks,max_phase_steps=max_phase,max_raw_entry_steps=max_total),
                 sort_keys=True,indent=2))
