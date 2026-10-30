#!/usr/bin/env python3
"""Independent small table of Fibonacci actions, plus temporal bit decoder."""
import json
from pathlib import Path

def require(ok,msg):
    if not ok: raise RuntimeError(msg)
def step(x,p):return ((x[0]+x[1])%p,x[0])
def orbit_count(p,stride):
    todo={(a,b) for a in range(p) for b in range(p)}-{(0,0)}
    total=0
    while todo:
        x=min(todo); o=set(); y=x
        while y not in o:
            o.add(y)
            for _ in range(stride):y=step(y,p)
        require(y==x,'not an orbit')
        require(not any(((-a)%p,(-b)%p) in o for a,b in o),'negative orbit')
        todo-=o;total+=1
    return total//2
rows=[]
for p,stride,expected in [(3,8,4),(5,4,4),(11,2,12),(19,2,20)]:
    v=(1,0);x=step(v,p);t=1
    while x!=v:x=step(x,p);t+=1
    bits=orbit_count(p,stride)
    require(bits==expected,'bits table')
    rows.append({'prime':p,'pisano_period':t,'preserving_stride':stride,'vertices':p*p,'orientation_bits':bits,'tournaments':2**bits})
w=(0,0,0,1,1)
triples=[tuple(w[(i+j)%5] for j in range(3)) for i in range(5)]
require(len(set(triples))==5,'three serial bits distinguish')
require(2**2<5,'two binary readings cannot suffice')
# g=Q^4 mod5 on each trace-history fibre of Q; phase histories have size5.
Q=(1,1,1,0)
def mul(a,b):return tuple(sum(a[2*i+k]*b[2*k+j] for k in range(2))%5 for i in range(2) for j in range(2))
states=[];x=(1,0,0,1)
for _ in range(20):states.append(x);x=mul(x,Q)
require(x==states[0] and len(set(states))==20,'Q cycle')
fibres={}
for j,x in enumerate(states):
    qx=mul(Q,x); key=((x[0]+x[3])%5,(qx[0]+qx[3])%5)
    fibres.setdefault(key,[]).append(j)
require(len(fibres)==4 and all(len(v)==5 for v in fibres.values()),'fibres')
beacon={}
for js in fibres.values():
    j0=min(js)
    for i in range(5):beacon[(j0+4*i)%20]=w[i]
for js in fibres.values():
    observed={tuple(beacon[(j+4*h)%20] for h in range(3)) for j in js}
    require(len(observed)==5,'actual clock decoder')
summary={'status':'PASS','prime_table':rows,'binary_word':list(w),'three_sample_windows':triples,'actual_Q_mod5_fibres':4,'sample_stride':4,'interpretation':'snapshot alphabet and temporal history are different information budgets'}
print(json.dumps(summary,indent=2))
