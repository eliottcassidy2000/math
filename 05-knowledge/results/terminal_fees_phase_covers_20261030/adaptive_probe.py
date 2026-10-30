#!/usr/bin/env python3
"""Discovery/control for a nondeterministic height-adaptive content observer."""
from collections import deque
from fractions import Fraction as F
from pathlib import Path
import json

def need(b,s):
    if not b:raise RuntimeError(str(s))

def run(E):
    q=5**E;D=8*q;N=5*q-5;cap=64
    def h(z):
        k=0
        while k<E and z%5==0:k+=1;z//=5
        return k
    heights=[h(z) for z in range(q)];moduli=[5**max(1,k) for k in heights]
    inv=[pow(2,-a,q) for a in range(cap+1)]
    states=[];offset=[]
    for z,r in enumerate(moduli):
        offset.append(len(states))
        states.extend((z,m) for m in range(1,r) if m%5)
    def index(z,m):return offset[z]+m-m//5-1
    need(len(states)==(16*E+84)*q//25,'state formula')
    def transitions(z,m,a):
        zz=(3*z*inv[a]+1-2*inv[a])%q;rr=moduli[zz]
        small=min(moduli[z],rr);mm=3*m*inv[a]%small
        for m2 in range(mm,rr,small):
            yield index(zz,m2),D*min(h(m2-1),heights[zz])-N*a
    start=index(0,1);V=[None]*len(states);V[start]=0;todo=deque([start]);inside={start}
    while todo:
        s=todo.popleft();inside.remove(s);z,m=states[s]
        for a in range(1,cap+1):
            for t,cost in transitions(z,m,a):
                new=V[s]+cost
                if V[t] is None or new>V[t]:
                    need(new<=0,('positive abstract path',E,s,a,t,new))
                    V[t]=new
                    if t not in inside:todo.append(t);inside.add(t)
    need(all(v is not None for v in V),'finite domain')
    # Separate verification pass with direct destination multiplier formula and layer indicators.
    edges=0
    for s,(z,m) in enumerate(states):
        for a in range(1,cap+1):
            zz=(3*z+pow(2,a,q)-2)*pow(pow(2,a,q),-1,q)%q
            rr=moduli[zz];small=min(moduli[z],rr)
            first=(3*m*pow(2,-a,small))%small
            for mm in range(first,rr,small):
                charge=sum(int(zz%(5**e)==0 and (mm-1)%(5**e)==0) for e in range(1,E+1))
                t=index(zz,mm)
                need(V[s]+D*charge-N*a<=V[t], 'independent edge inequality');edges+=1
    need(max(V)==0 and max(V)+D*E-65*N<=min(V),'all-letter tail')
    return dict(E=E,states=len(states),edges=edges,minimum=min(V),maximum=max(V),D=D,N=N,
                tail=True,status='PASS: all-word abstract enclosure at this precision')

rows=[run(E) for E in range(1,6)]
out=dict(status='VERIFIED finite-precision abstract certificates; uniform E remains OPEN',levels=rows)
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
