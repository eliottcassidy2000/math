#!/usr/bin/env python3
"""Separate finite controls: literal affine content, closed epochs and online F2."""
from itertools import product
from fractions import Fraction as F
from pathlib import Path
import json

def need(b,s):
    if not b: raise RuntimeError(str(s))

def height(z, E):
    h=0
    while h<E and z%5==0:
        h+=1;z//=5
    return h

def lam(e): return F(5,8)*(1-F(1,5**e))

def closed_epoch(m,h,L,r):
    if h==0 or m%5 not in (1,4): return F(0)
    total=F(0);last=None
    for e in range(1,min(r,h)+1):
        q=5**e;u=m%q;d=0
        while u!=1:
            u=u*3*pow(2,-1,q)%q;d+=1
            need(d<2*5**(e-1), 'low-phase order')
        total+=max(0,1+(L-d)//(2*5**(e-1)));last=d
    if h>r and L>=last:
        total+=(h-r)+(lam(h)-lam(r))*(L-last)
    return total

def check(word,E):
    q=5**E;P=Q=1;B=0;S=0;states=[];m=z=0
    initial=0
    for a in word:
        if a!=1: break
        initial+=1
    exact_initial=sum(initial//(2*5**j) for j in range(E))
    online=F(exact_initial);pending=False
    for i,a in enumerate(word):
        P,B,Q=3*P,3*B+Q,Q*2**a
        c=min(height(P-Q,E),height(B,E));S+=c
        m=P*pow(Q,-1,q)%q;z=(B-P+Q)*pow(Q,-1,q)%q;h=height(z,E)
        states.append((m,h))
        if i<initial: continue
        low=int(h>=1 and m%5==1)+int(h>=2 and m%25==1)
        reward=F(low)
        if a>1:
            pending=h>2 and m%5 in (1,4) and m%25!=1
            if h>2 and m%25==1: reward+=h-2
        elif pending:
            if m%25==1: reward+=h-2;pending=False
        elif h>2 and m%5 in (1,4): reward+=lam(h)-lam(2)
        online+=reward
    f1=f2=F(exact_initial)
    starts=[i for i,a in enumerate(word) if a>1]
    for j,i in enumerate(starts):
        end=starts[j+1] if j+1<len(starts) else len(word)
        mm,hh=states[i];L=end-i-1
        f1+=closed_epoch(mm,hh,L,1);f2+=closed_epoch(mm,hh,L,2)
    need(S<=f2<=f1,('envelope order',word,E,S,f2,f1))
    need(f2==online,('online/closed disagreement',word,E,f2,online))
    need(exact_initial<=lam(E)*initial,'initial all-ones payment')
    if m%25==1: need(not pending,'terminal flag')
    return dict(length=len(word),clock=sum(word),content=S,F1=str(f1),F2=str(f2),
                endpoint_multiplier_mod25=m%25,endpoint_defect_height=height(z,E))

def main():
    n=0
    for length in range(1,8):
        for word in product(range(1,5),repeat=length):
            check(word,6);n+=1
    R=(5,)+(1,)*8+(5,);W=R+(1,)*233+R;U=W+(1,)*8+W
    hostile=check(U,9)
    extra=[]
    for k in (0,1,9,10,49,50,271,781250):
        # The long case checks the exact initial payment formula without allocating its word.
        if k==781250:
            need(sum(k//(2*5**j) for j in range(9))==lam(9)*k,'sharp initial run')
        else:
            extra.append(check((1,)*k+(5,1,2,1,5),9))
    out=dict(status='FINITE-EXACT semantics controls',alphabet='1..4',lengths='1..7',
             words=n,precision=6,hostile_U_E9=hostile,initial_phase_controls=extra,
             sharp_initial_clock=781250)
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
