#!/usr/bin/env python3
"""Exact actual-word and independent first-arrival checks, no floating point."""
from fractions import Fraction as F
from itertools import groupby
from math import gcd
from pathlib import Path
import json

def need(b,msg):
    if not b: raise RuntimeError(msg)

def v5(n):
    need(n!=0,'finite valuation required');n=abs(n);e=0
    while n%5==0:n//=5;e+=1
    return e

def affine(w):
    p=q=1;b=0
    for a in w:p,q,b=3*p,2**a*q,3*b+q
    return p,q,b

def address(p,q,h):
    """j with p/q=(3/2)^j mod5^h, lifted five choices per digit."""
    if h==0:return []
    need(p*pow(q,-1,5)%5 in (1,4),'active quadratic coset')
    j=0 if p*pow(q,-1,5)%5==1 else 1
    ans=[(-j)%2];period=2
    for e in range(2,h+1):
        mod=5**e;target=p*pow(q,-1,mod)%mod;base=3*pow(2,-1,mod)%mod
        lifts=[j+k*period for k in range(5) if pow(base,j+k*period,mod)==target]
        need(len(lifts)==1,'unique lifted address');j=lifts[0];period*=5
        ans.append((-j)%period)
    return ans

def analyze(w):
    p=q=1;b=A=charge=0;epochs=[];initial=0
    for a in w:
        p,q,b=3*p,2**a*q,3*b+q;A+=a
        c=v5(gcd(abs(p-q),b));charge+=c
        if a>1:
            h=v5(b+q-p);active=p*pow(q,-1,5)%5 in (1,4)
            epochs.append(dict(a=a,h=h,k=v5(2**a-2),active=active,L=0,
                               offsets=address(p,q,h) if active else [],direct=c))
        elif epochs:epochs[-1]['L']+=1;epochs[-1]['direct']+=c
        else:initial+=1
    old=F(0);new=F(0);exact=0
    for x in epochs:
        a,h,k,active,L=(x[z] for z in ('a','h','k','active','L'))
        old-=F(5*a,8);new-=F(5*a,8)
        if active:
            old+=max(h,k)-F(5*L,8*5**h)
            startup=F(0);actual=0
            for e,r in enumerate(x['offsets'],1):
                t=2*5**(e-1)
                if L>=r:actual+=1+(L-r)//t;startup+=1-F(r,t)
            need(actual==x['direct'],'first-arrival identity vs literal gcd')
            new+=startup-F(5*L,8*5**h)
            x['startup']=str(startup)
        else:
            old-=F(5*L,8);new-=F(5*L,8);actual=0
            need(x['direct']==0,'inactive no charge')
        exact+=actual
    initial_charge=sum(initial//(2*5**(e-1)) for e in range(1,initial.bit_length()+1))
    need(exact+initial_charge==charge,'partition all source prefixes')
    need(charge<=F(5*A,8)+new<=F(5*A,8)+old,'bounds retain direction')
    return dict(length=len(w),clock=A,charge=charge,old_deficit=str(old),
                arrival_deficit=str(new),epochs=epochs)

def main():
    U=[3,1,3]+[1]*6+[3,1,3];w=U+[1]*231+U
    result=analyze(w)
    need(result['old_deficit']=='769/5000','positive witness')
    need(result['charge']==154 and result['clock']==271,'literal witness')
    bounds=[]
    # C((1), VU)=B_(VU)+Q_(VU)-P_(VU). Retain chronological VU.
    for k in range(len(U)+1):
        p,q,b=affine(U[:k]+U);bounds.append(v5(b+q-p))
    need(max(bounds)==5,'controller bounds every suffix phase')
    # Independent finite controls: every positive composition of clock<=15.
    count=0
    def compositions(A,prefix=()):
        if A==0:yield prefix
        else:
            for a in range(1,A+1):yield from compositions(A-a,prefix+(a,))
    rejected=0
    for A in range(1,16):
        for ww in compositions(A):
            r=analyze(ww);count+=1;rejected+=F(r['arrival_deficit'])>0
    family=[]
    for L in range(1250):
        r=analyze(U+[1]*L+U)
        need(1250*r['charge']<=781*r['clock'],'known sharp level-five bound')
        if L in (0,1,231,1249):family.append(dict(L=L,charge=r['charge'],
                old_deficit=r['old_deficit'],arrival_deficit=r['arrival_deficit']))
    out=dict(status='REFUTED old universal budget; PROVED arrival identity and bounded-controller reduction',
             witness=result,runs=[(a,len(list(g))) for a,g in groupby(w)],
             controller_phase_bounds=bounds,composition_controls=count,
             arrival_rejections_in_controls=rejected,family_controls=family)
    path=Path(__file__).with_name('collatz_audit.json');path.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('witness','runs')},indent=2))

if __name__=='__main__':main()
