#!/usr/bin/env python3
import json
import importlib.util
from pathlib import Path
spec=importlib.util.spec_from_file_location("inherited",Path(__file__).with_name("inherited.py"))
I=importlib.util.module_from_spec(spec)
spec.loader.exec_module(I)

def opposite(a):
    full=(1<<len(a))-1
    return [full ^ (1<<i) ^ r for i,r in enumerate(a)]

def join(a,b):
    # Both inputs rooted at 0; second is recorded in converse convention.
    b=opposite(b)
    n=len(a); mask=((1<<len(b))-1)<<n
    rows=[r|mask for r in a]+[r<<n for r in b]
    rows[0]^=1<<n
    rows[n]|=1
    return rows

def triple(a): return tuple(I.profile(a,0)[:3])
def value(x,y):
    h,e,u=x; k,f,v=y
    return h*k-e*f+u*v

def inverse_from_bank(N,bank):
    found=set(); cap=(2*N-1)//3
    for h in range(1,cap+1,2):
      for k in range(1,cap+1,2):
       for e in range(h+1):
        for f in range(k+1):
         if e==h and f==k: continue
         t=N-h*k+e*f
         if t<1: continue
         for u in range(h,t//k+1):
          if t%u: continue
          v=t//u
          x=(h,e,u); y=(k,f,v)
          if x in bank and y in bank: found.add((x,y))
    return found

def main():
    graphs=[a for n in range(1,6) for a in I.tournaments(n)]
    triples={tuple(a):triple(a) for a in graphs}
    bank={triples[tuple(a)]:a for a in graphs}
    literal=proper=0
    for a in graphs:
      for b in graphs:
        x=triples[tuple(a)]; y=triples[tuple(b)]
        h,e,u=x; k,f,v=y; N=value(x,y)
        if e<h or f<k:
          I.require(max(h,k)<N,'strict descent')
          I.require(3*h<=2*N-1 and 3*k<=2*N-1,'uniform quantitative descent')
          I.require(u*v<=N-1,'response product')
          if e<h: I.require(2*N>=(3*h+1)*k,'left port bound')
          if f<k: I.require(2*N>=(3*k+1)*h,'right port bound')
          proper+=1
        else: I.require(N==h*k,'degenerate product')
        if len(a)<=4 and len(b)<=4:
          I.require(I.profile(join(a,b),0)[0]==N,'literal join law')
          literal+=1
    controls=[]
    for N in range(3,42,2):
      direct={(x,y) for x in bank for y in bank if (x[1]<x[0] or y[1]<y[0]) and value(x,y)==N}
      inv=inverse_from_bank(N,bank)
      I.require(inv==direct,'bounded physical inverse')
      controls.append([N,len(inv)])
    # Same H, different proper scalar images.
    tt2=[2,0];tt3=[6,4,0];single=[0]
    I.require(triple(tt2)==(1,0,2) and triple(tt3)==(1,0,4),'sidecar controls')
    I.require(value(triple(tt2),triple(single))==3 and value(triple(tt3),triple(single))==5,'scalar replacement hostile')
    return dict(status='FINITE-EXACT',graphs=len(graphs),physical_triples=len(bank),pair_universe=len(graphs)**2,proper_pairs=proper,literal_pairs=literal,inverse_bounded_bank=controls,all_target_coverage='OPEN')
if __name__=='__main__': print(json.dumps(main(),indent=2,sort_keys=True))
