#!/usr/bin/env python3
"""Focused finite controls; all general statements are proved in CANONICAL_RETURNS.md."""
from itertools import product
from collections import deque
from check import C,D,one_attach,endpoint_state,strong,require

V=(1,2,0)

def run(word,v=V):
    for c in word:v=(C if c=='C' else D)(v)
    return v


def canonical(v):return v[2]==0 and v[1]==2*v[0]

def physical(word):
    x=(2,4,1);q=0
    for c in word:
        x,a=one_attach(x,q)
        if c=='C':q=a
    return x,q


def main():
    tested=0
    for n in range(9):
        for letters in product('CD',repeat=n):
            w=''.join(letters);v=run(w)
            predicted=(not w) or w=='CD' or w.endswith('DCD')
            require(canonical(v)==predicted,'return language: '+w)
            require(v[0]>0 and v[1]>0 and v[1]>=v[0]+v[2],'positive slack')
            tested+=1
    for n,expected in enumerate([{3},{5},{9},{15,17}]):
        hs={sum(run(''.join(w))[:2]) for w in product('CD',repeat=n)}
        require(hs==expected,'small depth boundary')
    for m in range(1,11):
        for w,f in [('DCC'+'D'*m+'CD',9+10*2**m),('DCCC'+'D'*m+'CD',15+19*2**m)]:
            require(run(w)==tuple(f*a for a in V),'parametric multiplier')
    for w,f in [('DCCDCD',29),('CCDDCD',29),('DCCDDCD',49),('DCCCDDCD',91)]:
        v=run(w);require(v==tuple(f*a for a in V),'specific multiplier')
        x,q=physical(w)
        require(strong(x),'literal source strength')
        require(endpoint_state(x,q)==v,'literal source endpoint state')
        print(w,'multiplier',f,'H',sum(v[:2]),'vertices',len(x))
    for k in range(1,6):
        m=7**k
        x=1;units=set()
        while x not in units:units.add(x);x=3*x%m
        require(x==1 and len(units)==6*7**(k-1),'unit generator')
        seen={3%m};queue=deque(seen)
        while queue:
            a=queue.popleft()
            for f in (3,91):
                b=a*f%m
                if b not in seen:seen.add(b);queue.append(b)
        require(len(seen)==m,'complete 7-power residues')
        print('modulus',m,'all',len(seen),'residues reached')
    print('PASS:',tested,'word controls; parametric returns; four literal sources; residues through7^5.')

if __name__=='__main__':main()
