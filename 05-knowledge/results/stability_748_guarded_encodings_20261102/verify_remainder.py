#!/usr/bin/env python3
"""Exact scalar controls for the unspent-remainder lemma, not a new order theorem."""
from fractions import Fraction as F
from itertools import product
import json,random

def need(ok,message):
    if not ok:raise RuntimeError(message)

def main():
    t=F(1,39);d=8*t;lam=t/8;checks=0
    def check(light,heavy):
        L=len(light);U=sum(w for w,S in light);Dl=sum(S for w,S in light)
        Dh=sum(S for w,S in heavy);D=Dl+Dh;Mh=sum(w for w,S in heavy);m=max(1,len(heavy))
        beta=max(F(0),Mh-F(m-1,4)-t)
        residue=U-t*(U//t)
        J=(Dl+U-2*t*L)**2/(4*L)+residue*(t-residue)/4+beta*Dh
        W=sum(max(w*S,w*S/2+S*S/4) for w,S in light+heavy)
        need(W>=t*D-lam*(d*L-4*U)+t*U/4+J,'retained aggregate invoice')
        xmin=max(F(0),min(D,2*t*L-U+2*beta*L))
        f=lambda x:(x+U-2*t*L)**2/(4*L)+beta*(D-x)
        need(f(Dl)>=f(xmin),'clamped minimum')
        need(0<=residue<t,'fractional last weight')
    for L in range(1,4):
        for weights in product([F(0),t/2,t],repeat=L):
            for masses in product([F(0),t,3*t],repeat=L):
                for heavy in [[],[(2*t,t)],[(F(1,4),2*t),(F(1,4),t)]]:
                    check(list(zip(weights,masses)),heavy);checks+=1
    rng=random.Random(314)
    for _ in range(1000):
        L=rng.randrange(1,20);m=rng.randrange(1,10)
        light=[(t*F(rng.randrange(101),100),F(rng.randrange(101),100)) for i in range(L)]
        heavy=[(t+(F(1,4)-t)*F(rng.randrange(1,101),100),F(rng.randrange(101),100)) for i in range(m)]
        check(light,heavy);checks+=1
    print(json.dumps({'status':'passed','scalar_group_controls':checks,'seed':314,
        'scope':'Checks the proved aggregate remainder lemma. Does not establish the proposed tail payment or certify any order below314.'},indent=2,sort_keys=True))

if __name__=='__main__':main()
