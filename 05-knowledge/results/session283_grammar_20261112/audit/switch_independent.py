#!/usr/bin/env python3
"""Direct matrix controls of the contraction-deck identities."""
from independent import graph,hp,need
from itertools import product
from math import comb
import json

def substitute(q,lengths,cyclic):
    verts=[(i,j) for i,l in enumerate(lengths) for j in range(l)]
    def arc(i,x,j,y):
        if i!=j:return q[i][j]
        if i in cyclic:return (y-x)%3==1
        return x<y
    return [[False if u==v else arc(*u,*v) for v in verts] for u in verts]

def analyze(q,sites,lengths):
    k=len(sites);deck={};out={}
    for mask in range(1<<k):
        selected={sites[j] for j in range(k) if mask>>j&1}
        ls=[1 if i in selected else l for i,l in enumerate(lengths)]
        deck[mask]=hp(substitute(q,ls,set()))
        out[mask]=hp(substitute(q,lengths,selected))
    for e,h in out.items():
        summands=[2**s.bit_count()*deck[s] for s in deck if s&e==s]
        need(h==sum(summands),'exact expansion')
        for r in range(1,5):
            short=sum(2**s.bit_count()*deck[s] for s in deck if s&e==s and s.bit_count()<r)
            need((h-short)%2**r==0,'exact low bits')
        bs=sum((deck[1<<j]-1)//2 for j in range(k) if e>>j&1)
        need((h-deck[0]-2*e.bit_count()**2-4*bs)%8==0,'mod8 formula')
    if k>=2:need(len(set(out.values()))<(max(out.values())-min(out.values()))//2+1,'cube gap')
    return deck,out

def main():
    families=outputs=0
    for n in range(1,4):
        for code in range(1<<comb(n,2)):
            q=graph(n,code)
            for mask in range(1<<n):
                sites=[i for i in range(n) if mask>>i&1]
                deck,out=analyze(q,sites,[3 if i in sites else 1 for i in range(n)])
                families+=1;outputs+=len(out)
    rows=[14,4,8,2];q=[[bool(rows[i]>>j&1) for j in range(4)] for i in range(4)]
    deck,out=analyze(q,[0,1],[3]*4)
    need(list(deck.values())==[2721,2721,51,51],'12vertex deck')
    need(list(out.values())==[2721,8163,2823,8469],'12vertex outputs')
    print(json.dumps({'status':'PASS','families':families,'outputs':outputs,'mixed_deck':deck,'mixed_outputs':out},sort_keys=True,indent=2))
if __name__=='__main__':main()
