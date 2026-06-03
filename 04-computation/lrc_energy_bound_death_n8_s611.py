#!/usr/bin/env python3
"""claudebox-S611: the resonance-energy bound proves the median config through n=7, collapses at
n=8 (the structural->computational boundary). See HYP-2165."""
import math, itertools, random
def ghat(m,n):
    d=1.0/n; return -math.sin(2*math.pi*m*d)/(math.pi*m) if m else (1-2*d)
def E_exact(V,M=3):
    n=len(V)+1; tot=0.0
    for m in itertools.product(range(-M,M+1),repeat=len(V)):
        if not any(m): continue
        if sum(a*b for a,b in zip(m,V))!=0: continue
        p=1.0
        for x in m: p*=abs(ghat(x,n))
        tot+=p
    return tot,(1-2/n)**(n-1)
if __name__=="__main__":
    random.seed(7)
    print('n  | main   | median E | frac E<main')
    for n in range(4,10):
        Es=[]; main=None; trials=80 if n<8 else 25
        for _ in range(trials):
            V=sorted(random.sample(range(1,5*n),n-1)); E,main=E_exact(V); Es.append(E)
        Es.sort(); print(f'{n:2d} | {main:.4f} | {Es[len(Es)//2]:.4f}   | {sum(1 for e in Es if e<main)/len(Es):.2f}')
