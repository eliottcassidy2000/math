#!/usr/bin/env python3
"""Exact optimal mixed-alphabet refresh congestion, retaining each coordinate."""
from itertools import product
from fractions import Fraction as F
from collections import defaultdict
from math import prod
import json
CHECKS=0
def need(x,m):
    global CHECKS
    CHECKS+=1
    if not x:raise RuntimeError(m)
def run(q,opt):
    d=len(q);N=prod(q);a=[F(n-1,n) for n in q];A=sum(a,F(0));c=[1+A-x for x in a]
    w=[x/sum(c) for x in c] if opt else [F(1,d)]*d
    states=list(product(*(range(n) for n in q)));loads=defaultdict(int);counts=defaultdict(int)
    for x,y in product(states,repeat=2):
        z=x;ell=sum(u!=v for u,v in zip(x,y))
        for i in range(d):
            if x[i]!=y[i]:
                v=z[:i]+(y[i],)+z[i+1:];loads[i,z,v]+=ell;counts[i,z,v]+=1;z=v
        need(z==y,'all full coordinates retained at endpoint')
    observed=[]
    for i in range(d):
        ratios={F(v,N*N)/(w[i]/(N*q[i])) for (j,x,y),v in loads.items() if j==i}
        need(ratios=={c[i]/w[i]},'exact coordinate load equals its cut lower bound')
        need(sum(v for (j,x,y),v in loads.items() if j==i)==N*N*a[i]*c[i],'coordinate-cut moment')
        observed.append(c[i]/w[i])
    if opt:need(max(observed)==d+(d-1)*A,'globally optimal coordinate probabilities')
    return {'alphabet':q,'weights':list(map(str,w)),'optimal_weights':opt,'states':N,
      'coordinate_congestion':list(map(str,observed)),'congestion':str(max(observed)),'gap':str(min(w))}
def main():
    alphabets=[(3,),(3,3),(2,3,4)]+[(2,)*r+(3,) for r in range(1,6)]
    rows=[run(q,opt) for q in alphabets for opt in (False,True)]
    return {'status':'FINITE-EXACT PASS; all-alphabet theorem proved separately','cases':rows,'checks':CHECKS}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
