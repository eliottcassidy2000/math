#!/usr/bin/env python3
"""Exact Hamming-chart balance, weighted Hall hostile and optimal route congestion."""
from fractions import Fraction as F
from math import comb
from collections import defaultdict
import json

CHECKS=0
def need(x,m):
    global CHECKS
    CHECKS+=1
    if not x:raise RuntimeError(m)
def main():
    rows=[]
    for r in range(1,9):
        N=1<<r;loads=defaultdict(int);distance2=0
        for k in range(r):
            left=[x for x in range(N) if x.bit_count()==k]
            right=[y for y in range(N) if y.bit_count()==k+1]
            mass=defaultdict(F)
            for x in left:
                for j in range(r):
                    y=x^(1<<j)
                    if y in right:mass[y]+=F(1,len(left)*(r-k))
            need(set(mass)==set(right),'complete adjacent-chart coverage')
            need(all(mass[y]==F(1,len(right)) for y in right),'exact target conditional marginal')
            need(comb(r,k)*(r-k)==comb(r,k+1)*(k+1),'reverse detailed balance')
        for x in range(N):
            for y in range(N):
                d=(x^y).bit_count();distance2+=d*d;z=x
                for j in range(r):
                    if (x^y)>>j&1:
                        w=z^(1<<j);loads[z,w]+=d;z=w
                need(z==y,'ordered coordinate route reaches actual endpoint')
        need(F(distance2,N*N)==F(r*(r+1),4),'exact second distance moment')
        congestion=max(F(v*2*r,N) for v in loads.values())
        need(congestion==F(r*(r+1),2),'achieves the universal lower bound')
        rows.append({'r':r,'states':N,'congestion':str(congestion),'E_distance_squared':str(F(distance2,N*N))})
    hostiles=[]
    for kappa in (1,2,3,10):
        supply=F(kappa,kappa+2);capacity=F(2,kappa+2)
        need((supply>capacity)==(kappa>2),'exact weighted Hall boundary')
        hostiles.append({'kappa':kappa,'source_singleton_mass':str(supply),'neighbor_pair_mass':str(capacity),
                         'direct_balanced_coupling_impossible':supply>capacity})
    return {'status':'FINITE-EXACT PASS; all-r proof separate','hypercubes':rows,'weighted_Hall_hostile':hostiles,'checks':CHECKS}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
