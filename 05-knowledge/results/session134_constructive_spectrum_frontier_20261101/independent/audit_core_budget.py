#!/usr/bin/env python3
"""Independent exact integer tests of the SCC-count capped core budget.
Uses direct component-size dynamic programming versus the (r,u,v) formula.
The floor values are an inherited cited input, not reproved by this audit.
"""
from math import comb
import json

def need(c,msg):
    if not c:raise RuntimeError(msg)
def floor_h(n):
    if n%3==0:return 3*5**(n//3-1)
    if n%3==1:return 5**(n//3)
    return 9*5**(n//3-1)
def formula(X,b):
    ans=0;arg=None
    for r in range(1,b+1):
        for u in range(2*r+1):
            z=3**u;v=0
            while z<=X:
                if v>=max(0,r-u) and r+2*u+3*v>ans:ans=r+2*u+3*v;arg=(r,u,v)
                z*=5;v+=1
    return ans,arg
def simplified(X,b):
    if not b:return 0
    best=0
    for u in range(2*b+1):
        a=3**u
        if a>X:break
        v=0
        while a*5<=X:a*=5;v+=1
        best=max(best,2*u+3*v+min(b,u+v))
    return best

def direct_table(Xmax,bmax):
    floors=[];n=3
    while floor_h(n)<=Xmax:floors.append((n,floor_h(n)));n+=1
    states={(0,0):1}
    for r in range(1,bmax+1):
        for (rr,order),h in list(states.items()):
            if rr!=r-1:continue
            for n,f in floors:
                val=h*f
                if val>Xmax:continue
                key=(r,order+n)
                if val<states.get(key,Xmax+1):states[key]=val
    return states

def bound(X,b):
    N,_=formula(X,b);M=N*(N-1)//2
    return sum(comb(M,j) for j in range(min(b,M)+1)),N

def main():
    states=direct_table(2000,8);cases=0
    for X in range(1,2001):
        for b in range(9):
            by_sizes=max((n for (r,n),h in states.items() if r<=b and h<=X),default=0)
            by_exponents,arg=formula(X,b)
            need(by_sizes==by_exponents==simplified(X,b),f'core budget mismatch {X,b,by_sizes,by_exponents,arg}')
            cases+=1
    rows=[];first=None
    for k in range(2,65):
        B,N=bound(2**k-1,6);interval=2**(k-2)
        if B<interval and first is None:first=k
        if k>=60:rows.append({'k':k,'J6':N,'B6':B,'odd_interval_size':interval,'strict_shortage':B<interval})
    need(first==64,'first sixarc dyadic counting gate')
    need(rows[-1]['J6']==88 and rows[-1]['B6']==4359919199804281268,'gate numeric witness')
    print(json.dumps({'status':'PASS','direct_component_knapsack_comparisons':cases,'formula_paths':['direct component-size product knapsack','r,u,v enumeration','optimized u-only maximum'],'X_range':[1,2000],'SCC_budget_range':[0,8],'six_arc_gate':rows,'first_certified_dyadic_shortage_k':first,'scope':'J_b exact for at most b nontrivial SCCs; only upper bound for feedback-arc budget b. Dyadic shortage is a counting-certificate failure, not location of the first missing value.'},indent=2,sort_keys=True))
if __name__=='__main__':main()
