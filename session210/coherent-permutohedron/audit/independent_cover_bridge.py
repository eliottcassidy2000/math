#!/usr/bin/env python3
"""Direct permutation check of the gate/permutohedron binomial bridge."""
from independent_carrier_checks import tournament,hp,gate,require
from itertools import permutations
from collections import Counter
from math import comb
import json
from pathlib import Path
ROOT=Path(__file__).parent

def adjacency_histogram(a):
    n=len(a)
    return Counter(sum(a[u]>>v&1 for u,v in zip(p,p[1:])) for p in permutations(range(n)))

def cover_by_allocations(a,q):
    # Independent bin assignment recursion using HP counts of induced bins.
    n=len(a);sub=[hp(a,m) for m in range(1<<n)];sums=[0]*q
    def rec(i):
        if i==n:
            v=1
            for s in sums:v*=sub[s]
            return v
        total=0
        for b in range(q):
            sums[b]|=1<<i;total+=rec(i+1);sums[b]^=1<<i
        return total
    return int(n==0) if q==0 else rec(0)

def binomial_cover(hist,n,q):
    return sum(c*comb(q+f,n) for f,c in hist.items() if q+f>=n)

def main():
    report={'orders':[]}
    for n in range(1,6):
        num=1<<(n*(n-1)//2)
        for bits in range(num):
            a=tournament(n,bits);h=adjacency_histogram(a);g=gate(a);hg=adjacency_histogram(g)
            require(all(h[f]==h[n-1-f] for f in range(n)),('reversal symmetry',n,bits))
            for q in range(0,5):
                p=binomial_cover(h,n,q)
                require(p==cover_by_allocations(a,q),('binomial transform',n,bits,q))
                if q:
                    lhs=binomial_cover(hg,n+2,q)
                    rhs=q*(binomial_cover(h,n,q+1)+q*p-binomial_cover(h,n,q-1))
                    require(lhs==rhs,('gate operator independent permutations',n,bits,q))
        report['orders'].append({'n':n,'labelled_tournaments':num,'q_values':[0,1,2,3,4]})
    report['status']='PASS'
    report['interpretation']='Bin assignments and permutation ascents give the same cover polynomial; gate operator checked without worker subset convolution.'
    (ROOT/'independent_cover_bridge.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
