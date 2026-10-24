#!/usr/bin/env python3
"""Independent finite checks for the all-prime positive-defect extraction.
Computes regular-matrix coefficients, not an unrestricted tournament census.
"""
from collections import Counter
from functools import lru_cache
from itertools import combinations
from math import comb
from pathlib import Path
import argparse,json

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def ordinary_coefficient(n,d):
    @lru_cache(None)
    def visit(caps):
        if not sum(caps):return 1
        if sum(caps)%d:return 0
        rows=sum(caps)//d
        if max(caps)>rows:return 0
        ans=0
        for selected in combinations([i for i,c in enumerate(caps) if c],d):
            c=list(caps)
            for i in selected:c[i]-=1
            ans+=visit(tuple(sorted(c)))
        return ans
    return visit((d,)*n)
def symmetric_coefficient(n,d):
    @lru_cache(None)
    def visit(caps):
        if not caps:return 1
        first,*tail=caps;ans=0
        for diagonal in (0,1):
            needed=first-diagonal
            if needed<0:continue
            for selected in combinations([i for i,c in enumerate(tail) if c],needed):
                c=tail[:]
                for i in selected:c[i]-=1
                if any(x>len(c) for x in c):continue
                ans+=visit(tuple(sorted(c)))
        return ans
    return visit((d,)*n)
def subset_key(D,p):return min(tuple(sorted((x+a)%p for x in D)) for a in range(p))
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',default='prime-normal.json');args=ap.parse_args();out=[]
    for p in (3,5,7):
        for d in range(1,(p+1)//2):
            subsets=list(combinations(range(p),d));keys={subset_key(D,p) for D in subsets}
            reflection_fixed=sum(subset_key(tuple(-x%p for x in D),p)==D for D in keys)
            B=comb(p,d);c=comb((p-1)//2,d//2)
            need(len(keys)==B//p,'nonconstant prime subset translation orbits')
            need(reflection_fixed==c,'negative-slope reversible necklace fixed count')
            A=ordinary_coefficient(p,d);S=symmetric_coefficient(p,d)
            N=(A+(p-1)**2*B)//(p*p);F=(S+(p-1)*c)//p
            nonrigid=(p-1)*B//p;nonrigid_sc=B//p+c
            rigid=(A-(p-1)*B)//(p*p);rigid_sc=(S-B-c)//p
            need(N==nonrigid+rigid and F==nonrigid_sc+rigid_sc,'ordinary and fixed error decomposition')
            need((rigid+rigid_sc)%2==0,'merged rigid integrality')
            out.append({'prime':p,'cross_degree':d,'labelled_regular_matrices':A,'symmetric_regular_matrices':S,'subset_necklaces':len(keys),'reversible_subset_necklaces':reflection_fixed,'ordinary_classes':N,'self_converse_classes':F,'ordinary_zero_error':nonrigid,'self_converse_zero_error':nonrigid_sc,'ordinary_positive_error':rigid,'self_converse_positive_error':rigid_sc,'merged_positive_error':(rigid+rigid_sc)//2,'error_per_positive_class':2*p*(p-1)})
    need([r['self_converse_positive_error'] for r in out if r['prime']==7]==[32,984,5016],'prime seven fixed positive sector values')
    need(any(r['reversible_subset_necklaces']<r['subset_necklaces'] for r in out),'not all necklaces reversible hostile')
    result={'status':'PASS','scope':'coefficient recursions and complete subset-necklace checks for primes three five seven; not p7 matrix enumeration','sectors':out}
    Path(args.output).write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps(result,indent=2,sort_keys=True))
if __name__=='__main__':main()
