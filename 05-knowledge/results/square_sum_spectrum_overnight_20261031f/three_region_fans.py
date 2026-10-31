#!/usr/bin/env python3
"""Exact dyadic coefficients for fixed left/center/right one-apex decorations."""
import json,sys
from fractions import Fraction as F
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'square_sum_spectrum_overnight_20261031e'))
from centered_fans import debt,count

def two(k):return F(1<<k) if k>=0 else F(1,1<<(-k))
def coefficients(L,S,R,D):
    U=sum((two(-l) for l in L),F(0));V=sum((two(-r) for r in R),F(0))
    A=sum(two(s) for s in S);B=sum(two(-s) for s in S)
    lead=two(D-1)*V*(1-U)
    middle=(1-U)*A/2-two(D-1)*V*B
    constant=1+sum(two(l-1) for l in L)-debt(L)-debt(S)-debt(R)
    return lead,middle,constant,U,V

def rows(L,S,R,D,a):
    K=2*a+D;losing=set(L)|{a+s for s in S}|{K-r for r in R}
    if a+min(S)-max(L,default=0)<2 or a+D-max(R)-max(S)<2:
        raise ValueError('separated cross gaps')
    out=[0]*(K+1)
    for i in range(1,K+1):
        for j in range(i+1,K+1):out[i]|=1<<j
        if i in losing:out[i]|=1
        else:out[0]|=1<<i
    return out

def subsets(items,mandatory=()):
    for mask in range(1<<len(items)):
        yield sorted(set(mandatory)|{x for j,x in enumerate(items) if mask>>j&1})

def main():
    controls=rigidity=0;examples=[]
    for L in subsets([2,3,4]):
        for R in subsets([1,2,3],[0]):
            U=sum((two(-l) for l in L),F(0));V=sum((two(-r) for r in R),F(0));z=V*(1-U)
            if not F(1,2)<z<2:raise RuntimeError('leading interval')
            if (z==1)!=(not L and R==[0]):raise RuntimeError('dyadic unit rigidity')
            rigidity+=1
            for S in subsets([1,2],[0]):
                for D in range(5):
                    coeff=coefficients(L,S,R,D)
                    a=max(2,max(L,default=0)+2,max(R)+max(S)-D+2)
                    for offset in (0,1):
                        n=a+offset;actual=count(rows(L,S,R,D,n),[0]);expected=coeff[0]*two(2*n)+coeff[1]*two(n)+coeff[2]
                        if expected!=actual:raise RuntimeError('three-region adjacency formula')
                        controls+=1
                    if L==[2] and R in [[0],[0,2]] and S==[0,1] and D==1:
                        examples.append({'left':L,'middle':S,'right_offsets':R,'D':D,'quadratic_coefficients':[str(x) for x in coeff[:3]],'relative_leading_factor':str(z)})
    approx=[]
    for r in range(1,9):
        R=list(range(0,2*r+1,2));V=sum(two(-s) for s in R);factor=F(3,4)*V
        if factor!=1-two(-2*r-2):raise RuntimeError('finite reciprocal train error')
        approx.append([r,str(factor),str(1-factor)])
    print(json.dumps({'status':'passed','native_adjacency_controls':controls,'finite_leading_rigidity_controls':rigidity,
                      'examples':examples,'approximate_reciprocal_train':approx,
                      'scope':'Fixed separated three-region grammar; growing boundary trains leave this grammar and require fresh overlap control.'},indent=2,sort_keys=True))

if __name__=='__main__':main()
