#!/usr/bin/env python3
"""Exact controls for attachments through a prefix of one shuffle chain."""
from fractions import Fraction
from itertools import product, combinations
from math import comb
from pathlib import Path
import json
from gold_terminal_chain import count,require,deficits
from gold_two_chain_fork import core, direct_flags, margin, upsets

def C(n,k):return comb(n,k) if 0<=k<=n else 0

def bound(m,n,j):
    r=m-j;den=(r+n)*(r+n-1);Q=n*n-(r+1)*n-r*(r-1)
    return min(Fraction(n-m,m+n),Fraction(r*(r-1),den),Fraction(-Q,den))

def nested(D,j):
    U=list(upsets(D))
    for W in product(U,repeat=j):
        if all(W[k+1]&W[k]==W[k+1] for k in range(j-1)):yield W

def build(m,n,H,D,W,p):
    base,_,e,f=core(m,n,H);base=list(base[:-2]);N=len(base)
    for v,x in enumerate(D):
        mask=x<<N
        for i,S in enumerate(W):
            if S>>v&1:mask|=1<<i
        base.append(mask)
    return tuple(base),(0,N+p),e,f

def interleave_weight(N,D,W,p,positions,af,pl):
    pred=[(1<<i)-1 for i in range(N)]
    for v,x in enumerate(D):
        mask=x<<N
        for t,S in zip(positions,W):
            if S>>v&1:mask|=1<<(t-1)
        if af:mask|=1
        pred.append(mask)
    if pl:pred[N+p]|=((1<<len(pred))-1)^(1<<(N+p))
    return count(tuple(pred))

def chain_table(m,n,j,L,h=0):
    T=m+n;N=T+h;r=m-j;A=B=Z=0
    for t in range(j,n+j+1):
        w=C(t-1,j-1)*C(N+L-t,L)
        A+=w*C(T-t-1,r-1);B+=w*C(T-t-2,r-1);Z+=w*C(T-t,r)
    return [[A,0],[B,Z-A-B]]

def main():
    actual=0;by_case=[]
    shapes=[(0,),(0,1),(0,0),(0,0,1),(0,1,3),(0,0,0)]
    for m,n,j in ((4,5,1),(5,5,2),(7,8,2),(8,8,3),(10,11,3)):
        local=0
        for H,D in product(((),(0,0)),shapes):
            for W in nested(D,j):
                for p in range(len(D)):
                    if any(x>>p&1 for x in D):continue
                    pred,ports,e,f=build(m,n,H,D,W,p);F=direct_flags(pred,ports,e,f)
                    require(margin(F)>=bound(m,n,j),('prefix margin',m,n,j,D,W,p))
                    for k,t in enumerate(F):
                        q=Fraction(t[0][0],sum(sum(t,[])))
                        require(Fraction(m-j,m-j+n)<=q<=Fraction(m,m+n),'prefix balance interval')
                        if k in (1,3):require(q<=Fraction(m-1,m+n-1),'a1-first upper bound')
                    actual+=1;local+=1
        by_case.append(dict(m=m,n=n,j=j,checks=local,universal_margin=str(bound(m,n,j))))
    fibre_checks=0;fibres=[]
    for m,n,j,D,W,p in ((5,5,2,(0,0,1),(7,4),2),(8,8,3,(0,0,1),(7,6,4),2),(8,8,3,(0,0,3),(7,5,4),2)):
        N=m+n+1
        for af,pl in product((False,True),repeat=2):
            vals=[]
            for t in range(j,n+j+1):
                prefixes=[(1,)+x for x in combinations(range(2,t),j-2)] if af else list(combinations(range(1,t),j-1))
                weights=[interleave_weight(N,D,W,p,x+(t,),af,pl) for x in prefixes]
                vals.append(Fraction(sum(weights),len(weights)));fibre_checks+=len(weights)
            require(all(x>=y for x,y in zip(vals,vals[1:])),'averaged prefix fibre decreasing')
            fibres.append(dict(m=m,n=n,j=j,a1_first=af,p_last=pl,weights=[str(x) for x in vals]))
    chain_checks=0;obstructions=[]
    for m,n,j in ((4,4,2),(5,6,2),(7,7,3),(5,5,2),(7,8,2)):
        for L in (1,2,5,11,17,30,93):
            D=tuple((1<<i)-1 for i in range(L));W=((1<<L)-1,)*j
            pred,ports,e,f=build(m,n,(),D,W,L-1)
            got=direct_flags(pred,ports,e,f)[0];expected=chain_table(m,n,j,L)
            require(got==expected,'chain collar concentration formula');chain_checks+=1
        for L in range(1,101):
            t=chain_table(m,n,j,L)
            if min(sum(deficits(t),[]))<0:
                obstructions.append(dict(m=m,n=n,j=j,first_failure_side_chain_length=L,table=t,deficits=deficits(t)));break
    feasible=0
    for m in range(2,101):
        for j in range(1,m):
            r=m-j
            require((bound(m,m,j)>=0)==(j*(m-1)<=r*r),'exact feasible prefix length')
            feasible+=1
    fib=[0,1]
    for _ in range(25):fib.append(fib[-1]+fib[-2])
    squares=[]
    for k in range(10):
        d,b,e=fib[2*k+1:2*k+4];m=d*e;j=d*d
        require((m-j)**2==j*(m-1),'Fibonacci-square prefix boundary')
        require(bound(m,m,j)==0,'Fibonacci boundary margin')
        squares.append(dict(k=k,m=m,n=m,j=j))
    out=dict(status='prefix attachment theorem proved; exact controls passed',actual_nested_upset_cases=actual,by_parameter=by_case,individual_fibre_counts=fibre_checks,averaged_fibres=fibres,actual_chain_concentration_checks=chain_checks,finite_strategy_obstructions=obstructions,prefix_feasibility_checks=feasible,fibonacci_square_boundary=squares)
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
