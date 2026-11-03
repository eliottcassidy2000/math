#!/usr/bin/env python3
"""Independent literal path/cut controls for UNIFORM_PLANCHEREL.md."""
from itertools import combinations, permutations
from math import comb, factorial
from fractions import Fraction
import json


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def templates(n):
    """Enumerate unordered path covers once, independently of subset DP."""
    out=[]
    for order in permutations(range(n)):
        for cuts in range(1 << (n-1)):
            blocks=[];start=0
            for i in range(n-1):
                if cuts >> i & 1:
                    blocks.append(order[start:i+1]);start=i+1
            blocks.append(order[start:])
            if blocks != sorted(blocks):
                continue
            arcs=0
            for block in blocks:
                for u,v in zip(block,block[1:]):
                    arcs |= 1 << (u*n+v)
            out.append((len(blocks),arcs))
    return out


def graph(n,bits,tournament=True):
    mask=0
    if tournament:
        for j,(u,v) in enumerate(combinations(range(n),2)):
            if bits >> j & 1:u,v=v,u
            mask |= 1 << (u*n+v)
    else:
        for j,(u,v) in enumerate(( (u,v) for u in range(n) for v in range(n) if u!=v )):
            if bits >> j & 1:mask |= 1 << (u*n+v)
    return mask


def strong(n,mask):
    for root in range(n):
        seen={root}; todo=[root]
        while todo:
            u=todo.pop()
            for v in range(n):
                if v not in seen and (mask >> (u*n+v)) & 1:
                    seen.add(v);todo.append(v)
        if len(seen)!=n:return False
    return True


def profile(n,mask,ts):
    p=[0]*(n+1)
    for k,arcs in ts:
        if arcs & mask == arcs:p[k]+=1
    return p


def fall(r,k):
    out=1
    for i in range(k):out*=r-i
    return out


def F(p,r):
    return sum(c*fall(r,k) for k,c in enumerate(p))


def check(p):
    n=len(p)-1
    delta=[0]+[k*(k+1)*(p[k+1] if k<n else 0)-(n-k)*p[k] for k in range(1,n+1)]
    require(all(d>=0 for d in delta),('negative cut defect',p,delta))
    for r in range(9):
        require(r*F(p,r+1)-(n+r)*F(p,r)==F(delta,r), ('polynomial defect',p,r))
    for r in range(1,9):
        require(F(p,r)>=p[1]*comb(n+r-1,r-1),('Hamilton bound',p,r))
        for s in range(1,r+1):
            require(F(p,s)*comb(n+r-1,r-1)<=F(p,r)*comb(n+s-1,s-1),('ratio bound',p,s,r))
        if r>=3:
            require(F(p,r-2)*(n+r-1)*(n+r-2)<=F(p,r)*(r-1)*(r-2),('two-slot ratio',p,r))
    if p[1]:
        tv=Fraction(p[1],2*F(p,3))
        require(tv<=Fraction(1,(n+1)*(n+2)),('S3 total variation',p))
    return delta


def main():
    report={'tournament_cases':0,'general_digraph_cases':0,'bidirected_sharpness_cases':0,'random_tournament_controls':[]}
    for n in range(1,6):
        ts=templates(n);sums=[0]*(n+1);max_ratio=Fraction(0)
        count=1<<(n*(n-1)//2)
        strong_H=0; strong_F3=0; strong_max_ratio=Fraction(0)
        for bits in range(count):
            mask=graph(n,bits);p=profile(n,mask,ts);check(p)
            if strong(n,mask):
                strong_H+=p[1];strong_F3+=F(p,3)
                strong_max_ratio=max(strong_max_ratio,Fraction(p[1],F(p,3)))
            report['tournament_cases']+=1
            sums=[x+y for x,y in zip(sums,p)]
            max_ratio=max(max_ratio,Fraction(p[1],F(p,3)))
        for k in range(1,n+1):
            expected=Fraction(factorial(n)*comb(n-1,k-1),factorial(k)*2**(n-k))
            require(Fraction(sums[k],count)==expected,('expected profile',n,k))
        require(F(sums,3)==(2*n*n+1)*sums[1],('expected F3/H',n))
        require(max_ratio>=Fraction(1,2*n*n+1),('quadratic sharpness existence',n))
        if n>=4:
            alpha=1-Fraction(n-1,2**(n-2))
            require(strong_H>=alpha*sums[1],('strong weighted union',n))
            require(strong_max_ratio>=alpha/Fraction(2*n*n+1),('strong quadratic existence',n))
        report['random_tournament_controls'].append({'n':n,'max_H_over_F3':str(max_ratio),'strong_max_H_over_F3':str(strong_max_ratio),'averaged_H_over_F3':str(Fraction(sums[1],F(sums,3)))})
        complete=graph(n,(1<<(n*(n-1)))-1,tournament=False)
        p=profile(n,complete,ts)
        require(not any(check(p)),('bidirected equality',n))
        report['bidirected_sharpness_cases']+=1
    n=3;ts=templates(n)
    for bits in range(1<<(n*(n-1))):
        check(profile(n,graph(n,bits,tournament=False),ts))
        report['general_digraph_cases']+=1
    report['status']='PASS'
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
