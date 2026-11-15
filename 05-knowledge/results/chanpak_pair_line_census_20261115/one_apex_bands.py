#!/usr/bin/env python3
"""Exact ranked-core bands, all one-apex strong presentations, no padding."""
from reproduce import stats
from itertools import combinations
def pair(m,mask,weights):
    h,q=1,1<<m
    zeros=mask^((1<<(m-1))-1)
    for d,a,b in weights:
        x=((mask<<d)&zeros).bit_count()
        h += x*a
        q += x*b
    return h,q
def main():
    controls=0
    print('n,words,distinct_H,distinct_pairs,distinct_HplusQ,min_HplusQ,max_HplusQ')
    for m in range(3,22):
        W=[(d,1<<max(1,d-1),(1<<(m-d-1))*(3**(d-1)+1)) for d in range(1,m-1)]
        pairs=set()
        for inner in range(1<<(m-3)):
            mask=1|(inner<<1)
            h,q=pair(m,mask,W)
            pairs.add((h,q))
            if m<=8:
                code=sum(1<<j for j,(u,v) in enumerate(combinations(range(m),2))
                         if v==m-1 and mask>>u&1)
                assert stats(m,code)[:2]==(h,q)
                controls+=1
        values={h+q for h,q in pairs}
        L=m-2
        w={d:a+b for d,a,b in W}
        low=1+(1<<m)+w[L]+sum(min(w[d],w[L-d]) for d in range(1,L))
        a,b=(m-1)//2,m//2
        high=2*3**(m-1)+7*2**(m-1)-2*(3**a*2**b+2**a*3**b)-3*(2**a+2**b)+5
        assert (min(values),max(values))==(low,high)
        assert high<2*3**(m-1)
        print(m,1<<(m-3),len({h for h,q in pairs}),len(pairs),len(values),min(values),max(values),sep=',')
    print('direct subset/cycle controls', controls)
if __name__=='__main__': main()
