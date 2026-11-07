#!/usr/bin/env python3
"""Finite probe only; deterministic random universes are not exhaustive claims."""
from random import Random
from time import monotonic
import json
from verify import H,tiles,require

def dpH(n,ones):
    incoming=[sum(1<<u for u in range(v+1,n)) for v in range(n)]
    for x,y in ones:
        incoming[y-1]^=1<<(x-1)
        incoming[x-1]^=1<<(y-1)
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n):dp[1<<v][v]=1
    for mask in range(1,1<<n):
        bs=mask
        while bs:
            vbit=bs&-bs;bs-=vbit;v=vbit.bit_length()-1
            prev=mask^vbit
            if not prev:continue
            incoming_prev=incoming[v]&prev
            total=0
            while incoming_prev:
                b=incoming_prev&-incoming_prev;incoming_prev-=b
                total+=dp[prev][b.bit_length()-1]
            dp[mask][v]=total
    return sum(dp[-1])

def main():
    start=monotonic();rng=Random(2222026);out={'seed':2222026,'universes':[]}
    for n in range(3,6):
        cells=tiles(n)
        for mask in range(1<<len(cells)):
            ones={p for i,p in enumerate(cells) if mask>>i&1}
            require(dpH(n,ones)==H(n,{p:1 for p in ones}))
    for n,limit in [(7,32768),(8,3000)]:
        cells=tiles(n);universe=range(1<<len(cells)) if n==7 else [0]+[rng.randrange(1,1<<len(cells)) for _ in range(limit)]
        row={'n':n,'sampling':'exhaustive' if n==7 else 'random with replacement plus zero','tested':0,'positive':0,'equal':0,'negative':0,'min_positive_gap':None}
        for mask in universe:
            ones={p for i,p in enumerate(cells) if mask>>i&1}
            a=dpH(n,ones);b=dpH(n+2,{(x+2,y) for x,y in ones});gap=b-a
            row['tested']+=1;row['positive' if gap>0 else 'equal' if gap==0 else 'negative']+=1
            if gap>0:row['min_positive_gap']=min(gap,row['min_positive_gap'] or gap)
            if gap<=0 and mask:
                row['hostile']={'mask':mask,'ones':sorted(ones),'old_H':a,'new_H':b};break
        out['universes'].append(row)
    out['elapsed_seconds']=round(monotonic()-start,3)
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
