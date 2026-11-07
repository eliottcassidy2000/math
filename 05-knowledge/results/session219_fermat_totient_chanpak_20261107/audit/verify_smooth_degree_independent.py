#!/usr/bin/env python3
import json
from pathlib import Path
N=20
legal=[(1,0),(0,2),(2,1),(0,3)]
relaxed=[(a,b) for a in range(N+1) for b in range(N+1) if (a,b) not in [(0,0),(0,1),(1,1)]]
def solve(gens):
    dp={(0,0):0}
    for a in range(N+1):
        for b in range(N+1):
            if a==b==0:continue
            vs=[dp[a-u,b-v]+1 for u,v in gens if u<=a and v<=b and (a-u,b-v) in dp]
            if vs:dp[a,b]=max(vs)
    return dp
lo,hi=solve(legal),solve(relaxed)
assert lo==hi
for (a,b),v in lo.items():
    expect=a if b==0 else a-1 if b==1 else a+b//2
    assert v==expect,(a,b,v,expect)
assert len(lo)==(N+1)**2-2
receipt={'status':'PASS','a_b_range':[0,N],'attained_cases':len(lo),'excluded_cases':[[0,1],[1,1]],'lower_dictionary':legal,'upper_relaxation':'all exponent pairs other than unit,7,21; may include nonexistent strong values','scope':'finite dynamic programming independent of all-parameter written proof'}
Path(__file__).with_name('smooth_degree_independent_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
