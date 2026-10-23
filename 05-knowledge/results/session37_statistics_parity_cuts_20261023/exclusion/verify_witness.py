#!/usr/bin/env python3
"""Independent permutation enumeration of det/per/H and cycle packings."""
from itertools import permutations, combinations
from collections import Counter
from pathlib import Path
import json,math

def adjacency(n,mask):
    a=[[0]*n for _ in range(n)]
    for b,(i,j) in enumerate(combinations(range(n),2)):
        a[i][j]=(mask>>b)&1;a[j][i]=1-a[i][j]
    return a

def verify(n,mask):
    a=adjacency(n,mask)
    pc=Counter();dc=Counter();h=0
    for p in permutations(range(n)):
        h+=all(a[p[i]][p[i+1]] for i in range(n-1))
        if all(i==p[i] or a[i][p[i]] for i in range(n)):
            fixed=sum(i==p[i] for i in range(n))
            inv=sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))
            pc[fixed]+=1;dc[fixed]+=(-1)**inv
    cyc={};tot=Counter()
    for k in range(3,n+1):
        for s in combinations(range(n),k):
            r=s[0];bits=sum(1<<i for i in s)
            count=sum(all(a[p[i]][p[i+1]] for i in range(k-2)) and a[r][p[0]] and a[p[-1]][r] for p in permutations(s[1:]))
            cyc[bits]=count;tot[k]+=count
    # Enumerate all vertex-disjoint cycle sets by smallest occupied vertex,
    # using optional unoccupied vertices as a deterministic branch.
    def packing(available):
        if not available:return Counter({():1})
        low=available&-available
        out=packing(available^low).copy()
        sub=available
        while sub:
            if sub&low and cyc.get(sub,0):
                for profile,ct in packing(available^sub).items():
                    out[tuple(sorted(profile+(sub.bit_count(),)))]+=cyc[sub]*ct
            sub=(sub-1)&available
        return out
    from functools import lru_cache
    packing=lru_cache(None)(packing)
    pack=packing((1<<n)-1)
    ocf=sum((2**len(prof))*ct for prof,ct in pack.items() if all(k%2 for k in prof))
    assert h==ocf,(h,ocf)
    p2=[0]*(n+1);d2=[0]*(n+1)
    for prof,ct in pack.items():
        k=n-sum(prof)
        p2[k]+=ct;d2[k]+=(-1)**sum(x-1 for x in prof)*ct
    assert p2==[pc[i] for i in range(n+1)]
    assert d2==[dc[i] for i in range(n+1)]
    return dict(n=n,mask=mask,outneighbor_masks=[sum(a[i][j]<<j for j in range(n)) for i in range(n)],hamiltonian=h,per_xIplusA=p2,det_xIplusA=d2,cycle_totals=dict(tot),packings={','.join(map(str,k)):v for k,v in sorted(pack.items())})

if __name__=='__main__':
    res=[verify(9,m) for m in [29759988540,18875161993]]
    assert res[0]['per_xIplusA']==res[1]['per_xIplusA']
    assert res[0]['det_xIplusA']==res[1]['det_xIplusA']
    assert res[0]['hamiltonian']!=res[1]['hamiltonian']
    Path(__file__).with_name('witness_verification.json').write_text(json.dumps(res,indent=2)+'\n')
    for r in res:print(json.dumps(r,sort_keys=True))
