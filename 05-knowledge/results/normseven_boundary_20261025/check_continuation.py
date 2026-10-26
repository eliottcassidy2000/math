#!/usr/bin/env python3
"""Finite exact continuation-depth tests, without importing producer code."""
from itertools import combinations, permutations, product
from pathlib import Path
import hashlib
import json


GATES=0


def need(ok, why):
    global GATES
    GATES+=1
    if not ok:
        raise RuntimeError(why)


def span(vectors):
    out={0}
    for v in vectors:
        out|={z^v for z in tuple(out)}
    return frozenset(out)


def carries(parts,n):
    e={}
    for chosen in product(*parts):
        z=0
        sums=[0]*n
        for v in chosen:
            z^=v
            for i in range(n):
                sums[i]+=(v>>i)&1
        if z in e:
            return None
        e[z]=tuple(s//2 for s in sums)
    return e if len(e)==1<<n else None


def requirements(e,z,words,n):
    out=[0]*n
    for y in words:
        for i in range(n):
            out[i]=max(out[i],e[z^y][i]-((z&y)>>i&1))
    return tuple(out)


def main():
    census={}
    one_step_failures={}
    n=4
    subspaces=sorted({span(v) for d in (1,2,3) for v in combinations(range(1,16),d)
                      if len(span(v))==1<<d},key=lambda s:(len(s),tuple(sorted(s))))
    for r in (2,3,4):
        count=0
        failures=0
        for parts in combinations(subspaces,r):
            size=1
            for c in parts:
                size*=len(c)
            if size!=16:
                continue
            e=carries(parts,n)
            if e is None:
                continue
            count+=1
            for z in range(16):
                full=requirements(e,z,range(16),n)
                bounded=requirements(e,z,[y for y in range(16) if y.bit_count()<=3],n)
                local=requirements(e,z,[0,1,2,4,8],n)
                need(full==bounded,'n-1 binary bank failed')
                if r==2:
                    need(full==local,'two-packet single-letter criterion failed')
                failures+=int(full!=local)
        census[r]=count
        one_step_failures[r]=failures
    need(census[2]==400,'all unordered complementary binary dimension-four pairs')
    need(one_step_failures[3]>0,'three-packet hostile must be present')

    def code(labels):
        return frozenset(sum(((a&v).bit_count()%2)<<i for i,v in enumerate(labels)) for a in range(8))
    codes=sorted({code((1,)+p) for p in permutations(range(2,8))},key=lambda s:tuple(sorted(s)))
    pairs=[pair for pair in combinations(codes,2) if pair[0]&pair[1]=={0}]
    need(len(pairs)==120,'complete Fano pair universe')
    for parts in pairs:
        for j in range(7):
            indices=[i for i in range(7) if i!=j]
            punctured=[frozenset(sum(((v>>i)&1)<<k for k,i in enumerate(indices)) for v in c) for c in parts]
            e=carries(punctured,6)
            need(e is not None,'punctured complement')
            for z in range(64):
                local=requirements(e,z,[0,1,2,4,8,16,32],6)
                need(local==tuple(1-(z>>i&1) for i in range(6)),'six-letter norm-seven criterion')

    hostiles=[]
    for n in range(2,13):
        parts=[frozenset((0,1))]+[frozenset((0,1|(1<<j))) for j in range(1,n)]
        e=carries(parts,n)
        k=n-2
        first_bad=n+1
        for y in range(1<<n):
            exponent=k+(y&1)
            z=(y&~1)|(exponent%2)
            coefficient=(exponent//2,)+(0,)*(n-1)
            belongs=all(a>=b for a,b in zip(coefficient,e[z]))
            if y.bit_count()<=k:
                need(belongs,'source failed before claimed distinguishing depth')
            if not belongs:
                first_bad=min(first_bad,y.bit_count())
            # Positive control: the conductor generator passes every reduced word.
            exponent=n-1+(y&1)
            z=(y&~1)|(exponent%2)
            coefficient=(exponent//2,)+(0,)*(n-1)
            need(all(a>=b for a,b in zip(coefficient,e[z])),'conductor positive control')
        need(first_bad==n-1,'sharp delayed failure')
        y=(1<<n)-2
        z=y|(k%2)
        need(e[z][0]-k//2==1,'precisely one missing b_0')
        hostiles.append({'coordinates':n,'source_root0_exponent':k,'first_failure_depth':first_bad,
                         'distinguishing_root_labels':list(range(1,n))})
    print(json.dumps({'status':'FINITE-EXACT','gates':GATES,
                     'dimension_four_unordered_decompositions':census,
                     'characters_with_one_step_strictly_weaker':one_step_failures,
                     'fano_pair_charts':840,'sharp_star_hostiles':hostiles,
                     'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2,sort_keys=True))


if __name__=='__main__':
    main()
