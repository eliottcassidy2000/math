#!/usr/bin/env python3
"""Independent C3 detector, truncated walk-series classes, interval flips."""
from independent_audit import mat,need
import itertools,json,random
from functools import lru_cache
from pathlib import Path
OUT=Path(__file__).resolve().parent

def cyclic(z,s):
    return all(sum(z[v][w]==1 for w in s)==1 for v in s)

def walk_coefficient(z,target):
    n=len(z)
    @lru_cache(None)
    def suffix(rem,last):
        if not any(rem):return 1
        out=0
        for j in range(n):
            if rem[j] and z[last][j]==1:
                nxt=list(rem);nxt[j]-=1;out+=suffix(tuple(nxt),j)
        return out
    ans=0
    for i in range(n):
        if target[i]:
            rem=list(target);rem[i]-=1;ans+=suffix(tuple(rem),i)
    return ans

def truncated_walks(z):
    n=len(z);out=[]
    @lru_cache(None)
    def ending(counts,last):
        if counts[last]==0:return 0
        counts=list(counts);counts[last]-=1
        if not any(counts):return 1
        return sum(ending(tuple(counts),v) for v in range(n) if counts[v] and z[v][last]==1)
    for c in itertools.product(range(3),repeat=n):out.append(sum(ending(c,j) for j in range(n)))
    return tuple(out)

def main():
    rows=[];rng=random.Random(360042)
    for n in range(1,6):
        groups={};detectors=intervals=0
        triples=list(itertools.combinations(range(n),3))
        for code in range(1<<(n*(n-1)//2)):
            z=mat(n,code);sig=tuple(cyclic(z,s) for s in triples);W=truncated_walks(z)
            need(sig not in groups or groups[sig]==W,('same C3 different truncated W',n,code));groups[sig]=W
            for s,isc in zip(triples,sig):
                c=tuple(n-1 if j in s else 1 for j in range(n));count=walk_coefficient(z,c)
                need(bool(count)==isc,('fullsupport triangle detector',n,code,s,count,isc));detectors+=1
            for mask in range(1,1<<n):
                inside=[j for j in range(n) if mask>>j&1];outside=[j for j in range(n) if not mask>>j&1]
                if len(inside)<2:continue
                if any(len({z[v][w] for v in inside})>1 for w in outside):continue
                flip=[[(-z[v][w] if v in inside and w in inside else z[v][w]) for w in range(n)] for v in range(n)]
                fsig=tuple(cyclic(flip,s) for s in triples)
                need(fsig==sig,('interval inversion C3',n,code,mask));intervals+=1
        row=dict(n=n,tournaments=1<<(n*(n-1)//2),C3_classes=len(groups),triple_detectors=detectors,interval_inversions=intervals)
        rows.append(row);print(row,flush=True)
    extra=[]
    for n,reps in ((6,64),(7,32),(8,16)):
        tests=0
        for _ in range(reps):
            z=mat(n,rng.getrandbits(n*(n-1)//2))
            for s in rng.sample(list(itertools.combinations(range(n),3)),5):
                c=tuple(n-1 if j in s else 1 for j in range(n));count=walk_coefficient(z,c)
                need(bool(count)==cyclic(z,s),('larger detector',n,s,count));tests+=1
        extra.append(dict(n=n,tournaments=reps,detectors=tests))
    report=dict(status='PASS',exhaustive=rows,larger_controls=extra,scope='Every truncated full walk coefficient with each exponent0..2 is constant on each C3 class; exact full-support triangle detector evaluated separately.')
    (OUT/'c3_walk_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)
if __name__=='__main__':main()
