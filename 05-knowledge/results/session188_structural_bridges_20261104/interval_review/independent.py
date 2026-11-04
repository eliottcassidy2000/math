#!/usr/bin/env python3
from itertools import combinations,permutations,product
from functools import lru_cache
from pathlib import Path
import json,sys,math
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'simultaneous_review'))
from independent import rows,literal_full,need


def converse(a):return tuple(sum(1<<v for v in range(len(a)) if v!=u and not a[u]>>v&1) for u in range(len(a)))


def substitute(q,blocks):
    ids=[];start=0
    for block in blocks:ids.append(list(range(start,start+len(block))));start+=len(block)
    output=[0]*start
    for b,block in enumerate(blocks):
        for i,u in enumerate(ids[b]):
            for j,v in enumerate(ids[b]):
                if block[i]>>j&1:output[u]|=1<<v
            for c,names in enumerate(ids):
                if q[b]>>c&1:
                    for v in names:output[u]|=1<<v
    return tuple(output)


@lru_cache(None)
def H(a):
    return sum(all(a[u]>>v&1 for u,v in zip(w,w[1:])) for w in permutations(range(len(a))))


@lru_cache(None)
def profile(a):
    n=len(a);full=(1<<n)-1;h=[0]*(full+1)
    for S in range(1,full+1):
        names=[u for u in range(n) if S>>u&1]
        h[S]=sum(all(a[u]>>v&1 for u,v in zip(w,w[1:])) for w in permutations(names))
    dp=[[0]*(full+1) for _ in range(n+1)];dp[0][0]=1
    for k in range(1,n+1):
        for S in range(1,full+1):
            B=S
            while B:
                dp[k][S]+=h[B]*dp[k-1][S^B];B=(B-1)&S
    return tuple(dp[k][full] for k in range(1,n+1))


def main():
    inversion=[]
    for k in range(1,7):
        pairs=list(combinations(range(k),2));preserved=0
        for code in range(1<<len(pairs)):
            delta={p:(code>>j)&1 for j,p in enumerate(pairs)}
            eta={p:delta[p]^delta.get((p[0]-1,p[1]),0)^delta.get((p[0],p[1]+1),0)^delta.get((p[0]-1,p[1]+1),0) for p in pairs}
            for i,j in pairs:
                need(sum(v for (a,b),v in eta.items() if a<=i and j<=b)%2==delta[i,j],'interval inversion')
            q=rows(k,code)
            transitive=sorted(x.bit_count() for x in q)==list(range(k))
            preserved+=transitive
        need(preserved==math.factorial(k),'k factorial')
        inversion.append(dict(k=k,patterns=1<<len(pairs),transitive_patterns=preserved))
    fixtures=[((0,),(0,),(0,)),((2,0),(0,),(0,)),((2,0),(2,0),(0,)),
              ((2,0),(2,0),(2,0)),((2,4,1),(0,),(2,0))]
    controls=0
    ordered=(6,4,0)
    for blocks in fixtures:
        source=substitute(ordered,blocks);base=profile(source)
        for qcode in range(8):
            quotient=tuple(rows(3,qcode));qh=H(quotient)
            for reverse_bits in range(8):
                changed=[converse(a) if reverse_bits>>j&1 else a for j,a in enumerate(blocks)]
                target=substitute(quotient,changed);got=profile(target)
                need(all(a>=b for a,b in zip(got,base)),'profile domination')
                need((got==base)==(qh==1),'equality iff transitive quotient')
                need(got[0]>=qh*math.prod(H(tuple(a)) for a in blocks),'Hamilton multiplier')
                controls+=1
    p=Path(__file__).resolve().parent.parent/'interval_modules'/'normal.json'
    receipt=json.loads(p.read_text());examples=[]
    for name,rec in receipt['full_contrast_records'].items():
        a,b=rec['original'],rec['changed'];n=len(a)
        pa,fa=literal_full(a);pb,fb=literal_full(b)
        difference={s:pa.get(s,0)-pb.get(s,0) for s in set(pa)|set(pb) if pa.get(s,0)!=pb.get(s,0)}
        need(rec['complete_contrast']==[[s&((1<<n)-1),s>>n,v] for s,v in sorted(difference.items())],'complete named contrast')
        need(rec['scalar_difference']==H(tuple(a))-H(tuple(b)),'scalar example')
        need(rec['full_degree']==max(s.bit_count() for s in difference),'quartic degree')
        for mask in range(len(fa)):
            need(sum(v for s,v in difference.items() if s&mask==s)==fa[mask]-fb[mask],'all exact named queries')
        examples.append(dict(name=name,degree=rec['full_degree'],all_queries=len(fa),nonzero_coefficients=len(difference)))
    rec=receipt['whole_converse_control'];a,b=rec['original'],rec['converse'];n=len(a)
    pa,fa=literal_full(a);pb,fb=literal_full(b)
    difference={s:pa.get(s,0)-pb.get(s,0) for s in set(pa)|set(pb) if pa.get(s,0)!=pb.get(s,0)}
    need(rec['complete_contrast']==[[s&31,s>>5,v] for s,v in sorted(difference.items())],'whole converse full contrast')
    need(all(s.bit_count()<=3 for s in difference),'whole converse quartics')
    return dict(status='PASS',interval_inversion=inversion,arbitrary_block_profile_controls=controls,
                exact_named_examples=examples,whole_converse_coefficients=len(difference))


if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
