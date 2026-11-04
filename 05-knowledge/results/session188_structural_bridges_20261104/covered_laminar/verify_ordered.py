#!/usr/bin/env python3
"""Independent exact controls for arbitrary-block order detection. No third-party imports."""
from itertools import combinations, permutations, product
from functools import lru_cache
from pathlib import Path
import json


def need(ok, why):
    if not ok: raise RuntimeError(why)


def tournament(n, code):
    a=[0]*n
    for b,(u,v) in enumerate(combinations(range(n),2)):
        if code>>b&1:a[u]|=1<<v
        else:a[v]|=1<<u
    return tuple(a)


def converse(a):
    return tuple(((1<<len(a))-1)^row^(1<<i) for i,row in enumerate(a))


def induced(a, names):
    return tuple(sum(1<<j for j,v in enumerate(names) if a[u]>>v&1) for u in names)


@lru_cache(None)
def hp_subsets(a):
    n=len(a);d=[[0]*n for _ in range(1<<n)];h=[0]*(1<<n);h[0]=1
    for v in range(n):d[1<<v][v]=1
    for s in range(1,1<<n):
        for v in range(n):
            z=d[s][v]
            if z:
                rem=a[v]&~s
                while rem:
                    b=rem&-rem;rem-=b;d[s|b][b.bit_length()-1]+=z
        h[s]=sum(d[s])
    return tuple(h)


def hp(a):return hp_subsets(tuple(a))[-1]


@lru_cache(None)
def profile(a):
    n=len(a);H=hp_subsets(a);full=(1<<n)-1;d=[0]*(1<<n);d[0]=1;out=[]
    for k in range(1,n+1):
        nxt=[0]*(1<<n)
        for s in range(1,full+1):
            t=s
            while t:
                nxt[s]+=H[t]*d[s^t];t=(t-1)&s
        out.append(nxt[full]);d=nxt
    return tuple(out)


def substitute(blocks, order):
    offsets=[];n=0
    for a in blocks:offsets.append(n);n+=len(a)
    out=[0]*n
    for a,o in zip(blocks,offsets):
        for i,row in enumerate(a):out[o+i]|=row<<o
    for i,j in combinations(order,2):
        mask=((1<<len(blocks[j]))-1)<<offsets[j]
        for u in range(offsets[i],offsets[i]+len(blocks[i])):out[u]|=mask
    return tuple(out),offsets


def extend(a,x,c):
    n=len(a);r=list(a)+[1<<(n+1),0]
    for i in range(n):
        if x>>i&1:r[i]|=1<<n
        else:r[n]|=1<<i
        if c>>i&1:r[i]|=1<<(n+1)
        else:r[n+1]|=1<<i
    return tuple(r)


def lift(a):
    """Fixed strong ambient quotient z->p->q->z."""
    n=len(a);r=list(a)+[1<<(n+1),(1<<n)-1]
    for i in range(n):r[i]|=1<<n
    return tuple(r)


def block_tests():
    specs=[(1,1),(2,1),(2,2),(3,2),(2,2,1),(2,1,2),(2,2,2),(3,1,2),(1,1,1,1),(2,1,1,1)]
    records=[];profile_checks=0;queries=0
    for sizes in specs:
        k=len(sizes)
        for variant in range(2):
            blocks=[tournament(s,(17*i+variant*5+1)%(1<<(s*(s-1)//2))) for i,s in enumerate(sizes)]
            replacements=[converse(a) if (i+variant)%2 else a for i,a in enumerate(blocks)]
            need(all(profile(a)==profile(b) for a,b in zip(blocks,replacements)),'equal replacement profiles')
            T,offsets=substitute(blocks,tuple(range(k)))
            for pi in permutations(range(k)):
                if pi==tuple(range(k)):continue
                U,_=substitute(replacements,pi)
                t=next(t for t in range(1,k) if set(pi[:t])!=set(range(t)))
                prefix=sum(((1<<sizes[i])-1)<<offsets[i] for i in range(t))
                X=extend(T,0,prefix);Y=extend(U,0,prefix)
                p,q=profile(X),profile(Y)
                need(all(b>=a for a,b in zip(p,q)),'complete profile domination')
                bound=2
                for a in blocks:bound*=hp(a)
                need(q[0]-p[0]>=bound,'strict local quantitative bound')
                # The old exterior q points to the module, p loses to it.
                m=len(T);a_mask=1<<(m+1);c_mask=prefix|a_mask
                f=hp(extend(lift(T),a_mask,c_mask));g=hp(extend(lift(U),a_mask,c_mask))
                need(g-f>=3*bound,'strong ambient quantitative bound')
                need(hp(T)==hp(U),'old scalar equality')
                records.append(dict(sizes=sizes,variant=variant,order=pi,cut=t,local_difference=q[0]-p[0],local_lower_bound=bound,ambient_difference=g-f))
                profile_checks+=len(p);queries+=4
    return dict(cases=len(records),profile_coordinates=profile_checks,literal_queries=queries,records=records)


def modules(a):
    n=len(a)
    return [m for m in range(1,1<<n) if m.bit_count()>=2 and all(a[v]&m in (0,m) for v in range(n) if not m>>v&1)]


def maximal(ms):return [m for m in ms if not any(m!=n and m&n==m for n in ms)]

def strong(a):
    n=len(a);full=(1<<n)-1
    for v in range(n):
        r=1<<v
        while True:
            s=r
            for u in range(n):
                if r>>u&1:s|=a[u]
            if s==r:break
            r=s
        if r!=full:return False
    return True


def toggle(a,ms):
    r=list(a)
    for m in ms:
        for u,v in combinations(range(len(a)),2):
            if m>>u&1 and m>>v&1:r[u]^=1<<v;r[v]^=1<<u
    return tuple(r)


def components(a):
    n=len(a);reach=[]
    for v in range(n):
        r=1<<v
        while True:
            s=r
            for u in range(n):
                if r>>u&1:s|=a[u]
            if s==r:break
            r=s
        reach.append(r)
    unseen=(1<<n)-1;cs=[]
    while unseen:
        v=(unseen&-unseen).bit_length()-1;m=sum(1<<u for u in range(n) if reach[v]>>u&1 and reach[u]>>v&1);cs.append(m);unseen^=m
    return sorted(cs,key=lambda m:-sum(bool(a[(m&-m).bit_length()-1]&n) for n in cs if n!=m))


def laminar_tests():
    root=Path(__file__).resolve().parent.parent/'simultaneous_response';count=0;queries=0;records=[];all_covered=0
    for n in range(2,6):
        bank=json.loads((root/f'full_response_n{n}.json').read_text());basis=bank['basis'];vectors=[r['vector'] for r in bank['records']]
        for code in range(1<<(n*(n-1)//2)):
            T=tournament(n,code);mods=modules(T)
            for z in range(1,1<<len(mods)):
                fs=[m for j,m in enumerate(mods) if z>>j&1]
                if any(a&b and a&b not in (a,b) for a,b in combinations(fs,2)):continue
                maxs=maximal(fs)
                if any(sum(maximal([a for a in fs if a!=B and a&B==a]))!=B for B in maxs):continue
                all_covered+=1
                candidates=[B for B in maxs if not strong(induced(T,[v for v in range(n) if B>>v&1]))]
                if not candidates:continue
                B=candidates[0];children=maximal([a for a in fs if a!=B and a&B==a]);reps=[(a&-a).bit_length()-1 for a in children];K=induced(T,reps);cs=components(K)
                need(len(cs)>=2,'nonstrong quotient')
                groups=[sum(children[i] for i in range(len(children)) if c>>i&1) for c in cs]
                prefix=groups[0];v=(B&-B).bit_length()-1
                outside=sum(1<<u for u in range(n) if not B>>u&1 and T[u]>>v&1)
                a=outside;c=outside|prefix;U=toggle(T,fs)
                need(all(g in modules(T) and g in modules(U) for g in groups),'persistent SCC groups')
                need(all(profile(induced(T,[v for v in range(n) if g>>v&1]))==profile(induced(U,[v for v in range(n) if g>>v&1])) for g in groups),'SCC-group profile equality')
                fa=hp(extend(T,a,c));fb=hp(extend(U,a,c));need(fb>fa,'covered laminar strict query')
                # Check every induced-subset H explicitly, independently of the bank.
                need(hp_subsets(T)==hp_subsets(U),'hereditary scalar equality')
                def evaluate(vec):return sum(z for z,(s,t) in zip(vec,basis) if s&a==s and t&c==t)
                ucode=sum(1<<b for b,(i,j) in enumerate(combinations(range(n),2)) if U[i]>>j&1)
                need(evaluate(vectors[code])==fa and evaluate(vectors[ucode])==fb,'complete vector query agreement')
                records.append(dict(n=n,source_code=code,modules=fs,maximal_module=B,ordered_groups=groups,query=[a,c],values=[fa,fb]))
                count+=1;queries+=2
    need(count==all_covered,'all covered families through five have a nonstrong maximal module')
    return dict(covered_families=all_covered,nonstrong_maximal_families=count,literal_queries=queries,records=records)


def main():return dict(status='PASS',ordered_blocks=block_tests(),covered_laminar=laminar_tests())
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,separators=(',',':')))
