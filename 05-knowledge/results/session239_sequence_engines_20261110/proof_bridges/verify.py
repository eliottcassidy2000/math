#!/usr/bin/env python3
"""Exact independent controls for value-refined fixed-permutation SCC calculus."""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import comb, factorial, gcd, lcm
import json


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def partitions(n, lo=1):
    if not n:
        yield ()
    for a in range(lo, n + 1):
        for tail in partitions(n-a, a):
            yield (a,) + tail


def permutation(parts):
    p=[]
    for d in parts:
        start=len(p)
        p.extend(start+(j+1)%d for j in range(d))
    return tuple(p)


def edges(n):
    return [(i,j) for i in range(n) for j in range(i+1,n)]


def relabel(mask, p):
    es=edges(len(p)); index={e:i for i,e in enumerate(es)}; out=0
    for k,(a,b) in enumerate(es):
        u,v=p[a],p[b]; bit=((mask>>k)&1) ^ (u>v)
        if bit: out |= 1<<index[min(u,v),max(u,v)]
    return out


def rows(n,mask):
    out=[0]*n
    for k,(i,j) in enumerate(edges(n)):
        a,b=(i,j) if (mask>>k)&1 else (j,i)
        out[a] |= 1<<b
    return out


def path_count(out):
    n=len(out)
    if not n: return 1
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1,1<<n):
        for v in range(n):
            if dp[mask][v]:
                rem=out[v]&~mask
                while rem:
                    bit=rem&-rem; rem-=bit; w=bit.bit_length()-1
                    dp[mask|bit][w]+=dp[mask][v]
    return sum(dp[-1])


def scc_count(out):
    n=len(out); reach=[out[i]|(1<<i) for i in range(n)]
    for k in range(n):
        for i in range(n):
            if (reach[i]>>k)&1: reach[i]|=reach[k]
    return len({tuple(j for j in range(n) if (reach[i]>>j)&1 and (reach[j]>>i)&1) for i in range(n)})


def fixed_masks(parts):
    """Enumerate one binary choice per non-reversing unordered-pair orbit."""
    p=permutation(parts); es=edges(len(p)); idx={e:k for k,e in enumerate(es)}
    unseen=set(range(len(es))); orbits=[]
    while unseen:
        first=min(unseen); a,b=es[first]; directed=[]
        while True:
            k=idx[min(a,b),max(a,b)]
            if directed and k==first:
                if (a,b)!=es[first]: return
                break
            directed.append((k,int(a<b))); unseen.discard(k)
            a,b=p[a],p[b]
        orbits.append(directed)
    for bits in product((0,1),repeat=len(orbits)):
        mask=0
        for flip,orb in zip(bits,orbits):
            for k,sign in orb: mask |= (sign^flip)<<k
        yield mask


def z(parts):
    c=Counter(parts); ans=1
    for d,m in c.items(): ans*=d**m*factorial(m)
    return ans


def subprofiles(parts):
    c=Counter(parts); ds=sorted(c)
    for ks in product(*(range(c[d]+1) for d in ds)):
        sub=tuple(d for d,k in zip(ds,ks) for _ in range(k))
        rest=tuple(d for d,k in zip(ds,ks) for _ in range(c[d]-k))
        mult=1
        for d,k in zip(ds,ks): mult*=comb(c[d],k)
        yield sub,rest,mult


def convolution(a,b):
    out=Counter()
    for h,x in a.items():
        for k,y in b.items(): out[h*k]+=x*y
    return out


def subtract(a,b,scale=1):
    for h,c in b.items(): a[h]-=scale*c
    return Counter({h:c for h,c in a.items() if c})


def scalar_counts(N):
    A=[]
    for n in range(N+1):
        val=Fraction()
        for part in partitions(n):
            if any(d%2==0 for d in part): continue
            e=(sum(gcd(i,j) for i in part for j in part)-len(part))//2
            val+=Fraction(2**e,z(part))
        require(val.denominator==1,'Davis integrality')
        A.append(int(val))
    S=[0]*(N+1); L=[2**comb(n,2) for n in range(N+1)]; LS=[0]*(N+1)
    for n in range(1,N+1):
        S[n]=A[n]-sum(S[k]*A[n-k] for k in range(1,n))
        LS[n]=L[n]-sum(comb(n,k)*LS[k]*L[n-k] for k in range(1,n))
    return A,S,LS


def main():
    F={():Counter({1:1})}; G={():Counter()}; output={'scope':'all labelled tournaments and all permutation cycle types through order 6'}
    class_tables={}; checks=0
    for n in range(1,7):
        facts={mask:(path_count(rows(n,mask)),scc_count(rows(n,mask))) for mask in range(1<<comb(n,2))}
        fixed_class=Counter(); fixed_strong=Counter()
        for part in partitions(n):
            p=permutation(part)
            direct=Counter(h for mask,(h,r) in facts.items() if relabel(mask,p)==mask)
            direct_s=Counter(h for mask,(h,r) in facts.items() if r==1 and relabel(mask,p)==mask)
            generated=list(fixed_masks(part))
            F[part]=Counter(facts[m][0] for m in generated)
            require(F[part]==direct,('fixed orientation generator',part))
            G[part]=F[part].copy()
            for sub,rest,mult in subprofiles(part):
                if not sub or not rest: continue
                G[part]=subtract(G[part],convolution(G[sub],F[rest]),mult)
            require(G[part]==direct_s,('SCC profile recurrence',part))
            for h,c in F[part].items(): fixed_class[h]+=Fraction(c,z(part))
            for h,c in G[part].items(): fixed_strong[h]+=Fraction(c,z(part))
            order=lcm(*part)
            require(all(h%order==0 for h in F[part]),('free action divisibility',part))
            for target in (1,3,5,7,9,15,21,25,35,49,63):
                cut=lambda c:Counter({h:v for h,v in c.items() if target%h==0})
                got=cut(F[part])
                for sub,rest,mult in subprofiles(part):
                    if sub and rest: got=subtract(got,cut(convolution(cut(G[sub]),cut(F[rest]))),mult)
                require(got==cut(G[part]),('divisor truncation',part,target)); checks+=1
        # Independent quotient: connected components under adjacent transpositions.
        unseen=set(facts); cls=Counter(); strong=Counter(); by_scc=Counter()
        swaps=[]
        for i in range(n-1):
            p=list(range(n)); p[i],p[i+1]=p[i+1],p[i]; swaps.append(p)
        while unseen:
            seed=min(unseen); unseen.remove(seed); stack=[seed]; h,r=facts[seed]
            cls[h]+=1; by_scc[r]+=1
            if r==1: strong[h]+=1
            while stack:
                mask=stack.pop()
                for p in swaps:
                    nxt=relabel(mask,p)
                    if nxt in unseen: unseen.remove(nxt); stack.append(nxt)
        require(cls==fixed_class,('weighted Burnside vs orbit BFS',n))
        require(strong==fixed_strong,('strong weighted Burnside vs orbit BFS',n))
        class_tables[n]={'classes':sum(cls.values()),'strong_classes':sum(strong.values()),'H_histogram':dict(sorted(cls.items())),'strong_H_histogram':dict(sorted(strong.items())),'SCC_histogram':dict(sorted(by_scc.items()))}
    A,S,LS=scalar_counts(20)
    for n,t in class_tables.items():
        require(A[n]==t['classes'] and S[n]==t['strong_classes'],'scalar specialization')
    # Hostile: ordinary polynomial multiplication adds H, rather than multiplying it.
    require(path_count(rows(3,2))==3,'C3 encoding')
    require(F[(3,3)][9]>0,'two cyclic orbits realize H9')
    require(9!=3+3,'hostile additive product')
    # Hostile: SCC words are ordered, not multisets. Triangle and singleton give 2 classes.
    require(class_tables[4]['H_histogram'][3]==2,'ordered SCC word hostile')
    require(F[(2,)]==Counter(),'reversing pair even-cycle hostile')
    # Exact cyclic minima for 3,5,7, independent of the order<=6 census.
    minima={}
    for d in (3,5,7):
        vals=Counter(path_count(rows(d,m)) for m in fixed_masks((d,)))
        minima[d]={'minimum':min(vals),'histogram':dict(sorted(vals.items()))}
    require([minima[d]['minimum'] for d in (3,5,7)]==[3,15,175],'inherited cyclic minima')
    for part,c in F.items():
        if not part or not c: continue
        lower=1
        for d in part: lower*=minima[d]['minimum'] if d>1 else 1
        require(min(c)>=lower,('cyclic profile floor',part))
    output.update(class_tables=class_tables,divisor_truncation_checks=checks,unlabelled=A,strong_unlabelled=S,strong_labelled=LS,cyclic_minima=minima,status='PASS')
    print(json.dumps(output,indent=2,sort_keys=True))

if __name__=='__main__': main()
