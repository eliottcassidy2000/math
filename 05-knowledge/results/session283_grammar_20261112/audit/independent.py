#!/usr/bin/env python3
"""Independent matrix/permutation checks; does not import producer algorithms."""
from itertools import permutations, combinations, product
from collections import Counter
from fractions import Fraction
from math import comb, factorial
import json

def need(p,s):
    if not p: raise RuntimeError(s)

def graph(n,k):
    a=[[False]*n for _ in range(n)]
    for b,(u,v) in enumerate(combinations(range(n),2)):
        a[u][v]=bool((k>>b)&1);a[v][u]=not a[u][v]
    return a

def hp(a):
    n=len(a);d={}
    for u in range(n):d[(1<<u,u)]=1
    for mask in range(1,1<<n):
        for u in range(n):
            c=d.get((mask,u),0)
            if c:
                for v in range(n):
                    if not (mask>>v&1) and a[u][v]:
                        key=(mask|1<<v,v);d[key]=d.get(key,0)+c
    return sum(d.get(((1<<n)-1,u),0) for u in range(n))

def cover(a):
    n=len(a);out=[0]*n
    for p in permutations(range(n)):
        bad=sum(not a[p[i]][p[i+1]] for i in range(n-1))
        for c in range(bad,n):out[c]+=comb(n-1-bad,c-bad)
    return out

def surj(m,s):
    if s==0:return int(m==0)
    return sum((-1)**j*comb(s,j)*(s-j)**m for j in range(s+1))

def weights(n,m):
    # Enumerate actual quotient run words, not the closed coefficient formula.
    out=[]
    for r in range(1,n+1):
        value=0
        for start in range(3):
            for length in range(max(1,3*r-2),3*r+3):
                counts=[0,0,0]
                for i in range(length):counts[(start+i)%3]+=1
                if counts[0]==r and counts[1]>0 and counts[2]>0:
                    value+=surj(m,counts[1])*surj(m,counts[2])
        out.append(value)
    return out

def caps(a,m):
    n=len(a);N=n+2*m;out=[[False]*N for _ in range(N)]
    blocks=[range(n),range(n,n+m),range(n+m,N)]
    for i in range(n):
        for j in range(n):out[i][j]=a[i][j]
    for block in blocks[1:]:
        for i in block:
            for j in block:out[i][j]=i<j
    for k,block in enumerate(blocks):
        for i in block:
            for j in blocks[(k+1)%3]:out[i][j]=True
    return out

def profile(a):
    n=len(a);p=Counter()
    for u,v in combinations(range(n),2):
        d=sum(a[u][w]!=a[v][w] for w in range(n) if w not in (u,v));p[d]+=1
    return p

def flip(a,u,v):
    b=[r[:] for r in a];b[u][v]=not b[u][v];b[v][u]=not b[v][u];return b

def rank(a):
    # Compute components, independently of the n-e0 formula under review.
    n=len(a);seen=set();num=0
    for u in range(n):
        if u in seen:continue
        num+=1;seen.add(u);stack=[u]
        while stack:
            v=stack.pop()
            for w in range(n):
                if w not in seen and all(a[v][x]==a[w][x] for x in range(n) if x not in (v,w)):
                    seen.add(w);stack.append(w)
    return num

def inflate(q,lens):
    vertex=[(i,j) for i,l in enumerate(lens) for j in range(l)]
    return [[(x<y if i==j else q[i][j]) for j,y in vertex] for i,x in vertex]

def main():
    counts=Counter();profiles={}
    for n in range(1,6):
        for code in range(1<<comb(n,2)):
            a=graph(n,code);p=profile(a);b=rank(a);need(b==n-p[0],'forest identity')
            counts['literal_rank_forests']+=1
            jumps=Counter()
            for u,v in combinations(range(n),2):
                db=rank(flip(a,u,v))-b;jumps[db]+=1
                need(abs(db)<=4,'jump bound');counts['literal_arc_reversals']+=1
            need(sum(k*v for k,v in jumps.items())==2*((n-2)*(n-b)-p[1]),'mean drift')
            if n<=4:
                f=cover(a);profiles.setdefault((n,tuple(f)),a)
                for m in (1,2):
                    w=weights(n,m);need(hp(caps(a,m))==sum(x*y for x,y in zip(f,w)),'physical cap')
                    counts['physical_caps']+=1
    for (n,f),a in profiles.items():
        B=factorial(n)*2**(n-1);L=(72*(n+1)*B-1).bit_length();m=n*L;w=weights(n,m)
        need(all(w[r]>B*sum(w[:r]) for r in range(n)),'digit inequalities')
        z=sum(x*y for x,y in zip(f,w));rec=[]
        for x in reversed(w):q,z=divmod(z,x);rec.append(q)
        need(tuple(rec[::-1])==f and z==0,'physical profile decode');counts['independent_decodings']+=1
    # Inherited H9 pair: independent permutation-cut profiles and whole extension menus.
    menus=[]
    for rows in [(14,28,24,16,1),(10,28,25,16,1)]:
        a=[[bool(rows[i]>>j&1) for j in range(5)] for i in range(5)]
        need(cover(a)==[9,66,186,240,120],'H9 complete profile')
        menu=set()
        for cut in range(32):
            x=[r+[not bool(cut>>i&1)] for i,r in enumerate(a)]
            x.append([bool(cut>>i&1) for i in range(5)]+[False])
            menu.add(hp(x))
        menus.append(sorted(menu))
    need(menus==[[9,15,17,19,23,25,27,29,33,37,41],[9,15,17,19,23,29,31,33,37,43]],'H9 growth menus')
    counts['independent_H9_attachment_counts']=64
    # Bulk recurrence checked on all weight vectors2..4 for the triangle core.
    q=[[False,True,False],[False,False,True],[True,False,False]]
    for lens in product(range(2,5),repeat=3):
        a=inflate(q,lens);n=len(a);b=3;t2=lens.count(2)
        actual=Counter(rank(flip(a,u,v))-b for u,v in combinations(range(n),2))
        expected=+Counter({0:n-b,2:2*b*b-b-t2,3:2*b*(n-2*b-1)+2*t2,4:comb(n-2*b,2)-(n-3*b+t2)})
        need(actual==expected,'universal bulk kernel');counts['weighted_bulk_families']+=1
    # Spectral expected rank from complete trajectory distributions.
    for n in (3,4):
        M=comb(n,2);N=n-2
        for code in range(1<<M):
            p=profile(graph(n,code));dist={code:1};den=1
            for t in range(4):
                actual=Fraction(sum(c*rank(graph(n,k)) for k,c in dist.items()),den)
                spectral=Fraction(n)
                for s in range(N+1):
                    Q=sum(c*sum((-1)**j*comb(d,j)*comb(N-d,s-j) for j in range(max(0,s-N+d),min(d,s)+1)) for d,c in p.items())
                    spectral-=Fraction(Q,2**N)*Fraction(M-4*s,M)**t
                need(actual==spectral,'spectral trajectory');counts['literal_trajectory_moments']+=1
                nxt=Counter()
                for k,c in dist.items():
                    for e in range(M):nxt[k^(1<<e)]+=c
                dist=nxt;den*=M
    print(json.dumps({'status':'PASS','independent_counts':dict(counts),'scope':'new paths/grammar/metagraph formulas; enumeration audited separately'},sort_keys=True,indent=2))
if __name__=='__main__':main()
