#!/usr/bin/env python3
"""Exact free-first-word + up to two Wave-20 readers, compared to capped holes."""
from pathlib import Path
from local_11 import I, J, S0, WORDS, g, mat, mm, mv
from collections import deque
BASE=Path(__file__).resolve().parent

def source_g(cap):
    a=[1]*(cap+1)
    for x in range(1,len(a)):
        z,b=divmod(x,2);a[x]=a[z]+(0 if b else z+(z&1))
    return a

def reader_rows():
    rows={I[:3]:()}
    for b in WORDS:
        rows.setdefault(mat(b)[:3],(b,))
    for b in WORDS:
        for c in WORDS:
            rows.setdefault(mm(mat(c),mat(b))[:3],(b,c))
    return rows

def free_bridge(misses):
    mx=max(misses); asc=source_g(2*mx); unresolved=set(misses); found={}
    for z in range(0,2*mx+1,2):
        h=asc[z]
        if h in unresolved:
            unresolved.remove(h);found[h]=(z,())
    counts=[len(found)]
    for n in (1,2):
        for (a,b,c),r in reader_rows().items():
            if len(r)!=n:continue
            for z in range((mx+1)//b+1):
                h=a*asc[z]+b*(z+1)+c*(1-(z&1))
                if h in unresolved:
                    unresolved.remove(h);found[h]=(z,r)
        counts.append(len(found))
    return found,sorted(unresolved),counts

def low_ports(cap=961):
    p=[tuple(map(int,l.split())) for l in (BASE/'ports13.txt').read_text().splitlines()]
    q=deque([S0]); seen={S0}; val={1}
    q=deque([(1,0,1)]);seen={(1,0,1)}
    while q:
        u=q.popleft()
        for m in p:
            v=mv(m,u)
            if v[0]<=cap and v not in seen:
                seen.add(v);q.append(v);val.add(v[0])
    return val

SMALL={
    67:((0,1,11,19,39,15),28),
    163:((0,1,3,23,39,79,31),108),
    173:((0,1,11,35,15,87,31),60),
}

def literal_small(N):
    r,U=SMALL[N]
    n=len(r)
    r=tuple(v|((1-(U>>i&1))<<n) for i,v in enumerate(r))+(U,)
    M=len(r);F=(1<<M)-1
    D=[[0]*M for _ in range(1<<M)]
    for i in range(M):D[1<<i][i]=1
    for j,W in enumerate(D):
        for i,val in enumerate(W):
            d=r[i]&~j
            while val and d:
                k=d&-d;d^=k
                D[j|k][k.bit_length()-1]+=val
    assert sum(D[F])==N

def main():
    misses=list(map(int,(BASE/'orbit_head_misses.txt').read_text().split()))
    assert len(misses)==408 and max(misses)==418429
    found,rest,ct=free_bridge(misses)
    assert ct==[207,299,363] and len(rest)==45 and max(rest)==139201
    top=low_ports()
    last=[n for n in rest if n<=961 and n not in top]
    assert last==[7,21,67,163,173]
    assert all(n<=200000 for n in rest)
    for z in SMALL:literal_small(z)
    print('PASS:',ct,'patched by free prefix; remaining',last,'after certified 13-port head')
    if len(__import__('sys').argv)>1:
        for h,(z,w) in sorted(found.items()):print(h,z,'/'.join(v or 'e' for v in w))
if __name__=='__main__':main()
