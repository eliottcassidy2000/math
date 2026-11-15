#!/usr/bin/env python3
from itertools import combinations
from collections import defaultdict

def rows(n,x):
    a=[set() for _ in range(n)]
    for k,(i,j) in enumerate(combinations(range(n),2)):
        if x>>k&1:a[i].add(j)
        else:a[j].add(i)
    return a

def covers(a):
    n=len(a); f=(1<<n)-1
    mask=[sum(1<<u for u in r) for r in a]
    d=[defaultdict(int) for _ in range(f+1)]
    for i in range(n):d[1<<i][i]=1
    for s in range(1,f+1):
        for i,x in d[s].items():
            m=mask[i]&~s
            while m:
                z=m&-m;m-=z
                d[s|z][z.bit_length()-1]+=x
    h=[1]+[sum(d[s].values())for s in range(1,f+1)]
    return h[f],sum(h[s]*h[f^s]for s in range(f+1)),d,h

def equiv(a):
    n=len(a)
    return [set(range(n))-{i}-a[i]for i in range(n)]

def join(a,b):
    n,m=len(a),len(b)
    return [r|set(range(n,n+m))for r in a]+[{u+n for u in r}for r in b]

def tt(n):
    return [set(range(i+1,n))for i in range(n)]

def subst(a,parts):
    starts=[];out=[]
    for p in parts:
        starts.append(len(out));o=len(out)
        out.extend({v+o for v in r}for r in p)
    for i,r in enumerate(a):
        for j in r:
            for k in range(starts[i],starts[i]+len(parts[i])):
                out[k].update(range(starts[j],starts[j]+len(parts[j])))
    return out

def close(a,O,x):
    n=len(a);b=join(a,x)
    for i in range(n):
        if i not in O:b[i].add(len(b))
    for i in range(n,len(b)):b[i].add(len(b))
    return b+[set(O)]

def finish():
    B=[{1},{2,4},{0,3,4},{0,1,7},{0,3,5,6},{0,1,2,3},{0,1,2,3,5},{0,1,2,4,5,6}]
    O={1,3,4}
    h,q,ends,hp=covers(B)
    _,_,starts,_=covers(equiv(B))
    full=(1<<len(B))-1
    M=sum(hp[s]*sum(starts[full^s][i] for i in O) for s in range(full))
    W=sum(sum(ends[s][i] for i in O)*sum(starts[full^s][i] for i in O)for s in range(1,full))
    assert (h,M,h-W)==(217,1024,-1)
    for seed in (rows(4,17),rows(5,548),rows(3,2),tt(1)):
        hh,qq,*_=covers(seed)
        for t in (0,1):
            assert covers(close(B,O,join(seed,tt(t))))[0]==1024*2**t*qq-hh

    # Derive this literal order-five kernel instead from the canonical shape
    # clos(TT2,both,TTA=>TTB). Ports are its two middle ordered slots.
    K=close(tt(2),{0,1},tt(2))
    a=rows(5,548)
    assert covers(a)[:2]==(11,100)
    for j in range(3):
        ap=[]
        for i in range(j):ap=join(ap,a)
        for t in ((1,2) if not j else (0,1)):
            mod=join(ap,tt(t))
            z=subst(K,[tt(1),tt(1),mod,tt(1),tt(1)])
            assert covers(z)[0]==6*100**j*2**t-11**j
    print('verified literal anti-cone pair, 8 varied module closes, and 6 H11 block hyperbolic witnesses')

if __name__=='__main__':finish()
