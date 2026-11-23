"""Exact minimum nonadaptive decoder within the 64 allowed one-extra-apex probes."""
from itertools import combinations
from collections import Counter

A8=(16,1,3,7,14)
A10=(20,1,2,7,14)

def apex(R,S):
    n=len(R)
    return tuple(r|(0 if S>>i&1 else 1<<n) for i,r in enumerate(R))+(S,)

def height(R):
    n=len(R)
    M=1<<n
    out=[[0]*n for _ in range(M)]
    inc=[sum(1<<i for i,r in enumerate(R) if r>>j&1) for j in range(n)]
    for u in range(n):out[1<<u][u]=1
    for s in range(1,M):
        if not s&(s-1):continue
        rem=s
        while rem:
            b=rem&-rem;rem-=b;z=s^b;v=b.bit_length()-1
            a=inc[v]&z
            if a:out[s][v]=sum(out[z][j] for j in range(n) if a>>j&1)
    return sum(out[-1])

def go():
    states=[apex(r,s) for r in (A8,A10) for s in range(32)]
    data=[[height(apex(r,m)) for r in states] for m in range(64)]
    blocks=[]
    for vals in data:
        blocks.append(sum(1<<i for i,(u,v) in enumerate(combinations(vals,2)) if u==v))
    hist=Counter();win=[]
    for k in range(1,5):
        for inds in combinations(range(64),k):
            bb=blocks[inds[0]]
            for j in inds[1:]:bb&=blocks[j]
            if not bb:win.append(inds)
            hist[k]+=not bb
        print(k,"successful sets",hist[k], "first",win[:3],flush=True)
        if win:break
    print("max distinct single-probe outputs",max(len(set(u)) for u in data))
    assert hist=={1:0,2:74}
    physical=[(a,b) for a,b in win if (a^b).bit_count()==1]
    print("successful physical pairs",physical)
    assert physical==[(2,34),(10,42),(23,55)]
    sorted_maps=[]
    for a,b in physical:
        D={}
        for i,(x,y) in enumerate(zip(data[a],data[b])):
            key=tuple(sorted((x,y)))
            D.setdefault(key,[]).append((8 if i<32 else 10,i%32,x,y))
        print("unoriented",a,b,len(D),"collisions",[(k,v) for k,v in D.items() if len(v)>1])
        sorted_maps.append(len(D))
    assert sorted_maps==[64,62,64]
    assert len(set(zip(data[2],data[34])))==64
    for key in ((2,34),(3,33),(5,31)):
        a,b=key; z=[v-u for u,v in zip(data[a],data[b])]
        print("masks",key,"readout range",min(data[a]+data[b]),max(data[a]+data[b]),
              "distinct signed deltas",len(set(z)),"gap overlap",
              sorted(set(z[:32])&set(z[32:])))
    print("chosen masks",(2,34),"example A8,0",data[2][0],data[34][0])
    # The three sets below describe everything lost when using just one scalar.
    p0=set(data[2]);p1=set(data[34])
    assert len(p0)<64 and len(p1)<64 and len(set(zip(data[2],data[34])))==64
    print("two-apex decoder and all 64 individual probes PASS")

if __name__=="__main__":go()
