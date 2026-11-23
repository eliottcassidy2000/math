#!/usr/bin/env python3
"""Independent integer verification of off-ray prime/cofactor module witnesses."""
from math import isqrt

def counts(rows):
    n=len(rows); N=1<<n
    dp=[[0]*n for _ in range(N)]
    for i in range(n): dp[1<<i][i]=1
    for s in range(1,N):
        for i,v in enumerate(dp[s]):
            if v:
                nxt=rows[i]&~s
                while nxt:
                    b=nxt&-nxt; nxt-=b
                    dp[s|b][b.bit_length()-1]+=v
    return [1]+[sum(row) for row in dp[1:]]

def star(core,tags,y):
    tape=tuple(tags>>i&1 for i in range(len(core)))+tuple(1-int(c) for c in bin(y)[2:][::-1]) if y else tuple(tags>>i&1 for i in range(len(core)))
    n=len(tape); c_mask=(1<<len(core))-1
    h=[((1<<n)-1)^((1<<(i+1))-1) for i in range(n)]+[0]
    for i,row in enumerate(core): h[i]=(h[i]&~c_mask)|row
    for i,b in enumerate(tape):
        if b: h[n]|=1<<i
        else: h[i]|=1<<n
    return h

def insert(host,m,module):
    n=len(host); k=len(module); total=n+k-1
    def label(i): return (i,) if i<m else tuple(range(m,m+k)) if i==m else (i+k-1,)
    out=[0]*total
    for i,row in enumerate(host):
        for j in range(n):
            if row>>j&1:
                for u in label(i):
                    for v in label(j):
                        out[u]|=1<<v
    for i,row in enumerate(module): out[m+i]|=row<<m
    return out

def flip(rows,i,j):
    out=list(rows);out[i]^=1<<j;out[j]^=1<<i
    return out

examples=[
    ((226,100,137,19,39,76,156,58),661,8848,354581),
    ((38,60,56,49,33,0),11,200,8011),
    ((22,60,56,49,32,1),19,256,10259),
    ((28,53,56,50,32,1),23,272,10903),
    ((14,60,56,48,33,1),23,296,11863),
    ((26,44,57,48,34,1),31,328,13151),
    ((22,44,56,49,34,1),37,376,15077),
]
p0=star((8,5,9,18,7),11,10)
c0=star((16,9,67,21,38,15,59),68,3)
def strong(rows):
    full=(1<<len(rows))-1
    for start in range(len(rows)):
        seen=1<<start
        while True:
            nxt=seen
            for i,row in enumerate(rows):
                if seen>>i&1: nxt|=row
            if nxt==seen: break
            seen=nxt
        if seen!=full: return False
    return True

assert all(strong(r) for r in (p0,c0,flip(c0,2,4),flip(c0,0,2)))
for module,D,Q,p in examples:
    cnt=counts(module); f=len(cnt)-1
    assert p==40*Q+D and all(p%d for d in range(2,isqrt(p)+1))
    assert cnt[-1]==D and sum(cnt[s]*cnt[f^s] for s in range(f+1))==Q
    parent=counts(insert(p0,5,module))[-1]
    center=counts(insert(c0,7,module))[-1]
    lower=counts(insert(flip(c0,2,4),7,module))[-1]
    upper=counts(insert(flip(c0,0,2),7,module))[-1]
    assert (parent,lower,center,upper)==(7*p,21*p-2*D,21*p,21*p+2*D)
    print(D,Q,p,parent,lower,center,upper)

word="OOIOIOIOI"
tag=sum((c=='O')<<i for i,c in enumerate(word))
core=tuple(((1<<len(word))-1)^((1<<(i+1))-1) for i in range(len(word)))
assert counts(star(core,tag,0))[-1]==287
base=star((16,9,67,21,38,15,59),68,1)
assert [counts(A)[-1] for A in (flip(base,2,4),base,flip(base,0,2))]==[859,861,863]
print("dyadic s=0 parent 287 with common-center fork 859,861,863 verified")
