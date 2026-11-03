#!/usr/bin/env python3
"""Independent subset-endpoint control for the interval-chain code."""
def build(bits):
    n=3; edges={(2,0)} # 0-based reversed edge of initial triangle
    for i,b in enumerate(bits,1):
        k=2*i+int(b)
        v=n+k
        edges.add((v,n-1))
        n=v+1
    return n,edges

def hamilton_count(n, backward):
    out=[]
    for i in range(n):
        row=0
        for j in range(n):
            if i!=j and ((i<j and (j,i) not in backward) or (i>j and (i,j) in backward)):
                row|=1<<j
        out.append(row)
    dp=[[0]*n for _ in range(1<<n)]
    for i in range(n): dp[1<<i][i]=1
    for s in range(1,1<<n):
        avail=((1<<n)-1)^s
        for i in range(n):
            if not dp[s][i]: continue
            z=out[i]&avail
            while z:
                l=z&-z; j=l.bit_length()-1
                dp[s|l][j]+=dp[s][i]; z-=l
    return sum(dp[-1])

def decode(h):
    q=h//3; out=[]
    while q>1:
        v=q-1; k=(v&-v).bit_length()-1
        i=len(out)+1; b=k-2*i
        if b not in (0,1): raise ValueError((h,k,i))
        out.append(str(b)); q//=1+(1<<k)
    return ''.join(out)

for bits in ('','0','1','00','01','10','11'):
    n, back=build(bits)
    product=3
    for i,b in enumerate(bits,1): product*=1+(1<<(2*i+int(b)))
    got=hamilton_count(n,back)
    if got!=product or decode(got)!=bits: raise AssertionError((bits,n,got,product,decode(got)))
    print(repr(bits),n,got)
print('self-converse controls')
for bits in ('', '0', '1'):
    ks=[2*i+int(b) for i,b in enumerate(bits,1)]
    chain=[1]+ks+list(reversed(ks))+[1]
    n=1; back=set()
    for k in chain:
        v=n+k
        back.add((v,n-1)); n=v+1
    expected=9
    for k in ks: expected *= (1+(1<<k))**2
    reflected={(n-1-j,n-1-i) for i,j in back}
    if reflected!=back: raise AssertionError(('reflection',bits,back,reflected))
    got=hamilton_count(n,back)
    if got!=expected: raise AssertionError(('self-converse',bits,got,expected))
    print(repr(bits),n,got)
from itertools import product
from math import isqrt
for m in range(9):
    for seq in product('01',repeat=m):
        bits=''.join(seq)
        ks=[2*i+int(b) for i,b in enumerate(bits,1)]
        chain=[1]+ks+list(reversed(ks))+[1]
        n=1; back=set(); P=1
        for k in ks:P*=1+(1<<k)
        for k in chain:
            v=n+k; back.add((v,n-1));n=v+1
        if {(n-1-j,n-1-i) for i,j in back}!=back: raise AssertionError(('reflection',bits))
        if n!=5+2*m*m+4*m+2*sum(map(int,bits)): raise AssertionError(('order',bits,n))
        h=9*P*P
        root=isqrt(h//9)
        if root!=P or decode(3*root)!=bits:raise AssertionError(('decode',bits))
print('PASS palindrome reflection and algebraic decoding for all words of length <=8')
