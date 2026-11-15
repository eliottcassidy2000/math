#!/usr/bin/env python3
"""Exact absolute host selection from a THM-5192 carrier bank.

Usage: python bank_absolute_host.py T [forms_orders3_5.txt]
Without a bank, use the order-five mask-zero carrier only.
"""
from itertools import combinations
from math import isqrt
import sys

def p(n):
    q=0
    while n:
        b,n=n&1,n>>1
        if not b: q+=n+(n&1)
    return q

def reverse(n,k):
    r=0
    for _ in range(k): r,n=2*r+(n&1),n>>1
    return r

def ceildiv(a,b): return -(-a//b)

def inverse(T,form):
    k,q,s,e,t,f,n=form
    assert k>0 and q>0 and min(s,e,t,f,n)>=0
    for a in range(1,isqrt(T//q)+1):
        pa=p(a)
        for x_is_small in (True,False):
            if x_is_small:
                c,phi,Q=q*a+s,e,T-(t*a+k*pa+n+f*(a&1))
            else:
                c,phi,Q=q*a+t,f,T-(s*a+k*pa+n+e*(a&1))
            lo=max(a+(not x_is_small),ceildiv(Q-phi,c+k))
            hi=Q//c
            if lo>hi: continue
            j=max(1,(isqrt(T)//a).bit_length()-1)
            K=1<<j
            Rmax=(c+k)*(K-1)+phi+k
            for w in range(lo//K,hi//K+1):
                L=(c*K+k*(K-1))*w+k*p(w)
                if w:
                    rlo=max(0,ceildiv(L-Q,k*w)); rhi=min(K-1,(L+Rmax-Q)//(k*w))
                    candidates=(K*w+reverse(r,j) for r in range(rlo,rhi+1))
                else:
                    candidates=range(lo,min(hi,K-1)+1)
                for b in candidates:
                    if lo<=b<=hi and c*b+k*p(b)+phi*(b&1)==Q:
                        yield (a,b) if x_is_small else (b,a)

def opposite(out):
    return [sum(1<<j for j,b in enumerate(out) if b>>i&1) for i in range(len(out))]

def dp(out):
    n=len(out)
    D=[[0]*n for _ in range(1<<n)]
    for i in range(n): D[1<<i][i]=1
    for S in range(1,1<<n):
        for i,z in enumerate(D[S]):
            if not z: continue
            A=out[i]&~S
            while A:
                a=A&-A; A-=a
                D[S|a][a.bit_length()-1]+=z
    return D

def H(out): return sum(dp(out)[-1])

def append_apex(out,mask):
    N=len(out)
    return [r if mask>>i&1 else r|(1<<N) for i,r in enumerate(out)]+[mask]

def compile_host(x,y,out,mask):
    """Literal outgoing adjacency rows, with apex last."""
    bx=bin(x)[3:]; by=bin(y)[3:][::-1]
    i0=1+len(bx)
    marks=[1]+list(map(int,bx))+[(mask>>i)&1 for i in range(len(out))]+[1-int(b) for b in by]+[0]
    V=len(marks)
    rows=[0]*V
    for i,j in combinations(range(V),2):
        if i0<=i<j<i0+len(out):
            forward=out[i-i0]>>(j-i0)&1
        else:
            forward=True
        if forward: rows[i]|=1<<j
        else: rows[j]|=1<<i
    return append_apex(rows,sum(mark<<i for i,mark in enumerate(marks)))

MASK_ZERO=(11,100,0,0,89,11,11)
MASK_ZERO_WITNESS=([24,1,3,6,14],0)

def connected(out):
    want=(1<<len(out))-1; seen=1
    while True:
        nxt=seen
        for i,r in enumerate(out):
            if seen>>i&1: nxt|=r
        if nxt==seen: return seen==want
        seen=nxt

def recover_core(form,order):
    """Recover a marked labelled core from transfer coordinates; order <=5."""
    if form==MASK_ZERO: return MASK_ZERO_WITNESS
    a,q,s,e,t,f,n=form
    pairs=list(combinations(range(order),2)); V=(1<<order)-1
    for bits in range(1<<len(pairs)):
        out=[0]*order
        for k,(i,j) in enumerate(pairs):
            if bits>>k&1: out[i]|=1<<j
            else: out[j]|=1<<i
        if not connected(out) or not connected(opposite(out)): continue
        end=dp(out)
        if sum(end[-1])!=a: continue
        hs=[1]+[sum(row) for row in end[1:]]
        if sum(hs[S]*hs[V^S] for S in range(V+1))!=q: continue
        start=dp(opposite(out))
        def reads(D):
            plain=[0]*order; signed=[0]*order
            for S in range(1,V+1):
                sig=1 if S.bit_count()%2 else -1
                for i,z in enumerate(D[S]):
                    w=z*hs[V^S]
                    plain[i]+=w; signed[i]+=sig*w
            return plain,signed
        R,E=reads(start); T,F=reads(end)
        def subset_sum(vals,mask):
            return sum(v for i,v in enumerate(vals) if mask>>i&1)
        for mask in range(V+1):
            if (subset_sum(R,mask),subset_sum(E,mask),subset_sum(T,V^mask),subset_sum(F,V^mask))!=(s,e,t,f): continue
            if H(append_apex(out,mask))==n: return out,mask
    raise ValueError("No core of the indicated order has these transfer coordinates")

def read_bank(path):
    with open(path) as ff:
        for line in ff:
            data=line.split('#',1)[0].strip()
            if data:
                parts=data.split('|')
                form=tuple(map(int,parts[0].split()))
                order=int(parts[1].split()[0])
                yield form,order

def find_host(T,bank):
    for form,order in bank:
        xy=next(inverse(T,form),None)
        if xy is not None:
            out,mark=recover_core(form,order)
            return xy,form,out,mark,compile_host(*xy,out,mark)
    return None

def replay(bank):
    assert H(MASK_ZERO_WITNESS[0])==11
    for x,y in [(1,1),(2,3),(7,6),(15,29),(24,17)]:
        F=100*x*y+89*x+11*(p(x)+p(y)+1+(x&1))
        assert (x,y) in set(inverse(F,MASK_ZERO))
        assert H(compile_host(x,y,*MASK_ZERO_WITNESS))==F
    if bank is not None:
        for T in (333,497,1105,8899):
            result=find_host(T,list(read_bank(bank)))
            assert result and H(result[-1])==T, (T,result)
    print("Literal tournament assembly and exact prescribed-count controls passed")

if __name__=="__main__":
    bank_path=sys.argv[2] if len(sys.argv)>2 else None
    if sys.argv[1]=="--replay":
        replay(bank_path)
    else:
        T=int(sys.argv[1])
        if T<1: raise ValueError("target must be positive")
        bank=list(read_bank(bank_path)) if bank_path else [(MASK_ZERO,5)]
        result=find_host(T,bank)
        if result:
            (x,y),form,core,mark,adj=result
            print("X,Y:",x,y,"form:",form,"core:",core,"mark:",mark)
            print("Order:",len(adj),"outgoing adjacency rows:",adj)
        else:
            print("No exact host in the specified carrier bank")
