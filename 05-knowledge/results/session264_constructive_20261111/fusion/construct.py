#!/usr/bin/env python3
"""Root-amalgam constructor and exact finite controls. Standard library only."""
from itertools import combinations, product
from functools import lru_cache
from fractions import Fraction
from pathlib import Path
import json

@lru_cache(None)
def hp(T):
    n=len(T); size=1<<n
    dp=[[0]*n for _ in range(size)]
    for v in range(n): dp[1<<v][v]=1
    H=[1]+[0]*(size-1)
    for S in range(1,size):
        H[S]=sum(dp[S])
        for v in range(n):
            if not dp[S][v]:continue
            rem=T[v]&~S
            while rem:
                bit=rem&-rem; rem-=bit
                dp[S|bit][bit.bit_length()-1]+=dp[S][v]
    return tuple(H)

def op(T):return tuple(sum(1<<j for j in range(len(T)) if T[j]>>i&1) for i in range(len(T)))

def rooted(T):
    """Root is 0; return weighted counts by actual path vertex length."""
    n=len(T);H=hp(T);full=(1<<n)-1;dp=[[0]*n for _ in range(1<<n)];dp[1][0]=1
    P=[0]*(n+1)
    for S in range(1,1<<n):
        P[S.bit_count()]+=sum(dp[S])*H[full^S]
        for v in range(n):
            if not dp[S][v]:continue
            rem=T[v]&~S
            while rem:
                bit=rem&-rem;rem-=bit;dp[S|bit][bit.bit_length()-1]+=dp[S][v]
    return P

def profile(T):
    Hs=hp(T);full=(1<<len(T))-1;R=full^1;U=rooted(T);V=rooted(op(T));h=Hs[R]
    Q=k=0;S=R
    while True:
        val=Hs[S]*Hs[R^S];Q+=val;k+=(-1)**S.bit_count()*val
        if not S:break
        S=(S-1)&R
    d=sum((-1)**(j+1)*x for j,x in enumerate(U));e=sum((-1)**(j+1)*x for j,x in enumerate(V))
    return dict(H=Hs[full],h=h,Q=Q,k=k,U=sum(U),V=sum(V),d=d,e=e)

def matrix(p):
    h=p['h']
    return [[h,p['U']-h,p['d']-h,p['H']-h],[0,p['Q'],0,p['V']-h],[0,0,p['k'],p['e']-h],[0,0,0,h]]

def mul(A,B):return [[sum(A[i][k]*B[k][j] for k in range(4)) for j in range(4)] for i in range(4)]

def from_matrix(M):
    h=M[0][0]
    return dict(H=M[0][3]+h,h=h,Q=M[1][1],k=M[2][2],U=M[0][1]+h,V=M[1][3]+h,d=M[0][2]+h,e=M[2][3]+h)

def fuse(A,B):
    """Identify roots 0, orient all A-root vertices toward all B-root vertices."""
    n=len(A);m=len(B);out=[0]*(n+m-1)
    for T,inds in [(A,list(range(n))),(B,[0]+list(range(n,n+m-1)))]:
        for i in range(len(T)):
            for j in range(len(T)):
                if T[i]>>j&1:out[inds[i]]|=1<<inds[j]
    for i in range(1,n):
        for j in range(n,n+m-1):out[i]|=1<<j
    return tuple(out)

def tournaments(n):
    edges=list(combinations(range(n),2))
    for mask in range(1<<len(edges)):
        out=[0]*n
        for j,(a,b) in enumerate(edges):
            if mask>>j&1:out[a]|=1<<b
            else:out[b]|=1<<a
        yield tuple(out)

def power_formula(p,j):
    h=p['h'];rho=Fraction(p['Q'],h);x=Fraction(p['H']-h,h);b=Fraction(p['U']-h,h);c=Fraction(p['V']-h,h)
    delta=Fraction(p['d']-h,h);eta=Fraction(p['e']-h,h)
    assert j>=1 and len(str(h))>=1 and rho>1
    val=h**j*(1+j*x+b*c*(rho**j-j*rho+j-1)/(rho-1)**2+(j-1)*delta*eta)
    assert val.denominator==1
    return int(val)

def main():
    bank=[T for n in range(1,4) for T in tournaments(n)]
    controls=[]
    for A,B in product(bank,repeat=2):
        P=profile(fuse(A,B));expect=from_matrix(mul(matrix(profile(A)),matrix(profile(B))))
        assert P==expect,(A,B,P,expect)
        controls.append(dict(A=A,B=B,result=P))
    triple_bank=[(0,),(2,0),(0,1),(2,4,1)]
    for A,B,C in product(triple_bank,repeat=3):
        assert fuse(fuse(A,B),C)==fuse(A,fuse(B,C))
    towers=[]
    for A in bank[1:]:
        T=(0,);p=profile(A);seq=[]
        for j in range(1,6):
            T=fuse(T,A);actual=profile(T)['H'];closed=power_formula(p,j)
            assert actual==closed
            seq.append(actual)
        h=p['h'];q=p['Q']
        # Characteristic (lambda-q)(lambda-h)^2; k starts at 1.
        for j in range(3,len(seq)):
            assert seq[j]==(q+2*h)*seq[j-1]-(2*q*h+h*h)*seq[j-2]+q*h*h*seq[j-3]
        towers.append(dict(rows=A,profile=p,H=seq))
    # Pair of sources indistinguishable without parity, but distinguishable by a context.
    lookup={};hostile=None
    for n in range(1,6):
        for A in tournaments(n):
            p=profile(A);key=tuple(p[k] for k in ('H','h','Q','U','V'))
            if key in lookup:
                B=lookup[key];pb=profile(B)
                if (p['d'],p['e'])!=(pb['d'],pb['e']):
                    for C in bank:
                        for side in ['left','right']:
                            X=fuse(C,A) if side=='left' else fuse(A,C)
                            Y=fuse(C,B) if side=='left' else fuse(B,C)
                            if profile(X)['H']!=profile(Y)['H']:
                                hostile=dict(A=A,B=B,profile_A=p,profile_B=pb,context=C,side=side,H_A=profile(X)['H'],H_B=profile(Y)['H']);break
                        if hostile:break
            else:lookup[key]=A
            if hostile:break
        if hostile:break
    # A concrete order-sensitive pair, found deterministically.
    switch=None
    for A,B in product(bank,repeat=2):
        ha=profile(fuse(A,B))['H'];hb=profile(fuse(B,A))['H']
        if ha!=hb:switch=dict(A=A,B=B,H_AB=ha,H_BA=hb);break
    A=(50,44,57,17,34,8); B=(42,52,57,18,33,8)
    pa,pb=profile(A),profile(B)
    assert pa==dict(H=29,h=9,Q=84,k=0,U=52,V=41,d=4,e=5)
    assert pb==dict(H=29,h=9,Q=84,k=0,U=50,V=43,d=4,e=5)
    AB=mul(matrix(pa),matrix(pb));BA=mul(matrix(pb),matrix(pa))
    N=[[AB[i][j]-BA[i][j] for j in range(4)] for i in range(4)]
    assert N==[[0,150,0,150],[0,0,0,150],[0,0,0,0],[0,0,0,0]]
    context_controls=[]
    small=[(0,),(2,0),(0,1),(2,4,1)]
    # Literal small side contexts: largest composed graph 15 vertices.
    for L,R in product(small,repeat=2):
        pl,pr=profile(L),profile(R)
        X=fuse(fuse(fuse(L,A),B),R);Y=fuse(fuse(fuse(L,B),A),R)
        # Read H independently from the subset Hamilton DP, not rooted profiles.
        x,y=hp(X)[-1],hp(Y)[-1]
        delta=150*(pl['h']*pr['V']+pl['U']*pr['h']-pl['h']*pr['h'])
        assert x-y==delta and delta>0
        context_controls.append(dict(L=L,R=R,order=len(X),H_AB=x,H_BA=y,difference=delta))
    ladders=[]
    for a,b in [(1,1),(2,2),(3,3),(4,4)]:
        word=[0]*b+[1]*a # 0=B, 1=A; bubble A left.
        values=[];words=[]
        while True:
            M=matrix(profile((0,)))
            for bit in word:M=mul(M,matrix(pa if bit else pb))
            values.append(M[0][3]+M[0][0]);words.append(''.join('A' if bit else 'B' for bit in word))
            index=next((i for i in range(len(word)-1) if word[i:i+2]==[0,1]),None)
            if index is None:break
            word[index:index+2]=[1,0]
        assert len(values)==a*b+1 and all(x<y for x,y in zip(values,values[1:]))
        ladders.append(dict(a=a,b=b,order=1+5*(a+b),words=words,counts=values,method='exact proved matrix composition, not full graph DP'))
    result=dict(sorting_trade=dict(A=A,B=B,state_A=pa,state_B=pb,commutator=N,contexts=context_controls,ladders=ladders),status='FINITE-EXACT controls for new fusion law, not a tournament census',pair_controls=len(controls),largest_pair_order=5,associativity_controls=64,pairs=controls,towers=towers,parity_hostile=hostile,order_switch=switch)
    dest=Path(__file__).with_name('receipt.json');dest.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('pairs','towers','sorting_trade')},indent=2))
if __name__=='__main__':main()
