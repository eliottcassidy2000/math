"""Exact finite controls for the matching-capacity and section theorems."""
from functools import lru_cache
from itertools import permutations
from math import comb, factorial
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent /
                      "session328_fermat_gluing_and_phase_repair_20261114"))
from profile_tools import extend, hp_all, strong

def ctx_coeff(R, k):
    """Number of literal quotient walks; module index is last."""
    b = len(R)-1
    @lru_cache(None)
    def rec(last, mask, rem):
        if not mask and not rem:
            return 1
        total = rec(b, mask, rem-1) if rem and (last<0 or R[last]>>b&1) else 0
        for j in range(b):
            if mask>>j&1 and (last<0 or R[last]>>j&1):
                total += rec(j, mask^(1<<j), rem)
        return total
    return rec(-1, (1<<b)-1, k)

def matching(O, I, R):
    if not O:
        return 0
    x,*rest=O
    return max([matching(rest,I,R)]+
               [1+matching(rest,[a for a in I if a!=y],R) for y in I if R[x]>>y&1])

def labelled(n):
    out=[()]
    for i in range(n):
        out=[extend(A,c) for A in out for c in range(1<<i)]
    return out

def trans(n):
    return tuple((1<<n)-(1<<(j+1)) for j in range(n))

def ordinal(U,W):
    return tuple(x|(((1<<len(W))-1)<<len(U)) for x in U)+tuple(x<<len(U) for x in W)

def cyclic(U,W,X):
    n,b,c=map(len,(U,W,X))
    return (tuple(x|(((1<<b)-1)<<n) for x in U)+
            tuple((y<<n)|(((1<<c)-1)<<(n+b)) for y in W)+
            tuple((z<<(n+b))|((1<<n)-1) for z in X))

def replace(R,r,U):
    """Different implementation: build arrows on each pair of literal vertices."""
    tags=[i for i in range(len(R)) if i!=r]+[r]*len(U)
    G=[0]*len(tags)
    for i,t in enumerate(tags):
        for j,s in enumerate(tags):
            if t!=s:
                a=R[t]>>s&1
            elif t==r:
                a=U[i-len(R)+1]>>(j-len(R)+1)&1
            else:
                a=0
            G[i]|=a<<j
    return tuple(G)

def coefficients(A,r):
    H=hp_all(A)
    Z=H[:]
    Z[0]=0
    ans=[Z[-1]]
    for i in range(2,r+1):
        N=[0]*len(H)
        for s in range(1,len(H)):
            a=(s-1)&s
            while a:
                N[s]+=H[a]*Z[s^a]
                a=(a-1)&s
        Z=N
        ans.append(Z[-1])
    return ans

def triangle():
    return (2,4,1)

def skeleton(i,j,k):
    return 3 if i==j==k else int(max(i,j,k)-min(i,j,k)==1)

def coeff_cyclic(W,X,r):
    a=coefficients(W,len(W)); b=coefficients(X,len(X))
    return [sum(e*f*skeleton(i,j,k) for j,e in enumerate(a,1)
                for k,f in enumerate(b,1)) for i in range(1,r+1)]

def main():
    for b in range(5):
        num=0
        for R in labelled(b+1):
            O=[i for i in range(b) if R[b]>>i&1]
            I=[i for i in range(b) if not R[b]>>i&1]
            nu=matching(O,I,R)
            for k in range(1,b//2+3):
                co=ctx_coeff(R,k)
                assert bool(co)==(k<=nu+1) and co%factorial(k-1)==0
            if len(O)==len(I):
                p=sum(all(R[x]>>y&1 for x,y in zip(O,p)) for p in permutations(I))
                assert ctx_coeff(R,1+b//2)==factorial(b//2)*p
            num+=1
        print('context census', b, num)

    R=(2,4,9,3)
    controlA=(128,129,99,71,15,27,51,124)
    controlB=(36,97,50,7,11,216,29,95)
    for U in (controlA,controlB):
        h,p,*_=coefficients(U,3)
        S=[hp_all(replace(R,r,U))[-1] for r in (0,1)]
        assert S==[5*h+2*p,5*h+3*p]
        print('two-coordinate small section', S)

    for U,r in [(controlA,3),(controlB,3),(extend(controlA,0),4)]:
        b=2*r-2
        W=trans(1)
        cases=[(W,trans(b-1)),(W,ordinal(triangle(),trans(b-4)))]+[
            (trans(m),trans(b-m)) for m in range(2,r)]
        K=[coeff_cyclic(V,X,r) for V,X in cases]
        P=coefficients(U,r)
        Ys=[]
        for row,(V,X) in zip(K,cases):
            G=cyclic(U,V,X)
            assert strong(G)
            val=hp_all(G)[-1]
            assert sum(a*b for a,b in zip(row,P))==val
            Ys.append(val)
        a,b,*_=K[0]; c,d,*_=K[1]
        det=a*d-b*c
        h=(d*Ys[0]-b*Ys[1]);p=(a*Ys[1]-c*Ys[0])
        assert h%det==0 and p%det==0
        rec=[h//det,p//det]
        for row,v in zip(K[2:],Ys[2:]):
            x=v-sum(a*b for a,b in zip(row,rec));y=row[len(rec)]
            assert y>0 and x%y==0
            rec.append(x//y)
        assert rec==P
        print('triangle section', len(U), r, 'P',P, 'outputs',Ys)
        if U==controlA: score=Ys
        if U==controlB: assert Ys[:2]==score[:2] and Ys[2]-score[2]==2880
    R=list(cyclic((0,),trans(2),trans(2)))
    R[1]^=1<<4
    R[4]^=1<<1
    T=[replace(tuple(R),0,U) for U in (controlA,controlB)]
    assert all(strong(G) for G in T)
    N=[hp_all(G)[-1] for G in T]
    assert N[1]-N[0]==1440
    print('one missing crossing', N)
    for d in range(1,8):
        PP=list(permutations(range(d)))
        for a in range(1<<(d-1)):
            M=[[(i<=j) or (i and j==0 and (a>>(i-1))&1) for j in range(d)] for i in range(d)]
            perm=sum(all(M[i][j] for i,j in enumerate(p)) for p in PP)
            assert perm==a+1
    print('binary frontier permanents d=1..7: PASS')
    print('matching sections: PASS')


if __name__ == "__main__":
    main()
