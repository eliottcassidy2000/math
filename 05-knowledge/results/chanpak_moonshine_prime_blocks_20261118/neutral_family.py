"""Literal same-height clone/spine trade, with its degree marker."""
from pathlib import Path
from itertools import combinations
import sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'moonshine_prime_blocks'))
from replay import dp
from atlas import flip
from collections import Counter

def build(k,t,m,labels=None):
    if labels is None:
        S=tuple(range(m));C=tuple(range(m,m+k+t+1)); V=(m+k+t+1,m+k+t+2);z=m+k+t+3
    else:S,C,V,z=labels
    N=len(S)+len(C)+3;A=[0]*N
    def arc(i,j):A[i]|=1<<j
    r=C[k]
    for U in (S,C,V):
        for i,j in combinations(U,2):arc(i,j)
    for i in S:
        for j in (*C,*V,z):
            if j==r and i in S[:2]:arc(j,i)
            else:arc(i,j)
    for i in C:
        for v in V:arc(i,v)
        arc(z,i)
    for v in V:arc(v,z)
    return tuple(A),(S,C,V,z)

def reads(A,L):return sum(dp(A)[-1])

def degree_marker(T):
    N=len(T);ds={r.bit_count() for r in T}
    assert max(ds)==N-2
    gap=next(j for j in range(N-2,-1,-1) if j not in ds)
    return N-gap-1

def run():
    for j,w in ((4,3),(4,4),(4,5),(5,5),(5,6),(6,6)):
        k=j-2
        for t in range(1,w-1):
            m=w-t+1
            T,(S,C,V,z)=build(k,t,m)
            h=3**j*2**w+1
            assert reads(T,None)==h and degree_marker(T)==m
            oldr=C[k]
            S2=V[::-1]+C[k+1:][::-1]
            C2=C[:k][::-1]+(z,)+S[2:][::-1]
            X,_=build(k,m-2,t+2,(S2,C2,(S[1],S[0]),oldr))
            assert X==tuple(((1<<len(T))-1)^row^(1<<a) for a,row in enumerate(T))
            if t>1:
                x=C[k+1];Q,_=build(k,t-1,m+1,(S+(x,),C[:k+1]+C[k+2:],V,z))
                arcs={(i,l) for i,l in combinations(range(len(T)),2) if (T[i]>>l&1)!=(Q[i]>>l&1)}
                assert arcs=={tuple(sorted((x,u))) for u in (*C[:k+1],z)}
                assert reads(Q,None)==h and degree_marker(Q)==m+1
                r=C[k]
                A=flip(T,r,x)
                f1=(r,*S,z,*C[:k],*C[k+1:],*V)
                assert all(T[p]>>q&1 and A[p]>>q&1 for p,q in zip(f1,f1[1:]))
                assert reads(A,None)==3*(h-1)//2+1
                # Reidentify x as earliest predecessor.
                order=(x,)+C[:k]
                A1,_=build(k+1,t-1,m,(S,order+(r,)+C[k+2:],V,z))
                f2=(r,*S,x,*C[:k],*C[k+2:],*V,z)
                assert all(A1[p]>>q&1 and Q[p]>>q&1 for p,q in zip(f2,f2[1:]))
                assert flip(A1,z,x)==Q
                y=C[0]
                low=flip(T,z,y)
                both=flip(A,z,y)
                common=(r,*S,*C[:k],*C[k+1:],*V,z)
                assert all(all(W[a]>>b&1 for a,b in zip(common,common[1:]))
                           for W in (T,A,low,both))
                other,_=build(k,t-1,m+1,(S+(y,),C[1:k]+(x,r)+C[k+2:],V,z))
                assert reads(low,None)==2*(h-1)//3+1 and both==other
                assert reads(both,None)==h
                before=Counter(v.bit_count() for v in T)
                after=Counter(v.bit_count() for v in Q)
                assert sum((before-after).values())>2
        print(j,w,'height',h,'nonisomorphic vertices',w-2,'local edit bundle',j)
if __name__=='__main__':run()
