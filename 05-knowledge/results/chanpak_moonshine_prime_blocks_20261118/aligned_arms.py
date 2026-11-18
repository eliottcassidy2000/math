"""Exact positional decoder for aligned shared-root arms."""
from fractions import Fraction
from math import isqrt
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'moonshine_prime_blocks'))
from clone_replay import TT
from frame_replay import add
from replay import profile,dp,transpose

A=(26,12,9,16,6)
B=(6,12,8,1)

def relabel(T,J):
    return tuple(sum((T[i]>>j&1)<<k for k,j in enumerate(J)) for i in J)
def arm(T,r):
    T=relabel(T,tuple(i for i in range(len(T)) if i!=r)+(r,))
    n=len(T)-1; F=(1<<n)-1
    Q=tuple(a&F for a in T[:-1]); E=dp(Q); S=dp(transpose(Q))
    g,B=profile(Q); H=[1]+[sum(e) for e in E[1:]]
    U,V,ub,vb=0,0,0,0
    for s in range(1,F+1):
        O=sum(x for i,x in enumerate(S[s]) if T[-1]>>i&1)*H[F^s]
        I=sum(x for i,x in enumerate(E[s]) if T[i]>>n&1)*H[F^s]
        sig=(-1)**s.bit_count()
        U+=O;V+=sig*O;ub+=I;vb+=sig*I
    return sum(dp(T)[-1]),g,B,U,V,ub,vb

def glue(A,B):
    # Both roots last. Reuse them as a single last vertex.
    n,m=len(A)-1,len(B)-1;F=(1<<n)-1;G=(1<<m)-1
    return (tuple((r&F)|(G<<n)|(((r>>n)&1)<<(n+m)) for r in A[:-1])
           +tuple(((r&G)<<n)|(((r>>m)&1)<<(n+m)) for r in B[:-1])
           +((A[-1]&F)|((B[-1]&G)<<n),))

def upper(k,I):
    return 1+sum(1<<max(1,j-i-1) for i in range(k) for j in range(i+1,k)
                 if I>>i&1 and not I>>j&1)

def decode(p,m):
    center=(576*8**m+84*m+61)//49
    a,r=divmod(7*(center-p),12)
    if r or a<0:return None
    q=a+1
    j=(q.bit_length()-1)//3
    return j if 8**j==q and j<=m else None

def run():
    assert arm(A,4)==(13,3,24,12,0,15,-1)
    assert arm(B,3)==(5,1,8,4,0,6,0)
    assert arm(glue(A,B),7)[0]==97 and arm(glue(B,A),7)[0]==85
    rows=((48157,0,20,402849),(48049,2,19,15509))
    for p,j,k,I in rows:
        assert all(p%d for d in range(2,isqrt(p)+1))
        assert decode(p,4)==j and upper(k,I)==7*p
        W=[B]*j+[A]+[B]*(4-j)
        T=W[0]
        for X in W[1:]:T=glue(T,X)
        assert sum(dp(T)[-1])==p
        print('aligned ordinary word then different sevenfold graph',p,7*p)
if __name__=='__main__':run()
