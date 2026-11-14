from functools import cache
from math import isqrt

def rev(n,b):return int(format(n,f'0{b}b')[::-1],2) if b else 0
def E(n):
    e=0
    while n:
        t=n//2
        if n%2:e+=(t+1)//2
        n=t
    return e
def H(n):return 1+n-n%2-2*E(n)
def prime(n):return n>=2 and all(n%i for i in range(2,isqrt(n)+1))
def run_dp(word):
    d={(-1,False):1}
    for b in map(int,word):
        q={}
        for (left,right),c in d.items():
            k=(b,right)
            q[k]=q.get(k,0)+c
            if right or b==1:
                k=(left,True)
                q[k]=q.get(k,0)+c
        d=q
    return sum(c for (left,right),c in d.items() if left!=1)

def pair(A,b,d=1):
    B=rev((1<<b)-A-d-2,b)
    n=(8*A+5)*(1<<b)+B
    N=H(n)
    return N,N+2*d,n,B

def inverse(N,d=1):
    b=1
    while (1<<b)**2<=N:
        l=1<<(b-1)
        r=(1<<b)-d-2
        while l<=r:
            m=(l+r)//2
            h=pair(m,b,d)[0]
            if h==N:yield (m,b)
            if h<N:l=m+1
            else:r=m-1
        b+=1

def apex_rows(n,L):
    t=[]
    for i in range(L):
        row=((1<<L)-1)^((1<<(i+1))-1)
        if not(n>>(L-1-i)&1):row|=1<<L
        t.append(row)
    t.append(rev(n,L))
    return t

def table(T):
    n=len(T)
    @cache
    def start(S,v):
        s=S^(1<<v)
        if not s:return 1
        return sum(start(s,j) for j in range(n) if s>>j&1 and T[v]>>j&1)
    Hs=[1]+[sum(start(S,v) for v in range(n) if S>>v&1) for S in range(1,1<<n)]
    return Hs,start

def rooted_profile(T,v):
    Hs,start=table(T)
    m=len(Hs)-1
    u=w=0
    for S in range(1,m+1):
        if S>>v&1:
            c=start(S,v)*Hs[m^S]
            u+=c
            w+=(-1)**S.bit_count()*c
    return Hs[-1],u,w

def ordered(A,B):
    p,q=len(A),len(B)
    suffix=((1<<q)-1)<<p
    return tuple(r|suffix for r in A)+tuple(r<<p for r in B)

def end_subset(B,U):
    T=tuple(sum(1<<j for j in range(len(B)) if B[j]>>i&1) for i in range(len(B)))
    return sum(rooted_profile(T,u)[1] for u in U)

def bridge(T,U,r):
    n=len(T)
    out=list(ordered(T,tuple(((1<<r)-1)^((1<<(i+1))-1) for i in range(r))))
    w=len(out)
    out=[a|((1<<w) if i not in U else 0) for i,a in enumerate(out)]
    out.append(sum(1<<i for i in U))
    return out

def F2(B):
    t,_=table(B)
    return sum(a*b for a,b in zip(t,t[::-1]))

if __name__=='__main__':
    for n in range(1,256):
        assert H(n)==run_dp(bin(n)[2:])==table(apex_rows(n,n.bit_length()))[0][-1]
    for d in (1,2,5,11,13,26):
        for b in range(3,10):
            last=-1
            for A in range(1<<(b-1),(1<<b)-d-1):
                N,Q,n,B=pair(A,b,d)
                assert H(n+(1<<(b+1)))==Q
                assert N>last
                assert (A,b) in tuple(inverse(N,d))
                last=N
    for b,A,N in ((9,434,1517507),(10,646,3354407)):
        L,U,n,B=pair(A,b)
        assert L==N
        assert run_dp(bin(n)[2:])==N and run_dp(bin(n+(1<<(b+1)))[2:])==N+2
        K=N if N%7==0 else N+2
        J=N+2 if K==N else N
        assert K%7==0 and prime(K//7) and prime(J)
    mid=(14,12,8,0)
    cores=((),(0,),(2,0),(2,4,1))
    for L in cores:
        for B in cores:
            T=ordered(ordered(L,mid),B)
            for UL in ((),tuple(range(len(L))),tuple(range(0,len(L),2))):
                RL=sum(rooted_profile(L,j)[1] for j in UL)
                for UB in ((),tuple(range(len(B))),tuple(range(0,len(B),2))):
                    QB=end_subset(B,UB)
                    fixed=UL+tuple(len(L)+4+j for j in UB)
                    P=fixed+tuple(len(L)+j for j in (0,2,3))
                    Q=fixed+tuple(len(L)+j for j in (0,1,3))
                    hl,hb=table(L)[0][-1],table(B)[0][-1]
                    for r in (0,1):
                        assert table(bridge(T,Q,r))[0][-1]-table(bridge(T,P,r))[0][-1]==2*(hl*((F2(B)<<r)-QB)+hb*RL)
    P=(24,1,3,6,14);Q=(18,8,3,5,14)
    assert rooted_profile(P,0)==(11,13,-1) and rooted_profile(Q,4)==(13,26,-2)
    for r in range(1,5):
        assert table(bridge(P,(0,),r+1))[0][-1]==11+(26<<r)
        assert table(bridge(Q,(4,),r))[0][-1]==13+(26<<r)
    print("balanced contact words: 2d switches, quadratic exact inverse, prime/7p pairs, 11/13 rigidity, and 11/13/26 literal bridge passed")
