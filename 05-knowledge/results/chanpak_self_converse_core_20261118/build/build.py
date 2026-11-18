"""Literal self-converse hosts and exact 7-adic dimension compiler."""
from math import isqrt
import argparse,json

D4=(2,4,9,3)
def TT(n):return tuple((1<<n)-(1<<(i+1)) for i in range(n))
def transpose(a):return tuple(sum((x>>i&1)<<j for j,x in enumerate(a)) for i in range(len(a)))
def ordinal(A,B):return tuple(r|(((1<<len(B))-1)<<len(A)) for r in A)+tuple(r<<len(A) for r in B)
def cyclic(A,B,C):
    n,m,k=len(A),len(B),len(C)
    return tuple(r|(((1<<m)-1)<<n) for r in A)+tuple((r<<n)|(((1<<k)-1)<<(n+m)) for r in B)+tuple((r<<(n+m))|((1<<n)-1) for r in C)
def host(T,j,m,t):
    X=T
    for _ in range(j):X=ordinal(X,D4)
    X=ordinal(X,TT(m))
    return cyclic(X,transpose(X),TT(t))
def htable(a):
    n=len(a);P=[[0]*n for _ in range(1<<n)]
    for i in range(n):P[1<<i][i]=1
    for s in range(1,len(P)):
        for i,z in enumerate(P[s]):
            b=a[i]&~s
            while z and b:
                v=b&-b;b-=v;P[s|v][v.bit_length()-1]+=z
    return [1]+[sum(q) for q in P[1:]]
def profile(a):
    h=htable(a);M=len(h)-1;P2=[0]*(M+1)
    for s in range(M+1):
        t=s
        while True:
            P2[s]+=h[t]*h[s^t]
            if not t:break
            t=(t-1)&s
    return h[M],P2[M],sum(h[s]*P2[s^M] for s in range(M+1))
def F(h,a,b):
    r=a-2*h;s=b-3*a+3*h
    return 5*h*h+6*h*r+7*r*r+4*r*s+2*s*s
def H1(h,a):return (a-h)**2+2*h*h
def phase(h,a,b):
    if b%7==0 or h%7==a%7==0:raise ValueError('source is outside exact phase class')
    for j in range(6):
        H,A,B=(h*pow(5,j,7))%7,(a*pow(32,j,7))%7,(b*pow(117,j,7))%7
        for m in range(6):
            C,D=A*pow(2,m,7)%7,B*pow(3,m,7)%7
            if H and C and (C,D)==(2*H%7,4*H%7):return j,m
            if H and not C and D==H:return j,m
            if not H and C==D:return j,m
    raise AssertionError('phase not found')
def evaluate(h,a,b,j,m,mod):
    return F(h*pow(5,j,mod)%mod,a*pow(32,j,mod)*pow(2,m,mod)%mod,b*pow(117,j,mod)*pow(3,m,mod)%mod)%mod
def lift(h,a,b,target,k):
    if target%7 or k<1:raise ValueError('choose a multiple of seven and positive precision')
    j,m=phase(h,a,b);mod=7
    for _ in range(1,k):
        candidates=[m+6*(mod//7)*q for q in range(7)];mod*=7
        m,=[z for z in candidates if evaluate(h,a,b,j,z,mod)==target%mod]
    return j,m
def check_sc(a,half,t):
    assert len(a)==2*half+t
    s=[i+half if i<half else i-half if i<2*half else 4*half+t-1-i for i in range(len(a))]
    assert all(a[i]>>j&1==a[s[j]]>>s[i]&1 for i in range(len(a)) for j in range(len(a)))
def test():
    assert profile(D4)==(5,32,117)
    for T in (TT(1),TT(2),TT(3),(2,4,1),D4,(8,5,9,50,7,23,63)):
        h,a,b=profile(T)
        Y=[]
        for m in range(3):
            U=host(T,0,m,1);check_sc(U,len(T)+m,1)
            v=htable(U)[-1];assert v==H1(h,2**m*a);Y.append(v)
        assert Y[2]-6*Y[1]+8*Y[0]==9*h*h
        assert Y[2]-3*Y[1]+2*Y[0]==6*a*a
        U=host(T,0,0,2);check_sc(U,len(T),2)
        assert htable(U)[-1]==F(h,a,b)
    for X in (TT(1),TT(2),TT(3),(2,4,1),D4):
        for Z,v in ((TT(1),3),(TT(2),5),(TT(3),1),((2,4,1),7)):
            assert htable(cyclic(X,transpose(X),Z))[-1]%8==v
    T=(8,5,9,50,7,23,63);H=F(*profile(T));assert H==37050461
    p=H//7;assert p==5292923 and all(p%i for i in range(2,isqrt(p)+1))
    print('literal prime quotient host',H,p,*host(T,0,0,2))
    from itertools import product
    tested=0
    for h,a,b in product(range(7),repeat=3):
        if b and (h or a):
            for k in (2,3,6):
                target=(h*6000+40*a+5*b)*7
                j,m=lift(h,a,b,target,k)
                assert evaluate(h,a,b,j,m,7**k)==target%(7**k)
                assert 0<=j<6 and 0<=m<6*7**(k-1)
                tested+=1
        else:
            for j,m in product(range(6),repeat=2):
                r=evaluate(h,a,b,j,m,7)
                assert r!=0 or h==a==b==0
    print('PASS literal mirror counts and anti-automorphism, all 20 mod-eight controls,',tested,'phase lifts plus excluded profiles')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--profile',nargs=3,type=int);p.add_argument('--target',type=int,default=7);p.add_argument('-k',type=int,default=2);a=p.parse_args()
    if a.profile:print(json.dumps(lift(*a.profile,a.target,a.k)))
    else:test()
