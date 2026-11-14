#!/usr/bin/env python3
"""Construct literal strong tournaments at the two changed-core exact windows."""
from replay import h, dual, blocks, join
from math import isqrt
import sys
# (offset, wing factors a,b, contacts U, outgoing rows of the core)
AFFINE=(
 (885,10,1,115,(2,36,25,19,35,77,31)),
 (887,10,1,103,(2,36,41,19,7,89,31)),
 (889,5,6,31,(2,16,3,7,45,15,63)),
 (891,10,1,100,(6,12,24,49,35,7,63)),
 (893,7,1,117,(6,24,42,49,69,19,47)),
 (895,5,3,27,(2,12,17,37,11,87,31)),
 (897,5,3,15,(2,36,9,19,7,93,31)),
 (899,5,2,64,(16,5,9,19,102,15,47)),
 (901,7,1,93,(6,24,42,81,37,11,55)),
 (903,7,3,53,(0,5,81,7,43,15,59)),
 (905,10,6,30,(0,5,25,19,35,15)),
)
QUADRATIC=(
 (-11,3,3,21,(0,1,3,7,15)),
 (-9,1,1,33,(2,4,1,7,47,79,31)),
 (-7,3,3,17,(0,1,3,7,15)),
 (-5,6,6,7,(0,1,3)),
 (-3,9,1,0,(0,1,3,7,15)),
 (-1,3,3,9,(2,4,9,3)),
 (1,3,3,6,(2,4,9,3)),
 (3,1,1,34,(2,32,11,19,7,29)),
 (5,3,3,0,(2,4,9,3)),
)
def g(x):
    z=1
    while x:
        x,d=divmod(x,2)
        if not d:z+=x+x%2
    return z

def trans(n):
    return tuple(((1<<n)-1)^((1<<(i+1))-1) for i in range(n))

def attach(r,U,x,y):
    left=bin(x)[2:]
    right=''.join('1' if v=='0' else '0' for v in bin(y)[2:][::-1])
    n=len(left)+len(r)+len(right)
    p=sum(int(v)<<i for i,v in enumerate(left))+ (U<<len(left))
    p+=sum(int(v)<<(i+len(left)+len(r)) for i,v in enumerate(right))
    R=join((trans(len(left)),r,trans(len(right))))
    return tuple(a | ((1-(p>>i&1))<<n) for i,a in enumerate(R))+(p,)

def table(r):
    n=len(r)
    dp=[[0]*n for _ in range(1<<n)]
    for i in range(n):dp[1<<i][i]=1
    for m,v in enumerate(dp):
        for i,s in enumerate(v):
            if s:
                for j in range(n):
                    if r[i]>>j&1 and not m>>j&1:dp[m|1<<j][j]+=s
    return dp

def profile(r,U):
    n=len(r); F=(1<<n)-1
    e=table(r); s=table(dual(r)); H=[1]+list(map(sum,e[1:]))
    P=Q=E=0
    for m in range(F+1):
        P+=sum(s[m][i] for i in range(n) if U>>i&1)*H[F^m]
        Q+=sum(e[m][i] for i in range(n) if not (U>>i&1))*H[F^m]
        E+=sum(e[m][i] for i in range(n) if U>>i&1)*sum(s[F^m][i] for i in range(n) if U>>i&1)
    S=sum(e[F][i] for i in range(n) if U>>i&1)
    return H[F],sum(v*H[F^m] for m,v in enumerate(H)),P,Q,S,E

def case(mode,z,t):
    d,a,b,U,r=z
    u=1<<t
    return attach(r,U,a*u,b if mode=='affine' else b*u-1)

def select(N):
    for mode,C,bank in [('affine',13790,AFFINE),('quadratic',288,QUADRATIC)]:
        for z in bank:
            d=z[0]
            if (N-d)%C:continue
            q=(N-d)//C
            if mode=='quadratic':
                if q<4:continue
                u=isqrt(q)
                if u*u!=q:continue
            else:
                u=q
            if u>=2 and not (u&(u-1)):
                yield mode,u.bit_length()-1,d,case(mode,z,u.bit_length()-1)

def main():
    for mode,bank in [('quadratic',QUADRATIC),('affine',AFFINE)]:
        for z in bank:
            d,a,b,U,r=z
            assert all(not (v>>i&1) for i,v in enumerate(r))
            assert all((v>>j&1)+(r[j]>>i&1)==1 for i,v in enumerate(r) for j in range(i+1,len(r)))
            hc,f,p,q,s,e=profile(r,U)
            if mode=='quadratic':
                assert (f*a*b,b*p+a*(q-f+hc))==(288,0)
                assert d==hc-e+hc*(g(a)+g(b-1)-a+a%2-2)
                C=288
            else:
                assert a*(f*b+q+hc)==13790
                assert d==hc+p*(b+1)+s*(b%2-1)-e+hc*(g(a)+g(b)-a+a%2-2)
                C=13790
            print(mode,d,hc,f,p,q,s,e)
            for t in (1,2):
                N=C*(1<<(t*(2 if mode=='quadratic' else 1)))+d
                R=case(mode,z,t)
                assert len(blocks(R))==1 and h(R)==N
                assert (mode,t,d,R) in list(select(N))
    print('PASS: all core profiles and 40 literal tournament targets across four windows')

if __name__=='__main__':
    if len(sys.argv)==1: main()
    else:
        ans=list(select(int(sys.argv[1])))
        for t in ans:print(*t)
        if not ans: print('Not in these two certified families.')
