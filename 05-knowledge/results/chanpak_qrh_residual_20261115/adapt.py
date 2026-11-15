#!/usr/bin/env python3
from math import isqrt
from complement import value,baseline,first_across,K

def adaptive(N):
    assert N>=10000
    # s=ceil sqrt(N/108), computed without floating point.
    s=isqrt(N//108)
    s+=(108*s*s<N)
    z=[]
    g=1
    while g*g<=s:
        delta=(-1-s)%g
        S=s+delta
        score=1+isqrt(S*(delta+1))+g+(S+g-1)//g
        z.append((score,g,S))
        g*=2
    score,g,S=min(z)
    def xy(k):
        x=(S+1)//2+k
        return x,S-x
    def V(k):
        return value(*xy(k))
    def low(k):
        return baseline(*xy(k))+45*S+45+K[-1]<=N
    cap=S-(S+1)//2-1
    assert V(0)[0]>N and low(cap)
    _,end=first_across(0,cap,low)
    i,_=first_across(0,end,lambda k:V(k)[0]<=N)
    _,j=first_across(0,end,lambda k:V(k)[-1]<=N)
    lo,hi=(xy(j),V(j)),(xy(i),V(i))
    assert lo[1][-1]<=N<hi[1][0]
    # Test the quantitative bound with a conservative absolute constant.
    assert max(N-lo[1][0],hi[1][-1]-N)<=2000*score
    return (score,g,S,s,lo,hi)

if __name__=='__main__':
    ns=0
    for r in range(3,260):
        S=(1<<(2*r))-1
        for j in (1, S-1, 108*S):
            N=108*S*S-j
            score,g,T,s,lo,hi=adaptive(N)
            assert s==S and score<=10*isqrt(S)
            ns+=1
    print('adaptive aligned targets verified:',ns)
    for N in (35099,10**6+1,10**12+1,10**36+1):
        print(N,adaptive(N))
