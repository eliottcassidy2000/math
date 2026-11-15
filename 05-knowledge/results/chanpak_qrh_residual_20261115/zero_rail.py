from math import isqrt
from complement import value,baseline,first_across,K,transfer

def zero_rail(N):
    s=isqrt(N//108)
    s+=(108*s*s<N)
    assert N>=10000
    score,g,S=min((1+isqrt((s+(-s)%g)*(1+(-s)%g))+g+(s+(-s)%g+g-1)//g,g,s+(-s)%g)
                 for g in (1<<j for j in range(1,s.bit_length())) if g*g<=s)
    x0=(S+1)//2
    x0+=(not x0&1)
    def xy(k):
        x=x0+2*k
        return x,S-x
    def V(k):
        return value(*xy(k))
    def low(k):
        return baseline(*xy(k))+45*S+45+K[-1]<=N
    cap=(S-1-x0)//2
    assert V(0)[0]>N and low(cap)
    _,end=first_across(0,cap,low)
    i,_=first_across(0,end,lambda k:V(k)[0]<=N)
    _,j=first_across(0,end,lambda k:V(k)[-1]<=N)
    lo,hi=(xy(j),V(j)),(xy(i),V(i))
    assert lo[1][-1]<=N<hi[1][0]
    assert max(N-lo[1][0],hi[1][-1]-N)<=4000*score
    return score,g,S,lo,hi

if __name__=='__main__':
    count=0
    for s in range(1,11):
        g=1<<s
        for b in range(1,g,2):
            assert transfer(g,b)[0]+transfer(g,g-b)[0]==g//2-1
            count+=1
    print('low rail cancellations',count)
    for N in range(20_000,50_000):
        zero_rail(N)
    for r in range(4,260):
        S=1<<(2*r)
        N=108*S*S-1
        score,g,T,*_=zero_rail(N)
        assert score<=10*isqrt(S)
    print('validated 30,000 consecutive targets and 256 aligned heights')
