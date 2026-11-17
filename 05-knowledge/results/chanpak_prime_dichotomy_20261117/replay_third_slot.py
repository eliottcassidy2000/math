from slot_profiles import profile, end_profiles
from math import isqrt

def tt(k):
    return tuple((1<<k)-(2<<i) for i in range(k))

def fixture(B):
    size=len(B)
    groups=[(0,),(1,),(2,),tuple(range(3,3+size)),
            (3+size,),(4+size,),(5+size,)]
    n=size+7
    a=[0]*n
    for k,G in enumerate(groups):
        for i in G:
            for J in groups[k+1:]:
                a[i]|=sum(1<<j for j in J)
            if k<4:a[-1]|=1<<i
            else:a[i]|=1<<(n-1)
    for i,r in enumerate(B):a[i+3]|=r<<3
    return a

def prime(n):return n>=2 and all(n%d for d in range(2,isqrt(n)+1))

for B in (tt(1),tt(2),tt(3),(2,4,1),(118,88,26,97,104,6,36)):
    order=len(B); f=(1<<order)-1
    counts,_=end_profiles(B)
    h=counts[f]
    P=sum(counts[m]*counts[f^m] for m in range(f+1))
    P3=0
    for m in range(f+1):
        c=f^m
        s=c
        while True:
            P3+=counts[m]*counts[s]*counts[c^s]
            if not s:break
            s=(s-1)&c
    F=fixture(B)
    g,_,W,_,E,_,K=profile(F,5)
    U=g-K
    p=3*W+2*g
    if U+E!=g or (g,W,p)!=(56*P-5*h,494*P3-80*P-7*h,1482*P3-128*P-31*h):
        raise ValueError((B,g,W,U,E,p))
    print(B,(h,P,P3),(g,W,p),prime(p))
