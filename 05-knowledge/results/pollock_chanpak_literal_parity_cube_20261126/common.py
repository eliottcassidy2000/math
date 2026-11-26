from itertools import product
from math import isqrt
def F(x): return (x**3-x)//6
def inverse(x):
    lo,hi=0,4
    while F(hi)<=x: hi*=2
    while hi-lo>1:
        m=(lo+hi)//2
        if F(m)<=x: lo=m
        else: hi=m
    return lo
def E(t,a,b,d):return F(t)+F(t+a+b)+F(t+a+d)+F(t+b+d)
kmax=26
for k in range(1,kmax+1):
    N=(38*3**k+22*2**k-18)//2
    pairs=[]; n=0
    for i in range(k):
        for j in range(i,k-i):
            if j>k-1-i-j:continue
            a,b,d=2**i,2**j,2**(k-1-i-j)
            t=1
            while E(t,a,b,d)<=N:
                n+=1
                rem=N-E(t,a,b,d); z=inverse(rem)
                if F(z)==rem:pairs.append((t,a,b,d,z))
                t+=1
    print(k,N,n,pairs,flush=True)
