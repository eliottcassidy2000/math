from itertools import permutations
from math import comb, gcd

A = [(2,0),(4,0),(6,0),(7,0),(8,0),(7,1),(8,1),(5,2),(7,2),(6,3),(7,5),(8,5),(8,6)]
B = [(2,0),(4,0),(5,0),(6,0),(7,0),(5,1),(6,1),(7,2),(8,2),(5,3),(6,4),(8,5)]

def need(ok):
    if not ok:
        raise RuntimeError("certificate failed")

def rows(backs):
    r = [((1<<9)-1)^((2<<i)-1) for i in range(9)]
    for i,j in backs:
        r[i] ^= 1<<j
        r[j] ^= 1<<i
    return r

def warp(r):
    n=len(r)
    return [a|(1<<n) for a in r]+[1<<(n+1),(1<<n)-1]

def subset(r):
    n=len(r)
    dp=[[0]*n for _ in range(1<<n)]
    for j in range(n): dp[1<<j][j]=1
    for s in range(1,1<<n):
        for j,a in enumerate(dp[s]):
            if a:
                t=r[j]&~s
                while t:
                    b=t&-t
                    dp[s|b][b.bit_length()-1]+=a
                    t^=b
    return [1]+[sum(z) for z in dp[1:]]

def covers(P, slots):
    a=P
    for i in range(slots-1):
        b=[]
        for s in range(len(P)):
            z=s
            acc=0
            while True:
                acc+=P[z]*a[s^z]
                if not z: break
                z=(z-1)&s
            b.append(acc)
        a=b
    return a[-1]

def strong(r):
    for i in range(len(r)):
        seen=1<<i
        for _ in r:
            old=seen
            for j,a in enumerate(r):
                if (old>>j)&1: seen|=a
        if seen!=(1<<len(r))-1: return False
    return True

for backs, expected, upper in [(A,302037,782283),(B,302277,782763)]:
    r=rows(backs)
    need(strong(r))
    P=subset(r)
    goal=(2231,36088,expected)
    need(tuple(covers(P,m) for m in (1,2,3))==goal)
    bad=[0,0,0]
    for perm in permutations(range(9)):
        b=sum(not ((r[i]>>j)&1) for i,j in zip(perm,perm[1:]))
        if b<3: bad[b]+=1
    need(tuple(sum(bad[b]*comb(8+m-b,m-1-b) for b in range(m))
               for m in (1,2,3))==goal)
    w=warp(r)
    pp=subset(w)
    need(gcd(*goal[:2])==1 and pp[-1]==38319)
    need(covers(pp,2)==4*goal[1]+2*(goal[2]-goal[0]))
    need(subset(warp(w))[-1]==upper==-goal[0]+5*goal[1]+2*goal[2])
    print(goal,tuple(bad),pp[-1],upper)
