"""Independent literal and weighted contact replays for Mersenne-offset rails."""
from math import isqrt

def rows(word):
    n = len(word)
    spine = [(1 << n)-(1 << (i+1)) + ((w=="I") << n)
             for i,w in enumerate(word)]
    return spine + [sum((w=="O") << i for i,w in enumerate(word))]

def strong(rows):
    for start in range(len(rows)):
        seen=1<<start
        for _ in rows:
            for i,out in enumerate(rows):
                if seen>>i&1: seen |= out
        if seen != (1<<len(rows))-1: return False
    return True

def literal(word):
    R=rows(word); n=len(R)
    assert strong(R)
    dp=[[0]*n for _ in range(1<<n)]
    for i in range(n): dp[1<<i][i]=1
    for mask,row in enumerate(dp):
        for i,paths in enumerate(row):
            if not paths: continue
            nexts=R[i]&~mask
            while nexts:
                bit=nexts&-nexts; nexts-=bit
                dp[mask|bit][bit.bit_length()-1]+=paths
    return sum(dp[-1])

def weighted(word):
    h,u,previous=1,0,0
    for w in word:
        bit=int(w=="O")
        if not bit: h+=u+previous
        u=2*u+bit; previous=bit
    return h

def contact(a,b,c):
    return "O"+"I"*a+"O"+"I"*b+"O"+"I"*c

if __name__=="__main__":
    for a in range(1,10):
        for b in range(1,10):
            for c in range(1,10):
                x,y,z=(1<<a,1<<b,1<<c)
                H=4*x*y*z+(2*y+1)*(z-x)+1-y
                assert weighted(contact(a,b,c))==H
                if a==c and b%2==0:
                    assert (isqrt(H)**2==H)==(a==b//2-2)
    for q,w in [(4,"OIII"),(5,"OIIOOI"),(7,contact(1,3,1)),
                 (15,contact(5,3,5))]:
        H=literal(w)
        assert H==weighted(w)==(1<<q)-7 and isqrt(H)**2==H
        print(q,len(w)+1,H,w)
