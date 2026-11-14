from itertools import product
from random import Random

def count(R):
    n = len(R)
    dp = [[0]*n for _ in range(1<<n)]
    for j in range(n): dp[1<<j][j] = 1
    for s in range(1,1<<n):
        for j in range(n):
            if not dp[s][j]: continue
            z = R[j] & ~s
            while z:
                b = z & -z
                dp[s|b][b.bit_length()-1] += dp[s][j]
                z ^= b
    return sum(dp[-1])

def paths(word):
    R = [sum(1<<j for j in range(i+1,len(word))) | ((1-int(b))<<len(word))
         for i,b in enumerate(word)]
    R.append(sum(int(b)<<i for i,b in enumerate(word)))
    return count(R)

def h(x):
    v=1
    while x:
        x,r=divmod(x,2)
        if not r: v+=x+(x%2)
    return v

def rc(w): return ''.join('10'[int(b)] for b in reversed(w))
def counts(x,y,b):
    return h(x)+h(y)-1+2*x*y+([x,y][b]+[x,y][b]%2)

def decode(N,d,b):
    for e,f in product(range(2),repeat=2):
        k=2*d+e-f
        lo=max(0,-k)
        lo+=(e-lo)%2
        hi=2*abs(d)+N+4
        hi+=(e-hi)%2
        v=(hi-lo)//2
        u=0
        while v>u:
            m=(v+u)//2
            a=lo+2*m
            if counts(a,a+k,b)>=N: v=m
            else: u=m+1
        x=lo+2*u
        if counts(x,x+k,b)==N:
            return (x,x+k)
    return None

rng = Random(1729)
checks=0
for n in range(1,9):
    for w in map(''.join,product('01',repeat=n)):
        got=paths(w)
        for i in range(len(w)):
            a,b,c=w[:i],int(w[i]),w[i+1:]
            x,y=int(a or '0',2),int(rc(c) or '0',2)
            assert counts(x,y,b)==got
            d=(y+1)//2-(x+1)//2
            a2,c2=decode(got,d,b)
            assert counts(a2,c2,b)==got
            checks+=1
for _ in range(250):
    x,y=rng.randrange(1<<50),rng.randrange(1<<50)
    expected=(y+y%2)-(x+x%2)
    assert counts(x,y,1)-counts(x,y,0)==expected
    for b in (0,1):
        assert decode(counts(x,y,b),expected//2,b)
assert paths('1100010')==paths('1101010')==73
print(f'PASS: {checks} literal pivot decompositions, 250 large target inversions.')
