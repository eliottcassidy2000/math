from collections import Counter
from itertools import product, combinations
from math import gcd, factorial
from functools import reduce
def dp(rows):
    n=len(rows);d=[[0]*n for _ in range(1<<n)]
    assert all(not m>>i&1 for i,m in enumerate(rows))
    assert all((rows[i]>>j&1)+(rows[j]>>i&1)==1 for i,j in combinations(range(n),2))
    for i in range(n):d[1<<i][i]=1
    for s in range(1,1<<n):
        for i,h in enumerate(d[s]):
            if not h:continue
            r=rows[i]&~s
            while r:
                t=r&-r;r^=t;d[s|t][t.bit_length()-1]+=h
    return d

def transpose(rows):
    return [sum(1<<i for i,w in enumerate(rows) if w>>j&1) for j in range(len(rows))]

def connected(rows):
    for R in (rows,transpose(rows)):
        s=1
        for _ in R:
            for i,w in enumerate(R):
                if s>>i&1:s|=w
        if s!=(1<<len(rows))-1:return False
    return True

def unused(rows):
    A,B=dp(rows),dp(transpose(rows));full=len(A)-1;n=len(rows)
    return sum(not any(A[s][i] and B[full^s][j] for s in range(1,full))
               for i in range(n) for j in range(n) if rows[i]>>j&1)

def covers(rows):
    d=dp(rows); h=[1]+[sum(r) for r in d[1:]]
    T=[0]*len(h);T[0]=1
    out=[]
    for k in range(1,len(rows)+1):
        U=[0]*len(h)
        for s in range(1,len(h)):
            if s.bit_count()<k:continue
            t=s
            while t:
                U[s]+=h[t]*T[s^t]
                t=(t-1)&s
        T=U
        out.append(T[-1])
    return out

def content(Q):
    return reduce(gcd,Q),reduce(gcd,(factorial(i)*u for i,u in enumerate(Q)))
def polynomial_eval(Q,k):
    from math import comb
    return sum(comb(k,i)*u for i,u in enumerate(Q,1))
def bcase(rows,Q):
    keys=content(Q)
    print('order',len(rows),'H',Q[0],'q2',Q[1],'raw/context',keys,'unused',unused(rows),'rows',rows)
    for p in (3,5,7):
        for u in (1,2,3):
            assert (keys[1]%(p**u)==0)==all(polynomial_eval(Q,i)%(p**u)==0 for i in range(1,p+1))
def wrap(D,I):
    n=len(D);b=n+1
    return [r|((i in I)<<b) for i,r in enumerate(D)]+[(1<<n)-1]+[sum(1<<i for i in range(b) if i not in I)]
for T in ((6,12,24,17,3),(16,5,9,3,14)): bcase(T,covers(T))
for T in ((56,117,129,22,100,76,141,59),(212,137,186,97,170,195,22,72)):
    Q=covers(T);bcase(T,Q)
    print('slots [1..6]:',[polynomial_eval(Q,i) for i in range(1,7)])
R=tuple(sum(1<<j for j in range(7) if (j-i)%7 in (1,2,3))|128 for i in range(7))+(0,)
for T,N in [(wrap(R,{7}),4151),((954,1004,305,1012,322,912,37,852,576,84),4151),
            ((2008,89,1019,1760,1128,131,1696,530,1786,50,678),16079)]:
    Q=covers(T);assert Q[0]==N and connected(T); bcase(T,Q)
# Reproduce the 128 rotationally invariant orientations on F7.
orb=[];seen=set()
for p in combinations(range(7),2):
    if p in seen:continue
    oriented=[p,tuple(2*v%7 for v in p),tuple(4*v%7 for v in p)]
    orb.append(oriented);seen.update(tuple(sorted(e)) for e in oriented)
assert len(orb)==7
hist=Counter()
for bits in product(range(2),repeat=7):
    o=[0]*7
    for w,t in zip(bits,orb):
        for a,b in t:
            if w:a,b=b,a
            o[a]|=1<<b
    H=sum(dp(o)[-1]);assert H%3==0;hist[H]+=1
print('Eisenstein rotation spectrum:', sorted(hist.items()))
# Residual jet law at both small primes, including monic formal sharp examples.
for p in (3,5,7):
    m=p**3
    def f(x):
        a=x
        for j in range(1,p+1):a*=x*x-j*j
        return a
    assert all(f(i)%m==0 for i in range(1,p+1))
    for a in range(p):
        for t in range(p):
            assert (f(a+p*t)-p*p*a*t*(t+1))%m==0
    assert f(p+1)%m
    print(p,'formal parity threshold',2*p+1,'f(p+1) mod p^3',f(p+1)%m)
