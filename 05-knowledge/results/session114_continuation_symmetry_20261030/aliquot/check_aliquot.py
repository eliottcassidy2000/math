"""Exact finite controls for ALIQUOT_CONTINUATION.md; standard library only."""
import json
from math import gcd, isqrt

gates=0
def check(condition, message):
    global gates
    gates += 1
    if not condition:
        raise RuntimeError(message)

def sigma(n):
    total=0
    for d in range(1,isqrt(n)+1):
        if n%d==0:
            total+=d
            if d*d!=n: total+=n//d
    return total

def kernel(n):
    while n%2==0: n//=2
    out=1
    p=3
    while p*p<=n:
        parity=0
        while n%p==0:
            n//=p
            parity^=1
        if parity: out*=p
        p+=2
    return out*n

def obs(n): return (sigma(n)-n)%2

def formula(n,m): return int(kernel(n)==kernel(m)) ^ ((n%2)*(m%2))

for n in range(1,251):
    for m in range(1,121):
        check(obs(n*m)==formula(n,m),f'formula {n},{m}')

# Complete separation witnesses for all small pairs, using odd multipliers.
for n in range(1,101):
    for h in range(n+1,101):
        if (n%2,kernel(n))==(h%2,kernel(h)):
            check(all(formula(n,m)==formula(h,m) for m in range(1,42,2)), 'equal classes')
        else:
            m=101 if n%2 != h%2 else kernel(n)
            check(obs(n*m)!=obs(h*m),f'separation {n},{h},{m}')

# Exact source-derived fixed-alphabet counts, including r=0 and r=1 collisions.
primes=[3,5,7,11]
counts=[]
for r in range(5):
    P=primes[:r]
    multipliers=[]
    for mask in range(1<<r):
        m=1
        for i,p in enumerate(P):
            if mask>>i&1: m*=p
        multipliers.append(m)
    sources=[m*f for m in multipliers for f in (1,2)]+[13,26]
    functions={tuple(obs(n*m) for m in multipliers) for n in sources}
    expected=2 if r==0 else (4 if r==1 else 2+(1<<(r+1)))
    check(len(functions)==expected,f'fixed alphabet {r}')
    counts.append({'r':r,'functions':len(functions)})

exponents=[3,5,7,11,13,17,19]
witnesses=[]
for p in exponents:
    n=1<<p
    y=(1<<p)-1
    K=kernel(y)
    check(sigma(n)-n==y,'successor')
    check((n%2,kernel(n))==(0,1),'common class')
    check(y%8==7 and isqrt(y)**2!=y and K>1,'nonsquare')
    check(isqrt(y*K)**2==y*K,'square completion')
    witnesses.append({'p':p,'source':n,'successor':y,'kernel':K})
for a in witnesses:
    for b in witnesses:
        if a==b: continue
        check(gcd(a['successor'],b['successor'])==1,'coprime successors')
        check(gcd(a['kernel'],b['kernel'])==1,'coprime kernels')
        # Independent parity control from square testing avoids factoring a huge product.
        for w,expected in ((a,0),(b,1)):
            z=w['successor']*a['kernel']
            actual=1-int(isqrt(z)**2==z)
            check(actual==expected,'context separates')

# Original 25/49 hostile is retained as a small direct divisor control.
check((25%2,kernel(25))==(49%2,kernel(49)),'inherited collision')
check(obs(sigma(25)-25)!=obs(sigma(49)-49),'aliquot separates')
print(json.dumps({'status':'PASS','gates':gates,'formula_universe':{'n':250,'m':120},
 'pair_universe':100,'fixed_alphabet_counts':counts,'witnesses':witnesses},indent=2))
