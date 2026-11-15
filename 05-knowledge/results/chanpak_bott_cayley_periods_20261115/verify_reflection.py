"""Replay exact words on physical addresses and formal residue states."""
from itertools import product
from math import gcd
from random import Random

def need(x):
    if not x: raise RuntimeError("identity failed")

def p(x):
    s=0
    while x:
        x,t=divmod(x,2)
        if not t:s+=x+(x&1)
    return s

def raw(x,y):
    return 432*x*y+129*x+258*y+15*(x&1)+30*(y&1)+45*(p(x)+p(y))+135

def state(x,y):
    return x+1,y,p(x)+p(y)-y+(y&1),y&1

def tap(s,t):
    a,b,k,r=s
    return 2*a-t,2*b+t,k+t*(a-b-2*r),t

def read(s):
    a,b,k,r=s
    return 432*a*b+129*(a-b)+45*k+21-30*r

def probes(s):
    a,b,k,r=s
    need(r==0)
    return (1728*a*b+1122*a-303*b+45*k-123,
            6912*a*b+3972*a-561*b+45*k-237)

def reflection(s):
    a,b,k,r=s
    return -b,-a,k,r

def signatures(s,p):
    s1=tap(s,0)
    s2=tap(s1,0)
    t=tap(tap(s,1),0)
    R=(read(s),read(s1),read(s2),read(t),read(tap(t,0)))
    return tuple(v%p for v in R),tuple(v%p for v in R+probes(s))

rng=Random(20261115)
for _ in range(400):
    x=2*rng.randrange(1,1<<160)+1
    y=2*rng.randrange(1,1<<160)
    s=state(x,y)
    need(probes(s)==(raw(2*x+1,2*y+1),raw(4*x+3,4*y+2)))
    T=x+y+1
    need(tuple(a-b for a,b in zip(probes(s),probes(reflection(s))))==(819*T,3411*T))
    for _ in range(30):
        t=rng.randrange(2)
        x,y=2*x+1-t,2*y+t
        need(reflection(tap(s,t))==tap(reflection(s),t))
        s=tap(s,t)
        need(s==state(x,y) and raw(x,y)==read(s)==read(reflection(s)))
for f in (7,11,13,43):
    seen={};full=set()
    for a,b,k in product(range(f),repeat=3):
        s=a,b,k,0
        key,u=signatures(s,f)
        prof=(a*b%f,(a-b)%f,k)
        need(seen.setdefault(key,prof)==prof)
        need(u not in full);full.add(u)
need(gcd(819,3411)==9)
print("400 terminal probe comparisons and 12000 physical appends; full field signatures at 7,11,13,43")
