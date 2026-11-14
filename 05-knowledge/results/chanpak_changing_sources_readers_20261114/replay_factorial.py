"""Subset-path DP controlling literal sharpness and the source change."""
from functools import cache
from math import factorial
from random import Random

def check(flag):
    if not flag:
        raise AssertionError('control failed')

def profile(rows):
    n = len(rows)
    full = (1 << n)-1
    for i in range(n):
        check((rows[i] >> i)&1 == 0 and rows[i]<=full)
        for j in range(i):
            check(((rows[i] >> j)&1)+((rows[j] >> i)&1)==1)

    @cache
    def trail(s,v):
        rem=s^(1<<v)
        return sum(trail(rem,w) for w in range(n) if rows[v]&rem&(1<<w)) if rem else 1

    @cache
    def p(s,k):
        if k==0: return int(s==0)
        if k==1: return sum(trail(s,v) for v in range(n) if s&(1<<v)) if s else 1
        z=s
        v=0
        while True:
            v += p(z,1)*p(s^z,k-1)
            if not z: return v
            z=(z-1)&s
    return lambda k:p(full,k)

def pair(a,b):
    # Singleton 0 beats A; A beats B; B beats singleton.
    return ((1<<(1+len(a)))-2,) + tuple(
            (x<<1)|(((1<<len(b))-1)<<(len(a)+1)) for x in a
        ) + tuple((x<<(len(a)+1))|1 for x in b)

def op(g,a,t):
    if t==0: return 0
    return t*(g(t)*a(t) + (a(t+1)-a(t))*(g(t+1)-g(t))
              + (a(t-1)-a(t))*(g(t-1)-g(t)))

def random_t(n,rng):
    out=[0]*n
    for x in range(n):
        for y in range(x):
            a,b=(x,y) if rng.randrange(2) else (y,x)
            out[a]|=1<<b
    return tuple(out)

controls = [
    ((436,225,18,279,450,348,141,300,70),(388,21,328,211,421,143,179,262,106),
      2269,37124,311991,311751,39393,807333,806853),
    ((388,37,352,55,775,849,411,574,138,335),(602,520,163,916,454,731,142,259,103,468),
      9359,176552,1696167,1695447,185911,4265735,4264295),
]
for A,B,h,q,u,v,w,x,y in controls:
    for source, p3, w2 in ((A,u,x),(B,v,y)):
        p=profile(source)
        check([p(t) for t in (1,2,3)]==[h,q,p3])
        bot=pair((0,),source)
        check(profile(bot)(1)==w)
        check(profile(pair((0,),bot))(1)==w2)
    check(x-y==2*factorial(2)*factorial(5)*(3 if len(A)%2==0 else 1))

rng=Random(6769)
for _ in range(48):
    a,b=(random_t(rng.randrange(1,5),rng) for _ in range(2))
    pa,pb,pc=profile(a),profile(b),profile(pair(a,b))
    for t in (1,2,3):
        check(pc(t)==op(pb,pa,t))
for r in range(1,22):
    for t in range(50):
        v=1
        for j in range(-r,r+1):
            v*=t-j
        check(v%factorial(2*r+1)==0)
        check((t*v)%(factorial(2*r+1)*(r+1))==0)
for (a,b,*_), (s,t,x,y) in zip(controls,[(288,66,8931,8933),(320,5,37689,37691)]):
    check(s.bit_count()==t.bit_count()==2)
    for B,M,H in ((a,s,x),(b,t,y)):
        n=len(B)
        glued=tuple(row|((1<<n) if not (M>>v)&1 else 0)
                    for v,row in enumerate(B))+(M,)
        check(profile(glued)(1)==H)
print('Literal sharp r=2 pairs (both parities), 144 arbitrary-module checks, and factorial arithmetic controls and selective reader bypass pass.')
