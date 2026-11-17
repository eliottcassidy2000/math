from itertools import product
from random import Random
from collections import defaultdict

def F(x):
    y,z=0,1
    for b in map(int,bin(x)[2:]):
        if not b: z+=y+(y%2)
        y=2*y+b
    return z
def rc(r,t):
    return ((1<<t)-1)^int(format(r,'0'+str(t)+'b')[::-1],2)
def H(u,w):
    return 91+5*(F(u)+F(w)-2)+142*(u+w)+256*u*w+2*(w%2)
def ck(x):
    if not x: raise ArithmeticError('certificate failed')
for t,p,q,delta in [(3,(5,2),(3,4),512),(2,(3,0),(1,2),522)]:
    for a in (1,2,3,4,17,256,2**80+5):
        X,Y=[tuple((a<<t)+i for i in x) for x in (p,q)]
        ck((H(*Y)-H(*X),H(*X[::-1])-H(*X),H(*Y[::-1])-H(*Y))==(delta,2,2))
for t in range(2,10):
    Y=1<<t
    groups=defaultdict(set)
    for r in range(1,Y,2):
        for s in range(0,Y,2):
            groups[r+s,rc(r,t)+rc(s,t)].add(256*r*s+5*(F(r)+F(s)))
    g=min(y-x for l in map(sorted,groups.values()) for x,y in zip(l,l[1:]))
    ck(g==(522 if t==2 else 512))
rng=Random(522)
for _ in range(2000):
    n=rng.randrange(2,40)
    r=rng.randrange(1<<n)|1;s=rng.randrange(1<<n)&~1
    A=(r^s)&rng.randrange(1<<n)&~1
    r1,s1=r^A,s^A
    aa=sum(((r>>i&1)-(s>>i&1))<<i for i in range(n) if A>>i&1)
    B=(r^s)&~A
    bb=sum(((r>>i&1)-(s>>i&1))<<i for i in range(n) if B>>i&1)
    f=F(r1)+F(s1)-F(r)-F(s)
    ck(r1*s1-r*s==aa*bb)
    ck(abs(f)<=4*abs(aa*bb))
print('Two physical moves, width-2..9 sharp controls, 2000 trade identities pass.')
